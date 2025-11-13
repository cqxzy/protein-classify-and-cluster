#!/usr/bin/env python3
"""
classify_and_cluster_proteins.py (phi/psi default)

Replicates the clustering/visualization logic behind Fig.2 and Extended Data Fig.4
from "Hallucination of closed repeat proteins containing central pockets" for your
own set of protein structures.

Pipeline
1) Read all .pdb files from an input folder.
2) Assign per-residue secondary structure. **Default: pure-Python phi/psi heuristic**
   (robust on WSL). Optional: mkdssp via --ss dssp/auto if your mkdssp works.
3) Summarize helix/beta/coil fractions and classify topology: all-alpha / all-beta / alpha-beta.
4) Compute an all-by-all TM-score matrix using TM-align (average the two TM-scores reported
   by TM-align per pair to mimic the paper's length-normalized averaging).
5) Perform hierarchical clustering and plot a dendrogram with leaf labels colored by topology.
6) Grid-search the number of clusters with AgglomerativeClustering over the precomputed
   TM-distance matrix (1 - TMscore) to minimize singleton clusters while keeping high
   mean intra-cluster TMscore and low std. Save labels & quality metrics.
7) Pick medoid representatives per cluster; optionally export a PyMOL batch script to color
   helix/sheet/loop for Extended Data Fig.4-like panels.

Outputs (written to --out):
- secondary_structure_summary.csv : per-structure H/E/C fractions and topology label
- tmatrix.npy                     : NxN symmetric matrix of averaged TM-scores
- labels.csv                      : file list (order used for matrices/plots)
- dendrogram.png                  : dendrogram with colored leaf labels
- clusters.csv                    : cluster_id for each structure at the selected K
- cluster_selection.json          : selection stats across K and the chosen K
- representatives.csv             : medoid representative per cluster
- color_by_ss.pml                 : optional PyMOL script to color helix/sheet/loop

Usage
------
python classify_and_cluster_proteins.py \
  --inp ./pdbs \
  --out ./out \
  --tmalign /path/to/TMalign \
  --ss phi_psi \
  --jobs 8

Notes
- Requires: Python 3.9+, Biopython, numpy, pandas, scipy, scikit-learn, matplotlib
- External tools: TM-align executable (path via --tmalign). **DSSP is optional** and only used
  if you pass --ss dssp/auto and have a working mkdssp.
- If you have .cif files, convert to .pdb first (e.g., with pdb-tools or PyMOL)
"""

from __future__ import annotations
import argparse
import concurrent.futures as cf
import json
import math
import re
import subprocess
import tempfile, shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from Bio.PDB import PDBParser, PPBuilder
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
#1105修改的下面两行，为了没有图形界面
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ------------------------------
# Secondary structure utilities
# ------------------------------

SS_HELIX_SET = {"H", "G", "I"}  # alpha, 3-10, pi
SS_BETA_SET  = {"E", "B"}       # beta strand/bridge


def run_dssp(pdb_path: Path) -> Dict[str, float]:
    """Run DSSP via mkdssp (subprocess) with modern/legacy CLI fallback.
    Only used when --ss dssp/auto. Copies to /tmp to avoid WSL path quirks.
    """
    if not pdb_path.is_file():
        raise RuntimeError(f"Input is not a file: {pdb_path}")
    if pdb_path.stat().st_size == 0:
        raise RuntimeError(f"Input PDB is empty: {pdb_path}")

    with tempfile.TemporaryDirectory() as tmpd:
        tmp_in = Path(tmpd) / pdb_path.name
        shutil.copy2(pdb_path, tmp_in)
        tmp_out = Path(tmpd) / (pdb_path.stem + ".dssp")

        dssp_bin = shutil.which("mkdssp") or "mkdssp"
        cmds = [
            [dssp_bin, "-i", str(tmp_in), "-o", str(tmp_out)],  # modern
            [dssp_bin, str(tmp_in), str(tmp_out)],              # legacy
        ]
        last_err = None
        for cmd in cmds:
            try:
                subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                last_err = None
                break
            except subprocess.CalledProcessError as e:
                last_err = e
        if last_err is not None:
            raise RuntimeError(f"mkdssp failed on {pdb_path.name}: {(last_err.stderr or last_err.stdout or '').strip()}")

        if not tmp_out.exists() or tmp_out.stat().st_size == 0:
            raise RuntimeError(f"mkdssp produced no output for {pdb_path.name}")

        # Parse DSSP output text
        helix = beta = total = 0
        started = False
        with open(tmp_out, "r", errors="ignore") as fh:
            for line in fh:
                if not started:
                    if line.lstrip().startswith('#') and 'RESIDUE AA STRUCTURE' in line:
                        started = True
                    continue
                if not line.strip():
                    continue
                ch = line[16] if len(line) > 17 else ' '
                if ch not in ('H','G','I','E','B','T','S',' '):
                    window = line[14:25]
                    ch = next((c for c in window if c in ('H','G','I','E','B','T','S',' ')), ' ')
                if ch in ('H','G','I'):
                    helix += 1
                elif ch in ('E','B'):
                    beta += 1
                total += 1

        if total == 0:
            return {"helix_frac": 0.0, "beta_frac": 0.0, "coil_frac": 1.0}

        helix_frac = helix / total
        beta_frac = beta / total
        coil_frac = max(0.0, 1.0 - helix_frac - beta_frac)
        return {"helix_frac": helix_frac, "beta_frac": beta_frac, "coil_frac": coil_frac}


def estimate_ss_by_phipsi(pdb_path: Path) -> Dict[str, float]:
    """Assign secondary structure via phi/psi heuristic (pure-Python).
    Regions (degrees):
      alpha-helix  : phi in [-120,-30], psi in [-80,  -5]
      beta-strand  : phi in [-180,-100], psi in [ 70, 180]
    Others are counted as coil. Multi-chain supported (aggregate residues).
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = structure[0]
    ppb = PPBuilder()

    helix = beta = total = 0

    for chain in model:
        polypeptides = ppb.build_peptides(chain)
        for pp in polypeptides:
            phipsi = pp.get_phi_psi_list()
            for phi, psi in phipsi:
                if phi is None or psi is None:
                    continue
                phi_deg = math.degrees(phi)
                psi_deg = math.degrees(psi)
                if (-120 <= phi_deg <= -30) and (-80 <= psi_deg <= -5):
                    helix += 1
                elif (-180 <= phi_deg <= -100) and (70 <= psi_deg <= 180):
                    beta += 1
                total += 1

    if total == 0:
        return {"helix_frac": 0.0, "beta_frac": 0.0, "coil_frac": 1.0}
    helix_frac = helix / total
    beta_frac = beta / total
    coil_frac = max(0.0, 1.0 - helix_frac - beta_frac)
    return {"helix_frac": helix_frac, "beta_frac": beta_frac, "coil_frac": coil_frac}


def classify_topology(helix_frac: float, beta_frac: float,
                      alpha_thr: float = 0.60,
                      beta_thr: float = 0.45,
                      cross_ceil: float = 0.15) -> str:
    """Heuristic topology classification similar to the paper's alpha/beta/mixed panels."""
    if helix_frac >= alpha_thr and beta_frac <= cross_ceil:
        return "alpha"
    if beta_frac >= beta_thr and helix_frac <= cross_ceil:
        return "beta"
    return "alpha_beta"

# ------------------------------
# TM-align utilities
# ------------------------------

TM_LINE_RE = re.compile(r"TM-score\s*=\s*([0-9]*\.?[0-9]+)")


def tmalign_pair(tmalign_bin: str, pdb_a: Path, pdb_b: Path) -> float:
    """Run TM-align on two PDBs and return the average of the two reported TM-scores."""
    cmd = [tmalign_bin, str(pdb_a), str(pdb_b)]
    try:
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"TM-align failed on {pdb_a.name} vs {pdb_b.name}:\n{e.output}")

    scores = [float(m.group(1)) for m in TM_LINE_RE.finditer(out)]
    if len(scores) < 2:
        raise RuntimeError(f"Could not parse two TM-scores from TM-align output for {pdb_a.name} vs {pdb_b.name}.")
    return float(np.mean(scores[:2]))


def tm_matrix(pdbs: List[Path], tmalign_bin: str, jobs: int = 4) -> np.ndarray:
    """Compute an NxN symmetric matrix of averaged TM-scores for PDB list."""
    n = len(pdbs)
    tm = np.eye(n, dtype=float)

    def work(pair):
        i, j = pair
        return i, j, tmalign_pair(tmalign_bin, pdbs[i], pdbs[j])

    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]

    with cf.ThreadPoolExecutor(max_workers=jobs) as ex:
        for i, j, s in ex.map(work, pairs):
            tm[i, j] = s
            tm[j, i] = s
    return tm

# ------------------------------
# Clustering utilities
# ------------------------------

@dataclass
class ClusterStats:
    k: int
    prop_singletons: float
    mean_intra_tm: float
    std_intra_tm: float


def evaluate_clustering(labels: np.ndarray, tm: np.ndarray) -> ClusterStats:
    unique, counts = np.unique(labels, return_counts=True)
    singletons = (counts == 1).sum()
    prop_single = singletons / len(unique)

    intra_scores = []
    for lab in unique:
        members = np.where(labels == lab)[0]
        if len(members) <= 1:
            continue
        sub = tm[np.ix_(members, members)]
        iu = np.triu_indices_from(sub, k=1)
        if iu[0].size:
            intra_scores.extend(sub[iu].tolist())
    if len(intra_scores) == 0:
        mean_tm = 0.0
        std_tm = 0.0
    else:
        mean_tm = float(np.mean(intra_scores))
        std_tm = float(np.std(intra_scores))
    return ClusterStats(k=len(unique), prop_singletons=prop_single, mean_intra_tm=mean_tm, std_intra_tm=std_tm)


def select_k_from_grid(tm: np.ndarray, k_min: int, k_max: int,
                       mean_thr: float = 0.85, std_ceil: float = 0.12,
                       linkage_method: str = 'average') -> Tuple[int, ClusterStats, Dict[int, ClusterStats]]:
    """Grid-search k; prefer minimal singleton proportion subject to mean/std constraints."""
    n = tm.shape[0]
    dist = 1.0 - tm

    all_stats: Dict[int, ClusterStats] = {}
    best_k = None
    best_tuple = None  # (obj, -mean, std, k)

    for k in range(k_min, min(k_max, n) + 1):
        ac = AgglomerativeClustering(n_clusters=k, metric='precomputed', linkage=linkage_method)
        labels = ac.fit_predict(dist)
        stats = evaluate_clustering(labels, tm)
        all_stats[k] = stats

        meets = (stats.mean_intra_tm >= mean_thr) and (stats.std_intra_tm <= std_ceil)
        obj = stats.prop_singletons
        tuple_key = (obj, -stats.mean_intra_tm, stats.std_intra_tm, k)

        if meets:
            if best_k is None or tuple_key < best_tuple:
                best_k, best_tuple = k, tuple_key
        else:
            if best_k is None:
                if (best_tuple is None) or (tuple_key < best_tuple):
                    best_k, best_tuple = k, tuple_key

    return best_k, all_stats[best_k], all_stats


def medoid_indices_per_cluster(labels: np.ndarray, tm: np.ndarray) -> Dict[int, int]:
    reps: Dict[int, int] = {}
    for lab in np.unique(labels):
        idx = np.where(labels == lab)[0]
        if len(idx) == 1:
            reps[int(lab)] = int(idx[0])
            continue
        sub = tm[np.ix_(idx, idx)]
        mean_to_others = sub.mean(axis=1)
        reps[int(lab)] = int(idx[int(np.argmax(mean_to_others))])
    return reps

# ------------------------------
# Plotting
# ------------------------------

COLOR_MAP = {
    "alpha": "#008080",       # teal (to match paper)
    "beta": "#FF00FF",        # magenta
    "alpha_beta": "#003366",  # dark blue
}


def plot_dendrogram(tm: np.ndarray, names: List[str], topo_labels: List[str], out_png: Path):
    dist = 1.0 - tm
    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method='average')

    fig, ax = plt.subplots(figsize=(max(10, len(names) * 0.25), 8))
    dendrogram(Z, labels=names, leaf_rotation=90, leaf_font_size=8, ax=ax)

    label_to_topo = dict(zip(names, topo_labels))
    for lbl in ax.get_xmajorticklabels():
        text = lbl.get_text()
        topo = label_to_topo.get(text, "alpha_beta")
        lbl.set_color(COLOR_MAP.get(topo, "black"))

    ax.set_ylabel("Relative structural dissimilarity (1 - TMscore)")
    ax.set_title("Hierarchical clustering of proteins (leaf labels colored by topology)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

# ------------------------------
# PyMOL coloring script (optional)
# ------------------------------

PYMOL_SCRIPT = r"""
# color_by_ss.pml : color helix/sheet/loop for all loaded objects
# Usage in PyMOL:
#   pymol -qc color_by_ss.pml -- pdb_dir=out/representatives
# Expects PDBs in 'out/representatives'.

from pymol import cmd
import os

# read argument
pdb_dir = cmd.get("pdb_dir") or "representatives"

# load all pdbs in folder
for f in os.listdir(pdb_dir):
    if f.lower().endswith('.pdb'):
        obj = os.path.splitext(f)[0]
        cmd.load(os.path.join(pdb_dir, f), obj)
        cmd.show("cartoon", obj)
        cmd.util.cbaw(obj)  # base colors, then override
        cmd.color("teal", f"({obj} and ss h+g+i)")
        cmd.color("magenta", f"({obj} and ss s+b)")
        cmd.color("marine", f"({obj} and not (ss h+g+i or ss s+b))")
        cmd.orient(obj)
        cmd.ray(1200, 900)
        cmd.png(os.path.join(pdb_dir, obj + "_ss.png"), dpi=300)

cmd.quit()
"""

# ------------------------------
# Main
# ------------------------------

def main():
    ap = argparse.ArgumentParser(description="Classify proteins by secondary structure and cluster by TM-score.")
    ap.add_argument('--inp', required=True, help='Folder containing input .pdb files (monomers).')
    ap.add_argument('--out', required=True, help='Output folder.')
    ap.add_argument('--tmalign', required=True, help='Path to TM-align executable.')
    ap.add_argument('--jobs', type=int, default=4, help='Parallel TM-align workers (default: 4).')
    ap.add_argument('--ss', choices=['auto','phi_psi','dssp'], default='phi_psi',
                    help='Secondary-structure assignment: phi_psi (pure Python), dssp (mkdssp), or auto (try dssp then fall back). Default: phi_psi.')
    ap.add_argument('--alpha_thr', type=float, default=0.60, help='Threshold for all-alpha classification (default: 0.60).')
    ap.add_argument('--beta_thr', type=float, default=0.45, help='Threshold for all-beta classification (default: 0.45).')
    ap.add_argument('--cross_ceil', type=float, default=0.15, help='Ceiling for cross-secondary when classifying (default: 0.15).')
    ap.add_argument('--k_min', type=int, default=20, help='Min clusters to try (default: 20).')
    ap.add_argument('--k_max', type=int, default=400, help='Max clusters to try (default: 400).')
    ap.add_argument('--mean_thr', type=float, default=0.85, help='Min mean intra-cluster TMscore (default: 0.85).')
    ap.add_argument('--std_ceil', type=float, default=0.12, help='Max std of intra-cluster TMscore (default: 0.12).')
    args = ap.parse_args()

    inp = Path(args.inp)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    pdbs = sorted([p for p in inp.glob('*.pdb')])
    if len(pdbs) < 2:
        raise SystemExit(f"Need at least 2 PDBs in {inp}")

    # 1) Secondary structure + topology
    ss_rows = []
    warn_log = []
    mode = getattr(args, 'ss', 'phi_psi')
    print(f"[SS mode] {mode}")

    for p in pdbs:
        if mode == 'phi_psi':
            fracs = estimate_ss_by_phipsi(p)
            method = 'phi_psi'
        elif mode == 'dssp':
            fracs = run_dssp(p)
            method = 'dssp'
        else:  # auto
            try:
                fracs = run_dssp(p)
                method = 'dssp'
            except Exception as e:
                warn = f"DSSP failed on {p.name}: {e}. Falling back to phi/psi heuristic."
                print(warn)
                warn_log.append(warn)
                fracs = estimate_ss_by_phipsi(p)
                method = 'phi_psi'

        topo = classify_topology(fracs['helix_frac'], fracs['beta_frac'],
                                 alpha_thr=args.alpha_thr, beta_thr=args.beta_thr, cross_ceil=args.cross_ceil)
        ss_rows.append({
            'name': p.stem,
            'path': str(p.resolve()),
            **fracs,
            'topology': topo,
            'ss_method': method,
        })

    ss_df = pd.DataFrame(ss_rows)
    if warn_log:
        with open(out / 'secondary_structure_warnings.log', 'w') as f:
            f.write('\n'.join(warn_log))
    ss_df.to_csv(out / 'secondary_structure_summary.csv', index=False)

    # 2) TM-score matrix
    tm = tm_matrix(pdbs, args.tmalign, jobs=args.jobs)
    np.save(out / 'tmatrix.npy', tm)

    names = [p.stem for p in pdbs]
    pd.DataFrame({'name': names, 'path': [str(p.resolve()) for p in pdbs]}).to_csv(out / 'labels.csv', index=False)

    # 3) Dendrogram
    topo_labels = [ss_df.set_index('name').loc[n, 'topology'] for n in names]
    plot_dendrogram(tm, names, topo_labels, out / 'dendrogram.png')

    # 4) Cluster selection via constraints
    n = len(pdbs)
    k_max = min(args.k_max, max(args.k_min, n))
    k_min = min(args.k_min, k_max)

    best_k, best_stats, all_stats = select_k_from_grid(tm, k_min, k_max, args.mean_thr, args.std_ceil)

    # Fit final clustering for best_k
    dist = 1.0 - tm
    ac = AgglomerativeClustering(n_clusters=best_k, metric='precomputed', linkage='average')
    labels = ac.fit_predict(dist)

    # Save cluster labels
    pd.DataFrame({'name': names, 'cluster_id': labels, 'topology': topo_labels}).to_csv(out / 'clusters.csv', index=False)

    # Save stats grid and choice
    stats_json = {
        'chosen_k': best_k,
        'chosen_stats': best_stats.__dict__,
        'all_stats': {int(k): s.__dict__ for k, s in all_stats.items()},
    }
    with open(out / 'cluster_selection.json', 'w') as f:
        json.dump(stats_json, f, indent=2)

    # 5) Representatives (medoids)
    reps = medoid_indices_per_cluster(labels, tm)
    rep_rows = []
    rep_dir = out / 'representatives'
    rep_dir.mkdir(exist_ok=True)
    for cl, idx in reps.items():
        name = names[idx]
        src = Path(ss_df.set_index('name').loc[name, 'path'])
        dst = rep_dir / f"cluster{cl:03d}_{name}.pdb"
        dst.write_bytes(Path(src).read_bytes())
        rep_rows.append({'cluster_id': cl, 'name': name, 'path': str(dst.resolve())})
    pd.DataFrame(rep_rows).to_csv(out / 'representatives.csv', index=False)

    # 6) Write optional PyMOL coloring script
    with open(out / 'color_by_ss.pml', 'w') as f:
        f.write(PYMOL_SCRIPT)

    print("Done.")
    for fn in [
        'secondary_structure_summary.csv', 'tmatrix.npy', 'dendrogram.png',
        'clusters.csv', 'cluster_selection.json', 'representatives.csv', 'color_by_ss.pml']:
        print(f"Wrote: {out / fn}")


if __name__ == '__main__':
    main()
