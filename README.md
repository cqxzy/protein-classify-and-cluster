# protein-classify-and-cluster
A Python pipeline for protein structure clustering (PDBs) using all-vs-all TMalign, hierarchical analysis, and secondary structure classification.
Also, you may use it to reproduce the protein structure clustering method (TMalign + Hierarchical) from the Nature paper "Hallucination of closed repeat proteins".


# **Protein Structure Clustering Pipeline User Guide**

## **1\. Introduction**

This Python script provides an automated workflow for performing secondary structure analysis and hierarchical clustering on protein PDB files within a directory. This workflow reproduces the core clustering method described in the *Nature* paper ("Hallucination of closed repeat proteins containing central pockets").
you can get pictures like:
![the example](images/your-image.png)
**Process Overview:**

1. **Secondary Structure Analysis:** Traverse all PDB entries, compute their α-helix and β-sheet content, and perform coarse classification ("alpha", "beta", "alpha\_beta")  
2. **Similarity Matrix:** Perform all-vs-all structural alignment using TMalign to generate a TM-score similarity matrix.  
3. **Automatic Clustering:** Automatically search for the optimal number of clusters (K), determining the best K value by minimizing the proportion of "orphan clusters".  
4. **Result Output:** Generate a list of cluster members, a list of representative structures, and a final dendrogram colored by secondary structure.

## **2. Dependencies and Installation**


### **2.1 External Programs**


This script relies on a critical external bioinformatics tool:

* **TMalign**: **(Required)** Used to compute structural similarity. You must provide the full path to its executable in the --tmalign parameter.  
* **mkdssp**: **(Optional)** Only required when the ss parameter is set to dssp. The script's default phi\_psi mode has no such dependency.

### **2.2 Python Environment**

It is recommended to use conda to create an isolated environment (for convenience in experimentation, I set up an environment to perform these tasks; however, in practice, you only need to have these packages installed to proceed directly. The step of setting up an environment can be skipped).

\# 1\. Create and activate the environment  
conda create \-n protein\_cluster python=3.10  
conda activate protein\_cluster

\# 2\. Install core packages  
conda install \-c anaconda pandas numpy scipy scikit-learn  
conda install \-c conda-forge matplotlib  
conda install \-c bioconda biopython

## **3. Command Line Usage**

### **3.1 Basic Commands**

\# 1\. (Optional) Create an output directory  
mkdir \-p ./results

\# 2\. Run the script  
python classify_and_cluster_proteins.py \
    --inp "./my_pdbs" \
    --out "./results" \
    --tmalign "/path/to/your/TMalign" \
    --jobs 16

For example, here is an example I have used before:  
python \~/proteincluster/classify_and_cluster_proteins.py \
--inp "./final_sample_for_dendrogram" \
--out "./output_final_super_dendrogram" \ 
--tmalign "$(which TMalign)" \  
--jobs 8

All that's needed is to modify the input and output paths. If tmalign has no other requirements, it can match my example. Of course, lines 378-391 in the code contain additional configurable values. You can either directly modify the values corresponding to "default" in the code or set them via the command line using the same format as above.
In Markdown format, the command line will contain numerous backslashes ("\"). When using it in practice, avoid copying the example above directly, as it may include extra characters.

### **3.2 Performance Warning: O(N²) Complexity**

**Note:** The core operation of this script is computing an all-vs-all TMalign matrix. Its computational complexity is O(N²), where N is the number of input PDBs.

* **N = 100**: ~4,950 comparisons (minutes)  
* **N = 1,000**: ~500,000 comparisons (hours)  
* **N = 10,000**: ~50,000,000 comparisons (days to weeks)  

**Do not** directly input tens of thousands of PDB files without preprocessing (e.g., using Foldseek). For large datasets, first perform rapid clustering (e.g., Foldseek) to extract representatives. Then use these **representatives** (N < 1500) as input for this script for detailed clustering and visualization.


### **3.3 foldseek Acceleration Workflow**

Using the TMalign method is time-consuming and generates numerous files. It is recommended to first use foldseek for rapid screening of many proteins. Employ foldseek's easycluster mode.

You can use this to determine how many cores to use. For laptops, set this to avoid excessive resource consumption:
THREADS=$(($(nproc)-2))

Next, run this command:

foldseek easy-cluster input pdbs sample/my _clusters tmp fs \

--tmscore-threshold 0.6 \

--alignment-type 1 \

-s 9 \

--threads $THREADS \

--gpu 1

These two parameters represent the sensitivity level and the filtering model used:
\--alignment-type 1 \\ This is similar to the model category used; 2 is faster with minimal difference in results
\--s 9 \\ Higher values increase sensitivity but slow down processing
\--gpu 1 \\ Based on personal experience, this option doesn't seem to offer significant improvement

For more details, refer to the official foldseek GitHub documentation.

Alternatively, visit for advice: https://zread.ai/steineggerlab/foldseek 


## **4. Parameter Description**

#### **Core Parameters**

* --inp: **Input directory containing all .pdb files.  
* --out: **Directory storing all output results (CSV, JSON, PNG files).  
* --tmalign: **Full path** to the TMalign executable.  
* --jobs:  
  * ** **Number of CPU cores used for parallel execution of TMalign and secondary structure analysis.  
  * **Default value:** 4.  

#### **Secondary Structure (SS) Parameters**

* --ss:  
  * **Description:** Method for calculating secondary structure.  
  * phi_psi: (default) Uses Biopython's built-in Phi/Psi angle calculation, which is **fast** and **does not require mkdssp**.  
  * dssp: Uses the mkdssp program for more accurate results, but is **slower** and has **additional dependencies**.  
* --alpha_thr, --beta_thr:  
  * **Description:** Percentage thresholds for determining "all-alpha" or "all-beta" topologies.  
  * **Default:** 0.8 (i.e., 80%).  
  
#### **Clustering Tuning Parameters**

* --k_min, --k_max:  
  * **Description:** Minimum and maximum range for automatically searching the "optimal K value".  
  * **Default:** k_min=2, k_max=16.  
* --mean_thr, --std_ceil:  
  * **Description:** Metrics used to filter the "optimal K value". Defines a "good" cluster (average TM-score within the cluster > mean_thr, and standard deviation of TM-scores within the cluster < std_ceil).  
  * **Default values:** 0.7, 0.1.

## **5\. Output File Analysis**

After successful script execution, you will find the following files in the \--out directory:

* dendrogram.jpg:
  **Core visualization result**. A hierarchical clustering dendrogram with leaf labels automatically colored based on secondary structure classification. However, this generated image does not include protein structure diagrams, differing from the reproduced paper "Hallucination of closed repeat proteins containing central pockets." You can manually remove text from the image and attach protein images. Additionally, when classifying too many proteins, the image becomes very large. Using 20 proteins for visualization is reasonable; more proteins result in an extremely long image.   
* secondary_structure_summary.csv:  
  Contains the secondary structure content (helix_frac, beta_frac) and rough classification (topology) for each PDB.  
* cluster_selection.json:  
Records statistics from the automatic K-value search and the final "best K-value" (chosen_k).  
* clusters.csv:  
The final clustering assignment table. Lists each protein and its assigned cluster_id.  
* representatives.csv:  
  **(Important)** Lists one "representative" protein (i.e., the medoid, or central point of the cluster) for each cluster_id. You can first run a round on a moderate number of proteins, then use the proteins listed in this file for another round. This approach can select a reasonable number of proteins that also yield visually appropriate images.


## **6. How to Modify and Extend Functionality**

### **6.1 Improving Secondary Structure Classification**

The current classification (`ss_summary` function) only considers the **proportion** of secondary structures (e.g., `helix_frac > 0.8`), not their **sequence**. To implement topology-based classification like "alpha-beta-alpha," the following is required:

1. **Positioning Functions:** Open classify_and_cluster_proteins.py and locate the helper functions _ss_from_dssp or _ss_from_phi_psi.  
2. **Modify Logic:**  
   * Ensure these functions return not only the frac (proportion) but also a **string representing the topological order**.  
   * For example, if DSSP yields CCCHHHHHHEEEECCCCHHHH, you can implement a run-length encoding function to convert it to H5E4H4 (ignoring C, coil).  
3. **Save topology:** In the ss_summary function, store this new topology string (H5E4H4) in the df\['topology'\] column.  
4. **Modify plotting:** In the plot_dendrogram function, modify the ss_map dictionary to recognize your new topology string and map it to specific colors (e.g., all H-prefixed strings as red, E-prefixed as blue, etc.).

### **6.2 Replacing/Modifying Similarity Metrics**

* **Location:** tm_matrix and tmalign_pair functions.
* **Approach:** The script currently hardcodes TMalign. You can modify the tmalign_pair function to support calling US-align or other alignment tools and parsing their outputs.

### **6.3 Replace/Modify Clustering Algorithm**

* **Location:** select_k_from_grid and run_final_clustering functions.
* **Approach:** The script currently hardcodes AgglomerativeClustering and requires a complex K-value search. You can replace it with HDBSCAN (a more modern clustering algorithm that doesn't require a predefined K value), which would simplify or eliminate the need for the "automatic K value selection" section.

Given the author's limited expertise and rushed development timeline, this code and related documentation inevitably contain omissions and shortcomings.

This work was carried out in the group of Professor Li Zhe at Southern University of Science and Technology.

Xing Zhenyan
If you want, you may contact by: zyxing05@gmail.com

November 13, 2025
Hong Kong
