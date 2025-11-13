# **蛋白质结构聚类流程 (Protein Structure Clustering Pipeline) \- 用户指南**

## **1\. 简介**

本 Python 脚本 提供了一个自动化的流程，用于对一个目录下的蛋白质 PDB 文件进行二级结构分析和层次聚类。本流程复现了 *Nature* 论文 ("Hallucination of closed repeat proteins containing central pockets”) 中的核心聚类方法。

**流程概览：**

1. **二级结构分析：** 遍历所有 PDB，计算其 α-螺旋 和 β-折叠 的含量，并进行粗略分类 ("alpha", "beta", "alpha\_beta")  
2. **相似性矩阵：** 使用 TMalign 对所有蛋白质进行全对全（all-vs-all）的结构比对，生成一个TM-score相似性矩阵 。  
3. **自动聚类：** 自动搜索最佳聚类数量（K），通过最小化“孤儿簇”的比例来确定最佳K值 。  
4. **结果输出：** 生成聚类成员列表、代表性结构列表，以及一个按二级结构着色的最终树状图。

## **2\. 依赖与安装**

### **2.1 外部程序**

此脚本依赖一个关键的外部生物信息学工具：

* **TMalign**：**（必需）** 用于计算结构相似性。你必须在 \--tmalign 参数中提供其可执行文件的完整路径。  
* **mkdssp**：**（可选）** 仅在 ss 参数被设置为 dssp 时才需要 。脚本默认的 phi\_psi 模式 无此依赖。

### **2.2 Python 环境**

推荐使用 conda 创建一个独立的环境（为了实验方便，我建立了环境来去做这些事情，但是实际上只需要有了这些包，就可以直接用了，搭建环境这一步可以直接跳过）

\# 1\. 创建并激活环境  
conda create \-n protein\_cluster python=3.10  
conda activate protein\_cluster

\# 2\. 安装核心包  
conda install \-c anaconda pandas numpy scipy scikit-learn  
conda install \-c conda-forge matplotlib  
conda install \-c bioconda biopython

## **3\. 命令行使用**

### **3.1 基本命令**

\# 1\. (可选) 创建输出目录  
mkdir \-p ./results

\# 2\. 运行脚本  
python classify\_and\_cluster\_proteins.py \\  
    \--inp "./my\_pdbs" \\  
    \--out "./results" \\  
    \--tmalign "/path/to/your/TMalign" \\  
    \--jobs 16

比如下面是我曾经使用过的示例：  
python \~/proteincluster/classify\_and\_cluster\_proteins.py \\  
\--inp "./final\_sample\_for\_dendrogram" \\  
\--out "./output\_final\_super\_dendrogram" \\  
\--tmalign "$(which TMalign)" \\  
\--jobs 8

需要的只是修改input，output的路径，如果tmalign没有别的需求，可以和我的示例一致。当然，在代码的378-391行有更多的可以设置的值，既可以直接在代码中修改”default”对应的值，也可以在命令行中进行设置，格式同上。

### **3.2 性能警告：O(N²) 复杂度**

**注意：** 本脚本的核心是计算一个全对全（all-vs-all）的 TMalign 矩阵 。其计算复杂度为 O(N²)，N 是输入 PDB 的数量。

* **N \= 100**：\~4950 次比较 (几分钟)  
* **N \= 1000**：\~500,000 次比较 (几小时)  
* **N \= 10,000**：\~50,000,000 次比较 (几天到几周)

**请勿**在未预处理（例如使用 Foldseek）的情况下，将上万个 PDB 文件直接作为输入。对于大规模数据集，应先进行快速聚类（如 Foldseek）提取代表，再将**代表**（N \< 1500）作为本脚本的输入，用于精细聚类和可视化。

### **3.3 foldseek加速流程**

使用TMalign方法耗时很长，生成文件很多，建议是先使用foldseek对很多蛋白进行快速筛选，使用foldseek的easycluster模式，

可以用这个先确定用几个核，笔记本电脑可以设置这个一面占用太多：THREADS=$(($(nproc)-2))

接下来就是运行这个命令：

foldseek easy-cluster input pdbs sample/my \_clusters tmp fs \\

\--tmscore-threshold 0.6 \\

\--alignment-type 1 \\

\-s 9 \\

\--threads $THREADS \\

\--gpu 1

其中，这两个参数代表了筛选的敏感度核使用的筛选的模型：

\--alignment-type 1 \\这个类似于使用的模型的类别，2会更快，效果差距不太大

\-s 9 \\这个数字越高越敏感，速度更慢

\--gpu 1经本人体感，似乎没有太大改进

更多的相关内容可以查询foldseek官方GitHub文档进行学习

也可以在：https://zread.ai/steineggerlab/foldseek这个网页上看看

## **4\. 参数说明**

#### **核心参数**

* \--inp：** **存放所有 .pdb 文件的输入目录。  
* \--out：** **存放所有输出结果（CSV, JSON, PNG文件）的目录。  
* \--tmalign：TMalign 可执行文件的**完整路径**。  
* \--jobs：  
  * ** **用于并行运行 TMalign 和二级结构分析的 CPU 核心数。  
  * **默认值：** 4 。

#### **二级结构 (SS) 参数**

* \--ss：  
  * **描述：** 计算二级结构的方法。  
  * phi\_psi：（默认）使用 biopython 内置的 Phi/Psi 角度计算，**速度快**且**无需 mkdssp**。  
  * dssp：使用 mkdssp 程序，结果更准，但**速度慢**且有**额外依赖**。  
* \--alpha\_thr, \--beta\_thr：  
  * **描述：** 判定为 "all-alpha" 或 "all-beta" 拓扑的含量百分比阈值。  
  * **默认值：** 0.8 (即 80%) 。

#### **聚类调优参数**

* \--k\_min, \--k\_max：  
  * **描述：** 自动搜索“最佳K值”时的最小和最大范围。  
  * **默认值：** k\_min=2, k\_max=16 。  
* \--mean\_thr, \--std\_ceil：  
  * **描述：** 用于筛选“最佳K值”的指标。用于定义一个“好”的簇（簇内平均TM-score \> mean\_thr，且簇内TM-score标准差 \< std\_ceil）。  
  * **默认值：** 0.7, 0.1 。

## **5\. 输出文件解析**

脚本运行成功后，你会在 \--out 目录中找到以下文件：

* dendrogram.jpg：  
  **核心可视化结果**。一个层次聚类树状图，叶子标签已根据二级结构分类自动着色。但是这个生成的图像并没有包含蛋白质的结构图，这与复现的”Hallucination of closed repeat proteins containing central pockets“这篇文章有所不同，可以手动删除图片中文字并附上蛋白质图片。此外，当用于分类蛋白过多是，图片也会较大，选用20个蛋白出图比较合理，较多的蛋白会导致图片非常长。  
* secondary\_structure\_summary.csv：  
  每个 PDB 的二级结构含量（helix\_frac, beta\_frac）及其粗略分类（topology）。  
* cluster\_selection.json：  
  记录了自动K值搜索的统计数据，以及最终选择的“最佳K值”（chosen\_k）。  
* clusters.csv：  
  最终的聚类分配表。列出了每个蛋白质及其被分配到的 cluster\_id 。  
* representatives.csv：  
  **（重要）** 为每一个 cluster\_id 列出了一个“代表性”蛋白质（即该簇的中心点，medoid）可以先对适量蛋白做一轮，然后根据这个文件中的蛋白，再做一轮，可以选出数量较为合理且图片也合理的数量的蛋白。

## **6\. 如何修改和扩展功能**

### **6.1 改进二级结构分类 **

当前的分类（ss\_summary 函数 ）只关心二级结构**含量**（例如 helix\_frac \> 0.8 ），不关心**顺序**。要实现 "alpha-beta-alpha" 这种基于拓扑的分类，需要：

1. **定位函数：** 打开 classify\_and\_cluster\_proteins.py，找到 \_ss\_from\_dssp 或 \_ss\_from\_phi\_psi 这两个辅助函数。  
2. **修改逻辑：**  
   * 让这两个函数不仅返回 frac（比例），还返回一个**代表拓扑顺序的字符串**。  
   * 例如，通过 DSSP 得到 CCCHHHHHHEEEECCCCHHHH，你可以编写一个压缩函数（run-length encoding）将其转换为 H5E4H4 (忽略C, coil)。  
3. **保存拓扑：** 在 ss\_summary 函数中 ，将这个新的拓扑字符串（H5E4H4）保存到 df\['topology'\] 列中。  
4. **修改绘图：** 在 plot\_dendrogram 函数中 ，修改 ss\_map 字典 ，使其能够识别你的新拓扑字符串并映射到特定颜色（例如，所有 H 开头的都为红色，E 开头的为蓝色等）。

### **6.2 替换/修改相似性指标**

* **位置：** tm\_matrix 和 tmalign\_pair 函数 。  
* **思路：** 脚本目前硬编码了 TMalign 。你可以修改 tmalign\_pair 函数，使其支持调用 US-align 或其他比对工具，并解析它们的输出。

### **6.3 替换/修改聚类算法**

* **位置：** select\_k\_from\_grid 和 run\_final\_clustering 函数 。  
* **思路：** 脚本目前硬编码了 AgglomerativeClustering 并需要一个复杂的K值搜索。你可以将其替换为 HDBSCAN（一种更现代的、不需要预设K值的聚类算法），这会使“自动K值选择”部分变得更简单或不再需要。

鉴于作者水平有限，且开发时间仓促，本代码及相关说明文档难免存在疏漏与不足之处。

幸桢炎 2025.11.13