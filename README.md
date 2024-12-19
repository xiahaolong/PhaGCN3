# PhaGCN3 

Update Log： (December 20,2024)

Several important updates: 

* **Memory Optimization**: In the previous version, PhaGCN2, we used a batch size of 1000 to save memory. This resulted in a lack of connectivity between batches in the network graph. In this update, we have adjusted the batch size to 100,000, which should be sufficient for the classification of most viral genome datasets. If with sufficient system memory, you could process millions of viral sequences in one batch.

* **Runtime Optimization**: Thanks to updates in [Python 3.13](https://www.python.org/downloads/release/python-3130/?featured_on=pythonbytes), most steps now utilize multithreading. We have implemented this in PhaGCN3, significantly increasing the speed of virus classification. On a machine with 64 CPUs and 256GB of RAM, classifying 60,000 viral sequences takes about 7.5 hours. Because classification requires constructing a network graph, we do not recommend processing an excessively large number of sequences (>100000) in one batch. As the number of virus increases, the time and memory taken will increases exponentially
* **Precision Optimization**: During previous testing, we identified isolated connected subgraphs in the network graph predicted as "_like".  The precision of these subgraphs was poor due to uncertainties inherent in the GCN graph.  Therefore, we extracted these clusters and introduced [Genomad](https://github.com/apcamargo/genomad/) for classification.  PhaGCN3 currently only assigns a cluster ID to these nodes，but the specific classification of this cluster is now provided by [Genomad](https://github.com/apcamargo/genomad/). Viruses with the same cluster ID exhibit high similarity.
* **Results Optimization**: We have introduced a confidence score for the classification results.  Results with a confidence score above 0.5 are considered high-confidence predictions.
* **Visualization Optimization**:We now support direct output of network graph visualizations, as shown in the image below. We have tested and confirmed that the current version supports visualizing network graphs with fewer than 70,000 nodes. If you need more flexible visualization options, we also provide a network source file compatible with [Gephi](https://gephi.org/) , located at **tmp/node.csv,tmp/edge.csv** in the **results** folder.

<img src="https://wenguang.oss-cn-hangzhou.aliyuncs.com/figure/image-20241218170530956.png" alt="image-20241218170530956" style="zoom: 50%;" />

PhaGCN3 是一个基于 GCN 的模型，它可以通过深度学习分类器学习物种掩码特征，用于新的病毒分类学分类。以下指南将帮助您快速完成安装和运行。

## 安装

### 1. 克隆存储库并创建 Conda 环境

确保您的系统已安装 Conda，并按照以下步骤操作：

```bash
# 克隆 VirusGCN 存储库
git clone https://github.com/xiahaolong/VirusGCN.git
cd VirusGCN

# 创建 Conda 环境
conda env create -f VirusGCN.yaml -n virusgcn

#进入conda环境
conda activate virusgcn
```

### 2. 安装 Python 3.13t（自由线程版本）

> [!TIP]
>
> 由于需要安装python3.13t版本，并保证该版本在环境中可运行，我们推荐以下安装步骤，用户也可以自行在virusgcn环境中安装python3.13t



#### **2.1 准备安装环境**

执行以下命令安装必要的依赖：

```bash
sudo apt-get upgrade
sudo apt-get install build-essential libssl-dev zlib1g-dev \
    libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm \
    libncurses5-dev libncursesw5-dev xz-utils tk-dev libffi-dev liblzma-dev
```

#### **2.2 下载并解压 Python 3.13**

```bash
wget https://www.python.org/ftp/python/3.13.0/Python-3.13.0.tgz
tar xzf Python-3.13.0.tgz
```

#### **2.3 使用脚本安装 Python 3.13以及需要的库**

运行提供的安装脚本：

```bash
sh install_python3.13.sh
```

### 3. 安装 Genomad 并准备数据库

#### **3.1 安装 Genomad**

```bash
curl -fsSL https://pixi.sh/install.sh | bash
source ~/.bashrc
pixi global install -c conda-forge -c bioconda genomad
genomad download-database .
```

#### **3.2 解压数据库**

```bash
cd database
tar -zxvf ALL_protein.tar.gz
cd ..
```

## 使用示例

以下是使用 VirusGCN 的完整示例：

### 输入文件

存储库中提供了 `contigs.fa` 示例文件，其中包含从大肠杆菌噬菌体模拟的重叠群。

### 运行 VirusGCN

在每次运行 VirusGCN 前，激活环境：

```bash
conda activate virusgcn
```

运行以下命令：

```bash
python run_Speed_up.py --contigs contigs.fa --outpath result
```

#### 参数说明

- `--contigs`：重叠群文件路径。
- `--clustering`：预测批次大小（默认值：100,000）。
- `--outpath`：输出结果路径。

### 输出文件

输出文件位于指定的 `outpath` 目录中，包括：

1. `final_prediction.csv`：最终预测结果。
2. `final_network.ntw`：网络结构文件。
3. `tmp/node.csv,tmp/edge.csv`: 可用于直接输入 [Gephi](https://gephi.org/) 或者[cytoscape](https://cytoscape.org/)的边文件及点文件
4. 额外文件（当出现“_like”分类误差时）：
   - `processed_test_nodes.csv`
   - `filtered_test_edges.csv`
   - `processed_test_nodes_taxonomy.tsv`（含 Genomad 注释结果）

## 绘制网络图

> [!TIP]
>
> 类似于gephi，当点数量多的时候，该步骤可能需要大量计算网络图拓扑或者绘制时间，经测试，60000个点大约需要3小时运行绘制时间，因此我们将该步骤设置为单独可选式运行



### 1. 创建绘图环境

```bash
conda env create -f draw.yaml -n draw
```

### 2. 运行绘图程序

确保绘图程序的输出目录与 VirusGCN 的输出目录一致,最终的图在您的outpath中：

```bash
conda activate draw
python draw_network.py --outpath result
```

