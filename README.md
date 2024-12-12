# VirusGCN

VirusGCN 是一个专为病毒基因组分析设计的工具，支持 GPU 模式运行。以下指南将帮助您快速完成安装和运行。

## 安装指南

### 1. 克隆存储库并创建 Conda 环境

确保您的系统已安装 CUDA，并按照以下步骤操作：

```bash
# 克隆 VirusGCN 存储库
git clone https://github.com/xiahaolong/VirusGCN.git
cd VirusGCN

# 创建 Conda 环境
conda env create -f VirusGCN.yaml -n virusgcn
```

### 2. 安装 Python 3.13t（自由线程版本）

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
sudo wget https://www.python.org/ftp/python/3.13.0/Python-3.13.0.tgz
sudo tar xzf Python-3.13.0.tgz
```

#### **2.3 使用脚本安装 Python 3.13**

运行提供的安装脚本：

```bash
./install_python3.13.sh
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
python run_Speed_up.py --contigs contigs.fa --len 8000 --outpath result
```

#### 参数说明

- `--contigs`：重叠群文件路径。
- `--len`：重叠群长度（默认值：8000 bp，支持的最短长度为 1700 bp）。
- `--clustering`：预测批次大小（默认值：100,000）。
- `--outpath`：输出结果路径。

### 输出文件

输出文件位于指定的 `outpath` 目录中，包括：

1. `final_prediction.csv`：最终预测结果。
2. `final_network.ntw`：网络结构文件。
3. 额外文件（当出现分类误差时）：
   - `processed_test_nodes.csv`
   - `filtered_test_edges.csv`
   - `processed_test_nodes_taxonomy.tsv`（含 Genomad 注释结果）

## 绘制网络图

### 1. 创建绘图环境

```bash
conda env create -f draw.yaml -n draw
```

### 2. 运行绘图程序

确保绘图程序的输出目录与 VirusGCN 的输出目录一致：

```bash
conda activate draw
python draw_network.py --outpath result
```

