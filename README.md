# PhaGCN3 

Update Log： (December 20,2024)

Several important updates: 

* **Memory Optimization**: In the previous version, PhaGCN2, we used a batch size of 1000 to save memory. This resulted in a lack of connectivity between batches in the network graph. In this update, we have adjusted the batch size to 100,000, which should be sufficient for the classification of most viral genome datasets. If with sufficient system memory, you could process millions of viral sequences in one batch.

* **Runtime Optimization**: Thanks to updates in [Python 3.13](https://www.python.org/downloads/release/python-3130/?featured_on=pythonbytes), most steps now utilize multithreading. We have implemented this in PhaGCN3, significantly increasing the speed of virus classification. On a machine with 64 CPUs and 256GB of RAM([AWS](https://aws.amazon.com/) g4ad.16xlarge instance), classifying 60,000 viral sequences takes about 7.5 hours. Because classification requires constructing a network graph, we do not recommend processing an excessively large number of sequences (>100000) in one batch. As the number of virus increases, the time and memory taken will increases exponentially
* **Precision Optimization**: During previous testing, we identified isolated connected subgraphs in the network graph predicted as "_like".  The precision of these subgraphs was poor due to uncertainties inherent in the GCN graph.  Therefore, we extracted these clusters and introduced [Genomad](https://github.com/apcamargo/genomad/) for classification.  PhaGCN3 currently only assigns a cluster ID to these nodes，but the specific classification of this cluster is now provided by [Genomad](https://github.com/apcamargo/genomad/). Viruses with the same cluster ID exhibit high similarity.
* **Results Optimization**: We have introduced a confidence score for the classification results.  Results with a confidence score above 0.5 are considered high-confidence predictions.
* **Visualization Optimization**:We now support direct output of network graph visualizations, as shown in the image below. We have tested and confirmed that the current version supports visualizing network graphs with fewer than 70,000 nodes. If you need more flexible visualization options, we also provide a network source file compatible with [Gephi](https://gephi.org/) , located at **tmp/node.csv,tmp/edge.csv** in the **results** folder.

<img src="https://wenguang.oss-cn-hangzhou.aliyuncs.com/figure/36.pdf" alt="Optimized Network Graph" style="zoom: 50%;" />

PhaGCN3 is a GCN-based model that uses deep learning classifiers to learn species mask features for novel virus taxonomy classification. The following guide will help you quickly install and run it.

## Installation

### 1. Clone the repository and create a Conda environment

Ensure that Conda is installed on your system, and follow these steps:

```bash
# Clone the PhaGCN3 repository
git clone https://github.com/xiahaolong/PhaGCN3.git
cd PhaGCN3

# Create Conda environment
conda env create -f PhaGCN3.yaml -n PhaGCN3

# Activate the Conda environment
conda activate PhaGCN3
```

### 2. Install Python 3.13t (free-thread version)

> [!TIP]
>
> Since Python 3.13t is required and should be runnable in the environment, we recommend the following installation steps, but users can also install Python 3.13t manually within the PhaGCN3 environment.



#### **2.1 Prepare the installation environment**

Run the following commands to install necessary dependencies:

```bash
sudo apt-get upgrade
sudo apt-get install build-essential libssl-dev zlib1g-dev \
    libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm \
    libncurses5-dev libncursesw5-dev xz-utils tk-dev libffi-dev liblzma-dev
```

#### **2.2 Download and extract Python 3.13**

```bash
wget https://www.python.org/ftp/python/3.13.0/Python-3.13.0.tgz
tar xzf Python-3.13.0.tgz
```

#### **2.3 Install Python 3.13 and required libraries using the script**

Run the provided installation script:

```bash
sh install_python3.13.sh
```

### 3. Install Genomad and prepare the database

#### **3.1 Install Genomad**

```bash
curl -fsSL https://pixi.sh/install.sh | bash
source ~/.bashrc
pixi global install -c conda-forge -c bioconda genomad
genomad download-database .
```

#### **3.2 Extract the database**

```bash
cd database
tar -zxvf ALL_protein.tar.gz
cd ..
```

## Example Usage

Here is a complete example of how to use PhaGCN3:

### Input Files

The repository provides an example file `contigs.fa`, which contains overlapping clusters simulated from E. coli bacteriophage.

### Run PhaGCN3

Before each run, activate the environment:

```bash
conda activate PhaGCN3
```

Run the following command:

```bash
python run_Speed_up.py --contigs contigs.fa --outpath result
```

#### Parameter Descriptions

- `--contigs`：Path to the contig file.
- `--clustering`：Batch size for prediction (default: 100,000).
- `--outpath`：Path to output the results.

### Output Files

The output files will be in the specified `outpath` directory, including:

1. `final_prediction.csv`：Final prediction results.
2. `final_network.ntw`：Network structure file.
3. `tmp/node.csv,tmp/edge.csv`: Edge and node files for direct input into [Gephi](https://gephi.org/) or[cytoscape](https://cytoscape.org/).
4. Additional files (when there is a "_like" classification error):
   - `processed_test_nodes.csv`
   - `filtered_test_edges.csv`
   - `subgraph_nodes.csv`
   - `processed_test_nodes_taxonomy.tsv`(including Genomad annotation results)

## Drawing Network Graph

> [!TIP]
>
> The draw_network.py has been updated. It uses the 16,000 sequences from ictv's vmr_msl39 as test sequences to generate the PhaGCN3 network graph. It only takes ten minutes to produce the graph, and the clustering effect is better than that of Gephi. In practice, it can successfully draw a network graph for 280,000 test sequences with an edge file size of 20 GB, which takes approximately 10 hours.



### 1. Create the drawing environment

```bash
conda env create -f draw.yaml -n draw
```

### 2. Run the drawing program

Make sure the output directory for the drawing program matches the output directory for PhaGCN3. The final graph will be in your `outpath`：

```bash
conda activate draw
python draw_network.py --outpath result
```

