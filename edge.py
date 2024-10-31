import os
import numpy as np
import pandas as pd
import networkx as nx
import joblib
import argparse
from concurrent.futures import ProcessPoolExecutor

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--outpath', type=str, default="result")
parser.add_argument('--n', type=int, default=0)
args = parser.parse_args()

out_f = f"{args.outpath}/out/"
contig_in = f"{args.outpath}/input/"
contig_out = f"{args.outpath}/single_contig/"
file_in_fn = f"{args.outpath}/single_contig/"
file_out_fn = f"{args.outpath}/all_proteins/"
Knowledge_graph = f"{args.outpath}/Cyber_data/"
all_protein_f = out_f + "all_translate_proteins.fa"
proteins_aa_fp = all_protein_f
diamond_out_fn = '{}.diamond.tab'.format(os.path.basename(proteins_aa_fp).rsplit('.', 1)[0])
diamond_out_fp = os.path.join(out_f, diamond_out_fn)
contig_abc_fp = out_f + diamond_out_fn + ".abc"
print("\n\n" + "{:-^80}".format("Calculating E-edges"))

# Loading database
database = "database/"
gene2genome = pd.read_csv(database + "ALL_gene_to_genomes.csv")
contig_id = gene2genome["contig_id"].values
contig_id = [item.replace(" ", "~") for item in contig_id]
gene2genome["contig_id"] = contig_id

# 创建蛋白质到参考基因组的映射字典
protein_to_ref = {protein: ref for protein, ref in zip(gene2genome["protein_id"].values, gene2genome["contig_id"].values)}

# 创建 contig 的唯一标识符集合
contig_set = list(set(gene2genome["contig_id"].values))
ID_to_ref = {i: ref for i, ref in enumerate(contig_set)}
ref_to_ID = {ref: i for i, ref in enumerate(contig_set)}

# 文件处理
fn = f"{args.outpath}/single_contig/"
contig_to_id = {}
file_list = os.listdir(fn)
file_list = sorted(file_list)
for file_n in file_list:
    name = file_n.split(".")[0]
    contig_to_id[name] = file_list.index(file_n)

# 记录每个 contig 的行 ID
id_to_contig = {value: key for key, value in contig_to_id.items()}

# 文件路径设置
fn = f"{args.outpath}/out/"

# 读取 blastp 结果文件
blastp = pd.read_csv(contig_abc_fp, sep=" ", names=["contigs", "ref", "e-value"])

# 读取基因与基因组映射的 CSV 文件
gene_to_genome = pd.read_csv(fn + "contig_gene_to_genome.csv", sep=",")

# 初始化 e_matrix
e_matrix = np.ones((len(contig_to_id), len(ref_to_ID.keys())))
joblib.dump(contig_to_id, open(f"{args.outpath}/out/contig_to_id.joblib", "wb"))
# 从 blastp 数据框中提取数据
blast_contigs = blastp["contigs"].values
blast_ref = blastp["ref"].values
blast_value = blastp["e-value"].values

# 更新 e_matrix 的函数
def update_e_matrix(i):
    contig_name = gene_to_genome[gene_to_genome["protein_id"] == blast_contigs[i]]["contig_id"].values[0]
    row_id = contig_to_id[contig_name]
    reference = protein_to_ref[blast_ref[i]]
    col_id = ref_to_ID[reference]
    e_value = float(blast_value[i])
    if e_value == 0:
        e_value = 1e-250
    return row_id, col_id, e_value

# 使用多核处理更新 e_matrix
with ProcessPoolExecutor() as executor:
    results = list(executor.map(update_e_matrix, range(len(blast_contigs))))

for row_id, col_id, e_value in results:
    if e_matrix[row_id][col_id] == 1:
        e_matrix[row_id][col_id] = e_value
    else:
        e_matrix[row_id][col_id] += e_value

# 计算 e_weight
e_weight = -np.log10(e_matrix) - 50
e_weight[e_weight < 1] = 0

# 生成 P-edges
print("\n\n" + "{:-^80}".format("Calculating P-edges"))

# 创建从参考名称到 ID 的映射
name_to_id = {}
reference_df = pd.read_csv("database/reference_name_id.csv")
tmp_ref = reference_df["name"].values
tmp_id = reference_df["idx"].values
for ref, idx in zip(tmp_ref, tmp_id):
    name_to_id[ref.replace(" ", "~")] = idx

# 读取网络边数据
edges = pd.read_csv(out_f + "network.ntw", sep=' ', names=["node1", "node2", "weight"])

# 读取合并的基因组文件和分类标签文件
merged_df = pd.read_csv(database + "ALL_genome_profile.csv", header=0, index_col=0)
Taxonomic_df = pd.read_csv(database + "taxonomic_label.csv")

# 合并数据框
merged_df = pd.merge(merged_df, Taxonomic_df, left_on="contig_id", right_on="contig_id", how="inner")

# 提取 contig_id 和 family 信息
contig_id = merged_df["contig_id"].values
family = merged_df["class"].values

# 创建 contig 到 family 的映射
contig_to_family = {name: fam for name, fam in zip(contig_id, family) if type(fam) != type(np.nan)}
joblib.dump(contig_to_family, open(f"{args.outpath}/out/contig_to_family.joblib", "wb"))
# 创建一个空的 NetworkX 图对象
G = nx.Graph()

# 将 P-edges 添加到图中
with open(out_f + "/network.ntw") as file_in:
    for line in file_in.readlines():
        tmp = line[:-1].split(" ")
        node1 = tmp[0]
        node2 = tmp[1]
        weight = float(tmp[2])

        if "~" in node1 and node1 not in name_to_id.keys():
            print(node1)
            print("ERROR")
            exit(1)
        if "~" in node2 and node2 not in name_to_id.keys():
            print(node2)
            print("ERROR")
            exit(1)

        G.add_edge(node1, node2, weight=1)

# 添加 E-edges 的函数
def add_edges_to_graph(row_id):
    contig_name = id_to_contig[row_id]
    cnt = 0
    sorted_idx = np.argsort(e_weight[row_id])
    for j in range(5):
        idx = sorted_idx[-j]
        if e_weight[row_id][idx] != 0:
            ref_name = ID_to_ref[idx]
            if ref_name in G.nodes():
                G.add_edge(contig_name, ref_name, weight=1)
                cnt += 1
    return cnt

# 使用多核处理添加 E-edges
with ProcessPoolExecutor() as executor:
    counts = list(executor.map(add_edges_to_graph, range(e_weight.shape[0])))

# 统计总数
total_count = sum(counts)

# 移除未压缩的节点
node_list = list(G.nodes())
for node in node_list:
    if "~" in node and node not in contig_to_family.keys():
        G.remove_node(node)
joblib.dump(G, open(f"{args.outpath}/out/graph.joblib", "wb"))
