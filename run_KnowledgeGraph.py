import os
import sys
import Bio
import logging
import argparse
import subprocess
import scipy as sp
import numpy as np
import pandas as pd
import pickle as pkl
import networkx as nx
import scipy.stats as stats
import scipy.sparse as sparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import joblib
from scipy.sparse import csr_matrix, save_npz

# Defined folder
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--outpath', type=str, default = "result")
parser.add_argument('--n', type=int, default = 0)
args = parser.parse_args()

out_f = f"{args.outpath}/out/"
contig_in = f"{args.outpath}/input/"
contig_out = f"{args.outpath}/single_contig/"
file_in_fn = f"{args.outpath}/single_contig/"
file_out_fn = f"{args.outpath}/all_proteins/"
Knowledge_graph = f"{args.outpath}/Cyber_data/"
################################################################################
############################  Check the folder #################################
################################################################################
if not os.path.exists(out_f):
    _ = os.makedirs(out_f)
else:
    print("folder {0} exist... cleaning dictionary".format(out_f))
    if os.listdir(out_f):
        try:
            _ = subprocess.check_call("rm -rf {0}".format(out_f), shell=True)
            _ = os.makedirs(out_f)
            print("Dictionary cleaned")
        except:
            print("Cannot clean your folder... permission denied")
            exit(1)

if not os.path.exists(file_in_fn):
    _ = os.makedirs(file_in_fn)
else:
    print("folder {0} exist... cleaning dictionary".format(file_in_fn))
    if os.listdir(file_in_fn):
        try:
            _ = subprocess.check_call("rm -rf {0}".format(file_in_fn), shell=True)
            _ = os.makedirs(file_in_fn)
            print("Dictionary cleaned")
        except:
            print("Cannot clean your folder... permission denied")
            exit(1)

if not os.path.exists(file_out_fn):
    _ = os.makedirs(file_out_fn)
else:
    print("folder {0} exist... cleaning dictionary".format(file_out_fn))
    if os.listdir(file_out_fn):
        try:
            _ = subprocess.check_call("rm -rf {0}".format(file_out_fn), shell=True)
            _ = os.makedirs(file_out_fn)
            print("Dictionary cleaned")
        except:
            print("Cannot clean your folder... permission denied")
            exit(1)



################################################################################
############################  Rename the files  ################################
################################################################################

# Split contigs files into single contig per file.
file_list = sorted(os.listdir(contig_in))
seq = []
old_file_id = 0
contig_id = 0
with open(f"{args.outpath}/name_list.csv",'w') as list_out:
    list_out.write("contig_name,idx\n")
    for file_n in file_list:
        for record in SeqIO.parse(contig_in+file_n, "fasta"):
            name = str(old_file_id) + "_" + str(contig_id)
            contig_id += 1
            list_out.write(record.id + "," + name + "\n")
            _ = SeqIO.write(record, contig_out+name+".fasta", "fasta")
        old_file_id += 1


################################################################################
###################### Translate contigs into 6 ORFs ###########################
################################################################################

def return_protien(contig):
    frame = Bio.Seq.translate(contig)
    protein_list = frame.split("*")
    sorted_list = sorted(protein_list, key=len, reverse=True)
    sorted_list = [item for item in sorted_list if len(item) > 10 ]
    return sorted_list


file_list = os.listdir(file_in_fn)
for file in file_list:
    old_file_name = file.rsplit(".", 1)[0]
    contig_id = int(old_file_name.split("_")[-1])
    label_id = int(old_file_name.split("_")[0])
    for record in SeqIO.parse(file_in_fn+file, "fasta"):
        protein_record = []
        contig = str(record.seq)
        # conver into protein
        frame1 = return_protien(contig)
        frame2 = return_protien(contig[1:])
        frame3 = return_protien(contig[2:])
        rev_contig = Bio.Seq.reverse_complement(contig)
        frame4 = return_protien(rev_contig)
        frame5 = return_protien(rev_contig[1:])
        frame6 = return_protien(rev_contig[2:])
        proteins = np.concatenate([frame1, frame2, frame3, frame4, frame5, frame6])
        for i in range(len(proteins)):
            rec = SeqRecord(Seq(proteins[i]), id=str(label_id)+ "_" + str(contig_id) + "_" + str(i), description="")
            protein_record.append(rec)
        _ = SeqIO.write(protein_record, file_out_fn+file, "fasta")

all_protein_f = out_f + "all_translate_proteins.fa"
subprocess.run(f"find {file_out_fn} -type f -exec cat {{}} + > {all_protein_f}", shell=True, check=True)


################################################################################
############################## Run diamond BLASTp  #############################
################################################################################

print("\n\n" + "{:-^80}".format("Diamond BLASTp"))
print("Creating Diamond database and running Diamond...")

def make_diamond_db(cpu: int):
    diamond_db_bp = "database/database.dmnd"
    aa_fp = "database/ALL_protein.fasta"

    make_diamond_cmd = ['diamond', 'makedb', '--threads', str(cpu), '--in', aa_fp, '-d', diamond_db_bp]
    print("Creating Diamond database...")
    res = subprocess.run(make_diamond_cmd, check=True, stdout=subprocess.PIPE)
    if res.returncode != 0:
        print('Error creating Diamond database')
        exit(1)
    diamond_db_fp = diamond_db_bp + '.dmnd'
    return diamond_db_fp


def run_diamond(aa_fp, db_fp, cpu: int, diamond_out_fn):
    # More sensitive as an option?
    diamond_cmd = ['diamond', 'blastp', '--threads', str(cpu), '--sensitive', '-d', db_fp, '-q', aa_fp,
                   '-o', diamond_out_fn]
    print("Running Diamond...")
    res = subprocess.run(diamond_cmd, check=True, stdout=subprocess.PIPE)
    if res.returncode != 0:
        print('Error running Diamond')
        exit(1)
    return diamond_out_fn

proteins_aa_fp = all_protein_f
db_fp = "database/database.dmnd"
diamond_out_fn = '{}.diamond.tab'.format(os.path.basename(proteins_aa_fp).rsplit('.', 1)[0])
diamond_out_fp = os.path.join(out_f, diamond_out_fn)
# Create database
_ = make_diamond_db(128)
# Run BLASTP
similarity_fp = run_diamond(proteins_aa_fp, db_fp, 128, diamond_out_fp)

# capture the query, referencde, e-value from diamond output
contig_abc_fp = out_f + diamond_out_fn + ".abc"
abc_fp = out_f+"merged.abc"
_ = subprocess.check_call("awk '$1!=$2 {{print $1,$2,$11}}' {0} > {1}".format(diamond_out_fp, contig_abc_fp), shell=True)
_ = subprocess.check_call("cat database/database.self-diamond.tab.abc {0} > {1}".format(contig_abc_fp, abc_fp), shell=True)

# Generating gene-to-genome.csv: protein_id, contig_id, keywords
blastp = pd.read_csv(contig_abc_fp, sep=' ', names=["contig", "ref", "e-value"])
protein_id = sorted(list(set(blastp["contig"].values)))
contig_id = [item.rsplit("_", 1)[0] for item in protein_id]
description = ["hypothetical protein" for item in protein_id]
gene2genome = pd.DataFrame({"protein_id": protein_id, "contig_id": contig_id ,"keywords": description})
gene2genome.to_csv(out_f+"contig_gene_to_genome.csv", index=None)



# Combining the gene-to-genomes files
_ = subprocess.check_call("cat database/ALL_gene_to_genomes.csv {0}contig_gene_to_genome.csv > {1}gene_to_genome.csv".format(out_f, out_f), shell=True)



# Run MCL
print("\n\n" + "{:-^80}".format("Protein clustering"))
print("Loading proteins...")
gene2genome_fp = out_f+"gene_to_genome.csv"
gene2genome_df = pd.read_csv(gene2genome_fp, sep=',', header=0)
def make_protein_clusters_mcl(abc_fp, out_p, inflation=2):   
    
    print("Running MCL...")  
    abc_fn = "merged"  
    mci_fn = '{}.mci'.format(abc_fn)  
    mci_fp = os.path.join(out_p, mci_fn)  
    mcxload_fn = '{}_mcxload.tab'.format(abc_fn)  
    mcxload_fp = os.path.join(out_p, mcxload_fn)  

    
    subprocess.check_call("mcxload -abc {0} --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o {1} -write-tab {2}".format(
        abc_fp, mci_fp, mcxload_fp), shell=True)  

    
    mcx_fp = '{}.mcx'.format(abc_fn)  
    mcx_fp_full = os.path.join(out_p, mcx_fp)  
    subprocess.check_call("mcx convert {0} {1}".format(mci_fp, mcx_fp_full), shell=True)
    mcl_clstr_fn = "{0}_mcl{1}.clusters".format(abc_fn, int(inflation*10))  
    mcl_clstr_fp = os.path.join(out_p, mcl_clstr_fn)  

    subprocess.check_call("mcl {0} -I {1} -use-tab {2} -o {3} -te 128".format(
        mcx_fp_full, inflation, mcxload_fp, mcl_clstr_fp), shell=True)  
    
    return mcl_clstr_fp

pc_overlap, pc_penalty, pc_haircut, pc_inflation = 0.8, 2.0, 0.1, 2.0
pcs_fp = make_protein_clusters_mcl(abc_fp, out_f, pc_inflation)



def load_mcl_clusters(fi):
    """
    Load given clusters file
    
    Args:
        fi (str): path to clusters file
        proteins_df (dataframe): A dataframe giving the protein and its contig.
    Returns: 
        tuple: dataframe proteins and dataframe clusters
    """
    # Read MCL
    with open(fi) as f:
        c = [line.rstrip("\n").split("\t") for line in f]
    c = [x for x in c if len(c) > 1]
    nb_clusters = len(c)
    formatter = "PC_{{:>0{}}}".format(int(round(np.log10(nb_clusters))+1))
    name = [formatter.format(str(i)) for i in range(nb_clusters)]
    size = [len(i) for i in c]
    clusters_df = pd.DataFrame({"size": size, "pc_id": name}).set_index("pc_id")
    return clusters_df, name, c


def build_clusters(fp, gene2genome):
    """
        Build clusters given clusters file

        Args:
            fp (str): filepath of clusters file
            gene2genome (dataframe): A dataframe giving the protein and its genome.
            mode (str): clustering method
        Returns:
            tuple: dataframe of proteins, clusters, profiles and contigs
        """
    # Read MCL
    clusters_df, name, c = load_mcl_clusters(fp)
    print("Using MCL to generate PCs.")
    # Assign each prot to its cluster
    gene2genome.set_index("protein_id", inplace=True)  # id, contig, keywords, cluster
    for prots, clust in zip(c, name):
        try:
            gene2genome.loc[prots, "cluster"] = clust
        except KeyError:
            prots_in = [p for p in prots if p in gene2genome.index]
            not_in = frozenset(prots) - frozenset(prots_in)
            print("{} protein(s) without contig: {}".format(len(not_in), not_in))
            gene2genome.loc[prots_in, "cluster"] = clust
    # Keys
    for clust, prots in gene2genome.groupby("cluster"):
        clusters_df.loc[clust, "annotated"] = prots.keywords.count()
        if prots.keywords.count():
            keys = ";".join(prots.keywords.dropna().values).split(";")
            key_count = {}
            for k in keys:
                k = k.strip()
                try:
                    key_count[k] += 1
                except KeyError:
                    key_count[k] = 1
            clusters_df.loc[clust, "keys"] = "; ".join(["{} ({})".format(x, y) for x, y in key_count.items()])
    gene2genome.reset_index(inplace=True)
    clusters_df.reset_index(inplace=True)
    profiles_df = gene2genome.loc[:, ["contig_id", "cluster"]].drop_duplicates()
    profiles_df.columns = ["contig_id", "pc_id"]
    contigs_df = pd.DataFrame(gene2genome.fillna(0).groupby("contig_id").count().protein_id)
    contigs_df.index.name = "contig_id"
    contigs_df.columns = ["proteins"]
    contigs_df.reset_index(inplace=True)
    return gene2genome, clusters_df, profiles_df, contigs_df

print("Building the cluster and profiles (this may take some time...)")

protein_df, clusters_df, profiles_df, contigs_df = build_clusters(pcs_fp, gene2genome_df)

print("Saving files")
dfs = [gene2genome_df, contigs_df, clusters_df]
names = ['proteins', 'contigs', 'pcs']
output_dir = out_f



for name, df in zip(names, dfs):
    fn = "Cyber_{}.csv".format(name)
    fp = os.path.join(output_dir, fn)
    index_id = name.strip('s') + '_id'
    if not os.path.exists(fp):
        df.set_index(index_id).to_csv(fp)
    else:
        print("File {} exists and will be used. Use -f to overwrite.".format(fn))

        
        
        
profiles_fn = "Cyber_profiles.csv"
profiles_fp = os.path.join(out_f, profiles_fn)
if not os.path.exists(profiles_fp):
    profiles_df.to_csv(profiles_fp, index=False)
else:
    print("File {} exists and will be used. Use -f to overwrite.".format(profiles_fn))



# Create P-edges
def build_pc_matrices(profiles, contigs, pcs):
    """
    Build the pc profiles matrices (shared & singletons) from dataframes.
    
    Args:
        profiles (dataframe): required fields are contig_id and pc_id.
        contigs (dataframe): contigs info, required fields are proteins, pos and id.
        pcs (dataframe): pcs info, required fields are pos and id.

    Returns:
        (tuple of np.ndarray): Shared PCs and singletons matrix as NumPy arrays.
    """
    pc_by_cont = profiles.groupby("contig_id").count().pc_id
    pc_by_cont = pd.merge(contigs.sort_values("pos").loc[:, ["pos", "contig_id", "proteins"]], pc_by_cont.to_frame(), how="left",
                          left_on="contig_id", right_on="contig_id").fillna(0)
    singletons = (pc_by_cont.proteins - pc_by_cont.pc_id).values

    # Matrix
    profiles.index.name = "pos"
    profiles.reset_index(inplace=True)
    profiles = pd.merge(profiles, pcs.loc[:, ["pc_id", "pos"]], left_on="pc_id", right_on="pc_id", how="inner",
                            suffixes=["", "_pc"])
    profiles = pd.merge(profiles, contigs.loc[:, ["contig_id", "pos"]], left_on="contig_id", right_on="contig_id", how="inner",
                            suffixes=["", "_contig"])
    profiles = profiles.loc[:, ["pos_contig", "pos_pc"]]
    row_indices, col_indices = zip(*profiles.values)  # 将结果解压为行和列的索引

    
    matrix = np.zeros((len(contigs), len(pcs)))
    for r, c in zip(row_indices, col_indices):
        matrix[r, c] = 1  

    return matrix, singletons  


def run_network_computation():
    cmd = f"python3.13t network_compute.py --outpath {args.outpath}"
    out = subprocess.check_call(cmd, shell=True)

    
    S = np.load(f"{args.outpath}/out/output_network.npz")
    return S


# Loding dataset
contigs_df = pd.read_csv(f"{args.outpath}/out/Cyber_contigs.csv")
clusters_df = pd.read_csv(f"{args.outpath}/out/Cyber_pcs.csv")
profiles_df = pd.read_csv(f"{args.outpath}/out/Cyber_profiles.csv")

# Replace names
contigs_csv_df = contigs_df.copy()
contigs_csv_df['contig_id'] = contigs_csv_df['contig_id'].str.replace(' ', '~')
print("Read {} entries from {}".format(len(contigs_csv_df), os.path.join(output_dir, '{}_contigs.csv'.format(name))))
contigs_csv_df.index.name = "pos"
contigs_csv_df.reset_index(inplace=True)

pcs_csv_df = clusters_df.copy()
profiles = profiles_df.copy()
profiles['contig_id'] = profiles['contig_id'].str.replace(' ', '~')  # ClusterONE can't handle spaces

# Filtering the PC profiles that appears only onc
before_filter = len(profiles)
cont_by_pc = profiles.groupby("pc_id").count().contig_id.reset_index()

# get the number of contigs for each pcs and add it to the dataframe
cont_by_pc.columns = ["pc_id", "nb_proteins"]
pcs_csv_df = pd.merge(pcs_csv_df, cont_by_pc, left_on="pc_id", right_on="pc_id", how="left")
pcs_csv_df.fillna({"nb_proteins": 0}, inplace=True)

# Drop the pcs that <= 1 contig from the profiles.
pcs_csv_df = pcs_csv_df[pcs_csv_df['nb_proteins'] > 1]  # .query("nb_contigs>1")
at_least_a_cont = cont_by_pc[cont_by_pc['nb_proteins'] > 1]  # cont_by_pc.query("nb_contigs>1")
profiles = profiles[profiles['pc_id'].isin(at_least_a_cont.pc_id)]
print("Read {} entries (dropped {} singletons) from {}".format(len(profiles), (before_filter - len(profiles)), profiles_fp))
pcs_csv_df = pcs_csv_df.reset_index(drop=True)
pcs_csv_df.index.name = "pos"
pcs_csv_df = pcs_csv_df.reset_index()

matrix, singletons = build_pc_matrices(profiles, contigs_csv_df, pcs_csv_df)
profiles_csv = {"matrix": matrix, "singletons": singletons}
merged_df = contigs_csv_df
merged_fp = os.path.join(output_dir, 'merged_df.csv')
merged_df.to_csv(merged_fp)

with open(f"{args.outpath}/out/matrix.npy", 'wb') as f:
    np.save(f, matrix)
with open(f"{args.outpath}/out/singletons.npy", 'wb') as f:
    np.save(f, singletons)




run_network_computation()
data = np.load(f"{args.outpath}/out/output_network.npz")
ntw = data['S']
#sparse_ntw = csr_matrix(ntw)
#save_npz("output_network_sparse.npz", sparse_ntw)

def to_clusterer(matrix, fi, contigs=None,names=None):
    """Save a network in a file ready for MCL and/or ClusterONE

    Args:
        matrix (scipy.sparse_matrix): network.
        fi (str): filename .
        names (pandas.dataframe): with the columns
            "pos":  (int) is the position in the matrix.
            "id": (str) column contain the id of the node.
            If None, self.contigs is used.

    Returns:
        str: filename
    """
    names = contigs if names is None else names
    names = names.set_index("pos").contig_id
    with open(fi, "wt") as f:
        matrix = sparse.dok_matrix(matrix)
        for r, c in zip(*matrix.nonzero()):
            f.write(" ".join([str(x) for x in (names[r], names[c], matrix[r, c])]))
            f.write("\n")
    print("Saving network in file {0} ({1} lines).".format(fi, matrix.getnnz()))
    return fi

fi = to_clusterer(ntw, out_f+"network.ntw", merged_df.copy())

print(os.getcwd())
cmd = f"python3.13t edge.py --outpath {args.outpath}"
out = subprocess.check_call(cmd, shell=True)
G = joblib.load(open(f"{args.outpath}/out/graph.joblib", "rb"))


name_to_id = {}
reference_df = pd.read_csv("database/reference_name_id.csv")
tmp_ref = reference_df["name"].values
tmp_id  = reference_df["idx"].values
for ref, idx in zip(tmp_ref,tmp_id):
    name_to_id[ref.replace(" ", "~")] = idx
    
test_to_id = {}
class_to_label = {56:56,62:62,174:174,199:199,66:66,117:117,29:29,157:157,212:212,38:38,9:9,106:106,153:153,46:46,98:98,154:154,79:79,129:129,101:101,108:108,204:204,203:203,124:124,191:191,188:188,1:1,10:10,180:180,131:131,42:42,58:58,159:159,173:173,175:175,20:20,166:166,120:120,88:88,151:151,134:134,67:67,27:27,2:2,90:90,0:0,182:182,208:208,183:183,13:13,57:57,21:21,47:47,95:95,186:186,43:43,91:91,181:181,63:63,12:12,8:8,197:197,209:209,195:195,99:99,135:135,143:143,3:3,64:64,69:69,105:105,19:19,111:111,142:142,190:190,150:150,50:50,68:68,28:28,94:94,17:17,44:44,130:130,24:24,81:81,107:107,78:78,83:83,5:5,41:41,176:176,104:104,102:102,72:72,210:210,140:140,74:74,185:185,86:86,163:163,200:200,70:70,97:97,155:155,7:7,138:138,198:198,145:145,128:128,127:127,39:39,16:16,217:217,184:184,15:15,169:169,126:126,100:100,139:139,61:61,14:14,137:137,11:11,164:164,122:122,193:193,60:60,85:85,213:213,149:149,35:35,172:172,82:82,156:156,45:45,114:114,77:77,206:206,215:215,152:152,161:161,30:30,113:113,84:84,53:53,89:89,168:168,48:48,144:144,211:211,148:148,54:54,109:109,37:37,125:125,189:189,205:205,87:87,165:165,93:93,52:52,22:22,115:115,192:192,49:49,202:202,71:71,214:214,103:103,4:4,160:160,34:34,55:55,179:179,110:110,133:133,201:201,158:158,170:170,123:123,96:96,132:132,40:40,59:59,75:75,162:162,65:65,33:33,136:136,36:36,141:141,32:32,76:76,167:167,92:92,171:171,6:6,23:23,146:146,73:73,194:194,216:216,177:177,207:207,51:51,80:80,196:196,31:31,178:178,116:116,147:147,119:119,118:118,187:187,112:112,121:121,18:18,25:25,26:26}
# Generating the Knowledge Graph
print("\n\n" + "{:-^80}".format("Generating Knowledge graph"))
contig_to_family = joblib.load(open(f"{args.outpath}/out/contig_to_family.joblib", "rb"))
mode = "testing"
if mode == "validation":
    test_mask = []
    label = []
    cnt = 0
    for node in G.nodes():
        try:
            label.append(class_to_label[contig_to_family[node]])
            cnt+=1
        except:
            if "_" in node:
                try:
                    class_ = int(node.split("_")[0])
                    label.append(class_)
                    test_mask.append(cnt)
                    test_to_id[node] = cnt
                    cnt+=1
                except:
                    print(node)
            else:
                print(node)
    pkl.dump(test_mask, open(f"{args.outpath}/Cyber_data/contig.mask", "wb" ) )
    pkl.dump(label, open(f"{args.outpath}/Cyber_data/contig.label", "wb" ) )
    adj = nx.adjacency_matrix(G)
    pkl.dump(adj, open(f"{args.outpath}/Cyber_data/contig.graph", "wb" ) )
    pkl.dump(test_to_id, open(f"{args.outpath}/Cyber_data/contig.dict", "wb" ) )



if mode == "testing":
    test_mask = []
    label = []
    cnt = 0
    for node in G.nodes():
        try:
            label.append(class_to_label[contig_to_family[node]])
            cnt+=1
        except:
            if "_" in node:
                try:
                    label.append(-1)
                    test_mask.append(cnt)
                    test_to_id[node] = cnt
                    cnt+=1
                except:
                    print(node)
            else:
                print(node)
    pkl.dump(test_mask, open(f"{args.outpath}/Cyber_data/contig.mask", "wb" ) )
    adj = nx.adjacency_matrix(G)
    pkl.dump(adj, open(f"{args.outpath}/Cyber_data/contig.graph", "wb" ) )
    pkl.dump(test_to_id, open(f"{args.outpath}/Cyber_data/contig.dict", "wb" ) )


# contructing feature map
fn = "database"
contig_feature = pkl.load(open(f"{args.outpath}/Cyber_data/contig.F",'rb'))
database_feature = pkl.load(open(fn+"/dataset_compressF",'rb'))
contig_to_id = joblib.load(open(f"{args.outpath}/out/contig_to_id.joblib", "rb"))


feature = []
for node in G.nodes():
    if "~" not in node:
        idx = contig_to_id[node]
        feature.append(contig_feature[idx])
    else:
        try:
            idx = int(name_to_id[node])
            feature.append(database_feature[idx])
        except:
            print(node,"error")

feature = np.array(feature)
if mode == "testing":
    pkl.dump(feature, open(f"{args.outpath}/Cyber_data/contig.feature", "wb" ) )
else:
    pkl.dump(feature, open(f"{args.outpath}/Cyber_data/contig.feature", "wb" ) )


# Graph check for each testing samples
cnt = 0
for node in G.nodes:
    flag = 0
    if "~" not in node:
        neighbor_label = []
        for edge in G.edges(node):
            neighbor = edge[1]
            if "~" in neighbor:
                neighbor_label.append(class_to_label[contig_to_family[neighbor]])
                flag =  1
        if flag == 0:
            label[test_to_id[node]] = -2
            #print(label[test_to_id[node]])
        if len(set(neighbor_label)) == 1:
            label[test_to_id[node]] = neighbor_label[0]
            cnt += 1



with open(f"{args.outpath}/network/phage_"+str(args.n)+".ntw","w") as out_f:
    for node in G.nodes:
        for edge in G.edges(node):
            neighbor = edge[1]
            out_f.write(node+","+neighbor+"\n")

pkl.dump(label, open(f"{args.outpath}/Cyber_data/contig.label", "wb" ) )
