import networkx as nx
import csv
import re
import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--contigs', type=str, default='contigs.fa')
parser.add_argument('--outpath', type=str, default="result")
args = parser.parse_args()

def extract_subgraphs_and_process_nodes(input_csv, output_edge_csv, output_node_csv, output_subgraph_node_csv):
    G = nx.Graph()

    with open(input_csv, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if len(row) >= 2:
                G.add_edge(row[0], row[1])

    connected_components = list(nx.connected_components(G))
    filtered_edges = []
    filtered_nodes = set()  
    subgraph_nodes = []
    
    for i, component in enumerate(connected_components, start=1):
        subgraph = G.subgraph(component)
        if all(edge[0].startswith("test_") and edge[1].startswith("test_") for edge in subgraph.edges()):
            filtered_edges.extend(subgraph.edges())
            filtered_nodes.update(subgraph.nodes())
            subgraph_nodes.extend([(f"subgraph{i}", re.sub(r'^test_\d+_', '', node)) for node in subgraph.nodes()])


    processed_nodes = {
        re.sub(r'^test_\d+_', '', node) for node in filtered_nodes
    }

    with open(output_edge_csv, 'w') as file:
        writer = csv.writer(file)
        writer.writerows(filtered_edges)

    with open(output_node_csv, 'w') as file:
        writer = csv.writer(file)
        for node in processed_nodes:
            writer.writerow([node])

    with open(output_subgraph_node_csv, 'w') as file:
        writer = csv.writer(file)
        subgraph_counter = 1  
        for component in connected_components:
            subgraph = G.subgraph(component)
            if all(edge[0].startswith("test_") and edge[1].startswith("test_") for edge in subgraph.edges()):
                for node in subgraph.nodes():
                    writer.writerow([f"subgraph{subgraph_counter}", re.sub(r'^test_\d+_', '', node)])
                subgraph_counter += 1  

    print(f"The edges of the subgraph that meet the conditions have been saved to {output_edge_csv}.")
    print(f"The processed subgraph nodes have been saved to {output_node_csv}.")
    print(f"The list of nodes for each subgraph has been saved to {output_subgraph_node_csv}.")

input_csv_path = f"{args.outpath}/final_network.ntw"  
output_edge_csv_path = f"{args.outpath}/filtered_test_edges.csv"  
output_node_csv_path = f"{args.outpath}/processed_test_nodes.csv"  
output_subgraph_node_csv_path = f"{args.outpath}/subgraph_nodes.csv" 

extract_subgraphs_and_process_nodes(input_csv_path, output_edge_csv_path, output_node_csv_path, output_subgraph_node_csv_path)

if os.path.exists(output_node_csv_path) and os.stat(output_node_csv_path).st_size == 0:
    print(f"{output_node_csv_path} is empty. Skipping genomad.")
    exit(0)  


def remove_matched_rows(nodes_file, predictions_file):
    with open(nodes_file, 'r') as file:
        nodes = set(line.strip() for line in file)

    remaining_rows = []

    with open(predictions_file, 'r') as file:
        reader = csv.reader(file)
        header = next(reader)  
        remaining_rows.append(header)  
        
        for row in reader:
            contig_name = row[0]  

            if contig_name not in nodes:
                remaining_rows.append(row)

    with open(predictions_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(remaining_rows)

    print(f"The file {predictions_file} has been successfully updated, with the matching rows deleted!")

nodes_file_path = f"{args.outpath}/processed_test_nodes.csv"  
predictions_file_path = f"{args.outpath}/final_prediction.csv"  

remove_matched_rows(nodes_file_path, predictions_file_path)

try:
    cmd = (
    f"seqkit grep -f {args.outpath}/processed_test_nodes.csv "
    f"{args.contigs} > {args.outpath}/processed_test_nodes.fasta"
)
    out = subprocess.check_call(cmd, shell=True)
    cmd1 = f"genomad annotate {args.outpath}/processed_test_nodes.fasta {args.outpath}/genomad_output genomad_db"
    out = subprocess.check_call(cmd1, shell=True)
    cmd2 = f"mv {args.outpath}/genomad_output/processed_test_nodes_annotate/processed_test_nodes_taxonomy.tsv {args.outpath}"
    out = subprocess.check_call(cmd2, shell=True)
except Exception as e:
    print(f"genomad error: {e}")




input_file = f"{args.outpath}/processed_test_nodes_taxonomy.tsv"


modifications = {
    "Viruses;Riboviria;Orthornavirae;Negarnaviricota;Ellioviricetes;Bunyavirales;Phenuiviridae": "Viruses;Riboviria;Orthornavirae;Negarnaviricota;Bunyaviricetes;Elliovirales;Phenuiviridae",
    "Viruses;Riboviria;Orthornavirae;Negarnaviricota;Ellioviricetes;Bunyavirales;Arenaviridae": "Viruses;Riboviria;Orthornavirae;Negarnaviricota;Bunyaviricetes;Elliovirales;Arenaviridae",
    "Viruses;Riboviria;Orthornavirae;Negarnaviricota;Ellioviricetes;Bunyavirales;Peribunyaviridae": "Viruses;Riboviria;Orthornavirae;Negarnaviricota;Bunyaviricetes;Elliovirales;Peribunyaviridae",
    "Viruses;Riboviria;Orthornavirae;Negarnaviricota;Ellioviricetes;Bunyavirales;Nairoviridae": "Viruses;Riboviria;Orthornavirae;Negarnaviricota;Bunyaviricetes;Elliovirales;Nairoviridae",
    "Viruses;Duplodnaviria;Heunggongvirae;Peploviricota;Herviviricetes;Herpesvirales;Herpesviridae": "Viruses;Duplodnaviria;Heunggongvirae;Peploviricota;Herviviricetes;Herpesvirales;",
    "Viruses;Varidnaviria;Bamfordvirae;Preplasmiviricota;Maveriviricetes;Priklausovirales;Lavidaviridae": "Viruses;Varidnaviria;Bamfordvirae;Preplasmiviricota;Maveriviricetes;Lavidavirales;"
}


with open(input_file, "r", encoding="utf-8") as file:
    lines = file.readlines()

with open(input_file, "w", encoding="utf-8") as file:
    for line in lines:
        columns = line.strip().split("\t")  
        for original, modified in modifications.items():
            columns = [re.sub(re.escape(original), modified, col) for col in columns]
        file.write("\t".join(columns) + "\n")
