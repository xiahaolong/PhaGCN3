import pandas as pd
import argparse
import os
import graph_tool.all as gt
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import rgb2hex
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from distinctipy import distinctipy
import colorsys
from collections import defaultdict

parser = argparse.ArgumentParser(description="Generate node.csv and edge.csv files.")
parser.add_argument('--outpath', type=str, default="result", help="Output directory path")
parser.add_argument('--degree_threshold', type=int, default=20, help='Minimum degree threshold for filtering nodes (default: 20)')
parser.add_argument('--top_n', type=int, default=20, help='Number of top taxa categories to visualize (default: 20)')
args = parser.parse_args()

if os.path.exists(f"{args.outpath}/tmp/edge.csv") and os.path.exists(f"{args.outpath}/tmp/node.csv"):
    print("Edge.csv and node.csv already exist. Skipping processing.")
else:
    os.makedirs(f"{args.outpath}/tmp", exist_ok=True)
    os.system(f"sed '1i source,target' {args.outpath}/final_network.ntw > {args.outpath}/tmp/edge.csv")
    
    os.system(f"cut -f 1 -d ',' {args.outpath}/final_network.ntw > {args.outpath}/tmp/line1.txt")
    os.system(f"cut -f 2 -d ',' {args.outpath}/final_network.ntw > {args.outpath}/tmp/line2.txt")
    os.system(f"cat {args.outpath}/tmp/line2.txt {args.outpath}/tmp/line1.txt | sort | uniq -c > {args.outpath}/tmp/node.csv")
    os.system(f"sed 's/^[ \t]*//' {args.outpath}/tmp/node.csv | awk -F' ' '{{print $1}}' > {args.outpath}/tmp/node_num.csv")
    os.system(f"sed 's/^[ \t]*//' {args.outpath}/tmp/node.csv | awk -F' ' '{{print $2}}' > {args.outpath}/tmp/node_id.csv")

    print(f"Output files are stored in {args.outpath}")

    edge_df = pd.read_csv(f'{args.outpath}/tmp/node_id.csv', header=None, names=['source'])
    edge_df['target'] = edge_df['source']
    edge_df['weight'] = 1
    edge_df.to_csv(f'{args.outpath}/tmp/node_1.csv', index=False)

    data = pd.read_csv(f'{args.outpath}/tmp/node_1.csv')

    taxonomic_label_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "database", "taxonomic_label.csv")
    final_prediction_path = f"{args.outpath}/final_prediction.csv"
    taxonomic_path_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "database", "taxonomic_path.csv")
    edge_num_path = f"{args.outpath}/tmp/node_num.csv"

    taxonomic_label = pd.read_csv(taxonomic_label_path)
    taxonomic_label.set_index("contig_id", inplace=True)

    pred_to_label = {199:"Belpaoviridae",103:"Bacilladnaviridae",137:"Closteroviridae",36:"Azeredovirinae",200:"Caulimoviridae",59:"Queuovirinae",127:"Pseudototiviridae",10:"Autographiviridae",194:"Potyviridae",112:"Nanoviridae",217:"Polydnaviriformidae",0:"Lipothrixviridae",38:"Bronfenbrennervirinae",126:"Fusagraviridae",33:"Zobellviridae",114:"Pecoviridae",179:"Fusariviridae",104:"Circoviridae",20:"Mesyanzhinovviridae",108:"Gandrviridae",115:"Naryaviridae",107:"Smacoviridae",2:"Alloherpesviridae",208:"Asfarviridae",31:"Vilmaviridae",125:"Chrysoviridae",66:"Alexandravirus",37:"Bclasvirinae",14:"Demerecviridae",170:"Peribunyaviridae",47:"Guarnerosvirinae",80:"Ignaciovirus",175:"Phenuiviridae",101:"Papillomaviridae",152:"Solspiviridae",207:"Iridoviridae",142:"Alphaflexiviridae",176:"Orthomyxoviridae",190:"Picornaviridae",145:"Flaviviridae",99:"Microviridae",164:"Nyamiviridae",72:"Casadabanvirus",7:"Hafunaviridae",13:"Chimalliviridae",62:"Stephanstirmvirinae",71:"Benedictvirus",94:"Veracruzvirus",118:"Kanorauviridae",193:"Solemoviridae",106:"Vilyaviridae",215:"Anelloviridae",18:"Herelleviridae",162:"Lispiviridae",148:"Mitoviridae",75:"Efquatrovirus",122:"Geplanaviridae",65:"Weiservirinae",149:"Atkinsviridae",166:"Rhabdoviridae",60:"Sepvirinae",189:"Marnaviridae",210:"Tectiviridae",213:"Nudiviridae",16:"Grimontviridae",214:"Nimaviridae",159:"Artoviridae",24:"Rountreeviridae",8:"Ackermannviridae",91:"Skunavirus",27:"Schitoviridae",30:"Umezonoviridae",48:"Guernseyvirinae",156:"Qinviridae",32:"Zierdtviridae",120:"Geminiviridae",51:"Joanripponvirinae",198:"Hepadnaviridae",163:"Mymonaviridae",82:"Kostyavirus",83:"Kroosvirus",181:"Partitiviridae",169:"Hantaviridae",12:"Chaseviridae",43:"Gclasvirinae",88:"Paclarkvirus",110:"Draupnirviridae",136:"Bromoviridae",143:"Betaflexiviridae",205:"Phycodnaviridae",76:"Fernvirus",63:"Tybeckvirinae",146:"Nodaviridae",109:"Ouroboviridae",64:"Vequintavirinae",105:"Endolinaviridae",42:"Eucampyvirinae",6:"Suoliviridae",81:"Korravirus",34:"Arquatrovirinae",78:"Gladiatorvirus",17:"Guelinviridae",121:"Genomoviridae",183:"Coronaviridae",61:"Skurskavirinae",128:"Spiciviridae",138:"Endornaviridae",165:"Paramyxoviridae",46:"Gracegardnervirinae",53:"Kutznervirinae",4:"Intestiviridae",11:"Casjensviridae",95:"Vividuovirus",206:"Mimiviridae",102:"Parvoviridae",124:"Botybirnaviridae",41:"Dolichocephalovirinae",23:"Pootjesviridae",157:"Aliusviridae",49:"Hendrixvirinae",70:"Beenievirus",178:"Curvulaviridae",25:"Saffermanviridae",216:"Fuselloviridae",67:"Andromedavirus",3:"Orthoherpesviridae",92:"Triavirus",58:"Pclasvirinae",86:"Mudcatvirus",184:"Mesoniviridae",39:"Ceeclamvirinae",26:"Salasmaviridae",151:"Fiersviridae",73:"Ceduovirus",150:"Duinviridae",158:"Chuviridae",173:"Arenaviridae",202:"Pseudoviridae",90:"Pbunavirus",15:"Drexlerviridae",9:"Aliceevansviridae",35:"Azeevirinae",56:"Nymbaxtervirinae",195:"Astroviridae",180:"Hypoviridae",123:"Pleolipoviridae",182:"Arteriviridae",129:"Artiviridae",55:"Nclasvirinae",68:"Attisvirus",119:"Mahapunaviridae",133:"Sedoreoviridae",161:"Filoviridae",211:"Adenoviridae",45:"Gorskivirinae",87:"Obolenskvirus",188:"Iflaviridae",187:"Dicistroviridae",79:"Gordonvirus",203:"Retroviridae",174:"Nairoviridae",201:"Metaviridae",116:"Adamaviridae",135:"Hepeviridae",28:"Stanwilliamsviridae",186:"Caliciviridae",212:"Baculoviridae",52:"Jondennisvirinae",168:"Fimoviridae",130:"Inseviridae",84:"Marthavirus",209:"Poxviridae",139:"Kitaviridae",155:"Botourmiaviridae",97:"Wizardvirus",160:"Bornaviridae",96:"Wbetavirus",153:"Blumeviridae",167:"Xinmoviridae",134:"Spinareoviridae",111:"Anicreviridae",132:"Orthototiviridae",5:"Steigviridae",113:"Redondoviridae",144:"Tymoviridae",191:"Polycipiviridae",141:"Virgaviridae",131:"Lebotiviridae",89:"Pahexavirus",85:"Montyvirus",172:"Tospoviridae",69:"Backyardiganvirus",54:"Mccleskeyvirinae",192:"Secoviridae",100:"Polyomaviridae",98:"Inoviridae",19:"Kyanoviridae",74:"Dhillonvirus",197:"Birnaviridae",50:"Jameshumphriesvirinae",21:"Orlajensenviridae",147:"Tombusviridae",1:"Rudiviridae",204:"Polymycoviridae",196:"Yadokariviridae",177:"Amalgaviridae",154:"Steitzviridae",140:"Togaviridae",29:"Straboviridae",44:"Gordonclarkvirinae",57:"Ounavirinae",117:"Kirkoviridae",171:"Phasmaviridae",22:"Peduoviridae",185:"Tobaniviridae",93:"Turbidovirus",40:"Deejayvirinae",77:"Fromanvirus"}

    final_prediction = pd.read_csv(final_prediction_path, usecols=["contig_name", "idx", "prediction", "prediction_score", "full_path"])
    final_prediction.set_index("contig_name", inplace=True)

    taxonomic_path = pd.read_csv(taxonomic_path_file, sep='\t', header=None, names=["family", "path"])

    edge_num = pd.read_csv(edge_num_path, usecols=[0], header=None)

    def get_family_and_realm(source):
        attribute = "" 
        realm = "" 

        if source.startswith("test"):
            contig_name = "_".join(source.split("_", 2)[2:])
            if contig_name in final_prediction.index:
                family = final_prediction.loc[contig_name, "prediction"]
                if "like" in family.lower() and family:
                    attribute = "like"
                elif family:
                    attribute = "test"
                if pd.notna(final_prediction.loc[contig_name, "full_path"]):
                    path = final_prediction.loc[contig_name, "full_path"]
                    realm = path.split(";")[1] if ";" in path else ""
                return family, attribute, realm
            else:
                return "", "", "" 
        else:
            if source in taxonomic_label.index:
                pred_key = taxonomic_label.loc[source, "family"]
                if isinstance(pred_key, pd.Series):
                    pred_key = pred_key.iloc[0]
                family = pred_to_label.get(pred_key, "")
                attribute = "database" if family else ""
            
                taxonomic_row = taxonomic_path[taxonomic_path["family"] == family]
                if not taxonomic_row.empty:
                    path = taxonomic_row.iloc[0]["path"]
                    realm = path.split(";")[1].split(",")[0] if ";" in path else ""
                return family, attribute, realm
            else:
                return "", "", ""

    node_data = pd.DataFrame({
        "id": data["source"],
        "Attribute": "", 
        "Family": "",
        "Realm": "",  
        "size": ""
    })
    for idx, row in node_data.iterrows():
        family, attribute, realm = get_family_and_realm(row["id"])
        node_data.at[idx, "Family"] = family
        node_data.at[idx, "Attribute"] = attribute
        node_data.at[idx, "Realm"] = realm
        node_data["size"] = edge_num[0].values
    output_file = f"{args.outpath}/tmp/node.csv"
    node_data.to_csv(output_file, index=False)

    print(f"Output file node.csv as {output_file}")

print("Reading edge and node data...")
edges_df = pd.read_csv(f"{args.outpath}/tmp/edge.csv", memory_map=True)
nodes_df = pd.read_csv(f"{args.outpath}/tmp/node.csv", sep=",")
print("Data read complete.")

if not {'source', 'target'}.issubset(edges_df.columns):
    raise ValueError("CSV file must contain 'source' and 'target' columns.")
else:
    print("CSV file contains 'source' and 'target' columns.")
    
edge_file = f"{args.outpath}/tmp/edge.csv"
num_edges = sum(1 for _ in open(edge_file)) - 1  
print(f"Detected {num_edges} edges, adjusting parameters...")
if num_edges < 2_000_000:
    C, gamma, p, theta, epsilon, max_iter = 2.0, 0.2, 2.0, 0.6, 0.01, 0
elif num_edges < 10_000_000:
    C, gamma, p, theta, epsilon, max_iter = 4.0, 1.0, 1.5, 0.6, 0.01, 0
elif num_edges < 50_000_000:
    C, gamma, p, theta, epsilon, max_iter = 6.0, 3.0, 1.0, 0.6, 0.01, 0
else:
    C, gamma, p, theta, epsilon, max_iter = 10.0, 5.0, 0.5, 0.6, 0.01, 0
print(f"Parameters set: C={C}, gamma={gamma}, p={p}, theta={theta}, epsilon={epsilon}, max_iter={max_iter}")

def filter_low_degree_nodes(edges_df, degree_threshold=3):
    while True:
        node_degrees = edges_df['source'].value_counts().add(edges_df['target'].value_counts(), fill_value=0)
        valid_nodes = set(node_degrees[node_degrees >= degree_threshold].index)
        before_count = len(edges_df)
        edges_df = edges_df[edges_df['source'].isin(valid_nodes) & edges_df['target'].isin(valid_nodes)]
        after_count = len(edges_df)
        if before_count == after_count:
            break
    return edges_df, valid_nodes

edges_df, valid_nodes = filter_low_degree_nodes(edges_df, degree_threshold=args.degree_threshold)
print(f"Remaining edges: {len(edges_df)}, nodes: {len(valid_nodes)}")

print("Computing color mapping for taxa...")
taxa = "Family"
top_n = args.top_n
default_color = "#c0c0c0"

counts = nodes_df[taxa].value_counts()
top_categories = counts.nlargest(top_n).index

candidate_colors = distinctipy.get_colors(top_n * 2)

def is_not_gray(color):
    h, s, v = colorsys.rgb_to_hsv(*color)
    return s > 0.2 and v > 0.2  

filtered_colors = [c for c in candidate_colors if is_not_gray(c)]
if len(filtered_colors) < top_n:
    raise ValueError(f"Not enough usable colors after filtering, got {len(filtered_colors)}.")

palette = [rgb2hex(c) for c in filtered_colors[:top_n]]
color_map = {category: palette[i] for i, category in enumerate(top_categories)}

nodes_df["Mapped Color"] = nodes_df[taxa].apply(lambda x: color_map.get(x, default_color))
node_colors = dict(zip(nodes_df["id"], nodes_df["Mapped Color"]))

print("Color mapping complete.")
print("Building the graph...")
G = gt.Graph(directed=False)
vertex_map = {}
color_prop = G.new_vertex_property("string")

for node_id in valid_nodes:
    v = G.add_vertex()
    vertex_map[node_id] = v
    color_prop[v] = node_colors.get(node_id, default_color)

edges_df['weight'] = 1  
weight_prop = G.new_edge_property("float")  

for _, row in edges_df.iterrows():
    e = G.add_edge(vertex_map[row['source']], vertex_map[row['target']])
    weight_prop[e] = row['weight']

G.ep["weight"] = weight_prop  
print("Graph construction complete.")

print("Computing connected components...")
comp, hist = gt.label_components(G)
largest_component = max(hist)
print(f"Connected component calculation complete. Size of the largest component: {largest_component}")

print("Computing initial layout...")
pos = gt.sfdp_layout(G, C=C, gamma=gamma, p=p, theta=theta, eweight=G.ep["weight"], epsilon=epsilon, max_iter=max_iter)

component_positions = {}
for v in G.vertices():
    c = comp[v]
    if c not in component_positions:
        component_positions[c] = []
    component_positions[c].append(pos[v])

component_centers = {c: np.mean(component_positions[c], axis=0) for c in component_positions}

print("Performing subgraph boundary collision detection and repulsion...")

subgraph_bounds = {}
for comp_id in np.unique(comp.a):
    sub_coords = np.array([pos[v] for v in G.vertices() if comp[v] == comp_id])
    min_x, min_y = sub_coords.min(axis=0)
    max_x, max_y = sub_coords.max(axis=0)
    subgraph_bounds[comp_id] = {
        "coords": sub_coords,
        "bbox": [min_x, max_x, min_y, max_y],
        "shift": np.array([0.0, 0.0])
    }

max_iter = 1000
for _ in range(max_iter):
    overlap = False
    for i, (id_i, box_i) in enumerate(subgraph_bounds.items()):
        for j, (id_j, box_j) in enumerate(subgraph_bounds.items()):
            if i >= j:
                continue
            min_x1, max_x1, min_y1, max_y1 = box_i["bbox"]
            min_x2, max_x2, min_y2, max_y2 = box_j["bbox"]

            if not (max_x1 < min_x2 or max_x2 < min_x1 or
                    max_y1 < min_y2 or max_y2 < min_y1):
                overlap = True
                center_i = np.mean(box_i["coords"], axis=0)
                center_j = np.mean(box_j["coords"], axis=0)
                direction = center_i - center_j
                norm = np.linalg.norm(direction)
                if norm == 0:
                    direction = np.array([1.0, 0.0])
                else:
                    direction = direction / norm

                push = 0.7  
                box_i["shift"] += push * direction
                box_j["shift"] -= push * direction

    for comp_id, box in subgraph_bounds.items():
        shift = box["shift"]
        if np.linalg.norm(shift) > 0:
            for v in G.vertices():
                if comp[v] == comp_id:
                    pos[v] = (pos[v][0] + shift[0], pos[v][1] + shift[1])
            new_coords = np.array([pos[v] for v in G.vertices() if comp[v] == comp_id])
            min_x, min_y = new_coords.min(axis=0)
            max_x, max_y = new_coords.max(axis=0)
            subgraph_bounds[comp_id]["coords"] = new_coords
            subgraph_bounds[comp_id]["bbox"] = [min_x, max_x, min_y, max_y]
            subgraph_bounds[comp_id]["shift"] = np.array([0.0, 0.0])
    
    if not overlap:
        break

print("Subgraph boundary repulsion complete.")

print("Drawing and saving network graph...")
output_path = f"{args.outpath}/optimized_network_graph_compact.pdf"
with PdfPages(output_path) as pdf:
    fig, ax = plt.subplots(figsize=(14, 12), dpi=500)
    
    gt.graph_draw(G, pos=pos, 
                  output_size=(5000, 5000), 
                  vertex_fill_color=color_prop,  
                  edge_pen_width=0.1,  
                  edge_color="lightgray",
                  vertex_size=0.5,
                  mplfig=ax)

    ax.set_title("Optimized Network Graph (Compact Layout)", fontsize=16)

    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color,
                                   markersize=8, label=category) for category, color in color_map.items()]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1),
              fontsize=8, title="Taxa Categories", title_fontsize=10, frameon=True)
    
    pdf.savefig(fig, bbox_inches="tight", dpi=500)
    plt.close(fig)

print(f"Output saved to '{output_path}'")

print("Computing color mapping for Attribute...")
attribute = "Attribute"  
top_attribute = 3  

attribute_counts = nodes_df[attribute].value_counts()
top_attribute_categories = attribute_counts.nlargest(top_attribute).index
attribute_palette = sns.color_palette("Set2", n_colors=len(top_attribute_categories))

attribute_color_map = {category: rgb2hex(attribute_palette[i]) for i, category in enumerate(top_attribute_categories)}

nodes_df["Attribute Color"] = nodes_df[attribute].apply(lambda x: attribute_color_map[x] if x in attribute_color_map else default_color)

attribute_node_colors = dict(zip(nodes_df["id"], nodes_df["Attribute Color"]))
print("Color mapping for Attribute complete.")

print("Drawing and saving network graph with Attribute coloring...")
output_path_attribute = f"{args.outpath}/optimized_network_graph_attribute_colored.pdf"
with PdfPages(output_path_attribute) as pdf:
    fig, ax = plt.subplots(figsize=(14, 12), dpi=500)

    attribute_color_prop = G.new_vertex_property("string")
    for node_id in valid_nodes:
        v = vertex_map[node_id]
        attribute_color_prop[v] = attribute_node_colors.get(node_id, default_color)

    gt.graph_draw(G, pos=pos, 
                  output_size=(5000, 5000), 
                  vertex_fill_color=attribute_color_prop,  
                  edge_pen_width=0.1,  
                  edge_color="lightgray",
                  vertex_size=0.5,
                  mplfig=ax)

    ax.set_title("Network Graph with Attribute Coloring", fontsize=16)

    attribute_legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color,
                                            markersize=8, label=category) for category, color in attribute_color_map.items()]
    ax.legend(handles=attribute_legend_elements, loc='upper left', bbox_to_anchor=(1, 1),
              fontsize=8, title="Attribute Categories", title_fontsize=10, frameon=True)

    pdf.savefig(fig, bbox_inches="tight", dpi=500)
    plt.close(fig)

print(f"Output saved to '{output_path_attribute}'")