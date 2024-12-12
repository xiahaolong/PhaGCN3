import argparse
import os
import pandas as pd
import igraph as ig
import matplotlib.pyplot as plt
from matplotlib.colors import rgb2hex
from matplotlib.backends.backend_pdf import PdfPages
from fa2 import ForceAtlas2
from scipy.sparse import coo_matrix
import seaborn as sns


parser = argparse.ArgumentParser(description="Generate node.csv and edge.csv files.")
parser.add_argument('--outpath', type=str, default="result", help="Output directory path")
args = parser.parse_args()

os.makedirs(f"{args.outpath}/tmp", exist_ok=True)
os.system(f"sed '1i source,target' {args.outpath}/final_network.ntw > {args.outpath}/tmp/final_network_modified.ntw")
os.system(f"cat {args.outpath}/tmp/final_network_modified.ntw > {args.outpath}/tmp/edge.csv")

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
        # print(contig_name)
        if contig_name in final_prediction.index:
            family = final_prediction.loc[contig_name, "prediction"]
            if "like" in family.lower() and family:
                attribute = "like"
            elif family:
                attribute = "test"
            if pd.notna(final_prediction.loc[contig_name, "full_path"]):
                path = final_prediction.loc[contig_name, "full_path"]
                realm = path.split(";")[1] if ";" in path else ""
                # print(realm)
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
                # print(path)
                realm = path.split(";")[1].split(",")[0] if ";" in path else ""
                # print(realm)
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

print(f"output file node.csv as {output_file}")

edges_df = pd.read_csv(f"{args.outpath}/tmp/edge.csv")
nodes_df = pd.read_csv(f"{args.outpath}/tmp/node.csv", sep=",")

if not {'source', 'target'}.issubset(edges_df.columns):
    raise ValueError("CSV file node 'source' and 'target' column.")

taxa = "Family"
top_n = 10 
default_color = "#c0c0c0"

counts = nodes_df[taxa].value_counts()
top_categories = counts.nlargest(top_n).index
categories = set(top_categories)

palette = sns.color_palette("tab10", n_colors=top_n)
color_map = {category: rgb2hex(palette[i]) for i, category in enumerate(top_categories)}

nodes_df['Mapped Color'] = nodes_df[taxa].map(color_map).fillna(default_color)
node_colors = {node_id: color for node_id, color in zip(nodes_df['id'], nodes_df['Mapped Color'])}

edges_df['weight'] = 1
edges_weighted = edges_df.groupby(['source', 'target'], as_index=False).sum()

G = ig.Graph()
G.add_vertices(list(set(edges_df['source']).union(set(edges_df['target']))))
edges = list(zip(edges_weighted['source'], edges_weighted['target']))
G.add_edges(edges)
G.es['weight'] = edges_weighted['weight']

node_indices = {v: i for i, v in enumerate(G.vs['name'])}
rows = [node_indices[edge[0]] for edge in edges]
cols = [node_indices[edge[1]] for edge in edges]
weights = G.es['weight']
matrix = coo_matrix((weights, (rows, cols)), shape=(len(G.vs), len(G.vs)))

forceatlas2 = ForceAtlas2(
    outboundAttractionDistribution=True,
    linLogMode=False,
    adjustSizes=False,
    edgeWeightInfluence=1.0,
    jitterTolerance=1.0,
    barnesHutOptimize=True,
    barnesHutTheta=1.2,
    scalingRatio=1.0,
    strongGravityMode=False,
    gravity=300.0
)
positions = forceatlas2.forceatlas2(matrix, iterations=5000)
layout = [(pos[0], pos[1]) for pos in positions]

G.vs['color'] = [node_colors.get(node['name'], default_color) for node in G.vs]

output = f"{args.outpath}/forceatlas2_network_graph_colored_2.pdf"

with PdfPages(output) as pdf:
    fig, ax = plt.subplots(figsize=(14, 12))
    ig.plot(
        G, layout=layout, target=ax, vertex_size=5,
        vertex_color=G.vs['color'], edge_color="gray",
        edge_width=[0.01 + 0.5 * w for w in G.es['weight']]
    )
    ax.set_title("ForceAtlas2 Network Graph (Top 10 Colored by Family)", fontsize=16)

    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color,
                                   markersize=8, label=category) for category, color in color_map.items()]
    ax.legend(
        handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1),
        fontsize=8, title="Top 10 Taxa Categories", title_fontsize=10, frameon=True
    )
    pdf.savefig(fig, bbox_inches="tight") 
    plt.close(fig)

print(f"'{output}'")

taxa_attribute = "Attribute"
default_color_attribute = "#c0c0c0"
categories_attribute = set(nodes_df[taxa_attribute])
counts_attribute = nodes_df[taxa_attribute].value_counts()

palette_attribute = sns.color_palette("Set2", n_colors=len(categories_attribute))
color_map_attribute = {category: rgb2hex(palette_attribute[i]) for i, category in enumerate(categories_attribute)}

nodes_df['Attribute Color'] = nodes_df[taxa_attribute].map(color_map_attribute).fillna(default_color_attribute)
node_colors_attribute = {node_id: color for node_id, color in zip(nodes_df['id'], nodes_df['Attribute Color'])}

G.vs['color_attribute'] = [node_colors_attribute.get(node['name'], default_color_attribute) for node in G.vs]

output_path_attribute = f"{args.outpath}/forceatlas2_network_graph_colored_attribute.pdf"

with PdfPages(output_path_attribute) as pdf:
    fig, ax = plt.subplots(figsize=(14, 12))
    ig.plot(
        G, layout=layout, target=ax, vertex_size=5,
        vertex_color=G.vs['color_attribute'], edge_color="gray",
        edge_width=[0.01 + 0.5 * w for w in G.es['weight']]
    )
    ax.set_title("ForceAtlas2 Network Graph (Colored by Attribute)", fontsize=16)

    legend_elements_attribute = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color,
                   markersize=8, label=category)
        for category, color in color_map_attribute.items()
    ]
    ax.legend(
        handles=legend_elements_attribute, loc='upper left', bbox_to_anchor=(1, 1),
        fontsize=8, title="Attribute Categories", title_fontsize=10, frameon=True
    )
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

print(f"'{output_path_attribute}'")