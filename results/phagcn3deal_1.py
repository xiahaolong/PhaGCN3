import pandas as pd

classification_map = {
    'r': 'Realm (-viria)',
    'sr':'Subrealm (-vira)',
    'k': 'Kingdom (-virae)',
    'sk': 'Subkingdom (-virites)',
    'p': 'Phylum (-viricota)',
    'sp': 'Subphylum (-viricotina)',
    'c': 'Class (-viricetes)',
    'sc': 'Subclass (-viricetidae)',
    'o': 'Order (-virales)',
    'so': 'Suborder (-virineae)',
    'f': 'Family (-viridae)',
    'sf': 'Subfamily (-virinae)',
    'g': 'Genus (-virus)',
    'sg': 'Subgenus (-virus)',
    's': 'Species (binomial)'
}

input_csv = "final_prediction_1.csv"
output_csv = "final_prediction_1_out.csv"

data = pd.read_csv(input_csv)

new_columns = []
for key, value in classification_map.items():
    new_columns.append(value)
    new_columns.append(f"{value.split(' ')[0]}_score")

output_data = pd.DataFrame(columns=["contig_name"] + new_columns)

rows = []  
for index, row in data.iterrows():
    contig_name = row["contig_name"]
    classifications = row[1:].dropna().tolist()
    
    new_row = {col: "NA" for col in output_data.columns}
    new_row["contig_name"] = contig_name
    
    for entry in classifications:
        parts = entry.split(';', 1)
        if len(parts) == 2:
            key, value = parts
            column_name = classification_map.get(key, None)
            if column_name:
                new_row[column_name] = value
    
    rows.append(new_row)  

output_data = pd.concat([output_data, pd.DataFrame(rows)], ignore_index=True)

output_data.to_csv(output_csv, index=False)