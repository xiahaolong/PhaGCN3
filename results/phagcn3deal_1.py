import pandas as pd

# 定义分类映射
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

# 原始 CSV 文件路径
input_csv = "final_prediction_1.csv"
# 输出 CSV 文件路径
output_csv = "final_prediction_1_out.csv"

# 读取原始文件
data = pd.read_csv(input_csv)

# 初始化新的 DataFrame
new_columns = []
for key, value in classification_map.items():
    new_columns.append(value)
    new_columns.append(f"{value.split(' ')[0]}_score")

# 创建一个空的 DataFrame
output_data = pd.DataFrame(columns=["contig_name"] + new_columns)

# 填充数据
rows = []  # 用于暂存新行
for index, row in data.iterrows():
    contig_name = row["contig_name"]
    classifications = row[1:].dropna().tolist()
    
    # 初始化新行
    new_row = {col: "NA" for col in output_data.columns}
    new_row["contig_name"] = contig_name
    
    for entry in classifications:
        # 提取分类标识和名称
        parts = entry.split(';', 1)
        if len(parts) == 2:
            key, value = parts
            column_name = classification_map.get(key, None)
            if column_name:
                new_row[column_name] = value
    
    rows.append(new_row)  # 将新行添加到暂存列表

# 将所有新行一次性添加到 DataFrame
output_data = pd.concat([output_data, pd.DataFrame(rows)], ignore_index=True)

# 保存到 CSV 文件
output_data.to_csv(output_csv, index=False)