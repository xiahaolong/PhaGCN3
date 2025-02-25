import pandas as pd

# 定义文件路径
input_csv = "final_prediction_1_out.csv"
input3_csv = "score.csv"
output_csv = "ICTVTaxoChallenge.csv"

# 读取文件
data = pd.read_csv(input_csv)
data2 = pd.read_csv(input3_csv)

# 遍历 input3.csv 的每一行
for index, row in data2.iterrows():
    contig_name = row["contig_name"]
    prediction = row["prediction"]
    prediction_score = row["prediction_score"]

    # 找到 output2.csv 中 contig_name 匹配的行
    matching_rows = data[data["contig_name"] == contig_name]

    if not matching_rows.empty:
        for idx in matching_rows.index:
            # 遍历分类列，找到 prediction 匹配的列
            for col in data.columns:
                if data.at[idx, col] == prediction:
                    score_col = f"{col.split(' ')[0]}_score"
                    if score_col in data.columns:
                        data.at[idx, score_col] = prediction_score

# 填充空单元格为 "NA"
data.fillna("NA", inplace=True)

# 保存更新后的数据到新的 CSV 文件
data.to_csv(output_csv, index=False)
