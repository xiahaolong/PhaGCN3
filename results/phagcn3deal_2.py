import pandas as pd

input_csv = "final_prediction_1_out.csv"
input3_csv = "score.csv"
output_csv = "ICTVTaxoChallenge.csv"

data = pd.read_csv(input_csv)
data2 = pd.read_csv(input3_csv)

for index, row in data2.iterrows():
    contig_name = row["contig_name"]
    prediction = row["prediction"]
    prediction_score = row["prediction_score"]

    matching_rows = data[data["contig_name"] == contig_name]

    if not matching_rows.empty:
        for idx in matching_rows.index:
            for col in data.columns:
                if data.at[idx, col] == prediction:
                    score_col = f"{col.split(' ')[0]}_score"
                    if score_col in data.columns:
                        data.at[idx, score_col] = prediction_score

data.fillna("NA", inplace=True)

data.to_csv(output_csv, index=False)
