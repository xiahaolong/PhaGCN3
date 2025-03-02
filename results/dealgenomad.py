import csv
import sys

def process_csv(input_file, output_file):
    with open(input_file, newline='', encoding='utf-8') as csv_in, \
         open(output_file, 'w', newline='', encoding='utf-8') as csv_out:
        reader = csv.reader(csv_in)
        writer = csv.writer(csv_out)
        
        for row in reader:
            last_nonempty = ""
            for cell in reversed(row):
                if cell.strip() != "":
                    last_nonempty = cell
                    break
            
            if len(row) < 2:
                row.extend([""] * (2 - len(row)))
            
            row[1] = last_nonempty
            
            writer.writerow(row)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: python dealgenomad.py processed_test_nodes_taxonomy.csv processed_test_nodes_taxonomy_output.csv")
        sys.exit(1)
        
    input_csv = sys.argv[1]
    output_csv = sys.argv[2]
    process_csv(input_csv, output_csv)
