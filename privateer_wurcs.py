# Extract WURCS code and sugar chain IDs for PDB files

import os
from privateer import privateer_core as pvt
import pandas as pd

directory = "/Users/lucyschofield/phd/xhpi/pdb_files"
output_csv_file_path = "WURCS_privateer_output.csv"

def get_sugar_id (totalWurcs_list, output_csv_file_path, file_name):
    ids = []
    wurcs_list = []
    
    for i in range(1, len(totalWurcs_list), 2):  # Start from the second line, skipping the first line
        ids.append(totalWurcs_list[i].strip())
        wurcs_list.append(totalWurcs_list[i + 1].strip())

    df = pd.DataFrame({"ID": ids, "WURCS": wurcs_list})
    df['TSChainId'] = df['ID'].apply(lambda x: x.split('_')[0][-1] if x.split('_')[0] else None)
    df['FileName'] = file_name
    df = df[['FileName', 'TSChainId', 'ID', 'WURCS']]
    df.to_csv(output_csv_file_path, mode='a', header=not os.path.exists(output_csv_file_path))


def get_wurcs (file_path, output_csv_file_path):
    file_name = os.path.splitext(os.path.basename(file_path))[0]
    totalWURCS = pvt.print_wurcs(file_path)
    totalWurcs_list = totalWURCS.splitlines()
    get_sugar_id(totalWurcs_list, output_csv_file_path, file_name)

if __name__ == "__main__":
    # List all files in the directory with a ".pdb" extension
    pdb_files = [f for f in os.listdir(directory) if f.endswith(".pdb")]
    for file in pdb_files:
        file_path = os.path.join(directory, file)
        get_wurcs(file_path, output_csv_file_path)
