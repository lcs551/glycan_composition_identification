import os
import signal
from privateer import privateer_core as pvt
import pandas as pd
from tqdm import tqdm

directory = "/Users/lcs551/phd/year_1/xhpi/glycan_composition_identification/data"
output_csv_file_path = "/Users/lcs551/phd/year_1/xhpi/glycan_composition_identification/data/delete_WURCS_privateer_output.csv"
error_output_directory = "/Users/lcs551/phd/year_1/xhpi/glycan_composition_identification/data/fail_outputs"

if not os.path.exists(error_output_directory):
    os.makedirs(error_output_directory)

class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException()

def get_sugar_id(totalWurcs_list, output_csv_file_path, file_name):
    ids = []
    wurcs_list = []

    for i in range(1, len(totalWurcs_list), 2):  # Start from the second line, skipping the first line
        ids.append(totalWurcs_list[i].strip())
        wurcs_list.append(totalWurcs_list[i + 1].strip())

    df = pd.DataFrame({"ID": ids, "WURCS": wurcs_list})
    df['TSChainId'] = df['ID'].apply(lambda x: x.split('_')[0][-1] if x.split('_')[0] else None)
    df['FileName'] = file_name
    df = df[['FileName', 'TSChainId', 'ID', 'WURCS']]
    
    try:
        df.to_csv(output_csv_file_path, mode='a', header=not os.path.exists(output_csv_file_path))
    except Exception as e:
        error_file_path = os.path.join(error_output_directory, f"fail_{file_name}.txt")
        with open(error_file_path, 'w') as file:
            file.write(f"CSV Write Error: {e}")

def get_wurcs(file_path, output_csv_file_path):
    file_name = os.path.splitext(os.path.basename(file_path))[0]

    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(600)  # 10 minutes

    try:
        totalWURCS = pvt.print_wurcs(file_path)
        signal.alarm(0)  # Disable the alarm
    except TimeoutException:
        error_file_path = os.path.join(error_output_directory, f"fail_{file_name}.txt")
        with open(error_file_path, 'w') as file:
            file.write("Timeout Error: Function call took too long")
        return
    except Exception as e:
        signal.alarm(0)  # Disable the alarm
        error_file_path = os.path.join(error_output_directory, f"fail_{file_name}.txt")
        with open(error_file_path, 'w') as file:
            file.write(f"{e}")
        return

    totalWurcs_list = totalWURCS.splitlines()
    if not totalWurcs_list:
        error_file_path = os.path.join(error_output_directory, f"fail_{file_name}.txt")
        with open(error_file_path, 'w') as file:
            file.write("Empty WURCS data")
        return

    try:
        get_sugar_id(totalWurcs_list, output_csv_file_path, file_name)
    except Exception as e:
        error_file_path = os.path.join(error_output_directory, f"fail_{file_name}.txt")
        with open(error_file_path, 'w') as file:
            file.write(f"{e}")

if __name__ == "__main__":
    # List all files in the directory with a ".pdb" extension
    pdb_files = [f for f in os.listdir(directory) if f.endswith(".pdb")]
    
    file_count = 0

    with tqdm(total=len(pdb_files), desc="Processing files", unit="file") as pbar:
        for file in pdb_files:
            file_count += 1
            file_name = file.replace(".pdb", "")
            file_path = os.path.join(directory, file)

            try:
                get_wurcs(file_path, output_csv_file_path)
            except Exception as e:
                error_file_path = os.path.join(error_output_directory, f"fail_{file_name}.txt")
                with open(error_file_path, 'w') as file:
                    file.write(f"{e}")
                break

            pbar.update(1)

    print(f"Processed {file_count} files")
