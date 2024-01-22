# Take .csv of WURCS and returns glycan tree type (using glycan_tree_type_identifier)
import argparse
import csv
from glycan_tree_type_identifier import check_type

def process_csv(input_csv, output_csv):
    with open(input_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        header = next(reader, None)  # Read the header if present

        # Add a new column 'Results' to the header
        header_with_results = header + ['Results']

        writer = csv.writer(outfile)
        writer.writerow(header_with_results)  # Write the modified header

        for row in reader:
            wurcs_code = row[0]
            user_wurcs = f'"{wurcs_code}"'
            result = check_type(user_wurcs)
            
            # Add the result to the row
            row_with_result = row + [result]
            
            writer.writerow(row_with_result)  # Write the row with the result

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="glycan_composition_identification",
        description="""Identify whether your glycan tree is high mannose, complex or hybrid."""
    )
    parser.add_argument(
        "-i",
        "--input_csv",
        help="Path to the input CSV file containing WURCS codes"
    )
    parser.add_argument(
        "-o",
        "--output_csv",
        help="Path to the output CSV file with added 'Results' column"
    )

    args = parser.parse_args()

    if args.input_csv and args.output_csv:
        process_csv(args.input_csv, args.output_csv)
    else:
        print("Please provide paths to the input and output CSV files using -i/--input_csv and -o/--output_csv options.")
