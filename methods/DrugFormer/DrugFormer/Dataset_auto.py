import argparse
import csv
import os
import pickle
import sys

from datasets import Dataset

csv.field_size_limit(500000)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a DrugFormer expression matrix to HuggingFace Dataset format"
    )
    parser.add_argument(
        "--input_file",
        type=str,
        required=True,
        help="Path to the DrugFormer input expression matrix"
    )
    parser.add_argument(
        "--dict_file",
        type=str,
        required=True,
        help="Path to token_gene_dictionary.pickle"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory to save the processed dataset"
    )
    parser.add_argument(
        "--max_len",
        type=int,
        default=2048,
        help="Maximum sequence length (default: 2048)"
    )
    parser.add_argument(
        "--pos_keyword",
        type=str,
        default="sensitive",
        help="Keyword corresponding to label 1"
    )
    parser.add_argument(
        "--neg_keyword",
        type=str,
        default="resistant",
        help="Keyword corresponding to label 0"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if not os.path.exists(args.dict_file):
        print(f"[Error] Dictionary file not found: {args.dict_file}")
        sys.exit(1)

    if not os.path.exists(args.input_file):
        print(f"[Error] Input file not found: {args.input_file}")
        sys.exit(1)

    with open(args.dict_file, "rb") as file:
        dictionary_genename_token_pair = pickle.load(file)

    labels_cell = []
    samples = []

    print(f"Processing input file: {args.input_file}")

    with open(args.input_file, "r", encoding="utf-8") as file:
        csv_reader = csv.reader(file)
        pattern = []
        header = True

        for row in csv_reader:
            if not row:
                continue

            row_data = row[0].split("\t")

            if header:
                for gene in row_data:
                    gene = gene.replace('"', "")
                    if gene in dictionary_genename_token_pair:
                        pattern.append(dictionary_genename_token_pair[gene])
                    else:
                        pattern.append(-99999)
                header = False

            else:
                if len(pattern) != len(row_data):
                    print(
                        f"[Error] Pattern length ({len(pattern)}) does not match "
                        f"row length ({len(row_data)})."
                    )
                    sys.exit(1)

                seq_pattern_order_id_EXPscore = []

                for i in range(len(row_data)):
                    if i == 0:
                        pass

                    elif i == 1:
                        if args.pos_keyword in row_data[i]:
                            labels_cell.append(1)
                        elif args.neg_keyword in row_data[i]:
                            labels_cell.append(0)

                    else:
                        if row_data[i] == "0":
                            pass
                        else:
                            if pattern[i] == -99999:
                                pass
                            else:
                                seq_pattern_order_id_EXPscore.append(
                                    (pattern[i], row_data[i])
                                )

                sorted_seq_pattern = sorted(
                    seq_pattern_order_id_EXPscore,
                    key=lambda x: x[1],
                    reverse=True
                )

                sample = [item[0] for item in sorted_seq_pattern]

                while len(sample) < args.max_len:
                    sample.append(0)

                assert len(sample) >= args.max_len
                sample = sample[:args.max_len]

                samples.append(sample)

    if len(samples) != len(labels_cell):
        print(
            f"[Error] Sample-label mismatch: "
            f"{len(samples)} samples and {len(labels_cell)} labels."
        )
        sys.exit(1)

    dt = {
        "input_ids": samples,
        "cell_label": labels_cell
    }

    os.makedirs(os.path.dirname(args.output_dir), exist_ok=True)

    print(f"Saving dataset with {len(samples)} samples to: {args.output_dir}")

    my_dataset = Dataset.from_dict(dt)
    my_dataset.save_to_disk(args.output_dir)

    print("Dataset saved successfully.")


if __name__ == "__main__":
    main()
