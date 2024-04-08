import argparse
import csv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_file_names",
        type=str,
        required=True,
        nargs="+",
        help="input file names",
    )
    parser.add_argument(
        "-o", "--output_file_name", type=str, required=True, help="output file name"
    )
    args = parser.parse_args()
    input_file_names = args.input_file_names
    output_file_name = args.output_file_name
    primer_pairs = set()
    records_to_be_written = []
    with open(input_file_names[0]) as f:
        l = next(f).strip().split("\t")
        primer_pairs.add((l[1], l[2]))
        records_to_be_written.append(l)
    for each_file in input_file_names:
        with open(each_file) as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                primer_pair = (row[1], row[2])
                if primer_pair in primer_pairs:
                    continue
                primer_pairs.add(primer_pair)
                records_to_be_written.append(row)
    with open(output_file_name, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(records_to_be_written)


if __name__ == "__main__":
    main()
