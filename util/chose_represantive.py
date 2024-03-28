"""
複数のTSVファイルを受け取る
第一カラムのIDは重複があるが、その中から１つ選び、それ以外の行を削除する
argparseとcsvをimportする
"""

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
    id_values = set()
    records_to_be_written = []
    with open(input_file_names[0]) as f:
        l = next(f).strip().split("\t")
        id_values.add(l[0])
        records_to_be_written.append(l)
    for each_file in input_file_names:
        with open(each_file) as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                id_value = row[0]
                if id_value in id_values:
                    continue
                id_values.add(id_value)
                records_to_be_written.append(row)
    with open(output_file_name, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(records_to_be_written)
