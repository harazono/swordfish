#! /usr/bin/env python3

# from primer3_result_parser import primer3_result_parser, extract_primer_pairs
from itertools import combinations
import sys
import argparse
import pprint
import json
import csv
pp = pprint.PrettyPrinter(indent=2)


def file_parser(input_file_name: str) -> list:
    csv.field_size_limit(int(1e6))
    # header付きのtsvを読み込み、最終カラムの文字列をjson.parseする。ヘッダーをキーとする辞書で返す。
    with open(input_file_name, "r") as f:
        raw_list = list(csv.DictReader(f, delimiter="\t"))
    for each_dict in raw_list:
        each_dict["blast hits"] = json.loads(each_dict["blast hits"])
    return raw_list

def analysis_combination_of_primers(input_data_as_dict: list) -> list:
    considered_pair_set = set()
    no_intersection_primer_pairs = []
    for each_pair in combinations(input_data_as_dict, 2):
        primer1, primer2 = each_pair
        primer1_id = primer1["primer id"]
        primer2_id = primer2["primer id"]
        if primer1_id == primer2_id or (primer1_id, primer2_id) in considered_pair_set:
            continue
        primer1_blast_hit = set(map(lambda x: ";".join(x), primer1["blast hits"]))
        primer2_blast_hit = set(map(lambda x: ";".join(x), primer2["blast hits"]))
        intersection = set(primer1_blast_hit).intersection(set(primer2_blast_hit))
        if any(list(map(lambda x: "human" in x, primer1_blast_hit | primer2_blast_hit))):
            continue
        # print_msg = f"{primer1['primer id']} and {primer2['primer id']} have {len(intersection)} BLAST hits in common\t{'/'.join(list(intersection))}"
        # print(print_msg)
        if len(intersection) == 0:
            primer1_blast_hit_count = len(primer1_blast_hit)
            primer2_blast_hit_count = len(primer2_blast_hit)
            total_count = primer1_blast_hit_count + primer2_blast_hit_count
            no_intersection_primer_pairs.append((total_count, primer1, primer2))

        considered_pair_set.add((primer1_id, primer2_id))
    return sorted(no_intersection_primer_pairs, key=lambda x: x[0])


def main():
    parser = argparse.ArgumentParser(description="survivor から2つ選び、BLAST Hit先が重複しないか検証するスクリプト")
    parser.add_argument("input", metavar="input", type=str, help="input file name")
    parser.add_argument("-o", metavar="output_prefix", type=str, default="final_result", help="output file name (default = final_result)")
    args = parser.parse_args()
    print(args, file=sys.stderr)
    input_data_as_dict = file_parser(args.input)
    print(f"{len(input_data_as_dict)} primer pairs to be considered")
    print(f"total number of combination is {len(list(combinations(input_data_as_dict, 2)))}")
    no_intersection_primer_pairs_sorted = analysis_combination_of_primers(input_data_as_dict)
    print(f"{len(no_intersection_primer_pairs_sorted)} primer pairs have no intersection")
    # print(json.dumps(no_intersection_primer_pairs_sorted, indent=2))


if __name__ == "__main__":
    main()