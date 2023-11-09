#! /usr/bin/env python3

import sys
import argparse
import pprint
import json
import csv
from Bio import SeqIO

pp = pprint.PrettyPrinter(indent=2)


class BlastResult:
    def __init__(
        self,
        qseqid,
        sseqid,
        sacc,
        qlen,
        qstart,
        qend,
        slen,
        sstart,
        send,
        qseq,
        sseq,
        evalue,
        length,
        staxid,
        staxids,
        ssciname,
        scomname,
    ):
        self.qseqid = qseqid
        self.sseqid = sseqid
        self.sacc = sacc
        self.qlen = int(qlen)
        self.qstart = int(qstart)
        self.qend = int(qend)
        self.slen = int(slen)
        self.sstart = int(sstart)
        self.send = int(send)
        self.qseq = qseq
        self.sseq = sseq
        self.evalue = float(evalue)
        self.length = int(length)
        self.staxid = staxid
        self.staxids = staxids
        self.ssciname = ssciname
        self.scomname = scomname

    def __str__(self):
        return "\t".join([f"{k}:{v}" for k, v in self.__dict__.items()])


overlook_taxon_ids = [
    4528,  # Oryza longistaminata
    4529,  # Oryza rufipogon
    4530,  # Oryza sativa
    39946,  # Oryza sativa Indica Group
    39947,  # Oryza sativa Japonica Group
    1080340,  # Oryza sativa Japonica Group x Oryza sativa Indica Group
    1736656,  # Oryza sativa tropical japonica subgroup
    1736657,  # Oryza sativa temperate japonica subgroup
    1736658,  # Oryza sativa aromatic subgroup
    1736659,  # Oryza sativa aus subgroup
    1771142,  # Oryza sativa indica subgroup
    2998809,  # Oryza sativa f. spontanea
    1050722,  # Oryza sativa Indica Group x Oryza sativa Japonica Group
    256318,  # metagenome
    408170,  # human gut metagenome
    408172,  # marine metagenome
    410658,  # soil metagenome
    410661,  # mouse gut metagenome
    412755,  # marine sediment metagenome
    496920,  # saltern metagenome
    527639,  # wastewater metagenome
    527640,  # microbial mat metagenome
    652676,  # hydrothermal vent metagenome
    717931,  # groundwater metagenome
    749906,  # gut metagenome
    749907,  # sediment metagenome
    797283,  # ant fungus garden metagenome
    904678,  # hypersaline lake metagenome
    1076179,  # bioreactor metagenome
    1911570,  # Mangrovicoccus ximenensis
    2509717,  # Ciceribacter ferrooxidans
    "N/A",
]


def main():
    parser = argparse.ArgumentParser(
        description="プライマー加工前のLR-tupleをBLAST検索した結果、どれだけ篩にかけられたかを確認するスクリプト"
    )
    parser.add_argument("fasta", metavar="fasta", type=str, help="primers.fa file name")
    parser.add_argument(
        "blast",
        metavar="blast",
        type=str,
        help="blast output file name. outfmt must be '6 qseqid sseqid sacc slen qstart qend sstart send qseq sseq evalue length staxid staxids ssciname scomname'",
    )
    args = parser.parse_args()
    print(args, file=sys.stderr)

    print(f"start reading {args.fasta}", file=sys.stderr)
    sequence_names = []
    with open(args.fasta, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence_names.append(record.id)
    print(
        f"finish reading {args.fasta}, {len(sequence_names)} records", file=sys.stderr
    )

    print(f"start reading {args.blast}", file=sys.stderr)
    blast_results = []
    with open(args.blast) as f:
        reader = csv.reader(f, delimiter="\t")
        for each_record in reader:
            try:
                each_record_Obj = BlastResult(*each_record)
            except TypeError:
                print(type(each_record), each_record, file=sys.stderr)
                continue  # ここでcontinueすると、不完全なBLAST hitを無視する結果になるのでは？

            taxon_id = -1
            try:
                taxon_id = int(each_record_Obj.staxid)
            except ValueError:
                taxon_id = each_record_Obj.staxid
            if taxon_id in overlook_taxon_ids:
                continue
            distance = -1
            # if each_record_Obj.qseqid.endswith("L") and each_record_Obj.qlen - each_record_Obj.qstart - each_record_Obj.length != 0:
            if each_record_Obj.qseqid.endswith("L"):
                distance = each_record_Obj.qlen - each_record_Obj.qend + 1
            elif each_record_Obj.qseqid.endswith("R"):
                distance = each_record_Obj.qstart
            else:
                print("Never reached", file=sys.stderr)
                sys.exit(1)
            # print(distance, args.offset, distance > args.offset, each_record_Obj, file = sys.stderr)
            if distance > 1:
                continue  # 救済
            # metagenomeという文字列がscomnameやscinameに入っていればcontinue
            if (
                each_record_Obj.scomname is not None
                and "metagenome" in each_record_Obj.scomname
            ):
                continue
            if (
                each_record_Obj.ssciname is not None
                and "metagenome" in each_record_Obj.ssciname
            ):
                continue
            blast_results.append(each_record_Obj)
    print(f"found {len(blast_results)} blast results", file=sys.stderr)

    lr_tuple_set = set(sequence_names)
    blasthit_results_set = set(blast_results)
    print(f"lr_tuple_set: {lr_tuple_set}")
    print(f"blasthit_results_set: {blasthit_results_set}")
    print(f"lr_tuple_set - blasthit_results_set: {lr_tuple_set - blasthit_results_set}")
