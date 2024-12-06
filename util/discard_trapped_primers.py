#! /usr/bin/env python3

# from primer3_result_parser import primer3_result_parser, extract_primer_pairs
# from itertools import combinations
import sys
import argparse
import pprint
import json
import csv
from Bio import SeqIO
from enum import Enum
import itertools
import gzip

pp = pprint.PrettyPrinter(indent=2)


class Direction(Enum):
    RIGHT = "R"
    LEFT = "L"


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
        # 向きを調べる
        # Lプライマーの時
        if qseqid.endswith("L"):
            # qstart < qendのとき
            if qstart < qend:
                self.direction = Direction.RIGHT
            # qstart > qendのとき
            else:
                self.direction = Direction.LEFT
        # Rプライマーの時
        elif qseqid.endswith("R"):
            # qstart > qendのとき
            if qstart > qend:
                self.direction = Direction.RIGHT
            # qstart < qendのとき
            else:
                self.direction = Direction.LEFT

    def __str__(self):
        return "\t".join([f"{k}:{v}" for k, v in self.__dict__.items()])

    def subject_info(self):
        return ";".join(
            [
                f"{k}:{v}"
                for k, v in self.__dict__.items()
                if k
                not in [
                    "qseq",
                    "sseq",
                    "qlen",
                    "slen",
                    "qstart",
                    "qend",
                    "qseq",
                    "sseq",
                    "evalue",
                    "length",
                ]
            ]
        )


def grep_input_record_name(primer_pairs_dict):
    retset = set()
    for k, v in primer_pairs_dict.items():
        for idx, each_pair in enumerate(v["Primer3_output"]):
            retset.add(f">{k}_{idx}_L")
            retset.add(f">{k}_{idx}_M")
            retset.add(f">{k}_{idx}_R")
    return retset


def blast_hits_string(hit_1, hit_2):
    string = f"{hit_1.qseqid}\t{hit_2.qseqid}\t{hit_1.sacc == hit_2.sacc}\t{hit_1.direction}\t{hit_2.direction}\t{hit_1.sstart}\t{hit_2.sstart}"
    if hit_1.sstart < hit_2.sstart:
        string += f"\thit_1:{hit_1.sstart} {hit_1.direction.value}\thit_2:{hit_2.sstart} {hit_2.direction.value}"
    else:
        string += f"\thit_2:{hit_2.sstart} {hit_2.direction.value}\thit_1:{hit_1.sstart} {hit_1.direction.value}"
    return string


def main():
    parser = argparse.ArgumentParser(
        description="compare input fasta file and blast result"
    )
    parser.add_argument("fasta", metavar="fasta", type=str, help="primers.fa file name")
    parser.add_argument(
        "blast",
        metavar="blast",
        type=str,
        help="blast output file name. outfmt must be '6 qseqid sseqid sacc slen qstart qend sstart send qseq sseq evalue length staxid staxids ssciname scomname'",
    )
    parser.add_argument(
        "primer3", metavar="primer3", type=str, help="primers.json file name"
    )
    parser.add_argument(
        "-o",
        metavar="output_prefix",
        type=str,
        default="final_result",
        help="output file name (default = final_result)",
    )
    parser.add_argument("--distance", type = int, default = 20000, help = "Distance considered to be far enough away")
    args = parser.parse_args()
    print("\n".join([f"{key}...{value}" for key, value in vars(args).items()]), file=sys.stderr)

    filename = args.fasta
    fasta_ids = set()
    print(f"start reading {filename}", file=sys.stderr)
    with open(filename) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            fasta_ids.add(record.id)
    print(f"found {len(fasta_ids)} fasta records", file=sys.stderr)
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
    # overlook_taxon_ids = [4530, 39947, 1080340, 1050722, 1736656, 1736657, 1736658, 1736659, 1771142, 2998809, "N/A"]
    # 4529はOryza rufipogon, 野生の稲
    # 4528はOryza longistaminata
    print(f"start reading {args.blast}", file=sys.stderr)
    blast_results = []

    with gzip.open(args.blast, "rt") as f:
        reader = csv.reader(f, delimiter="\t")
        for each_record in reader:
            try:
                each_record_Obj = BlastResult(*each_record)
            except TypeError:
                print(type(each_record), each_record, file=sys.stderr)
                continue
            taxon_id = -1
            try:
                taxon_id = int(each_record_Obj.staxid)
            except ValueError:
                taxon_id = each_record_Obj.staxid
            if taxon_id in overlook_taxon_ids:
                continue

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
            # each_record_Obj.qstart > 3ならば、3'に2塩基のミスマッチがあると考える
            if each_record_Obj.qstart > 3:
                continue
            blast_results.append(each_record_Obj)
    print(f"found {len(blast_results)} blast results", file=sys.stderr)

    # プライマーペアの読み込み
    primer3_info = None
    with open(args.primer3, "r") as f:
        primer3_info = json.load(f)

    # blastの結果をプライマーペアに紐付ける
    primer_blasthit_dict = {k: [] for k in fasta_ids}
    for each_hit in blast_results:
        primer_blasthit_dict[each_hit.qseqid].append(each_hit)

    """
    primerのペアを周回する
    ペアのヒット先を調べる
    ヒット先のsseqidが同一か調べる→同一でなければサバイバー扱いにする
    同一のものがあった場合、さらにLL, LR, RRで挟めるかを調べる→挟めないのならサバイバー扱いにする
    """
    blast_trapped_seq_ids = set()
    blast_trapped_primer_ids = set()
    salvation_reason = {
        "total": 0,
        "hit to different sequence": 0,
        "hit to same sequence, same direction": 0,
        "hit to same sequence, opposite direction, no intersection": 0,
        "hit to same sequence, opposite direction, far enough away": 0,
    }
    for each_primer_id, info in primer3_info.items():
        for i, c in enumerate(info["Primer3_output"]):
            seqname_L = each_primer_id + "_" + str(i) + "_L"
            seqname_R = each_primer_id + "_" + str(i) + "_R"
            # ヒットをまとめる
            blast_hits = []  # list of BlastHit objects
            # Lのヒット先
            blast_hits.extend(primer_blasthit_dict.get(seqname_L, None))
            # Rのヒット先
            blast_hits.extend(primer_blasthit_dict.get(seqname_R, None))
            # print(blast_hits)
            for hit_1, hit_2 in itertools.combinations(blast_hits, 2):
                salvation_reason["total"] += 1
                # print(blast_hits_string(hit_1, hit_2), file=sys.stderr)
                if hit_1.sacc != hit_2.sacc:
                    salvation_reason["hit to different sequence"] += 1
                    continue
                if hit_1.direction == hit_2.direction:
                    salvation_reason["hit to same sequence, same direction"] += 1
                    continue
                if hit_1.sstart < hit_2.sstart and hit_1.direction == Direction.LEFT:
                    salvation_reason[
                        "hit to same sequence, opposite direction, no intersection"
                    ] += 1
                    continue
                if hit_1.sstart > hit_2.sstart and hit_1.direction == Direction.RIGHT:
                    salvation_reason[
                        "hit to same sequence, opposite direction, no intersection"
                    ] += 1
                    continue
                if abs(hit_1.sstart - hit_2.sstart) > args.distance:
                    salvation_reason[
                        "hit to same sequence, opposite direction, far enough away"
                    ] += 1
                    continue
                blast_trapped_seq_ids.add((hit_1, hit_2))
                blast_trapped_primer_ids.add(each_primer_id)

    finalist = list(
        filter(lambda i: i[0] not in blast_trapped_primer_ids, primer3_info.items())
    )

    report_file = open(args.o + ".report", mode="w")
    finalist_tsv_file = open(args.o + ".finalist.tsv", mode="w")
    finalist_namelist_file = open(args.o + ".finalist_name.txt", mode="w")
    cross_reactive_file = open(args.o + ".cross_reactive_species.txt", mode="w")
    print(f"{args}", file=report_file)
    print(
        f"fasta_ids in {args.fasta}...{len(fasta_ids)}",
        file=report_file,
    )
    print(
        f"primer pairs {args.fasta}...{len(fasta_ids)/2}",
        file=report_file,
    )
    print(
        f"blast hits in {args.blast}...{len(blast_results)}",
        file=report_file,
    )
    print(
        f"number of lr-tuple...{len(primer3_info)}",
        file=report_file,
    )
    if len(primer3_info) != 0:
        average_of_primers_from_a_lr_tuple = round(
            len(fasta_ids) / len(primer3_info) / 2, 2
        )
    else:
        average_of_primers_from_a_lr_tuple = 0
    print(
        f"average number of primers from a lr-tuple...{average_of_primers_from_a_lr_tuple}",
        file=report_file,
    )
    tmpstr = "\n".join([f"{k}:{v}" for k, v in salvation_reason.items()])
    print(
        f"Breakdown of reasons for not treating as hit...\n{tmpstr}",
        file=report_file,
    )
    print(
        f"trapped primers...{len(blast_trapped_primer_ids)}",
        file=report_file,
    )
    print(
        f"number of primer candidates...{len(finalist)}",
        file=report_file,
    )

    print(
        "\t".join(
            [
                "primer id",
                "left primer",
                "right primer",
                "primer left Tm",
                "primer right Tm",
                "primer pair product Tm",
            ]
        ),
        file=finalist_tsv_file,
    )
    cross_reactive_species_list = ["\t".join([x[0].staxid, x[0].staxids, x[0].ssciname, x[0].scomname]) for x in blast_trapped_seq_ids]
    print("\n".join(cross_reactive_species_list), file = cross_reactive_file)

    for each_finalist in finalist:
        id = each_finalist[0]
        print(id, file=finalist_namelist_file)
        primers = each_finalist[1]["Primer3_output"]
        for each_primer in primers:
            finalist_info = each_primer
            """ 
            left_blast_hit = sequence_blasthit_dict[
                finalist_info["PRIMER_LEFT_SEQUENCE"]
            ]
            right_blast_hit = sequence_blasthit_dict[
                finalist_info["PRIMER_RIGHT_SEQUENCE"]
            ]
            blasthit_cnt = len(left_blast_hit) + len(right_blast_hit)
            """
            print(
                "\t".join(
                    [
                        str(x)
                        for x in [
                            id,
                            finalist_info["PRIMER_LEFT_SEQUENCE"],
                            finalist_info["PRIMER_RIGHT_SEQUENCE"],
                            finalist_info["PRIMER_LEFT_TM"],
                            finalist_info["PRIMER_RIGHT_TM"],
                            finalist_info["PRIMER_PAIR_PRODUCT_TM"],
                        ]
                    ]
                ),
                file=finalist_tsv_file,
            )
    print("Done", file=sys.stderr)


if __name__ == "__main__":
    main()
