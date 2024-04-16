import argparse
from Bio import SeqIO
import sys
import itertools
def parse_blast_output(blast_output_file):
    """
    BLASTのoutfmt 6形式の出力から、各プライマーのヒット情報を抽出します。
    """
    hits = {}
    with open(blast_output_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            query_id = parts[0]
            sseqid   = parts[1]
            sstart   = int(parts[7])
            send     = int(parts[8])
            if sseqid not in hits: 
                hits[sseqid] = []
            hits[sseqid].append((query_id, sstart, send))
    return hits

def extract_amplicons(fasta_file, primer_hits, max_length):
    """
    指定されたプライマーのヒット情報と最大アンプリコン長を基に、FASTAファイルからアンプリコンを抽出します。
    """
    reads = SeqIO.parse(fasta_file, "fasta")
    cnt = 1
    amplicons = []
    for record in reads:
        print(f"\r{cnt}", end="", file=sys.stderr)
        cnt += 1
        current_record_id = record.id
        if current_record_id not in primer_hits.keys():
            continue
        hits = primer_hits[current_record_id]
        if len(hits) == 0:
            continue
        for v in itertools.combinations(hits, 2):
            hit_1, hit_2 = v
            sequence = str(record.seq).upper()
            # hit_1が右向き（end - start > 0）で、hit_2が左向きの場合
            if hit_1[2] - hit_1[1] > 0 and hit_2[1] - hit_2[2] > 0 and hit_2[2] - hit_1[1] > max_length:
                # read.idとhit_1[1]とhit_2[2]の値を配列名とするfastaを生成してampliconsに追加。配列はsequenceの部分配列とする。
                amplicon_name = f"{record.id}|{hit_1[1]}|{hit_2[2]}"
                amplicon_sequence = sequence[hit_1[1]:hit_2[2]]
                amplicons.append((amplicon_name, hit_1[1], hit_2[2], amplicon_sequence))
            # hit_1が左向き（end - start < 0）で、hit_2が右向きの場合
            elif hit_1[2] - hit_1[1] < 0 and hit_2[1] - hit_2[2] > 0 and hit_1[2] - hit_2[1] > max_length:
                amplicon_name = f"{record.id}|{hit_2[1]}|{hit_1[2]}"
                amplicon_sequence = sequence[hit_2[1]:hit_1[2]]
                amplicons.append((amplicon_name, hit_2[1], hit_1[2], amplicon_sequence))
    return amplicons

def write_hits_to_file(hits, output_file):
    """
    抽出されたヒット情報をファイルに書き出します。
    """
    with open(output_file, 'w') as f:
        for sseqid, hit_list in hits.items():
            f.write(f"{sseqid}\n")
            for hit in hit_list:
                qid, start, end = hit
                f.write(f"\t{qid}\t{start}\t{end}\n")

def main(args):
    primer_hits = parse_blast_output(args.blast_output)
    write_hits_to_file(primer_hits, args.output_file + ".hits")
    amplicons = extract_amplicons(args.fasta_file, primer_hits, args.max_length)
    with open(args.output_file + ".fa", 'w') as f:
        for amplicon_info in amplicons:
            if len(amplicon_info[3]) != 0:
                f.write(f">{amplicon_info[0]}\n{amplicon_info[3]}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract amplicons based on primer hits from a BLAST output.")
    parser.add_argument("-f", "--fasta_file", required=True, help="Path to the input FASTA file.")
    parser.add_argument("-b", "--blast_output", required=True, help="Path to the BLAST output file in outfmt 6 format.")
    parser.add_argument("-m", "--max_length", type=int, default=2000, help="Maximum amplicon length (default is 2000).")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output file to write hits.")
    args = parser.parse_args()
    main(args)
