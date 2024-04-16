import argparse
from Bio import SeqIO
import sys
def parse_blast_output(blast_output_file):
    """
    BLASTのoutfmt 6形式の出力から、各プライマーのヒット情報を抽出します。
    """
    hits = {}
    with open(blast_output_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            query_id = parts[0]
            sstart = int(parts[7])
            send = int(parts[8])
            if query_id not in hits:
                hits[query_id] = []
            hits[query_id].append((min(sstart, send), max(sstart, send)))
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
        sequence = str(record.seq).upper()
        for primer_id, hits in primer_hits.items():
            for hit1 in hits:
                for hit2 in hits:
                    if hit1 == hit2:
                        continue
                    start, end = sorted([hit1[0], hit2[1]])
                    if 0 < end - start <= max_length:
                        amplicon = sequence[start-1:end]
                        amplicons.append([record.id, primer_id, str(start), str(end), amplicon])
    return amplicons

def write_hits_to_file(hits, output_file):
    """
    抽出されたヒット情報をファイルに書き出します。
    """
    with open(output_file, 'w') as f:
        for query_id, hit_list in hits.items():
            for hit in hit_list:
                start, end = hit
                f.write(f"Primer/Probe: {query_id}, Start: {start}, End: {end}\n")

def main(args):
    primer_hits = parse_blast_output(args.blast_output)
    write_hits_to_file(primer_hits, args.output_file + ".hits")
    amplicons = extract_amplicons(args.fasta_file, primer_hits, args.max_length)
    for amplicon_info in amplicons:
        print(f">{amplicon_info[0]}|{amplicon_info[1]}|{amplicon_info[2]}-{amplicon_info[3]}\n{amplicon_info[4]}", file=args.output_file + ".fa")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract amplicons based on primer hits from a BLAST output.")
    parser.add_argument("-f", "--fasta_file", required=True, help="Path to the input FASTA file.")
    parser.add_argument("-b", "--blast_output", required=True, help="Path to the BLAST output file in outfmt 6 format.")
    parser.add_argument("-m", "--max_length", type=int, default=2000, help="Maximum amplicon length (default is 2000).")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output file to write hits.")
    args = parser.parse_args()
    main(args)
