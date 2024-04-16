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
            sseqid = parts[1]
            record = dict(
                qseqid=parts[0],
                sseqid=parts[1],
                sacc=parts[2],
                qlen=parts[3],
                qstart=int(parts[4]),
                qend=int(parts[5]),
                slen=int(parts[6]),
                sstart=int(parts[7]),
                send=int(parts[8]),
                qseq=parts[9],
                sseq=parts[10],
                length=int(parts[11]),
                sstrand=parts[12],
            )
            if sseqid not in hits: 
                hits[sseqid] = []
            hits[sseqid].append(record)
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
            if hit_1.sstrand == hit_2.sstrand:
                continue
            if hit_1.send < hit_2.sstart:
                amplicon_name = f"{record.id}|{hit_1.qseqid}{hit_1.sstrand}|{hit_2.qseqid}{hit_2.sstrand}"
                assert hit_1.sstart < hit_2.send, "assertion failed"
                amplicon_sequence = sequence[hit_1.sstart:hit_2.send]
                amplicons.append((amplicon_name, amplicon_sequence))
            if hit_2.send < hit_1.sstart:
                amplicon_name = f"{record.id}|{hit_1.qseqid}{hit_1.sstrand}|{hit_2.qseqid}{hit_2.sstrand}"
                assert hit_2.sstart < hit_1.send, "assertion failed"
                amplicon_sequence = sequence[hit_2.sstart:hit_1.send]
                amplicons.append((amplicon_name, amplicon_sequence))
    return amplicons

def write_hits_to_file(hits, output_file):
    """
    抽出されたヒット情報をファイルに書き出します。
    """
    with open(output_file, 'w') as f:
        for sseqid, hit_list in hits.items():
            f.write(f"{sseqid}\n")
            for hit in hit_list:
                f.write("\t".join([f"{str(k)}:{str(v)}" for k,v in hit.items()]) + "\n")

def main(args):
    primer_hits = parse_blast_output(args.blast_output)
    write_hits_to_file(primer_hits, args.output_file + ".hits")
    amplicons = extract_amplicons(args.fasta_file, primer_hits, args.max_length)
    with open(args.output_file + ".fa", 'w') as f:
        for amplicon_info in amplicons:
            if len(amplicon_info[1]) != 0:
                f.write(f">{amplicon_info[0]}\n{amplicon_info[1]}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract amplicons based on primer hits from a BLAST output.")
    parser.add_argument("-f", "--fasta_file", required=True, help="Path to the input FASTA file.")
    parser.add_argument("-b", "--blast_output", required=True, help="Path to the BLAST output file in outfmt 6 format.")
    parser.add_argument("-m", "--max_length", type=int, default=2000, help="Maximum amplicon length (default is 2000).")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output file to write hits.")
    args = parser.parse_args()
    main(args)
