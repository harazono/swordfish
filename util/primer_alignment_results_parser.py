import sys

def parse_blast_output(blast_output_file):
    alignments = {}
    with open(blast_output_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            # フィールドの取得
            qseqid, sseqid, sacc, qlen, qstart, qend, slen, sstart, send, qseq, sseq, length = parts
            if sseqid not in alignments:
                alignments[sseqid] = []
            alignments[sseqid].append({
                'primer_probe': qseqid,
                'sstart': sstart,
                'send': send,
                'alignment_length': length,
                'qseq': qseq,
                'sseq': sseq
            })

    return alignments

def print_alignments(alignments):
    for sseqid, aligns in alignments.items():
        print(f"NGS Read ID: {sseqid}")
        for align in aligns:
            print(f"\tPrimer/Probe: {align['primer_probe']}, Start: {align['sstart']}, End: {align['send']}, Alignment Length: {align['alignment_length']}")
            print(f"\tQuery Sequence: {align['qseq']}")
            print(f"\tSubject Sequence: {align['sseq']}\n")


if __name__ == "__main__":
    blast_output_file = sys.argv[1]  # BLAST結果ファイルへのパス
    alignments = parse_blast_output(blast_output_file)
    print_alignments(alignments)
