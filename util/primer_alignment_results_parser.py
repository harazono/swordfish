import sys
from Bio import SeqIO

def extract_amplicons(fasta_file, max_length=2000):
    amplicons = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq).upper()
        # 2種類のプライマーの位置を見つける (ここでは仮のプライマーを使用)
        primer1_pos = sequence.find("PRIMER1_SEQUENCE")
        primer2_pos = sequence.find("PRIMER2_SEQUENCE")
        if primer1_pos != -1 and primer2_pos != -1 and 0 < primer2_pos - primer1_pos <= max_length:
            # プライマーに挟まれた領域を抽出
            amplicon = sequence[primer1_pos:primer2_pos+len("PRIMER2_SEQUENCE")]
            amplicons.append((record.id, amplicon))

    return amplicons

def main(fasta_file, max_length):
    amplicons = extract_amplicons(fasta_file, max_length)
    for amplicon_id, amplicon_seq in amplicons:
        print(f">{amplicon_id}\n{amplicon_seq}")


if __name__ == "__main__":
    fasta_file = sys.argv[1]  # FASTAファイルへのパス
    max_length = int(sys.argv[2]) if len(sys.argv) > 2 else 2000  # アンプリコンの最大長
    main(fasta_file, max_length)
