import argparse
from Bio import SeqIO
from collections import Counter
import subprocess
import os


def main():
    parser = argparse.ArgumentParser(description="This tool redesigns PCR primers by extracting and grouping the first and last 30 bases of sequences from a FASTA file. Within each group, it identifies frequently occurring k-mers in the central region of the sequences. Using these 30-base segments and the most common k-mer, it designs primers for PCR amplification.")
    parser.add_argument("fasta_file", type=str, help="FASTA File path")
    parser.add_argument(
        "primer3_config_template_path", type=str, help="primer3 config template path"
    )
    parser.add_argument("output_dir", type=str, help="output directory path")
    args = parser.parse_args()
    sequence_groups = group_sequences(args.fasta_file)
    print(f"Number of sequence groups: {len(sequence_groups)}")
    for group, sequences in sequence_groups.items():
        print(
            f"Group: {group}, sequences: {len(sequences)}, left: {sequences[0].seq[:30]}, right: {sequences[0].seq[-30:]}"
        )
        frequent_mers = find_frequent_mers(sequences, 50)
        l = sequences[0].seq[:30]
        m = frequent_mers.most_common(1)[0][0]
        r = sequences[0].seq[-30:]
        many_n = "N" * 20
        concatinated_sequence = f"{l}{many_n}{m}{many_n}{r}"
        # print(concatinated_sequence)
        config_file_name = os.path.join(args.output_dir, group + ".config")
        output_path = os.path.join(args.output_dir, group + ".primer3.out")
        run_primer3(
            concatinated_sequence,
            args.primer3_config_template_path,
            config_file_name,
            output_path,
        )


def group_sequences(fasta_file):
    groups = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        group_name = record.id.split("_")[0]
        if group_name not in groups:
            groups[group_name] = []
        groups[group_name].append(record)
    return groups


def find_frequent_mers(sequences, mer_size):
    mer_counts = Counter()
    for i, record in enumerate(sequences):
        print(f"processing {i + 1} in {len(sequences)}")
        sequence = str(record.seq)
        middle_sequence = sequence[30:-30]
        for i in range(len(middle_sequence) - mer_size + 1):
            mer = middle_sequence[i : i + mer_size]
            mer_counts[mer] += 1
    return mer_counts


def run_primer3(
    concatenated_sequence,
    primer3_template_path,
    primer3_config_path,
    primer3_output_path,
):
    if not os.path.exists("primer3_config"):
        os.makedirs("primer3_config")
    if not os.path.exists("primer3_results"):
        os.makedirs("primer3_results")

    with open(primer3_template_path, "r") as template_file:
        config_content = template_file.read()
    config_content += (
        f"SEQUENCE_ID=example\nSEQUENCE_TEMPLATE={concatenated_sequence}\n=\n"
    )
    with open("primer3_config/" + primer3_config_path, "w") as config_file:
        config_file.write(config_content)
    subprocess.run(
        [
            "primer3_core",
            "primer3_config/" + primer3_config_path,
            "--output",
            "primer3_results_with_amplicon/" + primer3_output_path,
        ]
    )


if __name__ == "__main__":
    main()
