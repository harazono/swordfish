#! /usr/bin/env python3

import sys
import argparse
import pprint
import json
import csv
from Bio import SeqIO
pp = pprint.PrettyPrinter(indent = 2)
from primer3_result_parser import primer3_result_parser, extract_primer_pairs


class BlastResult():
	def __init__(self, qseqid, sseqid, sacc, qlen, qstart, qend, slen, sstart, send, qseq, sseq, evalue, length, staxid, staxids, ssciname, scomname):
		self.qseqid   = qseqid
		self.sseqid   = sseqid
		self.sacc     = sacc
		self.qlen     = qlen
		self.qstart   = qstart
		self.qend     = qend
		self.slen     = slen
		self.sstart   = sstart
		self.send     = send
		self.qseq     = qseq
		self.sseq     = sseq
		self.evalue   = evalue
		self.length   = length
		self.staxid   = staxid
		self.staxids  = staxids
		self.ssciname = ssciname
		self.scomname = scomname


def grep_input_record_name(primer_pairs_dict):
	retset = set()
	for k, v in primer_pairs_dict.items():
		for idx, each_pair in enumerate(v["Primer3_output"]):
			retset.add(f">{k}_{idx}_L")
			retset.add(f">{k}_{idx}_M")
			retset.add(f">{k}_{idx}_R")
	return retset


def main():
	parser = argparse.ArgumentParser(description = "compare input fasta file and blast result")
	parser.add_argument("fasta",   metavar = "fasta",       type = str, help = "fasta file name")
	parser.add_argument("blast",   metavar = "blast",       type = str, help = "blast output file name. outfmt must be '6 qseqid sseqid sacc slen qstart qend sstart send qseq sseq evalue length staxid staxids ssciname scomname'")
	parser.add_argument("primer3", metavar = "primer3",     type = str, help = "primer3 output file (json)")
	parser.add_argument("-o",      metavar = "output_file", type = str, default = "sys.stdout", help = "output file name (default = sys.stdout)")
	parser.add_argument("--fasta", action='store_true',     help = "output as fasta")
	args = parser.parse_args()
	filename = args.fasta
	fasta_ids = set()
	output_file = open(args.o, "w") if args.o is not None else sys.stdout

	with open(filename) as handle:
		for record in SeqIO.parse(handle, "fasta"):
			fasta_ids.add(record.id)

	blast_results = []
	with open(args.blast) as f:
		reader = csv.reader(f, delimiter='\t')
		for each_record in reader:
			each_record_Obj = BlastResult(*each_record)
			#救済対象か判断。対象でなければappend。
			if int(each_record_Obj.staxid) not in [4530, 39947, 1080340, 1050722, 1736656, 1736657, 1736658, 1736659, 1771142, 2998809]:
				if "metagenome" not in each_record_Obj.scomname:
					if "metagenome" not in each_record_Obj.ssciname:
						blast_results.append(each_record_Obj)

	primer3_info = None
	with open (args.primer3, "r") as f:
		primer3_info = json.load(f)


	primer_blasthit_dict = {k: set() for k in fasta_ids}
	for each_hit in blast_results:
		primer_blasthit_dict[each_hit.qseqid].add(f"{each_hit.scomname}({each_hit.staxid}: {each_hit.ssciname})")
	print(primer_blasthit_dict, file = sys.stderr)
	#for i in primer_blasthit_dict:
	#	print(len(primer_blasthit_dict[i]), primer_blasthit_dict[i])


	blast_trapped_seq_ids = set()
	for each_hit in blast_results:
		blast_trapped_seq_ids.add(each_hit.qseqid)



	survivor = fasta_ids - blast_trapped_seq_ids 
	survivor_pair = set()
	for each_primer in survivor:
		pair = ""
		if each_primer[-1] == "L":
			pair = each_primer[0:-1] + "R"
		else:
			pair = each_primer[0:-1] + "L"
		if pair in survivor:
			survivor_pair.add(each_primer[0:-2])
		else:
			pass
			#print(f"{each_primer} does not have a partner")
	print(f"total count of input sequence                  : {len(fasta_ids)}"            , file = output_file)
	print(f"total count of blast hits                      : {len(blast_results)}"        , file = output_file)
	print(f"cardinarity of input which was trapped by blast: {len(blast_trapped_seq_ids)}", file = output_file)
	print(f"survivor                                       : {len(survivor)}"             , file = output_file)
	print(f"survivor pair                                  : {len(survivor_pair)}"        , file = output_file)

	print("\t".join(["primer id", "left primer", "right primer", "primer left Tm", "primer right Tm", "primer pair product Tm", "survived side", "trapped side", "blast hits"]))
	for each_survivor in survivor:
		survivor_info_raw = primer3_info[each_survivor.split("_")[0]]
		survivor_index    = int(each_survivor.split("_")[1])
		primer_side       = each_survivor.split("_")[2]
		survivor_info     = survivor_info_raw["Primer3_output"][survivor_index]
		primer_pair       = each_survivor[:-1] + "L" if primer_side == "R" else each_survivor[:-1] + "R"
		print(each_survivor, file = sys.stderr)
		blast_hits = primer_blasthit_dict[primer_pair]
		print("\t".join([str(x) for x in [each_survivor.split("_")[0], survivor_info["PRIMER_LEFT_SEQUENCE"], survivor_info["PRIMER_RIGHT_SEQUENCE"], survivor_info["PRIMER_LEFT_TM"], survivor_info["PRIMER_RIGHT_TM"], survivor_info["PRIMER_PAIR_PRODUCT_TM"], primer_side, "L" if primer_side == "R" else "R", blast_hits]]))
		#print(f"{primer_side} survived\n---\n")


	primer_pairs_dict = None
	with open(args.primer3, "r") as f:
		primer_pairs_dict = json.load(f)
	#pp.pprint(primer_pairs_dict)


	for each_survivor_pair in survivor_pair:
		seqid = each_survivor_pair[0:-2]
		index = each_survivor_pair[-1]
		print(seqid, index, file = output_file)
		primer_info = primer_pairs_dict[seqid]["Primer3_output"][int(index)]
		pp.pprint(primer_info, file = output_file)


if __name__ == '__main__':
	main()


"""
c863b6bb4258705996f0adf5ad3a0d89 4
{ 'PRIMER_INTERNAL': '59,22',
  'PRIMER_INTERNAL_GC_PERCENT': '59.091',
  'PRIMER_INTERNAL_HAIRPIN_TH': '41.51',
  'PRIMER_INTERNAL_PENALTY': '2.441787',
  'PRIMER_INTERNAL_SELF_ANY_TH': '14.73',
  'PRIMER_INTERNAL_SELF_END_TH': '10.07',
  'PRIMER_INTERNAL_SEQUENCE': 'CGACTAACCGCGCCGTTAAGGT',
  'PRIMER_INTERNAL_TM': '59.558',
  'PRIMER_LEFT': '3,18',
  'PRIMER_LEFT_END_STABILITY': '4.0100',
  'PRIMER_LEFT_GC_PERCENT': '55.556',
  'PRIMER_LEFT_HAIRPIN_TH': '44.95',
  'PRIMER_LEFT_PENALTY': '2.261888',
  'PRIMER_LEFT_SELF_ANY_TH': '7.38',
  'PRIMER_LEFT_SELF_END_TH': '0.00',
  'PRIMER_LEFT_SEQUENCE': 'ACGATGTCGGTGTCAAGC',
  'PRIMER_LEFT_TM': '57.738',
  'PRIMER_PAIR_COMPL_ANY_TH': '0.00',
  'PRIMER_PAIR_COMPL_END_TH': '1.26',
  'PRIMER_PAIR_PENALTY': '4.295782',
  'PRIMER_PAIR_PRODUCT_SIZE': '137',
  'PRIMER_PAIR_PRODUCT_TM': '74.6',
  'PRIMER_RIGHT': '139,18',
  'PRIMER_RIGHT_END_STABILITY': '5.2800',
  'PRIMER_RIGHT_GC_PERCENT': '61.111',
  'PRIMER_RIGHT_HAIRPIN_TH': '0.00',
  'PRIMER_RIGHT_PENALTY': '2.033894',
  'PRIMER_RIGHT_SELF_ANY_TH': '0.00',
  'PRIMER_RIGHT_SELF_END_TH': '0.00',
  'PRIMER_RIGHT_SEQUENCE': 'GCTCGATTCCATGACCGG',
  'PRIMER_RIGHT_TM': '57.966'}"""