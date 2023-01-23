#! /usr/bin/env python3

import argparse
import re
import pprint
import json
import sys
pp = pprint.PrettyPrinter(indent = 2)

"""
SEQUENCE_ID=c863b6bb425c1665bc2b29f5ad3a0d89
SEQUENCE_TEMPLATE=TAGACGATGTCGGTGTCAAGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTAACCGCGCCGTTAAGGTAGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTTCCGGTCATGGAATCGAGC
SEQUENCE_TARGET=28,86
PRIMER_TASK=pick_pcr_primers_and_hyb_probe
PRIMER_OPT_SIZE=18
PRIMER_MIN_SIZE=15
PRIMER_MAX_SIZE=21
PRIMER_PRODUCT_SIZE_RANGE=201-300 101-200
P3_FILE_FLAG=0
PRIMER_EXPLAIN_FLAG=1
PRIMER_LEFT_EXPLAIN=considered 35, too many Ns 7, low tm 22, ok 6
PRIMER_RIGHT_EXPLAIN=considered 33, too many Ns 5, low tm 23, high hairpin stability 1, ok 4
PRIMER_INTERNAL_EXPLAIN=considered 154, too many Ns 119, low tm 33, ok 2
PRIMER_PAIR_EXPLAIN=considered 36, unacceptable product size 30, ok 6
PRIMER_LEFT_NUM_RETURNED=5
PRIMER_RIGHT_NUM_RETURNED=5
PRIMER_INTERNAL_NUM_RETURNED=5
PRIMER_PAIR_NUM_RETURNED=5
"""
primer_pair_identifier = re.compile(r"_\d_")



class PrimerPairReport():
	def __init__(self, SEQUENCE_ID, SEQUENCE_TEMPLATE, SEQUENCE_TARGET, PRIMER_TASK, PRIMER_OPT_SIZE, PRIMER_MIN_SIZE, PRIMER_MAX_SIZE, PRIMER_PRODUCT_SIZE_RANGE, P3_FILE_FLAG, PRIMER_EXPLAIN_FLAG, PRIMER_LEFT_EXPLAIN, PRIMER_RIGHT_EXPLAIN, PRIMER_INTERNAL_EXPLAIN, PRIMER_PAIR_EXPLAIN, PRIMER_LEFT_NUM_RETURNED, PRIMER_RIGHT_NUM_RETURNED, PRIMER_INTERNAL_NUM_RETURNED, PRIMER_PAIR_NUM_RETURNED, PrimerPairs):
		self.SEQUENCE_ID                  = SEQUENCE_ID
		self.SEQUENCE_TEMPLATE            = SEQUENCE_TEMPLATE
		self.SEQUENCE_TARGET              = SEQUENCE_TARGET
		self.PRIMER_TASK                  = PRIMER_TASK
		self.PRIMER_OPT_SIZE              = PRIMER_OPT_SIZE
		self.PRIMER_MIN_SIZE              = PRIMER_MIN_SIZE
		self.PRIMER_MAX_SIZE              = PRIMER_MAX_SIZE
		self.PRIMER_PRODUCT_SIZE_RANGE    = PRIMER_PRODUCT_SIZE_RANGE
		self.P3_FILE_FLAG                 = P3_FILE_FLAG
		self.PRIMER_EXPLAIN_FLAG          = PRIMER_EXPLAIN_FLAG
		self.PRIMER_LEFT_EXPLAIN          = PRIMER_LEFT_EXPLAIN
		self.PRIMER_RIGHT_EXPLAIN         = PRIMER_RIGHT_EXPLAIN
		self.PRIMER_INTERNAL_EXPLAIN      = PRIMER_INTERNAL_EXPLAIN
		self.PRIMER_PAIR_EXPLAIN          = PRIMER_PAIR_EXPLAIN
		self.PRIMER_LEFT_NUM_RETURNED     = PRIMER_LEFT_NUM_RETURNED
		self.PRIMER_RIGHT_NUM_RETURNED    = PRIMER_RIGHT_NUM_RETURNED
		self.PRIMER_INTERNAL_NUM_RETURNED = PRIMER_INTERNAL_NUM_RETURNED
		self.PRIMER_PAIR_NUM_RETURNED     = PRIMER_PAIR_NUM_RETURNED
		self.Primer_Pairs                 = PrimerPairs
		

class PrimerPair():
	def __init__(self, PRIMER_PAIR_PENALTY, PRIMER_LEFT_PENALTY, PRIMER_RIGHT_PENALTY, PRIMER_INTERNAL_PENALTY, PRIMER_LEFT_SEQUENCE, PRIMER_RIGHT_SEQUENCE, PRIMER_INTERNAL_SEQUENCE, PRIMER_LEFT, PRIMER_RIGHT, PRIMER_INTERNAL, PRIMER_LEFT_TM, PRIMER_RIGHT_TM, PRIMER_INTERNAL_TM, PRIMER_LEFT_GC_PERCENT, PRIMER_RIGHT_GC_PERCENT, PRIMER_INTERNAL_GC_PERCENT, PRIMER_INTERNAL_SELF_ANY_TH, PRIMER_LEFT_SELF_ANY_TH, PRIMER_RIGHT_SELF_ANY_TH, PRIMER_INTERNAL_SELF_END_TH, PRIMER_LEFT_SELF_END_TH, PRIMER_RIGHT_SELF_END_TH, PRIMER_LEFT_HAIRPIN_TH, PRIMER_RIGHT_HAIRPIN_TH, PRIMER_INTERNAL_HAIRPIN_TH, PRIMER_LEFT_END_STABILITY, PRIMER_RIGHT_END_STABILITY, PRIMER_PAIR_COMPL_ANY_TH, PRIMER_PAIR_COMPL_END_TH, PRIMER_PAIR_PRODUCT_SIZE, PRIMER_PAIR_PRODUCT_TM):
		self.PRIMER_PAIR_PENALTY         = PRIMER_PAIR_PENALTY
		self.PRIMER_LEFT_PENALTY         = PRIMER_LEFT_PENALTY
		self.PRIMER_RIGHT_PENALTY        = PRIMER_RIGHT_PENALTY
		self.PRIMER_INTERNAL_PENALTY     = PRIMER_INTERNAL_PENALTY
		self.PRIMER_LEFT_SEQUENCE        = PRIMER_LEFT_SEQUENCE
		self.PRIMER_RIGHT_SEQUENCE       = PRIMER_RIGHT_SEQUENCE
		self.PRIMER_INTERNAL_SEQUENCE    = PRIMER_INTERNAL_SEQUENCE
		self.PRIMER_LEFT                 = PRIMER_LEFT
		self.PRIMER_RIGHT                = PRIMER_RIGHT
		self.PRIMER_INTERNAL             = PRIMER_INTERNAL
		self.PRIMER_LEFT_TM              = PRIMER_LEFT_TM
		self.PRIMER_RIGHT_TM             = PRIMER_RIGHT_TM
		self.PRIMER_INTERNAL_TM          = PRIMER_INTERNAL_TM
		self.PRIMER_LEFT_GC_PERCENT      = PRIMER_LEFT_GC_PERCENT
		self.PRIMER_RIGHT_GC_PERCENT     = PRIMER_RIGHT_GC_PERCENT
		self.PRIMER_INTERNAL_GC_PERCENT  = PRIMER_INTERNAL_GC_PERCENT
		self.PRIMER_INTERNAL_SELF_ANY_TH = PRIMER_INTERNAL_SELF_ANY_TH
		self.PRIMER_LEFT_SELF_ANY_TH     = PRIMER_LEFT_SELF_ANY_TH
		self.PRIMER_RIGHT_SELF_ANY_TH    = PRIMER_RIGHT_SELF_ANY_TH
		self.PRIMER_INTERNAL_SELF_END_TH = PRIMER_INTERNAL_SELF_END_TH
		self.PRIMER_LEFT_SELF_END_TH     = PRIMER_LEFT_SELF_END_TH
		self.PRIMER_RIGHT_SELF_END_TH    = PRIMER_RIGHT_SELF_END_TH
		self.PRIMER_LEFT_HAIRPIN_TH      = PRIMER_LEFT_HAIRPIN_TH
		self.PRIMER_RIGHT_HAIRPIN_TH     = PRIMER_RIGHT_HAIRPIN_TH
		self.PRIMER_INTERNAL_HAIRPIN_TH  = PRIMER_INTERNAL_HAIRPIN_TH
		self.PRIMER_LEFT_END_STABILITY   = PRIMER_LEFT_END_STABILITY
		self.PRIMER_RIGHT_END_STABILITY  = PRIMER_RIGHT_END_STABILITY
		self.PRIMER_PAIR_COMPL_ANY_TH    = PRIMER_PAIR_COMPL_ANY_TH
		self.PRIMER_PAIR_COMPL_END_TH    = PRIMER_PAIR_COMPL_END_TH
		self.PRIMER_PAIR_PRODUCT_SIZE    = PRIMER_PAIR_PRODUCT_SIZE
		self.PRIMER_PAIR_PRODUCT_TM      = PRIMER_PAIR_PRODUCT_TM


def fmt4fasta(primer_pairs_dict, trim_flag = False):
	retarray = []
	start_index = -19 if trim_flag else 0
	for k, v in primer_pairs_dict.items():
		for idx, each_pair in enumerate(v["Primer3_output"]):
			tmpstr = ""
			if each_pair['PRIMER_LEFT_SEQUENCE'] is not None:
				tmpstr += f">{k}_{idx}_L\n{each_pair['PRIMER_LEFT_SEQUENCE'][start_index:]}\n"
			if each_pair['PRIMER_INTERNAL_SEQUENCE'] is not None:
				tmpstr += f">{k}_{idx}_M\n{each_pair['PRIMER_INTERNAL_SEQUENCE'][start_index:]}\n"
			if each_pair['PRIMER_RIGHT_SEQUENCE'] is not None:
				tmpstr += f">{k}_{idx}_R\n{each_pair['PRIMER_RIGHT_SEQUENCE'][start_index:]}\n"
			retarray.append(tmpstr)
	return "\n".join(retarray)



def tempfmt(primer_pairs_dict):
	retarray = []
	start_index = 0
	for k, v in primer_pairs_dict.items():
		for idx, each_pair in enumerate(v["Primer3_output"]):
			tmpstr = ""
			if each_pair['PRIMER_LEFT_SEQUENCE'] is not None:
				tmpstr += f"{each_pair['PRIMER_LEFT_SEQUENCE'][start_index:]}\t"
			if each_pair['PRIMER_INTERNAL_SEQUENCE'] is not None:
				tmpstr += f"{each_pair['PRIMER_INTERNAL_SEQUENCE'][start_index:]}\t"
			if each_pair['PRIMER_RIGHT_SEQUENCE'] is not None:
				tmpstr += f"{each_pair['PRIMER_RIGHT_SEQUENCE'][start_index:]}"
			retarray.append(tmpstr)
	return "\n".join(retarray)



def primer3_result_parser(filename):
	retarray = []
	current_SEQUENCE_ID                  = ""
	current_SEQUENCE_TEMPLATE            = ""
	current_SEQUENCE_TARGET              = ""
	current_PRIMER_TASK                  = ""
	current_PRIMER_OPT_SIZE              = ""
	current_PRIMER_MIN_SIZE              = ""
	current_PRIMER_MAX_SIZE              = ""
	current_PRIMER_PRODUCT_SIZE_RANGE    = ""
	current_P3_FILE_FLAG                 = ""
	current_PRIMER_EXPLAIN_FLAG          = ""
	current_PRIMER_LEFT_EXPLAIN          = ""
	current_PRIMER_RIGHT_EXPLAIN         = ""
	current_PRIMER_INTERNAL_EXPLAIN      = ""
	current_PRIMER_PAIR_EXPLAIN          = ""
	current_PRIMER_LEFT_NUM_RETURNED     = ""
	current_PRIMER_RIGHT_NUM_RETURNED    = ""
	current_PRIMER_INTERNAL_NUM_RETURNED = ""
	current_PRIMER_PAIR_NUM_RETURNED     = ""
	current_PrimerPairs                  = {}
	current_match = None
	with open(filename, 'r') as f:
		for l in f:
			if l == "\n":
				continue
			if l.startswith("="):
				tmp_PrimerPair = PrimerPairReport(
					current_SEQUENCE_ID,
					current_SEQUENCE_TEMPLATE,
					current_SEQUENCE_TARGET,
					current_PRIMER_TASK,
					current_PRIMER_OPT_SIZE,
					current_PRIMER_MIN_SIZE,
					current_PRIMER_MAX_SIZE,
					current_PRIMER_PRODUCT_SIZE_RANGE,
					current_P3_FILE_FLAG,
					current_PRIMER_EXPLAIN_FLAG,
					current_PRIMER_LEFT_EXPLAIN,
					current_PRIMER_RIGHT_EXPLAIN,
					current_PRIMER_INTERNAL_EXPLAIN,
					current_PRIMER_PAIR_EXPLAIN,
					current_PRIMER_LEFT_NUM_RETURNED,
					current_PRIMER_RIGHT_NUM_RETURNED,
					current_PRIMER_INTERNAL_NUM_RETURNED,
					current_PRIMER_PAIR_NUM_RETURNED,
					current_PrimerPairs
				)
				retarray.append(tmp_PrimerPair)
				current_SEQUENCE_ID                  = None
				current_SEQUENCE_TEMPLATE            = None
				current_SEQUENCE_TARGET              = None
				current_PRIMER_TASK                  = None
				current_PRIMER_OPT_SIZE              = None
				current_PRIMER_MIN_SIZE              = None
				current_PRIMER_MAX_SIZE              = None
				current_PRIMER_PRODUCT_SIZE_RANGE    = None
				current_P3_FILE_FLAG                 = None
				current_PRIMER_EXPLAIN_FLAG          = None
				current_PRIMER_LEFT_EXPLAIN          = None
				current_PRIMER_RIGHT_EXPLAIN         = None
				current_PRIMER_INTERNAL_EXPLAIN      = None
				current_PRIMER_PAIR_EXPLAIN          = None
				current_PRIMER_LEFT_NUM_RETURNED     = None
				current_PRIMER_RIGHT_NUM_RETURNED    = None
				current_PRIMER_INTERNAL_NUM_RETURNED = None
				current_PRIMER_PAIR_NUM_RETURNED     = None
				current_PrimerPairs                  = {}


			else:
				leftterm, rightterm = l.strip().split("=")
				if leftterm == "SEQUENCE_ID":
					current_SEQUENCE_ID = rightterm
				elif leftterm == "SEQUENCE_TEMPLATE":
					current_SEQUENCE_TEMPLATE = rightterm
				elif leftterm == "SEQUENCE_TARGET":
					current_SEQUENCE_TARGET = rightterm
				elif leftterm == "PRIMER_TASK":
					current_PRIMER_TASK = rightterm
				elif leftterm == "PRIMER_OPT_SIZE":
					current_PRIMER_OPT_SIZE = rightterm
				elif leftterm == "PRIMER_MIN_SIZE":
					current_PRIMER_MIN_SIZE = rightterm
				elif leftterm == "PRIMER_MAX_SIZE":
					current_PRIMER_MAX_SIZE = rightterm
				elif leftterm == "PRIMER_PRODUCT_SIZE_RANGE":
					current_PRIMER_PRODUCT_SIZE_RANGE = rightterm
				elif leftterm == "P3_FILE_FLAG":
					current_P3_FILE_FLAG = rightterm
				elif leftterm == "PRIMER_EXPLAIN_FLAG":
					current_PRIMER_EXPLAIN_FLAG = rightterm
				elif leftterm == "PRIMER_LEFT_EXPLAIN":
					current_PRIMER_LEFT_EXPLAIN = rightterm
				elif leftterm == "PRIMER_RIGHT_EXPLAIN":
					current_PRIMER_RIGHT_EXPLAIN = rightterm
				elif leftterm == "PRIMER_INTERNAL_EXPLAIN":
					current_PRIMER_INTERNAL_EXPLAIN = rightterm
				elif leftterm == "PRIMER_PAIR_EXPLAIN":
					current_PRIMER_PAIR_EXPLAIN = rightterm
				elif leftterm == "PRIMER_LEFT_NUM_RETURNED":
					current_PRIMER_LEFT_NUM_RETURNED = rightterm
				elif leftterm == "PRIMER_RIGHT_NUM_RETURNED":
					current_PRIMER_RIGHT_NUM_RETURNED = rightterm
				elif leftterm == "PRIMER_INTERNAL_NUM_RETURNED":
					current_PRIMER_INTERNAL_NUM_RETURNED = rightterm
				elif leftterm == "PRIMER_PAIR_NUM_RETURNED":
					current_PRIMER_PAIR_NUM_RETURNED = rightterm
				else:
					current_PrimerPairs[leftterm] = rightterm
	return retarray

				
def extract_primer_pairs(primers_candidates):
	retdict = {}
	for each_primer in primers_candidates:
		if len(each_primer.Primer_Pairs) != 0:
			each_pairs = []
			i = 0
			while True:
				if each_primer.Primer_Pairs.get(f"PRIMER_PAIR_{i}_PENALTY") == None:
					break
				current_PRIMER_PAIR_PENALTY         = each_primer.Primer_Pairs.get(f"PRIMER_PAIR_{i}_PENALTY")
				current_PRIMER_LEFT_PENALTY         = each_primer.Primer_Pairs.get(f"PRIMER_LEFT_{i}_PENALTY")
				current_PRIMER_RIGHT_PENALTY        = each_primer.Primer_Pairs.get(f"PRIMER_RIGHT_{i}_PENALTY")
				current_PRIMER_INTERNAL_PENALTY     = each_primer.Primer_Pairs.get(f"PRIMER_INTERNAL_{i}_PENALTY")
				current_PRIMER_LEFT_SEQUENCE        = each_primer.Primer_Pairs.get(f"PRIMER_LEFT_{i}_SEQUENCE")
				current_PRIMER_RIGHT_SEQUENCE       = each_primer.Primer_Pairs.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
				current_PRIMER_INTERNAL_SEQUENCE    = each_primer.Primer_Pairs.get(f"PRIMER_INTERNAL_{i}_SEQUENCE")
				current_PRIMER_LEFT                 = each_primer.Primer_Pairs.get(f"PRIMER_LEFT_{i}")
				current_PRIMER_RIGHT                = each_primer.Primer_Pairs.get(f"PRIMER_RIGHT_{i}")
				current_PRIMER_INTERNAL             = each_primer.Primer_Pairs.get(f"PRIMER_INTERNAL_{i}")
				current_PRIMER_LEFT_TM              = each_primer.Primer_Pairs.get(f"PRIMER_LEFT_{i}_TM")
				current_PRIMER_RIGHT_TM             = each_primer.Primer_Pairs.get(f"PRIMER_RIGHT_{i}_TM")
				current_PRIMER_INTERNAL_TM          = each_primer.Primer_Pairs.get(f"PRIMER_INTERNAL_{i}_TM")
				current_PRIMER_LEFT_GC_PERCENT      = each_primer.Primer_Pairs.get(f"PRIMER_LEFT_{i}_GC_PERCENT")
				current_PRIMER_RIGHT_GC_PERCENT     = each_primer.Primer_Pairs.get(f"PRIMER_RIGHT_{i}_GC_PERCENT")
				current_PRIMER_INTERNAL_GC_PERCENT  = each_primer.Primer_Pairs.get(f"PRIMER_INTERNAL_{i}_GC_PERCENT")
				current_PRIMER_INTERNAL_SELF_ANY_TH = each_primer.Primer_Pairs.get(f"PRIMER_INTERNAL_{i}_SELF_ANY_TH")
				current_PRIMER_LEFT_SELF_ANY_TH     = each_primer.Primer_Pairs.get(f"PRIMER_LEFT_{i}_SELF_ANY_TH")
				current_PRIMER_RIGHT_SELF_ANY_TH    = each_primer.Primer_Pairs.get(f"PRIMER_RIGHT_{i}_SELF_ANY_TH")
				current_PRIMER_INTERNAL_SELF_END_TH = each_primer.Primer_Pairs.get(f"PRIMER_INTERNAL_{i}_SELF_END_TH")
				current_PRIMER_LEFT_SELF_END_TH     = each_primer.Primer_Pairs.get(f"PRIMER_LEFT_{i}_SELF_END_TH")
				current_PRIMER_RIGHT_SELF_END_TH    = each_primer.Primer_Pairs.get(f"PRIMER_RIGHT_{i}_SELF_END_TH")
				current_PRIMER_LEFT_HAIRPIN_TH      = each_primer.Primer_Pairs.get(f"PRIMER_LEFT_{i}_HAIRPIN_TH")
				current_PRIMER_RIGHT_HAIRPIN_TH     = each_primer.Primer_Pairs.get(f"PRIMER_RIGHT_{i}_HAIRPIN_TH")
				current_PRIMER_INTERNAL_HAIRPIN_TH  = each_primer.Primer_Pairs.get(f"PRIMER_INTERNAL_{i}_HAIRPIN_TH")
				current_PRIMER_LEFT_END_STABILITY   = each_primer.Primer_Pairs.get(f"PRIMER_LEFT_{i}_END_STABILITY")
				current_PRIMER_RIGHT_END_STABILITY  = each_primer.Primer_Pairs.get(f"PRIMER_RIGHT_{i}_END_STABILITY")
				current_PRIMER_PAIR_COMPL_ANY_TH    = each_primer.Primer_Pairs.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH")
				current_PRIMER_PAIR_COMPL_END_TH    = each_primer.Primer_Pairs.get(f"PRIMER_PAIR_{i}_COMPL_END_TH")
				current_PRIMER_PAIR_PRODUCT_SIZE    = each_primer.Primer_Pairs.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE")
				current_PRIMER_PAIR_PRODUCT_TM      = each_primer.Primer_Pairs.get(f"PRIMER_PAIR_{i}_PRODUCT_TM")
				tmp_PrimerPair = PrimerPair(current_PRIMER_PAIR_PENALTY, current_PRIMER_LEFT_PENALTY, current_PRIMER_RIGHT_PENALTY, current_PRIMER_INTERNAL_PENALTY, current_PRIMER_LEFT_SEQUENCE, current_PRIMER_RIGHT_SEQUENCE, current_PRIMER_INTERNAL_SEQUENCE, current_PRIMER_LEFT, current_PRIMER_RIGHT, current_PRIMER_INTERNAL, current_PRIMER_LEFT_TM, current_PRIMER_RIGHT_TM, current_PRIMER_INTERNAL_TM, current_PRIMER_LEFT_GC_PERCENT, current_PRIMER_RIGHT_GC_PERCENT, current_PRIMER_INTERNAL_GC_PERCENT, current_PRIMER_INTERNAL_SELF_ANY_TH, current_PRIMER_LEFT_SELF_ANY_TH, current_PRIMER_RIGHT_SELF_ANY_TH, current_PRIMER_INTERNAL_SELF_END_TH, current_PRIMER_LEFT_SELF_END_TH, current_PRIMER_RIGHT_SELF_END_TH, current_PRIMER_LEFT_HAIRPIN_TH, current_PRIMER_RIGHT_HAIRPIN_TH, current_PRIMER_INTERNAL_HAIRPIN_TH, current_PRIMER_LEFT_END_STABILITY, current_PRIMER_RIGHT_END_STABILITY, current_PRIMER_PAIR_COMPL_ANY_TH, current_PRIMER_PAIR_COMPL_END_TH, current_PRIMER_PAIR_PRODUCT_SIZE, current_PRIMER_PAIR_PRODUCT_TM)
				each_pairs.append(tmp_PrimerPair.__dict__)
				i = i + 1
			retdict[each_primer.SEQUENCE_ID] = {"Primer3_input": each_primer.__dict__, "Primer3_output": each_pairs}
	return retdict



def main():
	parser = argparse.ArgumentParser(description = "primer3 result parser")
	parser.add_argument("primer3_result", metavar = "primer3_result", type = str, help = "primer3 results file name")
	parser.add_argument("-o",    metavar = "output_file",    type = str, default = "sys.stdout", help = "output file name (default = sys.stdout)")
	parser.add_argument("--fasta", action='store_true', help = "output as fasta")
	parser.add_argument("--trimmed_fasta", action='store_true', help = "output as trimmed fasta")
	parser.add_argument("--tempfmt", action='store_true', help = "output as temporaly format")

	args = parser.parse_args()
	filename = args.primer3_result
	primers_candidates = primer3_result_parser(filename)
	primer_pairs_dict = extract_primer_pairs(primers_candidates)

	output_file = open(args.o, "w") if args.o is not None else sys.stdout
	if args.fasta:
		print(fmt4fasta(primer_pairs_dict), file = output_file)
	elif args.trimmed_fasta:
		print(fmt4fasta(primer_pairs_dict, trim_flag = True), file = output_file)
	elif args.tempfmt:
		print(tempfmt(primer_pairs_dict), file = output_file)
	else:
		print(json.dumps(primer_pairs_dict, indent = 2), file = output_file)
		

if __name__ == '__main__':
	main()