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
	parser.add_argument("fasta",     metavar = "fasta",         type = str, help = "primers.fa file name")
	parser.add_argument("blast",     metavar = "blast",         type = str, help = "blast output file name. outfmt must be '6 qseqid sseqid sacc slen qstart qend sstart send qseq sseq evalue length staxid staxids ssciname scomname'")
	parser.add_argument("primer3",   metavar = "primer3",       type = str, help = "primer3 output file (json)")
	parser.add_argument("--discard", metavar = "namelist",      type = str, nargs="+", help = "list of sequence names to be discarded. one name per line.")
	parser.add_argument("-o",        metavar = "output_prefix", type = str, default = "final_result", help = "output file name (default = final_result)")
	#parser.add_argument("--fasta", action='store_true', help = "output as fasta")
	args = parser.parse_args()
	filename = args.fasta
	fasta_ids = set()

	with open(filename) as handle:
		for record in SeqIO.parse(handle, "fasta"):
			fasta_ids.add(record.id)

	overlook_taxon_ids = [
		4528, #Oryza longistaminata
		4529, #Oryza rufipogon
		4530, #Oryza sativa
		39946, #Oryza sativa Indica Group
		39947, #Oryza sativa Japonica Group
		1080340, #Oryza sativa Japonica Group x Oryza sativa Indica Group
		1736656, #Oryza sativa tropical japonica subgroup
		1736657, #Oryza sativa temperate japonica subgroup
		1736658, #Oryza sativa aromatic subgroup
		1736659, #Oryza sativa aus subgroup
		1771142, #Oryza sativa indica subgroup
		2998809, #Oryza sativa f. spontanea
		1050722, #Oryza sativa Indica Group x Oryza sativa Japonica Group
		256318, #metagenome
		408170, #human gut metagenome
		408172, #marine metagenome
		410658, #soil metagenome
		410661, #mouse gut metagenome
		412755, #marine sediment metagenome
		496920, #saltern metagenome
		527639, #wastewater metagenome
		527640, #microbial mat metagenome
		652676, #hydrothermal vent metagenome
		717931, #groundwater metagenome
		749906, #gut metagenome
		749907, #sediment metagenome
		797283, #ant fungus garden metagenome
		904678, #hypersaline lake metagenome
		1076179, #bioreactor metagenome
		1911570, #Mangrovicoccus ximenensis
		2509717, #Ciceribacter ferrooxidans
	]
	#overlook_taxon_ids = [4530, 39947, 1080340, 1050722, 1736656, 1736657, 1736658, 1736659, 1771142, 2998809, "N/A"]
	#4529はOryza rufipogon, 野生の稲
	#4528はOryza longistaminata
	blast_results = []
	with open(args.blast) as f:
		reader = csv.reader(f, delimiter='\t')
		for each_record in reader:
			each_record_Obj = BlastResult(*each_record)
			taxon_id = -1
			try:
				taxon_id = int(each_record_Obj.staxid)
			except:
				taxon_id = each_record_Obj.staxid
			if taxon_id not in overlook_taxon_ids:
				blast_results.append(each_record_Obj)


	"""
			if taxon_id in overlook_taxon_ids:
				continue
			if "N/A" in each_record_Obj.scomname:
				continue
			if "N/A" in each_record_Obj.ssciname:
				continue
			if "metagenome" in each_record_Obj.scomname:
				continue
			if "metagenome" in each_record_Obj.ssciname:
				continue
			blast_results.append(each_record_Obj)
	"""

	primer3_info = None
	with open (args.primer3, "r") as f:
		primer3_info = json.load(f)


	primer_blasthit_dict = {k: set() for k in fasta_ids}
	for each_hit in blast_results:
		primer_blasthit_dict[each_hit.qseqid].add(f"{each_hit.sseqid};{each_hit.staxid};{each_hit.scomname};{each_hit.ssciname}")



	blast_trapped_seq_ids = set()
	for each_hit in blast_results:
		blast_trapped_seq_ids.add(each_hit.qseqid)

	discard_set = set()
	if args.discard is not None:
		discard_filenames = args.discard
		for discard_filename in discard_filenames:
			with open(discard_filename) as f:
				for line in f:
					discard_set.add(line.strip())



	survivor = fasta_ids - blast_trapped_seq_ids - discard_set
	survivor_pair = set()
	for each_primer in survivor:
		id, idx, side = each_primer.split("_")
		pair = ""
		if side == "L":
			pair = "R"
		else:
			pair = "L"
		#pair_full_name = id + "_" + idx + pair#間違い
		pair_full_name = id + "_" + idx + "_" + pair
		if pair_full_name in survivor:
			survivor_pair.add(id + "_" + idx)
		else:
			pass

	print(list(primer3_info.keys())[0], file = sys.stderr) #6d83378ab5107afd062baf2cca8e913
	print(list(survivor)[0], file = sys.stderr) #e116136dc273515db5cee535731c145_2_R
	if len(survivor_pair) > 0:
		print(list(survivor_pair)[0], file = sys.stderr) #62baf2cca8e91329bcaba327c863b6b_4

	#print(survivor, file = sys.stderr)#84db709cd45706d4ee535731c145dbe_4_L
	report_file            = open(args.o + ".report", mode = "w")
	survivor_tsv_file      = open(args.o + ".survivor.tsv", mode = "w")
	survivor_pair_tsv_file = open(args.o + ".survivor_pair.tsv", mode = "w")
	survivor_namelist_file = open(args.o + ".survivor_name.txt", mode = "w")

	print(f"total count of input sequence                  : {len(fasta_ids)}"            , file = report_file)
	print(f"total count of blast hits                      : {len(blast_results)}"        , file = report_file)
	print(f"cardinarity of input which was trapped by blast: {len(blast_trapped_seq_ids)}", file = report_file)
	print(f"survivor                                       : {len(survivor)}"             , file = report_file)
	print(f"survivor pair                                  : {len(survivor_pair)}"        , file = report_file)

	print("\t".join(["primer id", "left primer", "right primer", "primer left Tm", "primer right Tm", "primer pair product Tm", "survived side", "trapped side", "blast hits"]), file = survivor_tsv_file)
	for each_survivor in survivor:
		survivor_info_raw = primer3_info[each_survivor.split("_")[0]]
		survivor_index    = int(each_survivor.split("_")[1])
		primer_side       = each_survivor.split("_")[2]
		survivor_info     = survivor_info_raw["Primer3_output"][survivor_index]
		primer_pair       = each_survivor[:-1] + "L" if primer_side == "R" else each_survivor[:-1] + "R"
		blast_hits = primer_blasthit_dict[primer_pair]
		print("\t".join([str(x) for x in [each_survivor.split("_")[0], survivor_info["PRIMER_LEFT_SEQUENCE"], survivor_info["PRIMER_RIGHT_SEQUENCE"], survivor_info["PRIMER_LEFT_TM"], survivor_info["PRIMER_RIGHT_TM"], survivor_info["PRIMER_PAIR_PRODUCT_TM"], primer_side, "L" if primer_side == "R" else "R", blast_hits]]), file = survivor_tsv_file)
		print(each_survivor, file = survivor_namelist_file)


	print("\t".join(["primer id", "left primer", "right primer", "primer left Tm", "primer right Tm", "primer pair product Tm"]), file = survivor_pair_tsv_file)
	for each_survivor_pair in survivor_pair:
		survivor_info_raw = primer3_info[each_survivor_pair.split("_")[0]]
		survivor_index    = int(each_survivor_pair.split("_")[1])
		survivor_info     = survivor_info_raw["Primer3_output"][survivor_index]
		print("\t".join([str(x) for x in [each_survivor_pair, survivor_info["PRIMER_LEFT_SEQUENCE"], survivor_info["PRIMER_RIGHT_SEQUENCE"], survivor_info["PRIMER_LEFT_TM"], survivor_info["PRIMER_RIGHT_TM"], survivor_info["PRIMER_PAIR_PRODUCT_TM"]]]), file = survivor_pair_tsv_file)


if __name__ == '__main__':
	main()

