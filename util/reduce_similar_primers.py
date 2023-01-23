#! /usr/bin/env python3

import argparse
import re
import pprint
import csv
import sys
pp = pprint.PrettyPrinter(indent = 2)

species_name_identifier = re.compile(r"'([\w\s]+)\(")


class Primer():
	def __init__(self,primer_id,left_primer,right_primer,primer_left_Tm,primer_right_Tm,primer_pair_product_Tm,survived_side,trapped_side,blast_hits) -> None:
		self.primer_id              = primer_id
		self.left_primer            = left_primer
		self.right_primer           = right_primer
		self.primer_left_Tm         = primer_left_Tm
		self.primer_right_Tm        = primer_right_Tm
		self.primer_pair_product_Tm = primer_pair_product_Tm
		self.survived_side          = survived_side
		self.trapped_side           = trapped_side
		self.blast_hits             = set(re.findall(species_name_identifier, blast_hits))
		self.blast_hits_str         = ", ".join(self.blast_hits)
		#self.blast_hits             = json.loads(blast_hits.replace(":", ";").replace("'", '"').replace("{", "[").replace("}", "]"))
		# {'Mangrovicoccus ximenensis(1911570: Mangrovicoccus ximenensis)', 'N/A(N/A: N/A)', 'Ciceribacter ferrooxidans(2509717: Ciceribacter ferrooxidans)'}
	def __str__(self):
		return "\t".join([f"{k}" for k in self.__dict__.values()])
		#f"{self.left_primer}\t{self.right_primer}"


def main():
	parser = argparse.ArgumentParser(description = "reduce, merge similar primers")
	parser.add_argument("primer3_result", metavar = "primer3_result", type = str, help = "primer3 results file name")
	parser.add_argument("-s",    metavar = "size",    type = int, default = 15, help = "Length from 5' used to determine if they are considered the same")
	parser.add_argument("-d", action="store_true", help = "Print primer pairs whose wrong target does not have any intersection")
	args = parser.parse_args()
	filename = args.primer3_result
	primer_data = []
	species_redlist = ['Lolium rigidum',# ボウムギ
	'Anopheles albimanus', #メキシコの蚊
	'red palm weevil', # ヤシオオオサゾウムシ
	'Ciceribacter ferrooxidans' #ダイズとクロスしそうな配列を除く
	]
	primerid_redlist = [
		"2cca8e913286f2aed8705996f0aca7a",
		"2cca8e913286f2af61c1665bc2b29e9",
		"2cca8e913286f2af61c1665bc2b29e9",
		"949717c1a66bcb605706d73b94d5cc7",
		"525c5f0699af2d845706d73b94d5cc7",
		"cb32a3a44ca1bcab61c1665bc2b29e9",
		"2cca8e913286f2aed8705996f0aca7a",
		"cb32a3a44ca1bcaad8705996f0aca7a",
		"949717c1a66bcb615c1b5cee535731c",
		"525c5f0699af2d855c1b5cee535731c",
		"2cca8e913286f2aed8705996f0aca7a"
	]
	with open(filename, "r") as f:
		header = next(csv.reader(f, delimiter="\t"))
		reader = csv.reader(f, delimiter="\t")
		for line in reader:
			try:
				primer_data.append(Primer(*line))
			except Exception as e:
				#print(e, line)
				pass
	primer_merged_dict = {}
	for each_primer in primer_data:
		left_primer_subseq  = each_primer.left_primer[-1 * args.s:]
		right_primer_subseq = each_primer.right_primer[-1 * args.s:]
		key = (left_primer_subseq, right_primer_subseq)
		redlist_check_flag = False
		for i in species_redlist:
			if i in set(each_primer.blast_hits):
				redlist_check_flag = True
				break
		if each_primer.primer_id in primerid_redlist:
			redlist_check_flag == True
		if redlist_check_flag:
			continue
		if key in primer_merged_dict:
			primer_merged_dict[key].append(each_primer)
		else:
			primer_merged_dict[key] = [each_primer]
	print(len(primer_merged_dict), file = sys.stderr)



	if not args.d:
		print("似たプライマーを除去して残った数：", len(primer_merged_dict.keys()))
		print("primer group	primer id	left primer	right primer	primer left Tm	primer right Tm	primer pair product Tm	survived side	trapped side	blast hits")
		for k, v in primer_merged_dict.items():
			set_name = f"{k[0]}-{k[1]}"
			print_flag = True
			for each_primer in v:
				print(set_name, end = "\t") if print_flag else print(" ", end = "\t")
				print_flag = False
				print(each_primer)
			print()

	else:
		key_pair = {}
		safe_pair_cnt = 0
		for each_key_1 in primer_merged_dict.keys():
			for each_key_2 in primer_merged_dict.keys():
				if each_key_1 == each_key_2:
					continue
				blast_trap_in_set1 = set()
				for i in [x.blast_hits for x in primer_merged_dict[each_key_1]]:
					blast_trap_in_set1 |= i
				blast_trap_in_set2 = set()
				for i in [x.blast_hits for x in primer_merged_dict[each_key_2]]:
					blast_trap_in_set2 |= i
				primer_identical_flag = len(set([each_key_1[0], each_key_2[0], each_key_1[1], each_key_2[1]])) == 4
				offtarget_flag        = len(blast_trap_in_set1 & blast_trap_in_set2) == 0
				if primer_identical_flag and offtarget_flag:
					safe_pair_cnt = safe_pair_cnt + 1
					key_pair_str = f"{each_key_1[0]}-{each_key_1[1]}\t{each_key_2[0]}-{each_key_2[1]}"
					if key_pair_str not in key_pair:
						key_pair[key_pair_str] = (len(blast_trap_in_set1) + len(blast_trap_in_set2), blast_trap_in_set1, blast_trap_in_set2)
		print(safe_pair_cnt, file = sys.stderr)
		sorted_key_pair = sorted(key_pair.items(), key = lambda x :x[1][0])
		print("primer group	primer id	left primer	right primer	primer left Tm	primer right Tm	primer pair product Tm	survived side	trapped side	blast hits")
		for i in sorted_key_pair:
			key_1 = tuple(i[0].split("\t")[0].split("-"))
			key_2 = tuple(i[0].split("\t")[1].split("-"))
			primer_str_1 = primer_merged_dict[key_1]
			primer_str_2 = primer_merged_dict[key_2]
			print_flag = True
			for j in primer_str_1:
				print(i[0].split("\t")[0], end = "\t") if print_flag else print(" ", end = "\t")
				print_flag = False
				print(str(j))
			print_flag = True
			for j in primer_str_2:
				print(i[0].split("\t")[1], end = "\t") if print_flag else print(" ", end = "\t")
				print_flag = False
				print(str(j))
			print()


			#print(i[0], end = "\t")
			#print(i[1][0], end = "\t")
			#print(",".join(i[1][1]), end = "\t")
			#print(",".join(i[1][2]), end = "\n")

if __name__ == '__main__':
	main()