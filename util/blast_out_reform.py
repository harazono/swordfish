#! /usr/bin/env python3
import argparse
import pprint
from abc import ABC
from dataclasses import dataclass
import sys
pp = pprint.PrettyPrinter(indent = 2)

@dataclass
class Record(ABC):
	qseqid   : str
	sseqid   : str
	sacc     : str
	slen     : int
	qstart   : int
	qend     : int
	sstart   : int
	send     : int
	qseq     : str
	sseq     : str
	evalue   : str
	length   : int
	staxid   : str
	staxids  : str
	ssciname : str
	scomname : str

def blast_checker_1(primer_pairs, records):
	discarded_primer_pairs = set()
	for each_primer_pair in primer_pairs:
		subject_to_be_checked = list(filter(lambda x: x.seqid == each_primer_pair, records))
		roles_in_each_primer_pair = set()
		for i in subject_to_be_checked:
			roles_in_each_primer_pair.add(i.role)
		if len(roles_in_each_primer_pair) == 3:
			sseqids = set([x.sseqid for x in subject_to_be_checked])
			for each_sseqid in sseqids:
				role_L_coverage = [(x.sstart, x.send) for x in subject_to_be_checked if x.role == "L" and x.sseqid == each_sseqid]
				role_M_coverage = [(x.sstart, x.send) for x in subject_to_be_checked if x.role == "M" and x.sseqid == each_sseqid]
				role_R_coverage = [(x.sstart, x.send) for x in subject_to_be_checked if x.role == "R" and x.sseqid == each_sseqid]
				if role_L_coverage != [] and role_M_coverage != [] and role_R_coverage != []:
					# test all combination
					for l in role_L_coverage:
						for m in role_M_coverage:
							for r in role_R_coverage:
								if 0 < m[0] - l[1] < 76 and 0 < r[0] - m[1] < 76:
									discarded_primer_pairs.add(each_primer_pair)
	return discarded_primer_pairs

def blast_checker_2(records):
	return set([x.qseqid for x in records])

def main():
	parser = argparse.ArgumentParser(description = "blastn-short result reformer")
	parser.add_argument("-b", "--blast_result", required=True, nargs="+", type=str, metavar="Blast result", help="Blast results file (outfmt 6)")
	parser.add_argument("-p", "--primer_candidates",required=True, type=str, metavar="Primer candidates", help="Primer candidate name list file. One candidate in u128 in a line.")
	args = parser.parse_args()
	filenames = args.blast_result
	discarded_primer_pairs = set()

	primer_candidates_set = set()
	with open(args.primer_candidates) as f:
		for line in f:
			primer_candidates_set.add(int(line, 16))

	for each_db in filenames:
		records = []
		with open(each_db) as f:
			for line in f:
				elm = line.strip().split("\t")
				seqid = int(elm[0].split("-")[0], 16)
				role = elm[0].split("-")[1]
				qseqid, sseqid, sacc, slen, qstart, qend, sstart, send, qseq, sseq, evalue, length, staxid, staxids, ssciname, scomname = elm
				qseqid   = str(qseqid)
				sseqid   = str(sseqid)
				sacc     = str(sacc)
				slen     = int(slen)
				qstart   = int(qstart)
				qend     = int(qend)
				sstart   = int(sstart)
				send     = int(send)
				qseq     = str(qseq)
				sseq     = str(sseq)
				evalue   = str(evalue)
				length   = int(length)
				staxid   = str(staxid)
				staxids  = str(staxids)
				ssciname = str(ssciname)
				scomname = str(scomname)


				records.append(Record(qseqid,sseqid,sacc,slen,qstart,qend,sstart,send,qseq,sseq,evalue,length,staxid,staxids,ssciname,scomname))

			#discarded_primer_pairs = blast_checker_1(primer_pairs, records)
			discarded_primer_pairs = discarded_primer_pairs | blast_checker_2(records)



	remained = primer_candidates_set - discarded_primer_pairs
	#print(f"candidates: {len(primer_candidates_set)}", file = sys.stderr)
	#print(f"discarded:  {len(discarded_primer_pairs)}", file = sys.stderr)
	#print(f"remained:   {len(remained)}", file = sys.stderr)

	for i in remained:
		print(hex(i).replace("0x", ""))
	# for i in primer_candidates_set - discarded_primer_pairs:
	# 	print(f'{hex(i).replace("0x", "")}')
	# print(f"input file: {filename}", file = sys.stderr)
	# print(f"# of records in the inout file: {len(records)}", file = sys.stderr)
	# print(f"{len(primer_candidates_set)} - {len(discarded_primer_pairs)} = {len(primer_candidates_set - discarded_primer_pairs)}", file = sys.stderr)
	# print(f"{len(primer_candidates_set - discarded_primer_pairs)/len(primer_candidates_set)}", file = sys.stderr)



if __name__ == "__main__":
	main()



