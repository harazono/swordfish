#! /usr/bin/env python3

import requests
import time
import sys
from bs4 import BeautifulSoup
import json
import pprint
import argparse
import time
pp = pprint.PrettyPrinter(indent=4)

def main():
	parser = argparse.ArgumentParser(description = "NCBI record checker")
	parser.add_argument("input", type=str, metavar="file", help="input file")
	args = parser.parse_args()
	filename = args.input
	prefix = "https://www.ncbi.nlm.nih.gov/nuccore/"
	with open(filename, "r") as f:
		for line in f:
			url = prefix + line.strip()
			res = requests.get(prefix + line.strip())
			soup = BeautifulSoup(res.text, 'html.parser')
			title_text = soup.find('title').get_text()
			print(line.strip(), end = "\t")
			print(title_text, end = "\t")
			print(url)
			time.sleep(10)


if __name__ == "__main__":
	main()