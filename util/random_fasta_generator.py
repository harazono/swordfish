#! /usr/bin/env python3
import random


random_sequences = []
for i in range(10):
	random_str = "".join([random.choice(["A", "C", "G", "T"]) for x in range(80)])
	random_sequences.append(random_str)

for i in range(1, 201):
	num = str(i).zfill(3)
	print(f">dummy_sequence_{num} {i}th record")
	for j in range(5):
		r = round(random.random() * 100) % 10
		print(random_sequences[r])