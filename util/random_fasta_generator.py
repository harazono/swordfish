#! /usr/bin/env python3
import random


random_sequences = []
for i in range(10):
	random_str = "".join(["".join(random.sample(["A", "C", "G", "T"], 4)) for x in range(48)])
	random_sequences.append(random_str)

for i in range(1, 201):
	num = str(i).zfill(3)
	print(f">dummy_sequence_{num} {i}th record")
	for j in range(1):
		r = round(random.random() * 100) % 10
		print(random_sequences[r])
