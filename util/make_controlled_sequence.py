#! /usr/bin/env python3
import random
import itertools

sequence_pool = "ACGT"
sequence_permutations = ["".join(x) for x in list(itertools.permutations(sequence_pool, 4))]

random_sequences = []
for i in range(10):
    random_str = "".join(["".join(random.sample(sequence_permutations, 4)) for _ in range(20)])
    random_sequences.append(random_str)

for i in range(1, 201):
    num = str(i).zfill(3)
    print(f">dummy_sequence_{num} {i}th record")
    for j in range(1):
        r = round(random.random() * 100) % 10
        print(random_sequences[r])
