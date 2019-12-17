#!/usr/bin/env python
import argparse
import numpy as np
import sys
import os

# helices: G, H, I
# beta: E, B
# turn: T
# bend : S
# none: NonSeq

sec_type = {"G" : "helix", "H" : "helix", "I" : "helix", "T": "turn", "S" : "bend", "C": "misc", "NoSeq": "none", "*EOS*" : "beg", "*SOS*" : "end"}

def scoring(proposed, actual):
	# 1 if same
	# one if within same set (see above)
	if len(proposed) > len(actual):
		x = proposed
		y = actual[:len(proposed)]
	else:
		x = proposed[:len(actual)]
		y = actual
	score = 0
	for i, j in zip(x, y):
		if i == j:
			score += 1
		elif sec_type[i] == sec_type[j]:
			score += 0.5
		elif (sec_type[i] == "helix" and sec_type[j] == "turn") or (sec_type[j] == "helix" and sec_type[i] == "turn"):
			score += 0.25
	return score, max(len(proposed), len(actual))

def valid_file(parser, arg):
    if arg and not os.path.exists(arg):
        parser.error('The file doesn\'t exist: {}'.format(arg))
    else:
        return arg

    return args


lst = raw_input().split()

total_score = 0
total_len = 0
count = 0
# fi = sys.
# print("FILE", fi)
for line in sys.stdin.readlines():
	line = line.split()
	score, lenx = scoring(line, lst[count])
	total_score += score
	total_len += lenx
		# print(score)
	count += 1
if total_len == 0:
	total_len = 1
print("{}".format(100 * total_score / total_len))