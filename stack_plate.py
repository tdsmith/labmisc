#!/usr/bin/env python3

import sys

buf = sys.stdin.read().replace("\r\n", "\n").replace("\r", "\n").split("\n")
buf = [line.split('\t') for line in buf]

max_cols = 0
for line in buf:
    line[-1] = line[-1].strip()
    if len(line) > max_cols:
        max_cols = len(line)

rowlabels = "ABCDEFGH"

print("Well,Label")

for col in range(max_cols):
    for row, line in enumerate(buf):
        if col >= len(line):
            continue
        label = "{}{:02}".format(rowlabels[row], col+1)
        contents = line[col]
        if not contents:
            continue
        print("{},{}".format(label, contents))
