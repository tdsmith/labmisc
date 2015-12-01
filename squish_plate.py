#!/usr/bin/env python3
import sys

import numpy as np

buf = sys.stdin.read().replace("\r\n", "\n").replace("\r", "\n").split("\n")
buf = [line.split('\t') for line in buf]

max_rows = 0
max_cols = 0
for well, value in buf:
    row = ord(well[0]) - ord('A')
    col = int(well[1:]) - 1
    max_rows = max(max_rows, row)
    max_cols = max(max_cols, col)

plate = np.empty((max_rows+1, max_cols+1)) * np.nan

for well, value in buf:
    row = ord(well[0]) - ord('A')
    col = int(well[1:]) - 1
    plate[row, col] = float(value)

np.savetxt(sys.stdout.buffer, plate, delimiter='\t')
