#!/usr/bin/env python

import sys, xlrd
import numpy as np

def usage():
    print "Usage: assemble_samples.py layout.xlsx quantitation.xls > samples.tab"

def assemble_samples(layout_path, quant_path):
    A1 = np.array((3, 2))
    layout = xlrd.open_workbook(layout_path).sheet_by_name('Sheet1')
    quant = xlrd.open_workbook(quant_path).sheet_by_name('0')
    # make sure we're actually at A1
    if layout.cell_value(*(A1 + (0,-1))) != 'A' or layout.cell_value(*(A1 + (-1,0))) != 1:
        raise ValueError("A1 seems to be in the wrong place or the input files are swapped.")
    rows = 'ABCDEFGH'
    cols = np.arange(12)+1
    sample_index = {}
    for (i, row) in enumerate(rows):
        for (j, col) in enumerate(cols):
            value = layout.cell_value(*(A1+(i,j)))
            if value: sample_index['%s%02d' % (row, col)] = str(value)

    start_row = 1
    name_col = 1
    cq_col = 6
    cq_index = {}
    for row in range(96):
        name = quant.cell_value(start_row+row, name_col)
        value = quant.cell_value(start_row+row, cq_col) or 'nan'
        cq_index[name] = float(value)

    print 'Well\tSample\tCq\tTarget'
    wells = sorted(sample_index.keys())
    for well in wells:
        print '%s\t%s\t%f\t' % (well, sample_index[well], cq_index[well])

def main():
    if len(sys.argv) != 3:
        usage()
        sys.exit(1)
    try:
        assemble_samples(sys.argv[1], sys.argv[2])
    except Exception, e:
        usage()
        print e
        sys.exit(1)
    return 0

if __name__ == '__main__':
    main()
