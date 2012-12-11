#!/usr/bin/env python

import sys
from expression import expression

def reconstitute(filename, ref_gene, ref_sample):
    # assumes headings: Well Sample Cq Target
    # further assume that line ends may have been munged by Excel 2011
    f = open(filename, 'r')
    buf = f.read().replace('\r', '\n').split('\n')
    norm = lambda line: [line[3], line[1], float(line[2])]
    buf = [norm(line.split('\t')) for line in buf[1:] if len(line) > 1]
    # now each line in buf looks like: ['Gapdh', 'Sample1', 20.95]
    sample_d, ignored_genes, ignored_samples = expression(buf, ref_gene, ref_sample)
    if ignored_genes: print >> sys.stderr, 'ignored genes:', ignored_genes
    if ignored_samples: print >> sys.stderr, 'ignored samples:', ignored_samples
    print 'Target\tSample\tValues'
    for target in sorted(sample_d.keys()):
        for sample in sorted(sample_d[target].keys()):
            print '%s\t%s\t%s' % (target, sample, '\t'.join([str(i) for i in sample_d[target][sample]]))


def main():
    reconstitute(sys.argv[1], sys.argv[2], sys.argv[3])

if __name__ == '__main__':
    main()

