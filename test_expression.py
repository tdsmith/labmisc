#!/usr/bin/env python

import unittest
from expression import expression, average_cq, rank_genes
from types import *

mean = lambda seq: float(sum(seq))/len(seq)

class TestExpression(unittest.TestCase):
    def setUp(self):
        self.sample_list = [
                ['Exp1', 'Sample 1', 5.0],
                ['Exp1', 'Sample 2', 6.0],
                ['Ref1', 'Sample 1', 3.3],
                ['Ref1', 'Sample 2', 3.3]]

    def test_returns_tuple(self):
        r = expression(self.sample_list, 'Ref1', 'Sample 1')
        self.assertIsInstance(r, TupleType)
        self.assertEqual(len(r), 3)
        self.assertIsInstance(r[0], DictType)
        self.assertIsInstance(r[1], ListType)
        self.assertIsInstance(r[2], ListType)

    def test_does_not_skip_good_genes_or_samples(self):
        sample_d, ignored_genes, ignored_samples = expression(self.sample_list, 'Ref1', 'Sample 1')
        self.assertEqual(ignored_genes, [])
        self.assertEqual(ignored_samples, [])

    def test_skips_genes_with_no_references(self):
        sample_list = self.sample_list + [['Exp2', 'Sample 2', 10.0]]
        sample_d, ignored_genes, ignored_samples = expression(sample_list, 'Ref1', 'Sample 1')
        self.assertIn('Exp2', ignored_genes)

    def test_skips_samples_with_no_reference(self):
        sample_list = self.sample_list + [['Exp1', 'Sample 3', 10.0]]
        sample_d, ignored_genes, ignored_samples = expression(sample_list, 'Ref1', 'Sample 1')
        self.assertIn('Sample 3', ignored_samples)

    def test_right_ballpark(self):
        sample_d, ignored_genes, ignored_samples = expression(self.sample_list, 'Ref1', 'Sample 1')
        self.assertAlmostEqual(mean(sample_d['Exp1']['Sample 2']), 0.5)

    def test_silences_near_ntc(self):
        addme =  [['Exp2', 'Sample 1', 10.0], ['Exp2', 'Sample 2', 5.0]]
        sample_list = self.sample_list + addme
        sample_d, ignored_genes, ignored_samples = expression(sample_list, 'Ref1', 'Sample 1')
        self.assertAlmostEqual(mean(sample_d['Exp2']['Sample 2']), 2**5)
        
        sample_list = sample_list + [['Exp2', 'NTC', 11.0]]
        sample_d, ignored_genes, ignored_samples = expression(sample_list, 'Ref1', 'Sample 1')
        self.assertIn('Exp2', ignored_genes)


class TestAverage(unittest.TestCase):
    def test_averaging(self):
        self.assertAlmostEqual(average_cq((3.0,5.0)), 3.678071905112638)
        self.assertAlmostEqual(average_cq((3.0,4.0,5.0)), 3.7776075786635523)

class TestRankGenes(unittest.TestCase):
    def test_ranking(self):
        with open('vandesompele-2002-cq.txt') as f:
            buf = f.readlines()
        gene_names = buf[0][:-1].split('\t')[2:]
        # slicing to cut off newlines on each line and skip header row
        # yields ['Tissue', 'Sample', Cq0, Cq1, Cq2...]
        buf = [line[:-1].split('\t') for line in buf[1:] if len(line) > 2]
        d = {}
        for line in buf:
            for (index, name) in enumerate(gene_names):
                d.setdefault(line[0], []).append([name, line[1], float(line[2+index])])
        # e.g. d['Fib'] = [['ACTB', 'Fib1', 12.34], ...]
        key = {'Neuroblastoma': ['HPRT1', 'GAPD', 'SDHA', 'UBC', 'HMBS', 'YWHAZ', 'TBP', 'ACTB', 'RPL13A', 'B2M'],
               'Fib': ['HPRT1', 'GAPD', 'YWHAZ', 'UBC', 'ACTB', 'TBP', 'SDHA', 'RPL13A', 'B2M', 'HMBS'],
               'Leukocyte': ['UBC', 'YWHAZ', 'B2M', 'GAPD', 'RPL13A', 'TBP', 'SDHA', 'HPRT1', 'HMBS', 'ACTB'],
               'BM': ['UBC', 'RPL13A', 'YWHAZ', 'HPRT1', 'GAPD', 'SDHA', 'TBP', 'HMBS', 'B2M', 'ACTB'],
               'Pool': ['SDHA', 'GAPD', 'HMBS', 'HPRT1', 'TBP', 'UBC', 'RPL13A', 'YWHAZ', 'ACTB', 'B2M']}
        for tissue in d:
            ranked = rank_genes(d[tissue], gene_names, d[tissue][0][1])
            # no preference between the two best genes
            self.assertIn(ranked[0], key[tissue][:2])
            self.assertIn(ranked[1], key[tissue][:2])
            self.assertEqual(ranked[2:], key[tissue][2:])

if __name__ == '__main__':
    unittest.main()

