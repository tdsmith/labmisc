#!/usr/bin/env python

import unittest
from expression import *
from types import *
from sys import stdout
import pandas as pd

mean = lambda seq: float(sum(seq))/len(seq)

class TestExpression(unittest.TestCase):
    def setUp(self):
        self.sample_list = [
                ['Exp1', 'Sample 1', 5.0],
                ['Exp1', 'Sample 2', 6.0],
                ['Ref1', 'Sample 1', 3.3],
                ['Ref1', 'Sample 2', 3.3]]
        self.sample_frame = pd.DataFrame(self.sample_list, columns=['Target', 'Sample', 'Cq'])

    def test_hello(self):
        r = expression_frame(self.sample_frame, 'Ref1', 'Sample 1')
        self.sample_frame['RelExp'] = r

    def test_returns_tuple(self):
        r = expression(self.sample_frame, 'Ref1', 'Sample 1')
        self.assertIsInstance(r, TupleType)
        self.assertEqual(len(r), 3)
        self.assertIsInstance(r[0], DictType)
        self.assertIsInstance(r[1], ListType)
        self.assertIsInstance(r[2], ListType)

    def test_does_not_skip_good_genes_or_samples(self):
        sample_d, ignored_genes, ignored_samples = expression(self.sample_frame, 'Ref1', 'Sample 1')
        self.assertEqual(ignored_genes, [])
        self.assertEqual(ignored_samples, [])

    def test_skips_genes_with_no_references(self):
        sample_list = self.sample_list + [['Exp2', 'Sample 2', 10.0]]
        sample_frame = pd.DataFrame(sample_list, columns=['Target', 'Sample', 'Cq'])
        sample_d, ignored_genes, ignored_samples = expression(sample_frame, 'Ref1', 'Sample 1')
        self.assertIn('Exp2', ignored_genes)

    def test_skips_samples_with_no_reference(self):
        sample_list = self.sample_list + [['Exp1', 'Sample 3', 10.0]]
        sample_frame = pd.DataFrame(sample_list, columns=['Target', 'Sample', 'Cq'])
        sample_d, ignored_genes, ignored_samples = expression(sample_frame, 'Ref1', 'Sample 1')
        self.assertIn('Sample 3', ignored_samples)

    def test_right_ballpark(self):
        sample_d, ignored_genes, ignored_samples = expression(self.sample_frame, 'Ref1', 'Sample 1')
        self.assertAlmostEqual(mean(sample_d['Exp1']['Sample 2']), 0.5)

    def test_censor_background(self):
        addme =  [['Exp2', 'Sample 1', 10.0], ['Exp2', 'Sample 2', 5.0]]
        sample_list = self.sample_list + addme
        sample_list = sample_list + [['Exp2', 'NTC', 11.0]]
        sample_frame = pd.DataFrame(sample_list, columns=['Target', 'Sample', 'Cq'])
        censored = censor_frame_background(sample_frame, ntc_samples=['NTC'], margin=log2(10))

    def test_silences_near_ntc(self):
        addme =  [['Exp2', 'Sample 1', 10.0], ['Exp2', 'Sample 2', 5.0]]
        sample_list = self.sample_list + addme
        sample_frame = pd.DataFrame(sample_list, columns=['Target', 'Sample', 'Cq'])
        sample_d, ignored_genes, ignored_samples = expression(sample_frame, 'Ref1', 'Sample 1')
        self.assertAlmostEqual(mean(sample_d['Exp2']['Sample 2']), 2**5)
        
        sample_list = sample_list + [['Exp2', 'NTC', 11.0]]
        sample_frame = pd.DataFrame(sample_list, columns=['Target', 'Sample', 'Cq'])
        sample_d, ignored_genes, ignored_samples = expression(sample_frame, 'Ref1', 'Sample 1')
        self.assertIn('Exp2', ignored_genes)

class TestAverage(unittest.TestCase):
    def test_averaging(self):
        self.assertAlmostEqual(average_cq((3.0,5.0)), 3.678071905112638)
        self.assertAlmostEqual(average_cq((3.0,4.0,5.0)), 3.7776075786635523)

class TestRankGenes(unittest.TestCase):
    def setUp(self):
        with open('vandesompele-2002-cq.txt') as f:
            buf = f.readlines()
        self.gene_names = buf[0][:-1].split('\t')[2:]
        # slicing to cut off newlines on each line and skip header row
        # yields ['Tissue', 'Sample', Cq0, Cq1, Cq2...]
        buf = [line[:-1].split('\t') for line in buf[1:] if len(line) > 2]
        self.d = {}
        for line in buf:
            for (index, name) in enumerate(self.gene_names):
                self.d.setdefault(line[0], []).append([name, line[1], float(line[2+index])])
        for key in self.d:
            self.d[key] = pd.DataFrame(self.d[key], columns=['Target', 'Sample', 'Cq'])
        # e.g. d['Fib'] = [['ACTB', 'Fib1', 12.34], ...]

    def test_ranking(self):
        key = {'Neuroblastoma': ['HPRT1', 'GAPD', 'SDHA', 'UBC', 'HMBS', 'YWHAZ', 'TBP', 'ACTB', 'RPL13A', 'B2M'],
           'Fib': ['HPRT1', 'GAPD', 'YWHAZ', 'UBC', 'ACTB', 'TBP', 'SDHA', 'RPL13A', 'B2M', 'HMBS'],
           'Leukocyte': ['UBC', 'YWHAZ', 'B2M', 'GAPD', 'RPL13A', 'TBP', 'SDHA', 'HPRT1', 'HMBS', 'ACTB'],
           'BM': ['UBC', 'RPL13A', 'YWHAZ', 'HPRT1', 'GAPD', 'SDHA', 'TBP', 'HMBS', 'B2M', 'ACTB'],
           'Pool': ['SDHA', 'GAPD', 'HMBS', 'HPRT1', 'TBP', 'UBC', 'RPL13A', 'YWHAZ', 'ACTB', 'B2M']}

        for tissue in self.d:
            ranked = rank_genes(self.d[tissue], self.gene_names, self.d[tissue].ix[0, 'Sample'])
            # no preference between the two best genes
            self.assertIn(ranked[0], key[tissue][:2])
            self.assertIn(ranked[1], key[tissue][:2])
            self.assertEqual(ranked[2:], key[tissue][2:])

    def test_frame_ranking(self):
        key = {'Neuroblastoma': ['HPRT1', 'GAPD', 'SDHA', 'UBC', 'HMBS', 'YWHAZ', 'TBP', 'ACTB', 'RPL13A', 'B2M'],
           'Fib': ['HPRT1', 'GAPD', 'YWHAZ', 'UBC', 'ACTB', 'TBP', 'SDHA', 'RPL13A', 'B2M', 'HMBS'],
           'Leukocyte': ['UBC', 'YWHAZ', 'B2M', 'GAPD', 'RPL13A', 'TBP', 'SDHA', 'HPRT1', 'HMBS', 'ACTB'],
           'BM': ['UBC', 'RPL13A', 'YWHAZ', 'HPRT1', 'GAPD', 'SDHA', 'TBP', 'HMBS', 'B2M', 'ACTB'],
           'Pool': ['SDHA', 'GAPD', 'HMBS', 'HPRT1', 'TBP', 'UBC', 'RPL13A', 'YWHAZ', 'ACTB', 'B2M']}

        for tissue in self.d:
            ranked = rank_target_frame(self.d[tissue], self.gene_names, self.d[tissue].ix[0, 'Sample'])
            # no preference between the two best genes
            self.assertIn(ranked.ix[0, 'Target'], key[tissue][:2])
            self.assertIn(ranked.ix[1, 'Target'], key[tissue][:2])
            self.assertEqual(list(ranked.ix[2:, 'Target']), key[tissue][2:])

    def test_nfs(self):
        comparison_vs = {'Fib': (3, 0.109409), 'BM': (2, 0.113745), 'Leukocyte': (3, 0.116988),
                'Pool': (2, 0.202805), 'Neuroblastoma': (4, 0.138337)}
        for tissue in self.d:
            ranked = rank_genes(self.d[tissue], self.gene_names, self.d[tissue].ix[0, 'Sample'])
            nfs = calculate_all_nfs(self.d[tissue], ranked, self.d[tissue].ix[0, 'Sample'])
            vs = nf_v(nfs)
            test_v, test_value = comparison_vs[tissue]
            self.assertAlmostEqual(vs[test_v], test_value, places=5)
            # for v in sorted(vs.keys()):
                # if v > 3 and vs[v-1] < 0.15: break
                # print '%d: %f' % (v, vs[v]),
                # stdout.flush()

    def test_nfs_frame(self):
        comparison_vs = {'Fib': (3, 0.109409), 'BM': (2, 0.113745), 'Leukocyte': (3, 0.116988),
                'Pool': (2, 0.202805), 'Neuroblastoma': (4, 0.138337)}
        for tissue in self.d:
            print '\n%s' % tissue
            ranked = rank_genes(self.d[tissue], self.gene_names, self.d[tissue].ix[0, 'Sample'])
            nfs = calculate_all_nfs_frame(self.d[tissue], ranked, self.d[tissue].ix[0, 'Sample'])
            vs = nf_v_frame(nfs)
            test_v, test_value = comparison_vs[tissue]
            self.assertAlmostEqual(vs[test_v], test_value, places=5)

class TestExPd(unittest.TestCase):
    def test_validation(self):
        # empty frame
        a = pd.DataFrame()
        self.assertRaises(TypeError, validate_sample_frame, ['a', 'b', 'c'])
        self.assertRaises(TypeError, validate_sample_frame, None)
        self.assertRaises(ValueError, validate_sample_frame, a)

        # valid minimal frame
        b = pd.DataFrame({'foo': [0], 'Cq': [1.0], 'Sample': ['bar'], 'Target': ['baz']})
        self.assertIs(validate_sample_frame(b), True)
        
        # Integer Cq column should fail
        c = pd.DataFrame({'foo': [0], 'Cq': [1], 'Sample': ['bar'], 'Target': ['baz']})
        self.assertRaises(ValueError, validate_sample_frame, c)

        # Invalid sample types (Content field) should fail
        # d = pd.DataFrame({'foo': [0, 1], 'Cq': [1.0, 2.0], 'Sample': ['bar', 'bar'], 'Target': ['baz', 'blip'], 'Content': ['Unknown', 'NTC']})
        # e = pd.DataFrame({'foo': [0, 1], 'Cq': [1.0, 2.0], 'Sample': ['bar', 'bar'], 'Target': ['baz', 'blip'], 'Content': ['Unknown', 'blah']})
        # self.assertIs(validate_sample_frame(d), True)
        # self.assertRaises(ValueError, validate_sample_frame, e)


if __name__ == '__main__':
    unittest.main()

