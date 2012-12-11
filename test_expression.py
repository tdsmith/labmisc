#!/usr/bin/env python

import unittest
from expression import expression, average_cq
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

if __name__ == '__main__':
    unittest.main()

