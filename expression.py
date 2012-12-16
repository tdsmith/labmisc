#!/usr/bin/env python

from math import isnan, log
from numpy import mean, std
from types import *

log2 = lambda x: log(x)/log(2)

def average_cq(seq, efficiency=1.0):
    denominator = sum( [pow(2.0*efficiency, -Ci) for Ci in seq] )
    return log(len(seq)/denominator)/log(2.0*efficiency)

def expression(sample_list, ref_gene, ref_sample):
    # out, ignored_targets, ignored_samples = expression(sample_list, ref_gene, ref_sample)
    # sample_list has form [['target1', 'sample1', Cq1], ...]
    # replicates for experimental genes are considered independently
    # replicates for reference genes are averaged for referencing
    # returns (sample_dictionary, [ignored genes], [ignored samples])
    # sample_dictionary has the form s_d[target][sample] = [rel0, rel1] where
    # each rel_n is from a replicate in the input list and is relative to the
    # reference sample

    d = {} # d[target][sample] = [cq1, cq2, ...]
    for (target, sample, cq) in sample_list:
        assert type(cq) is FloatType
        if(isnan(cq)): continue # drop silently
        d.setdefault(target,{}).setdefault(sample,[]).append(cq)

    # censor anything too close to any NTCs we have
    margin = log(10)/log(2)
    for target in d:
        if not ('NTC' in d[target]): continue
        too_close = min(d[target]['NTC']) - margin
        for sample in d[target].keys():
            if sample == 'NTC': continue
            d[target][sample] = [i for i in d[target][sample] if i < too_close]
            if len(d[target][sample]) == 0: del d[target][sample]
        del d[target]['NTC']

    # make sure all genes have a reference sample
    ignored_targets = []
    for target in d.keys():
        if not (ref_sample in d[target]):
            ignored_targets.append(target)
            del d[target]

    # reference genes need to be defined for each sample. ID them first
    ignored_samples = []
    for target in d:
        for sample in d[target]:
            if not(sample in d[ref_gene]):
                ignored_samples.append(sample)

    # and now clean up
    for target in d:
        for ignoreme in ignored_samples:
            if ignoreme in d[target]:
                del d[target][ignoreme]
    
    out = {}
    # perform 2**delta-delta Cq expression comparison
    for target in d:
        for sample in d[target]:
            for cq in d[target][sample]:
                delta_ref = average_cq(d[ref_gene][sample]) - average_cq(d[ref_gene][ref_sample])
                delta_target = cq - average_cq(d[target][ref_sample])
                rel = pow(2, delta_ref - delta_target)
                out.setdefault(target, {}).setdefault(sample, []).append(rel)

    return (out, ignored_targets, ignored_samples)

def rank_genes(sample_list, ref_genes, ref_sample):
    # applies the Vandesompele et al. (2002) method to rank reference genes
    # in order of stability.
    # doi:10.1186/gb-2002-3-7-research0034

    ref_genes = set(ref_genes)

    # first, normalize all samples by each of the reference genes
    table = {}
    available_samples = set([sample[1] for sample in sample_list])
    for ref_gene in ref_genes:
        table[ref_gene], ignored_targets, ignored_samples = expression(sample_list, ref_gene, ref_sample)
        if set(ignored_targets).intersection(ref_genes):
            raise ValueError("Bad reference sample. One or more reference gene values are not available.")
        available_samples -= set(ignored_samples)

    worst = []
    while len(ref_genes) - len(worst) > 2:
        # calculate M for each gene
        M = []
        for test_gene in ref_genes:
            if test_gene in worst: continue
            Vs = []
            for ref_gene in ref_genes:
                if ref_gene == test_gene: continue
                if ref_gene in worst: continue
                # compute A_(ref_gene, test_gene)
                A = [log2(mean(table[ref_gene][test_gene][sample])) for sample in available_samples]
                Vs.append(std(A))
            M.append( (sum(Vs) / (len(ref_genes) - len(worst) - 1), test_gene) )
        worst.append(max(M)[1])
    best = ref_genes - set(worst) # whatever's left over
    worst.reverse()
    return list(best) + worst
