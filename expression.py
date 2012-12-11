#!/usr/bin/env python

from math import isnan, log
from types import *

def average_cq(seq, efficiency=1.0):
    denominator = sum( [pow(2.0*efficiency, -Ci) for Ci in seq] )
    return log(len(seq)/denominator)/log(2.0*efficiency)

def expression(sample_list, ref_gene, ref_sample):
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
