#!/usr/bin/env python

from math import isnan
from numpy import mean, std, power, asarray, log
from scipy.stats.mstats import gmean
from warnings import warn
from types import *
from itertools import repeat
import pandas as pd

log2 = lambda x: log(x)/log(2)

def make_sample_dict(sample_frame):
    d = {}
    g = sample_frame.groupby(['Target', 'Sample'])
    for (target, sample), h in g:
        d.setdefault(target, {})[sample] = list(h['Cq'])
    return d

def average_cq(seq, efficiency=1.0):
    # given a set of Cq values, return the Cq value that represents the
    # average expression level of the input.
    # The intent is to average the expression levels of the samples,
    # since the average of Cq values is not biologically meaningful.
    denominator = sum( [pow(2.0*efficiency, -Ci) for Ci in seq] )
    return log(len(seq)/denominator)/log(2.0*efficiency)

def validate_sample_frame(sample_frame):
    if not isinstance(sample_frame, pd.core.frame.DataFrame):
        raise TypeError("Expected a pandas DataFrame, received {}".format(type(sample_frame)))
    for col in ['Sample', 'Target', 'Cq']:
        if col not in sample_frame:
            raise ValueError("Missing column {} in sample frame".format(col))
    if sample_frame['Cq'].dtype.kind != 'f':
        raise ValueError("Expected Cq column to have float type; has type {} instead".format(str(sample_frame['Cq'].dtype)))
    return True

def censor_frame_background(sdf, ntc_samples, margin):
    ntcs = sdf.loc[ sdf['Sample'].apply(lambda x: x in ntc_samples), ]
    if ntcs.empty:
        return sdf
    g = ntcs.groupby('Target')
    min_ntcs = g['Cq'].min()
    # if a target has no NTC, min_ntcs.loc[sample] is NaN
    # we should retain all values from targets with no NTC
    # all comparisons with NaN are false
    # so we test for the "wrong" condition and invert the result
    censored = sdf.loc[ ~(sdf['Cq'] > (min_ntcs.loc[sdf['Target']] - margin)) ]
    return censored

def expression_frame(sdf, ref_target, ref_sample, ntc_samples=['NTC'], ntc_margin=log2(10)):
    # It might be more correct to replace asarray calls (to discard indexes)
    # with proper joins.
    censored = censor_frame_background(sdf, ntc_samples, ntc_margin)
  
    ref_target_df = censored.ix[censored['Target'] == ref_target, ['Sample', 'Cq']]
    ref_target_grouped = ref_target_df.groupby('Sample')
    ref_target_mean_by_sample = ref_target_grouped['Cq'].aggregate(average_cq)
    ref_target_mean_list = ref_target_mean_by_sample.ix[censored['Sample']]
    ref_target_delta = asarray(ref_target_mean_list - ref_target_mean_by_sample[ref_sample])

    ref_sample_df = censored.ix[censored['Sample'] == ref_sample, ['Target', 'Cq']]
    ref_sample_grouped = ref_sample_df.groupby('Target')
    ref_sample_mean_by_target = ref_sample_grouped['Cq'].aggregate(average_cq)
    ref_sample_delta = asarray(censored['Cq'] - asarray(ref_sample_mean_by_target.ix[censored['Target']]))

    rel_exp = pd.Series(
            power(2, ref_target_delta - ref_sample_delta),
            index = censored.index)

    return rel_exp
    
def expression_nf_frame(sample_frame, nf_n, ref_sample, ntc_samples=['NTC'], ntc_margin=log2(10)):
    censored = censor_frame_background(sample_frame, ntc_samples, ntc_margin)
    ref_sample_df = censored.ix[censored['Sample'] == ref_sample, ['Target', 'Cq']]
    ref_sample_cq = ref_sample_df.groupby('Target')['Cq'].aggregate(average_cq)

    delta = -censored['Cq'] + asarray(ref_sample_cq.ix[censored['Target']])
    rel = power(2, delta) / asarray(nf_n.ix[censored['Sample']])
    return rel

def collect_expression_frame(sample_frame, ref_targets, ref_sample):
    by_gene = {'Sample': sample_frame['Sample'], 'Target': sample_frame['Target']}
    for target in ref_targets:
        by_gene[target] = expression_frame(sample_frame, target, ref_sample)
    return pd.DataFrame(by_gene)

def rank_target_frame(sample_frame, ref_targets, ref_sample):
    table = collect_expression_frame(sample_frame, ref_targets, ref_sample)
    all_samples = sample_frame['Sample'].unique()
    t = table.groupby(['Sample', 'Target']).mean()
    logt = log2(t)
    ref_targets = set(ref_targets)

    worst = []
    worst_m = []
    while len(ref_targets) - len(worst) > 1:
        M = []
        for test_target in ref_targets:
            if test_target in worst: continue
            Vs = []
            for ref_target in ref_targets:
                if ref_target == test_target or ref_target in worst: continue
                A = logt.ix[zip(all_samples, repeat(test_target)), ref_target]
                Vs.append(A.std())
            M.append( (sum(Vs)/(len(ref_targets)-len(worst)-1), test_target) )
        worst.append(max(M)[1])
        worst_m.append(max(M)[0])
    best = ref_targets - set(worst)
    worst.reverse()
    worst_m.reverse()
    worst_m = [worst_m[0]] + worst_m
    return pd.DataFrame({'Target': list(best) + worst, 'M': worst_m}, columns=['Target', 'M'])

def calculate_all_nfs_frame(sdf, ranked_genes, ref_sample):
    # Returns a DataFrame, where rows represent samples and columns represent a number of reference genes.
    grouped = sdf.groupby(['Target', 'Sample'])['Cq'].aggregate(average_cq)
    samples = sdf['Sample'].unique()
    nfs = {}
    for i in xrange(1, len(ranked_genes)+1):
        nfs[i] = gmean([pow(2, -grouped.ix[zip(repeat(ref_gene), samples)] + grouped.ix[ref_gene, ref_sample]) for ref_gene in ranked_genes[:i]])
    return pd.DataFrame(nfs, index=samples)

"""
def nfs_from(sample_list, ref_genes, ref_sample):
    # Returns a dictionary nfs, where nfs[sample] = nf.
    nfs = {}
    d = make_sample_dict(sample_list) # d[ref_gene][target][sample] = [rel0, rel1]
    for sample in d[ref_genes[0]]:
        if sample == 'NTC': continue
        nfs[sample] = gmean([pow(2, -average_cq(d[ref_gene][sample]) + average_cq(d[ref_gene][ref_sample])) for ref_gene in ref_genes])
    return nfs
"""

def nf_v_frame(nfs):
    v = []
    if (nfs.columns != range(1, nfs.columns[-1]+1)).any():
        raise ValueError("Column names invalid in nf_v_frame")
    for i in nfs.columns[:-1]:
        v.append(std(log2(nfs[i]/nfs[i+1]), ddof=1))
    return pd.Series(v, index=nfs.columns[:-1])

def recommend_refset(sample_list, ref_genes, ref_sample):
    ranked_genes = rank_genes(sample_list, ref_genes, ref_sample)
    nfs = calculate_all_nfs(sample_list, ref_genes, ref_sample)
    vs = nf_v(nfs)
    rec = [(ranked_genes[0], 0)]
    for v in sorted(vs.index):
        if v > 3 and vs[v-1] < 0.15: break
        rec.append((ranked_genes[v-1], vs[v]))
    return rec

