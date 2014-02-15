#!/usr/bin/env python

from math import isnan
from numpy import mean, std, power, asarray, log
from scipy.stats.mstats import gmean
from warnings import warn
from types import *
from itertools import repeat
import pandas as pd

log2 = lambda x: log(x)/log(2)

def average_cq(seq, efficiency=1.0):
    # given a set of Cq values, return the Cq value that represents the
    # average expression level of the input.
    # The intent is to average the expression levels of the samples,
    # since the average of Cq values is not biologically meaningful.
    denominator = sum( [pow(2.0*efficiency, -Ci) for Ci in seq] )
    return log(len(seq)/denominator)/log(2.0*efficiency)

def make_sample_dict(sample_frame):
    d = {}
    g = sample_frame.groupby(['Target', 'Sample'])
    for (target, sample), h in g:
        d.setdefault(target, {})[sample] = list(h['Cq'])
    return d

def validate_sample_frame(sample_frame):
    if not isinstance(sample_frame, pd.core.frame.DataFrame):
        raise TypeError("Expected a pandas DataFrame, received {}".format(type(sample_frame)))
    for col in ['Sample', 'Target', 'Cq']:
        if col not in sample_frame:
            raise ValueError("Missing column {} in sample frame".format(col))
    if sample_frame['Cq'].dtype.kind != 'f':
        raise ValueError("Expected Cq column to have float type; has type {} instead".format(str(sample_frame['Cq'].dtype)))
    return True

def censor_background(d):
    # censor anything too close to any NTCs we have and remove NTCs from d
    # operates on d in place
    margin = log2(10)
    for target in d:
        if not ('NTC' in d[target]): continue
        too_close = min(d[target]['NTC']) - margin
        for sample in d[target].keys():
            if sample == 'NTC': continue
            d[target][sample] = [i for i in d[target][sample] if i < too_close]
            if len(d[target][sample]) == 0: del d[target][sample]
        del d[target]['NTC']

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

def censor_targets_missing_refsample(d, ref_sample):
    # all genes must have a reference sample
    ignored_targets = []
    for target in d.keys():
        if ref_sample not in d[target]:
            ignored_targets.append(target)
            del d[target]
    return ignored_targets

def expression(sample_list, ref_gene, ref_sample):
    # out, ignored_targets, ignored_samples = expression(sample_list, ref_gene, ref_sample)
    # sample_list has form [['target1', 'sample1', Cq1], ...]
    # replicates for experimental genes are considered independently
    # replicates for reference genes are averaged for referencing
    # returns (sample_dictionary, [ignored genes], [ignored samples])
    # sample_dictionary has the form s_d[target][sample] = [rel0, rel1] where
    # each rel_n is from a replicate in the input list and is relative to the
    # reference sample

    d = make_sample_dict(sample_list)
    censor_background(d)
    ignored_targets = censor_targets_missing_refsample(d, ref_sample)

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
    
def expression_nf(sample_list, nfs, ref_sample):
    d = make_sample_dict(sample_list)
    censor_background(d)
    ignored_targets = censor_targets_missing_refsample(d, ref_sample)

    # make sure all samples have a NF
    ignored_samples = []
    for target in d:
        for sample in d[target]:
            if sample not in nfs:
                ignored_samples.append(sample)
    for target in d:
        for ignoreme in ignored_samples:
            if ignoreme in d['target']:
                del d[target][ignoreme]

    out = {}
    for target in d:
        for sample in d[target]:
            for cq in d[target][sample]:
                delta = -cq + average_cq(d[target][ref_sample])
                rel = pow(2, delta) / nfs[sample]
                out.setdefault(target, {}).setdefault(sample, []).append(rel)

    return (out, ignored_targets, ignored_samples)


def collect_expression(sample_frame, ref_genes, ref_sample):
    # normalize all samples by each of the reference genes
    table = {}
    if not (type(ref_genes) is set): ref_genes = set(ref_genes)
    available_samples = set(sample_frame['Sample']) - {'NTC'}
    for ref_gene in ref_genes:
        table[ref_gene], ignored_targets, ignored_samples = expression(sample_frame, ref_gene, ref_sample)
        if set(ignored_targets).intersection(ref_genes):
            raise ValueError("Bad reference sample. One or more reference gene values are not available.")
        available_samples -= set(ignored_samples)
    return table, available_samples

def collect_expression_frame(sample_frame, ref_targets, ref_sample):
    by_gene = {'Sample': sample_frame['Sample'], 'Target': sample_frame['Target']}
    for target in ref_targets:
        by_gene[target] = expression_frame(sample_frame, target, ref_sample)
    return pd.DataFrame(by_gene)

def rank_genes(sample_list, ref_genes, ref_sample, with_m=False):
    # applies the Vandesompele et al. (2002) method to rank reference genes
    # in order of stability.
    # doi:10.1186/gb-2002-3-7-research0034

    ref_genes = set(ref_genes)

    table, available_samples = collect_expression(sample_list, ref_genes, ref_sample)

    worst = []
    worst_m = []
    while len(ref_genes) - len(worst) > 1:
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
        worst_m.append(max(M)[0])
    best = ref_genes - set(worst) # whatever's left over
    worst.reverse()
    worst_m.reverse()
    worst_m = [worst_m[0]] + worst_m
    if with_m:
        return list(best) + worst, worst_m
    else:
        return list(best) + worst

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

def calculate_all_nfs(sample_list, ranked_genes, ref_sample):
    # Given a list of sample results and an ordered list of reference genes
    # (as produced by rank_genes), calculate_all_nfs computes normalization factors
    # for all samples, using successively larger subsets of the reference gene
    # list. All computed normalization factors are returned. Use nf_v() to
    # determine how many reference genes to include in your experimental
    # samples.
    # Returns a dictionary nfs, where nfs[n_genes][sample] = nf, for
    # 2 <= n_genes <= len(ref_genes).

    # d[ref_gene][target][sample] = [rel0, rel1]
    nfs = {}
    d = make_sample_dict(sample_list)
    for i in range(1, len(ranked_genes)+1):
        nfs[i] = {}
        for sample in d[ranked_genes[0]]:
            nfs[i][sample] = gmean([pow(2, -average_cq(d[ref_gene][sample]) + average_cq(d[ref_gene][ref_sample])) for ref_gene in ranked_genes[:i]])
    return nfs

def calculate_all_nfs_frame(sdf, ranked_genes, ref_sample):
    # Returns a DataFrame, where rows represent samples and columns represent a number of reference genes.
    grouped = sdf.groupby(['Target', 'Sample'])['Cq'].aggregate(average_cq)
    samples = sdf['Sample'].unique()
    nfs = {}
    for i in xrange(1, len(ranked_genes)+1):
        nfs[i] = gmean([pow(2, -grouped.ix[zip(repeat(ref_gene), samples)] + grouped.ix[ref_gene, ref_sample]) for ref_gene in ranked_genes[:i]])
    return pd.DataFrame(nfs, index=samples)

def nfs_from(sample_list, ref_genes, ref_sample):
    # Returns a dictionary nfs, where nfs[sample] = nf.
    nfs = {}
    d = make_sample_dict(sample_list) # d[ref_gene][target][sample] = [rel0, rel1]
    for sample in d[ref_genes[0]]:
        if sample == 'NTC': continue
        nfs[sample] = gmean([pow(2, -average_cq(d[ref_gene][sample]) + average_cq(d[ref_gene][ref_sample])) for ref_gene in ref_genes])
    return nfs

def nf_v(nfs):
    v = {}
    for i in range(2, max(nfs.keys())+1):
        if not (i in nfs): break
        if not (i+1 in nfs): break
        v[i] = std([ log2(nfs[i][sample]/nfs[i+1][sample]) for sample in nfs[i] ], ddof=1)
    return v

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
    for v in sorted(vs.keys()):
        if v > 3 and vs[v-1] < 0.15: break
        rec.append((ranked_genes[v-1], vs[v]))
    return rec

