#!/usr/bin/env python

import sys
import re
import os

gatk, ngsep, tassel = sys.argv[1:]

ngsep_vars = {}
tassel_vars = {}

def avg_ab(tup):
    if sum(tup) > 0:
        return str(tup[1]/sum(tup))
    else:
        return 'NA'


def sum_dp(tup, an):
    if sum(tup) == 0:
        return 'NA'
    else:
        return str(sum(tup)/an)


def get_ab_dp(content, ab_index):
    ab_nonref = [0, 0]
    ab_ref = [0, 0]
    ab_homvar = [0, 0]
    counts = [0, 0, 0]
    try:
        for samp in content[9:]:
            if samp.startswith('0/0'):
                ab_ref[0] += [int(x) for x in samp.split(':')[ab_index].split(',')][0]
                ab_ref[1] += [int(x) for x in samp.split(':')[ab_index].split(',')][1]
                counts[0] += 1
            elif samp.startswith('0/1'):
                ab_nonref[0] += [int(x) for x in samp.split(':')[ab_index].split(',')][0]
                ab_nonref[1] += [int(x) for x in samp.split(':')[ab_index].split(',')][1]
                counts[1] += 1
            elif samp.startswith('1/1'):
                ab_homvar[0] += [int(x) for x in samp.split(':')[ab_index].split(',')][0]
                ab_homvar[1] += [int(x) for x in samp.split(':')[ab_index].split(',')][1]
                counts[2] += 1
        ab_stats = '\t'.join([avg_ab(x) for x in [ab_ref, ab_nonref, ab_homvar]])
        dp_stats = '\t'.join([sum_dp(x, counts[i]) for i, x in enumerate([ab_ref, ab_nonref, ab_homvar])])
        return f'{ab_stats}\t{dp_stats}'
    except:
        return 'NA\tNA\tNA\tNA\tNA\tNA'


def get_ab_dp_bsdp(content, ab_index):
    ab_nonref = [0, 0]
    ab_ref = [0, 0]
    ab_homvar = [0, 0]
    counts = [0, 0, 0]
    nuc_numbers = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    if content[4] not in nuc_numbers or content[3] not in nuc_numbers:
        return 'NA\tNA\tNA\tNA\tNA\tNA'
    alleles = [nuc_numbers[content[3]], nuc_numbers[content[4]]]
    try:
        for samp in content[9:]:
            if samp.startswith('0/0'):
                ab_ref[0] += [int(x) for i, x in enumerate(samp.split(':')[ab_index].split(',')) if i in alleles][0]
                ab_ref[1] += [int(x) for i, x in enumerate(samp.split(':')[ab_index].split(',')) if i in alleles][1]
                counts[0] += 1
            elif samp.startswith('0/1'):
                ab_nonref[0] += [int(x) for i, x in enumerate(samp.split(':')[ab_index].split(',')) if i in alleles][0]
                ab_nonref[1] += [int(x) for i, x in enumerate(samp.split(':')[ab_index].split(',')) if i in alleles][1]
                counts[1] += 1
            elif samp.startswith('1/1'):
                ab_homvar[0] += [int(x) for i, x in enumerate(samp.split(':')[ab_index].split(',')) if i in alleles][0]
                ab_homvar[1] += [int(x) for i, x in enumerate(samp.split(':')[ab_index].split(',')) if i in alleles][1]
                counts[2] += 1
        ab_stats = '\t'.join([avg_ab(x) for x in [ab_ref, ab_nonref, ab_homvar]])
        dp_stats = '\t'.join([sum_dp(x, counts[i]) for i, x in enumerate([ab_ref, ab_nonref, ab_homvar])])
        return f'{ab_stats}\t{dp_stats}'
    except:
        return 'NA\tNA\tNA\tNA\tNA\tNA'



with open(ngsep, 'r') as ng_handle:
    for line in ng_handle:
        if line.startswith('#'):
            continue
        content = line.strip().split('\t')
        if ',' in content[4]:
            continue
        ngsep_nhet = len(re.findall('\t0/1:', line))
        ngsep_nalt = len(re.findall('\t[01]/1:', line))
        ngsep_af = (ngsep_nhet + 2 * (ngsep_nalt - ngsep_nhet))/(2 * len(re.findall('\t\d/\d:', line)))
        var_id = f'{content[0]}:{content[1]}'
        stats = get_ab_dp_bsdp(content, 4)
        print(f'{var_id}\tNGSEP\t{ngsep_af}\t{ngsep_nhet}\t{ngsep_nalt}\t{stats}')

with open(tassel, 'r') as ts_handle:
    for line in ts_handle:
        if line.startswith('#'):
            continue
        content = line.strip().split('\t')
        if ',' in content[4]:
            continue
        ts_nhet = len(re.findall('\t0/1:', line))
        ts_nalt = len(re.findall('\t[01]/1:', line))
        ts_af = (ts_nhet + 2 * (ts_nalt - ts_nhet))/(2 * len(re.findall('\t\d/\d:', line)))
        var_id = f'{content[0]}:{content[1]}'.lower()
        stats = get_ab_dp(content, 1)
        print(f'{var_id}\tTASSEL\t{ts_af}\t{ts_nhet}\t{ts_nalt}\t{stats}')


with open(gatk, 'r') as ga_handle:
    for line in ga_handle:
        if line.startswith('#'):
            continue
        content = line.strip().split('\t')
        if ',' in content[4]:
            continue
        gatk_nhet = len(re.findall('\t0/1:', line))
        gatk_nalt = len(re.findall('\t[01]/1:', line))
        try:
            gatk_af = (gatk_nhet + 2 * (gatk_nalt - gatk_nhet))/(2 * len(re.findall('\t\d/\d:', line)))
        except:
            gatk_af = 'NA'
        var_id = f'{content[0]}:{content[1]}'.lower()
        stats = get_ab_dp(content, 1)
        print(f'{var_id}\tGATK\t{gatk_af}\t{gatk_nhet}\t{gatk_nalt}\t{stats}')
