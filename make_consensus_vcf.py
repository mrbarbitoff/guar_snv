#!/usr/bin/env python

import sys
import re
import os

gatk, ngsep, tassel = sys.argv[1:]

ngsep_vars = {}
tassel_vars = {}

with open(ngsep, 'r') as ng_handle:
    for line in ng_handle:
        if line.startswith('#'):
            continue
        content = line.strip().split('\t')
        if ',' in content[4]:
            continue
        ngsep_nhet = len(re.findall('\t0/1:', line))
        ngsep_nalt = len(re.findall('\t[01]/1:', line))
        ngsep_af = re.findall('AF=([^;^\t]+)', line)[0]
        var_id = f'{content[0]}:{content[1]}'
        ngsep_vars[var_id] = (ngsep_af, ngsep_nhet, ngsep_nalt)

with open(tassel, 'r') as ts_handle:
    for line in ts_handle:
        if line.startswith('#'):
            continue
        content = line.strip().split('\t')
        if ',' in content[4]:
            continue
        ts_nhet = len(re.findall('\t0/1:', line))
        ts_nalt = len(re.findall('\t[01]/1:', line))
        ts_af = (ts_nhet + 2 * ts_nalt)/len(re.findall('\t\d/\d:', line))
        var_id = f'{content[0]}:{content[1]}'.lower()
        tassel_vars[var_id] = (ts_af, ts_nhet, ts_nalt)


with open(gatk, 'r') as ga_handle:
    for line in ga_handle:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                print('##INFO=<ID=NGSEP_NHET,Number=A,Type=Integer,Description="HET genotype count (NGSEP)">')
                print('##INFO=<ID=NGSEP_NALT,Number=A,Type=Integer,Description="Non-reference genotype count (NGSEP)">')
                print('##INFO=<ID=NGSEP_AF,Number=A,Type=Float,Description="Allele Frequency (NGSEP)">')
                print('##INFO=<ID=TASSEL_NHET,Number=A,Type=Integer,Description="HET genotype count (TASSEL)">')
                print('##INFO=<ID=TASSEL_NALT,Number=A,Type=Integer,Description="Non-reference genotype count (TASSEL)">')
                print('##INFO=<ID=TASSEL_AF,Number=A,Type=Float,Description="Allele Frequency (TASSEL)">')
            print(line.strip())
            continue
        content = line.strip().split('\t')
        var_id = f'{content[0]}:{content[1]}'
        if var_id not in ngsep_vars and var_id not in tassel_vars:
#            print(var_id)
            continue
        n1, n2, n3 = ngsep_vars[var_id] if var_id in ngsep_vars else (0, 0, 0)
        t1, t2, t3 = tassel_vars[var_id] if var_id in tassel_vars else (0, 0, 0)
        content[7] = f'{content[7]};NGSEP_NHET={n2};NGSEP_NALT={n3};NGSEP_AF={n1};TASSEL_NHET={t2};TASSEL_NALT={t3};TASSEL_AF={t1}'
        print('\t'.join(content))
