#!/usr/bin/env python3

import sys
import os
import re
import pandas as pd

fasta, vcf, transcripts_match, gwas_folder, gm_file = sys.argv[1:]

variant_map = {}
with open(gm_file, 'r') as gm_handle:
    for line in gm_handle:
        if line.startswith('SNP'):
            continue
        content = line.strip().split('\t')
        variant_map[content[0]] = f'{content[1]}:{content[2]}'

gwas_results = {}
for tsv_file in [str(f) for f in os.listdir(gwas_folder) if os.path.isfile(os.path.join(gwas_folder, f))]:
    trait_id = tsv_file.replace('.tsv', '') 
    gwas_results[trait_id] = {}
    print(f'Reading {trait_id}')
    with open(f'{gwas_folder}/{tsv_file}', 'r') as tsv_handle:
        for line in tsv_handle:
            if not line.startswith('GUARVAR'):
                continue
            content = line.strip().split('\t')
            variant_id = variant_map[content[0]]
            gwas_results[trait_id][variant_id] = format(float(content[3]), 'f')


in_transcripts = {}
with open(transcripts_match, 'r') as tm_file:
    for line in tm_file:
        content = line.strip().split('\t')
        var_id = ':'.join(content[:5])
        if var_id not in in_transcripts:
            in_transcripts[var_id] = []
        in_transcripts[var_id].append(content[8])


sequences = {}
with open(fasta, 'r') as f_handle:
    for line in f_handle:
        if line.startswith('>'):
            seqname = line.strip().split('-')[0][1:]
        else:
            sequences[seqname] = (line.strip()[:150], line.strip()[151:])

data_list = {'Contig': [],
        'Position': [],
        'Upstream sequence': [],
        'Reference': [],
        'Alt': [],
        'Downstream sequence': [],
        'Transcripts': [],
        'AF': [],
        'Average depth': [],
        '# HET': [],
        'HET alt allele %': []}
for gwas_trait in gwas_results:
    data_list[gwas_trait] = []
names = []
with open(vcf, 'r') as vcf_handle:
    for line in vcf_handle:
        if line.startswith('#'):
            if not line.startswith('#CHROM'):
                continue
            content = line.strip().split('\t')
            for sth in content[9:]:
                data_list[sth] = []
                names.append(sth)
            continue
        content = line.strip().split('\t')
#        if len(content[4]) > 1 or len(content[3]) > 1:
#            continue
#        if len(re.findall('\t[01]/[01]:', line)) < (0.95 * (len(content) - 9)):
#            continue
        seq_id = f'{content[0]}:{int(content[1]) - 151}'
        if seq_id not in sequences:
            continue
        data_list['Contig'].append(content[0])
        data_list['Position'].append(content[1])
        data_list['Reference'].append(content[3])
        data_list['Alt'].append(content[4])
        data_list['Upstream sequence'].append(sequences[seq_id][0])
        data_list['Downstream sequence'].append(sequences[seq_id][1])
        if ':'.join(content[:5]) in in_transcripts:
            data_list['Transcripts'].append(';'.join(in_transcripts[':'.join(content[:5])]))
        else:
            data_list['Transcripts'].append('none')
        af_value = re.findall('AF=([^;]+)', content[7])[0]
        an_value = int(re.findall('AN=(\d+)', content[7])[0])
        dp_value = int(re.findall('DP=([^;^\t]+)', content[7])[0])
        data_list['AF'].append(af_value)
        data_list['Average depth'].append(dp_value/an_value * 2)
        where_ad = [i for i, x in enumerate(content[8].split(':')) if x == 'AD']
        if len(where_ad) > 0:
            ad_index = where_ad[0]
        else:
            ab_value = 'NA'
        total_dps = [0, 0]
        het_count = 0
        for i, fmt in enumerate(content[9:]):
            if fmt.startswith('0/0'):
                gt = f'{content[3]}/{content[3]}'
            elif fmt.startswith('0/1'):
                het_count += 1
                gt = f'{content[3]}/{content[4]}'
                allele_dps = [int(x) for x in fmt.split(':')[ad_index].split(',')]
                total_dps[0] += allele_dps[0]
                total_dps[1] += allele_dps[1]
            elif fmt.startswith('1/1'):
                 gt = f'{content[4]}/{content[4]}'
            else:
                gt = 'NA'
            data_list[names[i]].append(gt)
        if total_dps[1] <= 0:
            ab_value = 'NA'
        else:
            ab_value = total_dps[1]/sum(total_dps)
        data_list['# HET'].append(het_count)
        data_list['HET alt allele %'].append(ab_value)
        variant_id = f'{content[0]}:{content[1]}'
        for trait_id in gwas_results:
            if variant_id in gwas_results[trait_id]:
                data_list[trait_id].append(gwas_results[trait_id][variant_id])
            else:
                data_list[trait_id].append('NA')

df = pd.DataFrame(data_list)
df.to_excel('genotypes.xlsx')
