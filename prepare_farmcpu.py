#!/usr/bin/env python3

import sys

vcf, gd_file, gm_file = sys.argv[1:]

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


def num_gt(gt):
    if gt == '0/0' or gt == './.':
        return 0
    elif gt == '0/1':
        return 1
    else:
        return 2

#print('AB_HR\tAB_HET\tAB_HV\tDP_HR\tDP_HET\tDP_HV')
with open(vcf, 'r') as vcf_file:
    names = []
    genotypes = {}
    for line in vcf_file:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                content = line.strip().split('\t')
                for samp in content[9:]:
                    names.append(samp.replace('S', ''))
                    genotypes[samp.replace('S', '')] = {}
            continue
        if '.' in content[4]:
            continue
        content = line.strip().split('\t')
#        ab_nonref = [0, 0]
#        ab_ref = [0, 0]
#        ab_homvar = [0, 0]
#        counts = [0, 0, 0]
#        if content[8].split(':')[1] != 'AD':
#            continue
#        try:
#            for samp in content[9:]:
#                if samp.startswith('0/0'):
#                    ab_ref[0] += [int(x) for x in samp.split(':')[1].split(',')][0]
#                    ab_ref[1] += [int(x) for x in samp.split(':')[1].split(',')][1]
#                    counts[0] += 1
#                elif samp.startswith('0/1'):
#                    ab_nonref[0] += [int(x) for x in samp.split(':')[1].split(',')][0]
#                    ab_nonref[1] += [int(x) for x in samp.split(':')[1].split(',')][1]
#                    counts[1] += 1
#                elif samp.startswith('1/1'):
#                    ab_homvar[0] += [int(x) for x in samp.split(':')[1].split(',')][0]
#                    ab_homvar[1] += [int(x) for x in samp.split(':')[1].split(',')][1]
#                    counts[2] += 1
#        except:
#            continue
#        if avg_ab(ab_nonref) == 'NA' or (float(avg_ab(ab_nonref)) != 0 and float(avg_ab(ab_nonref)) < 0.3):
#            continue
        variant_id = ':'.join(content[:5])
        for i, samp in enumerate(content[9:]):
            genotype = num_gt(samp.split(':')[0])
            this_name = names[i]
            genotypes[this_name][variant_id] = str(genotype)
#        ab_stats = '\t'.join([avg_ab(x) for x in [ab_ref, ab_nonref, ab_homvar]])
#        dp_stats = '\t'.join([sum_dp(x, counts[i]) for i, x in enumerate([ab_ref, ab_nonref, ab_homvar])])
        
#        print(f'{ab_stats}\t{dp_stats}')
#        print(f'At variant {content[2]} ab_ref {ab_ref[1]/sum(ab_ref)} ab_nonref {ab_nonref[1]/sum(ab_nonref)}, ab_nonref {ab_nonref}')
#        print(f'{ab_ref[1]/sum(ab_ref)}\t{ab_nonref[1]/sum(ab_nonref)}')

 
gd_handle = open(gd_file, 'w')
gm_handle = open(gm_file, 'w')
variant_head = ''
vcounter = 1
print('SNP\tChromosome\tPosition', file=gm_handle)
for var_id in list(genotypes[names[0]].keys()):
#    print(var_id)
    chrom, pos = var_id.split(':')[:2]
    varname = f'GUARVAR{vcounter}'
    print(f'{varname}\t{chrom}\t{pos}', file=gm_handle)
    vcounter += 1
    variant_head += f'{varname}\t'
variant_head = variant_head.strip()

print(f'taxa\t{variant_head}', file=gd_handle)
for samp in genotypes:
    matrix_row = '\t'.join(list(genotypes[samp].values()))
    print(f'guar{samp}\t{matrix_row}', file=gd_handle)

gm_handle.close()
gd_handle.close()
