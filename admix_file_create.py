#!/usr/bin/env python3

import re
import sys

vcf, genome, out_vcf, out_genome = sys.argv[1:]

contig_map = {}
with open(vcf, 'r') as vcf_file, open(out_vcf, 'w') as out_file:
    previous_contig = ''
    contig_index = 0
    for line in vcf_file:
        if line.startswith('#'):
            if 'contig=' not in line:
                print(line.strip(), file=out_file)
            continue
        content = line.strip().split('\t')
        if len(content[4]) > 1:
            continue
        if content[0] != previous_contig:
            previous_contig = content[0]
            contig_index += 1
            contig_map[content[0]] = contig_index
        print(line.strip().replace(content[0], f'{contig_index}'), file=out_file)

with open(genome, 'r') as genome_file, open(out_genome, 'w') as out_genome_file:
    for line in genome_file:
        if line.startswith('>'):
            chr_name = line.strip().split()[0][1:]
            if chr_name not in contig_map:
                keep_chr = False
            else:
                keep_chr = True
                print(f'>{contig_map[chr_name]}', file=out_genome_file)
        else:
            if keep_chr:
                print(line.strip(), file=out_genome_file)

