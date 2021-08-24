#!/bin/bash

VCF=$1
PVAL_FOLDER=$2
GMFILE=$3

grep -vP '^#' $VCF | cut -f 1-2 | awk 'OFS="\t" { $2=$2-151 ; $3=$2+301 ; if ($2 < 0) $2=0; print }' > ${VCF%%.vcf}.300bp.bed
bedtools getfasta -fi reference.fa -bed ${VCF%%.vcf}.300bp.bed > ${VCF%%.vcf}.300bp.fasta
bedtools intersect -wo -a $VCF -b transcripts_alignments.bed | cut -f 1-5,176-182 > ${VCF%%.vcf}_matched_to_transcripts.tsv
./parse_vcf_fcpu.py ${VCF%%.vcf}.300bp.fasta $VCF ${VCF%%.vcf}_matched_to_transcripts.tsv $PVAL_FOLDER $GMFILE
