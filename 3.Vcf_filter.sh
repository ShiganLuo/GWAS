#!/bin/bash
#PBS -N Splendens_vcf_filter
#PBS -l nodes=10:ppn=12
#PBS -e /home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/err_vcffilter
#PBS -o /home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/out_vcffilter
#PBS -M 421298109@qq.com
#PBS -q cal
#PBS -S /bin/bash
# Kill script if any commands fail
set -e

vcftools=/home/zwc1988/Software/bin/vcftools
plink=/home/zwc1988/Software/bin/plink
raw_vcf=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/Splendens_complex_all.hardfiltered.vcf.gz
outdir=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling
bgzip=/home/zwc1988/Software/bin/bgzip
tabix=/home/zwc1988/Software/bin/tabix

cd ${outdir}

#$bgzip -c $raw_vcf
#$tabix -p vcf 


##########filter the original vcf file and prune the vcf file##############
$vcftools --gzvcf $raw_vcf --maf 0.03 --min-alleles 2 --max-alleles 2 --max-missing 0.5 --minQ 30 --recode --recode-INFO-all --stdout | $bgzip -c > Splendens_410ind_maf0.03miss0.5minQ30allele2.vcf.gz
#maf 0.05;?
#$plink --vcf Splendens_410ind_maf0.03miss0.5minQ30allele2.vcf.gz --recode --allow-extra-chr --out Splendens410_clean
#$plink --file Splendens410_clean --recode --make-bed --allow-extra-chr --out Splendens410_clean

#zcat Splendens_410ind_maf0.03miss0.5minQ30allele2.vcf.gz | sed 's/Chr/chr/g' | gzip > Splendens_410ind_maf0.03miss0.5minQ30allele2_chr.vcf.gz
#$plink --vcf Splendens_410ind_maf0.03miss0.5minQ30allele2_chr.vcf --recode --allow-extra-chr --out Splendens_410ind_maf0.03miss0.5minQ30allele2_chr
#$plink --file Splendens_410ind_maf0.03miss0.5minQ30allele2_chr --make-bed --allow-extra-chr --out Splendens_410ind_maf0.03miss0.5minQ30allele2_chr
#$plink --file Splendens_410ind_maf0.03miss0.5minQ30allele2_chr --indep-pairwise 50 10 0.3 --allow-extra-chr --out Splendens_410ind_pruned
#$plink --file Splendens_410ind_maf0.03miss0.5minQ30allele2_chr --extract Splendens_410ind_pruned.prune.in --allow-extra-chr --recode --out Splendens_410ind_pruned
#$plink --file Splendens_410ind_pruned --make-bed --allow-extra-chr --out Splendens_410ind_pruned
#################################################################################################





