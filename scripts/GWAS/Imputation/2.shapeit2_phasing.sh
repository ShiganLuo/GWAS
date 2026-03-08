#!/bin/bash
#PBS -N Splendens410_Imputation
#PBS -l nodes=6:ppn=12
#PBS -e /home/zwc1988/Data/Betta/3.Splendens_complex/Imputation/error
#PBS -o /home/zwc1988/Data/Betta/3.Splendens_complex/Imputation/out
#PBS -M 421298109@qq.com
#PBS -q cu
#PBS -S /bin/bash
# Kill script if any commands fail
set -e

angsd=/home/zwc1988/Software/bin/angsd
vcftools=/home/zwc1988/Software/bin/vcftools
shapeit2=/home/zwc1988/Software/bin/shapeit
beagle=/home/zwc1988/Software/beagle/beagle.5.1.jar
output=/home/zwc1988/Data/Betta/3.Splendens_complex/Imputation

cd ${output}

#for i in {1..21}
#do
#i=3
#vcf=/home/zwc1988/Data/Betta/3.Splendens_complex/Imputation/Splendens410_clean.vcf

#$vcftools --vcf ${vcf} --chr ${i} --recode --recode-INFO-all --out Splendens410_chr${i}

#$shapeit2 -V Splendens410_chr${i}.recode.vcf -T 24 -O Splendens410_filter_chr${i}_phased

#rm Splendens410_chr${i}.recode.vcf

#done




