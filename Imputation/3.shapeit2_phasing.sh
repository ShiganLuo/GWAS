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
output=/home/zwc1988/Data/Betta/3.Splendens_complex/Imputation/chr

cd ${output}
chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 19 21 22 24 
"NW_026577841.1" "NW_026577842.1" "NW_026577843.1" "NW_026577844.1"
"NW_026577845.1" "NW_026577846.1" "NW_026577847.1" "NW_026577848.1"
"NW_026577849.1" "NW_026577850.1" "NW_026577851.1" "NW_026577852.1"
"NW_026577853.1" "NW_026577854.1" "NC_026581.1")
vcf=/home/zwc1988/Data/Betta/3.Splendens_complex/Imputation/Splendens_615ind_maf0.02miss0.5minQ30allele2.numberChr.pass.vcf
a=$1
i=${chromosomes[${a}]}
# for i in ${chromosomes[@]}
# do
echo $i

# $vcftools --vcf ${vcf} --chr ${i} --recode --recode-INFO-all --out ${output}/Splendens615_chr${i}

# $shapeit2 -V ${output}/Splendens615_chr${i}.recode.vcf -T 24 --force -O ${output}/Splendens615_filter_chr${i}_phased
$shapeit2 -convert --input-haps ${output}/Splendens615_filter_chr${i}_phased --output-vcf ${output}/Splendens615_filter_chr${i}_phased.vcf

# rm ${output}/Splendens615_chr${i}.recode.vcf

# done




