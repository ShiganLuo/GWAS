#!/bin/bash
#PBS -N Splendens410_Imputation
#PBS -l nodes=10:ppn=12
#PBS -e /home/zwc1988/Data/Betta/3.Splendens_complex/Imputation/Imputation410.err
#PBS -o /home/zwc1988/Data/Betta/3.Splendens_complex/Imputation/Imputation410.out
#PBS -M 421298109@qq.com
#PBS -q cu
#PBS -S /bin/bash
# Kill script if any commands fail
set -e

angsd=/home/zwc1988/Software/bin/angsd
vcftools=/home/zwc1988/Software/bin/vcftools
beagle=/home/zwc1988/Software/beagle/beagle.5.1.jar
GATK=/home/zwc1988/Software/gatk/gatk/gatk-package-4.2.0.0-local.jar
picard=/home/zwc1988/Software/picard/picard.jar
plink=/home/zwc1988/Software/bin/plink
output=/home/zwc1988/Data/Betta/3.Splendens_complex/Imputation
Refseq=/home/zwc1988/Data/Betta/1.Betta_splendens/RefSeq/BS_chrRef.fa

cd ${output}

#file=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/Splendens410_clean
#$plink --file $file --allow-extra-chr --recode vcf --out ${output}/Splendens410_clean

vcf=${output}/Splendens410_clean.vcf

#for i in {3..21}
#do
#i=2
#java -Xmx200G -jar $beagle \
#     gt=${vcf} \
#     nthreads=48 \
#    chrom=${i} \
#     window=0.4 \
#     overlap=0.2 \
#     out=${output}/Splendens410_phased_chr${i}
#done

###################################################################
###########fix the header of the vcf file then gather###########
#for i in {1,2,3,4,5,6,7,8,9,10,.21}

#do
#i=2
#java -Xmx200G -jar ${picard} FixVcfHeader \
#       I=Splendens410_phased_chr${i}.vcf.gz \
#       O=Splendens410_phased_chr${i}.fixed.vcf.gz \
#       HEADER=Header.vcf
#done
########################################################################
#chr_list=/home/zwc1988/Data/Betta/3.Splendens_complex/Imputation/chrall.vcf.list

#java -Xmx100G -jar $GATK GatherVcfs -I ${chr_list} -O ${output}/Splendens410_phased.vcf.gz

tabix -p vcf Splendens410_phased.vcf.gz





