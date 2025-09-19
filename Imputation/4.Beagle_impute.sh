#!/bin/bash
#PBS -N Splendens615_Imputation
#PBS -l nodes=10:ppn=12
#PBS -e /home/zwc1988/Data/Betta/3.Splendens_complex/Imputation/log/Imputation615.err
#PBS -o /home/zwc1988/Data/Betta/3.Splendens_complex/Imputation/log/Imputation615.out
#PBS -M 2530320102@qq.com
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
RefSeq=/home/zwc1988/Data/Betta/3.Splendens_complex/genome/ncbi_dataset/data/GCF_900634795.4/GCF_900634795.4_fBetSpl5.4_genomic.fna


cd ${output}



output=/home/zwc1988/Data/Betta/3.Splendens_complex/Imputation/chr
chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 19 21 22 24 
"NW_026577841.1" "NW_026577842.1" "NW_026577843.1" "NW_026577844.1"
"NW_026577845.1" "NW_026577846.1" "NW_026577847.1" "NW_026577848.1"
"NW_026577849.1" "NW_026577850.1" "NW_026577851.1" "NW_026577852.1"
"NW_026577853.1" "NW_026577854.1" "NC_026581.1")
# a=$1
# i=${chromosomes[${a}]}
# echo $i
# vcf=${output}/Splendens615_filter_chr${i}_phased.vcf
# echo ${vcf}
# echo $i
###########Imputation###########
# java -Xmx200G -jar $beagle \
#     gt=${vcf} \
#     nthreads=48 \
#    chrom=${i} \
#     window=0.4 \
#     overlap=0.2 \
#     out=${output}/Splendens615_phased_chr${i}
# echo "#### ${i} done ####"

###################################################################
###########fix the header of the vcf file then gather###########
# for i in ${chromosomes[@]}
# do
# java -Xmx200G -jar ${picard} FixVcfHeader \
#       I=./chr/Splendens615_phased_chr${i}.vcf.gz \
#       O=./chr/Splendens615_phased_chr${i}.fixed.vcf.gz \
#       HEADER=Splendens615.header.vcf
# done
########################################################################
chr_list=/home/zwc1988/Data/Betta/3.Splendens_complex/Imputation/Splendens615.chrall.vcf.list

java -Xmx100G -jar $GATK GatherVcfs -I ${chr_list} -O Splendens615_phased.chr.vcf.gz

# tabix -p vcf Splendens615_phased.vcf.gz





