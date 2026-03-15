#!/bin/bash
#PBS -N Splendens_complex_GATK_filter
#PBS -l nodes=4:ppn=6
#PBS -e /home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/SNPcalling_filter.err
#PBS -o /home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/SNPcalling_filter.out
#PBS -M 421298109@qq.com
#PBS -q cu
#PBS -S /bin/bash
# Kill script if any commands fail
set -e

platypus=/home/zwc1988/Software/platypus/Platypus_0.8.1/Platypus.py
GATK=/home/zwc1988/Software/gatk/gatk/gatk-package-4.2.0.0-local.jar
RefSeq=/home/zwc1988/Data/Betta/1.Betta_splendens/RefSeq/BS_chrRef.fa
outdir=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling
bgzip=/home/zwc1988/Software/bin/bgzip
tabix=/home/zwc1988/Software/bin/tabix

cd ${outdir}

###过滤结果
#使用SelectVariants，选出SNP
#java -Xmx200G -jar $GATK SelectVariants \
#     -select-type SNP \
#     -V Splendens_complex_raw.vcf.gz \
#     -O Splendens_complex_snp.vcf.gz

#为SNP作硬过滤
#java -Xmx200G -jar $GATK VariantFiltration \
#     -V Splendens_complex_snp.vcf.gz \
#     --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || ReadPosRankSum < -8.0" \
#     --filter-name "SNP_Filter" \
#     -O Splendens_complex_snp.hardfilter.vcf.gz

#使用SelectVariants，选出Indel
#java -Xmx200G -jar $GATK SelectVariants \
#     -select-type INDEL \
#     -V Splendens_complex_raw.vcf.gz \
#     -O Splendens_complex_indel.vcf.gz

#为Indel作过滤
#java -Xmx200G -jar $GATK VariantFiltration \
#     -V Splendens_complex_indel.vcf.gz \
#     --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
#     --filter-name "INDEL_Filter" \
#     -O Splendens_complex_indel.hardfilter.vcf.gz

#重新合并过滤后的SNP和Indel
#java -Xmx200G -jar $GATK MergeVcfs \
#    -I Splendens_complex_snp.hardfilter.vcf.gz \
#    -I Splendens_complex_indel.hardfilter.vcf.gz \
#    -O Splendens_complex_all.hardfilter.vcf.gz

#java -Xmx80G -jar $GATK SelectVariants \
#	--select "vc.isNotFiltered()" \
#	-R $RefSeq \
#	-V Splendens_complex_all.hardfilter.vcf.gz \
#	-O Splendens_complex_all.hardfiltered.vcf.gz




