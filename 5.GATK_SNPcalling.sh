#!/bin/bash
#PBS -N Splendens_complex_GATK
#PBS -l nodes=8:ppn=12
#PBS -e /home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/SNPcalling_GATK.err
#PBS -o /home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/SNPcalling_GATK.out
#PBS -M 421298109@qq.com
#PBS -q cal
#PBS -S /bin/bash
# Kill script if any commands fail
set -e

platypus=/home/zwc1988/Software/platypus/Platypus_0.8.1/Platypus.py
GATK=/home/zwc1988/Software/gatk/gatk/gatk-package-4.2.0.0-local.jar
gatk=/home/zwc1988/Software/gatk/gatk/gatk
samtools=/home/zwc1988/Software/bin/samtools
bgzip=/home/zwc1988/Software/bin/bgzip
tabix=/home/zwc1988/Software/bin/tabix

RefSeq=/home/zwc1988/Data/Betta/3.Splendens_complex/genome/ncbi_dataset/data/GCF_900634795.4/GCF_900634795.4_fBetSpl5.4_genomic.fna
bamfiledir=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/Bam_file
#bamfilelist=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/Bamfile.list
Markupdir=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/MarkDup
GVCFdir=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/GVCF
outdir=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling

cd ${outdir}

#################1)build reference dict,can't cover exit data##########################################
# $gatk CreateSequenceDictionary \
#      -R $RefSeq \
#      -O /home/zwc1988/Data/Betta/3.Splendens_complex/genome/ncbi_dataset/data/GCF_900634795.4/GCF_900634795.4_fBetSpl5.4_genomic.dict
# samtools faidx $RefSeq
#################################################################################

# #########2)mark PCR duplicated reads###############################################
# INDS=($(for i in $bamfiledir/*_sort.bam; do echo $(basename ${i%_sort*}); done))
# INDregion=${INDS[@]:${markup_index}:${markup_n}}
# echo $INDregion
# for IND in ${INDregion[@]}
# do
# echo $IND
# BAM1=${bamfiledir}/${IND}_sort.bam
# java -jar $GATK MarkDuplicates \
#       -I $BAM1 \
#       -O ${Markupdir}/${IND}_MarkDup.bam \
#       -M ${Markupdir}/${IND}.metrics
# $samtools index ${Markupdir}/${IND}_MarkDup.bam
# # rm -f $BAM1
# done
#########################################################################################

###########################3)produce gVCF file#############################################
# markup_index=$1
# markup_n=$2
# INDS=($(for i in ${Markupdir}/*_MarkDup.bam; do echo $(basename ${i%_MarkDup*}); done))
# INDregion=${INDS[@]:${markup_index}:${markup_n}}
# echo $INDregion
# for IND in ${INDregion[@]}
# do
# java -jar $GATK HaplotypeCaller \
#     -R $RefSeq \
#     -I ${Markupdir}/${IND}_MarkDup.bam -O ${GVCFdir}/${IND}_gatk.g.vcf \
#     --emit-ref-confidence GVCF \
#     -stand-call-conf 30 && echo "**$IND gvcf done **"
# done
################################################################################

#################4)merge multiple gVCF file########################################
#find ${GVCF}/ -name "*.g.vcf" > input.list
#input_list=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/gvcf.input.list
#java -Xmx240G -jar $GATK CombineGVCFs \
#     --reference $RefSeq \
#     --variant ${input_list} \
#     --output ${outdir}/Combined/Splendens_complex_combined.g.vcf.gz 
#$tabix -p vcf Splendens_complex_combined.g.vcf.gz
#rm ${GVCF}/*.g.vcf
###########################################################################
##################4)split into chromosome gvcfs,can't cover exit data###################################
# chr_n=$1
input_file=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/gvcf_615ind.input.list
chr=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/chr.txt
mapfile -t chrId < ${chr}
mapfile -t input_list < ${input_file}
# # echo ${input_list[0]} ${input_list[1]}
for i in ${input_list[@]}
do
    V+="-V ${i} "
done
echo ${V}
# i=${chrId[${chr_n}]}
# echo $i
# echo ${outdir}/GVCF/chr_gz/${i}
# java -Xmx200g -jar $GATK GenomicsDBImport \
#     ${V} \
#     -genomicsdb-workspace-path ${outdir}/GVCF/chr_gz/${i} \
#     -intervals ${i} \

# export TILEDB_DISABLE_FILE_LOCKING=1
# java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx500g -jar $GATK GenotypeGVCFs \
#     -R ${RefSeq} \
#     -V gendb://${outdir}/GVCF/chr_gz/${i} \
#     -O ${outdir}/GVCF/chr_gz/Betta.${i}.vcf.gz
# echo "### ${i} done ###"
### all sites,can be ajusted in pixy
m=$1
n=$2
for((a=${m};a<=${n};a++))
do
i=${chrId[${a}]}
echo $i
java -Xmx200g -jar $GATK GenomicsDBImport \
    ${V} \
    -genomicsdb-workspace-path ${outdir}/GVCF/allSites/${i} \
    -L ${i} 
java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx500g -jar $GATK GenotypeGVCFs \
    -R ${RefSeq} \
    -V gendb://${outdir}/GVCF/allSites/${i} \
    -all-sites \
    -L ${i} \
    -O ${outdir}/GVCF/allSites/Betta.all.${i}.vcf.gz
echo "### ${i} done ###"
done
###########################################################################

#######################4.1)merge chr.gz#############################################################
# chr=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/chr.txt
# mapfile -t chrId < ${chr}
# for i in ${chrId[@]}
# do
#     I+="-I ${outdir}/GVCF/chr_gz/Betta.${i}.vcf.gz " 
# done
# echo "${I}"
# java -Xmx200g -jar $GATK MergeVcfs \
#     ${I} \
#     -O ${outdir}/GVCF/chr_gz/Betta.merge.vcf.gz 

  
###################################################################################
####################################################################################
# input_list=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/gvcf.input.list
# java -Xmx220g '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' -jar $GATK GenomicsDBImport \
#     --variant ${input_list} \
#     --genomicsdb-workspace-path ${outdir}/New_combined_Chr1 \
#     --intervals Chr1
###################################################################################
###################5)joint genotyping##############################################
# java -Xmx240g -jar $GATK GenotypeGVCFs \
#     -R $RefSeq \
#     -V ${outdir}/GVCF/chr_gz/Betta.merge.vcf.gz  \
#     -G StandardAnnotation \
#     -O ${outdir}/GVCF/chr_gz/Splendens_complex_raw.vcf.gz

#rm combined.g.vcf

#$bgzip -f Splendens_complex_raw.vcf
#$tabix -p vcf Splendens_complex_raw.vcf.gz
#################################################################################




