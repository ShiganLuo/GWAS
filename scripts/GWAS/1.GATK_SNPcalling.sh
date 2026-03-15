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
GATK=/home/zwc1988/Software/GATK4/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar
samtools=/home/zwc1988/Software/bin/samtools
bgzip=/home/zwc1988/Software/bin/bgzip
tabix=/home/zwc1988/Software/bin/tabix
RefSeq=/home/zwc1988/Data/Betta/1.Betta_splendens/RefSeq/BS_chrRef.fa
bamfiledir=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/Bam_file
#bamfilelist=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/Bamfile.list
GVCFdir=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/GVCF
outdir=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling

cd ${outdir}

#################build reference dict##########################################
#$gatk CreateSequenceDictionary \
#      -R $RefSeq \
#      -O /home/zwc1988/Data/Betta/1.Betta_splendens/RefSeq/BS_chrRef.dict
#################################################################################

#########mark PCR duplicated reads###############################################
# INDS=($(for i in $bamfiledir/*_sort.bam; do echo $(basename ${i%_sort*}); done))
# for IND in ${INDS[@]}
# do
# BAM1=${bamfiledir}/${IND}_sort.bam
# java -jar $GATK MarkDuplicates \
#       -I $BAM1 \
#       -O ${IND}_MarkDup.bam \
#       -M ${IND}.metrics
# $samtools index ${IND}_MarkDup.bam
# #rm -f $BAM1
# done
#########################################################################################

###########################produce gVCF file#############################################
#INDS=($(for i in $bamfiledir/*_MarkDup.bam; do echo $(basename ${i%_MarkDup*}); done))
#SplendensIND=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/Splendens410Ind.list
#INDS=`cat $SplendensIND | sed -n '401,410p' | cut -d "/" -f 9 | cut -d "_" -f 1`
#for IND in ${INDS[@]}
#do
#IND=BSG2
#java -jar $GATK HaplotypeCaller \
#     -R $RefSeq \
#     -I ${bamfiledir}/${IND}_MarkDup.bam -O ${GVCFdir}/${IND}_gatk.g.vcf \
#     --emit-ref-confidence GVCF \
#     -stand-call-conf 30 && echo "**$IND gvcf done **"
#done
################################################################################

#################merge multiple gVCF file########################################
#find ${GVCF}/ -name "*.g.vcf" > input.list
#input_list=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/gvcf.input.list
#java -Xmx240G -jar $GATK CombineGVCFs \
#     --reference $RefSeq \
#     --variant ${input_list} \
#     --output ${outdir}/Combined/Splendens_complex_combined.g.vcf.gz 
#$tabix -p vcf Splendens_complex_combined.g.vcf.gz
#rm ${GVCF}/*.g.vcf
###########################################################################

####################################################################################
#input_list=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/gvcf.input.list
#java -Xmx220g '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' -jar $GATK GenomicsDBImport \
#     --variant ${input_list} \
#     --genomicsdb-workspace-path ${outdir}/New_combined_Chr1 \
#     --intervals Chr1
###################################################################################

###################joint genotyping##############################################
#java -Xmx240g -jar $GATK GenotypeGVCFs \
#     -R $RefSeq \
#     -V ${outdir}/Combined/Splendens_complex_combined.g.vcf.gz \
#     -G StandardAnnotation \
#     -O ${outdir}/Splendens_complex_raw.vcf.gz

#rm combined.g.vcf

#$bgzip -f Splendens_complex_raw.vcf
#$tabix -p vcf Splendens_complex_raw.vcf.gz
#################################################################################




