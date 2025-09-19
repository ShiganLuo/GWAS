#!/bin/bash
#PBS -N Splendens_align
#PBS -l nodes=80:ppn=24
#PBS -j oe
#PBS -o /home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/log/wild21_36.log
#PBS -M 2530320102@qq.com
#PBS -q cu
#PBS -S /bin/bash
##Kill script if any commands fail
set -e

# a simple script to map lots of individuals
bwa=/home/zwc1988/Software/bin/bwa
samtools=/home/zwc1988/Software/bin/samtools

REF=/home/zwc1988/Data/Betta/1.Betta_splendens/RefSeq/updated_assembly/Siamese_fighting_fish_newcoord.fa
output=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/Bam_file

dir1=/home/zwc1988/Data/Betta/3.Splendens_complex/fq_file
file1=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/fq.Betta59ind.list


INDS=`cat $file1 | sed -n '1,59p'`

for IND in ${INDS[@]}
do
	# declare variables
	FORWARD=${dir1}/${IND}_1_clean.fq.gz
	REVERSE=${dir1}/${IND}_2_clean.fq.gz
	#such as:/home/zwc1988/Data/Betta/3.Splendens_complex/fq_file/BSS2_1_clean.fq.gz
	OUTPUT=${output}/${IND}_sort.bam
	RG="@RG\\tID:${IND}\\tLB:${IND}\\tPL:ILLUMINA\\tSM:${IND}"
	# then align and sort
	echo "Aligning $IND with bwa"
time	$bwa mem -M -t 24 -R ${RG} \
	$REF \
	$FORWARD \
	$REVERSE | \
time	$samtools view -bSu | \
time	$samtools sort -@ 24 -m 2000000000 -O bam -o $OUTPUT

done



