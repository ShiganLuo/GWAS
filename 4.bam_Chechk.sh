#!/bin/bash
samtools=/home/zwc1988/Software/bin/samtools
for file in /home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/Bam_file/*
do
echo $file
    $samtools quickcheck $file
done 