#!/bin/bash
indat=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/fq.all_ind.list
function sta(){
    dir=/home/zwc1988/Data/Betta
    out=/home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/fq_Allpath.txt
    ID=$1
    # path=($(find $dir -name $ID*fq.gz))
    path=($(find $dir -type f -regex ".*/${ID}_[^/]*fq.gz"))
    n=${#path[@]}
    echo ${path[@]}
    for i in ${path[@]};do
        size+="$(du -sh ${i})\t"
    done
    echo -e "${ID}\t${n}\t${size}" >> $out
}
export -f sta
cat $indat | parallel -j 10 sta {1}
# sta BSU32