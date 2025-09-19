#!/bin/bash

# 定义节点列表
nodes=("fat01" "cu01" "cu02" "cu03" "cu04" "cu05" "cu06" "cu07" "cu08" 
"cu09" "cu10" "cu11" "cu12" "cu13" "cu14" "cu15" "cu16" "cu17" "cu18" 
"cu19" "cu20" "cu21" "cu22" "cu23" "cu24" "cu25" "cu26" "cu27" "cu29" 
"cu30" "cu31" "cu32" "cu33" "cu34" "cu35" "cu36" "cu37" "cu38" "cu39" 
"cu40" "cu41" "cu42" "cu43" "cu44" "cu45" "cu46" "cu47" "cu48" "cu49" "cu50" 
"cu51" "cu52" "cu53" "cu54" "cu55" "cu56" "cu57" "cu58" "cu59" "cu60" "cu61" 
"cu62" "cu63" "cu64" "cu65" "cu66" "cu67" "cu68" "cu69" "cu70" "cu71" "cu72" 
"cu73" "cu74" "cu75" "cu76" "cu77" "cu78" "cu79" "cu80") 
# command="cd /home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling"  # 替换为你要执行的命令

#########mark PCR duplicated reads###############################################
# k=30
# for i in $(seq 575 2 615)
# do
#     j=$((i+1))
#     node=${nodes[k]}
#     echo "$i $j $k $node"
#     echo "Logging into $node and executing command..."
#     ssh $node "
#         cd /home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling;
#         nohup bash 5.GATK_SNPcalling.sh ${i} 2 > log/GATK_SNPcalling/markup/index_${i}_${j}.log 2>&1 &
#     "
#     #注意${i}后面的数字一定要和seq中间的数字对应上;j和i的关系也要和其对上
#     echo "nohup bash 5.GATK_SNPcalling.sh ${i} 2 > log/GATK_SNPcalling/markup/index_${i}_${j}.log 2>&1 &"
#     # ssh $node "pkill -f 5.GATK_SNPcalling.sh"
#     # ssh $node "pkill -f /home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/"
#     k=$((k+1))
# done

###########################produce gVCF file#############################################
# k=0
# for i in $(seq 556 1 614)
# do
#     j=$((i+0))
#     node=${nodes[k]}
#     echo "$i $j $k $node"
#     echo "Logging into $node and executing command..."
#     ssh $node "
#         cd /home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling;
#         nohup bash 5.GATK_SNPcalling.sh ${i} 1 > log/GVCF/index_${i}_${j}.log 2>&1 &
#     "
#     #注意${i}后面的数字一定要和seq中间的数字对应上;j和i的关系也要和其对上；另外要考虑最后一个节点是否够得上(k),应该比预期序列小,以及节点从哪里开始
#     echo "nohup bash 5.GATK_SNPcalling.sh ${i} 1 > log/GVCF/index_${i}_${j}.log 2>&1 &"
#     # ssh $node "pkill -f 5.GATK_SNPcalling.sh"
#     # ssh $node "pkill -f /home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling/"
#     k=$((k+1))
# done
# # ssh cu68 "
# #         cd /home/zwc1988/Data/Betta/3.Splendens_complex/SNPcalling;
# #         nohup bash 5.GATK_SNPcalling.sh 614 1 > log/GVCF/index_614_614.log 2>&1 &
# #     "
# # echo "nohup bash 5.GATK_SNPcalling.sh 614 1 > log/GVCF/index_614_614.log 2>&1 &"

##################4)split into chromosome gvcfs,can't cover exit data###################################
#36chr(0~35),less number ,switch node and i manually;nedd 200G at least,only can run at fat01
# for i in {30..35}
# do
#     nohup bash 5.GATK_SNPcalling.sh $i > log/GVCF/chr_gz/chr_${i}.log 2>&1 &
#     echo "nohup bash 5.GATK_SNPcalling.sh $i > log/GVCF/chr_gz/chr_${i}.log 2>&1 &"
# done
