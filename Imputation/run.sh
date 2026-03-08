#!/bin/bash
# 定义节点列表
nodes=("cu01" "cu02" "cu03" "cu04" "cu05" "cu06" "cu07" "cu08" 
"cu09" "cu10" "cu11" "cu12" "cu13" "cu14" "cu15" "cu16" "cu17" "cu18" 
"cu19" "cu20" "cu21" "cu22" "cu23" "cu24" "cu25" "cu26" "cu27" "cu29" 
"cu30" "cu31" "cu32" "cu33" "cu34" "cu35" "cu36" "cu37" "cu38" "cu39" 
"cu40" "cu41" "cu42" "cu43" "cu44" "cu45" "cu46" "cu47" "cu48" "cu49" "cu50" 
"cu51" "cu52" "cu53" "cu54" "cu55" "cu56" "cu57" "cu58" "cu59" "cu60" "cu61" 
"cu62" "cu63" "cu64" "cu65" "cu66" "cu67" "cu68" "cu69" "cu70" "cu71" "cu72" 
"cu73" "cu74" "cu75" "cu76" "cu77" "cu78" "cu79" "cu80") 
for i in {0..35}
do
echo $i
node=${nodes[$i]}
echo "Logging into $node and executing command..."
#########################################3.shapeit2_phasing.sh #########################################
# ssh $node "
#         cd /home/zwc1988/Data/Betta/3.Splendens_complex/Imputation;
#         nohup bash 3.shapeit2_phasing.sh ${i} > log/convert_${i}.log 2>&1 &
#     "
# echo "cd /home/zwc1988/Data/Betta/3.Splendens_complex/Imputation"
# echo "nohup bash 3.shapeit2_phasing.sh ${i} > log/convert_${i}.log 2>&1 &"
#########################################4.Beagle_impute.sh#########################################
ssh $node "
        cd /home/zwc1988/Data/Betta/3.Splendens_complex/Imputation;
        nohup bash 4.Beagle_impute.sh ${i} > log/Imputation_${i}.log 2>&1 &"
echo "cd /home/zwc1988/Data/Betta/3.Splendens_complex/Imputation"
echo "nohup bash 4.Beagle_impute.sh ${i} > log/Imputation_${i}.log 2>&1 &"        
done