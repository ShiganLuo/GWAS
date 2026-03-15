#? -O tmb_tsv=~{out_dir}/~{config.info.sample_name}_tmb.txt
#? -O c_tmb_tsv=~{out_dir}/~{config.info.sample_name}_cTMB.stat.tsv
#? -O snv_indel_tsv=~{out_dir}/~{config.info.sample_name}.anno.filter.add.tsv

oud=~{out_dir}
sample=~{config.info.sample_name}
somatic_snv_indel_anno_filter_tsv=~{inputs.snv_indel_tsv}
module_bin=~{config.dirs.module_bin}

python ${module_bin}/calculate_tmb_ctmb.py \
    --vep ~{inputs.vep_vcf} \
    --filter "${somatic_snv_indel_anno_filter_tsv}" \
    --absmaf "~{inputs.abs_maf}" \
    --sample ${sample} \
    --size 38 \
    --vaf 0.05 \
    --oud ${oud} \
    --thr 0.15
