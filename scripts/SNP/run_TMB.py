from burden.mutation_burden import compute_exonic_mutation_burden_from_vcf_dir
from burden.plot import plot_mutation_burden
from vep.vep_ann import VEP_Annotation
from burden.calculate_tmb_ctmb import calculate_vep
import pandas as pd
import os
import glob

def run_mutation_burden(
    vcf_dir:str,
    gtf_file:str,
    outprefix:str,
    sample_rename_func = None
):
    df = compute_exonic_mutation_burden_from_vcf_dir(
        vcf_dir=vcf_dir,
        gtf_file=gtf_file,
        qual_filter=30,
        chunksize=100000
    )
    def get_sample_id(sample):
        return sample.split("_")[1]   # GSMxxxxxxx
    def get_condition(sample):
        return sample.split("_")[1].split("-")[0]   # GSExxxxxx
    df = df.sort_values("sample").reset_index(drop=True)
    df.to_csv(outprefix + ".csv", index=False)
    plot_mutation_burden(df, outprefix + ".png", get_sample_id_func=get_sample_id, get_condition_func=get_condition,figsize=(10, 6))
    pass

def run_TMB_calc(
    vcf_dir:str,
    outdir:str,
    cds_size_mb:float = 34.0,
    **vep_kwargs
):
    os.makedirs(outdir, exist_ok=True)
    res_dict = {}
    analysis = VEP_Annotation(**vep_kwargs)
    vcf_files = glob.glob(os.path.join(vcf_dir, "*.vcf.gz"))
    for vcf_file in vcf_files:
        sample = os.path.basename(vcf_file).split(".vcf.gz")[0].split("_")[1]
        sample_outdir = os.path.join(outdir, sample)
        os.makedirs(sample_outdir, exist_ok=True)
        outfile = os.path.join(sample_outdir, f"{sample}.vep.vcf")
        # analysis.annotate_sv_vep(vcf_file, outfile)
        tmbcount, indelcount = calculate_vep(outfile,vaf=0.05, tumor_sample=sample)
        res_dict[sample] = {"TMBcount": tmbcount,"TMB" : tmbcount/cds_size_mb, "Indelcount": indelcount, "Indel(Mutation/Mb)": indelcount/cds_size_mb}
    res_df = pd.DataFrame.from_dict(res_dict, orient='index').reset_index().rename(columns={"index": "sample"})
    res_df.to_csv(os.path.join(outdir, "tmb_results.csv"), index=False)
    
def main():
    config1 = {
        "vcf_dir": "/data/pub/zhousha/20260207_Exome/output/Exome/gatk/mutect2-vcf/GRCm39",
        "gtf_file": "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf",
        "outprefix": "/data/pub/zhousha/20260207_Exome/output/Exome/results/TMB/exonic_mutation_burden"
    }
    run_mutation_burden(**config1)

if __name__ == "__main__":
    # main()
    vcf_dir = "/data/pub/zhousha/20260207_Exome/output/Exome/gatk/mutect2-vcf/GRCm39"
    outdir = "/data/pub/zhousha/20260207_Exome/output/Exome/results/TMB"
    vcf_files = run_TMB_calc(vcf_dir=vcf_dir, outdir=outdir)
