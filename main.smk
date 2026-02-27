shell.prefix("set -x; set -e;")
import os
import logging
import sys
EXECUTION_DIR = os.getcwd()
SNAKEFILE_FULL_PATH = workflow.snakefile
ROOT_DIR = os.path.dirname(SNAKEFILE_FULL_PATH)
logging.info(f"Execution directory: {EXECUTION_DIR}")
logging.info(f"Workflow Root directory: {ROOT_DIR}")

fqdir = config.get('indir', 'data')
outdir = config.get('outdir', 'output')
metadata = config.get('metadata')
log = config.get('log', 'workflow.log')
# paraser metadata and generate sample info for WES and WGS analysis
sys.path.append(f"{ROOT_DIR}/utils")
from meta_utils import MetadataVariantsUtils
metadataVariantsUtils = MetadataVariantsUtils(
    meta=metadata,
    outdir=outdir,
    fastq_dir=fqdir,
    log=log
)
samples_info_dict, pairs, workdir = metadataVariantsUtils.run()


# global variable to store the output files for WES analysis
outfiles = []

# function to get the output files for WES analysis
def get_WES_outfiles(samples_info_dict:dict):
    include: "subworkflow/WES.smk"
    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.layout == 'SE':
            outfiles.append(f"{outdir}/cutadapt/{sample_id}_1.fq.gz")
        elif sample_info.layout == 'PE':
            outfiles.append(f"{outdir}/cutadapt/{sample_id}_1.fq.gz")
            outfiles.append(f"{outdir}/cutadapt/{sample_id}_2.fq.gz")
        else:
            logging.warning(f"Unknown layout for sample {sample_id}: {sample_info['layout']}")
            continue
    return outfiles
get_WES_outfiles(samples_info_dict)
# function to get the output files for WGS analysis
def get_WGS_outfiles():

    return outfiles

rule all:
    input:
        outfiles