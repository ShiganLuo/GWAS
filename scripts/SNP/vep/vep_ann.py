import subprocess
import os
import shutil
import tempfile
import gzip
import logging
from pathlib import Path
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] [%(name)s] %(message)s')
logger = logging.getLogger(__name__)


class VEP_Annotation:
    def __init__(self, vep_cache_dir:str = "~/.vep", species:str = "mus_musculus", assembly:str = "GRCm39"):
        """
        初始化分析类
        :param vep_cache_dir: VEP 缓存的根目录
        :param species: 物种名称 (如 mus_musculus)
        :param assembly: 基因组版本 (如 GRCm39)
        """
        self.vep_cache_dir = str(Path(vep_cache_dir).expanduser().resolve()) # 支持 ~ 和绝对路径
        self.species = species
        self.assembly = assembly
        
        # 自动创建缓存根目录
        if not os.path.exists(self.vep_cache_dir):
            os.makedirs(self.vep_cache_dir, exist_ok=True)
            logger.info(f"Created VEP cache directory: {self.vep_cache_dir}")


    def _run_cmd(self, cmd:list):
        """
        执行外部命令，返回 stdout
        - 命令不存在：给出清晰提示
        - 命令执行失败：打印 stdout / stderr
        """
        cmd_str = " ".join(cmd)
        cmd_bin = cmd[0]

        logger.info(f"Running: {cmd_str}")

        # 1️⃣ 预检查：命令是否存在（比 FileNotFoundError 更友好）
        if shutil.which(cmd_bin) is None:
            logger.error(f"Command not found: '{cmd_bin}'")
            logger.error("Please make sure it is installed and in $PATH")
            raise RuntimeError(f"Command not found: {cmd_bin}")

        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )

            if result.stdout:
                logger.info(f"Command Output:\n{result.stdout}")

            return result.stdout

        except subprocess.CalledProcessError as e:
            logger.error(f"Command failed with return code {e.returncode}")
            logger.error(f"STDOUT:\n{e.stdout or '[empty]'}")
            logger.error(f"STDERR:\n{e.stderr or '[empty]'}")
            raise RuntimeError(
                f"Command execution failed: {cmd_str}"
            ) from e


    def vep_annotation_install(self):
        """安装类指定的 VEP 缓存"""
        cmd = [
            "vep_install", "-a", "cf", 
            "-s", self.species, 
            "-y", self.assembly, 
            "-c", self.vep_cache_dir
        ]
        self._run_cmd(cmd)


    def annotate_sv_vep(
            self, 
            in_vcf:str, 
            outfile:str,
            result_format:str = "vcf"
        ):
        """
        Function: Annotate SVs using VEP with comprehensive options for structural variants.
        Parameters:
        - in_vcf: input VCF file containing SVs
        - outfile: output file path for annotated results
        - result_format: output format (vcf or tab) --{result_format} will be passed to VEP
        """
        # 检查特定物种的缓存子目录是否存在
        species_cache = os.path.join(self.vep_cache_dir, self.species)
        logger.info(species_cache)
        if not os.path.exists(species_cache):
            logger.warning(f"Cache for {self.species} not found. Attempting install...")
            self.vep_annotation_install()

        cmd = [
            "vep", "-i", in_vcf, "-o", outfile,
            "--cache", "--dir_cache", self.vep_cache_dir,
            "--species", self.species,
            "--assembly", self.assembly,
            "--format", "vcf", f"--{result_format}", "--force_overwrite",
            "--everything", "--pick", "--per_gene", "--offline"
        ]
        self._run_cmd(cmd)
        logger.info(f"VEP annotation finished: {outfile}")


if __name__ == "__main__":
    analysis = VEP_Annotation(
        vep_cache_dir="~/.vep",
        species="mus_musculus",
        assembly="GRCm39"
    )
    vcf_file = "/data/pub/zhousha/20260207_Exome/output/Exome/gatk/mutect2-vcf/GRCm39/E14-1_TLSC-1.vcf.gz"
    out_vcf = "/data/pub/zhousha/20260207_Exome/output/Exome/results/vep/E14-1_TLSC-1.vep.vcf"
    analysis.annotate_sv_vep(vcf_file, out_vcf)