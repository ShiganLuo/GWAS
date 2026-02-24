#!/usr/bin/env python3
from utils import version
from pathlib import Path
import argparse
import logging


formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
logger.addHandler(console_handler)


# Pipeline Metadata
__version__ = version
__authors__ = "Shigan Luo"
__email__ = "2530320102@qq.com"
__home__ = Path(__file__).resolve()
_name = "variants"
_description = "call variants from WES/WGS pipeline"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=_description)
    parser.add_argument("--assay", choices=["WES", "WGS"], default="WES", help="Assay type driving reference selection.")
    parser.add_argument("--sample", required=True, help="Sample identifier.")
    parser.add_argument("--germline_vcf", required=True, help="Path to the filtered germline VCF.")
    parser.add_argument("--somatic_vcf", required=True, help="Path to the filtered somatic VCF.")
    parser.add_argument("--log", required=True, help="Log file path.")
    return parser.parse_args()


def main():
    args = parse_args()

    file_handler = logging.FileHandler(args.log)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    if not any(isinstance(h, logging.FileHandler) for h in logger.handlers):
        logger.addHandler(file_handler)
    else:
        logger.addHandler(file_handler)

    logger.info(f"variants pipeline version: {__version__}")
    logger.info(f"assay: {args.assay}; sample: {args.sample}")
    logger.info(f"germline VCF: {args.germline_vcf}")
    logger.info(f"somatic VCF: {args.somatic_vcf}")
    logger.info("placeholder: integrate downstream annotation/aggregation here")


if __name__ == "__main__":
    main()