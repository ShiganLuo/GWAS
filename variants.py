#!/usr/bin/env python3
from utils import version
from pathlib import Path
import sys
import logging
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.WARNING)
logger.addHandler(console_handler)
# Pipeline Metadata
__version__ = version
__authors__ = 'Shigan Luo'
__email__ = '2530320102@qq.com'
__home__  =  Path(__file__).resolve()
_name = 'variants'
_description = 'call variants from WES/WGS pipeline'


def main(
        log:str
):
    file_handler = logging.FileHandler(log)
    file_handler.setLevel(logging.DEBUG)
    logger.addHandler(file_handler)
    logger.info(__version__)
    pass        

if __name__ == "__main__":
    main("a.log")