import os
import re
import yaml
from typing import List
from snakemake.utils import min_version


# Snakemake version when Singularity support was added
min_version("4.2.0")


#======================================================
# Config files
#======================================================
configfile: "config.yaml"

with open('cluster.yaml', 'r') as fh:
    cluster_config = yaml.load(fh)


#======================================================
# Functions and Classes
#======================================================
class InvalidBarcode(Exception):
    __module__ = Exception.__module__


def barcode_parser(barcodes_string: str) -> List[str]:
    """Parses the barcodes string and ensures they follow correct format"""
    msg = "Barcode must be of the form BC01. That is, BC followed by 2 digits."
    regex = r'\bBC\d{2}\b'
    barcodes = barcodes_string.split()
    for barcode in barcodes:
        if not (len(barcode) == 4 and re.match(regex, barcode)):
            raise InvalidBarcode(barcode + '\n' + msg)
    return barcodes


#======================================================
# Global variables
#======================================================
RULES_DIR = 'rules'
MULTIPLEXED = config["multiplexed"]
if MULTIPLEXED:
    SAMPLES = barcode_parser(config["barcodes"])
else:
    SAMPLES = [config["sample_name"]]


#======================================================
# Rules
#======================================================
rule all:
    input:
        expand("docs/report_{sample}.html", sample=SAMPLES)


# the snakemake files that run the different parts of the pipeline
include: os.path.join(RULES_DIR, 'porechop.smk')
include: os.path.join(RULES_DIR, 'align.smk')
include: os.path.join(RULES_DIR, 'mykrobe.smk')
include: os.path.join(RULES_DIR, 'reports.smk')
