#!/usr/bin/env python3
"""Build Nextflow samplesheets from SRA metadata (SraRunTable.csv).

Writes:
  samplesheet_all.csv : one row per run, all runs (for the full analysis)
"""
import pandas as pd
import argparse

CONDITION_MAP = {'10nM E2 treatment for 24h': 'E2',
                'Control ethanol 24h': 'Control'}        
PREFIX = {"Control": "ctrl", "E2": "e2"}

parser = argparse.ArgumentParser()

parser.add_argument("--input", 
                    "-i",
                    type=str,
                    required=True,
                    help="path to sraRunTable")
parser.add_argument("--output",
                    "-o",
                    type=str,
                    default="samplesheet.csv",
                    required=True,
                    help="output file name")
parser.add_argument("--urls",
                    "-u",
                    required=True,
                    help='path to fastq download url')
parser.add_argument("--md5",
                    "-m",
                    required=True,
                    help="path to checksum file")

args = parser.parse_args()

runs = pd.read_csv(args.input)

# fill in treatment and fastq paths columns
sample_sheet = pd.DataFrame()
sample_sheet['fastq'] = runs['Run'].map(lambda x: f"fastq/{x}.fastq.gz")
sample_sheet['condition'] = runs['treatment'].map(CONDITION_MAP)

# create mapping from biosample name 
e2_biosamples = runs.loc[sample_sheet['condition'] == 'E2', 'BioSample'].unique()
ctrl_biosamples = runs.loc[sample_sheet['condition'] == 'Control', 'BioSample'].unique()
e2_names = {bs: f"e2_{i}" for i, bs in enumerate(e2_biosamples)}
ctrl_names = {bs: f"ctrl_{i}" for i, bs in enumerate(ctrl_biosamples)}
combined_names = {**e2_names, **ctrl_names}

# fill in sample column
sample_sheet['run_accession'] = runs['Run']
sample_sheet['sample'] = runs['BioSample'].map(combined_names)

urls = pd.read_csv(args.urls, header=None, names=["url"])
urls["run_accession"] = urls["url"].str.extract(r"(SRR\d+)\.fastq\.gz$", expand=False)
md5 = pd.read_csv(args.md5, sep="\t")  
sample_sheet = pd.merge(pd.merge(sample_sheet, urls, on='run_accession'), md5, on='run_accession')
sample_sheet.to_csv(args.output, index=False)

