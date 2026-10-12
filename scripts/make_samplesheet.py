#!/usr/bin/env python3
"""Build Nextflow samplesheets from SRA metadata (SraRunTable.csv).

Run it once to see the columns and candidate condition values, set CONDITION_MAP below, then run it again.

Writes:
  samplesheet_all.csv : one row per run, all runs (for the full analysis)
"""
import pandas as pd

CONDITION_MAP = {'10nM E2 treatment for 24h': 'E2',
                'Control ethanol 24h': 'Control'}        
PREFIX = {"Control": "ctrl", "E2": "e2"}

runs = pd.read_csv("SraRunTable.csv")

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
sample_sheet['sample'] = runs['BioSample'].map(combined_names)
sample_sheet.to_csv('samplesheet_all.csv', index=False)
print("Done")
