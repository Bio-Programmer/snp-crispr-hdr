# Title: crispr-to-motifbreakR.py
# Description: motifbreakR analysis script for processing snp-crispr-hdr.py output
# Author: Sophia Longo
# Date: 4 August 2021

import sys
import pandas as pd

# Get input files
mbr_input_filepath = sys.argv[1]
# mbr_input_filepath = '/Users/sophia/Desktop/MHC_ARVID/hdr/snp-crispr-hdr/data/motifbreakR_motif_hdr_snps_080321_wPval.tsv'
crispr_input_filepath = sys.argv[2]
# crispr_input_filepath = '/Users/sophia/Desktop/MHC_ARVID/hdr/snp-crispr-hdr/data/snp-crispr-hdr_output_073121.csv'
output_filepath = sys.argv[3]
# output_filepath = '/Users/sophia/Desktop/MHC_ARVID/hdr/snp-crispr-hdr/data/snp-crispr-hdr_motifbreakR_apm_080421.csv'

# (1) Pre-processing

# Create unique motif_ids for each motifbreakR entry
def generate_motif_ids(mbr_df):
    motif_ids_list = []
    for idx, motif_row in mbr_df.iterrows():
        motif_id = str(motif_row['SNP_id']) + '_' + str(motif_row['geneSymbol']) + '_' + str(motif_row['motifPos']) + '_' + str(motif_row['dataSource'])
        motif_ids_list.append(motif_id)
    mbr_df.insert(loc=0, column='motif_id', value=motif_ids_list)
    return mbr_df

# Make it so that there are proper geneSymbols for every motif entry
def clean_mbr_geneSymbols(mbr_df):
    missing_geneSymbols_dict = {}
    for idx, motif_row in mbr_df.iterrows():
        if str(motif_row['geneSymbol']) == 'nan':
            motif_id = motif_row['motif_id']
            missing_geneSymbols_dict[motif_id] = str(motif_row['providerId'])

    all_geneSymbols = []
    for idx, motif_row in mbr_df.iterrows():
        if str(motif_row['geneSymbol']) == 'nan':
            motif_id = motif_row['motif_id']
            all_geneSymbols.append(missing_geneSymbols_dict[motif_id])
        else:
            all_geneSymbols.append(motif_row['geneSymbol'])

    mbr_df['geneSymbol'] = all_geneSymbols
    return mbr_df

# 2. Collect motifbreakR data in the format of a dictionary
def assemble_mbr_dict(mbr_df):
    mbr_snps = list(set(mbr_df['SNP_id']))
    all_snps_mbr_dict = {}
    for snp in mbr_snps:
        snp_mbr_df = mbr_df.loc[mbr_df['SNP_id'] == snp]
        snp_mbr_tfs = list(set(snp_mbr_df['geneSymbol']))
        snp_mbr_sources = list(set(snp_mbr_df['dataSource']))
        snp_mbr_motifIds = list(set(snp_mbr_df['motif_id']))
        snp_tf_alleleDiff_list = []
        for tf in snp_mbr_tfs:
            snp_tf_mbr_df = snp_mbr_df.loc[snp_mbr_df['geneSymbol'] == tf]
            snp_tf_alleleDiff = list(snp_tf_mbr_df['alleleDiff'])
            snp_tf_avg_alleleDiff = sum(snp_tf_alleleDiff) / len(snp_tf_alleleDiff)
            tf_alleleDiff = tf + ';' + str(round(snp_tf_avg_alleleDiff, 2))
            snp_tf_alleleDiff_list.append(tf_alleleDiff)
        snp_tf_alleleDiff_str = '|'.join(snp_tf_alleleDiff_list)
        snp_mbr_sources_str = ';'.join(snp_mbr_sources)
        snp_mbr_motifIds_str = ';'.join(snp_mbr_motifIds)
        all_snps_mbr_dict[snp] = {'tf_alleleDiffAvg':snp_tf_alleleDiff_str,
                                  'data_sources':snp_mbr_sources_str,
                                  'snp_mbr_motifIds':snp_mbr_motifIds_str}
    return all_snps_mbr_dict

# 3. Leverage the motifbreakR dictionary of data to add motifbreakR data to the crispr_df
def attach_mbr_data(crispr_df, mbr_data_dict):
    allSNPs_tf_alleleDiffAvg_list = []
    allSNPs_tf_sources_list = []
    allSNPs_motifIDs_list = []
    for idx, snp_row in crispr_df.iterrows():
        snp = snp_row['snp']
        if snp in mbr_data_dict:
            allSNPs_tf_alleleDiffAvg_list.append(mbr_data_dict[snp]['tf_alleleDiffAvg'])
            allSNPs_tf_sources_list.append(mbr_data_dict[snp]['data_sources'])
        else:
            allSNPs_tf_alleleDiffAvg_list.append('noneFound')
            allSNPs_tf_sources_list.append('NA')
            allSNPs_motifIDs_list.append('NA')
    crispr_df['mbr_tf_alleDiffAvg'] = allSNPs_tf_alleleDiffAvg_list
    crispr_df['mbr_dataSources'] = allSNPs_tf_sources_list
    return crispr_df

def execute(mbr_input_filepath, crispr_input_filepath, output_filepath):
    mbr_input_df = pd.read_csv(mbr_input_filepath, sep='\t')
    crispr_input_df = pd.read_csv(crispr_input_filepath)
    mbr_clean = clean_mbr_geneSymbols(generate_motif_ids(mbr_input_df))
    print('Successfully processed your motifbreakR input file.')
    mbr_data_dict = assemble_mbr_dict(mbr_clean)
    print('Successfully extracted data from your motifbreakR input file.')
    crispr_mbr_df = attach_mbr_data(crispr_input_df, mbr_data_dict)
    crispr_mbr_df.to_csv(output_filepath, index=False)
    print('Congratulations! Successfully annotated your snp-crispr-hdr.py file with your motifbreakR data.')
    print('Your annotated file has been written out to: ', output_filepath)

execute(mbr_input_filepath, crispr_input_filepath, output_filepath)
