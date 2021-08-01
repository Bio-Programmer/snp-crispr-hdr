# snp-crispr-hdr.py
# Purpose: To obtain SNPs that fulfill CRISPR Cas9 and CRISPR Cas12 conditions for homology directed repair experiments
# Input: .csv containing a table with columns: mpra_experiment, snp, chr, pos, disease, eGenes
# Output: a table containing SNPs that satisfy 1 or more of the CRISPR Cas9 and Cas12 conditions
# Author: Sophia Longo, 31 July 2021

import sys
import pandas as pd
import urllib.request
import xmltodict
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import requests

#################################################
##### (0) User input and output file paths #####
###############################################
# mpra_input_path = '/Users/sophia/Desktop/MHC_ARVID/hdr/data/mpraHits_to_crisprFilter_073021.csv'
mpra_input_path = sys.argv[1] # locally: /Users/sophia/Desktop/MHC_ARVID/hdr/data/mpraHits_to_crisprFilter_073021.csv
# crispr_output_path = '/Users/sophia/Desktop/MHC_ARVID/hdr/snp-crispr-hdr_output_test.csv'
crispr_output_path = sys.argv[2]
mpra_input_df = pd.read_csv(mpra_input_path)

################################################################
##### (1) Obtain nucleotide sequences of surrounding SNPs #####
################################################################

def prompt_padding():
    while True:
        pad_len = int(input('How many base pairs would you like surrounding your SNPs? (25 recommended): '))
        if pad_len > 0:
            break
        continue
    return pad_len

# prepare to obtain sequences by sepcifying amount of paddding surrounding each SNP
def add_padding(mpra_df, pad_bp):
    mpra_df['pad_start'] = mpra_df['pos'].map(int) - pad_bp
    mpra_df['pad_end'] = mpra_df['pos'].map(int) + pad_bp
    mpra_df = mpra_df[['mpra_experiment', 'snp', 'chr', 'pos',
                       'pad_start', 'pad_end', 'MHC_eGenes', 'disease']]
    return(mpra_df)

# obtain sequences that you want through API calls to UCSC genome browser
# IMPORTANT: make sure all positions are in desired HG version (default is hg38)
# hgLiftOVer: Tool to convert to correct version: https://genome.ucsc.edu/cgi-bin/hgLiftOver
def get_sequence(row, chr_col='chr', start_col='pad_start', stop_col='pad_end', version='hg38'):
    # use ucsc's DAS server to scrape sequences
    url = 'http://genome.ucsc.edu/cgi-bin/das/{}/dna?segment={}:{},{}'.format(version,
                                                                              row[chr_col],
                                                                              row[start_col],
                                                                              row[stop_col])
    # read the data in from the desired URL
    contents = urllib.request.urlopen(url).read()
    # convert from XML to ordered dict and index down to part of ordered dict that contains sequence
    # sequence might go over multiple lines, so join the lines together
    sequence = ''.join(xmltodict.parse(contents)['DASDNA']['SEQUENCE']['DNA']['#text'].upper().split('\n'))
    return(sequence)

# obtain dataframe of desired sequences
def assemble_sequence_df(mpra_df, v_hg = 'hg38'):
    mpra_df['sequence'] = mpra_df.apply(lambda x:get_sequence(x, version=v_hg, chr_col='chr',
                                                              start_col='pad_start', stop_col='pad_end'), axis=1)
    return(mpra_df)

# demarcate SNP within sequence for readability
def star_snp_seqs(mpra_seq_df):
    seqs_starred_list = []
    for idx, snp_row in mpra_seq_df.iterrows():
        seq = snp_row['sequence']
        seq = seq[:25] + '*' + seq[25:]
        seq_starred = seq[:27] + '*' + seq[27:]
        seqs_starred_list.append(seq_starred)
    mpra_seq_df['sequence'] = seqs_starred_list
    return(mpra_seq_df)

def prompt_for_allele_info(rsid):
    print('SNP not found in GTEx: ', rsid)
    nucleotides = ['A', 'T', 'C', 'G']
    while True:
        ref = input('Input ' + rsid + ' ref allele: ').upper()
        if ref in nucleotides:
            break
        print('Invalid input, please try again.')
        continue
    while True:
        alt = input('Input ' + rsid + ' alternate allele: ').upper()
        if alt in nucleotides:
            break
        print('Invalid input, please try again.')
        continue
    while True:
        is_maf01 = input('Input ' + rsid + ' if minor allele frequency > 1% (TRUE/FALSE): ')
        if is_maf01 == 'TRUE' or is_maf01 == 'FALSE':
            break
        print('Invalid input, please try again.')
        continue
    return {'ref':ref, 'alt':alt, 'is_maf01':is_maf01}

# obtain allele info from GTEx (reference, alternate and whether minor allele frequency is > 1% (is_maf01))
def get_allele_info(rsid):
    # rsid = 'rs12153855'
    url = 'https://gtexportal.org/rest/v1/dataset/variant?format=json&snpId={}&datasetId=gtex_v8'.format(rsid)
    response = requests.get(url)
    if response:
        snp_json = response.json()
        variant_info = snp_json['variant']
        if len(variant_info) > 0:
            variant_info = snp_json['variant'][0]
            ref = variant_info['ref']
            alt = variant_info['alt']
            is_maf01 = variant_info['maf01']
        else:
            return prompt_for_allele_info(rsid)
    else:
        return prompt_for_allele_info(rsid)
    return {'ref':ref, 'alt':alt, 'is_maf01':is_maf01}

def assemble_allele_info(mpra_seq_df):
    snp_list = list(mpra_seq_df['snp'])
    ref_list = []
    alt_list = []
    is_maf01_list = []
    for snp in snp_list:
        allele_info = get_allele_info(snp)
        ref_list.append(allele_info['ref'])
        alt_list.append(allele_info['alt'])
        is_maf01_list.append(allele_info['is_maf01'])
    mpra_seq_df['ref'] = ref_list
    mpra_seq_df['alt'] = alt_list
    mpra_seq_df['is_maf01'] = is_maf01_list
    return mpra_seq_df

##############################################################################
##### (2) Apply CRISPR filters to MPRA dataframe that contains sequences #####
#############################################################################

def attach_reverse_complement_to_df(seq_df):
    rev_comps_list = []
    for idx, snp in seq_df.iterrows():
        seq = Seq(snp['sequence'])
        reverse_comp_seq = str(seq.reverse_complement())
        rev_comps_list.append(reverse_comp_seq)
    # add column next to seq column
    # seq_df = seq_df.insert(loc=seq_df.columns.get_loc('sequence')+1, column='reverse_complement', value=rev_comps_list)
    seq_df['reverse_complement'] = rev_comps_list
    return seq_df

def is_in_cas9_PAM(sequence, alt_allele_nt):
    # cas9 PAM (5' to 3') = NGG
    snp_index = sequence.index('*') + 1
    two_nt_5to3 = sequence[snp_index+2:snp_index+4]
    potential_cutReactivator = sequence[snp_index+4]
    if two_nt_5to3 == 'GG' and alt_allele_nt != 'G' and potential_cutReactivator != 'G':
        return True
    return False

def is_in_cas12_PAM(sequence, alt_allele_nt):
    # cas12 PAM (5' to 3') = VTTT
    snp_index = sequence.index('*') + 1
    snp_nt = sequence[snp_index]
    three_nt_5to3 = sequence[snp_index+2:snp_index+5]
    potential_cutReactivator = sequence[snp_index+5]
    if three_nt_5to3 == 'TTT' and snp_nt != 'T' and alt_allele_nt != 'T' and potential_cutReactivator != 'T':
        return True
    return False

def is_in_cas9_gRNA_range(snp_seq, gRNA_range=8):
    lower_nt_upstream_of_PAM = 1
    upper_nt_upstream_of_PAM = gRNA_range
    snp_index = snp_seq.index('*') + 1
    nt_frame = ''
    for i in range(lower_nt_upstream_of_PAM, upper_nt_upstream_of_PAM):
        nt_frame = snp_seq[snp_index+i:snp_index+3+i]
        if nt_frame[1:3] == 'GG':
            return True
    return False

def is_in_cas12_gRNA_range(snp_seq, lower_nt_downstreamOfPAM=3, upper_nt_downstreamOfPAM=16):
    snp_index = snp_seq.index('*') + 1
    nt_frame = ''
    for i in range(lower_nt_downstreamOfPAM, upper_nt_downstreamOfPAM):
        nt_frame = snp_seq[snp_index-4-i:snp_index-i]
        if nt_frame[0:3] == 'TTT':
            return True
    return False

def get_is_candidate_bools(seq_df):
    cas9_PAM_bool_list = []
    cas12_PAM_bool_list = []
    cas9_gRNA_bool_list = []
    cas12_gRNA_bool_list = []
    is_candidate_bool_list = []
    gc_content_list = []

    # obtain candidates for each condition
    for idx, snp in seq_df.iterrows():
        seq = snp['sequence']
        alt_allele = snp['alt']
        alt_allele_comp = str(Seq(alt_allele).reverse_complement())
        rev_com_seq = snp['reverse_complement']
        is_candidate = False
        # compute GC content
        gc_content = GC(Seq(seq))
        gc_content_list.append(gc_content)
        # PAM sequence filters
        if is_in_cas9_PAM(seq, alt_allele) == True or is_in_cas9_PAM(rev_com_seq, alt_allele_comp):
            cas9_PAM_bool_list.append(True)
            is_candidate = True
        else:
            cas9_PAM_bool_list.append(False)
        if is_in_cas12_PAM(seq, alt_allele) == True or is_in_cas12_PAM(rev_com_seq, alt_allele_comp):
            cas12_PAM_bool_list.append(True)
            is_candidate = True
        else:
            cas12_PAM_bool_list.append(False)
        # gRNA sequence filters
        if is_in_cas9_gRNA_range(seq) == True or is_in_cas9_gRNA_range(rev_com_seq):
            cas9_gRNA_bool_list.append(True)
            is_candidate = True
        else:
            cas9_gRNA_bool_list.append(False)
        if is_in_cas12_gRNA_range(seq) == True or is_in_cas12_gRNA_range(rev_com_seq):
            cas12_gRNA_bool_list.append(True)
            is_candidate = True
        else:
            cas12_gRNA_bool_list.append(False)
        # mark SNP as HDR candidate if it passes any of the above filters
        if is_candidate == True:
            is_candidate_bool_list.append(True)
        else:
            is_candidate_bool_list.append(False)

    all_conditions_bool_dict = {'cas9_PAM':cas9_PAM_bool_list, 'cas12_PAM':cas12_PAM_bool_list,
                                'cas9_gRNA':cas9_gRNA_bool_list, 'cas12_gRNA':cas12_gRNA_bool_list,
                                'is_candidate':is_candidate_bool_list, 'gc_content':gc_content_list}
    return all_conditions_bool_dict

def percent_gc_content_filter(hdr_df, lower, upper):
    hdr_df = hdr_df[hdr_df['percent_GC'] >= lower]
    hdr_df = hdr_df[hdr_df['percent_GC'] <= upper]
    return hdr_df

def prompt_gc_content_filter(hdr_df):
    print('You have narrowed down your candidates to: ', hdr_df.shape[0])
    while True:
        gc_filter = input('Do you wish to specify a % GC conent filter (Y/N): ')
        if gc_filter == 'Y':
            while True:
                lower = float(input('Please specify lower bounds (0-99): '))
                upper = float(input('Please specify upper bounds (1-100): '))
                if lower > upper:
                    print('Invalid specification: lower cannot be greater than upper.')
                    continue
                break
            hdr_df = percent_gc_content_filter(hdr_df, lower, upper)
            print('Table filtered successfully by % GC content.')
            break
        elif gc_filter == 'N':
            break
        else:
            print('Invalid input, please try again.')
            continue
    return hdr_df

def get_all_candidates(seq_df):
    seq_df = attach_reverse_complement_to_df(seq_df)
    original_num_candidates = seq_df.shape[0]
    all_conditions_bool_dict = get_is_candidate_bools(seq_df)
    # append candidate booleans to seq_df
    seq_df['percent_GC'] = all_conditions_bool_dict['gc_content']
    seq_df['is_cas9_PAM'] = all_conditions_bool_dict['cas9_PAM']
    seq_df['is_cas12_PAM'] = all_conditions_bool_dict['cas12_PAM']
    seq_df['is_cas9_gRNA'] = all_conditions_bool_dict['cas9_gRNA']
    seq_df['is_cas12_gRNA'] = all_conditions_bool_dict['cas12_gRNA']
    seq_df['is_candidate'] = all_conditions_bool_dict['is_candidate']
    # sort by GC content
    seq_df = seq_df.sort_values(by='percent_GC')
    # keep only SNPs that are candidates
    seq_df = seq_df[seq_df.is_candidate == True]
    seq_df = seq_df.drop(columns='is_candidate')
    seq_df = seq_df[['mpra_experiment', 'disease', 'snp', 'MHC_eGenes',
                     'chr', 'pos', 'pad_start', 'pad_end', 'ref', 'alt', 'is_maf01', 'sequence', 'reverse_complement',
                     'percent_GC', 'is_cas9_PAM', 'is_cas12_PAM', 'is_cas9_gRNA', 'is_cas12_gRNA']]
    print('Original number of candidates: ', original_num_candidates)
    return seq_df


##########################################################################
##### (3) Execute: Put all of the above steps and functions together #####
##########################################################################

def execute(mpra_input_df, crispr_output_path):
    print('Welcome to snp-crispr-hdr!')
    sequence_df = star_snp_seqs(assemble_sequence_df(add_padding(mpra_input_df, prompt_padding())))
    print('Success: Obtained sequences surrounding your SNPs.')
    sequence_df_wAllele_info = assemble_allele_info(sequence_df)
    print('Success: Obtained allele information from GTEx.')
    hdr_df = get_all_candidates(sequence_df_wAllele_info)
    print('Success: Applied CRISPR filter.')
    hdr_df = prompt_gc_content_filter(hdr_df)
    hdr_df.to_csv(crispr_output_path, index=False)
    print('Congratulations! HDR candidates table written out to: ', crispr_output_path)

execute(mpra_input_df, crispr_output_path)
