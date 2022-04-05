import pandas
import os
from Bio.Seq import Seq
from collections import Counter

blast_output_columns = "qseqid qseq qlen sseqid sseq slen sstrand pident length mismatch gapopen qstart qend sstart send evalue bitscore"

# blast

os.system("blastn -db ../../HMPv35/seqs_v35 -query ../MSK_ASV_ALL.fasta -out matched_hmp_asvs.blast -outfmt \"7 qseqid qseq qlen sseqid sseq slen sstrand pident length mismatch gapopen qstart qend sstart send evalue bitscore\" -perc_identity 100 -evalue 1e-10 -num_threads 12")

def parse_blast(
    hmp_asvs_dict: dict,   # keys are HMP ASV IDs, values are corresponding sequences
    user_asvs_dict: dict,  # keys are user-provided ASV IDs, values are corresponding sequences
    output_folder:str,     # output folder directory
    verbose:bool=True      # whether to print intermediate results to screen during data processing
)->pd.DataFrame:
    #------------------
    # read blast output
    #------------------
    df_blast = pd.read_csv('unfiltered_blast_output.txt', sep='\t', engine='python', header=None, comment='#', names=blast_output_columns.split(' '))
    
    #------------------------------------------------------------------------------------------------------------
    # filter blast output
    # we consider two sequences are matched when no mismatches are found for sequences outside the aligned region
    #------------------------------------------------------------------------------------------------------------
    rowindex2keep = [] # which rows of blast output to keep
    for idx in df_blast.index:
        
        # get alignment information for each pair of HMP and user provided ASV
        qseqid = df.loc[idx,'qseqid']   # query sequence id
        sseqid = df.loc[idx,'sseqid']   # subject sequence id
        sstrand = df.loc[idx,'sstrand'] # subject sequence strand
        qseq = df.loc[idx,'qseq']       # matched region query sequence
        sseq = df.loc[idx,'sseq']       # matched region from subject sequence
        
        # since we require percent identity = 100, the aligned sequences must be equal
        if qseq != sseq:
            raise ValueError("aligned sequences between %s and %s are not equal."%(qseqid, sseqid))

        # get full length sequences of HMP ASVs and user provided ASVs
        if sseqid in hmp_asvs_dict.keys():
            curr_hmp_asv = hmp_asvs_dict[sseqid]
            if sstrand == "minus":
                curr_hmp_asv = str(Seq(curr_hmp_asv).reverse_complement())
        else:
            raise ValueError("could not find %s in HMP asvs dictionary."%(sseqid))

        if qseqid in user_asvs_dict.keys():
            curr_user_asv = user_asvs_dict[qseqid]
        else:
            raise ValueError("could not find %s in user defined asvs dictionary."%(qseqid))

        # make sure that aligned sequences are part of the full-legnth sequences
        if sseq in curr_hmp_asv:
            curr_hmp_asv_split = curr_hmp_asv.split(sseq)
            if len(curr_hmp_asv_split) != 2:
                raise ValueError("%d unmatched regions of %s were found. should be 2."%(len(curr_hmp_asv_split), sseqid))
        else:
            raise ValueError("aligned sequence should be a substring of the full-length sequence of %s."%(sseqid))
        
        if qseq in curr_user_asv:
            curr_user_asv_split = curr_user_asv.split(qseq)
            if len(curr_user_asv_split) != 2:
                raise ValueError("%d unmatched regions of %s were found. should be 2."%(len(curr_user_asv_split), qseqid))
        else:
            raise ValueError("aligned sequence should be a substring of the full-length sequence of %s."%(qseqid))

        # check for any unmatched sequences outside the aligned region
        if curr_user_asv_split[0] != '' and curr_hmp_asv_split[0] != '' and curr_user_asv_split[0] != curr_hmp_asv_split[0]:
            continue
        if curr_user_asv_split[1] != '' and curr_hmp_asv_split[1] != '' and curr_user_asv_split[1] != curr_hmp_asv_split[1]:
            continue

        rowindex2keep.append(idx)
        
    # filter blast output
    df_blast_filtered = df_blast.loc[rowindex2keep]

    # pick the match with highest overlap length if multiple matches are found
    user_asvs_w_multiple_matches = [k for k,v in dict(Counter(df_blast_filtered.qseqid)).items() if v>1]
    rowindex2remove = []
    for curr_user_asv in user_asvs_w_multiple_matches:
        rowindex2remove.extend(list(df_blast_filtered[df_blast_filtered.qseqid==user_asv].sort_values('length', ascending=False).index)[1:])
    df_blast_filtered = df_blast_filtered.loc[[idx for idx in df_blast_filtered.index if idx not in rowindex2remove]]

    # print number of matched sequences
    if verbose:
        print("%d unique user provided sequences are matched to HMP sequences."%(len(set(df_blast_filtered.qseqid))))
        
    # save data to file
    df_blast_filtered.to_csv("filtered_blast_output.txt", sep="\t", index=False)
    
    return df_blast_filtered
    
def merge_hmp_user_asv_count_tables(
    df_filtered_blast_output:pd.Dataframe, # filtered blast output
    filepath_user_asv_count_table:pd.DataFrame, #
    filepath_hmp_asv_count_table:pd.DataFrame, #
    is_user_asv_count_table_melted:bool,
    is_hmp_asv_count_table_melted:bool
):
    # get mapping between HMP and MSK
    df_blast_filtered = pd.read_csv("filtered_blast_results.csv")
    df_asv_mapping = df_blast_filtered[['qseqid','sseqid']]
    df_asv_mapping.columns = ['ASV','HMP_ASV']

    # transform MSK ASV to HMP ASV
    df_msk_count = pd.read_csv("../../MSKCC_alloHCT/tblcounts_asv_melt.csv")
    df_msk_count = pd.merge(df_msk_count, df_asv_mapping, left_on='ASV', right_on='ASV', how='left')
    df_msk_count.loc[df_msk_count.HMP_ASV.isnull(), 'HMP_ASV'] = ['unmapped_'+asv for asv in df_msk_count.loc[df_msk_count.HMP_ASV.isnull(), 'ASV']]

    # read HMP ASV table
    df_hmp_count = pd.read_csv("../../HMPv35/feature_table.txt", sep="\t", low_memory=False, index_col=0)
    df_hmp_count_stacked = df_hmp_count.T.stack().reset_index()
    df_hmp_count_stacked.columns = ['SampleID','HMP_ASV','Count']
    df_hmp_count_stacked = df_hmp_count_stacked[df_hmp_count_stacked.Count>0]

    # concat the two tables
    df_join = pd.concat([df_msk_count[['SampleID','HMP_ASV','Count']], df_hmp_count_stacked], axis=0)
    df_join.to_csv("tblcounts_asv_melt_hmp_msk_joined.csv", index=False)
    df_join.head()

