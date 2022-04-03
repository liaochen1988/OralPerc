import pandas as pd
import os
from Bio.Seq import Seq
from collections import Counter

def convert_strings_to_fasta(
    query_sequence_table_path:str,  # path of query sequences listed in a two-column table (ASV ID, ASV sequence)
    query_sequence_fasta_path:str   # path of query sequences in fasta format
)->None:
    df_query = pd.read_csv(query_sequence_table_path)
    df_query.columns = ['USER_ASV','Sequence']
    fasta_query = open(query_sequence_fasta_path, "w")
    for asv_id, asv_seq in zip(df_query.USER_ASV, df_query.Sequence):
        fasta_query.write(">%s\n%s\n"%(asv_id, asv_seq))
    fasta_query.close()
    return None

def run_blast(
    query_sequence_fasta_path:str,       # path of fasta files of query sequences
    blast_db_path:str="../data/HMPv35/", # path of source database that query sequences are blasted to
    evalue_cutoff:float=1e-10,           # evalue cutoff
    num_of_threads:int=1                 # number of threads for blast
)->None:
    blast_output_columns = "qseqid qseq qlen sseqid sseq slen sstrand pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    blast_command = "blastn -db %s -query %s -out unfiltered_blast_output.txt -outfmt \"7 %s\" -perc_identity 100 -evalue %2.2e -num_threads %d"(blast_db_path,query_sequence_fasta_path,blast_output_columns,evalue_cutoff,num_of_threads)
    ok = os.system(blast_command)
    if ok!=0: # program fails
        raise ValueError("command <%s> fails."%(blast_command))
    return None

def parse_blast(
    hmp_asvs_dict: dict,   # keys are HMP ASV IDs, values are corresponding sequences
    user_asvs_dict: dict,  # keys are user-provided ASV IDs, values are corresponding sequences
    output_folder_path:str,     # output folder directory
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
    df_blast_filtered:pd.Dataframe,                 # filtered blast output
    user_asv_count_table_path:str,              # file path of count table of user provided ASVs
    hmp_asv_count_table_path:str,               # file path of count table of HMP ASVs
    is_user_asv_count_table_melted:bool=False,      # is user ASV table melted or expanded?
    is_hmp_asv_count_table_melted:bool=False,       # is HMP ASV table melted or expanded
    verbose:bool=True                               # whether to print intermediate results during running
)->pd.DataFrame:
    # read user provided ASV count table (sample by ASV)
    df_user_asv_count = pd.read_csv(user_asv_count_table_path)
    if not is_user_asv_count_table_melted:
        # convert the table to the melted form
        df_user_asv_count = df_user_asv_count.stack().reset_index()

    # rename columns
    df_user_asv_count.columns = ['SampleID','USER_ASV','Count']

    # rename sample ID to avoid overlaps with HMP samples
    df_user_asv_count.Sample = ['USER__'+sid for sid in df_user_asv_count.Sample]

    # remove ASVs with zero count
    df_user_asv_count = df_user_asv_count[df_user_asv_count.Count>0]
    if verbose:
        print("user provided ASV table contains %d unique ASVs in %d unique samples."%(len(set(df_user_asv_count.ASV)),len(set(df_user_asv_count.Sample))))

    # select and rename filtered blast output
    df2_blast_filtered = deepcopy(df_blast_filtered[['qseqid','sseqid']])
    df2_blast_filtered.columns = ['USER_ASV','HMP_ASV']

    # append a column of matched HMP ASVs
    df_user_asv_count = pd.merge(df_user_asv_count, df2_blast_filtered, left_on='USER_ASV', right_on='USER_ASV', how='left')

    # rename user provided ASVs that are not matched to any HMP sequence by "UMP__"(unmapped) + asv id
    df_user_asv_count.loc[df_user_asv_count.HMP_ASV.isnull(), 'HMP_ASV'] = ['UMP__'+asv for asv in df_user_asv_count.loc[df_user_asv_count.HMP_ASV.isnull(), 'ASV']]

    # read HMP ASV table
    df_hmp_asv_count = pd.read_csv("../data/HMPv35/feature_table.txt", sep="\t", low_memory=False, index_col=0) # ASV by sample
    df_hmp_asv_count = df_hmp_asv_count.T.stack().reset_index()
    df_hmp_asv_count.columns = ['SampleID','HMP_ASV','Count']
    df_hmp_asv_count = df_hmp_asv_count[df_hmp_asv_count.Count>0]

    # concatenate the two tables
    df_joined_asv_count = pd.concat([df_user_asv_count[['SampleID','HMP_ASV','Count']], df_hmp_asv_count], axis=0)

    # save data to file
    df_joined_asv_count.to_csv("tblASVcounts_hmp_user_joined.csv", index=False)

    return df_joined_asv_count

def build_per_asv_random_forest_model(
    df_source_meta:pd.DataFrame,    # meta data of source database
    df_source_relab:pd.DataFrame,   # ASV relative abundance of source database
    num_of_threads:int=1            # number of treads
)->pd.DataFrame:
    
    # grids of relative abundance for each ASV
    # we use random forest model to predict probability of body sites an ASV belongs to at each relative abundance
    abundance_grids = np.linspace(0,1,1001)

    # build random forest model for each source ASV
    def build_rf(SourceASV):
        clf = RandomForestClassifier(n_estimators=10000, random_state=0)
        clf.fit(np.array(df_source_relab[SourceASV]).reshape(-1,1),  df_source_meta['env'])
        probs_on_grids = np.transpose(clf.predict_proba(abundance_grids.reshape(-1,1)))
        
        probs_int = []
        for bs in list(set(df_source_meta.env)):
            idx = list(clf.classes_).index(bs)
            probs_int.append([SourceASV, bs, simpson(probs_on_grids[idx,:], abundance_grids)])
        df_probs_int = pd.DataFrame(probs_int, columns=['SourceASV','BodySite','Probability'])
        return df_probs_int
    random_forest_models = Parallel(n_jobs=-1)(delayed(build_rf)(hmp_asv) for hmp_asv in mapped_hmp_asvs)
    df_rf = pd.concat(random_forest_models)
    df_rf = pd.merge(df_asv_mapping, df_rf, left_on='HMP_ASV', right_on='HMP_ASV', how='left').drop_duplicates()
        return None
    
    
def estimate_oral_bacteria_fraction():
    return None

