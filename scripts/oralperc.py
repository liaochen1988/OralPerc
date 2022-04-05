import pandas as pd
import numpy as np
import sys, getopt, os
from Bio import SeqIO
from collections import Counter
from copy import deepcopy

# blast output columns
blast_output_columns = "qseqid qseq qlen sseqid sseq slen sstrand pident length mismatch gapopen qstart qend sstart send evalue bitscore"

def blast_against_source(
    query_sequence_fasta:str,   # path of fasta files of query sequences
    source:str,                 # name of source
    source_dir:str,            # path of source database that query sequences are blasted to
    evalue_cutoff:float,        # evalue cutoff
    num_of_threads:int,         # number of threads for blast
    output_dir:str              # directory of output folder
)->None:
    blast_command = "blastn -db %s/%s -query %s -out %s/unfiltered_blast_against_source_%s.txt -outfmt \"7 %s\" -perc_identity 100 -evalue %2.2e -num_threads %d"%(
        source_dir, source, query_sequence_fasta, output_dir, source, blast_output_columns, evalue_cutoff, num_of_threads)
    ok = os.system(blast_command)
    if ok!=0: # program fails
        raise RuntimeError("command <%s> fails."%(blast_command))
    return None

def parse_blast(
    asv_dict_source:dict,     # keys are source ASV IDs, values are sequences
    asv_dict_query:dict,      # keys are query ASV IDs, values are sequences
    source:str,               # source ID
    output_dir:str            # output directory
)->pd.DataFrame:
    #------------------
    # read blast output
    #------------------
    df_blast = pd.read_csv('%s/unfiltered_blast_against_source_%s.txt'%(output_dir, source), sep='\t', engine='python', header=None, comment='#', names=blast_output_columns.split(' '))

    #------------------------------------------------------------------------------------------------------------
    # filter blast output
    # we consider two sequences are matched when no mismatches are found for sequences outside the aligned region
    #------------------------------------------------------------------------------------------------------------
    rowindex2keep = [] # which rows of blast output to keep
    for idx in df_blast.index:

        # get alignment information for each pair of HMP and user provided ASV
        qseqid = df_blast.loc[idx,'qseqid']   # query sequence id
        sseqid = df_blast.loc[idx,'sseqid']   # subject sequence id
        sstrand = df_blast.loc[idx,'sstrand'] # subject sequence strand
        qseq = df_blast.loc[idx,'qseq']       # matched region query sequence
        sseq = df_blast.loc[idx,'sseq']       # matched region from subject sequence

        # since we require percent identity = 100, the aligned sequences must be equal
        if qseq != sseq:
            raise RuntimeError("aligned sequences between %s and %s are not equal."%(qseqid, sseqid))

        # get full length sequences of source and query ASVs
        if sseqid in asv_dict_source.keys():
            curr_source_asv = asv_dict_source[sseqid]
            if sstrand == "minus":
                curr_source_asv = str(Seq(curr_source_asv).reverse_complement())
        else:
            raise RuntimeError("could not find %s in source asv dictionary."%(sseqid))

        if qseqid in asv_dict_query.keys():
            curr_user_asv = asv_dict_query[qseqid]
        else:
            raise RuntimeError("could not find %s in query asv dictionary."%(qseqid))

        # make sure that aligned sequences are part of the full-legnth sequences
        if sseq in curr_source_asv:
            curr_source_asv_split = curr_source_asv.split(sseq)
            if len(curr_source_asv_split) != 2:
                raise RuntimeError("%d unmatched regions of %s were found. should be 2."%(len(curr_source_asv_split), sseqid))
        else:
            raise RuntimeError("aligned sequence should be a substring of the full-length sequence of %s."%(sseqid))

        if qseq in curr_user_asv:
            curr_user_asv_split = curr_user_asv.split(qseq)
            if len(curr_user_asv_split) != 2:
                raise RuntimeError("%d unmatched regions of %s were found. should be 2."%(len(curr_user_asv_split), qseqid))
        else:
            raise RuntimeError("aligned sequence should be a substring of the full-length sequence of %s."%(qseqid))

        # check for any unmatched sequences outside the aligned region
        if curr_user_asv_split[0] != '' and curr_source_asv_split[0] != '' and curr_user_asv_split[0] != curr_source_asv_split[0]:
            continue
        if curr_user_asv_split[1] != '' and curr_source_asv_split[1] != '' and curr_user_asv_split[1] != curr_source_asv_split[1]:
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
    print("%d unique query sequences are matched to source sequences."%(len(set(df_blast_filtered.qseqid))))

    # save data to file
    df_blast_filtered.to_csv("%s/filtered_blast_output_for_source_%s.txt"%(output_dir, source), sep="\t", index=False)

    return df_blast_filtered

def merge_query_source_count_tables(
    df_blast_filtered:pd.DataFrame,     # filtered blast output
    df_query_asv_count:pd.DataFrame,    # count table of query ASVs
    df_source_asv_count:pd.DataFrame,   # count table of source ASVs
    source:str,                         # source ID
    output_dir:str                      # outpur directory
)->pd.DataFrame:

    # Important Notes:
    # Both query and source count tables are passed in long format (not wide format).

    # rename query sample ID to avoid overlaps with source sample IDs
    df_query_asv_count = df_query_asv_count[['SampleID','ASV','Count']].rename({"ASV":"QueryASV"}, axis=1)

    # select and rename filtered blast output
    df2_blast_filtered = deepcopy(df_blast_filtered[['qseqid','sseqid']])
    df2_blast_filtered.columns = ['QueryASV','SourceASV']

    # append a column of matched source ASVs in query count table
    df_query_asv_count = pd.merge(df_query_asv_count, df2_blast_filtered, left_on='QueryASV', right_on='QueryASV', how='left')

    # rename query ASVs that are not matched to any source ASV sequence by adding prefix "Unmapped__"
    df_query_asv_count.loc[df_query_asv_count.SourceASV.isnull(), 'SourceASV'] = ['Unmapped__'+asv for asv in df_query_asv_count.loc[df_query_asv_count.SourceASV.isnull(), 'QueryASV']]

    # concatenate the two tables
    df_source_asv_count = df_source_asv_count.rename({'ASV':'SourceASV'}, axis=1)
    df_joined_asv_count = pd.concat([df_query_asv_count[['SampleID','SourceASV','Count']], df_source_asv_count[['SampleID','SourceASV','Count']]], axis=0)
    df_joined_asv_count = df_joined_asv_count[df_joined_asv_count.Count > 0] # remove ASVs of zero count
    df_joined_asv_count.to_csv("%s/joined_asv_count_for_source_%s.csv"%(output_dir, source), sep="\t", index=False)

    return df_joined_asv_count

def build_per_asv_random_forest_model(
    df_source_sample_meta:pd.DataFrame,     # meta data of source
    df_source_asv_count:pd.DataFrame,       # relative abundance of source ASVs (samples by ASVs)
    source:str,                             # source ID
    source_dir:str,                         # path of source dir
    num_of_threads:int                      # number of treads
)->pd.DataFrame:

    # do not recreate the wheel: check if output file has already been generated
    if os.path.exists("%s/rf_bodysite_prob_for_source_%s.txt"%(output_dir, source)):
        df_bodysite_rf_cla_source_asvs = pd.read_csv("%s/rf_bodysite_prob_for_source_%s.txt"%(output_dir, source), sep="\t")
        return df_bodysite_rf_cla_source_asvs
    else:
        # convert count to relative abundance
        # note: df_source_asv_count is passed in long format
        df_source_asv_relab = deepcopy(pd.pivot_table(df_source_asv_count, index="SampleID", columns="ASV", values="Count", aggfunc=np.sum).fillna(0))
        df_source_asv_relab = df_source_asv_relab.div(df_source_asv_relab.sum(axis=1), axis=0)
        df_source_asv_relab = df_source_asv_relab.loc[list(df_source_sample_meta.SampleID)]

    # set up grids of relative abundance for each source ASV
    asv_abundance_grids = np.linspace(0,1,1001)

    # build random forest model for each source ASV
    def build_rf(source_asv):
        # train RF model
        clf = RandomForestClassifier(n_estimators=10000, random_state=0)
        clf.fit(np.array(df_source_relab[source_asv]).reshape(-1,1),  df_source_meta['Environment'])
        probs_on_grids = np.transpose(clf.predict_proba(asv_abundance_grids.reshape(-1,1)))

        # get weights
        nbins = list(asv_abundance_grids-0.0005)+[1.0005]
        curr_asv_pos_dist = list(df_source_relab.loc[df_source_relab[source_asv]>0, source_asv]) # distribution of non-zero relative abundance of the current ASV
        if len(curr_asv_pos_dist)<=1:
            weights = np.array([1]*len(asv_abundance_grids))
        else:
            weights, _ = np.histogram(curr_asv_pos_dist, bins=nbins)
            weights = np.array(weights/simpson(weights, asv_abundance_grids))

        # compute probability of each body site as the integral over abundance grids
        probs_of_body_site = []
        for bs in list(set(df_source_meta.Environment)):
            idx = list(clf.classes_).index(bs)
            integral = simpson(probs_on_grids[idx,:]*weights, asv_abundance_grids)
            probs_of_body_site.append([source_asv, bs, integral])
        df_probs_of_body_site = pd.DataFrame(probs_of_body_site, columns=['SourceASV','BodySite','Prob'])
        return df_probs_of_body_site

    # run in parallel for all source ASVs
    bodysite_probs = Parallel(n_jobs=num_of_threads)(delayed(build_rf)(asv) for asv in df_source_asv_relab.index)
    df_bodysite_rf_cla_source_asvs = pd.concat(bodysite_probs, ignore_index=True)
    df_bodysite_rf_cla_source_asvs.to_csv("%s/rf_bodysite_prob_for_source_%s.txt"%(output_dir, source), sep="\t", index=False)

    return df_bodysite_rf_cla_source_asvs

def predict_oral_percentage(
    df_bodysite_rf_cla_source_asvs:pd.DataFrame,    # body site classification for source ASVs
    df_joined_asv_count:pd.DataFrame,               # joined (both query and source) asv count table
    num_of_threads:int,                             # number of treads
    prevalence_cutoff:float                         # minimum prevalence to classify body site origin from sources
):
    # convert count to relative abundance
    df_joined_asv_relab = deepcopy(pd.pivot_table(df_joined_asv_count, index="SourceASV", columns="SampleID", values="Count", aggfunc=np.sum).fillna(0))
    df_joined_asv_relab = df_joined_asv_relab.div(df_joined_asv_relab.sum(axis=0), axis=1)
    df_query_asv_relab = df_joined_asv_relab[[sid for sid in df_joined_asv_relab.columns if sid.startswith("QUERY__")]].T      # samples by ASVs
    df_source_asv_relab = df_joined_asv_relab[[sid for sid in df_joined_asv_relab.columns if sid.startswith("SOURCE__")]].T    # samples by ASVs

    # get all environment categories
    all_unique_envs = list(set(df_bodysite_rf_cla_source_asvs.BodySite))

    def identify_source_percentages(sid):
        df_curr_sample = deepcopy(df_query_asv_relab.loc[sid])
        df_curr_sample = df_curr_sample[df_curr_sample>0].to_frame()
        res = [sid.lstrip("QUERY__")] + [0]*(len(all_unique_envs)+1)
        for query_asv, query_relab in zip(df_curr_sample.index, df_curr_sample[sid]):
            if query_asv in list(df_bodysite_rf_cla_source_asvs.SourceASV):
                prevalence = len(df_source_asv_relab[df_source_asv_relab[query_asv]>0.0])/len(df_source_asv_relab)
                if prevalence < prevalence_cutoff:
                    # too few samples contains this ASV; set it as unknown source
                    res[-1] += relab
                else:
                    # sufficient number of samples contain this ASV; quantify its body site origin
                    probs_curr_asv = df_bodysite_rf_cla_source_asvs.loc[df_bodysite_rf_cla_source_asvs.SourceASV==query_asv, ['BodySite','Prob']].set_index('BodySite').T
                    for k,bs in enumerate(all_unique_envs):
                        res[k+1] += relab*list(probs_curr_asv[bs])[0]
            else:
                res[-1] += query_relab
        return res

    source_percentages = Parallel(n_jobs=num_of_threads)(delayed(identify_source_percentages)(sid) for sid in df_query_asv_relab.index)
    df_source_percentages = pd.DataFrame(source_proportion, columns=['SampleID']+all_unique_envs+['Unknown']).set_index('SampleID').stack().reset_index()
    df_source_percentages.columns = ['SampleID', 'Environment','Percentage']

    return df_source_percentages

#------------------------
# Main program
#------------------------
if __name__ == "__main__":

    #----------------------------
    # Get command line parameters
    #----------------------------

    # parse command line arguments
    # --usage
    # --query_feature_table <required>: a txt table with ASVs on the rows and samples on the columns
    # --query_asv_sequences <required>: a fasta file of ASV sequences
    # --source_dir: directory of sources. Each source subdirectory has at least three data files: sample_meta, feature_table, and asv_sequences.
    # --which_source: the source will be used against all samples. For heterogenous sources, use source_mapping_file instead. If both which_source and
    #                 source_mapping_file exist, which_source is prefered.
    # --source_mapping_file: heterogenous sources exist and specify which sources are used against samples
    # --inference_method [default RF]: choose between FEAST and RF (Random Forest)
    # --max_feast_runs [default 5]: maximum number of samples run by feast simutaneously (used only when method == feast)
    # --evalue_cutoff [default 1e-10]: evalue cutoff
    # --num_of_threads [default 1]: number of threads
    # --output_dir [default './output']: directory of all output files

    # default parameters
    source_dir = None
    which_source = None
    source_mapping_file = None
    inference_method = "RF"
    max_feast_runs = 5
    evalue_cutoff = 1e-10
    num_of_threads = 5
    prevalence_cutoff = 0.052
    output_dir = "./output"

    # parse command line
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "u:f:a:d:w:s:i:m:e:n:p:o:",["usage","query_feature_table=", "query_asv_sequences=", "source_dir=",
                                                                                 "which_source", "source_mapping_file", "inference_method", "max_feast_runs",
                                                                                 "evalue_cutoff","num_of_threads", "prevalence_cutoff", "output_dir"])
    except getopt.GetoptError as e:
        print('oralperc.py -u -f <query_feature_table> -a <query_asv_sequences> -d <source_dir> -w [which_source] -s [source_mapping_file] '\
              '-i [inference_method] -m [max_feast_runs] -e [evalue_cutoff] -n [num_of_threads] -p [prevalence_cutoff] -o [output_dir]')
        print('>>>> ERROR: %s' % str(e))
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-u':
            print('oralperc.py -u -f <query_feature_table> -a <query_asv_sequences> -d <source_dir> -w [which_source] -s [source_mapping_file] '\
                  '-i [inference_method] -m [max_feast_runs] -e [evalue_cutoff] -n [num_of_threads] -p [prevalence_cutoff] -o [output_dir]')
            sys.exit(0)
        elif opt in ("-f", "--query_feature_table"):
            query_feature_table = arg
        elif opt in ("-a", "--query_asv_sequences"):
            query_asv_sequences = arg
        elif opt in ("-d", "--source_dir"):
            source_dir = arg
        elif opt in ("-w", "which_source"):
            which_source = arg
        elif opt in ("-s", "source_mapping_file"):
            source_mapping_file = arg
        elif opt in ("-i", "--inference_method"):
            inference_method = arg
        elif opt in ("-m", "--max_feast_runs"):
            max_feast_runs = int(arg)
        elif opt in ("-e", "--evalue_cutoff"):
            evalue_cutoff = float(arg)
        elif opt in ("-n", "--num_of_threads"):
            num_of_threads = int(arg)
        elif opt in ("-p", "--prevalence_cutoff"):
            prevalence_cutoff = float(arg)
        elif opt in ("-o", "--output_dir"):
            output_dir = arg

    # check for required files
    if query_feature_table is None:
        raise RuntimeError('query_feature_table was not provided in argument list.')
    if query_asv_sequences is None:
        raise RuntimeError('query_asv_sequences was not provided in argument list.')
    if source_dir is None:
        raise RuntimeError('source_dir was not provided in argument list.')
    if which_source is None and source_mapping_file is None:
        raise RuntimeError('both which_source and source_mapping_file were not provided in argument list.')
    if not os.path.isdir(output_dir):
        os.system("mkdir %s"%(output_dir))

    #----------------------
    # Load query microbiome
    #----------------------

    # read ASV sequences
    obj_query_asv_sequences = SeqIO.parse(open(query_asv_sequences), 'fasta')

    # read feature table
    compression = None
    if query_feature_table.endswith(".gz"):
        compression = "gzip"
    df_query_feature_table = pd.read_csv(query_feature_table, compression=compression, sep="\t")
    if ("ASV" not in list(df_query_feature_table.columns)) or ("Count" not in list(df_query_feature_table.columns)) or ("SampleID" not in list(df_query_feature_table.columns)):
        df_query_feature_table = pd.read_csv(query_feature_table, compression=compression, sep="\t", index_col=0).stack().reset_index()
        df_query_feature_table.columns = ['ASV', 'SampleID', 'Count']
    else:
        df_query_feature_table = df_query_feature_table[['ASV', 'SampleID', 'Count']]

    #-----------------------------------------
    # Check if all required source files exist
    #-----------------------------------------

    # read or create source mapping file
    if which_source is not None:
        df_heterogenous_sources = pd.DataFrame(set(df_query_feature_table.SampleID), columns=['SampleID'])
        if os.path.exists("%s/%s"%(source_dir, which_source)):
            df_heterogenous_sources['SourceID'] = which_source
            df_heterogenous_sources['SourceDir'] = "%s/%s"%(source_dir, which_source)
        else:
            raise RuntimeError("could not find data folder for source %s."%(which_source))
    elif source_mapping_file is not None:
        df_heterogenous_sources = pd.read_csv(source_mapping_file, sep="\t")
        if ("SampleID" not in list(df_heterogenous_sources.columns) or "SourceID" not in list(df_heterogenous_sources.columns)) or "SourceDir" not in list(df_heterogenous_sources.columns):
            raise RuntimeError("could not find SampleID or SourceID or SourceDir in source mapping file.")
        else:
            df_heterogenous_sources = df_heterogenous_sources[['SampleID', 'SourceID', 'SourceDir']]

    dict_source_dir = {} # keys: source id, values: path to data files of each source
    df_unique_sources = df_heterogenous_sources[['SourceID', 'SourceDir']].drop_duplicates()
    for source, source_dir in zip(df_unique_sources.SourceID, df_unique_sources.SourceDir):
        # source folder exists?
        if not os.path.exists(source_dir):
            raise RuntimeError("could not find data folder for source %s."%(source))
        # sample metadata exists?
        if not os.path.exists("%s/sample_meta.txt"%(source_dir)):
            raise RuntimeError("could not find sample metadata for source %s."%(source))
        # feature table exists?
        if (not os.path.exists("%s/feature_table.txt"%(source_dir))) and (not os.path.exists("%s/feature_table.txt.gz"%(source_dir))):
            raise RuntimeError("could not find feature table for source %s."%(source))
        # asv sequences exist?
        if not os.path.exists("%s/asv_sequences.fasta"%(source_dir)):
            raise RuntimeError("could not find sequence file in fasta format for source %s."%(source))
        # blast database file exist? (check .nhr file as a marker)
        # if not, create blast database files
        if not os.path.exists("%s/%s.nhr"%(source_dir, source)):
            print("making blast db for source %s >>>>"%(source))
            #print("makeblastdb -in %s/%s/asv_sequences.fasta -title %s -dbtype nucl -out %s/%s/"%(source_dir, source, source, source_dir, source))
            ok = os.system("makeblastdb -in %s/asv_sequences.fasta -title %s -dbtype nucl -out %s/%s"%(source_dir, source, source_dir, source))
            if ok != 0:
                raise RuntimeError("makeblastdb fails against fasta file %s/asv_sequences.fasta."%(source_dir))

        dict_source_dir[source] = source_dir
    print("found %d unique sources: %s >>>>"%(len(dict_source_dir), (';').join(dict_source_dir.keys())))

    #--------------------------------------------------------------
    # quantify oral bacterial fraction for query microbiome samples
    #--------------------------------------------------------------
    if inference_method in ["FEAST", "RF"]:
        df_oralperc_output_summary = None
    else:
        raise RuntimeError("unknown inference method: %s."%(inference_method))

    for source, source_dir in dict_source_dir.items():
        print("working on samples matched to source %s >>>>"%(source))

        # read feature table of source microbiome
        df_source_feature_table = pd.read_csv("%s/feature_table.txt"%(source_dir), sep="\t")
        if ("ASV" not in list(df_source_feature_table.columns)) or ("Count" not in list(df_source_feature_table.columns)) or ("SampleID" not in list(df_source_feature_table.columns)):
            if os.path.exists("%s/feature_table.txt.gz"%(source_dir)):
                df_source_feature_table = pd.read_csv("%s/feature_table.txt.gz"%(source_dir), compression="gzip", sep="\t", index_col=0).stack().reset_index()
            elif os.path.exists("%s/feature_table.txt"%(source_dir)):
                df_source_feature_table = pd.read_csv("%s/feature_table.txt"%(source_dir), sep="\t", index_col=0).stack().reset_index()
            df_source_feature_table.columns = ['ASV', 'SampleID', 'Count']
        else:
            df_source_feature_table = df_source_feature_table[['ASV', 'SampleID', 'Count']]

        # read sample metadata of source microbiome
        df_source_sample_meta = pd.read_csv("%s/sample_meta.txt"%(source_dir), sep="\t")
        if ("SampleID" not in list(df_source_sample_meta.columns)) or ("Environment" not in list(df_source_sample_meta.columns)):
            raise RuntimeError("could not find SampleID or Environment in sample metadata file of source %s."%(source))

        # check for consistency of SampleIDs between feature table and sample metadata
        common_samples = set(df_source_sample_meta.SampleID).intersection(set(df_source_feature_table.SampleID))
        df_source_sample_meta = df_source_sample_meta[df_source_sample_meta.SampleID.isin(common_samples)]
        df_source_feature_table = df_source_feature_table[df_source_feature_table.SampleID.isin(common_samples)]

        # find query microbiome samples that use the current source
        query_samples_using_curr_source = list(set(df_heterogenous_sources[df_heterogenous_sources.SourceID==source].SampleID))
        df_query_feature_table_curr_source = deepcopy(df_query_feature_table[df_query_feature_table.SampleID.isin(query_samples_using_curr_source)])
        df_query_feature_table_curr_source = df_query_feature_table_curr_source[df_query_feature_table_curr_source.Count > 0] # remove query ASVs that are zero

        # read source asv sequences
        asv_sequences_curr_source = SeqIO.parse(open("%s/asv_sequences.fasta"%(source_dir)), 'fasta')
        asv_dict_curr_source = {}
        for fasta in asv_sequences_curr_source:
            asv_dict_curr_source[fasta.id] = str(fasta.seq)

        # write a temporary fasta file for current query ASVs
        asv_dict_curr_query = {}
        with open("%s/current_query_asv_sequences.fasta"%(output_dir), 'w') as out_file:
            for fasta in obj_query_asv_sequences:
                asv, seq = fasta.id, str(fasta.seq)
                if asv in list(df_query_feature_table_curr_source.ASV):
                    asv_dict_curr_query[asv] = seq
                    out_file.write(">%s\n%s\n"%(asv, seq))

        # blast against source: the blast output file is saved in the output dir
        blast_against_source(
            query_sequence_fasta = "%s/current_query_asv_sequences.fasta"%(output_dir),
            source = source,
            source_dir = dict_source_dir[source],
            evalue_cutoff = evalue_cutoff,
            num_of_threads = num_of_threads,
            output_dir = output_dir
        )

        # remove temporary query sequence fasta file
        os.system("rm %s/current_query_asv_sequences.fasta"%(output_dir))

        # parse and filter blast output
        df_blast_filtered = parse_blast(
            asv_dict_source = asv_dict_curr_source,
            asv_dict_query = asv_dict_curr_query,
            source = source,
            output_dir = output_dir
        )

        # match sequences between sources and samples
        df_joined_asv_count = merge_query_source_count_tables(
            df_blast_filtered = df_blast_filtered,
            df_query_asv_count = df_query_feature_table_curr_source,
            df_source_asv_count = df_source_feature_table,
            source = source,
            output_dir = output_dir
        )

        # predict oral percentage in query samples
        if inference_method == "FEAST":
            # generate FEAST input files

            # ASV table
            if not os.path.exists("%s/otu_table_%s.txt"%(output_dir, source)):
                df_joined_asv_relab = pd.pivot_table(df_joined_asv_count, index='SourceASV', columns='SampleID', values='Count', aggfunc=np.sum).fillna(0).astype(int)
                df_joined_asv_relab.to_csv("%s/otu_table_%s.txt"%(output_dir, source), sep="\t")

            # Sample metdata and R script
            query_samples = list(df_joined_asv_relab.columns)
            source_samples = list(df_joined_asv_relab.columns)
            output_files_to_be_produced = []
            for qsid in query_samples:
                output_file_path = "%s/%s_source_contributions_matrix.txt"%(output_dir, qsid)
                if os.path.exists(output_file_path):
                    continue

                # generate metadata
                metadata =  [[qsid, "Query", "Sink", 1]]
                for ssid, env in zip(df_source_sample_meta.SampleID, df_source_sample_meta.Environment):
                    if ssid in source_samples:
                        metadata.append([ssid, 'SOURCE__' + env, 'Source', 1])
                df_feast_metadata = pd.DataFrame(metadata, columns=['SampleID','Env','SourceSink','id'])
                df_feast_metadata.to_csv("%s/metadata_%s.txt"%(output_dir, qsid), sep="\t", index=False)

                # generate R scripts
                rscript = open("%s/run_feast_%s.r"%(output_dir, qsid),"w")
                rscript.write("suppressPackageStartupMessages(library(FEAST))\n"\
                              "metadata <- Load_metadata(metadata_path = \"%s/metadata_%s.txt\")\n"\
                              "otus <- Load_CountMatrix(CountMatrix_path = \"%s/otu_table_%s.txt\")\n"\
                              "FEAST_output <- FEAST(C = otus, metadata = metadata, COVERAGE = 1000,"\
                              "different_sources_flag = 0, dir_path = \"%s\",outfile=\"%s\")"%(
                                  output_dir, qsid, output_dir, source, output_dir, qsid)
                              )
                rscript.close()

                # Run R scripts
                os.system("Rscript %s/%s &"%(output_dir, "run_feast_%s.r"%(qsid)))
                output_files_to_be_produced.append("%s/%s_source_contributions_matrix.txt"%(output_dir, qsid))

                # Depending on the source sample size, run FEAST can be very computationally intensive.
                # Limit the number of parallel jobs to max_feast_runs
                # wait until at least one job finishes
                while len(output_files_to_be_produced) >= max_feast_runs:
                    for ofile in output_files_to_be_produced:
                        if os.path.exists(output_file_path):
                            output_files_to_be_produced.remove(output_file_path)

        elif inference_method == "RF":

            # generate random forest models per source ASV
            df_bodysite_rf_cla_source_asvs = build_per_asv_random_forest_model(
                df_source_sample_meta = df_source_sample_meta,
                df_source_asv_count = df_source_feature_table,
                source = source,
                source_dir = source_dir,
                num_of_threads = num_of_threads
            )

            # quantify oral bacterial fraction
            df_source_percentages = predict_oral_percentage(
                df_bodysite_rf_cla_source_asvs = df_bodysite_rf_cla_source_asvs,
                df_joined_asv_count = df_joined_asv_count,
                num_of_threads = num_of_threads,
                prevalence_cutoff = prevalence_cutoff
            )

            if df_oralperc_output_summary is None:
                df_oralperc_output_summary = deepcopy(df_source_percentages)
            else:
                df_oralperc_output_summary = pd.concat([df_oralperc_output_summary, df_source_percentages], axis=0)

    #-------------------------------------------------
    # Generate a summary table if using FEAST and save
    #-------------------------------------------------
    if inference_method == "FEAST":
        for sid in df_query_feature_table.columns:
            feast_ofile = "%s/%s_source_contributions_matrix.txt"%(output_dir, sid)
            if os.path.exists(feast_ofile):
                df_res = pd.read_csv(feast_ofile, sep="\t").T.reset_index()
                df_res.columns = ['Environment', 'Percentage']
                df_res['SampleID'] = sid
                if df_oralperc_output_summary is None:
                    df_oralperc_output_summary = deepcopy(df_res)
                else:
                    df_oralperc_output_summary = pd.concat([df_oralperc_output_summary, df_res], axis=0)
        df_oralperc_output_summary = df_oralperc_output_summary[df_oralperc_output_summary.Percentage > 0.0]
        df_oralperc_output_summary = pd.pivot_table(df_oralperc_output_summary, index='SampleID', columns='Environment', values='Percentage', aggfunc=np.sum).fillna(0)
        df_oralperc_output_summary.to_csv("%s/Feast_oral_percentage_summary.txt"%(output_dir), sep="\t")
    elif inference_method == "RF":
        df_oralperc_output_summary = pd.pivot_table(df_oralperc_output_summary, index="SampleID", columns="Environment", values="Percentage")
        df_oralperc_output_summary.to_csv("%s/Random_Forest_oral_percentage_summary.txt"%(source_dir), sep="\t")

