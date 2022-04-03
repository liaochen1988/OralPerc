import pandas as pd
import numpy as np
import sys, getopt, os
from Bio.Seq import Seq
from collections import Counter

def blast_against_source(
    query_sequence_fasta:str,   # path of fasta files of query sequences
    source:str,                 # name of source
    source_path:str,            # path of source database that query sequences are blasted to
    evalue_cutoff:float,        # evalue cutoff
    num_of_threads:int,         # number of threads for blast
    output_dir:str              # directory of output folder
)->None:
    blast_output_columns = "qseqid qseq qlen sseqid sseq slen sstrand pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    blast_command = "blastn -db %s -query %s -out %s/unfiltered_blast_against_source_%s.txt -outfmt \"7 %s\" -perc_identity 100 -evalue %2.2e -num_of_threads %d"(
        source, query_sequence_fasta, output_dir, source, blast_output_columns, evalue_cutoff, num_of_threads)
    ok = os.system(blast_command)
    if ok!=0: # program fails
        raise RuntimeError("command <%s> fails."%(blast_command))
    return None

def parse_blast(
    asv_dict_curr_source: dict,   # keys are HMP ASV IDs, values are corresponding sequences
    asv_dict_curr_query: dict,  # keys are user-provided ASV IDs, values are corresponding sequences
    output_dir:str,     # output dir directory
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
            raise RuntimeError("aligned sequences between %s and %s are not equal."%(qseqid, sseqid))

        # get full length sequences of HMP ASVs and user provided ASVs
        if sseqid in asv_dict_curr_source.keys():
            curr_source_asv = asv_dict_curr_source[sseqid]
            if sstrand == "minus":
                curr_source_asv = str(Seq(curr_source_asv).reverse_complement())
        else:
            raise RuntimeError("could not find %s in HMP asvs dictionary."%(sseqid))

        if qseqid in asv_dict_curr_query.keys():
            curr_user_asv = asv_dict_curr_query[qseqid]
        else:
            raise RuntimeError("could not find %s in user defined asvs dictionary."%(qseqid))

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
    if verbose:
        print("%d unique user provided sequences are matched to HMP sequences."%(len(set(df_blast_filtered.qseqid))))

    # save data to file
    df_blast_filtered.to_csv("filtered_blast_output.txt", sep="\t", index=False)

    return df_blast_filtered

def merge_query_source_count_tables(
    df_blast_filtered:pd.DataFrame,     # filtered blast output
    df_query_asv_count:pd.DataFrame,    # count table of query ASVs
    df_source_asv_count:pd.DataFrame,   # count table of source ASVs
    SourceID:str                        # Source ID
)->pd.DataFrame:

    # Important Notes:
    # Both query and source count tables are passed in long format (not wide format).

    # rename query sample ID to avoid overlaps with source sample IDs
    df_query_asv_count = df_query_asv_count[['SampleID','QueryASV','Count']]
    df_query_asv_count.SampleID = ['QUERY__'+sid for sid in df_query_asv_count.SampleID]

    # select and rename filtered blast output
    df2_blast_filtered = deepcopy(df_blast_filtered[['qseqid','sseqid']])
    df2_blast_filtered.columns = ['QueryASV','source_asv']

    # append a column of matched source ASVs in query count table
    df_query_asv_count = pd.merge(df_query_asv_count, df2_blast_filtered, left_on='QueryASV', right_on='QueryASV', how='left')

    # rename query ASVs that are not matched to any source ASV sequence by adding prefix "Unmapped__"
    df_query_asv_count.loc[df_query_asv_count.source_asv.isnull(), 'source_asv'] = ['Unmapped__'+asv for asv in df_query_asv_count.loc[df_query_asv_count.source_asv.isnull(), 'QueryASV']]

    # concatenate the two tables
    df_source_asv_count = df_source_asv_count.rename({'ASV':'source_asv'}, axis=1)
    df_source_asv_count.SampleID = ['SOURCE__'+sid for sid in df_source_asv_count.SampleID]
    df_joined_asv_count = pd.concat([df_query_asv_count[['SampleID','source_asv','Count']], df_source_asv_count[['SampleID','source_asv','Count']]], axis=0)
    df_joined_asv_count = df_joined_asv_count[df_joined_asv_count.Count > 0] # remove ASVs of zero count
    df_joined_asv_count.to_csv("joined_asv_count_for_%s.csv"%(SourceID), sep="\t", index=False)

    return df_joined_asv_count

def build_per_asv_random_forest_model(
    df_source_sample_meta:pd.DataFrame,     # meta data of source
    df_source_asv_count:pd.DataFrame,       # relative abundance of source ASVs (samples by ASVs)
    output_dir:str,                         # path of output dir
    num_of_threads:int,                     # number of treads
    source:str                              # Source ID
)->pd.DataFrame:

    # do not recreate the wheel: check if output file has already been generated
    if os.path.exists("%s/rf_bodysite_prob_%s.txt"%(output_dir, source)):
        df_bodysite_rf_cla_source_asvs = pd.read_csv("%s/rf_bodysite_prob_%s.txt"%(output_dir, source), sep="\t")
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
        df_probs_of_body_site = pd.DataFrame(probs_of_body_site, columns=['SourceAsv','BodySite','Prob'])
        return df_probs_of_body_site

    # run in parallel for all source ASVs
    bodysite_probs = Parallel(n_jobs=num_of_threads)(delayed(build_rf)(asv) for asv in df_source_asv_relab.index)
    df_bodysite_rf_cla_source_asvs = pd.concat(bodysite_probs, ignore_index=True)
    df_bodysite_rf_cla_source_asvs.to_csv("%s/rf_bodysite_prob_%s.txt"%(output_dir, source), sep="\t", index=False)

    return df_bodysite_rf_cla_source_asvs

def quantify_oral_bacteria_fraction(
    df_bodysite_rf_cla_source_asvs:pd.DataFrame,    # body site classification for source ASVs
    df_joined_asv_count:pd.DataFrame,               # joined (both query and source) asv count table
    num_of_threads:int,                             # number of treads
    prevalence_cutoff:float                         # minimum prevalence to classify body site origin from sources
):
    # convert count to relative abundance
    df_joined_asv_relab = deepcopy(pd.pivot_table(df_joined_asv_count, index="ASV", columns="SampleID", values="Count", aggfunc=np.sum).fillna(0))
    df_joined_asv_relab = df_joined_asv_relab.div(df_joined_asv_relab.sum(axis=0), axis=1)
    df_query_asv_relab = df_joined_asv_relab[[sid for sid in df_joined_asv_relab.columns if sid.startswith("QUERY__")]].T      # samples by ASVs
    df_source_asv_relab = df_joined_asv_relab[[sid for sid in df_joined_asv_relab.columns if sid.startswith("SOURCE__")]].T    # samples by ASVs

    # get all environment categories
    all_unique_envs = list(set(df_bodysite_rf_cla_source_asvs.BodySite))

    def identify_source_proportions(sid):
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

    source_proportions = Parallel(n_jobs=num_of_threads)(delayed(identify_source_proportions)(sid) for sid in df_query_asv_relab.index)
    df_source_proportions = pd.DataFrame(source_proportion, columns=['SampleID']+all_unique_envs+['Unknown']).set_index('SampleID').stack().reset_index()
    df_source_proportions.columns = ['SampleID', 'Environment','Proportion']

    return df_source_proportions

#------------------------
# Main program
#------------------------
if __name__ == "__main__":

    #----------------------------
    # Get command line parameters
    #----------------------------

    # parse command line arguments
    # --help
    # --feature_table <required>: a txt table with ASVs on the rows and samples on the columns
    # --asv_sequences <required>: a fasta file of ASV sequences
    # --HMP_source [default HMP_16S_v35]: choose between HMP_16S_v13, HMP_16S_v35 and HMP_16S_v69 (this parameter is used only when sample_meta does not specify the sources explicitly)
    # --customized_source_dir [default None]: directory of potential microbial sources
    #                                       : each source directory should contain 3 files: sample_meta, feature_table, and asv_sequences
    # --inference_method [default RF]: choose between FEAST and RF (Random Forest)
    # --max_feast_runs [default 5]: maximum number of samples run by feast simutaneously (used only when method == feast)
    # --evalue_cutoff [default 1e-10]: evalue cutoff
    # --num_of_threads [default 1]: number of threads
    # --output_dir [default './output']: directory of all output files

    # default parameters
    HMP_source = "HMP_16S_v35"
    customized_source = None
    customized_source_dir = None
    inference_method = "RF"
    max_feast_runs = 5
    evalue_cutoff = 1e-10
    num_of_threads = 5
    prevalence_cutoff = 0.052
    output_dir = "output"

    # parse command line
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hf:a:pcdimbnvo",["help","feature_table=", "asv_sequences=", "HMP_source",
                                                                 "customized_source", "customized_source_dir", "inference_method",
                                                                 "max_feast_runs", "evalue_cutoff", "num_of_threads", "prevalence_cutoff", "output_dir"])
    except getopt.GetoptError:
        print('oral_bacteria_detector.py -h -f <feature_table> -a <asv_sequences> -p [HMP_source] -c [customized_source] -d [customized_source_dir]
               -i [inference_method] -m [max_feast_runs] -b [evalue_cutoff] -n [num_of_threads] -v [prevalence_cutoff] -o [output_dir]')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('oral_bacteria_detector.py -h -f <feature_table> -a <asv_sequences> -p [HMP_source] -c [customized_source] -d [customized_source_dir]
                   -i [inference_method] -m [max_feast_runs] -b [evalue_cutoff] -n [num_of_threads] -v [prevalence_cutoff] -o [output_dir]')
            sys.exit(0)
        elif opt in ("-f", "--feature_table"):
            feature_table = arg
        elif opt in ("-a", "--asv_sequences"):
            asv_sequences = arg
        elif opt in ("-p", "--HMP_source"):
            HMP_source = arg
        elif opt in ("-c", "--customized_source"):
            customized_source = arg
        elif opt in ("-d", "customized_source_dir"):
            customized_source_dir = arg
        elif opt in ("-i", "--inference_method"):
            inference_method = arg
        elif opt in ("-m", "--max_feast_runs"):
            max_feast_runs = arg
        elif opt in ("-b", "--evalue_cutoff"):
            evalue_cutoff = arg
        elif opt in ("-n", "--num_of_threads"):
            num_of_threads = arg
        elif opt in ("-v", "--prevalence_cutoff"):
            prevalence_cutoff = arg
        elif opt in ("-o", "--output_dir"):
            output_dir = arg

    #----------------------
    # Load query microbiome
    #----------------------

    # read ASV sequences
    query_asv_sequences = SeqIO.parse(open(asv_sequences), 'fasta')

    # read feature table
    df_query_feature_table = pd.read_csv(feature_table, sep="\t")
    if list(df_query_feature_table.columns).sort() != ['ASV','Count','SampleID']:
        # convert to long format
        df_query_feature_table = df_query_feature_table.set_index('ASV').stack().reset_index()
        df_query_feature_table.columns = ['ASV', 'SampleID', 'Count']

    #-----------------------------------------
    # Check if all required source files exist
    #-----------------------------------------
    dict_source_path = {}
    if customized_source is None:
        dict_source_path[HMP_source] = "../HMP_sources/%s"%(HMP_source)
    else:
        if customized_source_dir is None:
            raise RuntimeError("customized_source_dir is not given.")
        df_customized_sources = pd.read_csv(customized_source, sep="\t")
        if list(df_customized_sources.columns).sort() != ['SampleID','SourceID']:
            raise RuntimeError("only two columns are allowed in customized source file: SampleID and SourceID.")
        else:
            unique_sources = set(df_customized_sources.SourceID)
            for source in unique_sources:
                if not source.startswith("HMP_"):
                    # User provided source

                    # source folder exists?
                    if os.path.exists("%s/%s"%(customized_source_dir, source)):
                        raise RuntimeError("could not find data folder for source %s."%(source))
                    # sample metadata exists?
                    if not os.path.exists("%s/%s/sample_meta.txt"%(customized_source_dir, source)):
                        raise RuntimeError("could not find sample metadata for source %s."%(source))
                    # feature table exists?
                    if not os.path.exists("%s/%s/feature_table.txt"%(customized_source_dir, source)):
                        raise RuntimeError("could not find feature table for source %s."%(source))
                    # asv sequences exist?
                    if not os.path_exists("%s/%s/asv_sequences.fasta"%(customized_source_dir, source)):
                        raise RuntimeError("could not find sequence file in fasta format for source %s."%(source))
                    # blast database file exist? (check .nhr file)
                    # if not, create blast database files
                    if not os.path_exists("%s/%s/%s.nhr"%(customized_source_dir, source, source)):
                        ok = os.system("makeblastdb -in %s/%s/asv_sequences.fasta -title %s -dbtype nucl -out %s/%s/"%(customized_source_dir, source, source, customized_source_dir, source))
                        if ok != 0:
                            raise RuntimeError("makeblastdb fails against fasta file %s/%s/asv_sequences.fasta."%(customized_source_dir, source))
                        else:
                            dict_source_path[source] = "%s/%s/"%(customized_source_dir, source)
                else:
                    # HMP source
                    dict_source_path[source] = "../HMP_sources/%s"%(source)

    #--------------------------------------------------------------
    # quantify oral bacterial fraction for query microbiome samples
    #--------------------------------------------------------------
    if inference_method == "RF":
        df_rf_summary = None
    for source, source_path in dict_source_path.items():

        # read feature table of source microbiome
        df_source_feature_table = pd.read_csv("%s/%s/feature_table.txt"%(customized_source_dir, source), sep="\t")
        if list(df_source_feature_table.columns).sort() != ['ASV','Count','SampleID']:
           # convert to long format
           df_source_feature_table = df_source_feature_table.set_index('ASV').stack().reset_index()
           df_source_feature_table.columns = ['ASV', 'SampleID', 'Count']

        # read sample metadata of source microbiome
        df_source_sample_meta = pd.read_csv("%s/%s/sample_meta.txt"%(customized_source_dir, source), sep="\t")
        if "SampleID" not in list(df_source_sample_meta.columns):
            raise RuntimeError("SampleID is not a column name in sample metadata file of %s."%(source))
        if "Environment" not in list(df_source_sample_meta.columns):
            raise RuntimeError("Environment is not a column name in sample metadata file of %s."%(source))

        # check for consistency of SampleIDs between feature table and sample metadata
        common_samples = set(df_source_sample_meta.SampleID).intersection(set(df_source_feature_table.SampleID))
        df_source_sample_meta = df_source_sample_meta[df_source_sample_meta.SampleID.isin(common_samples)]
        df_source_feature_table = df_source_feature_table[df_source_feature_table.SampleID.isin(common_samples)]

        # find query microbiome samples that use the current source
        query_samples_using_curr_source = list(set(df_customized_sources[df_customized_sources.SourceID==source].SampleID))
        df_query_feature_table_curr_source = deepcopy(df_query_feature_table[df_query_feature_table.SampleID.isin(query_samples_using_curr_source)])
        df_query_feature_table_curr_source = df_query_feature_table_curr_source.loc[~(df_query_feature_table_curr_source==0).all(axis=1)] # remove source ASVs that are zero for all samples

        # read source asv sequences
        asv_sequences_curr_source = SeqIO.parse(open("%s/%s/asv_sequences.fasta"%(customized_source_dir, source)), 'fasta')
        asv_dict_curr_source = {}
        for fasta in asv_sequences_curr_source:
            asv_dict_curr_source[fasta.id] = str(fasta.seq)

        # write a temporary fasta file for current query ASVs
        asv_dict_curr_query = {}
        with open("%s/current_query_asv_sequences.fasta"%(output_dir)) as out_file:
            for fasta in query_asv_sequences:
                asv, seq = fasta.id, str(fasta.seq)
                if asv in list(df_query_feature_table_curr_source.index):
                    asv_dict_curr_query[asv] = seq
                    out_file.write(">%s\n%s\n"%(asv, seq))

        # blast against source: the blast output file is saved in the output dir
        blast_against_source(
            query_sequence_fasta = "%s/current_query_asv_sequences.fasta"%(output_dir),
            source = source,
            source_path = dict_source_path[source],
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
            output_dir = output_dir
        )

        # match sequences between sources and samples
        df_joined_asv_count = merge_query_source_count_tables(
            df_blast_filtered = df_blast_filtered,
            df_query_asv_count = df_query_feature_table_curr_source,
            df_source_asv_count = df_source_feature_table,
            SourceID = source
        )

        # quantify oral fraction
        if inference_method == "FEAST":
            # generate FEAST input files

            # OTU table
            df_joined_asv_relab = pd.pivot_table(df_joined_asv_count, index='SampleID', columns='source_asv', values='Count', aggfunc=np.sum).fillna(0).astype(int)
            df_joined_asv_relab.to_csv("otu_table_%s.txt"%(source), sep="\t")

            # Sample metdata and R script
            query_samples = [sample_id for sample_id in df_joined_asv_relab.index if sample_id.startswith('QUERY__')]
            source_samples = [sample_id for sample_id in df_joined_asv_relab.index if sample_id.startswith('SOURCE__')]
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
                df_feast_metadata.to_csv("metadata_%s.txt"%(qsid), sep="\t", index=False)

                # generate R scripts
                rscript = open("run_feast_%s.r"%(qsid),"w")
                rscript.write("library(FEAST)\n"\
                              "metadata <- Load_metadata(metadata_path = \"metadata_%s.txt\")\n"\
                              "otus <- Load_CountMatrix(CountMatrix_path = \"otu_table.txt\")\n"\
                              "FEAST_output <- FEAST(C = otus, metadata = metadata, COVERAGE = 1000, different_sources_flag = 0, dir_path = \"%s\", outfile=\"%s\")"%(qsid, output_dir, qsid))
                rscript.close()

                # Run R scripts
                os.system("Rscript %s &"%("run_feast_%s.r"%(qsid)))
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
                output_dir = output_dir,
                num_of_threads = num_of_threads,
                source = source
            )

            # quantify oral bacterial fraction
            df_source_proportions = quantify_oral_bacteria_fraction(
                df_bodysite_rf_cla_source_asvs = df_bodysite_rf_cla_source_asvs,
                df_joined_asv_count = df_joined_asv_count,
                num_of_threads = num_of_threads,
                prevalence_cutoff = prevalence_cutoff
            )

            if df_rf_summary is None:
                df_rf_summary = deepcopy(df_source_proportions)
            else:
                df_rf_summary = pd.concat([df_rf_summary, df_source_proportions], axis=0)

        else:
            raise RuntimeError("unknown inference method: %s."%(inference_method))

    #-------------------------------------------------
    # Generate a summary table if using FEAST and save
    #-------------------------------------------------
    if inference_method == "FEAST":
        df_feast_summary = None
        for sid in df_query_feature_table.columns:
            feast_ofile = "%s/%s_source_contributions_matrix.txt"%(output_dir, sid)
            if os.path.exists(feast_ofile):
                df_res = pd.read_csv(feast_ofile, sep="\t").T.reset_index()
                df_res.columns = ['Environment', 'Fraction']
                df_res['SampleID'] = sid
                if df_feast_summary is None:
                    df_feast_summary = deepcopy(df_res)
                else:
                    df_feast_summary = pd.concat([df_feast_summary, df_res], axis=0)
        df_feast_summary = df_feast_summary[df_feast_summary.Fraction > 0.0]
        df_feast_summary = pd.pivot_table(df_feast_summary, index='SampleID', columns='Environment', values='Fraction', aggfunc=np.sum).fillna(0)
        df_feast_summary.to_csv("Feast_oral_fraction_summary.txt", sep="\t")
    elif inference_method == "RF":
        df_rf_summary = pd.pivot_table(df_rf_summary, index="SampleID", columns="Environment", values="Proportion")
        df_rf_summary.to_csv("Random_Forest_oral_fraction_summary.txt", sep="\t")

