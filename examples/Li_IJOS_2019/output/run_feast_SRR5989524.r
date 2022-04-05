suppressPackageStartupMessages(library(FEAST))
metadata <- Load_metadata(metadata_path = "/Users/liaoc/Projects/OralPerc/examples/Li_IJOS_2019/output/metadata_SRR5989524.txt")
otus <- Load_CountMatrix(CountMatrix_path = "/Users/liaoc/Projects/OralPerc/examples/Li_IJOS_2019/output/otu_table_saliva.txt")
FEAST_output <- FEAST(C = otus, metadata = metadata, COVERAGE = 1000,different_sources_flag = 0, dir_path = "/Users/liaoc/Projects/OralPerc/examples/Li_IJOS_2019/output",outfile="SRR5989524")