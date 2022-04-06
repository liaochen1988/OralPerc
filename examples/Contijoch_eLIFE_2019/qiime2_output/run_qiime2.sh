#---------------------------
# import pairend fastq files
#---------------------------
#qiime tools import \
#  --type 'SampleData[SequencesWithQuality]' \
#  --input-path manifest.tsv \
#  --output-path single-end-demux.qza \
#  --input-format SingleEndFastqManifestPhred33V2

#qiime demux summarize \
#  --i-data single-end-demux.qza \
#  --o-visualization single-end-demux.qzv

#qiime tools view single-end-demux.qzv

#--------------
# DADA2 denoise
#--------------
#qiime dada2 denoise-single \
#  --i-demultiplexed-seqs single-end-demux.qza \
#  --p-trunc-len 250 \
#  --p-n-threads 36 \
#  --o-table table.qza \
#  --o-representative-sequences rep-seqs.qza \
#  --o-denoising-stats denoising-stats.qza

#qiime tools export \
#  --input-path table.qza \
#  --output-path feature-table
#biom convert -i feature-table/feature-table.biom -o feature-table/feature-table.from_biom.txt --to-tsv

#qiime tools export \
#   --input-path rep-seqs.qza \
#   --output-path asv-sequences

#qiime feature-table summarize \
#  --i-table table.qza \
#  --o-visualization table.qzv \
#  --m-sample-metadata-file sample-metadata.tsv
#qiime tools view table.qzv

#qiime feature-table tabulate-seqs \
#  --i-data rep-seqs.qza \
#  --o-visualization rep-seqs.qzv
#qiime tools view rep-seqs.qzv

#qiime metadata tabulate \
#  --m-input-file denoising-stats.qza \
#  --o-visualization denoising-stats.qzv
#qiime tools view denoising-stats.qzv

#----------------------
# phylogenetic analysis
#----------------------
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

qiime tools export \
  --input-path unrooted-tree.qza \
  --output-path exported-tree

#---------------------
# rarefaction analysis
#---------------------
#qiime diversity alpha-rarefaction \
#  --i-table table.qza \
#  --i-phylogeny rooted-tree.qza \
#  --p-max-depth 4000 \
#  --m-metadata-file sample-metadata.tsv \
#  --o-visualization alpha-rarefaction.qzv
#qiime tools view alpha-rarefaction.qzv

#---------------------
# taxonomic assignment
#---------------------
qiime feature-classifier classify-sklearn \
  --i-classifier ../../../../silva-138-99-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --p-n-jobs 32 \
  --p-confidence 0.8 \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
