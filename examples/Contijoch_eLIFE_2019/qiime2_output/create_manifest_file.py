import pandas as pd
import os

df_meta = pd.read_csv("sample-metadata.tsv", sep="\t", skiprows=[1])
manifest = []
for acc in df_meta['sample-id']:
    manifest.append([acc, "$PWD/%s.fastq"%(acc)])
df_manifest = pd.DataFrame(manifest, columns=['sample-id','absolute-filepath'])
df_manifest.to_csv("manifest.tsv", sep="\t", index=False)

