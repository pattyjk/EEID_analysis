## PICRUSt2

```
#convert qiime2 artifacts to non-qiime artifacts
#load qiime2
source activate qiime2-2022.2

#export taxonomy
qiime tools export --input-path EEID_16S_1234_Taxonomywplasmid.qza --output-path taxonomy

#export rep set
qiime tools export --input-path Repset_merged_noPl_EEID_AqAdult_TempDose_Exp_1234.qza --output-path rep_set

#export otu table
qiime tools export --input-path Filtlowsamp_CopCor_True_abund_estimate_contamfilt_samplefilt_OTU-table_wplasmid_EEID_AqAdult_TempDose_Exp_merged1234.qza --output-path otu_table
biom convert -i otu_table/feature-table.biom -o otu_table.txt --to-tsv

#load the thingy
source activate picrust2

#run pipeline with and without 16S copy number normalization
picrust2_pipeline.py -s rep_set/dna-sequences.fasta -i otu_table/feature-table.biom -o picrust_out_norm
picrust2_pipeline.py -s rep_set/dna-sequences.fasta -i otu_table/feature-table.biom -o picrust_out_notnorm --skip_norm

#for both
63 of 3140 ASVs were above the max NSTI cut-off of 2.0 and were removed from the downstream analyses.
```
