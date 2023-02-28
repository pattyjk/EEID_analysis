## Calculation of 16S rRNA copy number based on PICRUSt normalization
```
#load package
library(ggplot2)

#read in OTU tables

no_norm_tab<-read.delim("picrust_out_notnorm/KO_metagenome_out/pred_metagenome_unstrat.tsv", row.names = 1, header=T)
norm_tab<-read.delim("picrust_out_norm/KO_metagenome_out/pred_metagenome_unstrat.tsv", row.names=1, header=T)

#calculate the rRNA copy number per ASV
div_tab<-no_norm_tab/norm_tab

#get averages per column, e.g. per sample
asv_means<-as.data.frame(colMeans(na.omit(div_tab)))
asv_means$SampleID<-row.names(asv_means)
names(asv_means)<-c("Weighted_operon", "SampleID")


#read in metadata
meta<-read.delim("Data/Metadata_NSFEEID_16SRuns_1234Merged_plusExpData.txt", header=T)

#append 16S operon number to metadata
meta<-merge(meta, asv_means, by="SampleID")

ggplot(meta, aes(Dose2, Weighted_operon))+
  geom_boxplot()+
  facet_wrap(~Temperature2)

ggplot(meta, aes(as.numeric(Temperature2), Weighted_operon))+
  geom_point()+
  facet_wrap(~Dose)

ggplot(meta, aes(TimeWeek, Weighted_operon))+
  geom_point()+
  facet_wrap(~Dose)

```
