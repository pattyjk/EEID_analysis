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
picrust2_pipeline.py -s rep_set/dna-sequences.fasta -i otu_table/feature-table.biom -o picrust_out_norm #this is normalized
picrust2_pipeline.py -s rep_set/dna-sequences.fasta -i otu_table/feature-table.biom -o picrust_out_notnorm --skip_norm
```

## Bring data into R
```
setwd("Github/EEID_analysis")

#read in KO IDs
full_kegg<-read.delim("./full_kegg.txt", header = T)

#read metadata
meta<-read.delim("Data/Metadata_NSFEEID_16SRuns_1234Merged_plusExpData.txt", header=T)

#read in KO table
ko_table<-read.delim("picrust_out_norm/KO_metagenome_out/pred_metagenome_unstrat.tsv", header=T, row.names=1)

#read in original ASV table
asv.tbl<-read.delim("./otu_table.txt", header=T, row.names=1)

#load vegan
library(vegan)

#calculate Bray-Curtis distance matrix for each table (KO/ASV)
asv.dis<-vegdist(t(asv.tbl), method='bray')
ko.dis<-vegdist(t(ko_table), method='bray')

#calculate a mantel test to see if the patterns differ between the methods
mantel(asv.dis, ko.dis, method="pearson", permutations=999, strata = NULL,
       na.rm = FALSE, parallel = getOption("mc.cores"))
```
### Mantel test output
Call:
  mantel(xdis = asv.dis, ydis = ko.dis, method = "pearson", permutations = 999,      strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores")) 

Mantel statistic r: 0.483 
Significance: 0.001 

Upper quantiles of permutations (null model):
  90%    95%  97.5%    99% 
  0.0145 0.0183 0.0215 0.0238 
Permutation: free
Number of permutations: 999

They matrices are mathmatically different with only ~48% correlation between the two. Plot the KO data as a PCoA

```
#calculate PCoA with BRay-Curtis
ko_pcoa<-capscale(t(ko_table)  ~ 1, distance='bray')


#pull out x/y coordinates
ko.scores<-scores(ko_pcoa)

#grab only sample coordinates, write to data frame
ko.coords<-as.data.frame(ko.scores$sites)


#create sample names as a column
ko.coords$SampleID<-row.names(ko.coords)

#map back meta data
ko.coords<-merge(ko.coords, meta, by.x='SampleID', by.y='SampleID')

#calculate percent variation explained for first two axis
100*round(ko_pcoa$CA$eig[1]/sum(ko_pcoa$CA$eig), 3)
#48.6
100*round(ko_pcoa$CA$eig[2]/sum(ko_pcoa$CA$eig), 3)
#24.8

#plot PCoA
ggplot(ko.coords, aes(MDS1, MDS2, colour=Dose2))+
  geom_point(aes(size=2))+
  theme_bw()+
  xlab("PC1- 48.6%")+
  ylab("PC2- 24.8%")

#what's causing the arch?
par(mfrow=c(2,2))
plot(ko.coords$Temperature2, ko.coords$MDS1)
plot(ko.coords$Temperature2, ko.coords$MDS2)
plot(ko.coords$Dose2, ko.coords$MDS1)
plot(ko.coords$Dose2, ko.coords$MDS2)

cor.test(as.numeric(ko.coords$Dose2), ko.coords$MDS2)
cor.test(as.numeric(ko.coords$Dose2), ko.coords$MDS1)
cor.test(as.numeric(ko.coords$Temperature2), ko.coords$MDS1)
cor.test(as.numeric(ko.coords$Temperature2), ko.coords$MDS1)

#temperature is significant for MDS1/MDS2

#calculate PCoA with BRay-Curtis
ko_pcoa<-decorana(t(ko_table))

#pull out x/y coordinates
ko.scores<-scores(ko_pcoa)

#grab only sample coordinates, write to data frame
ko.coords<-as.data.frame(ko.scores)


#create sample names as a column
ko.coords$SampleID<-row.names(ko.coords)

#map back meta data
ko.coords<-merge(ko.coords, meta, by.x='SampleID', by.y='SampleID')

#plot PCoA
ggplot(ko.coords, aes(DCA1, DCA2, colour=Dose2))+
  geom_point(aes(size=2))+
  theme_bw()+
  xlab("DCA1")+
  ylab("DCA2")

##analyze KO's
library(reshape2)

#reshape data
ko_table2<-ko_table
ko_table2$KO<-row.names(ko_table2)
ko_m<-melt(ko_table2)
str(ko_m)
names(ko_m)<-c("KO", "SampleID", "Abun")

#add metadata
ko_m<-merge(ko_m, meta, by='SampleID')

ggplot(ko_m, aes(x = SampleID, y = KO, fill = Abun)) +
  geom_tile()
```

## Look at secondary metabolites
````
#add kegg to table
ko_m<-merge(ko_m, full_kegg, by.x='KO', by.y='KO')

#filter out non-secondary metabolite KO
second_metab<-ko_m[which(ko_m$Level2 == "Biosynthesis of other secondary metabolites"),]

#check if filtering looks good, cause I paranoid like that
dim(second_metab)

#plot all secondary metabolites
ggplot(second_metab, aes(Dose2, Abun))+
  geom_boxplot()+
  facet_wrap(~Temperature)+
  scale_y_log10() 
#dose increase secondary metabolite procution, except at 6C

#split  table by Temperature
meta_split<-split(second_metab,f = second_metab$Temperature)

##calculate correlation between dose/abundance
set.seed(515)
cor.test(meta_split$T14$Abun, as.numeric(meta_split$T14$Dose2))
#no correlation

cor.test(meta_split$T6$Abun, as.numeric(meta_split$T6$Dose2))
#signficant negative correlation

cor.test(meta_split$T22$Abun, as.numeric(meta_split$T22$Dose2))
#significatn positive correlation

cor.test(second_metab$Abun, as.numeric(second_metab$Temperature2))
#signficatn negative correlation with temperature
```
