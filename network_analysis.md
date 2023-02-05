## Networking with CoNet/Cytoscape
---
### Splitting ASV tables

```
setwd("./Github/EEID_analysis/")
library(tidyverse)
#read ASV table
asv.tbl<-read.delim("otu_table.txt", header=T, row.names=1)

#read metadata
eeid.meta<-read.delim("Data/Metadata_NSFEEID_16SRuns_1234Merged_plusExpData.txt", header=T)

#Temperatures = 6, 14, 22
#doses = "Control" "5x10.3"  "5x10.4"  "5x10.5"  "5x10.6"  "C"

##subset each ASV table for each dose, ignore dose
#split map
meta_split<-split(eeid.meta, eeid.meta$Temperature)

#make temp-specifc tables, temperature=6
asv.temp6<-asv.tbl[,names(asv.tbl) %in% meta_split$T6$SampleID]
asv.temp6$OTU<-row.names(asv.tbl)
asv.temp6<-asv.temp6 %>% relocate(OTU)
write.table(asv.temp6, 'split_asv_tables/asv.tmp6.txt', quote=F, sep="\t", row.names = F)

#make temp-specifc tables, temperature=14
asv.temp14<-asv.tbl[,names(asv.tbl) %in% meta_split$T14$SampleID]
asv.temp14$OTU<-row.names(asv.tbl)
asv.temp14<-asv.temp14 %>% relocate(OTU)
write.table(asv.temp14, 'split_asv_tables/asv.tmp14.txt', quote=F, sep="\t", row.names = F)

#make temp-specifc tables, temperature=22
asv.temp22<-asv.tbl[,names(asv.tbl) %in% meta_split$T22$SampleID]
asv.temp22$OTU<-row.names(asv.tbl)
asv.temp22<-asv.temp22 %>% relocate(OTU)
write.table(asv.temp22, 'split_asv_tables/asv.tmp22.txt', quote=F, sep="\t", row.names = F)


#subset each ASV table for each dose, ignore temperature
#split map
meta_split<-split(eeid.meta, eeid.meta$Dose)

asv.dose106<-asv.tbl[,names(asv.tbl) %in% meta_split$d5x10.6$ SampleID]
asv.dose106$OTU<-row.names(asv.tbl)
asv.dose106<-asv.dose106 %>% relocate(OTU)
write.table(asv.dose106, 'split_asv_tables/asv.dose106.txt', quote=F, sep="\t", row.names = F)

asv.dose103<-asv.tbl[,names(asv.tbl) %in% meta_split$d5x10.3$SampleID]
asv.dose103$OTU<-row.names(asv.tbl)
asv.dose103<-asv.dose103 %>% relocate(OTU)
write.table(asv.dose103, 'split_asv_tables/asv.dose103.txt', quote=F, sep="\t", row.names = F)

asv.dose104<-asv.tbl[,names(asv.tbl) %in% meta_split$d5x10.4$SampleID]
asv.dose104$OTU<-row.names(asv.tbl)
asv.dose104<-asv.dose104 %>% relocate(OTU)
write.table(asv.dose104, 'split_asv_tables/asv.dose104.txt', quote=F, sep="\t", row.names = F)

asv.dose105<-asv.tbl[,names(asv.tbl) %in% meta_split$d5x10.5$SampleID]
asv.dose105$OTU<-row.names(asv.tbl)
asv.dose105<-asv.dose105 %>% relocate(OTU)
write.table(asv.dose105, 'split_asv_tables/asv.dose105.txt', quote=F, sep="\t", row.names = F)

asv.dosectrl<-asv.tbl[,names(asv.tbl) %in% meta_split$dControl$SampleID]
asv.dosectrl$OTU<-row.names(asv.tbl)
asv.dosectrl<-asv.dosectrl %>% relocate(OTU)
write.table(asv.dosectrl, 'split_asv_tables/asv.dosectrl.txt', quote=F, sep="\t", row.names = F)

#subset table for each dose/temp combo
```
Networks run in CoNet app in Cytoscape. Spearman correlations, p<0.05, Benjamini-Hochberg correct p-value, OTU abundance >2 

### Adding metadata to Cytoscape tables
```
test<-read.csv("Network_node_tables/dose103-t14.csv", header=T)
setwd("EEID_analysis")

#change to nodes folder
setwd("./Network_node_tables/")

#list files
setwd("/Users/patty/OneDrive/Documents/Github/EEID_analysis/Network_node_tables/")
nodes<-list.files("./")

#read in all files
node_tables<-sapply(nodes, function(x) read.csv(x, header=T))

#read in tax/inhibiton file
tax_inib<-read.delim("/Users/patty/OneDrive/Documents/Github/EEID_analysis/taxa_calls_inhib.txt", header=T)

#split taxonomy 
library(tidyverse)
tax_inhib<-separate(tax_inib, Taxon, sep = ";", into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

#add taxonomy to all files
node_tables<-lapply(node_tables, function(x) merge(x, tax_inhib, by.x='Label', by.y='Feature.ID', all.x=T, all.y=F))

#write them back to a file
dir.create("nodes_with_tax_inhib")
setwd("nodes_with_tax_inhib/")
for(i in names(node_tables)){
  write.table(node_tables[[i]], paste0(i,".txt"), sep="\t", row.names=F, quote=F)
}
setwd("/Users/patty/OneDrive/Documents/Github/EEID_analysis/")

#get node stats
library(plyr)
library(ggplot2)
library(dplyr)

#set working directory to node folder
setwd("Network_node_tables/")

#list files of interest
bsg <- dir(pattern = "\\.csv$")

#add file names to list of files
names(bsg) <- basename(bsg)

#read in the data and catenate
dbsg <- ldply(bsg, read.csv)

#read in sample names
samps<-read.delim("/Users/patty/OneDrive/Documents/Github/EEID_analysis/network_samps.txt", header=T)

#strip '.csv
samps$Name<-gsub("\\.csv", "", samps$Name)
dbsg$.id<-gsub("\\.csv", "", dbsg$.id)
samps$Name<-gsub("\\-", "", samps$Name)
dbsg$.id<-gsub("\\-", "", dbsg$.id)
dbsg$Name<-dbsg$.id

samps$Name<-dbsg$Name

samps$Name<-as.factor(samps$Name)
dbsg$Name<-as.factor(dbsg$Name)

dbsg$Row <- seq.int(nrow(dbsg))
dbsg$Name2 <-paste(dbsg$Row, dbsg$Name)

samps$Row <- seq.int(nrow(samps))
samps$Name2 <-paste(samps$Row, samps$Name)


#add sample metadata (dose & temp) to node table
dbsg2<-merge(dbsg, samps, by="Name2")

ggplot(dbsg2[,-which(dbsg2$Dose == N/A)], aes(Dose, posdegree))+
  geom_boxplot()+
  facet_wrap(~Temperature)

ggplot(dbsg2, aes(Dose, degree))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~Temperature)

```
