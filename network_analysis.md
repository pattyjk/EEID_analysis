## Networking with CoNet/Cytoscape

R

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

```
#using CoNet/Cytoscape
#calculate full network
java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input C:\Users\patty\OneDrive\Documents\Github\EEID_analysis\otu_table.txt --measure2 supp --scoremergestrategy mean --correlnonrandp none --multicorr none --minsupport 2 --randroutine none --minetmiestimator mi.shrink --kernelwidth 0.25 --nantreatmentparam 1 --networkmergestrategy
```
