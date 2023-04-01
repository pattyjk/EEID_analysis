## Networking with CoNet/Cytoscape

### Splitting ASV tables

```
library(tidyverse)
#read ASV table
asv.tbl<-read.delim("otu_table.txt", header=T, row.names=1)

#read metadata
eeid.meta<-read.delim("Data/Metadata_NSFEEID_16SRuns_1234Merged_plusExpData.txt", header=T)

#Temperatures = 6, 14, 22
#doses = "Control" "5x10.3"  "5x10.4"  "5x10.5"  "5x10.6"  "C"
#TimeWeekCat comparisons = A + B (T0 and 1 week post exposure)

##subset each ASV table for each dose, ignore dose
#split map
meta_split<-split(eeid.meta, eeid.meta$Temperature)

#make temp-specifc tables, temperature=6, ignores dose/time
asv.temp6<-asv.tbl[,names(asv.tbl) %in% meta_split$T6$SampleID]
asv.temp6$OTU<-row.names(asv.tbl)
asv.temp6<-asv.temp6 %>% relocate(OTU)
write.table(asv.temp6, 'split_asv_tables/asv.tmp6.txt', quote=F, sep="\t", row.names = F)

#make temp-specifc tables, temperature=14, ignores dose/time
asv.temp14<-asv.tbl[,names(asv.tbl) %in% meta_split$T14$SampleID]
asv.temp14$OTU<-row.names(asv.tbl)
asv.temp14<-asv.temp14 %>% relocate(OTU)
write.table(asv.temp14, 'split_asv_tables/asv.tmp14.txt', quote=F, sep="\t", row.names = F)

#make temp-specifc tables, temperature=22, ignores dose/time
asv.temp22<-asv.tbl[,names(asv.tbl) %in% meta_split$T22$SampleID]
asv.temp22$OTU<-row.names(asv.tbl)
asv.temp22<-asv.temp22 %>% relocate(OTU)
write.table(asv.temp22, 'split_asv_tables/asv.tmp22.txt', quote=F, sep="\t", row.names = F)


#subset each ASV table for each dose, ignore temperature/time
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
```

## Subset table for each dose/temp/time combo
```
setwd("Github/EEID_analysis/")
library(tidyverse)
#read ASV table
asv.tbl<-read.delim("otu_table.txt", header=T, row.names=1)

#read metadata
eeid.meta<-read.delim("Data/Metadata_NSFEEID_16SRuns_1234Merged_plusExpData.txt", header=T)

#Temperatures = 6, 14, 22
#doses = "Control" "5x10.3"  "5x10.4"  "5x10.5"  "5x10.6"  "C"
#TimeWeekCat comparisons = A + B (T0 and 1 week post exposure)


#subset table for each dose/temp/time combo

#split meta by time/temperature/dose
meta_split2<-split(eeid.meta, list(eeid.meta$TimeWeekCat,  eeid.meta$Dose, eeid.meta$Temperature))

#control dose
asv.dosectrl.6B<-asv.tbl[, names(asv.tbl) %in% meta_split2$B.dControl.T6$SampleID]
asv.dosectrl.6A<-asv.tbl[, names(asv.tbl) %in% meta_split2$A.dControl.T6$SampleID]
asv.dosectrl.14A<-asv.tbl[, names(asv.tbl) %in% meta_split2$A.dControl.T14$SampleID]
asv.dosectrl.14B<-asv.tbl[, names(asv.tbl) %in% meta_split2$B.dControl.T14$SampleID]
asv.dosectrl.22A<-asv.tbl[, names(asv.tbl) %in% meta_split2$A.dControl.T22$SampleID]
asv.dosectrl.22B<-asv.tbl[, names(asv.tbl) %in% meta_split2$B.dControl.T22$SampleID]

#create column for ASV names
asv.dosectrl.6A$OTU<-row.names(asv.dosectrl.6A)
asv.dosectrl.6B$OTU<-row.names(asv.dosectrl.6B)
asv.dosectrl.14A$OTU<-row.names(asv.dosectrl.14A)
asv.dosectrl.14B$OTU<-row.names(asv.dosectrl.14B)
asv.dosectrl.22A$OTU<-row.names(asv.dosectrl.22A)
asv.dosectrl.22B$OTU<-row.names(asv.dosectrl.22B)

#move ASV name to its own column
asv.dosectrl.6A<-asv.dosectrl.6A %>% relocate(OTU)
asv.dosectrl.6B<-asv.dosectrl.6B %>% relocate(OTU)
asv.dosectrl.14A<-asv.dosectrl.14A %>% relocate(OTU)
asv.dosectrl.14B<-asv.dosectrl.14B %>% relocate(OTU)
asv.dosectrl.22A<-asv.dosectrl.22A %>% relocate(OTU)
asv.dosectrl.22B<-asv.dosectrl.22B %>% relocate(OTU)

#write data as a table
write.table(asv.dosectrl.6A, 'split_asv_tables/asv.dosectrl.6A.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dosectrl.6B, 'split_asv_tables/asv.dosectrl.6B.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dosectrl.14A, 'split_asv_tables/asv.dosectrl.14A.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dosectrl.14B, 'split_asv_tables/asv.dosectrl.14B.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dosectrl.22A, 'split_asv_tables/asv.dosectrl.22A.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dosectrl.22B, 'split_asv_tables/asv.dosectrl.22B.txt', quote=F, sep="\t", row.names = F)

#dose 10^3
asv.dose103.6B<-asv.tbl[, names(asv.tbl) %in% meta_split2$B.d5x10.3.T6$SampleID]
asv.dose103.6A<-asv.tbl[, names(asv.tbl) %in% meta_split2$A.d5x10.3.T6$SampleID]
asv.dose103.14A<-asv.tbl[, names(asv.tbl) %in% meta_split2$A.d5x10.3.T14$SampleID]
asv.dose103.14B<-asv.tbl[, names(asv.tbl) %in% meta_split2$B.d5x10.3.T14$SampleID]
asv.dose103.22A<-asv.tbl[, names(asv.tbl) %in% meta_split2$A.d5x10.3.T22$SampleID]
asv.dose103.22B<-asv.tbl[, names(asv.tbl) %in% meta_split2$B.d5x10.3.T22$SampleID]

#create column for ASV names
asv.dose103.6A$OTU<-row.names(asv.dose103.6A)
asv.dose103.6B$OTU<-row.names(asv.dose103.6B)
asv.dose103.14A$OTU<-row.names(asv.dose103.14A)
asv.dose103.14B$OTU<-row.names(asv.dose103.14B)
asv.dose103.22A$OTU<-row.names(asv.dose103.22A)
asv.dose103.22B$OTU<-row.names(asv.dose103.22B)

#move ASV name to its own column
asv.dose103.6A<-asv.dose103.6A %>% relocate(OTU)
asv.dose103.6B<-asv.dose103.6B %>% relocate(OTU)
asv.dose103.14A<-asv.dose103.14A %>% relocate(OTU)
asv.dose103.14B<-asv.dose103.14B %>% relocate(OTU)
asv.dose103.22A<-asv.dose103.22A %>% relocate(OTU)
asv.dose103.22B<-asv.dose103.22B %>% relocate(OTU)

#write data as a table
write.table(asv.dose103.6A, 'split_asv_tables/asv.dose103.6A.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose103.6B, 'split_asv_tables/asv.dose103.6B.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose103.14A, 'split_asv_tables/asv.dose103.14A.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose103.14B, 'split_asv_tables/asv.dose103.14B.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose103.22A, 'split_asv_tables/asv.dose103.22A.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose103.22B, 'split_asv_tables/asv.dose103.22B.txt', quote=F, sep="\t", row.names = F)

#dose 10^4
asv.dose104.6B<-asv.tbl[, names(asv.tbl) %in% meta_split2$B.d5x10.4.T6$SampleID]
asv.dose104.6A<-asv.tbl[, names(asv.tbl) %in% meta_split2$A.d5x10.4.T6$SampleID]
asv.dose104.14A<-asv.tbl[, names(asv.tbl) %in% meta_split2$A.d5x10.4.T14$SampleID]
asv.dose104.14B<-asv.tbl[, names(asv.tbl) %in% meta_split2$B.d5x10.4.T14$SampleID]
asv.dose104.22A<-asv.tbl[, names(asv.tbl) %in% meta_split2$A.d5x10.4.T22$SampleID]
asv.dose104.22B<-asv.tbl[, names(asv.tbl) %in% meta_split2$B.d5x10.4.T22$SampleID]

#create column for ASV names
asv.dose104.6A$OTU<-row.names(asv.dose104.6A)
asv.dose104.6B$OTU<-row.names(asv.dose104.6B)
asv.dose104.14A$OTU<-row.names(asv.dose104.14A)
asv.dose104.14B$OTU<-row.names(asv.dose104.14B)
asv.dose104.22A$OTU<-row.names(asv.dose104.22A)
asv.dose104.22B$OTU<-row.names(asv.dose104.22B)

#move ASV name to its own column
asv.dose104.6A<-asv.dose104.6A %>% relocate(OTU)
asv.dose104.6B<-asv.dose104.6B %>% relocate(OTU)
asv.dose104.14A<-asv.dose104.14A %>% relocate(OTU)
asv.dose104.14B<-asv.dose104.14B %>% relocate(OTU)
asv.dose104.22A<-asv.dose104.22A %>% relocate(OTU)
asv.dose104.22B<-asv.dose104.22B %>% relocate(OTU)

#write data as a table
write.table(asv.dose104.6A, 'split_asv_tables/asv.dose104.6A.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose104.6B, 'split_asv_tables/asv.dose104.6B.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose104.14A, 'split_asv_tables/asv.dose104.14A.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose104.14B, 'split_asv_tables/asv.dose104.14B.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose104.22A, 'split_asv_tables/asv.dose104.22A.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose104.22B, 'split_asv_tables/asv.dose104.22B.txt', quote=F, sep="\t", row.names = F)

#dose 10^5
asv.dose105.6B<-asv.tbl[, names(asv.tbl) %in% meta_split2$B.d5x10.5.T6$SampleID]
asv.dose105.6A<-asv.tbl[, names(asv.tbl) %in% meta_split2$A.d5x10.5.T6$SampleID]
asv.dose105.14A<-asv.tbl[, names(asv.tbl) %in% meta_split2$A.d5x10.5.T14$SampleID]
asv.dose105.14B<-asv.tbl[, names(asv.tbl) %in% meta_split2$B.d5x10.5.T14$SampleID]
asv.dose105.22A<-asv.tbl[, names(asv.tbl) %in% meta_split2$A.d5x10.5.T22$SampleID]
asv.dose105.22B<-asv.tbl[, names(asv.tbl) %in% meta_split2$B.d5x10.5.T22$SampleID]

#create column for ASV names
asv.dose105.6A$OTU<-row.names(asv.dose105.6A)
asv.dose105.6B$OTU<-row.names(asv.dose105.6B)
asv.dose105.14A$OTU<-row.names(asv.dose105.14A)
asv.dose105.14B$OTU<-row.names(asv.dose105.14B)
asv.dose105.22A$OTU<-row.names(asv.dose105.22A)
asv.dose105.22B$OTU<-row.names(asv.dose105.22B)

#move ASV name to its own column
asv.dose105.6A<-asv.dose105.6A %>% relocate(OTU)
asv.dose105.6B<-asv.dose105.6B %>% relocate(OTU)
asv.dose105.14A<-asv.dose105.14A %>% relocate(OTU)
asv.dose105.14B<-asv.dose105.14B %>% relocate(OTU)
asv.dose105.22A<-asv.dose105.22A %>% relocate(OTU)
asv.dose105.22B<-asv.dose105.22B %>% relocate(OTU)

#write data as a table
write.table(asv.dose105.6A, 'split_asv_tables/asv.dose105.6A.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose105.6B, 'split_asv_tables/asv.dose105.6B.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose105.14A, 'split_asv_tables/asv.dose105.14A.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose105.14B, 'split_asv_tables/asv.dose105.14B.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose105.22A, 'split_asv_tables/asv.dose105.22A.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose105.22B, 'split_asv_tables/asv.dose105.22B.txt', quote=F, sep="\t", row.names = F)

#dose 10^6
asv.dose106.6B<-asv.tbl[, names(asv.tbl) %in% meta_split2$B.d5x10.6.T6$SampleID]
asv.dose106.6A<-asv.tbl[, names(asv.tbl) %in% meta_split2$A.d5x10.6.T6$SampleID]
asv.dose106.14A<-asv.tbl[, names(asv.tbl) %in% meta_split2$A.d5x10.6.T14$SampleID]
asv.dose106.14B<-asv.tbl[, names(asv.tbl) %in% meta_split2$B.d5x10.6.T14$SampleID]
asv.dose106.22A<-asv.tbl[, names(asv.tbl) %in% meta_split2$A.d5x10.6.T22$SampleID]
asv.dose106.22B<-asv.tbl[, names(asv.tbl) %in% meta_split2$B.d5x10.6.T22$SampleID]

#create column for ASV names
asv.dose106.6A$OTU<-row.names(asv.dose106.6A)
asv.dose106.6B$OTU<-row.names(asv.dose106.6B)
asv.dose106.14A$OTU<-row.names(asv.dose106.14A)
asv.dose106.14B$OTU<-row.names(asv.dose106.14B)
asv.dose106.22A$OTU<-row.names(asv.dose106.22A)
asv.dose106.22B$OTU<-row.names(asv.dose106.22B)

#move ASV name to its own column
asv.dose106.6A<-asv.dose106.6A %>% relocate(OTU)
asv.dose106.6B<-asv.dose106.6B %>% relocate(OTU)
asv.dose106.14A<-asv.dose106.14A %>% relocate(OTU)
asv.dose106.14B<-asv.dose106.14B %>% relocate(OTU)
asv.dose106.22A<-asv.dose106.22A %>% relocate(OTU)
asv.dose106.22B<-asv.dose106.22B %>% relocate(OTU)

#write data as a table
write.table(asv.dose106.6A, 'split_asv_tables/asv.dose106.6A.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose106.6B, 'split_asv_tables/asv.dose106.6B.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose106.14A, 'split_asv_tables/asv.dose106.14A.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose106.14B, 'split_asv_tables/asv.dose106.14B.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose106.22A, 'split_asv_tables/asv.dose106.22A.txt', quote=F, sep="\t", row.names = F)
write.table(asv.dose106.22B, 'split_asv_tables/asv.dose106.22B.txt', quote=F, sep="\t", row.names = F)
```

## Getting networks for 10^3 through time
```
setwd("./Github/EEID_analysis/")
library(tidyverse)
#read ASV table
asv.tbl<-read.delim("otu_table.txt", header=T, row.names=1)

#read metadata
eeid.meta<-read.delim("Data/Metadata_NSFEEID_16SRuns_1234Merged_plusExpData.txt", header=T)

#need dose of 10^3 only.
#subset meta data
meta_103<-eeid.meta[which(eeid.meta$Dose=='d5x10.3'),]
meta_103s<-split(meta_103, list(meta_103$TimeWeekCat, meta_103$Temperature))

dose103.T6.A<-asv.tbl[,names(asv.tbl) %in% meta_103s$A.T6$SampleID]
dose103.T6.B<-asv.tbl[,names(asv.tbl) %in% meta_103s$B.T6$SampleID]
dose103.T6.C<-asv.tbl[,names(asv.tbl) %in% meta_103s$C.T6$SampleID]
dose103.T6.D<-asv.tbl[,names(asv.tbl) %in% meta_103s$D.T6$SampleID]
dose103.T6.E<-asv.tbl[,names(asv.tbl) %in% meta_103s$E.T6$SampleID]
dose103.T6.F<-asv.tbl[,names(asv.tbl) %in% meta_103s$F.T6$SampleID]
dose103.T6.G<-asv.tbl[,names(asv.tbl) %in% meta_103s$G.T6$SampleID]
dose103.T6.H<-asv.tbl[,names(asv.tbl) %in% meta_103s$H.T6$SampleID]
#I->L don't have enough samples


dose103.T6.A$OTU<-row.names(dose103.T6.A)
dose103.T6.B$OTU<-row.names(dose103.T6.B)
dose103.T6.C$OTU<-row.names(dose103.T6.C)
dose103.T6.D$OTU<-row.names(dose103.T6.D)
dose103.T6.E$OTU<-row.names(dose103.T6.E)
dose103.T6.F$OTU<-row.names(dose103.T6.F)
dose103.T6.G$OTU<-row.names(dose103.T6.G)
dose103.T6.H$OTU<-row.names(dose103.T6.H)


dose103.T6.A<- dose103.T6.A %>% relocate(OTU)
dose103.T6.B<- dose103.T6.B %>% relocate(OTU)
dose103.T6.C<- dose103.T6.C %>% relocate(OTU)
dose103.T6.D<- dose103.T6.D %>% relocate(OTU)
dose103.T6.E<- dose103.T6.E %>% relocate(OTU)
dose103.T6.F<- dose103.T6.F %>% relocate(OTU)
dose103.T6.G<- dose103.T6.G %>% relocate(OTU)
dose103.T6.H<- dose103.T6.H %>% relocate(OTU)

  
write.table(dose103.T6.A, 'split_asv_tables_time/dose1063.T6.A.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T6.B, 'split_asv_tables_time/dose1063.T6.B.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T6.C, 'split_asv_tables_time/dose1063.T6.C.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T6.D, 'split_asv_tables_time/dose1063.T6.D.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T6.E, 'split_asv_tables_time/dose1063.T6.E.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T6.F, 'split_asv_tables_time/dose1063.T6.F.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T6.G, 'split_asv_tables_time/dose1063.T6.G.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T6.H, 'split_asv_tables_time/dose1063.T6.H.txt', quote=F, sep="\t", row.names = F)

#for 14C
dose103.T14.A<-asv.tbl[,names(asv.tbl) %in% meta_103s$A.T14$SampleID]
dose103.T14.B<-asv.tbl[,names(asv.tbl) %in% meta_103s$B.T14$SampleID]
dose103.T14.C<-asv.tbl[,names(asv.tbl) %in% meta_103s$C.T14$SampleID]
dose103.T14.D<-asv.tbl[,names(asv.tbl) %in% meta_103s$D.T14$SampleID]
dose103.T14.E<-asv.tbl[,names(asv.tbl) %in% meta_103s$E.T14$SampleID]
dose103.T14.F<-asv.tbl[,names(asv.tbl) %in% meta_103s$F.T14$SampleID]
dose103.T14.G<-asv.tbl[,names(asv.tbl) %in% meta_103s$G.T14$SampleID]
dose103.T14.H<-asv.tbl[,names(asv.tbl) %in% meta_103s$H.T14$SampleID]
#I->L don't have enough samples


dose103.T14.A$OTU<-row.names(dose103.T14.A)
dose103.T14.B$OTU<-row.names(dose103.T14.B)
dose103.T14.C$OTU<-row.names(dose103.T14.C)
dose103.T14.D$OTU<-row.names(dose103.T14.D)
dose103.T14.E$OTU<-row.names(dose103.T14.E)
dose103.T14.F$OTU<-row.names(dose103.T14.F)
dose103.T14.G$OTU<-row.names(dose103.T14.G)
dose103.T14.H$OTU<-row.names(dose103.T14.H)


dose103.T14.A<- dose103.T14.A %>% relocate(OTU)
dose103.T14.B<- dose103.T14.B %>% relocate(OTU)
dose103.T14.C<- dose103.T14.C %>% relocate(OTU)
dose103.T14.D<- dose103.T14.D %>% relocate(OTU)
dose103.T14.E<- dose103.T14.E %>% relocate(OTU)
dose103.T14.F<- dose103.T14.F %>% relocate(OTU)
dose103.T14.G<- dose103.T14.G %>% relocate(OTU)
dose103.T14.H<- dose103.T14.H %>% relocate(OTU)


write.table(dose103.T14.A, 'split_asv_tables_time/dose1063.T14.A.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T14.B, 'split_asv_tables_time/dose1063.T14.B.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T14.C, 'split_asv_tables_time/dose1063.T14.C.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T14.D, 'split_asv_tables_time/dose1063.T14.D.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T14.E, 'split_asv_tables_time/dose1063.T14.E.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T14.F, 'split_asv_tables_time/dose1063.T14.F.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T14.G, 'split_asv_tables_time/dose1063.T14.G.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T14.H, 'split_asv_tables_time/dose1063.T14.H.txt', quote=F, sep="\t", row.names = F)



dose103.T22.A<-asv.tbl[,names(asv.tbl) %in% meta_103s$A.T22$SampleID]
dose103.T22.B<-asv.tbl[,names(asv.tbl) %in% meta_103s$B.T22$SampleID]
dose103.T22.C<-asv.tbl[,names(asv.tbl) %in% meta_103s$C.T22$SampleID]
dose103.T22.D<-asv.tbl[,names(asv.tbl) %in% meta_103s$D.T22$SampleID]
dose103.T22.E<-asv.tbl[,names(asv.tbl) %in% meta_103s$E.T22$SampleID]
dose103.T22.F<-asv.tbl[,names(asv.tbl) %in% meta_103s$F.T22$SampleID]
dose103.T22.G<-asv.tbl[,names(asv.tbl) %in% meta_103s$G.T22$SampleID]
dose103.T22.H<-asv.tbl[,names(asv.tbl) %in% meta_103s$H.T22$SampleID]
#I->L don't have enough samples


dose103.T22.A$OTU<-row.names(dose103.T22.A)
dose103.T22.B$OTU<-row.names(dose103.T22.B)
dose103.T22.C$OTU<-row.names(dose103.T22.C)
dose103.T22.D$OTU<-row.names(dose103.T22.D)
dose103.T22.E$OTU<-row.names(dose103.T22.E)
dose103.T22.F$OTU<-row.names(dose103.T22.F)
dose103.T22.G$OTU<-row.names(dose103.T22.G)
dose103.T22.H$OTU<-row.names(dose103.T22.H)


dose103.T22.A<- dose103.T22.A %>% relocate(OTU)
dose103.T22.B<- dose103.T22.B %>% relocate(OTU)
dose103.T22.C<- dose103.T22.C %>% relocate(OTU)
dose103.T22.D<- dose103.T22.D %>% relocate(OTU)
dose103.T22.E<- dose103.T22.E %>% relocate(OTU)
dose103.T22.F<- dose103.T22.F %>% relocate(OTU)
dose103.T22.G<- dose103.T22.G %>% relocate(OTU)
dose103.T22.H<- dose103.T22.H %>% relocate(OTU)


write.table(dose103.T22.A, 'split_asv_tables_time/dose1063.T22.A.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T22.B, 'split_asv_tables_time/dose1063.T22.B.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T22.C, 'split_asv_tables_time/dose1063.T22.C.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T22.D, 'split_asv_tables_time/dose1063.T22.D.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T22.E, 'split_asv_tables_time/dose1063.T22.E.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T22.F, 'split_asv_tables_time/dose1063.T22.F.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T22.G, 'split_asv_tables_time/dose1063.T22.G.txt', quote=F, sep="\t", row.names = F)
write.table(dose103.T22.H, 'split_asv_tables_time/dose1063.T22.H.txt', quote=F, sep="\t", row.names = F)
```

Networks run in CoNet app in Cytoscape. Spearman correlations, p<0.05, Benjamini-Hochberg correct p-value, OTU abundance >2, Spearman rho>0.7

## Plot network statistics
```
#read in change in network statistics
stats<-read.delim("network_stats_delta.txt", header=T)

#load ggplot
library(ggplot2)

#plot each metric
ggplot(stats, aes(as.numeric(Dose), No_nodes))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Temp)+
  ylab("Change in # nodes (T1-T0)")+
  xlab("Dose")+
  geom_smooth()

ggplot(stats, aes(as.numeric(Dose), No_edges))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Temp)+
  ylab("Change in # edges (T1-T0)")+
  xlab("Dose")+
  geom_smooth()

ggplot(stats, aes(as.numeric(Dose), Avg_num_neighbors))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Temp)+
  ylab("Change in avg # edges (T1-T0)")+
  xlab("Dose")+
  geom_smooth()

ggplot(stats, aes(as.numeric(Dose), Clus_coeffi))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Temp)+
  ylab("Change in Clustering Coefficient (T1-T0)")+
  xlab("Dose")+
  geom_smooth()
```

## Network statistics through time at dose = 10^3
```
#read in timecourse network stats
stats<-read.delim("network_stats_timecourse.txt", header=T)

ggplot(stats, aes(Time, No_nodes))+
  facet_wrap(~Temp)+
  geom_point()+
  theme_bw()+
  ylab("Number of network nodes")+
  xlab("Time")+
  geom_smooth()

ggplot(stats, aes(Time, No_edges))+
  facet_wrap(~Temp)+
  geom_point()+
  theme_bw()+
  ylab("Number of network edges")+
  xlab("Time")+
  geom_smooth()

ggplot(stats, aes(Time, Avg_num_neighbors))+
  facet_wrap(~Temp)+
  geom_point()+
  theme_bw()+
  ylab("Average number of network neighbors")+
  xlab("Time")+
  geom_smooth()

ggplot(stats, aes(Time, Clus_coeffi))+
  facet_wrap(~Temp)+
  geom_point()+
  theme_bw()+
  ylab("Clustering Coefficient")+
  xlab("Time")+
  geom_smooth()
```


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
