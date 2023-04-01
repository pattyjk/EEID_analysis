## Analyze taxonomic information of the edges
```

#list files & read them into R
edge_tables<- list.files("./Network_edge_tables/", full.names=T) %>% 
  lapply(read.csv, header=T) 

#read in tax/inhibiton file
tax_inib<-read.delim("./taxa_calls_inhib.txt", header=T)

#read in metadata
meta<-read.delim("./Data/Metadata_NSFEEID_16SRuns_1234Merged_plusExpData.txt", header=T)

#bind all tables together
library(data.table)
edge_tab_comb<-rbindlist(edge_tables, fill=T, use.names = T, idcol="Condition")

#filter out negative correlations
edge_tab_comb<-edge_tab_comb[-which(edge_tab_comb$interactionType =="mutualExclusion"),]

#split interaction column
library(tidyverse)
edge_tab_comb<-separate(edge_tab_comb, col='Label', into = c("Taxa1", "Taxa2"), sep = "->")

#add taxonomy for each cooccuring taxa
edge_tab_comb<-merge(edge_tab_comb, tax_inib, by.x='Taxa1', by.y = 'Feature.ID', all.y=F, all.x=T)

#rename column names
names(edge_tab_comb)<-c("Taxa1", "Condition", "cooc_method", "EdgeBetweenness", "interaction","interactionType", "Taxa2","name", "selected", "shared.interaction" ,"shared.name", "weight", "Taxon1", "Confidence1" , "database_id1", "Inhibitory1")

edge_tab_comb<-merge(edge_tab_comb, tax_inib, by.x='Taxa2', by.y = 'Feature.ID', all.y=F, all.x=T)

names(edge_tab_comb)<-c("Taxa2", "Taxa1", "Condition", "cooc_method", "EdgeBetweenness", "interaction", "interactionType", "name", "selected", "shared.interaction" ,"shared.name", "weight", "Taxon1", "Confidence1" , "database_id1", "Inhibitory1", "Taxon2", "Confidence2", "database_id2", "Inhibitory2")

#split taxonomy columns
edge_tab_comb<-separate(edge_tab_comb, Taxon1, sep=";", into = c("Kingdom1", "Phylum1", "Class1", "Order1", "Family1", "Genus1", "Species1"))

edge_tab_comb<-separate(edge_tab_comb, Taxon2, sep=";", into = c("Kingdom2", "Phylum2", "Class2", "Order2", "Family2", "Genus2", "Species2"))
## errors fine due to lack of taxonomy for some ASVs

#create columne that identifies mismatches between orders (at the class level, 41004/41539 edges are within the same class) for all interactions
edge_tab_comb$classmatch <- edge_tab_comb$Class1 == edge_tab_comb$Class2
edge_tab_comb$classmatch<-gsub("TRUE", 1, edge_tab_comb$classmatch)
edge_tab_comb$classmatch<-gsub("FALSE", 0, edge_tab_comb$classmatch)
edge_tab_comb[is.na(edge_tab_comb)] <- 0
edge_tab_comb$classmatch<-as.numeric(edge_tab_comb$classmatch)


#get percentage of interactions for each category that are within a order
library(plyr)
edge_class_sum<-ddply(edge_tab_comb, c("Condition"), summarize, sum=sum(classmatch), n=length(classmatch), per_match=sum/n)

#rename samples 
edge_class_sum$Condition<-gsub("./", "", edge_class_sum$Condition)
edge_class_sum$Condition<-gsub('"', "", edge_class_sum$Condition)

#read in list of samps for edges
edge_samps<-read.delim("./edge_samps.txt", header=T)
edge_class_sum$Condition<-as.character(edge_class_sum$Condition)

#add the metadata to summary
edge_class_sum<-merge(edge_class_sum, edge_samps, by.x="Condition", by.y="Conditon")

#plot full data
library(ggplot2)
ggplot(edge_class_sum, aes(Dose, per_match, color=Time))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Temp)+
  ylab("Percent of edges within the same order")+
  xlab("Dose")
  
#split data by time
edge_class_sum_spl<-split(edge_class_sum, edge_class_sum$Time)

time_a<-edge_class_sum_spl$A
time_b<-edge_class_sum_spl$B

library(dplyr)
time_a<-arrange(time_a, Dose, Temp)
time_b<-arrange(time_b, Dose, Temp)

#subtract
delta<-time_a[,c(1,5:7)]
delta$change<-time_b$per_match-time_a$per_match

#plot delta
ggplot(delta, aes(Dose, change))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Temp)+
  ylab("Percent of edges within the same order (T1-T0)")+
  xlab("Dose")+
  geom_smooth()
```

