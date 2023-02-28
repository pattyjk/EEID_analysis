## Analyze taxonomic information of the edges
```
setwd("./Github/EEID_analysis//Network_edge_tables/")

#list files
edges<-list.files("./")

#read in all files
edge_tables<-sapply(edges, function(x) read.csv(x, header=T))

#read in tax/inhibiton file
tax_inib<-read.delim("/Users/patty/OneDrive/Documents/Github/EEID_analysis/taxa_calls_inhib.txt", header=T)

#read in metadata
meta<-read.delim("/Users/patty/OneDrive/Documents/Github/EEID_analysis/Data/Metadata_NSFEEID_16SRuns_1234Merged_plusExpData.txt", header=T)

#bind all tables together
library(data.table)
edge_tab_comb<-rbindlist(edge_tables, fill=T, use.names = T, idcol="Condition")

#split interaction column
library(tidyverse)
edge_tab_comb<-separate(edge_tab_comb, col='Label', into = c("Taxa1", "Taxa2"), sep = "->")

#add taxonomy for each cooccuring taxa
edge_tab_comb<-merge(edge_tab_comb, tax_inib, by.x='Taxa1', by.y = 'Feature.ID', all.y=F, all.x=T)
names(edge_tab_comb)<-c("Taxa1", "Condition", "cooc_method", "interaction","interactionType", "Taxa2","name", "pval", "selected", "shared.interaction" ,"shared.name", "weight", "EdgeBetweenness", "Taxon1", "Confidence1" , "database_id1" , "Inhibitory1")

edge_tab_comb<-merge(edge_tab_comb, tax_inib, by.x='Taxa2', by.y = 'Feature.ID', all.y=F, all.x=T)
names(edge_tab_comb)<-c("Taxa2", "Taxa1", "Condition", "cooc_method", "interaction", "interactionType", "name", "pval", "selected", "shared.interaction", "shared.name","weight", "EdgeBetweenness", "Taxon1", "Confidence1", "database_id1", "Inhibitory1", "Taxon2", "Confidence2", "database_id2", "Inhibitory2")

#split taxonomy columns
edge_tab_comb<-separate(edge_tab_comb, Taxon1, into = c("Kingdom1", "Phylum1", "Class1", "Order1", "Family1", "Genus1", "Species1"))
edge_tab_comb<-separate(edge_tab_comb, Taxon2, into = c("Kingdom2", "Phylum2", "Class2", "Order2", "Family2", "Genus2", "Species2"))
## errors fince due to lack of taxonomy

#create columne that identifies mismatches between orders (at the class level, 41004/41539 edges are within the same class) for all interactions
edge_tab_comb$classmatch <- edge_tab_comb$Order1 == edge_tab_comb$Order2
edge_tab_comb$classmatch<-gsub("TRUE", 1, edge_tab_comb$classmatch)
edge_tab_comb$classmatch<-gsub("FALSE", 0, edge_tab_comb$classmatch)
edge_tab_comb[is.na(edge_tab_comb)] <- 0
edge_tab_comb$classmatch<-as.numeric(edge_tab_comb$classmatch)


#get percentage of interactions for each category that are within a order
library(plyr)
edge_class_sum<-ddply(edge_tab_comb, c("Condition"), summarize, sum=sum(classmatch), n=length(classmatch), per_match=sum/n)

#read in list of samps for edges
edge_samps<-read.delim("/Users/patty/OneDrive/Documents/Github/EEID_analysis/edge_samps.txt", header=T)

#add the metadata to summary
edge_class_sum<-merge(edge_class_sum, edge_samps, by='Condition')

#plot it
library(ggplot2)
ggplot(edge_class_sum, aes(Dose, per_match, fill=Temp))+
  geom_bar(stat='identity')+
  theme_bw()+
  facet_wrap(~Temp)+
  ylab("Percent of edges within the same order")+
  xlab("Condition")+
  coord_flip()
```
