## Random forest of KOs
Notes: most common response to temp/dose is metabolism...wicked cool

```
library("randomForest")
library("plyr")
library("rfUtilities")
library("caret")
setwd("/Users/patty/OneDrive/Documents/Github/EEID_analysis/")

#read in mapping file
meta<-read.delim("Data/Metadata_NSFEEID_16SRuns_1234Merged_plusExpData.txt", header=T)

#read in KO table
s16<-read.delim("./picrust_out_norm/KO_metagenome_out/pred_metagenome_unstrat.tsv", header=T, row.names=1)
dim(s16)
#7301 KOs by 499 samples


#how many nonzero counts?
otu_nonzero_counts<-apply(s16, 1, function(y) sum(length(which(y > 0))))
hist(otu_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of KOs", xlab="Number of Non-Zero Values")

#remove rare KO
remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

otu_table_rare_removed <- remove_rare(table=s16, cutoff_pro=0.1)
dim(otu_table_rare_removed)
#6220 OTUs by 499 samples

#scale data
otu_table_scaled <- scale(otu_table_rare_removed, center = TRUE, scale = TRUE)  

#fix dataframe to add metadata
otu_table_scaled_treatment <- data.frame(t(otu_table_scaled))
otu_table_scaled_treatment$SampleID<-row.names(otu_table_scaled_treatment)

#select only dose/temp columns from metadata
meta_sub<-as.data.frame(meta[,c(1,14)])
otu_table_scaled_treatment<-merge(otu_table_scaled_treatment, meta_sub, by="SampleID")
otu_table_scaled_treatment<-otu_table_scaled_treatment[,-1]

#set DoseTemp to factor if character, otherwise model won't run son
otu_table_scaled_treatment$DoseTemp<-as.factor(otu_table_scaled_treatment$DoseTemp)

#run model
set.seed(515)
RF_treatment_classify<-randomForest(x=otu_table_scaled_treatment[,1:(ncol(otu_table_scaled_treatment)-1)], y=otu_table_scaled_treatment[ , ncol(otu_table_scaled_treatment)] , ntree=100, importance=TRUE, proximities=TRUE)

#permutation test to test for signficance of model
RF_treatment_classify_sig<-rf.significance(x=RF_treatment_classify, xdata=otu_table_scaled_treatment[,1:(ncol(otu_table_scaled_treatment)-1)], nperm=100, ntree=100)

#identifying important features
RF_state_classify_imp <- as.data.frame(RF_treatment_classify$importance)
RF_state_classify_imp$features <- rownames( RF_state_classify_imp)
RF_state_classify_imp_sorted <- arrange( RF_state_classify_imp  , desc(MeanDecreaseAccuracy)  )
barplot(RF_state_classify_imp_sorted$MeanDecreaseAccuracy, ylab="Mean Decrease in Accuracy (Variable Importance)", main="RF Classification Variable Importance Distribution")

#top 50 features
barplot(RF_state_classify_imp_sorted[1:50,"MeanDecreaseAccuracy"], las=2, names.arg=RF_state_classify_imp_sorted[1:50,"features"] , ylab="Mean Decrease in Accuracy (Variable Importance)", main="Classification RF")  

#write RF to file becauase it takes a long time to run and ain't nobody got time for that
save.image("C:/Users/patty/OneDrive/Documents/Github/EEID_analysis/random_forest_KO.RData")
```

## Plot mean abundance of top 50 most important KOs
```
top50_feat<-as.data.frame(RF_state_classify_imp_sorted$features[1:50])
names(top50_feat)<-c("KO")
str(top50_feat)

#read in KO info, extract taxonomy
tax<-read.delim("full_kegg.txt", header=T)

#change original KO table to relative abundance
library(vegan)
s16<-decostand(s16, method = 'total')

#check to make sure rel abund cause I got paranoia about these thangs
rowSums(s16)

#extract top 50 from OTU table
s16.top50<-s16[rownames(s16) %in% top50_feat$KO,]
dim(s16.top50)  
#50 by 499
s16.top50$KO<-row.names(s16.top50)

#get mean relative abundance of each OTU in the 4 cats
library(reshape)
top50_m<-melt(s16.top50)
names(top50_m)<-c("KO", "SampleID", 'Rel_abund')
library(plyr)
top50_m<-merge(top50_m, meta, by='SampleID')
top50_sum<-ddply(top50_m, c('KO', "DoseTemp"), summarize, mean=mean(Rel_abund), sd=sd(Rel_abund), n=length(Rel_abund), se=sd/n)
top50_sum<-merge(top50_sum, tax, by='KO')

```

## plot the data
```
#define a color pallet n=50
colour_tax<-rainbow(50, s=1, v=1)[sample(1:50,50)]

#plot it
library(ggplot2)
ggplot(top50_sum, aes(Level1, mean, fill=Level2))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=colour_tax)+
  facet_wrap(~DoseTemp, ncol = 4)+
  theme_bw()+
  coord_flip()+
  ggtitle("Top 50 KOs in distinguishing categories")+
  guides(fill=guide_legend(ncol=1))+
  ylab("Mean Relative Abundance")+
  xlab("")+
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se, width=0.2), stat="identity")
```