## Random forest of KOs

```
library("randomForest")
library("plyr")
library("rfUtilities")
library("caret")
library(ggplot2)

#read in mapping file
meta<-read.delim("Data/Metadata_NSFEEID_16SRuns_1234Merged_plusExpData.txt", header=T)

#read in KO table
s16<-read.delim("./picrust_out_norm/KO_metagenome_out/pred_metagenome_unstrat.tsv", header=T, row.names=1)
dim(s16)
#7301 KOs by 499 samples


#how many nonzero counts?
otu_nonzero_counts<-apply(s16, 1, function(y) sum(length(which(y > 0))))

#plot histogram of the frequency of counts
hist(otu_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of KOs", xlab="Number of Non-Zero Values")

#write function to remove rare KO
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

#use function to remove rare KO
otu_table_rare_removed <- remove_rare(table=s16, cutoff_pro=0.1)
dim(otu_table_rare_removed)
#6220 KOs by 499 samples, removed ~1000 KO that were rare (<1 percent)

#remove time points that aren't A/B (0 & 7 days post exposure)
meta_AB<-meta[which(meta$TimeWeekCat == "A" | meta$TimeWeekCat == "B"),]
meta_AB<-meta_AB[,1:2]

otu_table_rare_removed <-t(otu_table_rare_removed)
otu_table_rare_removed2 <-otu_table_rare_removed
rows_otu <-as.data.frame(row.names(otu_table_rare_removed2))
names(rows_otu)<-c("SampleID")
otu_table_rare_removed2<-cbind(otu_table_rare_removed2, rows_otu)
otu_table_rare_removed3<-merge(otu_table_rare_removed2, meta_AB, by='SampleID', all.x=F)
row.names(otu_table_rare_removed3)<-otu_table_rare_removed3$SampleID
otu_table_rare_removed3<-otu_table_rare_removed3[,-1]
otu_table_rare_removed3$BarcodeSequence=NULL
otu_table_rare_removed3<-as.data.frame(t(otu_table_rare_removed3))

#scale data
otu_table_scaled <- scale(otu_table_rare_removed3, center = TRUE, scale = TRUE)  

#fix data frame to add metadata
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

#top 100 features
barplot(RF_state_classify_imp_sorted[1:100,"MeanDecreaseAccuracy"], las=2, names.arg=RF_state_classify_imp_sorted[1:100,"features"] , ylab="Mean Decrease in Accuracy (Variable Importance)", main="Classification RF")  

#write RF to file because it takes a long time to run and ain't nobody got time for that
save.image("C:/Users/patty/OneDrive/Documents/Github/EEID_analysis/random_forest_KO.RData")
```

## Plot mean abundance of top 50 most important KOs
```
top50_feat<-as.data.frame(RF_state_classify_imp_sorted$features[1:30])
names(top50_feat)<-c("KO")
str(top50_feat)

#read in KO info, extract taxonomy
tax<-read.delim("full_kegg.txt", header=T)

#change original KO table to relative abundance
library(vegan)
s16<-decostand(otu_table_rare_removed3, method = 'total')

#check to make sure rel abund cause I got paranoia about these thangs
rowSums(s16)

#extract top 50 from KO table
s16.top50<-s16[rownames(s16) %in% top50_feat$KO,]
dim(s16.top50)  
#100 by 252
s16.top50$KO<-row.names(s16.top50)

#get mean relative abundance of each OTU in the 4 cats
library(reshape2)
library(plyr)
top50_m<-melt(s16.top50)
names(top50_m)<-c("KO", "SampleID", 'Rel_abund')
top50_m<-merge(top50_m, meta, by='SampleID', all.y=F)
top50_sum<-ddply(top50_m, c('KO', "Temperature", "TimeWeekCat", "Dose"), summarize, mean=mean(Rel_abund), sd=sd(Rel_abund), n=length(Rel_abund), se=sd/n)
top50_sum<-merge(top50_sum, tax, by='KO')

```

## plot the data
```
#define a color pallet n=50
colour_tax<-rainbow(50, s=1, v=1)[sample(1:50,50)]

#calculate delta
top50_sum_split<-split(top50_sum, top50_sum$TimeWeekCat)
sum_split<-top50_sum_split$A[,1:4]
sum_split$delta<-top50_sum_split$B$mean - top50_sum_split$A$mean
sum_split<-merge(sum_split, tax, by='KO')

meta2<-as.data.frame(meta$Dose2, meta$Dose)
meta2$Dose<-row.names(meta2)
names(meta2)<-c("Dose2", "Dose")
sum_split<-merge(sum_split, meta2, by='Dose')

sum_split2<-split(sum_split, sum_split$Temperature)

#define a pallet
pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")

#plot it
library(ggplot2)
ggplot(sum_split2$T14, aes(as.numeric(Dose2), delta, color=Level2))+
  geom_point()+
  geom_line()+
  facet_wrap(~Product)+
  scale_color_manual(values=pal)+
  theme_bw()+
  ylab("Change in Mean Relative Abundance (T1-T0)")+
  ggtitle("Top 30 KOs in distinguishing categories- T14")+
  guides(fill=guide_legend(ncol=1))+
  xlab("Dose")

ggplot(sum_split2$T6, aes(as.numeric(Dose2), delta, color=Level2))+
  geom_point()+
  scale_color_manual(values=pal)+
  geom_line()+
  facet_wrap(~Product)+
  theme_bw()+
  ylab("Change in Mean Relative Abundance (T1-T0)")+
  ggtitle("Top 30 KOs in distinguishing categories- T6")+
  guides(fill=guide_legend(ncol=1))+
  xlab("Dose")

ggplot(sum_split2$T22, aes(as.numeric(Dose2), delta, color=Level2))+
  geom_point()+
  geom_line()+
  scale_color_manual(values=pal)+
  facet_wrap(~Product)+
  theme_bw()+
  ylab("Change in Mean Relative Abundance (T1-T0)")+
  ggtitle("Top 30 KOs in distinguishing categories- T2")+
  guides(fill=guide_legend(ncol=1))+
  xlab("Dose")
  
###
rf_ko<-read.delim("~/Documents/GitHub/EEID_analysis/RF_ko_delta.txt", header=T)
View(rf_ko)
library(dplyr)
library(plyr)
library(reshape2)
library(ggplot2)

#calculate correlations between dose/delta abundance
rf_split<-split(rf_ko, rf_ko$Temperature)

#write function to do it
func <- function(xx)
{
  return(data.frame(COR = cor(xx$Dose, xx$delta)))
}

#calculate it
cor14<-ddply(rf_split$T14, .(KO), func)
cor22<-ddply(rf_split$T22, .(KO), func)
  cor6<-ddply(rf_split$T6, .(KO), func)

#fix headers and bind
names(cor14)<-c("KO", "Cor14")
names(cor22)<-c("KO", "Cor22")
names(cor6)<-c("KO", "Cor6")
cors<-merge(cor14, cor22, by='KO')
cors<-merge(cors, cor6, by='KO')

#reshape data
cors_m<-melt(cors)

#add metadata
full_kegg<-read.delim("~/Documents/GitHub/EEID_analysis/full_kegg.txt", header=T)
cors_m<-merge(cors_m, full_kegg, by='KO')

#reorder table by correlation size (ascending)
cors_m<-cors_m[order(cors_m$value),]


#define a pallett
pal2<-c("#000000","#114477",  "#117744", "#AAAA44",  "#774411",  "#771122", "#41AB5D", "#252525", "#525252")
pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")

pal3<-c("#800000FF", "#000000", "#5b8fa8ff", "#725663ff", "#adb17dff", "#ffb547ff", '#d6d6ceff', '#b1746fff', "#d49464ff")
library(RColorBrewer)

cors_m$variable<-gsub("Cor22", "22", cors_m$variable)
cors_m$variable<-gsub("Cor14", "14", cors_m$variable)
cors_m$variable<-gsub("Cor6", "6", cors_m$variable)
ggplot(cors_m, aes(KO, value, color=Level1))+
  geom_point(size=2.6)+
  theme_bw()+
  coord_flip()+
  ylab("")+
  facet_wrap(~variable, scales='free')+ 
  scale_color_brewer(n="KEGG Level 1", palette = "Paired")+
  geom_hline(yintercept=0, color='black', size=1)+
  ylab("Spearman Rho")
