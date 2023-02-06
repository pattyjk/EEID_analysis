## Neutral models of beta diversity
```
#models run according to pattyjk/beta_modeling repository

#read in neutral model output
model.out<-read.delim("neutral_model_outpu.txt", header=T)

#load libraries
library(ggplot2)

#plot sloan neutral models
ggplot(model.out, aes(Dose, Sloan_neutral_null_dev))+
  geom_point()+
  theme_bw()+
  xlab("Bd Dose")+
  scale_x_log10()+
  ylab("Abundance Null Deviation")+
  facet_wrap(~Temp)

#plot Tucker neutral models
ggplot(model.out, aes(Dose, Tucket_neutral_null_dev))+
  geom_point()+
  theme_bw()+
  xlab("Bd Dose")+
  scale_x_log10()+
  ylab("Abundance Null Deviation")+
  facet_wrap(~Temp)

#R2 values higher in control than infected samples as well
```
