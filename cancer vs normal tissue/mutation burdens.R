library(tidyverse)
library(viridis)

mutation_burdens<-read.csv("mutation burdens.csv")

ggplot(mutation_burdens,aes(colour=sample,fill=sample,x=mb))+facet_wrap(~paste(cancer,", ",mutations,sep=""),nrow=3,scales="free_x")+labs(x="mutation burden",colour="",fill="")+guides(fill="none",colour="none")+geom_density(alpha=0.2)+scale_x_continuous(trans="log10")+theme_classic()+
  scale_fill_manual(values=c("black",viridis(n=1,begin=0.8,end=0.8,option="G")))+scale_colour_manual(values=c("black",viridis(n=1,begin=0.8,end=0.8,option="G")))
