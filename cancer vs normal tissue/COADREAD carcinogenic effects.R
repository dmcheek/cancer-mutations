library(viridis)
library(tidyverse)
library(ggrepel)

CRC_stats<-read_csv("CRC gene summary.csv")


my_transx <- scales::trans_new('custom',
                               transform = function(x) {
                                 log(10^(-3)+x)
                               },
                               inverse = function(x) {
                                 exp(x)-10^(-3)
                               })
my_transy <- scales::trans_new('custom',
                               transform = function(x) {
                                 log(10^(-1)+x)
                               },
                               inverse = function(x) {
                                 exp(x)-10^(-1)
                               })

ggplot(CRC_stats,aes(x=freq_normal,y=freq_cancer))+
  geom_point(size=2)+
 geom_errorbar(alpha=0.1,width=0,aes(ymin=freq_cancer_l,ymax=freq_cancer_h))+
 geom_errorbarh(alpha=0.1,height=0,aes(xmin=freq_normal_l,xmax=freq_normal_h))+
  geom_text_repel(aes(label=Gene))+
  scale_y_continuous(trans=my_transy,breaks=c(0,0.1,0.4,1),labels=c(0,0.1,0.4,1),limits=c(0,1))+
  scale_x_continuous(trans=my_transx,breaks=c(0,0.001,0.004,0.01,0.1,1),labels=c(0,0.001,0.004,0.01,0.1,1),limits=c(0,0.025))+
  theme_classic()+
  labs(
    x="fraction of colon crypts with mutation",y="fraction of CRCs with mutation")

ggplot(CRC_stats,aes(y=reorder(Gene,ce_l),x=pmin(ce,2800),xmin=ce_l,xmax=pmin(ce_h,3000)))+geom_point(aes(shape=Gene=="SMAD4",size=Gene=="SMAD4"),show.legend=F)+scale_size_manual(values=c(1,4))+scale_shape_manual(values=c(15,62))+geom_errorbarh(height=0)+theme_classic()+scale_x_continuous(trans="log10")+labs(y="",x="carcinogenic effect estimate")

