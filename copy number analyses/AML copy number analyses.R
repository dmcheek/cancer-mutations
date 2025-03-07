library(tidyverse)
library(readxl)
library(viridis)
library(ggrepel)
library(weights)
library(effectsize)
library(ggpubr)
library(data.table)

bloodCNAburdens<-read_csv("blood SCNA burdens.csv")
ggplot(bloodCNAburdens,aes(x=CNAb,y=f,colour=group))+geom_line()+scale_y_continuous(trans="log10")+theme_classic()+labs(x="number of SNCAs per sample",y="fraction of samples",colour="")+scale_colour_manual(values=c(viridis(n=1,begin=0.2,end=0.2,option="mako"),viridis(n=1,begin=0.8,end=0.8,option="mako")))

AMLCHf<-read_csv("AML CH SCNAs.csv")
expand.grid(y=seq(0,sqrt(0.07),length.out=1000)^2,x=seq(0,sqrt(0.0003),length.out=1000)^2)%>%
  mutate(w=log10(y*(1-x)/(x*(1-y))))%>%filter(!is.na(w),x<0.03,w<5)->theory

expand.grid(causal_effect=seq(0,5,1),freq=seq(0,sqrt(0.0003),length.out=1000)^2)%>%
  mutate(
    freq_in_cancer=freq*10^(causal_effect)/(freq*10^(causal_effect)+1-freq))%>%filter(freq_in_cancer<=0.07)->
  theory_lines

ggplot(
  filter(AMLCHf),aes(x=CH_freq,y=AML_freq))+
  geom_raster(data=filter(theory,w<=5),aes(x=x,y=y,fill=w),show.legend=F)+#scale_fill_gradient2(low="skyblue",high="salmon2",mid="white",breaks=seq(-2,3,1),limits=c(-2,3))+
  #scale_fill_viridis(end=0.85,option="H",breaks=seq(-2,3,1),limits=c(-3,3))+
  scale_fill_gradient2(low=viridis(n=1,begin=0.2,end=0.2,option="G",alpha=0.3),high=viridis(n=1,begin=0.8,end=0.8,option="G",alpha=0.6),mid="white")+
  geom_point()+
  geom_line(data=filter(theory_lines),aes(x=freq,y=freq_in_cancer,group=factor(causal_effect)),alpha=0.4,linetype="dotted")+
  geom_errorbar(data=filter(AMLCHf,annotation1!=""),alpha=0.1,width=0,aes(ymin=AML_freq_l,ymax=AML_freq_h))+
  geom_errorbarh(data=filter(AMLCHf,annotation1!=""),alpha=0.1,height=0,aes(xmin=CH_freq_l,xmax=CH_freq_h))+
  geom_text_repel(aes(label=annotation1))+
  scale_x_continuous(trans="sqrt",breaks=c(0,0.00001,0.0001,0.0004,0.001))+scale_y_continuous(trans="sqrt",breaks=c(0,0.01,0.04,0.1))+
  theme_classic()+
  labs(caption="",
    x="frequency in normal blood",y="frequency in AML")
