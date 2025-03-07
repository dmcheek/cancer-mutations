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




AMLCHff<-filter(AMLCHf,AML_n>0)
AMLCHff$mutation_frequency_rank<-(rank(AMLCHff$AML_freq)+rank(AMLCHff$CH_freq))/2

AMLCHfff<-AMLCHff%>%slice_max(order_by=mutation_frequency_rank,n=20)

cor.test(AMLCHfff$mean_age,AMLCHfff$AML_freq,method="spearman")

wtd.cor(rank(AMLCHfff$CH_freq),rank(AMLCHfff$mean_age),weight=AMLCHfff$AML_freq)
spCH<-cor.test(AMLCHfff$CH_freq,AMLCHfff$mean_age,method="spearman")
ggplot(filter(AMLCHfff),aes(x=CH_freq,y=mean_age,label=annotation))+geom_smooth(method="lm",se=F,colour=viridis(n=1,begin=0.8,end=0.8,option="mako",alpha=0.5))+geom_text_repel()+geom_point()+theme_classic()+guides(size="none")+labs()+labs(x="frequency in normal blood",y="mean age in AML",colour="frequency\nin AML",caption=paste("rho=",unname(round(spCH$estimate,2)),", p=",signif(spCH$p.value,2),sep=""))+scale_colour_viridis(end=0.95,option="F",trans="log10",direction=-1)+scale_x_continuous(trans="sqrt",breaks=c(0,0.0001,0.0004,0.001))#+scale_y_continuous(trans="sqrt",breaks=c(0,0.01,0.04,0.1))

wtd.cor(rank(AMLCHfff$OR),rank(AMLCHfff$mean_age),weight=AMLCHfff$AML_freq)
spOR<-cor.test(AMLCHfff$OR,AMLCHfff$mean_age,method="spearman")
ggplot(AMLCHfff,aes(x=OR,y=mean_age,label=annotation))+geom_smooth(method="lm",se=F,colour=viridis(n=1,begin=0.8,end=0.8,option="mako",alpha=0.5))+geom_text_repel()+geom_point()+scale_shape_manual()+theme_classic()+guides(size="none")+labs()+labs(x="estimated carcinogenic effect",y="mean age in AML",colour="frequency\nin AML",caption=paste("rho=",unname(round(spOR$estimate,2)),", p=",signif(spOR$p.value,2),sep=""))+scale_colour_viridis(end=0.95,option="F",trans="log10",direction=-1)+scale_x_continuous(trans="log10",breaks=10^(1:5))

filtparplotw<-read_csv("AML SCNA filtering parameters.csv")
ggplot(filtparplotw,aes(x=upperCNA,y=topalts,fill=rho))+facet_wrap(~CHOR)+geom_tile()+scale_fill_gradient2(mid="white",high=viridis(n=1,begin=0.8,end=0.8,option="G"),low=viridis(n=1,begin=0.2,end=0.2,option="G"),midpoint=0)+theme_classic()+labs(x="upper bound on number of SCNAs per sample",y="top n most common alterations")

