library(tidyverse)
library(viridis)
library(ggrepel)
library(metR)

hazards<-expand.grid(w=10^seq(-1,3,0.01),x=c(10^seq(-3,0,0.5)))%>%mutate(R=1-x+x*w,y=x*w/(1-x+x*w))
ggplot(hazards,aes(x=w,y=R,colour=x,group=factor(x)))+geom_line()+scale_y_continuous(trans="log10")+scale_x_continuous(trans="log10")+
  scale_color_gradient(trans="log10",low=viridis(n=1,begin=0.8,end=0.8,option="G",alpha=0.3),high="black")+
  theme_classic()+labs(x="carcinogenic effect, w",y="hazard ratio for cancer in tissue given mutation",colour="mutation\nfrequency\nin tissue")


expand.grid(x=seq(0.001,0.999,0.002),y=seq(0.001,0.999,0.002))%>%
  mutate(w=log10(y*(1-x)/(x*(1-y))))->theory
expand.grid(causal_effect=seq(-2,2,1),freq=c(seq(0,10^(-5),10^(-6)),10^seq(-5,0,0.001)))%>%
  mutate(rate_fold_increase=freq*10^causal_effect+1-freq,
         freq_in_cancer=freq*10^causal_effect/(freq*10^causal_effect+1-freq))->
  theory_lines
ggplot(filter(theory,abs(w)<=2))+geom_raster(aes(x=x,y=y,fill=w))+
  scale_fill_gradient2(low=viridis(n=1,begin=0.2,end=0.2,option="G",alpha=0.3),high=viridis(n=1,begin=0.8,end=0.8,option="G",alpha=0.6),mid="white",limits=c(-2,2),breaks=seq(-3,3,1),labels=10^seq(-3,3,1))+
  theme_classic()+#geom_point(x=1/10,y=2/3,size=2)+
  geom_line(data=filter(theory_lines,freq<=0.999&freq>=0.001&freq_in_cancer>=0.001&freq_in_cancer<=0.999),aes(x=freq,y=freq_in_cancer,group=factor(causal_effect)),alpha=0.3,linetype="dotted")+
  labs(fill="",
       y="mutation frequency among cancers, y",x="mutation frequency in normal tissue, x")+theme(legend.position = "bottom", legend.key.width = unit(1.2, "cm"))



my_transx <- scales::trans_new('custom',
                               transform = function(x) {
                                 log(10^(-3)+x)
                               },
                               inverse = function(x) {
                                 exp(x)-10^(-3)
                               })
my_transy <- scales::trans_new('custom',
                               transform = function(x) {
                                 log(10^(-2)+x)
                               },
                               inverse = function(x) {
                                 exp(x)-10^(-2)
                               })

expand.grid(y=my_transy$inverse(seq(-5,0.01,0.01)),x=my_transx$inverse(seq(-7,0.01,0.01)))%>%
  mutate(w=log10(y*(1-x)/(x*(1-y))))%>%filter(!is.na(w))->theory

expand.grid(causal_effect=seq(-2,3,1),freq=my_transx$inverse(seq(-7,0,0.01)))%>%
  mutate(
    freq_in_cancer=freq*10^(causal_effect)/(freq*10^(causal_effect)+1-freq))->
  theory_lines
expand.grid(causal_effect=seq(-2,3,1),freq_in_cancer=my_transx$inverse(seq(-7,0,0.01)))%>%
  mutate(
    freq=freq_in_cancer/(10^(causal_effect)*(1-freq_in_cancer)+freq_in_cancer))->
  theory_lines2
rbind(theory_lines,theory_lines2)->theory_lines3
theory_lines3<-theory_lines3[!duplicated(paste(theory_lines3$causal_effect,theory_lines3$freq,theory_lines3$freq_in_cancer)),]


YF3<-read.csv("ESSC gene summary.csv")

YF3$frac_C_samples_jitter<-YF3$frac_C_samples*exp(rnorm(n=nrow(YF3),mean=0,sd=0.03))

ggplot(
  filter(YF3),aes(x=cellfrac_N,y=frac_C_samples_jitter))+
  geom_raster(data=filter(theory,w<=3&w>=-2),aes(x=x,y=y,fill=w),show.legend=F)+#scale_fill_gradient2(low="skyblue",high="salmon2",mid="white",breaks=seq(-2,3,1),limits=c(-2,3))+
  scale_fill_gradient2(low=viridis(n=1,begin=0.2,end=0.2,option="G",alpha=0.3),high=viridis(n=1,begin=0.8,end=0.8,option="G",alpha=0.6),mid="white")+
  geom_point(size=2,aes(alpha=m_C+m_N>=5),show.legend=F)+
  geom_line(data=filter(theory_lines3),aes(x=freq,y=freq_in_cancer,group=factor(causal_effect)),alpha=0.4,linetype="dotted")+
  geom_errorbar(alpha=0.1,width=0,aes(ymin=frac_C_samples_l,ymax=frac_C_samples_h))+
  geom_errorbarh(alpha=0.1,height=0,aes(xmin=cellfrac_N_l,xmax=cellfrac_N_h))+
  geom_text_repel(aes(label=Gene,alpha=m_C+m_N>=5),show.legend=F)+
  scale_alpha_manual(values=c(0.3,1))+
  scale_y_continuous(trans=my_transy,breaks=c(0,0.01,0.1,1),labels=c(0,0.01,0.1,1),limits=c(0,1))+
  scale_x_continuous(trans=my_transx,breaks=c(0,0.01,0.1,1),labels=c(0,0.01,0.1,1),limits=c(0,1))+
  theme_classic()+
  labs(x="mutation frequency in normal esophageal epithelium",y="mutation frequency among ESSCs",fill="estimated\ncausal\neffect")






my_transx <- scales::trans_new('custom',
                               transform = function(x) {
                                 log(10^(-4.5)+x)
                               },
                               inverse = function(x) {
                                 exp(x)-10^(-4.5)
                               })
my_transy <- scales::trans_new('custom',
                               transform = function(x) {
                                 log(10^(-2)+x)
                               },
                               inverse = function(x) {
                                 exp(x)-10^(-2)
                               })

expand.grid(y=my_transy$inverse(seq(-5,0.01,0.01)),x=my_transx$inverse(seq(-12,0.01,0.01)))%>%
  mutate(w=log10(y*(1-x)/(x*(1-y))))%>%filter(!is.na(w))->theory

expand.grid(causal_effect=seq(0,4,1),freq=my_transx$inverse(seq(-7,0,0.01)))%>%
  mutate(
    freq_in_cancer=freq*10^(causal_effect)/(freq*10^(causal_effect)+1-freq))->
  theory_lines
expand.grid(causal_effect=seq(0,4,1),freq_in_cancer=my_transx$inverse(seq(-12,0,0.01)))%>%
  mutate(
    freq=freq_in_cancer/(10^(causal_effect)*(1-freq_in_cancer)+freq_in_cancer))->
  theory_lines2
rbind(theory_lines,theory_lines2)->theory_lines3
theory_lines3<-theory_lines3[!duplicated(paste(theory_lines3$causal_effect,theory_lines3$freq,theory_lines3$freq_in_cancer)),]


AMLCH_overview <- read.csv("AML normal blood statistics.csv")



ggplot(
  AMLCH_overview,aes(x=cellfrac_N,y=AMLo_f))+
  geom_raster(data=filter(theory,w<=5&w>=-1),aes(
    #alpha=(w<4)+(w>=4)*(5-w)^2,
    x=x,y=y,fill=w),show.legend=F)+#scale_fill_gradient2(low="skyblue",high="salmon2",mid="white",breaks=seq(-2,5,1),limits=c(-2,5))+
  #scale_fill_viridis(end=0.95,option="H",breaks=seq(-2,5,1),limits=c(-5,5))+
  scale_fill_gradient2(low=viridis(n=1,begin=0.2,end=0.2,option="G",alpha=0.3),high=viridis(n=1,begin=0.8,end=0.8,option="G",alpha=0.6),mid="white",midpoint=0)+
  geom_point(size=2)+
  geom_line(data=filter(theory_lines3),aes(x=freq,y=freq_in_cancer,group=factor(causal_effect)),alpha=0.4,linetype="dotted")+
  geom_errorbar(alpha=0.1,width=0,aes(ymin=AMLo_f_l,ymax=AMLo_f_h))+
  geom_errorbarh(alpha=0.1,height=0,aes(xmin=cellfrac_N_l,xmax=cellfrac_N_h))+
  guides(colour = guide_legend(reverse=T))+
  geom_text_repel(aes(label=Gene))+
  scale_y_continuous(trans=my_transy,breaks=c(0,0.01,0.1,1),labels=c(0,0.01,0.1,1),limits=c(0.01,0.5))+
  scale_x_continuous(trans=my_transx,breaks=c(0,0.0001,0.001,0.01,0.1,1),labels=c(0,0.0001,0.001,0.01,0.1,1),limits=c(0,0.03))+
  theme_classic()+
  labs(
    x="mutation frequency in normal blood",y="mutation frequency among AMLs",fill="estimated\ncausal\neffect")




my_transy <- scales::trans_new('custom',
                               transform = function(x) {
                                 log(10^(-3)+x)
                               },
                               inverse = function(x) {
                                 exp(x)-10^(-3)
                               })
normalisation<-mean(AMLCH_overview$AMLy_f/(AMLCH_overview$L*(1-AMLCH_overview$AMLy_f)))/mean(AMLCH_overview$ce_y)

expand.grid(y=my_transy$inverse(seq(-7,0.01,0.01)),x=10^seq(2,4,0.01))%>%
  mutate(w=log10(y/(x*(1-y)*normalisation)))%>%filter(!is.na(w))->theory

expand.grid(w=seq(0,4,1),L=10^seq(1,4,0.01))%>%
  mutate(
    freq_in_cancer=normalisation*L*10^w/(1+normalisation*L*10^w))->
  theory_lines

ggplot(
  AMLCH_overview,aes(x=L,y=AMLy_f))+
  geom_raster(data=filter(theory,w<=4&w>=-1),aes(x=x,y=y,fill=w),show.legend=F)+
  scale_fill_gradient2(low=viridis(n=1,begin=0.2,end=0.2,option="G",alpha=0.3),high=viridis(n=1,begin=0.8,end=0.8,option="G",alpha=0.6),mid="white",midpoint=0)+
  geom_point(size=2)+
  geom_line(data=filter(theory_lines),aes(x=L,y=freq_in_cancer,group=factor(w)),alpha=0.4,linetype="dotted")+
  geom_errorbar(alpha=0.1,width=0,aes(ymin=AMLy_f_l,ymax=AMLy_f_h))+
  guides(colour = guide_legend(reverse=T))+
  geom_text_repel(aes(label=Gene))+
  scale_y_continuous(trans=my_transy,breaks=c(0,0.001,0.01,0.1,1),labels=c(0,0.001,0.01,0.1,1),limits=c(0.001,0.1))+
  scale_x_continuous(trans="log10",limits=c(250,7000)#,breaks=c(0,0.0001,0.001,0.01,0.1,1),labels=c(0,0.0001,0.001,0.01,0.1,1)
                     )+
  theme_classic()+
  labs(
    x="gene length",y="mutation frequency in young AML group",fill="estimated\ncausal\neffect")

