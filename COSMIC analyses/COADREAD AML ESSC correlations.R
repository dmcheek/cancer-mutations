library(tidyverse)
library(viridis)

together<-read_csv("COADREAD AML ESSC correlations.csv")

ggplot(together,aes(x=factor(var1),y=factor(var2)))+facet_wrap(~cancer,ncol=3)+scale_size_continuous(range=c(2,20),breaks=c(0.1,0.2,0.3,0.4))+geom_point(aes(fill=factor(sign(rho)),shape=factor(sign(rho)),size=abs(rho)))+theme_minimal()+labs(x="",y="",shape="correlation\nsign",fill="correlation\nsign",size="correlation\nmagnitude")+scale_y_discrete(limits=c("average\ncarcinogenic\neffect","mutation\nburden","driver\nburden"))+geom_text(aes(label=pstar))+scale_shape_manual(values=c(24,25),breaks=c(1,-1))+scale_fill_manual(breaks=c(1,-1),values=c("white",viridis(n=1,alpha=1,begin=0.8,end=0.8,option="mako")))+
  guides(size = guide_legend( 
    override.aes=list(shape = 24)),
    shape = guide_legend(override.aes = list(size = 10)))
