library(expm)
library(tidyverse)
library(viridis)
library(data.table)
library(spatstat)
library(ggforce)
setwd("~/HMS Dropbox/Naxerova_Lab/Dave/Causation/to publish/figures/figure 4")




  expand.grid(t=0:50,s=seq(0,0.1,0.005))%>%mutate(f=exp(s*t))->expocloneplot

ggplot(expocloneplot,aes(x=t,y=f,colour=s,group=factor(s)))+geom_line()+scale_y_continuous(trans="log10")+theme_classic()+labs(x="years after mutation arrival",y="expected size of mutant\nclone in normal tissue",colour="selection in\nnormal tissue")+ scale_colour_viridis(option="mako")
ggsave("growthcurves.pdf",width=3.5,height=2.8)

Ages<-tibble(t=0:80)

  expand.grid(t=0:80,s=seq(0,0.1,0.005))%>%mutate(f=exp(s*t))%>%left_join(Ages)->cloneplot

cloneplot%>%group_by(s)%>%mutate(cumsumf=cumsum(f))%>%summarise(meanf=mean(f),mean_age=weighted.mean(t,cumsumf*t^3))->cloneplot2
cloneplot2<-cloneplot2%>%mutate(rel_mean_age=mean_age-cloneplot2$mean_age[cloneplot2$s==0])
ggplot(filter(mutate(cloneplot2,ymin=-2,ymax=0,xmin=0,xmax=1000)),aes(x=meanf,y=rel_mean_age,colour=s,group=factor(s)))+geom_point(show.legend=F)+scale_x_continuous(trans="log10",limits=c(1,NA))+geom_hline(yintercept=0,linetype="dotted")+theme_classic()+
  labs(title="null model:\ncarcinogenic effect=1",x="frequency of mutation z normalized\nby neutral mutations in cancers",y="mean age for mutation z relative\nto neutral mutations in cancers",colour="causal\neffect = 1\n\nselection in\nnormal tissue"
  )+scale_y_continuous()+scale_colour_viridis(option="G")
ggsave("null_age_bias.pdf",width=3,height=3)


ggplot(tibble())+geom_line()+scale_x_continuous(trans="log10",limits=c(1,300))+geom_hline(yintercept=0,linetype="dotted")+theme_classic()+
  labs(x="frequency of mutation z normalised\nby neutral mutations in cancers",y="mean age for mutation z relative\nto neutral mutations in cancers",colour="growth curve\nin normal tissue"
  )+scale_y_continuous(limits=c(-5,5))
ggsave("null_age_bias_blank.pdf",width=3.2,height=3)



####figure S4AB

######S4A

Ages<-tibble(t=0:50)                 
a<-10^3                             
rbind(
  expand.grid(model="exponential model",t=0:80,s=seq(0,0.2,0.005))%>%mutate(f=exp(s*t)),
  expand.grid(model="logistic",t=0:80,s=seq(0,0.2,0.005))%>%mutate(f=a/(1+(a-1)*exp(-s*t*a/(a-1)))),  #####N*exp(-b*exp(-c*t)), where exp(-b)=1/N,b=log(N)
  expand.grid(model="Gompertz",t=0:80,s=seq(0,0.2,0.005))%>%mutate(f=exp(log(a)*(1-exp(-s*t/log(a))))),
  expand.grid(model="neutral",t=0:80,s=seq(0,0.2,0.005))%>%mutate(f=1),
  expand.grid(model="polynomial",t=0:80,s=seq(0,0.2,0.005))%>%mutate(f=(s*t/3+1)^3)
)%>%left_join(Ages)->cloneplot
ggplot(filter(cloneplot,model!="neutral",t<=50),aes(x=t,y=f,colour=s,group=factor(s)))+facet_wrap(~model,nrow=1)+geom_line()+theme_classic()+labs(x="years after mutation arrival",y="expected size of mutant\nclone in normal tissue",colour="initial\ngrowth\nrate"
)+scale_colour_viridis(option="mako")+scale_y_continuous(trans="log10")
ggsave("S4Agrowthcurves.pdf",width=9,height=3.5)

cloneplot%>%group_by(model,s)%>%mutate(cumsumf=cumsum(f))%>%summarise(meanf=mean(f),mean_age=weighted.mean(t,cumsumf*t^3))->cloneplot2
cloneplot2<-cloneplot2%>%mutate(rel_mean_age=mean_age-cloneplot2$mean_age[cloneplot2$model=="neutral"])
ggplot(filter(cloneplot2,model!="neutral"),aes(x=meanf,y=rel_mean_age,colour=s,group=factor(s)))+facet_wrap(~model,nrow=1,scales="free")+geom_point(show.legend=F)+scale_x_continuous(trans="log10")+geom_hline(yintercept=0,linetype="dotted")+theme_classic()+labs(x="frequency of mutation z normalized by neutral mutations in cancers",y="mean age for mutation z relative\nto neutral mutations in cancers",colour="selective\neffect",title="null hypothesis of no carcinogenicity"#,caption="young-age bias: incompatible with null hypothesis"
)+scale_colour_viridis(option="mako")+scale_y_continuous()#+geom_ribbon(data=tibble(x=c(1,1000),ymin=-2,ymax=0),mapping=aes(ymin=ymin,ymax=ymax,x=x,y=x),alpha=0.2,fill="firebrick3")
ggsave("S4Bnull_age_bias.pdf",width=9,height=3)



