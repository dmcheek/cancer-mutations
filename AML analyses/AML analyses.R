library(tidyverse)
library(ggrepel)
library(viridis)

AMLCH_overview <- read.csv("AML normal blood statistics.csv")

ggplot(AMLCH_overview,aes(x=reorder(Gene,AML_mean_age),y=AML_mean_age,ymin=AML_mean_age_l,ymax=AML_mean_age_h,colour=AML_f))+geom_point()+geom_errorbar(width=0)+theme_minimal()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+labs(x="",y="mean age in AML",colour="frequency\nin AML")+scale_colour_gradient(high="black",low=viridis(n=1,begin=0.8,end=0.8,option="mako",alpha=0.5),limits=c(0,NA))

transx <- scales::trans_new('custom',
                               transform = function(x) {
                                 log(10^(-3)+x)
                               },
                               inverse = function(x) {
                                 exp(x)-10^(-3)
                               })
transy <- scales::trans_new('custom',
                               transform = function(x) {
                                 log(10^(-1)+x)
                               },
                               inverse = function(x) {
                                 exp(x)-10^(-1)
                               })
log10p <- scales::trans_new('custom',
                            transform = function(x) {
                              log(10+x)
                            },
                            inverse = function(x) {
                              exp(x)-10
                            })
p<-cor.test(AMLCH_overview$AMLy_f,AMLCH_overview$AMLo_f,method="pearson")
sp<-cor.test(AMLCH_overview$AMLy_f,AMLCH_overview$AMLo_f,method="spearman")
ggplot(AMLCH_overview,aes(x=AMLy_f,y=AMLo_f,label=Gene))+
 geom_point()+geom_text_repel()+scale_x_continuous(trans="log10",breaks=c(0,0.001,0.003,0.01,0.03,0.1,0.3),labels=c(0,0.001,0.003,0.01,0.03,0.1,0.3))+scale_y_continuous(trans="log10",breaks=c(0,0.01,0.03,0.1,0.3),labels=c(0,0.01,0.03,0.1,0.3))+labs(x="mutation frequency in young-onset AML",y="mutation frequency in adult-onset AML",caption=paste("rho=",signif(unname(sp$estimate),2),", p=",signif(sp$p.value,2),#,"\nr=",signif(unname(p$estimate),2),", p=",signif(p$p.value,2)
                                                                                                                                                                                                                                                                                                                                                                                                           sep=""))+theme_classic()


transx <- scales::trans_new('custom',
                            transform = function(x) {
                              log(10^(-3.5)+x)
                            },
                            inverse = function(x) {
                              exp(x)-10^(-3.5)
                            })
transy <- scales::trans_new('custom',
                            transform = function(x) {
                              log(10^(-1.5)+x)
                            },
                            inverse = function(x) {
                              exp(x)-10^(-1.5)
                            })

p<-cor.test(AMLCH_overview$ce_o,AMLCH_overview$ce_y,method="pearson")
sp<-cor.test(AMLCH_overview$ce_o,AMLCH_overview$ce_y,method="spearman")
ggplot(AMLCH_overview,aes(x=ce_y,y=pmin(ce_o,30000),label=Gene))+geom_smooth(data=filter(AMLCH_overview,Gene!="NPM1"),method="lm",se=F,colour=viridis(n=1,begin=0.8,end=0.8,option="mako",alpha=0.5))+geom_point(aes(shape=Gene=="NPM1"),show.legend=F)+geom_text_repel()+scale_x_continuous(trans="log10",breaks=c(0,10,100,1000),labels=c(0,10,100,1000))+scale_y_continuous(trans="log10",breaks=c(0,1,10,100,1000,10000))+
  labs(x="estimated relative carcinogenic effect for young-onset AML",y="estimated carcinogenic effect for adult-onset AML",caption=paste("rho=",signif(unname(sp$estimate),2),", p=",signif(sp$p.value,2),#"\nr=",signif(unname(p$estimate),2),", p=",signif(p$p.value,2),
                                                                                                                      sep=""))+theme_classic()


spagef<-cor.test(AMLCH_overview$AML_mean_age,AMLCH_overview$AML_f,method="spearman")
ggplot(AMLCH_overview,aes(x=AML_f,y=AML_mean_age,ymin=AML_mean_age_l,ymax=AML_mean_age_h,label=Gene))+
  geom_point()+geom_text_repel()+#geom_errorbar(width=0,alpha=0.2)+
  theme_classic()+labs(x="mutation frequency in AML",y="mutations' mean age in AML",caption=paste("rho=",signif(spagef$estimate,2),", p=",signif(spagef$p.value,2)),title="AML patients of all ages")


ggplot(AMLCH_overview,aes(x=possel,y=AML_mean_age))+geom_boxplot(outlier.shape=NA)+geom_jitter()+stat_compare_means(comparisons=list(c("FALSE","TRUE")))+theme_classic()+labs(y="mean age in AML",x="dN/dS>1: positive selection in normal blood")


transx <- scales::trans_new('custom',
                            transform = function(x) {
                              log(10^(-4.5)+x)
                            },
                            inverse = function(x) {
                              exp(x)-10^(-4.5)
                            })

AMLCH_overview$signif_sel_normal_blood<-AMLCH_overview$possel
AMLCH_overview$signif_sel_normal_blood_text<-"q<0.05"
AMLCH_overview$signif_sel_normal_blood_text[!AMLCH_overview$signif_sel_normal_blood]<-"N.S."
sp<-cor.test(AMLCH_overview$cellfrac_N/AMLCH_overview$L,AMLCH_overview$AML_mean_age,method="spearman")
ggplot(AMLCH_overview,aes(x=1000*cellfrac_N/L,y=AML_mean_age,label=Gene))+geom_smooth(method="lm",se=F,colour=viridis(n=1,begin=0.8,end=0.8,option="mako",alpha=0.5))+geom_point(aes(shape=signif_sel_normal_blood_text))+scale_shape_manual(values=c(16,8))+geom_text_repel()+scale_x_continuous(trans=transx,breaks=c(0,0.0001,0.001,0.01),labels=c(0,0.0001,0.001,0.01))+
  labs(y="mean age in AML",x="mutations per kilobase per normal blood cell",colour="dN/dS>1\nin normal\nblood:\n-log10(q)",shape="dN/dS>1\nin normal\nblood",caption=paste("rho=",signif(unname(sp$estimate),2),", p=",signif(sp$p.value,2),sep=""))+theme_classic()#+scale_colour_viridis(direction=-1,option="G",breaks=c(0,3,6),limits=c(0,6),labels=c(0,3,">6"))
 

ggplot(AMLCH_overview,aes(x=pmin(ce_o,30000),y=AML_mean_age,label=Gene))+geom_point()+scale_shape_manual(values=c(16,8))+geom_text_repel()+scale_x_continuous(trans="log10")+geom_smooth(method="lm",se=F)+
  labs(y="mean age in AML",x="carcinogenic effect estimate: old-onset AML",colour="dN/dS>1\nin normal\nblood:\n-log10(q)",caption=paste("rho=",signif(unname(sp$estimate),2),", p=",signif(sp$p.value,2),sep=""))+theme_classic()#+scale_colour_viridis(direction=-1,option="G",breaks=c(0,3,6),limits=c(0,6),labels=c(0,3,">6"))
ggplot(AMLCH_overview,aes(x=ce_y,y=AML_mean_age,label=Gene))+geom_point()+scale_shape_manual(values=c(16,8))+geom_text_repel()+scale_x_continuous(trans="log10")+geom_smooth(method="lm",se=F)+
  labs(y="mean age in AML",x="carcinogenic effect estimate: young-onset AML",colour="dN/dS>1\nin normal\nblood:\n-log10(q)",caption=paste("rho=",signif(unname(sp$estimate),2),", p=",signif(sp$p.value,2),sep=""))+theme_classic()#+scale_colour_viridis(direction=-1,option="G",breaks=c(0,3,6),limits=c(0,6),labels=c(0,3,">6"))

spceyage<-cor.test(AMLCH_overview$AML_mean_age,AMLCH_overview$ce_y,method="spearman")
spceoage<-cor.test(AMLCH_overview$AML_mean_age,AMLCH_overview$ce_o,method="spearman")
a2<-rbind(transmute(AMLCH_overview,Gene=Gene,cellfrac_N=cellfrac_N,L=L,ce=ce_o,ce_l=pmin(ce_o,ce_y),ce_h=pmax(ce_o,ce_y),AML_mean_age=AML_mean_age,
                    `carcinogenic effect estimate`=paste("\nadult-onset AML\nrho=",signif(spceoage$estimate,2),", p=",signif(spceoage$p.value,2)),Gene=Gene),transmute(AMLCH_overview,Gene=Gene,ce=ce_y,ce_l=pmin(ce_o,ce_y),ce_h=pmax(ce_o,ce_y),cellfrac_N=cellfrac_N,L=L,AML_mean_age=AML_mean_age,
                    `carcinogenic effect estimate`=paste("\nyoung-onset AML\nrho=",signif(spceyage$estimate,2),", p=",signif(spceyage$p.value,2)),Gene=""))
spceage<-cor.test(a2$AML_mean_age,a2$ce,method="spearman")
ggplot(a2,aes(y=AML_mean_age,x=pmin(ce,30000),xmin=ce_l,xmax=pmin(ce_h,30000),colour=`carcinogenic effect estimate`,label=Gene))+geom_smooth(data=filter(a2,ce<Inf),method="lm",se=F)+geom_errorbarh(colour="gray")+geom_point(aes(shape=Gene=="NPM1"),show.legend=F)+geom_text_repel()+theme_classic()+scale_x_continuous(trans="log10")+labs(x="carcinogenic effect estimate",y="mean age in AML",colour="carcinogenic\neffect estimate")+scale_colour_manual(values=c(viridis(n=1,begin=0.2,end=0.2,option="G",alpha=1),viridis(n=1,begin=0.8,end=0.8,option="G",alpha=1)))


sp<-cor.test(AMLCH_overview$AML_mean_age,rank(AMLCH_overview$ce_o)+rank(AMLCH_overview$ce_y),method="spearman")
ggplot(AMLCH_overview,aes(x=(rank(ce_y)+rank(ce_o))/2,y=AML_mean_age,label=Gene,xmin=pmin(rank(ce_o),rank(ce_y)),xmax=pmax(rank(ce_o),rank(ce_y))))+geom_errorbarh(alpha=0.3)+geom_point()+geom_text_repel()+geom_smooth(method="lm",se=F,colour="black")+
  labs(y="mean age in AML",x="ranking carcinogenic effect estimate",colour="dN/dS>1\nin normal\nblood:\n-log10(q)",caption=paste("rho=",signif(unname(sp$estimate),2),", p=",signif(sp$p.value,2),sep=""))+theme_classic()

rbind(transmute(AMLCH_overview,Gene=Gene,group="combined",ce=(rank(ce_y)+rank(ce_o))/2,AML_mean_age=AML_mean_age,xmin=pmin(rank(ce_y),rank(ce_o)),xmax=pmax(rank(ce_y),rank(ce_o))),
      transmute(AMLCH_overview,Gene="",group="young",ce=rank(ce_y),AML_mean_age=AML_mean_age,xmin=pmin(rank(ce_y),rank(ce_o)),xmax=pmax(rank(ce_y),rank(ce_o))),
      transmute(AMLCH_overview,Gene="",group="old",ce=rank(ce_o),AML_mean_age=AML_mean_age,xmin=pmin(rank(ce_y),rank(ce_o)),xmax=pmax(rank(ce_y),rank(ce_o)))
)->a3
ggplot(a3,aes(x=ce,y=AML_mean_age,label=Gene,colour=group,xmin=xmin,xmax=xmax,shape=group))+scale_colour_manual(values=c(viridis(n=1,begin=0.5,end=0.5,option="mako"),viridis(n=1,begin=0,end=0,option="mako"),viridis(n=1,begin=1,end=1,option="mako")))+geom_smooth(method="lm",se=F,show.legend=F)+geom_errorbarh(height=0,colour="gray")+geom_point()+geom_text_repel(show.legend=F,colour="black")+
  labs(y="mean age in AML",x="rank carcinogenic effect estimate",colour="age group for\ncarcinogenic\neffect",shape="age group for\ncarcinogenic\neffect")+theme_classic()
ggplot(filter(a3,group=="combined"),aes(x=ce,y=AML_mean_age,label=Gene,xmin=xmin,xmax=xmax))+scale_colour_manual(values=c(viridis(n=1,begin=0.5,end=0.5,option="mako"),viridis(n=1,begin=0,end=0,option="mako"),viridis(n=1,begin=1,end=1,option="mako")))+geom_smooth(method="lm",se=F,show.legend=F)+geom_errorbarh(height=0,colour="gray")+geom_point()+geom_text_repel(show.legend=F,colour="black")+
  labs(y="mean age in AML",x="rank carcinogenic effect estimate",colour="age group for\ncarcinogenic\neffect",shape="age group for\ncarcinogenic\neffect")+theme_classic()

AMLcors<-tibble(variable=c("mutation frequency;\nAML","mutation density;\nAML","carcinogenic effect;\nyoung","carcinogenic effect;\nold","carcinogenic effect;\ncombined","mutation frequency;\nnormal blood","mutation density;\nnormal blood"),
                rho=c(
                  cor(AMLCH_overview$AML_mean_age,AMLCH_overview$AML_f,method="spearman"),
                  cor(AMLCH_overview$AML_mean_age,AMLCH_overview$AML_f/AMLCH_overview$L,method="spearman"),
                  cor(AMLCH_overview$AML_mean_age,AMLCH_overview$ce_y,method="spearman"),
                  cor(AMLCH_overview$AML_mean_age,AMLCH_overview$ce_o,method="spearman"),
                  cor(AMLCH_overview$AML_mean_age,rank(AMLCH_overview$ce_o)+rank(AMLCH_overview$ce_y),method="spearman"),
                  cor(AMLCH_overview$AML_mean_age,AMLCH_overview$cellfrac_N,method="spearman"),
                  cor(AMLCH_overview$AML_mean_age,AMLCH_overview$cellfrac_N/AMLCH_overview$L,method="spearman")
                  ),
                p=c(
                  cor.test(AMLCH_overview$AML_mean_age,AMLCH_overview$AML_f,method="spearman")$p.value,
                  cor.test(AMLCH_overview$AML_mean_age,AMLCH_overview$AML_f/AMLCH_overview$L,method="spearman")$p.value,
                  cor.test(AMLCH_overview$AML_mean_age,AMLCH_overview$ce_y,method="spearman")$p.value,
                  cor.test(AMLCH_overview$AML_mean_age,AMLCH_overview$ce_o,method="spearman")$p.value,
                  cor.test(AMLCH_overview$AML_mean_age,rank(AMLCH_overview$ce_o)+rank(AMLCH_overview$ce_y),method="spearman")$p.value,
                  cor.test(AMLCH_overview$AML_mean_age,AMLCH_overview$cellfrac_N,method="spearman")$p.value,
                  cor.test(AMLCH_overview$AML_mean_age,AMLCH_overview$cellfrac_N/AMLCH_overview$L,method="spearman")$p.value
                ))%>%mutate(signif="*")
AMLcors$signif<- -0.03*sign(AMLcors$rho)
AMLcors$signif[AMLcors$p>0.05]<-NA
ggplot(AMLcors,aes(x=reorder(variable,rho),fill=-log10(p)))+geom_point(mapping=aes(y=signif),shape=8)+geom_col(mapping=aes(y=rho),colour="black")+theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+labs(x="",y="correlation with mean age in AML")+scale_fill_gradient(low="white",high=viridis(n=1,begin=0.8,end=0.8,option="mako",alpha=1),limits=c(0,NA))




ggplot(AMLCH_overview,aes(y=reorder(Gene,diff_rank_ce),xmin=diff_rank_ce_l,xmax=diff_rank_ce_h,x=diff_rank_ce))+geom_errorbarh(height=0)+geom_point()+theme_classic()+labs(y="",x="old vs young difference of rank carcinogenic effect")


agedataboxes<-read.csv("AMLdataageboxes.csv")
agedataboxes$group[agedataboxes$group==" normal blood"]<-" normal\nblood"
agedataboxes$group[agedataboxes$group=="young AML"]<-"young\nAML"

ggplot(agedataboxes,aes(x=paste(data," (",n,")",sep=""),y=Age,fill=n))+geom_violin(outlier.shape=NA)+facet_wrap(~group,scales="free_x",nrow=1)+theme_classic()+
  scale_fill_gradient(low="white",high=viridis(n=1,begin=0,end=0,option="viridis"),limits=c(0,NA))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+labs(y="age",fill="number of\nsamples",x="")


cor_bounds<-read.csv("AMLrobustness.csv")
ggplot(cor_bounds,aes(x=young_upper,y=correlation,colour=comparison))+geom_line()+facet_wrap(~old_lower,nrow=1)+theme_classic()+labs(x="upper age limit for young AML group",subtitle="lower age limit for old AML group",colour="old vs young AML")+theme(plot.subtitle = element_text(hjust = 0.5))+scale_colour_manual(values=c(viridis(n=1,begin=0.2,end=0.2,option="mako"),viridis(n=1,begin=0.8,end=0.8,option="mako")))


