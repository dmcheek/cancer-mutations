library(tidyverse)
library(viridis)
library(effectsize)
library(ggrepel)
library(data.table)

samplesf<-read.csv("SCNAs_vs_gene_mutations.csv")

samplesf%>%filter(!any_alt|alt)%>%group_by(TCGA,Gene,TOG,CNA_type)%>%summarise(n_alt=sum(alt),n_WT=sum(!alt),
                                                                               alt_age_bias=
                                                                                 rank_biserial(Age[alt],Age[!alt])[[1]],
                                                                               alt_age_bias_l=
                                                                                 rank_biserial(Age[alt],Age[!alt])[[3]],
                                                                               alt_age_bias_h=
                                                                                 rank_biserial(Age[alt],Age[!alt])[[4]],
                                                                               f_alt=mean(alt))->alt_age0



samplesf%>%
      filter(
        !any_alt|alt
      )%>%filter(mut==alt)%>%group_by(CNA_type,TCGA,Gene,TOG)%>%summarise(n_alt_mut=sum(alt),n_WT=sum(!alt),
                                         alt_mut_age_bias=rank_biserial(Age[alt],Age[!alt])[[1]],
                                         alt_mut_age_bias_l=rank_biserial(Age[alt],Age[!alt])[[3]],
                                         alt_mut_age_bias_h=rank_biserial(Age[alt],Age[!alt])[[4]],
                                         f_alt_mut=mean(alt))->alt_mut_age
    
    
samplesf%>%filter(
      mut==F
    )%>%
      filter(
        !any_alt|alt
      )%>%group_by(CNA_type,TCGA,Gene,TOG)%>%summarise(n_alt_no_mut=sum(alt),n_alt_WT=sum(!alt),
                                         alt_age_bias_no_mut=rank_biserial(Age[alt],Age[!alt])[[1]],
                                         alt_age_bias_no_mut_l=rank_biserial(Age[alt],Age[!alt])[[3]],
                                         alt_age_bias_no_mut_h=rank_biserial(Age[alt],Age[!alt])[[4]],
                                         f_alt_no_mut=mean(alt))->alt_age

samplesf%>%filter(
      any_alt==F
    )%>%group_by(CNA_type,TCGA,Gene,TOG)%>%summarise(n_mut_no_alt=sum(mut),n_mut_WT=sum(!mut),
                                       mut_age_bias_no_alt=rank_biserial(Age[mut],Age[!mut])[[1]],
                                       mut_age_bias_no_alt_l=rank_biserial(Age[mut],Age[!mut])[[3]],
                                       mut_age_bias_no_alt_h=rank_biserial(Age[mut],Age[!mut])[[4]],
                                       f_mut_no_alt=mean(mut)
    )->mut_age
    

samplesf%>%group_by(Gene,TOG,TCGA,CNA_type)%>%summarise(n_mut_no_alt=sum(mut&!any_alt),n_alt_no_mut=sum(alt&!mut))%>%filter(n_mut_no_alt>4&n_alt_no_mut>4)%>%
  left_join(samplesf)%>%group_by(Gene,TOG,TCGA,CNA_type)%>%summarise(
  age_diff=median(Age[alt&!mut])-median(Age[mut&!alt]),p=wilcox.test(Age[alt&!mut],Age[mut&!alt])$p.value)->
  gene_specific



    age_biases<-inner_join(alt_age,mut_age)%>%left_join(alt_mut_age)


age_biases%>%filter(CNA_type=="loss",
              #  pmin(n_mut_no_alt,n_alt_no_mut)>=10
                
              pmin(f_mut_no_alt,f_alt_no_mut)>=0.05
                ,TOG=="tsg"
       
)->tsgloss

sptsg<-cor.test(tsgloss$mut_age_bias_no_alt,tsgloss$alt_age_bias_no_mut,method="spearman")
petsg<-cor.test(tsgloss$mut_age_bias_no_alt,tsgloss$alt_age_bias_no_mut,method="pearson")


age_biases%>%filter(CNA_type=="gain",
               #  pmin(n_mut_no_alt,n_alt_no_mut)>=10
               pmin(f_mut_no_alt,f_alt_no_mut)>=0.05
                
                ,TOG=="oncogene"
)->oncogain
sponco<-cor.test(oncogain$mut_age_bias_no_alt,oncogain$alt_age_bias_no_mut,method="spearman")
peonco<-cor.test(oncogain$mut_age_bias_no_alt,oncogain$alt_age_bias_no_mut,method="pearson")


CNAvsSNVsmallindel<-rbind(mutate(tsgloss,TOG=" tumour suppressor genes"),mutate(oncogain,TOG="oncogenes"))

ggplot(CNAvsSNVsmallindel,aes(y=mut_age_bias_no_alt,x=alt_age_bias_no_mut,label=paste(Gene),size=(f_alt_no_mut+f_mut_no_alt)/2))+
  geom_smooth(data=filter(CNAvsSNVsmallindel,TOG==" tumour suppressor genes"),method="lm",se=F,color=viridis(n=1,begin=0.8,end=0.8,option="mako",alpha=0.5))+facet_wrap(~TOG,scales="free_x")+
  geom_point(aes(shape=TCGA,colour=TCGA))+
  scale_shape_manual(values=1:length(unique(CNAvsSNVsmallindel$TCGA)))+
  scale_colour_viridis(option="viridis",discrete=T)+
  geom_text_repel(aes(colour=TCGA),show.legend=F)+
  theme_classic()+
  labs(caption=paste("rho=",signif(sptsg$estimate,2),", p=",signif(sptsg$p.value,2),"; rho=",signif(sponco$estimate,2),", p=",signif(sponco$p.value,2),sep=""),
         shape="",colour="",y="age bias of SNVs/small indels without CNAs",x="age bias of deletions (left) and amplifications (right) without SNVs/small indels")+guides(size="none")





###linear models
lm_tsg0<-summary(lm(alt_age_bias_no_mut~mut_age_bias_no_alt+TCGA-1,data=tsgloss))[["coefficients"]]
lm_tsg<-mutate(as_tibble(lm_tsg0),coef=rownames(lm_tsg0))%>%mutate(cih=Estimate+1.96*`Std. Error`,cil=Estimate-1.96*`Std. Error`)
lm_tsg$coef[1]<-"age bias of SNVs/small indels"
lm_tsg$coef[-1]<-substr(lm_tsg$coef[-1],start=5,stop=str_length(lm_tsg$coef[-1]))
ggplot(lm_tsg,aes(x=coef,y=Estimate,ymin=cil,ymax=cih))+geom_point()+geom_errorbar(width=0)+theme_classic()+geom_hline(yintercept=0,linetype="dashed")+
  labs(title="linear regression among tumour suppressor genes",y="estimated effect on age bias of deletion",x="")+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

lm_onco0<-summary(lm(alt_age_bias_no_mut~mut_age_bias_no_alt+TCGA-1,data=oncogain))[["coefficients"]]
lm_onco<-mutate(as_tibble(lm_onco0),coef=rownames(lm_onco0))%>%mutate(cih=Estimate+1.96*`Std. Error`,cil=Estimate-1.96*`Std. Error`)
lm_onco$coef[1]<-"age bias of SNVs/small indels"
lm_onco$coef[-1]<-substr(lm_onco$coef[-1],start=5,stop=str_length(lm_onco$coef[-1]))
ggplot(lm_onco,aes(x=coef,y=Estimate,ymin=cil,ymax=cih))+geom_point()+geom_errorbar(width=0)+theme_classic()+geom_hline(yintercept=0,linetype="dashed")+
  labs(title="linear regression among oncogenes",y="estimated effect on age bias of amplification",x="")+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))








#filtering parameters





mbf<-c(250,500,1000,Inf)


rbindlist(
  lapply(1:length(mbf),function(i){
  
samplesff<-filter(samplesf,mb<=mbf[i])

samplesff%>%filter(
  mut==F
)%>%
  filter(
    !any_alt|alt
  )%>%group_by(CNA_type,TCGA,Gene,TOG)%>%summarise(n_alt_no_mut=sum(alt),n_alt_WT=sum(!alt),
                                                   alt_age_bias_no_mut=
                                                     rank_biserial(Age[alt],Age[!alt])[[1]],
                                                   f_alt_no_mut=mean(alt))->alt_age

samplesff%>%filter(
  any_alt==F
)%>%group_by(CNA_type,TCGA,Gene,TOG)%>%summarise(n_mut_no_alt=sum(mut),n_mut_WT=sum(!mut),
                                                 mut_age_bias_no_alt=
                                                   rank_biserial(Age[mut],Age[!mut])[[1]],
                                                 f_mut_no_alt=mean(mut)
)->mut_age



age_biases<-inner_join(alt_age,mut_age)


age_biases%>%filter(CNA_type=="loss",
                  TOG=="tsg"
)->tsgloss



age_biases%>%filter(CNA_type=="gain",
                    TOG=="oncogene"
)->oncogain


rbindlist(lapply(seq(0,0.15,0.01),function(m){
  
  tsglossf<-tsgloss%>%filter(pmin(f_alt_no_mut,f_mut_no_alt)>m)
  oncogainf<-oncogain%>%filter(pmin(f_alt_no_mut,f_mut_no_alt)>m)
  
tibble(mbf=mbf[i],
         m=m,
       
cor_type=rep(c("Spearman","Pearson"),2),
alt_type=c("\noncogene mutation\nvs amplification","\noncogene mutation\nvs amplification","\nTSG mutation\nvs deletion","\nTSG mutation\nvs deletion"),
n_genes=c(nrow(oncogainf),nrow(oncogainf),nrow(tsglossf),nrow(tsglossf)),
         cor=c(cor(oncogainf$mut_age_bias_no_alt,oncogainf$alt_age_bias_no_mut,method="spearman"),
               pearson_onco=cor(oncogainf$mut_age_bias_no_alt,oncogainf$alt_age_bias_no_mut),
               cor(tsglossf$mut_age_bias_no_alt,tsglossf$alt_age_bias_no_mut,method="spearman"),
               cor(tsglossf$mut_age_bias_no_alt,tsglossf$alt_age_bias_no_mut)),
corp=c(cor.test(oncogainf$mut_age_bias_no_alt,oncogainf$alt_age_bias_no_mut,method="spearman")$p.value,
       pearson_onco=cor.test(oncogainf$mut_age_bias_no_alt,oncogainf$alt_age_bias_no_mut)$p.value,
       cor.test(tsglossf$mut_age_bias_no_alt,tsglossf$alt_age_bias_no_mut,method="spearman")$p.value,
       cor.test(tsglossf$mut_age_bias_no_alt,tsglossf$alt_age_bias_no_mut)$p.value)
)
  
}))




}))->corfilt




corfilt$mbf[corfilt$mbf==Inf]<-"no bound"
corfilt$mbf[corfilt$mbf==500]<-" 500"
corfilt$mbf[corfilt$mbf==250]<-" 250"

corfilt$corsig<-"p<0.05"

corfilt$corsig[corfilt$corp>=0.05]<-"N.S."

ggplot(filter(corfilt,cor_type=="Spearman",n_genes>=20),aes(x=m,y=cor,colour=alt_type))+geom_point(aes(shape=corsig))+geom_line()+facet_wrap(~mbf,nrow=1,scales="free_x")+labs(subtitle="upper bound on number of SNVs/small indels per sample",y="age bias correlation",x="lower bound on fraction of samples altered per gene-cancer type pair",colour="",linetype="",shape="")+theme_classic()+scale_x_continuous()+theme(plot.subtitle = element_text(hjust = 0.5))+
  scale_shape_manual(values=c(8,16),breaks=c("p<0.05","N.S."))+scale_colour_manual(values=c(viridis(n=1,begin=0.2,end=0.2,option="G"),viridis(n=1,begin=0.8,end=0.8,option="G")),
                                                                                   breaks=c("\nTSG mutation\nvs deletion","\noncogene mutation\nvs amplification"))



corfilt$alt_type1<-"TSGs"
corfilt$alt_type1[grepl("oncogene",corfilt$alt_type)]<-"oncogenes"
ggplot(corfilt,aes(x=m,y=n_genes,colour=alt_type1))+geom_line()+facet_wrap(~mbf,nrow=1)+labs(y="number of gene-cancer type pairs",x="lower bound on fraction of samples altered per gene-cancer type pair",colour="",linetype="",subtitle="upper bound on number of SNVs/small indels per sample")+theme_minimal()+
  scale_x_continuous(limits=c(0,0.1),breaks=seq(0,0.1,0.02))+scale_y_continuous(limits=c(0,NA))+theme(plot.subtitle = element_text(hjust = 0.5))+
  scale_colour_manual(values=c(viridis(n=1,begin=0.2,end=0.2,option="G"),viridis(n=1,begin=0.8,end=0.8,option="G")),
                      breaks=c("TSGs","oncogenes"))


