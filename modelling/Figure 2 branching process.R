setwd("~/HMS Dropbox/Naxerova_Lab/Dave/Causation/to publish/figures/Figure 2")
library(ggh4x)
library(tidyverse)
library(viridis)
library(data.table)




####here we provide code for figures 2 and S2
####these results are for a conceptual illustration only, so don't take the precise numerical values seriously

lifespan<-80 #how long are people living? (years)
N<-10^5 #stem cell pool size
v<-10^(-8) #mutation rate
timesteps<-52 #timesteps per year 
rbindlist(
lapply(seq(-16,-13,0.1),function(l){

r<-(1:lifespan)^3*10^l

par<-expand.grid(s=seq(0,0.15,0.01),w=10^seq(0,5,1),t=1:lifespan) #parameters: s is selective effect, w is carcinogenic effect, t is age

sapply(1:nrow(par),function(i){
  
  a<-(par$s[i]+1)/timesteps #division rate of mutant cells
  b<-1/timesteps            #death rate of mutant cells
  r1<-r*par$w[i]/timesteps#rate of carcinogenesis per mutant cell
  theta<-v*N/timesteps
  x<-par$t[i]    #age
  P<-1           #probability that no cancer with mutation arises between ages x and t, given no mutation in tissue at age x
  Q<-1     #probability that no cancer with mutation arises between ages x and t, given one cell with mutation in tissue at age x
  R<-1           #probability that no cancer without mutation arises between ages x and t
  while (x>0){
    for (y in 1:timesteps){
    P<-(1-theta)*P+theta*Q*P
    Q<-(1-a-b-r1[x])*Q+a*Q^2+b
    R<-R-r[x]*N/timesteps
    }
    x<-x-1
  }
  c(P,R) #probability that no cancer (with,without) mutation arises by age t, given no mutation in tissue at birth
})->PR
PR[1,]->par$no_mutant_cancer_yet
PR[2,]->par$no_wildtype_cancer_yet

par2<-par%>%mutate(no_cancer_yet=no_wildtype_cancer_yet*no_mutant_cancer_yet)
sw<-par2%>%group_by(s,w)%>%reframe(t=1:(lifespan-1),
                                   cancer_rate0=no_cancer_yet[-length(no_cancer_yet)]-no_cancer_yet[-1],
                                   mut_f0=(no_mutant_cancer_yet[-length(no_mutant_cancer_yet)]-no_mutant_cancer_yet[-1])/(no_mutant_cancer_yet[-length(no_mutant_cancer_yet)]-no_mutant_cancer_yet[-1]+no_wildtype_cancer_yet[-length(no_wildtype_cancer_yet)]-no_wildtype_cancer_yet[-1]))
                                   

mutate(sw,l=l)
})
)->swl

sw<-swl%>%group_by(s,w,t)%>%summarise(cancer_rate=mean(cancer_rate0),mut_f=weighted.mean(mut_f0,cancer_rate0),t=t[1]) #average over the different values for l to obtain the population wide cancer rates and mutation frequencies

sw$normal_tissue_mut_f<-sw$w^(-1)*sw$mut_f/(1-sw$mut_f+sw$w^(-1)*sw$mut_f) #calculate the mutation's frequency in normal tissue



my_trans5 <- scales::trans_new('custom',
                               transform = function(x) {
                                 log(10^(-5)+x)
                               },
                               inverse = function(x) {
                                 exp(x)-10^(-5)
                               })
ggplot(data=filter(sw,w==1))+geom_line(mapping=aes(x=t,y=normal_tissue_mut_f,colour=s,group=factor(s)),size=1.5)+
  scale_y_continuous(trans=my_trans5,breaks=c(0,10^seq(-5,0,1)),labels=c(0,10^seq(-5,0,1)),limits=c(0,0.11))+
  scale_x_continuous()+
  labs(linetype="w",
       colour="selective\neffect",
       y="expected fraction of normal tissue cells with mutation z",
       x="age")+
  guides(colour = "colorbar")+
  theme_classic()+
  # theme(legend.position = c(0.1, 0.9))+
  #scale_colour_viridis(option="G",discrete=F,breaks=c(0,0.05,0.1,0.15))+
  scale_color_gradient(high=viridis(n=1,begin=0.8,end=0.8,option="G",alpha=0.3),low="black",breaks=c(0,0.05,0.1,0.15))+
  scale_linetype_discrete(breaks=c(100,1))
ggsave("single rainbow.pdf",width=5,height=4.2)


ggplot(data=filter(sw,is.element(w,c(1,100))))+geom_line(mapping=aes(x=t,y=mut_f,colour=s,group=paste(s,w),linetype=factor(w)),size=1.5)+
  scale_y_continuous(trans=my_trans5,breaks=c(0,10^seq(-5,0,1)),labels=c(0,10^seq(-5,0,1)),limits=c(0,0.11))+
  scale_x_continuous()+
  labs(linetype="carcinogenic\neffect",
       colour="selective\neffect",
       y="expected fraction of cancers with mutation z",
       x="cancer initiation age")+
  guides(colour = "none")+
  theme_classic()+
  # theme(legend.position = c(0.1, 0.9))+
  #scale_colour_viridis(option="G",discrete=F,breaks=c(0,0.05,0.1,0.15))+
  scale_color_gradient(high=viridis(n=1,begin=0.8,end=0.8,option="G",alpha=0.3),low="black",breaks=c(0,0.05,0.1,0.15))+
  scale_linetype_manual(breaks=c(100,1),values=c("dotted","solid"))+
  theme(legend.key.width = unit(1.5, "cm"))

ggsave("double rainbow.pdf",width=5,height=4.2)




mean_freqs<-sw%>%filter(is.element(s,seq(0,0.15,0.03)))%>%group_by(s,w)%>%summarise(mean_freq=weighted.mean(mut_f,cancer_rate)) #mutation frequency averaged across all ages
mean_freqs$dNdS<-(mean_freqs$mean_freq/(1-mean_freqs$mean_freq))/(mean_freqs$mean_freq[mean_freqs$s==0&mean_freqs$w==1]/(1-mean_freqs$mean_freq[mean_freqs$s==0&mean_freqs$w==1]))
mean_freqs$dNdSlabel<-format(signif(mean_freqs$dNdS,2),scientific=T)

ggplot(mean_freqs,aes(x=factor(s),y=factor(w),fill=mean_freq))+geom_tile()+geom_text(aes(label=as.character(signif(mean_freq,2))#,color=(dNdS > 1000)),show.legend=F
))+
  #scale_color_manual(values=c("white","black"))+
  #scale_fill_viridis(option="B",trans="log10")+
  scale_fill_gradient(low="white",high=viridis(n=1,begin=0.8,end=0.8,option="G",alpha=1),trans="log10")+
  theme_classic()+labs(x="selective effect in normal tissue",y="carcinogenic effect",fill="frequency\nin cancers")

ggplot(mean_freqs,aes(x=factor(s),y=factor(w),fill=dNdS))+geom_tile()+geom_text(aes(label=dNdSlabel)#,color=(dNdS > 1000)),show.legend=F
                                                                                    )+
  #scale_color_manual(values=c("white","black"))+
  #scale_fill_viridis(option="B",trans="log10")+
  scale_fill_gradient(low="white",high=viridis(n=1,begin=0.8,end=0.8,option="G",alpha=1),trans="log10")+
  theme_classic()+labs(x="selective effect in normal tissue",y="carcinogenic effect",fill="frequency\nrelative to\nneutral\nmutations\nin cancers")
ggsave("dNdS.pdf",width=5.5,height=4)
mean_ages<-sw%>%filter(is.element(s,seq(0,0.15,0.03)))%>%group_by(s,w)%>%summarise(mean_age=weighted.mean(t,mut_f*cancer_rate))
ggplot(mean_ages,aes(x=factor(s),y=factor(w)))+geom_tile(aes(fill=mean_age))+
  #scale_fill_viridis(option="H",trans="reverse")+
  scale_fill_gradient2(high=viridis(n=1,begin=0.8,end=0.8,option="mako",alpha=1),low=viridis(n=1,begin=0,end=0,option="viridis",alpha=0.7),mid="white",midpoint=mean_ages$mean_age[mean_ages$s==0&mean_ages$w==1])+
  geom_text(aes(label=signif(mean_age,3)
                #,color=(mean_age < 70&mean_age>61)),show.legend=F
                ))+
                scale_color_manual(values=c("white","black"))+
  theme_classic()+labs(x="selective effect in normal tissue",y="carcinogenic effect",fill="mean age\nin cancers")+ guides(fill = guide_colorbar(reverse = T))
ggsave("meanage.pdf",width=5.5,height=4)


#We also need parameter exploration for figure S2
r<-(1:lifespan)^3*10^(-14)
par<-expand.grid(R=10^seq(0,1,1),theta=10^seq(-4,0,4),b=c(1,10),s=seq(0,0.15,0.03),w=10^seq(0,5,1),t=1:lifespan)
#parameters are as above, but with the introduction of v as the mutation rate per cell and we scrap the l to keep matters simple


    sapply(1:nrow(par),function(i){
      
      theta<-par$theta[i]/timesteps
      b<-par$b[i]/timesteps
      a<-b+par$s[i]/timesteps
      r1<-par$R[i]*r*par$w[i]/timesteps
      
      
      x<-par$t[i]    #age
      P<-1           #probability that no cancer with mutation arises between ages x and t, given no mutation in tissue at age x
      Q<-1     #probability that no cancer with mutation arises between ages x and t, given one cell with mutation in tissue at age x
      while (x>0){
        for (y in 1:timesteps){
          P<-(1-theta)*P+theta*Q*P
          Q<-(1-a-b-r1[x])*Q+a*Q^2+b
        }
        x<-x-1
      }
      P #probability that no cancer (with,without) mutation arises by age t, given no mutation in tissue at birth
    })->par$no_mutant_cancer_yet

    sw<-par%>%group_by(s,w,b,R,theta)%>%reframe(t=1:(lifespan-1),r=r[-lifespan],r1=r[-lifespan]*w[1],mutant_cancer_rate=(no_mutant_cancer_yet[-length(no_mutant_cancer_yet)]-no_mutant_cancer_yet[-1])/no_mutant_cancer_yet[-length(no_mutant_cancer_yet)])
    
    
R1<-paste("r_0=",unique(sw$R),sep="")
R1<-c("r_0=t^3*10^-14","r_0=t^3*10^-13")
R1<-factor(R1,R1)
b1<-paste("a_0=b_0=",unique(sw$b),sep="")
b1<-factor(b1,b1)
theta1<-paste("N*v=",unique(sw$theta),sep="")
theta1<-factor(theta1,theta1)
sw1<-sw%>%left_join(tibble(theta=unique(sw$theta),theta1=theta1))%>%left_join(tibble(b=unique(sw$b),b1=b1))%>%left_join(tibble(R=unique(sw$R),R1=R1))

mean_freqs<-sw1%>%group_by(s,w,theta1,b1,R1)%>%summarise(f=sum(mutant_cancer_rate))%>%group_by(theta1,b1,R1)%>%reframe(s=s,w=w,dNdS=f/f[s==0&w==1])

ggplot(mean_freqs,aes(x=factor(s),y=factor(w),fill=dNdS))+geom_tile()+facet_nested(b1~R1+theta1)+
  #scale_fill_viridis(option="G",trans="log10",direction=1)+
  scale_fill_gradient(trans="log10",high=viridis(n=1,begin=0.8,end=0.8,option="mako",alpha=1),low="white")+
  theme_classic()+labs(x="selective effect in normal tissue",y="carcinogenic effect",fill="mutation\nfrequency\nrelative to\nneutral\nmutations\nin cancers")+scale_x_discrete(breaks=c(0,0.15))+scale_y_discrete(breaks=c(1,10^5))
ggsave("dNdSS2.pdf",width=8,height=4)

mean_ages1<-sw1%>%group_by(s,w,b1,theta1,R1)%>%summarise(mean_age=weighted.mean(t,mutant_cancer_rate))%>%group_by(theta1,b1,R1)%>%reframe(s=s,w=w,rel_mean_age=mean_age-mean_age[s==0&w==1])

ggplot(mean_ages1,aes(x=factor(s),y=factor(w),fill=rel_mean_age))+geom_tile()+facet_nested(b1~R1+theta1)+
  #scale_fill_viridis(option="H",trans="reverse",direction=1,breaks=c(40,50,60,70))+
  scale_fill_gradient2(low=viridis(n=1,begin=0.8,end=0.8,option="mako",alpha=1),high=viridis(n=1,begin=0,end=0,option="viridis",alpha=0.7),mid="white",midpoint=0,trans="reverse")+
   theme_classic()+labs(x="selective effect in normal tissue",y="carcinogenic effect",fill="mutation\nmean age\nrelative to\nneutral\nmutations\nin cancers")+scale_x_discrete(breaks=c(0,0.15))+scale_y_discrete(breaks=c(1,10^5))
ggsave("meanageS2.pdf",width=8,height=4)





