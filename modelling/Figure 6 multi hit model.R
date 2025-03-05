library(data.table)
library(tidyverse)
library(expm)
library(ggridges)
library(viridis)

###parameter values
v<-0.003    #mutation rate
r0<-10^(-9) #baseline cancer rate
m<-10       #count up to m mutations
alpha<-2    #each mutation's log carcinogenic effect is sampled from a gamma distribution with shape and scale parameters alpha and beta
beta<-1     
lifespan<-80 #look at lineage for 80 years

####begin with Figure S6A
runs<-5000   #sample carcinogenic effects for m*runs mutations
logw<-rgamma(n=m*runs,shape=alpha,rate=beta) #log carcinogenic effects of mutations

rbindlist(lapply(1:runs,function(l){
  
  logw<-logw[((l-1)*m+1):(l*m)]
  
 
  Q<-matrix(data=0,nrow=m+1,ncol=m+1)
  for (i in 1:m){
    Q[i,i]<--v-r0*exp(c(0,cumsum(logw)))[i]
    Q[i,i+1]<-v
  }
  Q[m+1,m+1]<--v-r0*exp(c(0,cumsum(logw)))[m+1]
  
  rbindlist(lapply(seq(20,80,10),function(t){
    
    tibble(run=l,t=t,k=0:m,log_w=c(0,logw), cumw=r0*exp(cumsum(log_w)),fn=c(c(1,rep(0,m))%*%expm(Q*t)))
  }))
  
}
))->mutations


mutations%>%group_by(t,k)%>%summarise(fn=mean(fn),cr=mean(fn*cumw))->crk
crk%>%group_by(t)%>%reframe(k=k,fn=fn/sum(fn),cr=cr/sum(cr))->k
k1<-k%>%group_by(k)%>%summarise(maxcr=max(cr))%>%filter(maxcr>0.001)%>%left_join(k)
k1<-k
colnames(k1)[is.element(colnames(k1),c("fn","cr"))]<-c(" normal tissue","cancer")
k2<-k1%>%pivot_longer(cols=c(" normal tissue","cancer"),names_to="nc",values_to="frequency")
ggplot(k2,aes(x=t,y=frequency,fill=factor(k)))+facet_wrap(~nc,nrow=2)+geom_col(position="stack")+theme_classic()+scale_fill_viridis(direction=-1,option="G",discrete=T)+labs(x="age",y="simulated frequency",fill="number of\ncarcinogenic\nmutations")

#next for figures S6BCD

m<-20 #count up to m mutations

par<-expand.grid(v=c(0.001,0.003,0.01))
             
rbindlist(lapply(1:2000,function(l){
  
  logw<-rgamma(n=m,shape=alpha,rate=beta)
  rbindlist(lapply(1:nrow(par),function(b){
    
     v<-par$v[b]
    
    Q<-matrix(data=0,nrow=m+1,ncol=m+1)
    for (i in 1:m){
      Q[i,i]<--v-r0*exp(c(0,cumsum(logw)))[i]
      Q[i,i+1]<-v
    }
    Q[m+1,m+1]<--v-r0*exp(c(0,cumsum(logw)))[m+1]
    
    rbindlist(lapply(seq(10,80,10),function(t){
      
      tibble(run=par$run[l],v=v,t=t,k=0:m,cumlogw=cumsum(c(0,logw)),fn=c(c(1,rep(0,m))%*%expm(Q*t)))
    }))
    
    
    
  }))

  
  
  
}
))->mutations

mutations$fn[mutations$fn>1]<-1 #ridding some numerical instability
mutations1<-mutations%>%group_by(v,t,k)%>%summarise(fn=mean(fn),cr=mean(1-exp(-r0*exp(cumlogw)*fn)))

mutations1%>%summarise(mean_k=weighted.mean(k,cr))->mean_k
ggplot(mean_k,aes(x=t,y=mean_k,colour=factor(v)))+geom_line()+theme_classic()+scale_colour_viridis(breaks=sort(unique(mutations$v),decreasing=T),direction=-1,option="G",discrete=T)+labs(x="age",y="expected number of mutations in cancer",colour="mutation\nrate")+theme(plot.subtitle = element_text(hjust = 0.5))

inc<-mutations1%>%group_by(v,t)%>%summarise(incidence=sum(cr))
ggplot(inc,aes(x=t,y=incidence,colour=factor(v)))+geom_line()+theme_classic()+scale_colour_viridis(breaks=sort(unique(inc$v),decreasing=T),direction=-1,option="G",discrete=T)+labs(x="age",y="cancer rate",colour="mutation\nrate")+scale_y_continuous(trans="log10")+theme(plot.subtitle = element_text(hjust = 0.5))


meanlogw<-mutations%>%filter(k>0)%>%group_by(v,t)%>%summarise(empiricalmeanlogw=weighted.mean(cumlogw/k,1-exp(-r0*exp(cumlogw)*fn)))
ggplot(meanlogw,aes(x=t,y=empiricalmeanlogw,colour=factor(v)))+geom_line()+theme_classic()+scale_colour_viridis(breaks=sort(unique(inc$v),decreasing=T),direction=-1,option="G",discrete=T)+labs(x="age",y="mean log carcinogenic effect of mutations in cancers",colour="mutation\nrate")+theme(plot.subtitle = element_text(hjust = 0.5))



#Figures 6B and C
runs<-20
m<-10
logws<-matrix(rgamma(n=m*runs,shape=alpha,rate=beta),nrow=runs)
#######
par<-expand.grid(logw=seq(0,20,0.5),k=1:m)



rbindlist(lapply(1:nrow(par),function(i){
  
  
  logw<-par$logw[i]
  k<-par$k[i]
  
  logwseq<-expand.grid(posi=1:k,logwseq=1:runs)
  
  rbindlist(lapply(1:nrow(logwseq),function(l){
    
    posi<-logwseq$posi[l]
    logws2<-c(logws[logwseq$logwseq[l],1:k])
    logws2[posi]<-logw
    
    
    cr<-r0*exp(c(0,cumsum(logws2)))
    Q<-matrix(data=0,nrow=k+1,ncol=k+1)
    for (j in 1:k){
      Q[j,j]<--v-cr[j]
      Q[j,j+1]<-v
    }
    Q[k+1,k+1]<--v-cr[k+1]
    
    tibble(run=l,logw=logw,k=k,cr=sapply(seq(20,80,20),function(t0){(c(c(1,rep(0,k[1]))%*%expm(Q*t0))*cr)[k+1]}),t=seq(20,80,20))
    
    
  }))%>%group_by(k,t,logw)%>%summarise(cr=mean(cr))
}))%>%mutate(cr=cr*dgamma(x=logw,shape=alpha,rate=beta))->ridgeplots

mmax<-4
ridgeplots%>%group_by(k,logw)%>%summarise(cr=mean(cr))->ridgeplotsk
ridgeplotsk<-ridgeplotsk%>%group_by(k)%>%reframe(logw=logw,dens=cr/sum(cr),mean_logw=weighted.mean(logw,cr))%>%filter(k<mmax+1)
ridgeplotst<-ridgeplots%>%group_by(t,logw)%>%summarise(cr=mean(cr))%>%group_by(t)%>%reframe(logw=logw,dens=cr/sum(cr),mean_logw=weighted.mean(logw,cr))

ggplot(ridgeplotsk)+geom_ridgeline(aes(y=k,group=factor(k),height=20*dens,fill=mean_logw,x=logw))+coord_flip()+theme_classic()+scale_fill_viridis(option="G",begin=0.3,end=1,limits=c(min(ridgeplotst$mean_logw,ridgeplotsk$mean_logw),max(ridgeplotst$mean_logw,ridgeplotsk$mean_logw)))+labs(x="log carcinogenic effect of mutations",fill="mean log\ncarcinogenic\neffect",y="number of carcinogenic\nmutations per cancer")+scale_y_continuous(breaks=seq(1.3,mmax+0.3,1),labels=1:mmax)#+theme(legend.position="bottom")

ggplot(filter(ridgeplotst))+geom_ridgeline(aes(y=t,group=factor(t),height=500*dens,fill=mean_logw,x=logw))+coord_flip()+theme_classic()+scale_y_continuous(breaks=seq(20,80,20)+4,labels=seq(20,80,20))+scale_fill_viridis(option="G",begin=0.3,end=1,limits=c(min(ridgeplotst$mean_logw,ridgeplotsk$mean_logw),max(ridgeplotst$mean_logw,ridgeplotsk$mean_logw)))+labs(x="log carcinogenic effect of mutations",fill="mean log\ncarcinogenic\neffect",y="cancer initiation age\n")




          






