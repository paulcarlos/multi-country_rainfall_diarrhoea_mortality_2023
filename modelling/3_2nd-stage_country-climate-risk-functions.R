library(splines)
library(dlnm)
library(mixmeta)
library(lubridate)
library(tsModel)
library(ggplot2)
library(cowplot)
library(classInt)


# call Wald function
source("D:/R/functions/wald_test.R")

# load data
load("data/mod_out/model-out_1st-stage.rda")

# meta-analysis all locations
mtot1 <- mixmeta(formula=coef1~clim,random=~1|ctry,data=mtab,S=vcov1)
#summary(mtot1)
#summary(mtot1)$i2[1]
#summary(mtot1)$q$pvalue[1]
#fwald(mtot1,"clim")

# blup 
blup1 <- blup(mtot1,vcov=T) #with vcov
names(blup1) <- mtab$id

# remove duplicates from sub-climate
r.blup <- !duplicated(do.call(rbind,lapply(blup1,"[","blup")))
gtab <- mtab[r.blup,]

# create country names ####
ctry <- c("arg"="Argentina","bra"="Brazil","cri"="Costa Rica","ecu"="Ecuador","ind"="India",
          "per"="Peru","phl"="Philippines","tha"="Thailand","zaf"="South Africa")

# select country and climate
loc1 <- substr(gtab$id,1,5)

# dataframe 
coln0 <- c("ctry","clim","cen.val","cen.pct",paste0("pk.",c("rain","pct","mean","ci")),
           paste0("p5.",c("rain","mean","ci")),paste0("p95.",c("rain","mean","ci")))
rr <- data.frame(matrix(NA,nrow=length(loc1),ncol=length(coln0),dimnames=list(1:length(loc1),coln0))) # store RRs for 5th and 95th pctile
coln1 <- c("loc","rain","value","lower","upper")
plt <- data.frame(matrix(0,nrow=0,ncol=length(coln1),dimnames=list(NULL,coln1))) # store RRs
coln2 <- c("loc",paste0("p",sprintf("%02d",c(0,5,25,50,75,95,99))))
rain <- data.frame(matrix(0,nrow=length(loc1),ncol=length(coln2),dimnames=list(seq(loc1),coln2))); rain$loc=loc1
#cen <- data.frame("loc"=loc1,"cen.val"=NA,"cen.pct"=NA,"peak.val"=NA,"peak.pct"=NA)

# loop to populate dataframes
byx <- 10
for (i in seq(loc1)) {
  stab <- mtab[grep(loc1[i],mtab$id),] # locations
  # get precipitation
  if (nrow(stab)>1) {
    slist <- dlist[stab$id]
    day1 <- unlist(lapply(slist,function(z)as.character(z$date)))
    var <- unlist(lapply(slist,"[",paste0("pr",lag1)))
    ag1 <- aggregate(var~day1,FUN="mean")
    var1 <- ag1$var
  } else {
    slist <- dlist[[stab$id]]
    var1 <- slist[,paste0("pr",lag1)]
  }
  # centre
  pred1 <- blup1[grep(loc1[i],names(blup1))]
  pred1 <- pred1[[1]]
  predvar <- quantile(var1,0:90/100,na.rm=T)
  argvar1 <- list(x=predvar,fun="ns",knots=quantile(var1,rkn,na.rm=T),Bound=range(var1,na.rm=T))
  bvar <- do.call(onebasis,argvar1)
  cval <- round(predvar[which.min(bvar%*%pred1$blup)])
  
  # model
  argvar <- list(x=var1,fun="ns",knots=quantile(var1,rkn,na.rm=T))
  bvar <- do.call(onebasis,argvar)
  cp <- crosspred(bvar,coef=pred1$blup,vcov=pred1$vcov,model.link="log",by=byx,cen=cval)
  
  # dataframe populate
  df1 <- data.frame("loc"=loc1[i],"rain"=cp$predvar,"value"=cp$allRRfit,"lower"=cp$allRRlow,"upper"=cp$allRRhigh)
  plt <- rbind(plt,df1)
  rain[i,-1] <- round(quantile(var1,as.numeric(substr(colnames(rain)[-1],2,3))/100),-1)
  rr$ctry[i] <- ctry[substr(loc1[i],1,3)]; rr$clim[i] <- gtab$clim2[i]
  
  # center value
  rr$cen.val[i] <- cval; rr$cen.pct[i] <- round(ecdf(var1)(cval)*100) # centering value
  
  # peak % change
  pkval <- cp$predvar[which.max(cp$allRRfit)]
  rr$pk.rain[i] <- pkval; rr$pk.pct[i] <- round(ecdf(var1)(pkval)*100)
  rr$pk.mean[i] <- round((cp$allRRfit[as.character(pkval)]-1)*100,1)
  rr$pk.ci[i] <- paste0(round((cp$allRRlow[as.character(pkval)]-1)*100,1),"; ",round((cp$allRRhigh[as.character(pkval)]-1)*100,1))
  
  # extract 5th and 95th percentile risks
  p5 <- round(quantile(var1,0.05),-1); p95 <- round(quantile(var1,0.95),-1)
  rr$p5.rain[i] <- p5
  rr$p5.mean[i] <- round((cp$allRRfit[as.character(p5)]-1)*100,1)
  rr$p5.ci[i] <- paste0(round((cp$allRRlow[as.character(p5)]-1)*100,1),"; ",round((cp$allRRhigh[as.character(p5)]-1)*100,1))
  rr$p95.rain[i] <- p95
  rr$p95.mean[i] <- round((cp$allRRfit[as.character(p95)]-1)*100,1)
  rr$p95.ci[i] <- paste0(round((cp$allRRlow[as.character(p95)]-1)*100,1),"; ",round((cp$allRRhigh[as.character(p95)]-1)*100,1))
}
rm(i,slist,day1,var,var1,ag1,pred1,predvar,argvar1,bvar,argvar,cp,p5,p95,df1,stab,pkval)

# save table of RRs
write.csv(rr,"tables/suppl_tab_clim-ctry_rr_centering.csv",row.names=T)
saveRDS(plt,"data/rr_dataframe/rr_values_country-clim.rds")


# individual ggplots
y1 <- c(0.9,1.4) # y axis boundaries
fill1 <- "gray70" #confidence intervals
col1 <- "black" # mean RR
#lab <- c(5,"",50,"",95,"")
mar <- c(0.1,0.03,0.1,0.1)
msize <- 15
xtext <- 8
ytext <- 10
atitle <- 13.5
nbrk <- 5
lab1 <- 2
vcol <- "darkblue"
yint <- 1.35
xax <- "28-day rain (mm)"
yax <- "RR (95%CI)"

# arg.C
arg.C <- plt[plt$loc=="arg.C",]
arg.C.val <- unlist(rain[rain$loc=="arg.C",c("p05","p50","p95")])
p.arg.C <- ggplot(data=arg.C,aes(x=rain,y=value,ymin=lower,ymax=upper)) +
  geom_hline(yintercept=1,linewidth=1) +
  geom_ribbon(fill=fill1,alpha=0.5) + 
  geom_line(color=col1,linewidth=1) +
  #scale_x_continuous(breaks=unlist(rain[rain$loc=="arg.C",-c(1,2)]),labels=lab) +
  scale_x_continuous(breaks=round(classIntervals(seq(min(arg.C$rain),max(arg.C$rain),10),
                                                 n=nbrk,style="equal")[[2]],-1)) +
  scale_y_continuous(n.breaks=6) +
  geom_vline(xintercept=arg.C.val,
             color=vcol,lwd=0.8,alpha=0.5) +
  geom_label(aes(x=arg.C.val[1],y=yint),label="5th",size=lab1) +
  geom_label(aes(x=arg.C.val[2],y=yint),label="50th",size=lab1) +
  geom_label(aes(x=arg.C.val[3],y=yint),label="95th",size=lab1) +
  labs(title=paste0(ctry[substr("arg.C",1,3)]," (n=",length(grep("arg.C",mtab$id)),")"),
       x=xax,y="") +
  coord_cartesian(ylim=y1,xlim=unlist(rain[rain$loc=="arg.C",c("p00","p99")])) +
  theme_bw() +
  theme(plot.title=element_text(size=msize,hjust=0.5,face="bold"),
        panel.grid.minor.x=element_blank(),
        plot.margin=unit(mar,"cm"),
        axis.text.x=element_text(size=xtext),
        axis.text.y=element_text(size=ytext),
        axis.title=element_text(size=atitle))
#p.arg.C

# bra.A 
bra.A <- plt[plt$loc=="bra.A",]
bra.A.val <- unlist(rain[rain$loc=="bra.A",c("p05","p50","p95")])
p.bra.A <- ggplot(data=bra.A,aes(x=rain,y=value,ymin=lower,ymax=upper)) +
  geom_hline(yintercept=1,linewidth=1) +
  geom_ribbon(fill=fill1,alpha=0.5) + 
  geom_line(color=col1,linewidth=1) +
  #scale_x_continuous(breaks=unlist(rain[rain$loc=="bra.A",-c(1,2)]),labels=lab) +
  scale_x_continuous(breaks=round(classIntervals(seq(min(bra.A$rain),max(bra.A$rain),10),
                                                 n=nbrk,style="equal")[[2]],-1)) +
  scale_y_continuous(n.breaks=6) +
  geom_vline(xintercept=bra.A.val,
             color=vcol,lwd=0.8,alpha=0.5) +
  geom_label(aes(x=bra.A.val[1],y=yint),label="5th",size=lab1) +
  geom_label(aes(x=bra.A.val[2],y=yint),label="50th",size=lab1) +
  geom_label(aes(x=bra.A.val[3],y=yint),label="95th",size=lab1) +
  labs(title=paste0(ctry[substr("bra.A",1,3)]," (n=",length(grep("bra.A",mtab$id)),")"),
       x="",y=yax) +
  coord_cartesian(ylim=y1,xlim=unlist(rain[rain$loc=="bra.A",c("p00","p99")])) +
  theme_bw() +
  theme(plot.title=element_text(size=msize,hjust=0.5,face="bold"),
        panel.grid.minor.x=element_blank(),
        plot.margin=unit(mar,"cm"),
        axis.text.x=element_text(size=xtext),
        axis.text.y=element_text(size=ytext),
        axis.title=element_text(size=atitle))
#p.bra.A

# bra.B 
bra.B <- plt[plt$loc=="bra.B",]
bra.B.val <- unlist(rain[rain$loc=="bra.B",c("p05","p50","p95")])
p.bra.B <- ggplot(data=bra.B,aes(x=rain,y=value,ymin=lower,ymax=upper)) +
  geom_hline(yintercept=1,linewidth=1) +
  geom_ribbon(fill=fill1,alpha=0.5) + 
  geom_line(color=col1,linewidth=1) +
  #scale_x_continuous(breaks=unlist(rain[rain$loc=="bra.B",-c(1,2)]),labels=lab) +
  scale_x_continuous(breaks=round(classIntervals(seq(min(bra.B$rain),max(bra.B$rain),10),
                                                 n=nbrk,style="equal")[[2]],-1)) +
  scale_y_continuous(n.breaks=6) +
  geom_vline(xintercept=bra.B.val,
             color=vcol,lwd=0.8,alpha=0.5) +
  geom_label(aes(x=bra.B.val[1],y=yint),label="5th",size=lab1) +
  geom_label(aes(x=bra.B.val[2],y=yint),label="50th",size=lab1) +
  geom_label(aes(x=bra.B.val[3],y=yint),label="95th",size=lab1) +
  labs(title=paste0(ctry[substr("bra.B",1,3)]," (n=",length(grep("bra.B",mtab$id)),")"),
       x="",y=yax) +
  coord_cartesian(ylim=y1,xlim=unlist(rain[rain$loc=="bra.B",c("p00","p99")])) +
  theme_bw() +
  theme(plot.title=element_text(size=msize,hjust=0.5,face="bold"),
        panel.grid.minor.x=element_blank(),
        plot.margin=unit(mar,"cm"),
        axis.text.x=element_text(size=xtext),
        axis.text.y=element_text(size=ytext),
        axis.title=element_text(size=atitle))
#p.bra.B

# bra.C 
bra.C <- plt[plt$loc=="bra.C",]
bra.C.val <- unlist(rain[rain$loc=="bra.C",c("p05","p50","p95")])
p.bra.C <- ggplot(data=bra.C,aes(x=rain,y=value,ymin=lower,ymax=upper)) +
  geom_hline(yintercept=1,linewidth=1) +
  geom_ribbon(fill=fill1,alpha=0.5) + 
  geom_line(color=col1,linewidth=1) +
  #scale_x_continuous(breaks=unlist(rain[rain$loc=="bra.C",-c(1,2)]),labels=lab) +
  scale_x_continuous(breaks=round(classIntervals(seq(min(bra.C$rain),max(bra.C$rain),10),
                                                 n=nbrk,style="equal")[[2]],-1)) +
  scale_y_continuous(n.breaks=6) +
  geom_vline(xintercept=bra.C.val,
             color=vcol,lwd=0.8,alpha=0.5) +
  geom_label(aes(x=bra.C.val[1],y=yint),label="5th",size=lab1) +
  geom_label(aes(x=bra.C.val[2],y=yint),label="50th",size=lab1) +
  geom_label(aes(x=bra.C.val[3],y=yint),label="95th",size=lab1) +
  labs(title=paste0(ctry[substr("bra.C",1,3)]," (n=",length(grep("bra.C",mtab$id)),")"),
       x="",y=yax) +
  coord_cartesian(ylim=y1,xlim=unlist(rain[rain$loc=="bra.C",c("p00","p99")])) +
  theme_bw() +
  theme(plot.title=element_text(size=msize,hjust=0.5,face="bold"),
        panel.grid.minor.x=element_blank(),
        plot.margin=unit(mar,"cm"),
        axis.text.x=element_text(size=xtext),
        axis.text.y=element_text(size=ytext),
        axis.title=element_text(size=atitle))
#p.bra.C

# cri.A 
cri.A <- plt[plt$loc=="cri.A",]
cri.A.val <- unlist(rain[rain$loc=="cri.A",c("p05","p50","p95")])
p.cri.A <- ggplot(data=cri.A,aes(x=rain,y=value,ymin=lower,ymax=upper)) +
  geom_hline(yintercept=1,linewidth=1) +
  geom_ribbon(fill=fill1,alpha=0.5) + 
  geom_line(color=col1,linewidth=1) +
  #scale_x_continuous(breaks=unlist(rain[rain$loc=="cri.A",-c(1,2)]),labels=lab) +
  scale_x_continuous(breaks=round(classIntervals(seq(min(cri.A$rain),max(cri.A$rain),10),
                                                 n=nbrk,style="equal")[[2]],-1)) +
  scale_y_continuous(n.breaks=6) +
  geom_vline(xintercept=cri.A.val,
             color=vcol,lwd=0.8,alpha=0.5) +
  geom_label(aes(x=cri.A.val[1],y=yint),label="5th",size=lab1) +
  geom_label(aes(x=cri.A.val[2],y=yint),label="50th",size=lab1) +
  geom_label(aes(x=cri.A.val[3],y=yint),label="95th",size=lab1) +
  labs(title=paste0(ctry[substr("cri.A",1,3)]," (n=",length(grep("cri.A",mtab$id)),")"),
       x="",y="") +
  coord_cartesian(ylim=y1,xlim=unlist(rain[rain$loc=="cri.A",c("p00","p99")])) +
  theme_bw() +
  theme(plot.title=element_text(size=msize,hjust=0.5,face="bold"),
        panel.grid.minor.x=element_blank(),
        plot.margin=unit(mar,"cm"),
        axis.text.x=element_text(size=xtext),
        axis.text.y=element_text(size=ytext),
        axis.title=element_text(size=atitle))
#p.cri.A

# ind.A 
ind.A <- plt[plt$loc=="ind.A",]
ind.A.val <- unlist(rain[rain$loc=="ind.A",c("p05","p50","p95")])
p.ind.A <- ggplot(data=ind.A,aes(x=rain,y=value,ymin=lower,ymax=upper)) +
  geom_hline(yintercept=1,linewidth=1) +
  geom_ribbon(fill=fill1,alpha=0.5) + 
  geom_line(color=col1,linewidth=1) +
  #scale_x_continuous(breaks=unlist(rain[rain$loc=="ind.A",-c(1,2)]),labels=lab) +
  scale_x_continuous(breaks=round(classIntervals(seq(min(ind.A$rain),max(ind.A$rain),10),
                                                 n=nbrk,style="equal")[[2]],-1)) +
  scale_y_continuous(n.breaks=6) +
  geom_vline(xintercept=ind.A.val,
             color=vcol,lwd=0.8,alpha=0.5) +
  geom_label(aes(x=ind.A.val[1],y=yint),label="5th",size=lab1) +
  geom_label(aes(x=ind.A.val[2],y=yint),label="50th",size=lab1) +
  geom_label(aes(x=ind.A.val[3],y=yint),label="95th",size=lab1) +
  labs(title=paste0(ctry[substr("ind.A",1,3)]," (n=",length(grep("ind.A",mtab$id)),")"),
       x="",y="") +
  coord_cartesian(ylim=y1,xlim=unlist(rain[rain$loc=="ind.A",c("p00","p99")])) +
  theme_bw() +
  theme(plot.title=element_text(size=msize,hjust=0.5,face="bold"),
        panel.grid.minor.x=element_blank(),
        plot.margin=unit(mar,"cm"),
        axis.text.x=element_text(size=xtext),
        axis.text.y=element_text(size=ytext),
        axis.title=element_text(size=atitle))
#p.ind.A

# ind.B 
ind.B <- plt[plt$loc=="ind.B",]
ind.B.val <- unlist(rain[rain$loc=="ind.B",c("p05","p50","p95")])
p.ind.B <- ggplot(data=ind.B,aes(x=rain,y=value,ymin=lower,ymax=upper)) +
  geom_hline(yintercept=1,linewidth=1) +
  geom_ribbon(fill=fill1,alpha=0.5) + 
  geom_line(color=col1,linewidth=1) +
  #scale_x_continuous(breaks=unlist(rain[rain$loc=="ind.B",-c(1,2)]),labels=lab) +
  scale_x_continuous(breaks=round(classIntervals(seq(min(ind.B$rain),max(ind.B$rain),10),
                                                 n=nbrk,style="equal")[[2]],-1)) +
  scale_y_continuous(n.breaks=6) +
  geom_vline(xintercept=ind.B.val,
             color=vcol,lwd=0.8,alpha=0.5) +
  geom_label(aes(x=ind.B.val[1],y=yint),label="5th",size=lab1) +
  geom_label(aes(x=ind.B.val[2],y=yint),label="50th",size=lab1) +
  geom_label(aes(x=ind.B.val[3],y=yint),label="95th",size=lab1) +
  labs(title=paste0(ctry[substr("ind.B",1,3)]," (n=",length(grep("ind.B",mtab$id)),")"),
       x="",y="") +
  coord_cartesian(ylim=y1,xlim=unlist(rain[rain$loc=="ind.B",c("p00","p99")])) +
  theme_bw() +
  theme(plot.title=element_text(size=msize,hjust=0.5,face="bold"),
        panel.grid.minor.x=element_blank(),
        plot.margin=unit(mar,"cm"),
        axis.text.x=element_text(size=xtext),
        axis.text.y=element_text(size=ytext),
        axis.title=element_text(size=atitle))
#p.ind.B

# ind.C 
ind.C <- plt[plt$loc=="ind.C",]
ind.C.val <- unlist(rain[rain$loc=="ind.C",c("p05","p50","p95")])
p.ind.C <- ggplot(data=ind.C,aes(x=rain,y=value,ymin=lower,ymax=upper)) +
  geom_hline(yintercept=1,linewidth=1) +
  geom_ribbon(fill=fill1,alpha=0.5) + 
  geom_line(color=col1,linewidth=1) +
  #scale_x_continuous(breaks=unlist(rain[rain$loc=="ind.C",-c(1,2)]),labels=lab) +
  scale_x_continuous(breaks=round(classIntervals(seq(min(ind.C$rain),max(ind.C$rain),10),
                                                 n=nbrk,style="equal")[[2]],-1)) +
  scale_y_continuous(n.breaks=6) +
  geom_vline(xintercept=ind.C.val,
             color=vcol,lwd=0.8,alpha=0.5) +
  geom_label(aes(x=ind.C.val[1],y=yint),label="5th",size=lab1) +
  geom_label(aes(x=ind.C.val[2],y=yint),label="50th",size=lab1) +
  geom_label(aes(x=ind.C.val[3],y=yint),label="95th",size=lab1) +
  labs(title=paste0(ctry[substr("ind.C",1,3)]," (n=",length(grep("ind.C",mtab$id)),")"),
       x="",y="") +
  coord_cartesian(ylim=y1,xlim=unlist(rain[rain$loc=="ind.C",c("p00","p99")])) +
  theme_bw() +
  theme(plot.title=element_text(size=msize,hjust=0.5,face="bold"),
        panel.grid.minor.x=element_blank(),
        plot.margin=unit(mar,"cm"),
        axis.text.x=element_text(size=xtext),
        axis.text.y=element_text(size=ytext),
        axis.title=element_text(size=atitle))
#p.ind.C

# per.A 
per.A <- plt[plt$loc=="per.A",]
per.A.val <- unlist(rain[rain$loc=="per.A",c("p05","p50","p95")])
p.per.A <- ggplot(data=per.A,aes(x=rain,y=value,ymin=lower,ymax=upper)) +
  geom_hline(yintercept=1,linewidth=1) +
  geom_ribbon(fill=fill1,alpha=0.5) + 
  geom_line(color=col1,linewidth=1) +
  #scale_x_continuous(breaks=unlist(rain[rain$loc=="per.A",-c(1,2)]),labels=lab) +
  scale_x_continuous(breaks=round(classIntervals(seq(min(per.A$rain),max(per.A$rain),10),
                                                 n=nbrk,style="equal")[[2]],-1)) +
  scale_y_continuous(n.breaks=6) +
  geom_vline(xintercept=per.A.val,
             color=vcol,lwd=0.8,alpha=0.5) +
  geom_label(aes(x=per.A.val[1],y=yint),label="5th",size=lab1) +
  geom_label(aes(x=per.A.val[2],y=yint),label="50th",size=lab1) +
  geom_label(aes(x=per.A.val[3],y=yint),label="95th",size=lab1) +
  labs(title=paste0(ctry[substr("per.A",1,3)]," (n=",length(grep("per.A",mtab$id)),")"),
       x="",y="") +
  coord_cartesian(ylim=y1,xlim=unlist(rain[rain$loc=="per.A",c("p00","p99")])) +
  theme_bw() +
  theme(plot.title=element_text(size=msize,hjust=0.5,face="bold"),
        panel.grid.minor.x=element_blank(),
        plot.margin=unit(mar,"cm"),
        axis.text.x=element_text(size=xtext),
        axis.text.y=element_text(size=ytext),
        axis.title=element_text(size=atitle))
#p.per.A

# per.B 
per.B <- plt[plt$loc=="per.B",]
per.B.val <- unlist(rain[rain$loc=="per.B",c("p05","p50","p95")])
p.per.B <- ggplot(data=per.B,aes(x=rain,y=value,ymin=lower,ymax=upper)) +
  geom_hline(yintercept=1,linewidth=1) +
  geom_ribbon(fill=fill1,alpha=0.5) + 
  geom_line(color=col1,linewidth=1) +
  #scale_x_continuous(breaks=unlist(rain[rain$loc=="per.B",-c(1,2)]),labels=lab) +
  scale_x_continuous(breaks=round(classIntervals(seq(min(per.B$rain),max(per.B$rain),10),
                                                 n=nbrk,style="equal")[[2]],-1)) +
  scale_y_continuous(n.breaks=6) +
  geom_vline(xintercept=per.B.val,
             color=vcol,lwd=0.8,alpha=0.5) +
  geom_label(aes(x=per.B.val[1],y=yint),label="5th",size=lab1) +
  geom_label(aes(x=per.B.val[2],y=yint),label="50th",size=lab1) +
  geom_label(aes(x=per.B.val[3],y=yint),label="95th",size=lab1) +
  labs(title=paste0(ctry[substr("per.B",1,3)]," (n=",length(grep("per.B",mtab$id)),")"),
       x="",y="") +
  coord_cartesian(ylim=y1,xlim=unlist(rain[rain$loc=="per.B",c("p00","p99")])) +
  theme_bw() +
  theme(plot.title=element_text(size=msize,hjust=0.5,face="bold"),
        panel.grid.minor.x=element_blank(),
        plot.margin=unit(mar,"cm"),
        axis.text.x=element_text(size=xtext),
        axis.text.y=element_text(size=ytext),
        axis.title=element_text(size=atitle))
#p.per.B

# phl.A 
phl.A <- plt[plt$loc=="phl.A",]
phl.A.val <- unlist(rain[rain$loc=="phl.A",c("p05","p50","p95")])
p.phl.A <- ggplot(data=phl.A,aes(x=rain,y=value,ymin=lower,ymax=upper)) +
  geom_hline(yintercept=1,linewidth=1) +
  geom_ribbon(fill=fill1,alpha=0.5) + 
  geom_line(color=col1,linewidth=1) +
  #scale_x_continuous(breaks=unlist(rain[rain$loc=="phl.A",-c(1,2)]),labels=lab) +
  scale_x_continuous(breaks=round(classIntervals(seq(min(phl.A$rain),max(phl.A$rain),10),
                                                 n=nbrk,style="equal")[[2]],-1)) +
  scale_y_continuous(n.breaks=6) +
  geom_vline(xintercept=phl.A.val,
             color=vcol,lwd=0.8,alpha=0.5) +
  geom_label(aes(x=phl.A.val[1],y=yint),label="5th",size=lab1) +
  geom_label(aes(x=phl.A.val[2],y=yint),label="50th",size=lab1) +
  geom_label(aes(x=phl.A.val[3],y=yint),label="95th",size=lab1) +
  labs(title=paste0(ctry[substr("phl.A",1,3)]," (n=",length(grep("phl.A",mtab$id)),")"),
       x="",y="") +
  coord_cartesian(ylim=y1,xlim=unlist(rain[rain$loc=="phl.A",c("p00","p99")])) +
  theme_bw() +
  theme(plot.title=element_text(size=msize,hjust=0.5,face="bold"),
        panel.grid.minor.x=element_blank(),
        plot.margin=unit(mar,"cm"),
        axis.text.x=element_text(size=xtext),
        axis.text.y=element_text(size=ytext),
        axis.title=element_text(size=atitle))
#p.phl.A

# tha.A 
tha.A <- plt[plt$loc=="tha.A",]
tha.A.val <- unlist(rain[rain$loc=="tha.A",c("p05","p50","p95")])
p.tha.A <- ggplot(data=tha.A,aes(x=rain,y=value,ymin=lower,ymax=upper)) +
  geom_hline(yintercept=1,linewidth=1) +
  geom_ribbon(fill=fill1,alpha=0.5) + 
  geom_line(color=col1,linewidth=1) +
  #scale_x_continuous(breaks=unlist(rain[rain$loc=="tha.A",-c(1,2)]),labels=lab) +
  scale_x_continuous(breaks=round(classIntervals(seq(min(tha.A$rain),max(tha.A$rain),10),
                                                 n=nbrk,style="equal")[[2]],-1)) +
  scale_y_continuous(n.breaks=6) +
  geom_vline(xintercept=tha.A.val,
             color=vcol,lwd=0.8,alpha=0.5) +
  geom_label(aes(x=tha.A.val[1],y=yint),label="5th",size=lab1) +
  geom_label(aes(x=tha.A.val[2],y=yint),label="50th",size=lab1) +
  geom_label(aes(x=tha.A.val[3],y=yint),label="95th",size=lab1) +
  labs(title=paste0(ctry[substr("tha.A",1,3)]," (n=",length(grep("tha.A",mtab$id)),")"),
       x="",y="") +
  coord_cartesian(ylim=y1,xlim=unlist(rain[rain$loc=="tha.A",c("p00","p99")])) +
  theme_bw() +
  theme(plot.title=element_text(size=msize,hjust=0.5,face="bold"),
        panel.grid.minor.x=element_blank(),
        plot.margin=unit(mar,"cm"),
        axis.text.x=element_text(size=xtext),
        axis.text.y=element_text(size=ytext),
        axis.title=element_text(size=atitle))
#p.tha.A

# zaf.B 
zaf.B <- plt[plt$loc=="zaf.B",]
zaf.B.val <- unlist(rain[rain$loc=="zaf.B",c("p05","p50","p95")])
p.zaf.B <- ggplot(data=zaf.B,aes(x=rain,y=value,ymin=lower,ymax=upper)) +
  geom_hline(yintercept=1,linewidth=1) +
  geom_ribbon(fill=fill1,alpha=0.5) + 
  geom_line(color=col1,linewidth=1) +
  #scale_x_continuous(breaks=unlist(rain[rain$loc=="zaf.B",-c(1,2)]),labels=lab) +
  scale_x_continuous(breaks=round(classIntervals(seq(min(zaf.B$rain),max(zaf.B$rain),10),
                                                 n=nbrk,style="equal")[[2]],-1)) +
  scale_y_continuous(n.breaks=6) +
  geom_vline(xintercept=zaf.B.val,
             color=vcol,lwd=0.8,alpha=0.5) +
  geom_label(aes(x=zaf.B.val[1],y=yint),label="5th",size=lab1) +
  geom_label(aes(x=zaf.B.val[2],y=yint),label="50th",size=lab1) +
  geom_label(aes(x=zaf.B.val[3],y=yint),label="95th",size=lab1) +
  labs(title=paste0(ctry[substr("zaf.B",1,3)]," (n=",length(grep("zaf.B",mtab$id)),")"),
       x="",y="") +
  coord_cartesian(ylim=y1,xlim=unlist(rain[rain$loc=="zaf.B",c("p00","p99")])) +
  theme_bw() +
  theme(plot.title=element_text(size=msize,hjust=0.5,face="bold"),
        panel.grid.minor.x=element_blank(),
        plot.margin=unit(mar,"cm"),
        axis.text.x=element_text(size=xtext),
        axis.text.y=element_text(size=ytext),
        axis.title=element_text(size=atitle))
#p.zaf.B

# zaf.C 
zaf.C <- plt[plt$loc=="zaf.C",]
zaf.C.val <- unlist(rain[rain$loc=="zaf.C",c("p05","p50","p95")])
p.zaf.C <- ggplot(data=zaf.C,aes(x=rain,y=value,ymin=lower,ymax=upper)) +
  geom_hline(yintercept=1,linewidth=1) +
  geom_ribbon(fill=fill1,alpha=0.5) + 
  geom_line(color=col1,linewidth=1) +
  #scale_x_continuous(breaks=unlist(rain[rain$loc=="zaf.C",-c(1,2)]),labels=lab) +
  scale_x_continuous(breaks=round(classIntervals(seq(min(zaf.C$rain),max(zaf.C$rain),10),
                                                 n=nbrk,style="equal")[[2]],-1)) +
  scale_y_continuous(n.breaks=6) +
  geom_vline(xintercept=zaf.C.val,
             color=vcol,lwd=0.8,alpha=0.5) +
  geom_label(aes(x=zaf.C.val[1],y=yint),label="5th",size=lab1) +
  geom_label(aes(x=zaf.C.val[2],y=yint),label="50th",size=lab1) +
  geom_label(aes(x=zaf.C.val[3],y=yint),label="95th",size=lab1) +
  labs(title=paste0(ctry[substr("zaf.C",1,3)]," (n=",length(grep("zaf.C",mtab$id)),")"),
       x="",y="") +
  coord_cartesian(ylim=y1,xlim=unlist(rain[rain$loc=="zaf.C",c("p00","p99")])) +
  theme_bw() +
  theme(plot.title=element_text(size=msize,hjust=0.5,face="bold"),
        panel.grid.minor.x=element_blank(),
        plot.margin=unit(mar,"cm"),
        axis.text.x=element_text(size=xtext),
        axis.text.y=element_text(size=ytext),
        axis.title=element_text(size=atitle))
#p.zaf.C

# empty plots
bplot <- ggplot() + theme_void()

# compile by rows
row1 <- plot_grid(p.bra.A, p.ind.A, p.per.A, p.cri.A, p.phl.A, p.tha.A, ncol=6)
row2 <- plot_grid(p.bra.B, p.ind.B, p.per.B, p.zaf.B, bplot, bplot, ncol=6)
row3 <- plot_grid(p.bra.C, p.ind.C, p.arg.C, p.zaf.C, bplot, bplot, ncol=6)
mainplot <- plot_grid(row1,row2,row3,ncol=1,labels=c("A","B","C"),
                      label_size=25,hjust=c(-0.2,-0.2,-0.2),vjust=c(1,1,1))
#mainplot
save_plot("figures/main_fig_risk-functions_ctry_clim.png",mainplot,base_asp=3,base_height=8,base_width=15)
