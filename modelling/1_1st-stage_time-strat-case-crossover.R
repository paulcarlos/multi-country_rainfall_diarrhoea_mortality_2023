library(splines)
library(dlnm)
library(mixmeta)
library(gnm)
library(lubridate)
library(tsModel)

# ERA5-Land complete data list #####
arg <- readRDS("data/list_dataset/arg_list_clim-cat_era5-land.rds")
bra <- readRDS("data/list_dataset/bra_list_clim-cat_era5-land.rds")
cri <- readRDS("data/list_dataset/cri_list_clim-cat_era5-land.rds") 
ecu <- readRDS("data/list_dataset/ecu_list_clim-cat_era5-land.rds")
ind <- readRDS("data/list_dataset/ind_list_clim-cat_era5-land.rds")
per <- readRDS("data/list_dataset/per_list_clim-cat_era5-land.rds")
phl <- readRDS("data/list_dataset/phl_list_clim-cat_era5-land.rds")
tha <- readRDS("data/list_dataset/tha_list_clim-cat_era5-land.rds")
zaf <- readRDS("data/list_dataset/zaf_list_clim-cat_era5-land.rds")
dlist <- c(arg,bra,cri,ecu,ind,per,phl,tha,zaf)
rm(arg,bra,cri,ecu,ind,per,phl,tha,zaf)

# metadata ####
mtab <- read.csv("vars_meta-analysis_country-subclim.csv",stringsAsFactors=F)
mtab <- mtab[mtab$id%in%names(dlist),]
rownames(mtab) <- 1:nrow(mtab)

# model specs ####
rkn <- c(0.5,0.9)
tkn <- c(0.33,0.67)
lagnk <- 3
tlag <- 28 # temperature lag

# storage for coefficients and vcov matrices
coef1 <- matrix(data=NA,nrow=nrow(mtab),ncol=length(rkn)+1,dimnames=list(mtab$id))
vcov1 <- vector("list",nrow(mtab)); names(vcov1) <- mtab$id

# loop first stage ####
lag1 <- 28
strat <- "str_3mo"
par(mfrow=c(6,5),mar=c(3,2,3,2),oma=c(2,3,1,1))
for (i in 1:nrow(mtab)) {
  cat(i," ")
  # get data
  sdat <- dlist[[mtab$id[i]]]
  sdat$stratum <- sdat[,strat]
  
  # rainfall 
  pr2 <- sdat[,paste0("pr",lag1)]
  ob1 <- onebasis(pr2,fun="ns",knots=quantile(pr2,rkn,na.rm=TRUE))
  
  # temperature
  t2 <- sdat$t2m
  cb1 <- crossbasis(t2,lag=tlag,argvar=list(fun="ns",knots=quantile(t2,tkn)),arglag=list(fun="ns",knots=logknots(tlag,nk=lagnk)))
  
  # model +as.factor(month)
  mod1 <- gnm(mort~ob1+cb1+as.factor(month),data=sdat,family=quasipoisson,eliminate=stratum,na.action=na.exclude)
  cp <- crosspred(ob1,mod1,cen=median(pr2,na.rm=TRUE)) # rain
  
  # adjust for centering
  predvar <- quantile(pr2,0:90/100,na.rm=T)
  argvar1 <- list(x=predvar,fun="ns",knots=quantile(pr2,rkn,na.rm=T),Bound=range(pr2,na.rm=T))
  bvar <- do.call(onebasis,argvar1)
  cen <- round(predvar[which.min((bvar%*%cp$coefficients))])
  cp <- crosspred(ob1,mod1,cen=cen,by=10)
  
  # plot 
  plot(cp,ylim=c(0.8,1.5),main=loc[i])
  abline(v=cen,col="darkgray",lty=2)
  
  # save coefficients and vcov
  coef1[i,] <- coef(cp)
  vcov1[[i]] <- vcov(cp)
}
rm(i,sdat,pr2,ob1,t2,cb1,mod1,cp,predvar,argvar1,bvar,cen)

# meta-analysis all locations
mtot1 <- mixmeta(formula=coef1~clim,random=~1|ctry,data=mtab,S=vcov1)
summary(mtot1)

# save
save(dlist,mtab,rkn,tkn,lagnk,tlag,lag1,coef1,vcov1,file="data/mod_out/model-out_1st-stage.rda")
