library(ncdf4)
library(lubridate)
library(matrixStats)
library(weathermetrics)

nproj <- "+proj=longlat +datum=WGS84 +no_defs"
tz1 <- "Asia/Manila" # time zone

# reference table
rtab <- read.csv("phl_adm3_names.csv",stringsAsFactors=F)
kp <- read.csv("kp_class_beck2018.csv",sep=";",stringsAsFactors=FALSE)

# matrix of phil adm3
mlist <- readRDS("phl_matrix_adm3.rds")

# matrix of population
pop <- readRDS("un-wpp_pop-density_era5-land-aligned_avg_2005-2020_phl.rds")


###################### total rainfall ########################
fn1 <- list.files("ecmwf/era5-land/prcp/phl/",pattern="_ph",full.names=TRUE)
for (i in seq(fn1)) {
  cat("\n",i,"\n")
  # get netcdf files
  nc <- nc_open(fn1[i])
  #nc
  # dimensions
  lon <- as.vector(ncvar_get(nc,varid="longitude"))
  lat <- as.vector(ncvar_get(nc,varid="latitude"))
  tim <- as.vector(ncvar_get(nc,varid="time"))
  # create time
  rtim <- as.POSIXct(gsub("hours since ","",nc$dim$time$units),tz="UTC")
  ctim <- rtim+hours(tim)
  # create data frame
  df1 <- data.frame(matrix(NA,nrow=length(ctim),ncol=length(mlist)+1))
  colnames(df1) <- c("hour",names(mlist))
  df1$hour <- ctim
  # get vars
  var1 <- ncvar_get(nc,varid="tp")
  dimnames(var1) <- list(lon,lat,tim)
  # by location
  for (j in seq(mlist)) {
    mat1 <- mlist[[j]]
    dimnames(mat1) <- list(lon,lat)
    mat2 <- as.data.frame(as.table(mat1))
    mat3 <- mat2[mat2$Freq>0,]
    if (nrow(mat3)>1) {
      df2 <- data.frame(matrix(NA,nrow=length(ctim),ncol=nrow(mat3)))
      wt1 <- NA
      for (k in 1:nrow(mat3)) {
        sval1 <- var1[as.character(mat3$Var1[k]),as.character(mat3$Var2[k]),] # ACCUMULATED FROM Day 1 0H to Day 2 0H (25 hours)
        sval2 <- tsModel::Lag(sval1,1) 
        sval3 <- sval1-sval2 # subtract accumulation
        sval4 <- sval3
        sval4[hour(ctim)==1] <- sval1[hour(ctim)==1] - sval3[hour(ctim)==0] # 1H is corrected
        sval4[2] <- sval1[2]
        df2[,k] <- sval4 
        wt1[k] <- pop[as.character(mat3$Var1[k]),as.character(mat3$Var2[k])]
      }
      wt1[is.na(wt1)] <- 0
      df1[,1+j] <- apply(df2,1,weighted.mean,w=wt1,na.rm=T)
    } else {
      sval1 <- var1[as.character(mat3$Var1[k]),as.character(mat3$Var2[k]),]
      sval2 <- tsModel::Lag(sval1,1)
      sval3 <- sval1-sval2 # subtract accumulation
      sval4 <- sval3
      sval4[hour(ctim)==1] <- sval1[hour(ctim)==1] - sval3[hour(ctim)==0]
      sval4[2] <- sval1[2]
      df1[,1+j] <- sval4
    }
    #weightedMean(var1[,,1],w=mat2,na.rm=TRUE)
    #df1[,1+j] <- apply(var1,3,function(z)weightedMean(x=z,w=mat2,na.rm=TRUE))
    cat(paste0(i,"-",j)," ")
  }
  saveRDS(df1,paste0("data/extracted_era5land/phl/pr_hourly_era5land_",min(year(ctim)),"_",max(year(ctim)),"_phl_adm3.rds"))
  # remove files to reduce memory
  nc_close(nc)
  rm(nc,var1)
  gc()
}
rm(i,lon,lat,tim,rtim,ctim,df1,j,mat1,mat2,mat3,k,wt1,df2)

# convert hourly into daily
pr1 <- readRDS("data/extracted_era5land/phl/pr_hourly_era5land_2005_2012_phl_adm3.rds")
pr2 <- readRDS("data/extracted_era5land/phl/pr_hourly_era5land_2013_2020_phl_adm3.rds")
pr <- rbind(pr1,pr2)
rm(pr1,pr2); gc()
cnum <- match(names(which(is.na(colSums(pr[,-1])))),colnames(pr))
rtab$adm3.name[rtab$adm3.code%in%colnames(pr)[cnum]]
rtab$mort.dia[rtab$adm3.code%in%colnames(pr)[cnum]]
pr <- pr[,-cnum]
pr$date <- as.Date(pr$hour,tz=tz1)
ag1 <- aggregate(.~date,data=pr[,-1],FUN="sum")
saveRDS(ag1,"data/extracted_era5land/phl/pr_daily_era5land_2005_2020_phl_adm3.rds")

rm(pr1,pr2,pr,cnum,ag1)
gc()

# convert to clim cat
ag1 <- readRDS("data/extracted_era5land/phl/pr_daily_era5land_2005_2020_phl_adm3.rds")
#clim <- c("Af","Am","Aw")
clim <- names(clist)
df <- data.frame("date"=ag1$date)
for (i in clim) {
  loc1 <- intersect(colnames(ag1),clist[[i]])
  sdat <- ag1[,loc1]
  wt1 <- rtab$mort.dia[rtab$adm3.code%in%loc1]
  val <- apply(sdat,1,weighted.mean,w=wt1)
  df[,i] <- val*1000 #from meter to mm
}
rm(i,loc1,sdat,wt1,val)
plot(df$Cwb)
saveRDS(df,"data/extracted_era5land/phl/pr_daily_era5land_2005_2020_phl_clim-cat.rds")
rm(ag1,clim,df)
gc()
