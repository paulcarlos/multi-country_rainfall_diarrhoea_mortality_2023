library(sf)
library(raster)

ctry <- "phl"
lvl <- paste0(c("adm","ADM_"),2)

# load data
shp <- st_read(paste0("gadm41_",toupper(ctry),".gpkg"), layer=paste0("ADM_",lvl[2]),quiet=T) # country shapefile from GADM 4.1
pop <- readRDS(grep("pop-density",list.files(paste0("data/matrix_popd/",ctry),full.names=T),value=T)) # get population density matrix
clim <- readRDS(grep("kg-cat",list.files(paste0("data/matrix_popd/",ctry),full.names=T),value=T)) # get climate classification matrix
kg <- read.csv("kp_class_beck2018.csv",sep=";",stringsAsFactors=FALSE) # climate category
mdat <- read.csv("phl_adm3_names.csv"),stringsAsFactors=F) #metadata
adm <- readRDS(grep(paste0("era5land_",lvl[1]),list.files(paste0("data/matrix_loc/",ctry),full.names=T),value=T))  #matrices of sub-locations

# convert into raster
lon <- as.numeric(dimnames(pop)[[1]])
lat <- as.numeric(dimnames(pop)[[2]])
rpop <- raster(t(pop),xmn=min(lon),xmx=max(lon),ymn=min(lat),ymx=max(lat),crs=crs(shp))
#plot(rpop)
rclim <- raster(t(clim),xmn=min(lon),xmx=max(lon),ymn=min(lat),ymx=max(lat),crs=crs(shp))
#plot(rclim)

# loop
mdat$kg.era5 <- NA
for (i in seq(adm)) {
  # sample plot
  #plot(rclim)
  #plot(shp[i,"geom"],add=T)
  # get name of sub-location
  nm <- names(adm)[i] 
  # select location matrix
  loc <- adm[[nm]]
  #dimnames(loc) <- dimnames(clim)
  # multiply with pop and clim matrices to focus on location values only
  m.pop <- pop*loc
  m.clim <- clim*loc
  m.clim[is.na(m.clim)] <- m.pop[is.na(m.pop)] <- 0 #replace NAs with 0s
  # convert into dataframe
  df.loc <- as.data.frame.table(loc,stringsAsFactors = F)
  df.loc$clim <- as.data.frame.table(m.clim,stringsAsFactors = F)[,"Freq"]
  df.loc$pop <- as.data.frame.table(m.pop,stringsAsFactors = F)[,"Freq"]
  df.clim <- df.loc[df.loc$Freq==1,]
  ag <- aggregate(pop~clim,df.clim,"sum")
  ag1 <- ag[!ag$clim==0,]
  if (nrow(ag1>0)) {
    mdat$kg.era5[mdat[,paste0(lvl[1],".id.gadm")]==nm] <- kg$abbr_kp[kg$val==ag1$clim[ag1$pop==max(ag1$pop,na.rm=T)]]
  } else {
    buf <- st_buffer(shp[shp$GID_2==nm,],dist=10000)
    #plot(rclim)
    #plot(rclim,xlim=c(-83,-80),ylim=c(-3,-1))
    #plot(buf$geom,add=T)
    ext.clim <- extract(rclim,buf,df=TRUE,cellnumbers=TRUE,weights=TRUE)
    ext.pop <- extract(rpop,buf,df=TRUE,cellnumbers=TRUE,weights=TRUE)
    ext.clim$pop <- ext.pop$layer
    rval <- which(ext.clim$layer>0) # rows with climate categories
    sext <- ext.clim[rval,]
    mdat$kg.era5[mdat[,paste0(lvl[1],".id.gadm")]==nm] <- kg$abbr_kp[kg$val==sext$layer[sext$pop==max(sext$pop)]]
  }
  cat(i," ")
}
rm(i,loc,m.pop,m.clim,df.loc,ag)
sum(is.na(mdat$kg.era5))
aggregate(dia.mort~kg.era5,mdat,"sum")

# update metadata
write.csv(mdat,"phl_adm3_names.csv",row.names=F)
