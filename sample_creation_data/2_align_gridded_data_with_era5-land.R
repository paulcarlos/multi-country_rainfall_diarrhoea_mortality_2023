library(ncdf4)
library(sf)
library(raster)

ctry <- "phl" # country
vnam <- "t2m" # target variable, t2m = temperature

# upload shapefile
shp <- st_read(paste0("gadm41_",toupper(ctry),".gpkg"),layer="ADM_ADM_0",quiet=T)  # using GADM 4.1

# upload sample ERA5-Land raster
nc <- nc_open(list.files(paste0("ecmwf/era5-land/",vnam,"/",ctry,"/"),full.names=T)[1]) # get file names of downloaded NC files
lon <- as.vector(ncvar_get(nc,varid="longitude"))
lat <- as.vector(ncvar_get(nc,varid="latitude"))
mat <- ncvar_get(nc,varid=vnam,start=c(1,1,1),count=c(length(lon),length(lat),1)) # get 1 hour sample raster
dimnames(mat) <- list(lon,lat)
mat1 <- t(mat)
ras <- raster(mat1,xmx=max(lon),xmn=min(lon),ymx=max(lat),ymn=min(lat),crs=crs(shp))
plot(ras); plot(shp$geom,add=T)
nc_close(nc); rm(nc)



################ population matrix era5-land aligned ##############
nc <- nc_open("gpw_v4_population_density_adjusted_rev11_2pt5_min.nc")  # CIESIN GPW data
#nc #check details

# get dimensions
#ncvar_get(nc,varid="raster") #years: 1=2000; 2=2005; 3=2010; 4=2015; 5=2020
#yr <- 1:5
lon.pop <- ncvar_get(nc,varid="longitude")
lat.pop <- ncvar_get(nc,varid="latitude")

# limit longitude and latitude
lon.dif <- which(lon.pop>=floor(min(as.numeric(dimnames(mat1)[[2]]))) & lon.pop<=ceiling(max(as.numeric(dimnames(mat1)[[2]]))))
lat.dif <- which(lat.pop>=floor(min(as.numeric(dimnames(mat1)[[1]]))) & lat.pop<=ceiling(max(as.numeric(dimnames(mat1)[[1]]))))

# get matrix based on longitude & latitude from ERA5
yr1 <- seq(2000,2020,5)
yrstart <- 2005
yrend <- 2020 
pop <- ncvar_get(nc,varid="UN WPP-Adjusted Population Density, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes",
                 start=c(min(lon.dif),min(lat.dif),which(yr1==yrstart)),
                 count=c(length(lon.dif),length(lat.dif),which(yr1==yrend)-which(yr1==yrstart)+1))
pop.avg <- apply(pop,1:2,FUN="mean") #take average for year 2010,2015,2020
ras.pop <- raster(t(pop.avg),xmn=min(lon.pop[lon.dif]),xmx=max(lon.pop[lon.dif]),
                  ymn=min(lat.pop[lat.dif]),ymx=max(lat.pop[lat.dif]),crs=crs(shp))
plot(ras.pop)
plot(shp$geom,add=T)

# re-sample raster
popd <- projectRaster(ras.pop,ras,crs=crs(shp),method="bilinear")
plot(popd)
plot(shp$geom,add=T)

# create matrix
mat.pop <- t(as.matrix(popd)) #re-transposed
dimnames(mat.pop) <- list(lon,lat)
saveRDS(mat.pop,file=paste0("data/matrix_popd/",ctry,"/un-wpp_pop-density_era5-land-aligned_avg_",
                            yrstart,"-",yrend,"_",ctry,".rds"))

# remove
nc_close(nc)
#rm(list=setdiff(ls(),c("mat1","ras")))



################ Koppen-Geiger classification matrix era5-land aligned ##############

# Koppen-Geiger data
kp1 <- raster("Beck_KG_V1_present_0p0083.tif")  # Beck et ak 2018 data
kptab <- read.csv("kp_class_beck2018.csv",sep=";",stringsAsFactors = FALSE)

# re-project
ras.kp <- projectRaster(kp1,ras,crs=crs(shp),method="ngb")
plot(ras.kp)
plot(shp$geom,add=T)

# create matrix
mat.kp <- t(as.matrix(ras.kp))
dimnames(mat.kp) <- list(lon,lat)
saveRDS(mat.kp,file=paste0("data/matrix_popd/",ctry,"/beck_kg-cat_era5-land-aligned_",ctry,".rds"))
