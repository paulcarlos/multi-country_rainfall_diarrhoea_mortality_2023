library(sf)
library(ncdf4)
library(raster)

ctry <- "phl"  #phl=Philippines
vnam <- "t2m"  #t2m=temperature
lvl <- paste0(c("adm","ADM_"),1)  #administrative level

# load shapefile
shp <- st_read(paste0("gadm41_",toupper(ctry),".gpkg"), layer=paste0("ADM_",lvl[2]),quiet=T)

# load sample ERA5-Land
nc <- nc_open(list.files(paste0("ecmwf/era5-land/",vnam,"/",ctry,"/"),full.names=T)[1])
lon <- as.vector(ncvar_get(nc,varid="longitude"))
lat <- as.vector(ncvar_get(nc,varid="latitude"))
mat <- ncvar_get(nc,varid=vnam,start=c(1,1,1),count=c(length(lon),length(lat),1))  # matrix for 1 hour data
dimnames(mat) <- list(lon,lat)
mat1 <- t(mat)
ras <- raster(mat1,xmx=max(lon),xmn=min(lon),ymx=max(lat),ymn=min(lat),crs=crs(shp))
#plot(ras); plot(shp$geom,add=T)
nc_close(nc); rm(nc)

# select grids using extract convert into full matrix with 1s and 0s
bmat <- matrix(0,nrow=dim(mat1)[1],ncol=dim(mat1)[2],dimnames=dimnames(mat1)) #transposed to fit the rasterised version
mlist <- list() #blank list
for (i in 1:nrow(shp)) { #loop per district
  cat(i," ")
  ext <- extract(ras,shp[i,],df=TRUE,cellnumbers=TRUE,weights=TRUE)
  xy <- data.frame("rowras"=rowFromCell(ras,cell=ext$cell),"colras"=colFromCell(ras,cell=ext$cell))
  smat <- bmat
  for (j in 1:nrow(xy)) {
    smat[xy$rowras[j],xy$colras[j]] <- 1
  }
  mlist[[shp$GID_1[i]]] <- t(smat)
}
rm(i,j,ext,xy,smat)
saveRDS(mlist,paste0("data/matrix_loc/",ctry,"/matrix_",ctry,"_grids_era5land_",lvl[1],".rds"))
