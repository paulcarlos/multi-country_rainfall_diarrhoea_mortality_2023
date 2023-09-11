# activate libraries
library(ecmwfr)
library(lubridate)
library(sf)
library(zoo)

ctry <- "phl" # Philippines

# get extent
shp <- st_read(paste0("gadm41_",toupper(ctry),".gpkg"),"ADM_ADM_0",quiet=T)
ext <- st_bbox(shp)
ext[c("ymin","xmin")] <- floor(ext[c("ymin","xmin")]) #reduce min coordinates
ext[c("ymax","xmax")] <- ceiling(ext[c("ymax","xmax")]) #increase max coordinates

# set the key / authorisation to download from CDS
UID <- "" # add own UID from the CDS website: https://cds.climate.copernicus.eu/
wf_set_key(user=UID,
           key=wf_get_key(user=UID,service="cds"),
           service="cds")

# create repetitive variables
hrs <- paste0(sprintf("%02d",0:23),":00")
coor <- paste0(unname(ext[c("ymax","xmin","ymin","xmax")]),collapse="/") #rearrange

# loop download 
#vname <- "t2m"; var <- "2m_temperature"  # for temperature
vname <- "prcp"; var <- "total_precipitation" # for rainfall
pname <- paste0("ecmwf/era5-land/",vname,"/",ctry,"/") # path name for saving data
yrs <- 2005:2020  # select 1 year before and after of actual targeted years
mons <- sprintf("%02d",1:12) # months in two digits
for (i in yrs) {
  for (j in mons) {
    dnum <- sprintf("%02d",1:days_in_month(as.yearmon(paste0(i,"-",j)))) #create days
    fname <- paste0(vname,"_hourly_era5land_",
                    "lon_",ext["xmin"],"_",ext["xmax"],"_",
                    "lat_",ext["ymin"],"_",ext["ymax"],"_",
                    i,"-",j,"_",ctry,".nc") #filename
    rlist <- list(
      "dataset_short_name" = "reanalysis-era5-land",
      "variable" = var,
      "year" = i,
      "month" = j,
      "day" = dnum,
      "time" = hrs,
      "area" = coor,
      "format" = "netcdf",
      "target" = fname
    )
    wf_request(user=UID, request=rlist,transfer=T,path=pname)
  }
}
#rm(i,j,dnum,fname,rlist)
