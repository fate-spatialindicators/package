#-------------------------------------------------------------------------------------------------------------------------
# This script pulls JPL / ERDDAP satellite data of SST for the coordinates of the west coast trawl survey. It adds a spatial buffer around each point and averages the sst data for that area and date.
# Author: Jordan Watson (jordan.watson@noaa.gov)
# Date created: 07/19/2019
#-------------------------------------------------------------------------------------------------------------------------

# Load packages
library(tidyverse)
library(ncdf4)
library(RCurl)

devtools::install_github("nwfsc-assess/nwfscSurvey")

library(nwfscSurvey)
#catch = PullCatch.fn(SurveyName = "NWFSC.Combo")
#saveRDS(catch, "nwfsc_catch.rds")
haul = PullHaul.fn(SurveyName="NWFSC.Combo")
saveRDS(haul, "nwfsc_haul.rds")

haul <- readRDS("Data/nwfsc_haul.rds") %>% 
  mutate(year=lubridate::year(lubridate::parse_date_time(date_yyyymmdd,orders="ymd")),
         date=lubridate::parse_date_time(date_yyyymmdd,orders="$Y-%m-%d"),
         lat=round(latitude_dd,2),
         lon=round(longitude_dd,2))

#ggplot2::ggplot(haul, aes(longitude_dd, latitude_dd,col=temperature_at_gear_c_der)) + 
#  geom_point(size=0.1) + facet_wrap(~year)

hauldat <- haul %>% 
  distinct(date,lat,lon)


#  Create a function for querying different lat-long ranges. The script will create a polygon +/- 0.01 degrees around each point. 
mysstgridfun <- function(mydate,mylat,mylon){
  #mylat=38.08
  #mylon=-123.12
  #mydate="2003-07-29"
  #  Establish the spatial buffer around points
  buffer=0.01
  minlat=mylat-buffer
  maxlat=mylat+buffer
  minlon=mylon-buffer
  maxlon=mylon+buffer
  
  #  Download SST data from jplMURSST41 as netcdf
  x <- getBinaryURL(paste0("https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41.nc?analysed_sst[(",mydate,"T09:00:00Z):1:(",mydate,"T09:00:00Z)][(",minlat,"):1:(",maxlat,")][(",minlon,"):1:(",maxlon,")]"))
  
  #  Convert and open the netcdf file
  tmpSST <- tempfile(pattern="xwB", fileext=".nc")
  writeBin(object=x, con=tmpSST)
  nc <- nc_open(tmpSST)
  
  #  Extract and average the temperature records.
  tempneg <- expand.grid(lon=ncvar_get(nc, varid = "longitude"),
                         lat=ncvar_get(nc, varid = "latitude")) %>% 
                         {bind_cols(.,data.frame(sst=as.vector(ncvar_get(nc,"analysed_sst"))))} %>% 
    mutate(sst=ifelse(sst<(-2),-2,sst)) %>% 
    summarise(sst=mean(sst),
              date=mydate,
              lat=mylat,
              lon=mylon)
  
  nc_close(nc)
  
  return(tempneg)
}


system.time({
  for(i in 1:nrow(hauldat)){
    tryCatch({ # Don't let error messages stop the loop.
      print(i)
      
      #  Use the above function to query three separate lat/lon rectangles
      tempdat <- mysstgridfun(hauldat$date[i],hauldat$lat[i],hauldat$lon[i])
      
      #  Rename this data frame, which will become important for binding in the next iteration of the loop
      if(i==1){
        basket <- tempdat
      } else {
        # Bind successive years
        basket <- basket %>% 
          bind_rows(tempdat)
      }
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
})


saveRDS(basket,file="Data/hauldata_temperatures.RDS")

#  Join the environmental data back onto the haul dataset.
haul %>% 
  left_join(basket)
