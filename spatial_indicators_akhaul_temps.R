#-------------------------------------------------------------------------------------------------------------------------
# This script pulls JPL / ERDDAP satellite data of SST for the coordinates of the Alaska Bottom Trawl survey. It adds a spatial buffer around each point and averages the sst data for that area and date.
# Author: Jordan Watson (jordan.watson@noaa.gov)
# Date created: 07/19/2019
#-------------------------------------------------------------------------------------------------------------------------

# Load packages
library(tidyverse)
library(ncdf4)
library(RCurl)
library(lubridate)

#  Read in the haul data, find distinct sets of date-coordinates, add a year column
haul <- readRDS("Data/AK_BTS.rds")
hauldat <- haul %>% 
  dplyr::select(YEAR,LATITUDE,LONGITUDE,DATE) %>% 
  distinct() %>% 
  rename_all(tolower) %>% 
  mutate(year=lubridate::year(lubridate::parse_date_time(date,orders="ymd")))
  
#  Create a function for querying different lat-long ranges. The script will create a polygon +/- a buffer (0.01 degrees is the default) around each point. 
mysstgridfun <- function(mydate,mylat,mylon){
  #mylat=38.08 #dummy value for troubleshooting
  #mylon=-153.12 #dummy value for troubleshooting
  #mydate="2005-06-03" #dummy value for troubleshooting

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
    summarise(sstsd=round(sd(sst),3),
              sstn=length(sst),
              sst=round(mean(sst),2),
              date=mydate,
              lat=mylat,
              lon=mylon)
  
  nc_close(nc)
  
  return(tempneg)
}

# The JPL MUR dataset starts on 2002-06-01 so lots of the AK haul data were before this time. 
# Createa a version of the hauldata that contains only dates since the MUR dataset came online.
hauldatnew <- hauldat %>% filter(year>"2002-06-02")

#  Run a loop that pulls the satellite data for each day and set of lat-lon boundaries
for(i in 1:nrow(hauldatnew)){
    tryCatch({ # Don't let error messages stop the loop.
      print(i)
      
      #  Run the function that we created above.
      tempdat <- mysstgridfun(hauldatnew$date[i],hauldatnew$lat[i],hauldatnew$lon[i])
      
      #  Bind all dates together.
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

#  Save this new dataset.
saveRDS(basket,file="Data/hauldata_temperatures_AK_new.RDS")


#-----------------------------------------------------------
#  Pull temperature data for older survey data (prior to 2002-06-02) using
#  https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdPH2ssta1day.html
#  This dataset is at a much coarser spatial scale (about 0.04 degrees instead of 0.01)
#-----------------------------------------------------------

#  Create an older version of the data using this sst dataset back to 1981
#  This dataset discontinues in 2012. I ran it to include all the dates for which there is 
#  overlap with the JPL MUR dataset. So for the years from 2002-2012, we will have two sets of temperature data to compare.
hauldatold <- haul %>% 
  dplyr::select(YEAR,LATITUDE,LONGITUDE,DATE) %>% 
  distinct() %>% 
  rename_all(tolower) %>% 
  mutate(year=lubridate::year(lubridate::parse_date_time(date,orders="ymd"))) %>% 
  filter(date>"1981-10-31" & date<"2012-12-31")


#  Create a function for querying different lat-long ranges. Because this dataset is at such a coarser spatial scale, I did not 
#  put a spatial buffer around each point - by default there is already a much greater spatial buffer. So I have the buffer set to 0. 
mysstgridfun <- function(mydate,mylat,mylon){
  
  #mylat=59.66785 #dummy value for troubleshooting
  #mylon=-177.1424 #dummy value for troubleshooting
  #mydate="2012-07-25" #dummy value for troubleshooting
  
  #  Establish the spatial buffer around points
  buffer=0.2
  minlat=mylat-buffer
  maxlat=mylat+buffer
  minlon=mylon-buffer
  maxlon=mylon+buffer
  
  #  Download SST data from jplMURSST41 as netcdf
  x <- getBinaryURL(paste0("https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdPH2ssta1day.nc?sea_surface_temperature[(",mydate,"T12:00:00Z):1:(",mydate,"T12:00:00Z)][(",maxlat,"):1:(",minlat,")][(",minlon,"):1:(",maxlon,")]"))

  #  Convert and open the netcdf file
  tmpSST <- tempfile(pattern="xwB", fileext=".nc")
  writeBin(object=x, con=tmpSST)
  nc <- nc_open(tmpSST)
  
  #  Extract and average the temperature records.
  tempneg <- expand.grid(lon=ncvar_get(nc, varid = "longitude"),
                         lat=ncvar_get(nc, varid = "latitude")) %>% 
                         {bind_cols(.,data.frame(sst=as.vector(ncvar_get(nc,"sea_surface_temperature"))))} %>% 
    mutate(sst=ifelse(sst<(-2),-2,sst)) %>% 
    summarise(sst=round(mean(sst),2),
              date=mydate,
              lat=mylat,
              lon=mylon)
  
  nc_close(nc)
  
  return(tempneg)
}




#  This takes a few hours to run on almost 25000 data records.
for(i in 1:nrow(hauldatold)){
  tryCatch({ # Don't let error messages stop the loop.
    print(i)
    
    #  Use the above function to query three separate lat/lon rectangles
    tempdat <- mysstgridfun(hauldatold$date[i],hauldatold$lat[i],hauldatold$lon[i])
    
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


saveRDS(basket,file="Data/hauldata_temperatures_AK_old.RDS")

new <- readRDS("Data/hauldata_temperatures_AK_new.RDS") %>% 
  rename(sstMUR=sst)
old <- readRDS("Data/hauldata_temperatures_AK_old.RDS") %>% 
  rename(sstPH2=sst)

sstdata <- full_join(new,old)

minidat <- tail(hauldatold %>% arrange(date),100)

#  This takes a few hours to run on almost 25000 data records.
for(i in 1:nrow(minidat)){
  tryCatch({ # Don't let error messages stop the loop.
    print(i)
    
    #  Use the above function to query three separate lat/lon rectangles
    tempdat <- mysstgridfun(minidat$date[i],minidat$lat[i],minidat$lon[i])
    
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


basket05 <- basket
basket1 <- basket
basket2 <- basket





#  Join the environmental data back onto the haul dataset.
saveRDS(haul %>% 
          left_join(basket),
        file="Data/nwfsc_haul_with_sst.RDS")

test %>% 
  filter(!is.na(temperature_at_surface_c_der)) %>% 
  mutate(diff=sst-temperature_at_surface_c_der) %>% 
  ggplot(aes(sst,diff)) + 
  geom_point() + 
  facet_wrap(~year)

test %>% 
  filter(!is.na(temperature_at_surface_c_der)) %>% 
  mutate(diff=sst-temperature_at_surface_c_der) %>% 
  ggplot(aes(sst,diff)) + 
  geom_point() + 
  facet_wrap(~year+vessel)