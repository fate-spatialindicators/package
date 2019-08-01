library(tidyverse)
library(rgdal)
library(maps)
library(PBSmapping)
library(marmap)
library(viridis)

#  Load stat area shapefile.
statdat <- readOGR(dsn="Data",layer="adfg_stat_areas_simple") %>% 
  fortify(simple)

#-------------------------------------------------------------------
#  Plot basic example
statdat %>% 
  ggplot(aes(long,lat,group=group)) + 
  geom_polygon(fill=NA,color="black") + 
  theme_bw()

#  Determine whether longitudes are positive or negative
summary(statdat$long)

#  Load world map with only Alaska
world <- map_data("world2",region="USA") %>% 
  filter(subregion=="Alaska")

#  Plot stat areas with only Alaska
ggplot() + 
  geom_polygon(data=statdat,aes(long,lat,group=group),fill=NA,color="black") + 
  geom_polygon(data=world,aes(long-360,lat,group=group),fill="grey") + 
  theme_bw()
#-------------------------------------------------------------------


#  Created some latitude / longitude boundaries based on extent of stat areas grid
xmin <- -195
xmax <- -125
ymin <- 45
ymax <- 70

#  Load world map centered on the Pacific and convert longitudes to negative so they align with stat areas shapefile
myworld <- map_data("world2") %>% 
  mutate(long=long-360)

#  Rename world map columns to align with PBSmapping style.
names(myworld) <- c("X","Y","PID","POS","region","subregion")

#  Clip world map to lat/lon boundaries
myworld <- clipPolys(myworld, xlim=c(xmin,xmax),ylim=c(ymin,ymax), keepExtra=TRUE)

#  Now we should have a map that crosses the dateline and includes Asia. 
ggplot() + 
  geom_polygon(data=statdat,aes(long,lat,group=group),fill=NA,color="black") + 
  geom_polygon(data=myworld,
               aes(x=X,y=Y,group=factor(PID)),color="white") + 
  coord_map(xlim=c(xmin,xmax),ylim=c(ymin,ymax)) + 
  theme_bw()


#  Query NOAA bathymetry data as two objects - one on either side of the 180 line.
mymap1 <- getNOAA.bathy(lon1 = 180-(-xmin-180), lon2 = 180, lat1 = ymin, lat2 = ymax, resolution = 5) %>% 
  fortify.bathy() %>% 
  mutate(x=x-360)
mymap2 <- getNOAA.bathy(lon1 = -180, lon2 = xmax, lat1 = ymin, lat2 = ymax, resolution = 5) %>% 
  fortify.bathy()

#  Combine two bathy objects and convert to bathy object type
newmap <- bind_rows(mymap2,mymap1) %>% 
  as.bathy()

#  The "autoplot" function allows us to plot bathymetry objects within ggplot. 
#  In this case, we want to plot our bathymetry data as a nice pretty raster. 
#  We could also use geom=c("tile") which should give us contour lines instead. 
p1 <- autoplot(newmap,geom=c("raster"),coast=FALSE)+   
  geom_polygon(data=myworld,aes(x=X,y=Y,group=factor(PID))) + 
  theme_bw() + 
  scale_fill_continuous(guide = FALSE)

#  Plot bathy data with stat area grid.
png("Intro_Talk/AK_Map_with_grid_only.png",width=9.75,height=6.5,units="in",res=200)
p1 + 
  geom_polygon(data=statdat,aes(long,lat,group=group),fill=NA,color="grey") + 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0))  + 
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        plot.background = element_rect(fill="black"))
dev.off()

#  Read in fish-ticket data in order to identify the stat areas that are fished.
data <- readRDS("Data/AKFIN_comp_ft_trimmed_1991_2018.rds")

#  Create a list of the stat areas.
statareas <- (readOGR(dsn="Data",layer="adfg_stat_areas_simple"))@data
statvec <- data.frame(stat6=statareas$STAT_AREA,id=row.names(statareas))

#  Tally the number of fish ticket records in each stat area in each year.
#  Because pollock and salmon will dominate, remove them from this example
datasub <- data %>% 
  filter(SPECIES_NAME!="pollock, walleye" & FISHERY_DESCRIPTION!="SALMON") %>% 
  dplyr::select(year,stat6) %>% 
  group_by(stat6,year) %>% 
  summarise(n=n()) %>% 
  inner_join(statvec)

#  Create another version where NAs are filled with zeroes.
datasub.j <- statdat %>% 
  left_join(datasub) %>% 
  mutate(n=replace_na(n,0))

png("Intro_Talk/AK_Map_with_grid_no_bathy.png",width=9.75,height=6.5,units="in",res=200)
ggplot()+
  geom_polygon(data=myworld,aes(x=X,y=Y,group=factor(PID))) +
  theme_bw() +
  geom_polygon(data=datasub.j %>% filter(year==2018),aes(long,lat,group=group,fill=(n)),color="grey50") + 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0))  + 
  scale_fill_viridis(guide = FALSE) + 
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        plot.background = element_rect(fill="black"))
dev.off()

png("Intro_Talk/AK_Map_with_grid_and_bathy.png",width=9.75,height=6.5,units="in",res=200)
p1 +
  geom_polygon(data=datasub.j %>% filter(year==2018),aes(long,lat,group=group,fill=(n)),color="grey50") + 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0))  + 
  scale_fill_viridis(guide = FALSE) + 
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        plot.background = element_rect(fill="black"))
dev.off()