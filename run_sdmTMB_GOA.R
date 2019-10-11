# Simplified code example for model fitting and prediction with sdmTMB

#devtools::install_github("pbs-assess/sdmTMB")
#library(INLA) # if we want tools to make other meshes
library(ggplot2)
library(raster)
library(dplyr)
library(sdmTMB)
library(sp)

###########################################################################################################
# options 

# specify # of knots for mesh
n_knots = 450

# specify species to model
species = c("Dover sole","arrowtooth flounder", "Pacific halibut",
            "walleye pollock", "rex sole", "English sole","sablefish","Pacific cod",
            "spiny dogfish","longnose skate","big skate", "Pacific ocean perch")

###########################################################################################################
# Prepare data and fit models

# Load combined catch and haul data 
data <- readRDS(paste0(here::here(),"/data/AK/AK_BTS/AK_BTS.rds"))

# filter to GOA survey, remove tows with 0 bottom depth, and drop 2001 year when survey was incomplete
data <- data %>% filter(SURVEY == "GOA", BOTTOM_DEPTH > 0, YEAR != 2001)

# read in the grid cell data from the survey design (one may choose to pre-specify which hauls are in which cells)
#grid_cells = read.csv(paste0(here::here(),"/data/AK/AK_BTS/survey_grids/grid_GOA.csv"))

# project to UTM
coordinates(data) <- c("LONGITUDE", "LATITUDE")
proj4string(data) <- CRS("+proj=longlat +datum=WGS84")
data <- as.data.frame(spTransform(data, CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")))

colnames(data) = tolower(colnames(data))

# mesh for GOA, which you can pass to make_spde or let it create a default mesh as below
#meshbuilder() # or play with this shiny interactive tool to make a good mesh to pass to make_spde
#Coord <- cbind(data$X, data$Y)	
#bnd = inla.nonconvex.hull(Coord, convex=200000)		# boundary of the region of interest
#mesh = inla.mesh.2d(loc=Coord, boundary = bnd, max.edge=c(3000000,5000000), cutoff = 50000, offset=c(100000,100000))

# rescale coordinates
data$X <- data$longitude / 10000
data$Y <- data$latitude / 10000

# center year
data$year_centered = data$year - mean(unique(data$year)) # set intercept to mean, or other value to roughly center year

# fit same model structure to each species 
for(spp in 1:length(species)) {
  data_sub = data %>% dplyr::filter(common_name == species[spp]) %>% 
    dplyr::filter(latitude > min(latitude[which(cpue>0)]),
                  latitude <= max(latitude[which(cpue>0)]),
                  longitude > min(longitude[which(cpue>0)]),
                  longitude < max(longitude[which(cpue>0)]))

  c_spde <- make_spde(data_sub$X, data_sub$Y, n_knots = n_knots) 
  plot_spde(c_spde)
  
  # center and scale depth
  data_sub$log_depth_scaled = scale(log(data_sub$bottom_depth))[,1]
  data_sub$log_depth_scaled2 = data_sub$log_depth_scaled ^ 2
  
      density_model <- sdmTMB(formula = cpue ~ log_depth_scaled + log_depth_scaled2 + year_centered,
                              data = data_sub,
                              time = "year", spde = c_spde, anisotropy = TRUE,
                              silent = TRUE, spatial_trend = TRUE, spatial_only = FALSE,
                              family = tweedie(link = "log")
                              #control = sdmTMBcontrol(step.min = 0.01, step.max = 1)
                              )
      
      saveRDS(density_model, file=paste0("output/AK/",species[spp],"_density.rds"))
}

###########################################################################################################
# Prepare prediction data frame (with columns for time, coordinates, covariates)

# Starting with prediction grid from the broader region
Predict_data <- readRDS(paste0(here::here(), "/data/AK/AK_BTS/Predict_data_AK.Rds"))
Predict_data <- Predict_data %>% filter(GOA == 1) %>% select(LONG, LAT, BOTTOM_DEPTH)

# make erroneous depth estimates nearshore equal to minimum depth
Predict_data$BOTTOM_DEPTH <- ifelse(Predict_data$BOTTOM_DEPTH < min(data$BOTTOM_DEPTH), 
                                    min(data$bottom_depth),
                                    Predict_data$BOTTOM_DEPTH)

# make prediction data frame with same covariates as went into the model fit
Predict_data$log_depth_scaled = scale(log(Predict_data$BOTTOM_DEPTH))
Predict_data$log_depth_scaled2 = Predict_data$log_depth_scaled ^ 2
Predict_data$X = Predict_data$LONG / 10000
Predict_data$Y = Predict_data$LAT / 10000
Predict_data = select(Predict_data, X, Y, log_depth_scaled, log_depth_scaled2)
#Predict_data %>% rename(log_depth_scaled = log_depth_scaled.V1, log_depth_scaled2 = log_depth_scaled2.V1)

# make replicates for each year in data
Predict_data_years = Predict_data
Predict_data_years$year = min(unique(data$year))
for(i in 2:length(unique(data$year))) {
  Predict_data$year = sort(unique(data$year))[i]
  Predict_data$year_centered = sort(unique(data$year))[i] - mean(unique(data$year))
  Predict_data_years = rbind(Predict_data_years, Predict_data)
}

saveRDS(Predict_data_years, file=paste0("data/AK/AK_BTS/GOA_predict_data.rds")) # save prediction grid

###########################################################################################################
# Make plots from predictions, model fit (you can wrap this in a loop to perform by species)

# predict from model to full prediction domain (space and time)
p = predict(density_model, newdata=Predict_data_years, se_fit = FALSE)

# plotting functions
plot_map_raster <- function(dat, column = "omega_s") {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    scale_fill_viridis_c() +
    xlab("Longitude") +
    ylab("Latitude")
}

# make plot of predictions from full model (all fixed + random effects)
plot_map_raster(p, "est") +
  facet_wrap(~year) +
  coord_fixed()
