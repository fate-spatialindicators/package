# Simplified code example for model fitting and prediction with sdmTMB

#devtools::install_github("pbs-assess/sdmTMB")
#library(INLA) # if we want tools to make other meshes
library(ggplot2)
library(raster)
library(dplyr)
library(sdmTMB)
library(sp)

## if using Microsoft R Open:
setMKLthreads(parallel::detectCores(logical = FALSE) - 1)
# getMKLthreads()

#########################################################################################################
# options 

# specify # of knots for mesh (aim for ~0.002-0.003 * square kilometers of survey area)
n_knots = 750
unit_scale = 1000
# specify species to model
species = c("smooth_dogfish","spiny_dogfish", "winter_skate",
            "summer_flounder", "atlantic_cod", "american_plaice",
            "red_hake","silver_hake", "haddock",
            "cusk", "str_sea_bass", "sea_raven")[c(1,4,5)]

###########################################################################################################
# Prepare data and fit models

# Load combined catch and haul data 
data <- readRDS("data/NE/NE_BTS/NE_BTS.rds") %>% 
  dplyr::mutate(common_name = dplyr::case_when(SVSPP == "13" ~ "smooth_dogfish",
                                               SVSPP == "15" ~ "spiny_dogfish",
                                               SVSPP == "23" ~ "winter_skate",
                                               SVSPP == "103" ~ "summer_flounder",
                                               SVSPP == "73" ~ "atlantic_cod",
                                               SVSPP == "102" ~ "american_plaice",
                                               SVSPP == "77" ~ "red_hake",
                                               SVSPP == "72" ~ "silver_hake",
                                               SVSPP == "74" ~ "haddock",
                                               SVSPP == "84" ~ "cusk",
                                               SVSPP == "141" ~ "str_sea_bass",
                                               SVSPP == "163" ~ "sea_raven",
                                               TRUE ~ NA_character_)) %>% 
  dplyr::rename(LONGITUDE = LON,
                LATITUDE = LAT) %>% 
  dplyr::select(-TOW,
                -LENGTH,
                -NUMLEN,
                -CATCHSEX) %>%
  dplyr::distinct(.keep_all = TRUE)

# project to UTM
coordinates(data) <- c("LONGITUDE", "LATITUDE")
proj4string(data) <- sp::CRS("+proj=longlat +datum=WGS84")
# Geospatial calculations are executed under the North America Albers Equal Area Conic projection (SRID: 102008).
wg_crs <- sp::CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")
data <- as.data.frame(sp::spTransform(data, wg_crs))

colnames(data) = tolower(colnames(data))

# rescale coordinates
data$X <- data$longitude / unit_scale
data$Y <- data$latitude / unit_scale
year_range <- range(data$year, na.rm = TRUE)

spp = 1
sea = 1
seasons <- c("FALL", "SPRING")
# fit same model structure to each species 
for(spp in 1:length(species)) {
  for(sea in 1:length(seasons)) {
    
    ## Alternative approach uses 99% kernal density estimate
    kde <- data %>%
      dplyr::filter(common_name == species[spp]) %>%
      dplyr::select(longitude, latitude) %>%
      sp::SpatialPoints() %>%
      adehabitatHR::kernelUD() %>%
      adehabitatHR::getverticeshr(percent = 99) %>%
      sf::st_as_sf(coords = c("longitude", "latitude")) %>%
      sf::st_set_crs(wg_crs)
    
    ## taken from https://github.com/r-spatial/sf/issues/231#issuecomment-290817623
    sfc_as_cols <- function(x, names = c("x","y")) {
      stopifnot(inherits(x,"sf") && inherits(sf::st_geometry(x),"sfc_POINT"))
      ret <- sf::st_coordinates(x)
      ret <- tibble::as_tibble(ret)
      stopifnot(length(names) == ncol(ret))
      x <- x[ , !names(x) %in% names]
      ret <- setNames(ret,names)
      dplyr::bind_cols(x,ret)
    }
    
    station_dat <- data %>%
      dplyr::select(c("cruise6", "stratum", "station", "svvessel",
                      "year", "season", "est_towdate", "depth", 
                      "surftemp", "surfsalin", "bottemp", "botsalin", 
                      # "log_depth_scaled", "log_depth_scaled2", 
                      "longitude", "latitude", "X", "Y")) %>%
      dplyr::distinct(.keep_all = TRUE)
    
    ## Select svspp and join station data
    data_spp <- data %>%
      dplyr::filter(common_name == species[spp]) %>%
      right_join(station_dat, by = c("cruise6", "stratum", "station", "svvessel", 
                                     "year", "season", "est_towdate", "depth", 
                                     "surftemp", "surfsalin", "bottemp", "botsalin",
                                     "longitude", "latitude", "X", "Y")) %>%
      dplyr::filter(!is.na(depth),
                    year >= 1967,
                    season == season[sea]) %>%
      mutate(biomass = ifelse(is.na(biomass),
                              0,
                              biomass),
             abundance = ifelse(is.na(abundance) | biomass == 0,
                                0,
                                abundance),
             present = ifelse(abundance == 0,
                              0,
                              1),
             ## Successful tows were for 30 min at a speed of 3.5 knots before 2009 and 20 min at a speed of 3.0 knots since 2009,
             ## with area swept changing from 0.038 to 0.024 km2 (Politis et al., 2014; ASMFC, 2015; Li et al., 2018).
             density = case_when(abundance == 0 ~ 0,
                                 abundance != 0 & year >= 2009 ~ abundance/(0.024),
                                 abundance != 0 & year < 2009 ~ abundance/(0.038),
                                 TRUE ~ NA_real_),
             log_depth_scaled = scale(log(depth), scale = TRUE, center = TRUE),
             log_depth_scaled2 = log_depth_scaled^2)
    
    # create one mesh for all species, without subsetting by the geographic range of positive observations
    c_spde <- make_spde(data_spp$X, data_spp$Y, n_knots = n_knots) 
    
    # # filter by species and only include range of coordinates with positive observations over the timeseries
    # data_sub <-  data_spp %>%
    #   # dplyr::filter(common_name == species[spp]) %>%
    #   dplyr::filter(latitude >= min(latitude[which(density > 0)]),
    #                 latitude <= max(latitude[which(density > 0)]),
    #                 longitude >= min(longitude[which(density > 0)]),
    #                 longitude <= max(longitude[which(density > 0)])) %>% 
    #   dplyr::select(year,
    #                 X,
    #                 Y,
    #                 longitude,
    #                 latitude,
    #                 log_depth_scaled,
    #                 log_depth_scaled2,
    #                 depth,
    #                 density,
    #                 biomass,
    #                 abundance,
    #                 present)
    
    data_sub <- data_spp %>% 
      sf::st_as_sf(coords = c("longitude", "latitude")) %>%
      sf::st_set_crs(wg_crs) %>%
      sf::st_intersection(kde) %>%
      sfc_as_cols(names = c("longitude", "latitude")) %>%
      sf::st_drop_geometry() %>%
      dplyr::select(-id,
                    -area) %>%
      mutate(year = as.integer(year)) %>%
      arrange(year) %>% 
      select(year, 
             X, Y, 
             depth,
             density,
             log_depth_scaled,
             log_depth_scaled2,
             longitude, 
             latitude)
    
    # if subset survey area by species, then make new mesh for each species
    #c_spde <- make_spde(data_sub$X, data_sub$Y, n_knots = n_knots) 
    # plot_spde(c_spde)

    # start_time <- Sys.time()
    density_model <- sdmTMB(formula = density ~ 0 + as.factor(year),
                            time_varying = ~ 0 + log_depth_scaled + log_depth_scaled2,
                            data = data_sub,
                            time = "year", 
                            spde = c_spde, 
                            reml = TRUE,
                            anisotropy = TRUE,
                            silent = FALSE,
                            family = tweedie(link = "log"),
                            control = sdmTMBcontrol(step.min = 0.01, step.max = 1))
    # end_time <- Sys.time()

    saveRDS(density_model, 
            file = sprintf("output/NE/%s_%s_density_depth_varying.rds", species[spp], seasons[sea]))
    
    # start_time <- Sys.time()
    density_model_2 <- sdmTMB(formula = density ~ 0 + as.factor(year),
                            data = data_sub,
                            time = "year", 
                            spde = c_spde, 
                            reml = TRUE,
                            anisotropy = TRUE,
                            silent = FALSE,
                            family = tweedie(link = "log"),
                            control = sdmTMBcontrol(step.min = 0.01, step.max = 1))
    # end_time <- Sys.time()
    
    saveRDS(density_model2, 
            file = sprintf("output/NE/%s_%s_density_no_covar.rds", species[spp], seasons[sea]))
    # 
    # tt <- difftime(end_time, start_time)
    # time_difference <- sprintf("%s %s", round(tt[[1]], 0), attr(tt, "units"))
    # 
    # 
    # glue::glue("Hi!
    #             It looks like {tolower(seasons[sea])} {gsub('_', ' ', species[spp])} has finished
    #             running and has been saved. This model took about
    #             {time_difference} to run.
    #             ")

  }
}

###########################################################################################################
# Prepare prediction data frame (with columns for time, coordinates, covariates)
# 
## create prediction grid
neus_grid <- read.csv(here::here("data/NE/nes_lon_lat.csv"))
xlims <- c(-77, -65)
ylims <- c(35, 45)

# 200 m isobath layer
nesbath <- marmap::getNOAA.bathy(lon1 = xlims[1], lon2 = xlims[2],
                                 lat1 = ylims[1], lat2 = ylims[2],
                                 resolution = 6,
                                 keep = TRUE) %>%
  marmap::as.raster()

neus_sp <- neus_grid %>%
  dplyr::mutate(z = raster::extract(nesbath, .) * -1) %>%
  dplyr::filter(z > 0) %>% 
  dplyr::rename(longitude = x,
                latitude = y)
# # project to UTM
sp::coordinates(neus_sp) <- c("longitude", "latitude")
sp::proj4string(neus_sp) <- sp::CRS("+proj=longlat +datum=WGS84")
neus_sp <- as.data.frame(sp::spTransform(neus_sp, wg_crs))

# clip to NEFSC survey area
newdat <- neus_sp %>%
  dplyr::mutate(log_depth_scaled = base::scale(log(z), scale = TRUE, center = TRUE)[,1],
                log_depth_scaled2 = log_depth_scaled^2,
                depth = z,
                X = longitude / unit_scale,
                Y = latitude / unit_scale) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(year = list(1967:2015)) %>%
  tidyr::unnest(year) %>%
  dplyr::select(year,
                X, Y,
                depth,
                log_depth_scaled, 
                log_depth_scaled2,
                longitude,
                latitude)

saveRDS(newdat, file = paste0("data/NE/NE_BTS/predict_data.rds")) # save prediction grid

# ggplot() +
#   geom_point(data = newdat %>% filter(year %in% c(2000:2015)), aes(x = X, y = Y, color = depth)) +
#   geom_point(data = data_spp %>% filter(year %in% c(2000:2015)), aes(x = X, y = Y), color = "white") +
#   facet_wrap(~year)
# 
# 
# newdat_kde <- newdat %>%
#   sf::st_as_sf(coords = c("longitude", "latitude")) %>%
#   sf::st_set_crs(wg_crs) %>%
#   sf::st_intersection(kde) %>%
#   sfc_as_cols(names = c("longitude", "latitude")) %>%
#   sf::st_drop_geometry() %>%
#   dplyr::select(-id,
#                 -area) %>%
#   mutate(year = as.integer(year)) %>%
#   arrange(year)




# # Starting with prediction grid from the broader region
# Predict_data <- readRDS(paste0(here::here(), "/data/AK/AK_BTS/Predict_data_AK.Rds"))
# Predict_data <- Predict_data %>% filter(GOA == 1) %>% select(LONG, LAT, BOTTOM_DEPTH)
# 
# # make erroneous depth estimates nearshore equal to minimum depth
# Predict_data$BOTTOM_DEPTH <- ifelse(Predict_data$BOTTOM_DEPTH < min(data$bottom_depth), 
#                                     min(data$bottom_depth),
#                                     Predict_data$BOTTOM_DEPTH)
# 
# # make prediction data frame with same covariates as went into the model fit
# Predict_data$log_depth_scaled = scale(log(Predict_data$BOTTOM_DEPTH))
# Predict_data$log_depth_scaled2 = Predict_data$log_depth_scaled ^ 2
# Predict_data$X = Predict_data$LONG / 10000
# Predict_data$Y = Predict_data$LAT / 10000
# Predict_data = select(Predict_data, X, Y, log_depth_scaled, log_depth_scaled2)
# #Predict_data %>% rename(log_depth_scaled = log_depth_scaled.V1, log_depth_scaled2 = log_depth_scaled2.V1)
# 
# # make replicates for each year in data
# Predict_data_years = Predict_data
# Predict_data_years$year = min(unique(data$year))
# for(i in 2:length(unique(data$year))) {
#   Predict_data$year = sort(unique(data$year))[i]
#   Predict_data_years = rbind(Predict_data_years, Predict_data)
# }
# 
# saveRDS(Predict_data_years, file=paste0("data/AK/AK_BTS/GOA_predict_data.rds")) # save prediction grid

###########################################################################################################
# Make plots from predictions, model fit (you can wrap this in a loop to perform by species)

# predict from model to full prediction domain (space and time)
# p = predict(density_model, newdata = newdat, se_fit = FALSE)
# 
# # # plotting functions
# plot_map_raster <- function(dat, column = "omega_s") {
# ggplot(p %>% filter(year %in% c(2000:2015)), aes_string("X", "Y", color = "est")) +
#   geom_point() +
#   scale_color_viridis_c() +
#   xlab("Longitude") +
#   ylab("Latitude")
# }
# # 
# # # make plot of predictions from full model (all fixed + random effects)
# plot_map_raster(p, "est") +
#  facet_wrap(~year) +
#  coord_fixed()
