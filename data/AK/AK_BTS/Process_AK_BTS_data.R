## Get AK groundfish bottom trawl survey data for 3 primary surveys
# Result is cleaned cpue (kg/km^2) by haul, with zeros included

# Note: EBS includes all species, whereas GOA and AI are only a subset
# here, we will filter the data to only include species represented in all regions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)

# start from flat files exported from AFSC database
setwd("C:/Users/lewis.barnett/Work/AFSC/Data/AK_BTS/data-raw/")
EBS = read.csv("cpue_EBSshelf_all_spp.csv", stringsAsFactors = FALSE)
GOA = read.csv("cpue_GOA_selected_spp.csv", stringsAsFactors = FALSE) # CPUE is (num or kg / km^2)
AI = read.csv("cpue_AI_selected_spp.csv", stringsAsFactors = FALSE)

# standardize columns among surveys and combine
EBS <- select(EBS, -STATIONID)
data <- rbind(EBS,GOA,AI)

# filter to species represented in all surveys/regions
#species_common = Reduce(intersect, list(GOA$SPECIES_CODE, AI$SPECIES_CODE, EBS$SPECIES_CODE))
#data <- filter(data, SPECIES_CODE %in% species_common)

# join haul data to get coordinates, depth, bottom and surface temperature
haul <- read.csv("haul.csv", stringsAsFactors = FALSE)
haul <- cbind(haul, 
             geosphere::midPoint(cbind(haul$START_LONGITUDE, haul$START_LATITUDE), 
                                     cbind(haul$END_LONGITUDE, haul$END_LATITUDE))) # get haul midpoints
haul$DATE <- as.Date(haul$START_TIME, format = "%d-%b-%y")
haul$MONTH <- lubridate::month(haul$DATE)
haul$DAY <- lubridate::day(haul$DATE)
haul <- haul %>% select(HAULJOIN, GEAR_DEPTH, SURFACE_TEMPERATURE, GEAR_TEMPERATURE, LATITUDE = lat, LONGITUDE = lon, 
                        DATE, DAY, MONTH)
data <- inner_join(data, haul)

# join species names
species_codes =  read.csv("species.csv", stringsAsFactors = FALSE)
species_codes = select(species_codes, -YEAR_ADDED)
data <- inner_join(data, species_codes)

# select and rename columns, dropping rows with mising depths
data <- data %>% select(YEAR, SURVEY, BOTTOM_DEPTH = GEAR_DEPTH, SURFACE_TEMPERATURE, GEAR_TEMPERATURE, 
                        CPUE = WGTCPUE, LATITUDE, LONGITUDE, DATE, DAY, MONTH, SPECIES_NAME, COMMON_NAME, SPECIES_CODE) %>%
                 tidyr::drop_na(BOTTOM_DEPTH, LATITUDE, LONGITUDE)

write.csv(data, "C:/Users/lewis.barnett/Work/AFSC/Data/AK_BTS/data/AK_BTS.csv")
saveRDS(data, "C:/Users/lewis.barnett/Work/AFSC/Data/AK_BTS/data/AK_BTS.rds")

# CAVEATS
# 1) species included are those requested for assessments and may not be the most data-rich
# 2) some species have been lumped/split over the years, which could cause problems (e.g., northern, southern rock sole)

# FOR PREDICTIONS
# 1) note that temperature is missing for many observations