# Notes from Lewis:
#
# The first step would be to fit spatiotemporal models to the BC trawl survey
# data for these species: "Dover sole", "arrowtooth flounder", "rex sole",
# "English sole","sablefish","Pacific cod", "spiny dogfish","longnose
# skate","big skate", "Pacific ocean perch" (along with P halibut if observed
# often enough; I am assuming you guys don't get enough pollock to include
# them).  As for the modeling approach, I started by just fitting models to the
# survey data with only depth, depth^2 and year as fixed effects, as shown in
# this stripped-down script I made to help Scott get started fitting models to
# the NE data (see lines 69-73 for arguments to sdmTMB):
# https://github.com/fate-spatialindicators/spatialindicators/blob/master/run_sdmTMB_GOA.R
#
# The second step would be to use the predictions from the above models to make
# plots of COG and inertia in a similar way as Jordan did for the commercial
# data.  See some notes on this in the messages below, but I can link to some
# cleaned up plotting code on github. We should
# also extract the COG estimates directly from sdmTMB for comparison.
#
# Lastly, if it is not too much work, computing the same indicators on the
# commercial data for the same set of species would be great.

library(dplyr)
library(sdmTMB)
library(future)
library(ggplot2)
library(purrr)

plan(multisession, workers = floor(availableCores() / 2))
# plan(sequential)
dir.create("data/BC/", showWarnings = FALSE)
dir.create("data/BC/generated/", showWarnings = FALSE)
dir.create("figures/BC/", showWarnings = FALSE)

spp <- tolower(c(
  "Dover sole", "arrowtooth flounder", "rex sole",
  "English sole", "sablefish", "Pacific cod", "north pacific spiny dogfish",
  "longnose skate", "big skate", "Pacific ocean perch", "pacific halibut"
))

.file <- "data/BC/bc-synoptic-tows.rds"
if (!file.exists(.file)) {
  d <- list()
  for (i in seq_along(spp)) {
    .spp <- gsub(" ", "-", gsub("\\/", "-", tolower(spp[i])))
    # d[[i]] <- gfdata::get_survey_sets(.spp, ssid = c(1, 3, 4, 16),
    #   join_sample_ids = TRUE)
    d[[i]] <- readRDS(paste0("../../gfs/report/data-cache/", .spp, ".rds"))
    d[[i]] <- d[[i]]$survey_sets
  }
  d <- dplyr::bind_rows(d)
  d <- dplyr::filter(d, survey_series_id %in% c(1, 3, 4, 16)) %>%
    select(year,
      survey = survey_abbrev, species = species_common_name, longitude,
      latitude, depth_m, density_kgpm2
    )
  saveRDS(d, file = .file)
} else {
  d <- readRDS(.file)
}

.file <- "data/BC/bc-synoptic-grids.rds"
if (!file.exists(.file)) {
  syn_grid <- gfplot::synoptic_grid %>% # utm 9
    select(-cell_area, -survey_domain_year)
  saveRDS(syn_grid, file = .file)
} else {
  syn_grid <- readRDS(.file)
}

convert2utm <- function(x, utm_zone = 9) {
  x <- dplyr::rename(x, X = longitude, Y = latitude)
  attr(x, "projection") <- "LL"
  attr(x, "zone") <- utm_zone
  suppressMessages(PBSmapping::convUL(x))
}
dat <- convert2utm(d)

expand_years <- function(grid, tow_dat) {
  original_time <- sort(unique(tow_dat$year))
  nd <- do.call(
    "rbind",
    replicate(length(original_time), grid, simplify = FALSE)
  )
  nd[["year"]] <- rep(original_time, each = nrow(grid))
  nd
}

scale_dat <- function(grid, tow_dat) {
  tow_dat <- dplyr::filter(tow_dat, !is.na(depth_m))
  .mean <- mean(log(tow_dat$depth_m), na.rm = TRUE)
  .sd <- sd(log(tow_dat$depth_m), na.rm = TRUE)
  grid$log_depth_scaled <- (log(grid$depth) - .mean) / .sd
  tow_dat$log_depth_scaled <- (log(tow_dat$depth_m) - .mean) / .sd
  grid$log_depth_scaled2 <- grid$log_depth_scaled^2
  tow_dat$log_depth_scaled2 <- tow_dat$log_depth_scaled^2
  list(grid = grid, tow_dat = tow_dat)
}

fit_bc_survey <- function(.survey, .species, knots = 200L, silent = TRUE, ...) {
  sp_dat <- filter(dat, species == .species, survey %in% .survey) %>%
    mutate(density_kgp100m2 = density_kgpm2 * 100 * 100) # computational stability?
  survey_grid <- filter(syn_grid, survey %in% .survey) %>%
    expand_years(sp_dat) %>%
    select(year, X, Y, depth)
  if (nrow(sp_dat) > 1L) {
    scaled <- scale_dat(grid = survey_grid, tow_dat = sp_dat)
    sp_dat <- scaled$tow_dat
    survey_grid <- scaled$grid
    spde <- make_spde(sp_dat$X, sp_dat$Y, n_knots = knots)
    density_model <- tryCatch(sdmTMB(
      formula = density_kgpm2 ~ 0 + as.factor(year),
      time_varying = ~ 0 + log_depth_scaled + log_depth_scaled2,
      data = sp_dat,
      time = "year",
      spde = spde,
      reml = TRUE,
      anisotropy = TRUE,
      family = tweedie(link = "log"),
      silent = silent
    ), error = function(e) NA)
  } else {
    density_model <- NA
  }
  list(model = density_model, grid = survey_grid, tow_dat = sp_dat)
}

to_fit <- tidyr::expand_grid(
  .survey = unique(syn_grid$survey),
  .species = unique(dat$species)
)
to_fit <- filter(to_fit, !.survey %in% "SYN WCHG")
to_fit$.group <- to_fit$.survey
to_fit <- mutate(to_fit, .group =
    ifelse(.survey %in% c("SYN QCS", "SYN HS"), "SYN QCS SYN HS", .group))
sp <- split(to_fit$.species, paste(to_fit$.group, to_fit$.species)) %>% map(1)
sv <- split(to_fit$.survey, paste(to_fit$.group, to_fit$.species))
kn <- map(sv, ~ ifelse(.x == c("SYN QCS", "SYN HS"), 200, 125)[1])

.f <- "data/BC/generated/bc-survey-fits2.rds"
if (!file.exists(.f)) {
  # bc_fits <- pmap(list(sv, sp, kn), fit_bc_survey)
  bc_fits <- furrr::future_pmap(list(sv, sp, kn), fit_bc_survey)
  saveRDS(bc_fits, file = .f)
} else {
  bc_fits <- readRDS(.f)
}

.f1 <- "data/BC/generated/bc-survey-predictions2.rds"
.f2 <- "data/BC/generated/bc-survey-cogs2.rds"
if (!file.exists(.f1) || !file.exists(.f2)) {
  predictions <- furrr::future_map(bc_fits, function(.x) {
    tryCatch(predict(.x$model, newdata = .x$grid, return_tmb_object = TRUE),
      error = function(e) NA
    )
  })
  saveRDS(predictions, file = .f1)
  cogs <- furrr::future_map(predictions, function(.x) {
    tryCatch(get_cog(.x, bias_correct = FALSE), error = function(e) NA)
  })
  saveRDS(cogs, file = .f2)
} else {
  predictions <- readRDS(.f1)
  cogs <- readRDS(.f2)
}

prediction_df <- purrr::map_df(seq_along(sv), function(i) {
  .dat <- if (!is.na(predictions[[i]])[[1]]) predictions[[i]]$data else NA
  data.frame(
    survey = paste(sv[[i]], collapse = ", "), species = sp[[i]], .dat,
    stringsAsFactors = FALSE
  )
}) %>%
  select(-density_kgpm2)
saveRDS(prediction_df, file = "data/BC/generated/bc-survey-predictions2-df.rds")

cogs_df <- purrr::map_df(seq_along(sv), ~ data.frame(
  survey = paste(sv[[.]], collapse = ", "),
  species = sp[[.]],
  cogs[[.]], stringsAsFactors = FALSE
)) %>%
  as_tibble()

stopifnot(nrow(filter(cogs_df, is.na(est))) == 0)

g <- ggplot(filter(cogs_df, coord == "X")) +
  geom_ribbon(aes(x = year, ymin = lwr, ymax = upr), alpha = 0.4) +
  geom_line(aes(x = year, y = est)) +
  facet_grid(cols = vars(species), rows = vars(survey), scales = "free_y") +
  ggsidekick::theme_sleek()
ggsave("figures/BC/cog-surveys-X.pdf", width = 20, height = 8)

g <- ggplot(filter(cogs_df, coord == "Y")) +
  geom_ribbon(aes(x = year, ymin = lwr, ymax = upr), alpha = 0.4) +
  geom_line(aes(x = year, y = est)) +
  facet_grid(cols = vars(species), rows = vars(survey), scales = "free_y") +
  ggsidekick::theme_sleek()
ggsave("figures/BC/cog-surveys-Y.pdf", width = 20, height = 8)

cogs_wide_est <- reshape2::dcast(cogs_df, survey + species + year ~ coord, value.var = "est")
cogs_wide_lwr <- reshape2::dcast(cogs_df, survey + species + year ~ coord, value.var = "lwr") %>%
  mutate(X_lwr = X, Y_lwr = Y)
cogs_wide_upr <- reshape2::dcast(cogs_df, survey + species + year ~ coord, value.var = "upr") %>%
  mutate(X_upr = X, Y_upr = Y)
cogs_wide <- cogs_wide_est %>% 
  left_join(select(cogs_wide_lwr, survey, species, year, X_lwr, Y_lwr)) %>%
  left_join(select(cogs_wide_upr, survey, species, year, X_upr, Y_upr))

g <- ggplot(cogs_wide, aes(X, Y, color = year)) + geom_path(alpha = 0.7) +
  geom_point() +
  geom_segment(aes(x = X_lwr, xend = X_upr, y = Y, yend = Y), lwd = 0.2) +
  geom_segment(aes(x = X, xend = X, y = Y_lwr, yend = Y_upr), lwd = 0.2) +
  facet_wrap(survey ~ species, ncol = 11, scales = "free") +
  scale_color_viridis_c() +
  ggsidekick::theme_sleek()
ggsave("figures/BC/cog-surveys-XY.pdf", width = 19, height = 8)

cogs_wide <- mutate(cogs_wide, diameter_x = X_upr - X_lwr, diameter_y = Y_upr - Y_lwr)
g <- ggplot(cogs_wide, aes(X, Y, color = year, fill = year)) + 
  geom_path() +
  ggforce::geom_ellipse(aes(x0 = X, y0 = Y, a = diameter_x/2, b = diameter_y/2, angle = 0), alpha = 0.1) +
  geom_point() +
  facet_wrap(survey ~ species, ncol = 11, scales = "free") +
  scale_color_viridis_c(option = "C") +
  scale_fill_viridis_c(option = "C") +
  ggsidekick::theme_sleek()
ggsave("figures/BC/cog-surveys-XY-ellipse.pdf", width = 19, height = 8)

# columns produced by applying cgi(): year      xcg      ycg         I      Imax     Imin       Iso    xaxe1    yaxe1    xaxe2    yaxe2
#
# additional columns: type ("survey", "commercial"), response ("density", "cpue", "catch"), gear ("bottom_trawl", "longline"…), species (species common name, or the species group targeted by a commercial fishery sector, or simply "groundfish" if general), region ("BC,"WC","AK","NE"), subregion ("GOA"….), gic, xcg_model, ycg_model (COG estimates from sdmTMB),  xcg_model_lwr (lower bound of 95% CI for COG), ycg_model_upr, xcg_model_se, ycg_model_se
#
# One other note on standardization:  Let's try to be consistent with units, where catch will be in kg, density and cpue as kg/square km, and coordinates/distances in units of km (note that when you project geographic coordinates to utm, they will be in meters).