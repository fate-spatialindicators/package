library(dplyr)
library(ggplot2)
source("Spatial_indicators_functions_Woillez2009_modified.R")

# created by `01_run_sdmTMB_BC.R`:
d <- readRDS("data/BC/generated/bc-survey-predictions2-df.rds") %>%
  mutate(density = exp(est))

mycgi <- d %>%
  group_by(year, survey, species) %>%
  do((cgi(x = .$X, y = .$Y, z = .$density) %>% data.frame())) %>%
  data.frame() %>%
  as_tibble()

survey_order <- c("SYN WCVI", "SYN QCS", "SYN HS", "SYN WCHG")
mycgi$survey <- forcats::fct_relevel(mycgi$survey, survey_order)

mean_data <- mycgi %>%
  filter(!is.na(year)) %>%
  group_by(survey, species) %>%
  summarise(mean_xcg = mean(xcg), mean_ycg = mean(ycg))

all_data <- mycgi %>%
  filter(!is.na(year)) %>%
  mutate(survey_species = paste0(survey, "-", species))

g <- ggplot() +
  geom_path(data = all_data, aes(xaxe1, yaxe1, color = year, group = as.factor(year))) +
  geom_path(data = all_data, aes(xaxe2, yaxe2, color = year, group = as.factor(year))) +
  geom_vline(data = mean_data, aes(xintercept = mean_xcg), linetype = 2) +
  geom_hline(data = mean_data, aes(yintercept = mean_ycg), linetype = 2) +
  facet_wrap(survey ~ species, ncol = 11, scales = "free") +
  scale_color_viridis_c() +
  ggsidekick::theme_sleek() +
  theme(legend.position = "none")
ggsave("figures/BC/survey-mcgi.pdf", width = 26, height = 12)

mycgifun <- function(mycgi) {
  dplyr::bind_cols(
    mycgi %>%
      dplyr::select(year, xaxe1, xaxe2, survey_species) %>%
      tidyr::gather(xaxis, xval, -c(year, survey_species)) %>%
      dplyr::select(-xaxis),
    mycgi %>% dplyr::select(year, yaxe1, yaxe2, survey_species) %>%
      tidyr::gather(yaxis, yval, -c(year, survey_species)) %>%
      dplyr::select(-yaxis, -year)
  )
}

g <- all_data %>%
  mycgifun() %>%
  tidyr::separate(survey_species, into = c("survey", "species"), sep = "-") %>%
  mutate(survey = forcats::fct_relevel(survey, survey_order)) %>%
  ggplot(aes(xval, yval, fill = year, color = year, group = as.factor(year))) +
  ggforce::geom_mark_ellipse(expand = unit(0, "mm"), alpha = 0.1) +
  facet_wrap(survey ~ species, ncol = 11, scales = "free") +
  scale_fill_viridis_c() +
  scale_colour_viridis_c() +
  ggsidekick::theme_sleek() +
  theme(legend.position = "none")
ggsave("figures/BC/survey-mcgi-ellipse.pdf", width = 26, height = 12)
