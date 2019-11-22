library(sdmTMB)
library(ggplot2)
library(viridis)
library(tidyverse)
library(fpc)
library(gridExtra)
library(sp)
library(broom)
library(ggforce) # for plotting ellipses

Predict_data_years = readRDS("data/AK/AK_BTS/GOA_predict_data.rds") # save prediction grid

species = sort(c("Dover sole","arrowtooth flounder", "Pacific halibut",
            "walleye pollock", "rex sole", "English sole","sablefish","Pacific cod",
            "spiny dogfish","longnose skate","big skate", "Pacific ocean perch"))

anisotropy_plots = list()
qq_plots = list()
residuals_plots = list()
prediction_plots = list()
spatiotemporal_plots = list()
intercept_plots = list()
intercept_cluster_plots = list()
COG_plots_N = list()
COG_plots_E = list()
mycgi = list()
mycgi_cross_plots = list()
mycgi_ellipse_plots = list()
gic_plots = list()

# plotting functions
plot_map_point <- function(dat, column = "omega_s") {
  ggplot(dat, aes_string("X", "Y", colour = column)) +
    geom_point(size=0.1) +
    xlab("Longitude") +
    ylab("Latitude")
}
plot_map_raster <- function(dat, column = "omega_s") {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    scale_fill_viridis_c() +
    xlab("Longitude") +
    ylab("Latitude")
}

# functions to calculate COGs and inertia in 2 dimensions following Jordan's approach
source("Spatial_indicators_functions_Woillez2009_modified.R")
mycgifun <- function(mycgi){
  return(bind_cols(mycgi %>% 
                     dplyr::select(year,xaxe1,xaxe2) %>% 
                     gather(xaxis,xval,-c(year)) %>% 
                     dplyr::select(-xaxis),
                   mycgi %>% dplyr::select(year,yaxe1,yaxe2) %>% 
                     gather(yaxis,yval,-year) %>% 
                     dplyr::select(-yaxis,-year)))
}

# loop over species
for(spp in 1:length(species)) {
  d = readRDS(paste0("output/AK/", species[spp],"_density.rds"))
  # below 2 lines necessary for models fit with older sdmTMB versions
  #d$tmb_data$weights_i = rep(1, length(d$tmb_data$y_i))
  #d$tmb_data$calc_quadratic_range = as.integer(FALSE)
  p = predict(d, newdata=Predict_data_years)
  
  # calculate COGs and inertia in 2 dimensions
  mycgi[[spp]] <- p %>% 
    group_by(year) %>% 
    do((cgi(x=.$X,y=.$Y,z=exp(.$est)) %>% data.frame)) %>% 
    data.frame
  
  mycgi_ellipse_plots[[spp]] <- mycgifun(mycgi[[spp]]) %>% 
    ggplot(aes(xval,yval,fill=factor(year),color=factor(year))) +
    geom_mark_ellipse(expand = unit(0, "mm"),alpha=0.1) +
    scale_y_continuous(expand = expand_scale(mult = .15)) +
    scale_x_continuous(expand = expand_scale(mult = .15)) +
    theme_bw() + 
    theme(legend.position = "none", axis.title=element_blank()) +
    ggtitle(species[spp])
  
  mycgi_cross_plots[[spp]] <- ggplot() + 
    geom_path(data=mycgi[[spp]],aes(xaxe1,yaxe1)) + 
    geom_path(data=mycgi[[spp]],aes(xaxe2,yaxe2)) +
    facet_wrap(~year,ncol=5) + 
    geom_vline(xintercept = (mycgi[[spp]] %>% summarise(mean(xcg)))[[1]],linetype=2) + 
    geom_hline(yintercept = (mycgi[[spp]] %>% summarise(mean(ycg)))[[1]],linetype=2) +
    xlab("Eastings (10s km)") +
    ylab("Northings (10s km)") +
    ggtitle(paste0(species[spp],"_COG_Inertia"))
  
  # global index of collocation, comparing all years to individual years
  if(!(spp %in% c(9:12))){
    gic_plots[[spp]] <- p %>% 
      group_by(year) %>% 
      summarise(gic=gic(x1=X,
                        y1=Y,
                        z1=exp(est),
                        x2=p$X, 
                        y2=p$Y,
                        z2=exp(p$est))) %>%
      ggplot(aes(year,gic)) + 
      geom_line() +
      theme(axis.title=element_blank(), axis.text.x = element_blank()) +
      ggtitle(species[spp])
  }
  else{
    gic_plots[[spp]] <- p %>% 
      group_by(year) %>% 
      summarise(gic=gic(x1=X,
                        y1=Y,
                        z1=exp(est),
                        x2=p$X, # These funky lines allow us to refer to all years of the data and not just the current group_by year
                        y2=p$Y,
                        z2=exp(p$est))) %>%
      ggplot(aes(year,gic)) + 
      geom_line() +
      theme(axis.title=element_blank()) +
      scale_x_continuous(breaks=seq(1985, 2015, 10)) + 
      ggtitle(species[spp])
  }
  
  # make timeseries plots of COG and 95% CI from model estimates
  p_All = predict(d, newdata=Predict_data_years, return_tmb_object = TRUE)  
  COG = get_cog(p_All)
  # add x axis labels and legend only to specific panels to be combined later
  if(!(spp %in% c(9:12))){
    COG_plots_N[[spp]] = ggplot(filter(COG, coord == "Y"), aes(year, est)) +
      geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey70") +
      geom_line(color = "black") +
      labs(title = species[spp], tag = LETTERS[spp]) +
      theme(axis.title = element_blank(), title = element_text(size = rel(0.9)),
            axis.text.x = element_blank(), plot.margin = unit(c(0,0,1,3), "pt"),
            legend.position = "none")
  }
  else{
    COG_plots_N[[spp]] = ggplot(filter(COG, coord == "Y"), aes(year, est)) +
      geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey70") +
      geom_line(color = "black") +
      labs(title = species[spp], tag = LETTERS[spp]) +
      theme(axis.title = element_blank(), title = element_text(size = rel(0.9)),
            plot.margin = unit(c(0,0,0,3), "pt"),
            legend.position = "none") +
      scale_x_continuous(breaks=seq(1985, 2015, 10))
  }
  
  if(!(spp %in% c(9:12))){
    COG_plots_E[[spp]] = ggplot(filter(COG, coord == "X"), aes(year, est)) +
      geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey70") +
      geom_line(color = "black") +
      labs(title = species[spp], tag = LETTERS[spp]) +
      theme(axis.title = element_blank(), title = element_text(size = rel(0.9)),
            axis.text.x = element_blank(), plot.margin = unit(c(0,0,1,3), "pt"),
            legend.position = "none")
  }
  else{
    COG_plots_E[[spp]] = ggplot(filter(COG, coord == "X"), aes(year, est)) +
      geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey70") +
      geom_line(color = "black") +
      labs(title = species[spp], tag = LETTERS[spp]) +
      theme(axis.title = element_blank(), title = element_text(size = rel(0.9)),
            plot.margin = unit(c(0,0,0,3), "pt"),
            legend.position = "none") +
      scale_x_continuous(breaks=seq(1985, 2015, 10))
  }
  
## model checking
# check anisotropy
anisotropy_plots[[spp]] = plot_anisotropy(d) +
  ggtitle(paste0(species[spp],"_aniso"))
## residuals
data = select(d$data, X, Y, cpue, year)
data$residuals = residuals(d)
# qq plots
qq_plots[[spp]] = qqnorm(data$residuals)
# spatial residuals
residuals_plots[[spp]] = plot_map_point(data, "residuals") + facet_wrap(~year) + geom_point(size=0.05, alpha=0.1) +
  coord_fixed() + scale_color_gradient2() + ggtitle(species[spp])
# check convergence
#sd = as.data.frame(summary(TMB::sdreport(d$tmb_obj)))
#sink(file = "output/AK/sdreport.txt", append = TRUE)
#print(species[spp])
#print(d$sd_report)
#sink()
# check whether AR1 assumption is supported in models where fields are not IID, printing estimate and 95%CI for AR1 param
#print("AR1 estimate")
#print(sd$Estimate[row.names(sd) == "ar1_phi"])
#print("AR1 parameter 95% CI")
#print(sd$Estimate[row.names(sd) == "ar1_phi"] +
# c(-2, 2) * sd$`Std. Error`[row.names(sd) == "ar1_phi"])

# make plot of predictions from full model (all fixed + random effects)
prediction_plots[[spp]] = plot_map_raster(p, "est") +
  facet_wrap(~year) +
  coord_fixed() +
  ggtitle(paste0(species[spp],"_predicted_density"))
# make plot of predictions from spatiotemporal random effects
spatiotemporal_plots[[spp]] = plot_map_raster(p, "epsilon_st") +
  facet_wrap(~year) +
  coord_fixed() +
  ggtitle(paste0(species[spp],"_ST"))

# TO DO: make plot of average biomass, rather than spatial intercept as below
#intercept_plots[[spp]] = plot_map_raster(dplyr::filter(p,year==min(Predict_data_years$year)), "omega_s") +
#  theme(plot.title = element_blank(),
#        axis.title.x = element_text(margin = margin(t = -20)),
#        axis.title.y = element_blank(),
#        axis.text = element_blank(),
#        legend.key.width = unit(0.2,"cm"),
#        legend.title = element_blank(),
#        legend.position = c(0.1,0.8)) +
#  labs(x = "intercept")

}


# save results
save.image(file = "output/AK/plotdata_all.Rdata")

# plot COG and inertia as crosses and ellipses separately
ggsave(filename = "figures/AK/AK_BTS/crosses.pdf",
       plot = marrangeGrob(grobs = mycgi_cross_plots, nrow = 1, ncol = 1),
       width = 7, height = 7, units = c("in"))
ggsave(filename = "figures/AK/AK_BTS/ellipses.pdf",
       plot = arrangeGrob(grobs = mycgi_ellipse_plots, ncol = 4, bottom = "Eastings (10s km)",
                          left = grid::textGrob("Northings (10s km)", rot = 90, vjust = 0.2)),
       width = 7, height = 9, units = c("in"))
# plot legend separately for now
pdf(file = "figures/AK/AK_BTS/ellipses_legend.pdf", width = 1, height = 4)
plot_legend <- mycgifun(mycgi[[1]]) %>% 
  ggplot(aes(xval,yval,fill=factor(year), color = factor(year))) +
  geom_mark_ellipse(expand = unit(0, "mm"),alpha=0.1)
legend <- cowplot::get_legend(plot_legend)
grid::grid.newpage()
grid::grid.draw(legend)
dev.off()
# plot GIC by year
ggsave(filename = "figures/AK/AK_BTS/GIC.pdf",
       plot = arrangeGrob(grobs = gic_plots, ncol = 4, bottom = "Year",
                          left = grid::textGrob("Global index of collocation", rot = 90, vjust = 0.2)),
       width = 9, height = 7, units = c("in"))

# COG timeseries plots from model output
ggsave(filename = "figures/AK/AK_BTS/COG_model_est_N.pdf",
       plot = arrangeGrob(grobs = COG_plots_N, ncol = 4, bottom = "Year",
                          left = grid::textGrob("COG Northings (10s km)", rot = 90, vjust = 0.2)),
       width = 12, height = 8, units = c("in"))
ggsave(filename = "figures/AK/AK_BTS/COG_model_est_E.pdf",
       plot = arrangeGrob(grobs = COG_plots_E, ncol = 4, bottom = "Year",
                          left = grid::textGrob("COG Eastings (10s km)", rot = 90, vjust = 0.2)),
       width = 12, height = 8, units = c("in"))

# plot of predictions from full model (all fixed + random effects)
ggsave(filename = "figures/AK/AK_BTS/predicted_density_maps.pdf",
       plot = marrangeGrob(prediction_plots, nrow = 1, ncol = 1),
       width = 7, height = 9, units = c("in"))

# plot only spatiotemporal random effects
ggsave(filename = "figures/AK/AK_BTS/st_maps.pdf",
       plot = marrangeGrob(spatiotemporal_plots, nrow = 1, ncol = 1),
       width = 7, height = 9, units = c("in"))

## model checking plots
# anisotropy plots for each species
ggsave(filename = "figures/AK/AK_BTS/anisotropy.pdf",
       plot = arrangeGrob(grobs = anisotropy_plots, ncol = 4),
       width = 12, height = 8, units = c("in"))
# qq norm plots for each species
pdf(file = "figures/AK/AK_BTS/qqnorm.pdf", width = 12, height = 12)
par(mfrow=c(3,4))
for(spp in 1:length(species)) {
  plot(qq_plots[[spp]], pch = ".", main = paste0(species[spp],"_qq")); abline(a = 0, b = 1)
}
dev.off()
# residuals maps
ggsave(filename = "figures/AK/AK_BTS/residuals_maps.pdf",
       plot = marrangeGrob(residuals_plots, nrow = 1, ncol = 1),
       width = 7, height = 9, units = c("in"))
