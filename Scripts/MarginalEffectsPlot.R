## ---------------------------
##
## Script name: Marginal Effects
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-10-20
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes: Estimates the fixed effects panel model of growth and
##   plots the marginal effects of climate variables
##
## ---------------------------


## packages
librarian::shelf(sjPlot, lme4, patchwork, tidyverse, ggeffects,
                 fixest, clubSandwich, here, lmtest, sandwich, marginaleffects)

theme_set(
  theme_bw(base_size = 15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)


## data
rwldat <- read_csv(here("Data", "Primary_data", "paneldat_RWL.csv"))
climdat <- read_csv(here("Data","Primary_data", "climatewindows.csv"))

paneldat <- rwldat %>% 
  left_join(climdat)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The panel model 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod <-  feols(log(value+0.01) ~ tmax * ppt + year | tree_id,
                 data= paneldat, cluster = ~ plot_id_needle)

summary(fe_mod)
modelsummary(fe_mod)

tab_model(fe_mod, digits=4, show.ci = F, show.se = T)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Margins plot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Compute temperature marginal effects
mtemp <- slopes(fe_mod, 
              variables = "tmax",
              by = "ppt",
              newdata = datagrid(ppt = seq(min(paneldat$ppt), max(paneldat$ppt), length.out = 100)))

# the plot
marg_temp = ggplot(mtemp, aes(x = ppt*1000, y = estimate, color)) +
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(color="#91376F") +
  ylim(-.03, .08)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color="#91376F", fill="#91376F") +
  labs(x = "Precipitation (mm)", 
       y = "Marginal effect of temperature")+
  scale_x_continuous(breaks = seq(0, 3800, by = 500))


# Compute marginal effects of precip
mprecip <- slopes(fe_mod, 
                variables = "ppt",
                by = "tmax",
                newdata = datagrid(tmax = seq(min(paneldat$tmax), max(paneldat$tmax), length.out = 100)))

# The plot
marg_precip = ggplot(mprecip, aes(x = tmax, y = estimate)) +
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(color= "#303077") +
  ylim(-.2, .45)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color="#303077", fill = "#303077") +
  labs(x = "Temperature (°C)", 
       y = "Marginal effect of precipitation")



marg_temp + marg_precip +   plot_annotation(tag_levels = "A")
 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Adding density distribution panels to both plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Add labels A and B to the density plots only
ppt_density <- ggplot(paneldat, aes(x = ppt*1000)) +
  geom_density(fill = "#91376F", alpha = 0.6) +
  scale_x_continuous(breaks = seq(0, 3800, by = 500)) +
  scale_y_continuous(
    breaks = function(x) pretty(x, n = 3),
    labels = function(x) format(x, scientific = FALSE)  # Add this line
  ) +
  ggtitle("A") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = margin(5, 0, -5, 0),  # Negative bottom margin to reduce gap
    panel.border = element_blank(),
    axis.line.x = element_line(color = "gray70"),
    axis.line.y = element_line(color = "gray70"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0, vjust = 1),
    plot.title.position = "plot"
  ) +
  labs(y = "Density")


# Adjust the temperature margins plot
marg_temp <- marg_temp +
  theme(plot.margin = margin(0, 0, 0, 0))

# Combine with temperature marginal effects plot
combined_temp_plot <- ppt_density / marg_temp +
  plot_layout(heights = c(2, 5), guides = "collect") +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5))) &
  theme(plot.margin = margin(5, 5, 5, 5))

# Add label B to the temperature density plot
tmax_density <- ggplot(paneldat, aes(x = tmax)) +
  geom_density(fill = "#303077", alpha = 0.6) +
  scale_y_continuous(breaks = function(x) pretty(x, n = 3)) +
  ggtitle("B") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = margin(5, 0, -5, 0),  # Negative bottom margin to reduce gap
    panel.border = element_blank(),
    axis.line.x = element_line(color = "gray70"),
    axis.line.y = element_line(color = "gray70"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0, vjust = 1),
    plot.title.position = "plot"
  ) +
  labs(y = "Density")

# Adjust the precipitation margins plot
marg_precip <- marg_precip +
  theme(plot.margin = margin(0, 0, 0, 0))

# Combine with precipitation marginal effects plot
combined_precip_plot <- tmax_density / marg_precip +
  plot_layout(heights = c(2, 5), guides = "collect") +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5))) &
  theme(plot.margin = margin(5, 5, 5, 5))

# Display the combined plots
# combined_temp_plot
# combined_precip_plot

# Set top margin of margin plots to negative to reduce gap
marg_temp <- marg_temp +
  theme(
    plot.margin = margin(-5, 0, 5, 0)  # Negative top margin to reduce gap
  )

# Set top margin of precipitation plot to negative to reduce gap
marg_precip <- marg_precip +
  theme(
    plot.margin = margin(-5, 0, 5, 0)  # Negative top margin to reduce gap
  )

# Combine with minimal spacing between panels
combined_temp_plot <- ppt_density / marg_temp +
  plot_layout(heights = c(1.5, 5)) & 
  theme(plot.margin = margin(0, 5, 0, 5))

combined_precip_plot <- tmax_density / marg_precip +
  plot_layout(heights = c(1.5, 5)) & 
  theme(plot.margin = margin(0, 5, 0, 5))

# Wrap each combined plot without additional theming
wrapped_temp_plot <- combined_temp_plot
wrapped_precip_plot <- combined_precip_plot

# Combine the wrapped plots
final_figure <- wrapped_temp_plot | wrapped_precip_plot

# Display the final combined figure
final_figure

