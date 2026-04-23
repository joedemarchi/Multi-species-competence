#Q1 Plots

# load libraries

# set working directory
setwd("C:/Users/demar/OneDrive/Documents/UTK/Research/Projects/multi_sp_competence")


# Load required packages
library(dplyr)
library(ggplot2)
library(scales)
library(readr)
library(data.table)

library(patchwork)
library(cowplot)
#load plots for all 6 
library(png)
library(grid)
library(ggplotify)


abundance = read.csv("amphibianTN_field/abundance.csv") #Site and species-level abundance
#abundance$species = factor(abundance$species, levels = species_levels)

abundance_wh <- subset(abundance, site == 15) # Subset for WH site)
abundance_wh <- subset(abundance_wh, species %in% c("novi", "racl", "raca"))

mean_abundance <- abundance_wh %>%
  group_by(species) %>%
  summarize(mean_abundance = mean(abund, na.rm = TRUE),
            sd_abundance   = sd(abund, na.rm = TRUE),
            n              = sum(!is.na(abund)),
            se_abundance   = sd_abundance / sqrt(n),
            min_abundance  = min(abund, na.rm = TRUE),
            max_abundance  = max(abund, na.rm = TRUE),
            .groups = "drop"
  )

mean_abundance


movement_distance <- data.frame(
  species = c("N.viridescens", "R.clamitans", "R.catesbeianus"),
  movement_distance = c(10, 40, 30)*0.5 # movement distances in m^3 per 6 days in a pond 0.5 m deep
)

#--------------------------------------------------
# Add in field data
#-------------------------------------------------

# Species and site labels ----

labels = c("hych"="H. chysoscelis","novi"="N. viridescens","pscr"="P. crucifer","psfe"="P. feriarum","raca"="R. catesbeianus","racl"="R. clamitans")
species_levels = names(labels)

# Set working directory ----
setwd("C:/Users/demar/OneDrive/Documents/UTK/Research/Projects/multi_sp_competence/amphibianTN_field")
# Load datasets ----



## Raw data ----

### Site & Survey level environmental data ----
# Read in dataset
pond_data = read.csv("pond_data.csv", stringsAsFactors = FALSE) #Pond sampling level data (may have multiple rows for same survey when multiple water quality measurements (ysi) were taken for a single site)

# Read in field descriptions
comments_df <- read.table("pond_data_comments.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Reattach comments
for (i in seq_len(nrow(comments_df))) {
  col_name <- comments_df$column[i]
  col_comment <- comments_df$comment[i]
  if (col_name %in% names(pond_data)) {
    comment(pond_data[[col_name]]) <- col_comment
  }
}

str(pond_data)

pond_data_lim = pond_data[pond_data$site_id==2 | pond_data$site_id==7 | pond_data$site_id==10 | pond_data$site_id==15,]


### Area ----
pond_area = read.csv("pond_area_df.csv") #Sampling-level data of pond area for Lotus and William Hastie
site_area = read.csv("site_area_df.csv") #Area of full site sampled (land and water)

### Water temperature ----
temperature = read.csv("temperature_dat.csv") #Site-level water temperature 

### Abundance ----
abundance = read.csv("abundance.csv") #Site and species-level abundance
abundance$species = factor(abundance$species, levels = species_levels)


## Interpolated data ----
### Water temperature ----
temperature_gam = read.csv("temp_gam_df.csv") #Site-level water temperature interpolated for the full year

### Abundance ----
abundance_gam = read.csv("abund_gam_df.csv") #Site and species-level abundance interpolated for the full year
abundance_gam$species = factor(abundance_gam$species, levels = species_levels)


abundance$species = factor(abundance$species, levels = species_levels)

# mean_abundance <- abundance_wh %>%
#   group_by(species) %>%
#   summarize(mean_abundance = mean(abund, na.rm = TRUE),
#             sd_abundance   = sd(abund, na.rm = TRUE),
#             n              = sum(!is.na(abund)),
#             se_abundance   = sd_abundance / sqrt(n),
#             min_abundance  = min(abund, na.rm = TRUE),
#             max_abundance  = max(abund, na.rm = TRUE),
#             .groups = "drop"
#   )

### Calculate tadpole
#load data
tadpole_abundance<-read.csv("amphibianTN_field/pipe_sampling.csv")

# subset jsut for william hastie
tad_abund <- subset(tadpole_abundance, subsite_id %like% "WH")
head(tad_abund)
#find mean abundance for each species
# Individual species abundance

# mean tadpoles per pipe
mean_count <- mean(tad_abund$tadpole_count, na.rm = TRUE)

# pipe area and pond area
pipe_area <- 0.8  # m^2
pond_radius <- 15 # m
pond_area <- pi * pond_radius^2

# by night
night_density <- tadpole_abundance %>%
  group_by(survey_id) %>% 
  summarize(
    mean_count = mean(tadpole_count),
    density = mean_count / 0.8
  )

# Night abundance estimate
night_density$abund_est <- night_density$density*pond_area

mean(night_density$abund_est) # approx 1000 tadpoles


### boot strap it to check
boot_abund <- replicate(5000, {
  m <- mean(sample(night_density$density, replace = TRUE))
  m * pond_area
})

quantile(boot_abund, c(0.025, 0.5, 0.975)) # ~2000

## Figure out what proportion of total tadpoles are bullfrogs vs green frogs
# adult density
N_g <- 238
N_b <- 23

#clutch size
c_g <- 5000
c_b <- 15000
# surval of egg to tadpole
s = 0.10 # survival rate or eggs

#total average tadpoles in pond

total_tads <- 2000

## find out proportion of tadpoles

## Relative contribution (egg × clutch × survival)
T_g_rel <- N_g * c_g * s
T_b_rel <- N_b * c_b * s

## Species proportions
p_g <- T_g_rel / (T_g_rel + T_b_rel)
p_b <- T_b_rel / (T_g_rel + T_b_rel)

p_g
p_b

## Apply proportions to total observed tadpoles
greenfrog_tads  <- p_g * total_tads
bullfrog_tads   <- p_b * total_tads

greenfrog_tads # 3400
bullfrog_tads # 1000

## Bootstrap to get CI
B <- 5000 # sims
# bootstrap storage
p_g_boot <- numeric(B)
p_b_boot <- numeric(B)
green_tads_boot <- numeric(B)
bull_tads_boot  <- numeric(B)

for (i in 1:B) {
  
  # Bootstrap adult numbers (±10%)
  N_g_i <- rnorm(1, mean = 238, sd = 0.10 * 238)
  N_b_i <- rnorm(1, mean = 23,  sd = 0.10 * 23)
  
  # Bootstrap clutch sizes (SD = 10% of mean)
  c_g_i <- rnorm(1, mean = 5000,  sd = 0.10 * 5000)
  c_b_i <- rnorm(1, mean = 15000, sd = 0.10 * 15000)
  
  # Bootstrap survival (uniform or beta)
  s_i <- runif(1, 0.10, 0.20)   # survival between 10%–20%
  
  # Bootstrap uncertainty in total tadpole abundance
  total_tads_i <- rnorm(1, mean = 2000, sd = 0.10 * 2000)
  
  # Compute relative contributions
  T_g_rel_i <- N_g_i * c_g_i * s_i
  T_b_rel_i <- N_b_i * c_b_i * s_i
  
  # Species proportions
  p_g_boot[i] <- T_g_rel_i / (T_g_rel_i + T_b_rel_i)
  p_b_boot[i] <- 1 - p_g_boot[i]
  
  # Apply to total abundance
  green_tads_boot[i] <- p_g_boot[i] * total_tads_i
  bull_tads_boot[i]  <- p_b_boot[i] * total_tads_i
}

# Summaries
quantile(p_g_boot, c(0.025, 0.5, 0.975))
quantile(p_b_boot, c(0.025, 0.5, 0.975))

quantile(green_tads_boot, c(0.025, 0.5, 0.975))
quantile(bull_tads_boot, c(0.025, 0.5, 0.975))
bull_ci <- quantile(bull_tads_boot,  c(0.025, 0.5, 0.975))
green_ci <- quantile(green_tads_boot, c(0.025, 0.5, 0.975))

# 
# density (tadpoles per m^2)
density_est <- mean_count / pipe_area

# abundance
abundance_est <- density_est * pond_area

abundance_est


#------------------------------------
#Bullfrog Adult
#-------------------------------------

load("bull_r0_df.RData")
# set persitence values

per_L <- calculate_pond_volume(diameter = 22, avg_depth = 0.5)

per_m3 <- calculate_pond_volume(diameter = 22, avg_depth = 0.5)/1000


bf <- mean_abundance %>% 
  filter(species == "raca") %>%
  mutate(density_m3 = mean_abundance/per_m3,
         max_density_m3 = max_abundance/per_m3,
         min_density_m3 = min_abundance/per_m3)


bf$movement_distance <- movement_distance %>%
  filter(species == "R.catesbeianus") %>%
  pull(movement_distance)


obs <- bf %>% slice(1)  # one bullfrog observation

bull_point <- bull_df %>%
  mutate(
    dist = sqrt(
      (movement - obs$movement_distance)^2 +
        (host_m3 - obs$density_m3)^2
    )
  ) %>%
  slice_min(dist, n = 1)

bull_point$R0      #identifies the nearest cell to the red dot
bull_point        # 0.43

# percent of plot greater than 1: 

pct_r01 <- mean(bull_df$R0 > 1, na.rm = TRUE) * 100
#79.57778


## Sanity check
ipm_chk <- ipm_dists
ipm_chk[["R.catesbeianus"]]$movement_distance <- bull_point$movement
ipm_chk[["R.catesbeianus"]]$df_density        <- bull_point$density  # NOT host_m3

r_direct <- compute_R0_stage(
  ipm_chk,
  stage="adult",
  adult_name="R.catesbeianus",
  tadpole_name="Raca.tadpole",
  pond_diameter=22,
  lower=-15, upper=18, bins=150
)

c(grid=bull_point$R0, direct=r_direct)

ipm_dists[["R.catesbeianus"]][c("movement_distance","df_density")]
bull_point[c("movement","density","R0")]


#plot the data
bull_r0_plt <- ggplot(bull_df, aes(x = movement, y = host_m3)) +
  geom_tile(aes(fill = R0)) +
  geom_contour(aes(z = R0), breaks = 1, colour = "white", linewidth = 0.6) +
  scale_fill_viridis_c(
    name = expression(R[0]),
    limits = range(bull_df$R0, na.rm = TRUE)
  ) +
  guides(fill = guide_colorbar(barheight = unit(100, "pt"), barwidth = unit(20, "pt"))) +
  coord_fixed(ratio = (diff(range(bull_df$movement)) / diff(range(bull_df$host_m3)))) +
  theme_minimal(base_size = 16) +
  labs(title = expression(italic("R. catesbeianus") ~ "adult"),
       x = NULL, y = NULL) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),  # ensure no axis titles
    axis.text  = element_text(size = 16),  # hide tick labels
    axis.ticks = element_blank(),  # hide tick marks
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 0)),
    legend.title = element_text(size = 18)
    #plot.margin = margin(8, 8, 8, 8)  # leaves room for shared outer labels
  ) 

bull_r0_plot <- bull_r0_plt +
  # vertical line at observed movement distance
  geom_vline(data = bf, aes(xintercept = movement_distance),
             linetype = "dashed", linewidth = 0.8, color = "tomato2") +
  # horizontal line at observed density (use the column that matches your y-axis units)
  geom_hline(data = bf, aes(yintercept = density_m3),
             linetype = "dashed", linewidth = 0.8, color = "tomato2") +
  # point at the intersection
  geom_point(data = bf, aes(x = movement_distance, y = density_m3),
             shape = 21, fill = "white", stroke = 1.2, size = 3, color = "tomato2") +
  #point at max abundance
  geom_point(data = bf, aes(x = movement_distance, y = max_density_m3),
             shape = 24, fill = "white", stroke = 1, size = 1.5, color = "tomato2") +
  # point at min abundance
  geom_point(data = bf, aes(x = movement_distance, y = min_density_m3),
             shape = 25, fill = "white", stroke = 1, size = 1.5, color = "tomato2") +
  annotate(
    "text",
    x = 50,
    y = 0.045,  # slight vertical offset
    label = paste0("R[0] == ", round(r_direct, 2)),
    parse = TRUE,
    color = "tomato2",
    size = 5,
    fontface = "bold"
  )



print(bull_r0_plot)
ggsave("plots/bull_r0_plot.png", plot = bull_r0_plot, width = 6, height = 5, dpi = 300)


#------------------------------------
#Bullfrog Tadpole
#-------------------------------------

load("btad_r0_df.RData")
# create tad density

per_L <- calculate_pond_volume(diameter = 22, avg_depth = 0.5)

per_m3 <- calculate_pond_volume(diameter = 22, avg_depth = 0.5)/1000

#tad_dens <- (2700/per_m3) # for green frogs
tad_dens <- (450/per_m3) # for bullfrogs

pct_r01_btad <- mean(tad_df$R0 > 1, na.rm = TRUE) * 100

# observed tadpole density
tad_dens <- 450 / per_m3 # convert to m^3
tad_move <- 6   #

tad_point <- tad_df %>%
  mutate(
    dist = sqrt(
      (movement - tad_move)^2 +
        (host_m3 - tad_dens)^2
    )
  ) %>%
  slice_min(dist, n = 1)

tad_point$R0   # R0 at the grid cell closest to that tadpole point
tad_point # 0.722

## Sanity check
ipm_chk <- ipm_dists
ipm_chk[["Raca.tadpole"]]$movement_distance <- tad_point$movement
ipm_chk[["Raca.tadpole"]]$df_density        <- tad_point$density  # NOT host_m3

r_direct <- compute_R0_stage(
  ipm_chk,
  stage="tadpole",
  adult_name="R.catesbeianus",
  tadpole_name="Raca.tadpole",
  pond_diameter=22,
  lower=-15, upper=18, bins=150
)

c(grid=tad_point$R0, direct=r_direct)

ipm_dists[["Raca.tadpole"]][c("movement_distance","df_density")]
tad_point[c("movement","density","R0")]

#plot the data
tad_r0_plt <- ggplot(tad_df, aes(x = movement, y = host_m3)) +
  geom_tile(aes(fill = R0)) +
  geom_contour(aes(z = R0), breaks = 1, colour = "white", linewidth = 0.6) +
  scale_fill_viridis_c(
    name = expression(R[0]),
    limits = range(tad_df$R0, na.rm = TRUE)
  ) +
  guides(fill = guide_colorbar(barheight = unit(100, "pt"), barwidth = unit(20, "pt"))) +
  coord_fixed(ratio = (diff(range(tad_df$movement)) / diff(range(tad_df$host_m3)))) +
  labs(title = expression(italic("R. catesbeianus") ~ "tadpole"),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),  # ensure no axis titles
    axis.text  = element_text(size = 16),  # hide tick labels
    axis.ticks = element_blank(),  # hide tick marks
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 0)),
    legend.title = element_text(size = 18)
    #plot.margin = margin(8, 8, 8, 8),  # leaves room for shared outer labels
  ) +
  geom_vline(aes(xintercept = 6*0.5), linetype = "dashed", linewidth = 0.8, color = "tomato2") +
  geom_hline(aes(yintercept = tad_dens), linetype = "dashed", linewidth = 0.8, color = "tomato2") +
  geom_point(aes(x = 6*0.5, y = tad_dens),
             shape = 21, fill = "white", stroke = 1.2, size = 3, color = "tomato2", inherit.aes = FALSE) +
  #point at max abundance
  geom_point(aes(x = 6*0.5, y = unname((bull_ci[3])/per_m3)),
             shape = 24, fill = "white", stroke = 1, size = 1.5, color = "tomato2", inherit.aes = FALSE) +
  # point at min abundance
  geom_point(aes(x = 6*0.5, y = unname((bull_ci[1])/per_m3)),
             shape = 25, fill = "white", stroke = 1, size = 1.5, color = "tomato2", inherit.aes = FALSE) +
  annotate(
    "text",
    x = 8.5,
    y = 3,  # slight vertical offset
    label = paste0("R[0] == ", round(r_direct, 2)),
    parse = TRUE,
    color = "tomato2",
    size = 5,
    fontface = "bold"
  )



print(tad_r0_plt)

ggsave("plots/tad_r0_plot.png", plot = tad_r0_plt, width = 6, height = 5, dpi = 300)

#------------------------------------
#Green frog Adult
#-------------------------------------

load("green_r0_df.RData")


per_L <- calculate_pond_volume(diameter = 22, avg_depth = 0.5)

per_m3 <- calculate_pond_volume(diameter = 22, avg_depth = 0.5)/1000

gf <- mean_abundance %>% 
  filter(species == "racl") %>%
  mutate(density_m3 = mean_abundance/per_m3,
         max_density_m3 = max_abundance/per_m3,
         min_density_m3 = min_abundance/per_m3)


gf$movement_distance <- movement_distance %>%
  filter(species == "R.clamitans") %>%
  pull(movement_distance)


obs <- gf %>% slice(1)  # one bullfrog observation

green_point <- var_df %>%
  mutate(
    dist = sqrt(
      (movement - obs$movement_distance)^2 +
        (host_m3 - obs$density_m3)^2
    )
  ) %>%
  slice_min(dist, n = 1)

green_point$R0      #identifies the nearest cell to the red dot
green_point        # 0.15

pct_r01_green <- mean(var_df$R0 > 1, na.rm = TRUE) * 100
pct_r01_green
#40.71

## Sanity check
ipm_chk <- ipm_dists
ipm_chk[["R.clamitans"]]$movement_distance <- green_point$movement
ipm_chk[["R.clamitans"]]$df_density        <- green_point$density  # NOT host_m3

r_direct <- compute_R0_stage(
  ipm_chk,
  stage="adult",
  adult_name="R.clamitans",
  tadpole_name="Raca.tadpole",
  pond_diameter=22,
  lower=-15, upper=18, bins=150
)

c(grid=green_point$R0, direct=r_direct)

ipm_dists[["R.clamitans"]][c("movement_distance","df_density")]
green_point[c("movement","density","R0")]


#plot the data
r0_plt <- ggplot(var_df, aes(x = movement, y = host_m3)) +
  geom_tile(aes(fill = R0)) +
  geom_contour(aes(z = R0), breaks = 1, colour = "white", linewidth = 0.6) +
  scale_fill_viridis_c(
    name = expression(R[0]),
    limits = range(var_df$R0, na.rm = TRUE)
  ) +
  guides(fill = guide_colorbar(barheight = unit(100, "pt"), barwidth = unit(20, "pt"))) +
  coord_fixed(ratio = (diff(range(var_df$movement)) / diff(range(var_df$host_m3)))) +
  theme_minimal(base_size = 16) +
  labs(title = expression(italic("R. clamitans") ~ "adult"),
       x = NULL, y = NULL) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),  # ensure no axis titles
    axis.text  = element_text(size = 16),  # hide tick labels
    axis.ticks = element_blank(),  # hide tick marks
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 0)),
    legend.title = element_text(size = 18)
    #plot.margin = margin(8, 8, 8, 8)  # leaves room for shared outer labels
  )

green_r0_plot <- r0_plt +
  geom_vline(data = gf, aes(xintercept = movement_distance),
             linetype = "dashed", linewidth = 0.8, color = "tomato2") +
  # horizontal line at observed density (use the column that matches your y-axis units)
  geom_hline(data = gf, aes(yintercept = density_m3),
             linetype = "dashed", linewidth = 0.8, color = "tomato2") +
  # point at the intersection
  geom_point(data = gf, aes(x = movement_distance, y = density_m3),
             shape = 21, fill = "white", stroke = 1.2, size = 3, color = "tomato2") +
  #point at max abundance
  geom_point(data = gf, aes(x = movement_distance, y = max_density_m3),
             shape = 24, fill = "white", stroke = 1, size = 1.5, color = "tomato2") +
  # point at min abundance
  geom_point(data = gf, aes(x = movement_distance, y = min_density_m3),
             shape = 25, fill = "white", stroke = 1, size = 1.5, color = "tomato2") +
  # Min–max abundance range
  annotate(
    "text",
    x = 60,
    y = 0.8,  # slight vertical offset
    label = paste0("R[0] == ", round(r_direct, 2)),
    parse = TRUE,
    color = "tomato2",
    size = 5,
    fontface = "bold"
  )

print(green_r0_plot)

ggsave("plots/green_r0_plot.png", plot = green_r0_plot, width = 6, height = 5, dpi = 300)


#------------------------------------
#Greenfrog Tadpole
#-------------------------------------

load("gtad_r0_df.RData")

## Find Host density per liter or pond volume

per_L <- calculate_pond_volume(diameter = 22, avg_depth = 0.5)

per_m3 <- calculate_pond_volume(diameter = 22, avg_depth = 0.5)/1000
# create tad density
gtad_dens <- (1560/per_m3) # for green frogs
#tad_dens <- (420/per_m3) # for bullfrogs

pct_r01_gtad <- mean(gtad_df$R0 > 1, na.rm = TRUE) * 100
pct_r01_gtad
#12.515

gtad_dens <- 1560/per_m3
gtad_move <- 3   # or some other value you choose

gtad_point <- gtad_df %>%
  mutate(
    dist = sqrt(
      (movement - gtad_move)^2 +
        (host_m3 - gtad_dens)^2
    )
  ) %>%
  slice_min(dist, n = 1)

gtad_point$R0   # R0 at the grid cell closest to that tadpole point
gtad_point 

## Sanity check
ipm_chk <- ipm_dists
ipm_chk[["Racl.tadpole"]]$movement_distance <- gtad_point$movement
ipm_chk[["Racl.tadpole"]]$df_density        <- gtad_point$density  # NOT host_m3

r_direct <- compute_R0_stage(
  ipm_chk,
  stage="tadpole",
  adult_name="R.clamitans",
  tadpole_name="Racl.tadpole",
  pond_diameter=22,
  lower=-15, upper=18, bins=150
)

c(grid=gtad_point$R0, direct=r_direct)

ipm_dists[["Racl.tadpole"]][c("movement_distance","df_density")]
gtad_point[c("movement","density","R0")]

#plot the data
gtad_r0_plt <- ggplot(gtad_df, aes(x = movement, y = host_m3)) +
  geom_tile(aes(fill = R0)) +
  geom_contour(aes(z = R0), breaks = 1, colour = "white", linewidth = 0.6) +
  scale_fill_viridis_c(
    name = expression(R[0]),
    #limits = range(tad_df$R0, na.rm = TRUE)
  ) +
  #scale_x_continuous(name = expression("Movement distance (" *m^3* ~ "/" ~ 6 ~ plain("days") * ")"), labels = label_number(accuracy = 1)) +
  #scale_y_continuous(name = expression("Host density ("*m^3*")"), labels = label_number(big.mark = ",")) +
  guides(fill = guide_colorbar(barheight = unit(100, "pt"), barwidth = unit(20, "pt"))) +
  coord_fixed(ratio = (diff(range(gtad_df$movement)) / diff(range(gtad_df$host_m3)))) +
  labs(title = expression(italic("R. clamitans") ~ "tadpole"),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),  # ensure no axis titles
    axis.text  = element_text(size = 16),  # hide tick labels
    axis.ticks = element_blank(),  # hide tick marks
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 0)),
    legend.title = element_text(size = 18)
    #plot.margin = margin(8, 8, 8, 8),  # leaves room for shared outer labels
  ) +
  geom_vline(aes(xintercept = 3*0.5), linetype = "dashed", linewidth = 0.8, color = "tomato2") +
  geom_hline(aes(yintercept = gtad_dens), linetype = "dashed", linewidth = 0.8, color = "tomato2") +
  geom_point(aes(x = 3*0.5, y = gtad_dens),
             shape = 21, fill = "white", stroke = 1.2, size = 3, color = "tomato2") +
  geom_point(aes(x = 3*0.5, y = unname((green_ci[3])/per_m3)),
             shape = 24, fill = "white", stroke = 1, size = 1.5, color = "tomato2", inherit.aes = FALSE) +
  # point at min abundance
  geom_point(aes(x = 3*0.5, y = unname((green_ci[1])/per_m3)),
             shape = 25, fill = "white", stroke = 1, size = 1.5, color = "tomato2", inherit.aes = FALSE) +
  annotate(
    "text",
    x = 5,
    y = 7,  # slight vertical offset
    label = paste0("R[0] == ", round(r_direct, 2)),
    parse = TRUE,
    color = "tomato2",
    size = 5,
    fontface = "bold"
  )


print(gtad_r0_plt)

ggsave("plots/gtad_r0_plot.png", plot = gtad_r0_plt, width = 6, height = 5, dpi = 300)



#------------------------------------
#Newt
#-------------------------------------
load("newt_r0_df.RData")


per_L <- calculate_pond_volume(diameter = 22, avg_depth = 0.5)

per_m3 <- calculate_pond_volume(diameter = 22, avg_depth = 0.5)/1000


# set persitence values

nf <- mean_abundance %>% 
  filter(species == "novi") %>%
  mutate(density_m3 = mean_abundance/per_m3,
         max_density_m3 = max_abundance/per_m3,
         min_density_m3 = min_abundance/per_m3)


nf$movement_distance <- movement_distance %>%
  filter(species == "N.viridescens") %>%
  pull(movement_distance)

### Find R0 of observed density

obs <- nf %>% slice(1)  # one bullfrog observation

newt_point <- newt_df %>%
  mutate(
    dist = sqrt(
      (movement - obs$movement_distance)^2 +
        (host_m3 - obs$density_m3)^2
    )
  ) %>%
  slice_min(dist, n = 1)

newt_point$R0      #identifies the nearest cell to the red dot
newt_point        # 0.43

## Sanity check
ipm_chk <- ipm_dists
ipm_chk[["N.viridescens"]]$movement_distance <- newt_point$movement
ipm_chk[["N.viridescens"]]$df_density <- newt_point$density  # NOT host_m3

r_direct <- compute_R0_stage(
  ipm_chk,
  stage="adult",
  adult_name="N.viridescens",
  tadpole_name="Racl.tadpole",
  pond_diameter=22,
  lower=-15, upper=18, bins=150
)

c(grid=newt_point$R0, direct=r_direct)

ipm_dists[["N.viridescens"]][c("movement_distance","df_density")]
newt_point[c("movement","density","R0")]


#plot the data
newt_r0_plt <- ggplot(newt_df, aes(x = movement, y = host_m3)) +
  geom_tile(aes(fill = R0), linewidth = 0.15) +
  geom_contour(aes(z = R0), breaks = 1, colour = "white", linewidth = 0.6) +
  scale_fill_viridis_c(
    name = expression(R[0]),
    #option = "magma",
    trans = "log10",
    limits = range(newt_df$R0, na.rm = TRUE)
  ) +
  guides(fill = guide_colorbar(barheight = unit(100, "pt"), barwidth = unit(20, "pt"))) +
  coord_fixed(ratio = (diff(range(newt_df$movement)) / diff(range(newt_df$host_m3)))) +
  theme_minimal(base_size = 16) +
  labs(title = expression(italic("N. viridescens") ~ "adult"),
       x = NULL, y = NULL) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),  # ensure no axis titles
    axis.text  = element_text(size = 16),  # hide tick labels
    #axis.ticks = element_text(size = 10),  # hide tick marks
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 0)),
    legend.title = element_text(size = 18)
    #plot.margin = margin(8, 8, 8, 8)  # leaves room for shared outer labels
  )

newt_r0_plot <- newt_r0_plt +
  # vertical line at observed movement distance
  geom_vline(data = nf, aes(xintercept = movement_distance),
             linetype = "dashed", linewidth = 0.8, color = "tomato2") +
  # horizontal line at observed density (use the column that matches your y-axis units)
  geom_hline(data = nf, aes(yintercept = density_m3),
             linetype = "dashed", linewidth = 0.8, color = "tomato2") +
  # point at the intersection
  geom_point(data = nf, aes(x = movement_distance, y = density_m3),
             shape = 21, fill = "white", stroke = 1.2, size = 3, color = "tomato2") +
  # point at min abundance
  geom_point(data = nf, aes(x = movement_distance, y = min_density_m3),
             shape = 25, fill = "white", stroke = 1, size = 1.5, color = "tomato2") +
  annotate(
    "text",
    x = 12,
    y = 4,  # slight vertical offset
    label = paste0("R[0] == ", round(r_direct, 2)),
    parse = TRUE,
    color = "tomato2",
    size = 5,
    fontface = "bold"
  )






print(newt_r0_plot)

ggsave("plots/newt_r0min_plot.png", plot = newt_r0_plot, width = 6, height = 5, dpi = 300)

#------------------------------------
# Combine figures
#-------------------------------------

bull_r0_plot <- as.ggplot(
  rasterGrob(readPNG("plots/bull_r0_plot.png"),
             interpolate = TRUE)
)

newt_r0_plot <- as.ggplot(
  rasterGrob(readPNG("plots/newt_r0min_plot.png"),
             interpolate = TRUE)
)

tad_r0_plot <- as.ggplot(
  rasterGrob(readPNG("plots/tad_r0_plot.png"),
             interpolate = TRUE)
)
green_r0_plot <- as.ggplot(
  rasterGrob(readPNG("plots/green_r0_plot.png"),
             interpolate = TRUE)
)


gtad_r0_plot <- as.ggplot(
  rasterGrob(readPNG("plots/gtad_r0_plot.png"),
             interpolate = TRUE)
)
#### Make grid

# universal theme
pub_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    axis.ticks = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      size = 12,
      face = "plain",
      margin = margin(b = 4)
    ),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    plot.margin = margin(2, 2, 2, 2)
  )

strip_margins <- function(p) {
  p + theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2))
}

# strip_margins <- function(p) {
#   p + theme(plot.margin = margin(0, 0, 0, 0))
# }

newt_r0_plot  <- strip_margins(newt_r0_plot)
bull_r0_plot  <- strip_margins(bull_r0_plot)
green_r0_plot <- strip_margins(green_r0_plot)
tad_r0_plot   <- strip_margins(tad_r0_plot)
gtad_r0_plot <- strip_margins(gtad_r0_plot)

top_row <- plot_grid(
  newt_r0_plot, bull_r0_plot, tad_r0_plot,
  ncol = 3,
  labels = c("A", "B", "C"),
  label_size = 12,
  label_fontface = "bold",
  label_x = 0.03,
  label_y = 0.93,
  hjust = 0,
  vjust = 1
)

bottom_row <- plot_grid(
  NULL, green_r0_plot, gtad_r0_plot, NULL,
  ncol = 4,
  rel_widths = c(0.5, 1, 1, 0.5),
  labels = c("", "D", "E", ""),
  label_size = 12,
  label_fontface = "bold",
  label_x = 0.03,
  label_y = 0.93,
  hjust = 0,
  vjust = 1
)


pgrid <- plot_grid(
  top_row, bottom_row,
  ncol = 1,
  rel_heights = c(1, 1)
)

# pgrid <- cowplot::plot_grid(
#   newt_r0_plot, bull_r0_plot, tad_r0_plot,
#   green_r0_plot, gtad_r0_plot,
#   ncol = 3,
#   align = "hv",
#   axis = "tblr",
#   labels = c("A", "B", "C", "D", "E"),
#   label_size = 14,
#   label_fontface = "bold"
# )

##Make the plot
left_margin   <- 0.09
bottom_margin <- 0.09
right_margin  <- 0.02
top_margin    <- 0.02

final <- ggdraw() +
  draw_plot(
    pgrid,
    x = left_margin,
    y = bottom_margin,
    width  = 1 - left_margin - right_margin,
    height = 1 - bottom_margin - top_margin
  ) +
  annotate(
    "text",
    x = left_margin + (1 - left_margin - right_margin) / 2,
    y = bottom_margin * 0.7,
    size = 5.8,
    hjust = 0.5,
    label = "paste('Movement distance (', frac(m^3, 6~plain('days')), ')')",
    parse = TRUE
  ) +
  annotate(
    "text",
    x = left_margin * 0.65,
    y = bottom_margin + (1 - bottom_margin - top_margin) / 2,
    angle = 90,
    size = 5.8,
    hjust = 0.5,
    label = "paste('Host density (per ', m^3, ')')",
    parse = TRUE
  )

## Save the plot
print(final)

ggsave(
  "plots/R0_all_final.png",
  final,
  width = 12, height = 8,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
  
)

#-----------------------------------------
# Community R0
#----------------------------------------

#--------------------------------------
#Community Bullfrog
#----------------------------------------

load("multi_bf_r0_df.RData")

# observed counts
N_A <- 23
N_T <- 450

# pond volume
per_m3 <- calculate_pond_volume(diameter = 22, avg_depth = 0.5) / 1000

# observed total density (y-axis value)
obs_total_density_m3 <- (N_A + N_T) / per_m3
max_density_m3 <- 720/per_m3
min_density_m3 <- 310/per_m3

# observed movements (same units as x-axis)
m_A <- 30   # adult
m_T <- 6   # tadpole

### Paper graph

pct_r01_ball <- mean(bf_df$R0 > 1, na.rm = TRUE) * 100
pct_r01_ball
#55.815

## compute R0 at exact combined R0 point
# observed counts in pond (DFE)
N_A <- 23
N_T <- 450

# observed movement distances (your units)
m_A <- 30
m_T <- 6

# make a copy so you don't overwrite ipm_dists
ipm_obs <- ipm_dists

# set observed densities by stage (total individuals in pond)
ipm_obs[["R.catesbeianus"]]$df_density <- N_A
ipm_obs[["Raca.tadpole"]]$df_density  <- N_T

# set observed movements by stage
ipm_obs[["R.catesbeianus"]]$movement_distance <- m_A
ipm_obs[["Raca.tadpole"]]$movement_distance  <- m_T

# compute exact combined R0 at that point
R0_obs <- compute_bullfrog_2stage_R0_current(
  ipm_d = ipm_obs,
  adult_name = "R.catesbeianus",
  tadpole_name = "Raca.tadpole",
  pond_diameter = 22,
  lower = -15, upper = 18, bins = 150
)

R0_obs

## Find R0 along obs density value:

R0_target <- as.numeric(R0_obs)
y_obs <- obs_total_density_m3
x_obs <- 1

# choose the closest y level in the grid
y0 <- bf_df %>%
  mutate(dy = abs(host_m3_total - y_obs)) %>%
  arrange(dy) %>%
  slice(1) %>%
  pull(host_m3_total)

# within that y0 slice, find meff that makes R0 closest to R0_obs
pt_slice <- bf_df %>%
  filter(host_m3_total == y0) %>%
  mutate(dR = abs(R0 - R0_target)) %>%
  arrange(dR) %>%
  slice(1)

pt_slice %>% select(k, host_m3_total, R0, dR)


bf_r0_plt <- ggplot(bf_df, aes(x = k, y = host_m3_total)) +
  geom_tile(aes(fill = R0)) +
  geom_contour(aes(z = R0), breaks = 1, colour = "white", linewidth = 0.6) +
  scale_fill_viridis_c(
    name = expression(R[0]),
    limits = range(bf_df$R0, na.rm = TRUE)
  ) +
  guides(fill = guide_colorbar(barheight = unit(100, "pt"),
                               barwidth  = unit(20,  "pt"))) +
  # coord_fixed(
  #   ratio = diff(range(bf_df$meff)) / diff(range(bf_df$host_m3_total))
  # ) +
  theme_minimal(base_size = 16) +
  labs(
    title = "<i>R. catesbeianus</i><br>(adult + tadpole)",
    x = NULL,
    y = NULL
  ) +
  theme(
    panel.grid  = element_blank(),
    axis.title  = element_blank(),
    axis.text   = element_text(size = 16),
    axis.ticks  = element_blank(),
    plot.title   = ggtext::element_markdown(
      hjust = 0.5,
      margin = margin(b = 0),
      lineheight = 1.1
    ),
    legend.title = element_text(size = 16)
  )


bf_r0_plt <- bf_r0_plt +
  geom_vline(
    xintercept = pt_slice$k,
    linetype   = "dashed",
    linewidth  = 0.8,
    color      = "tomato2"
  ) +
  geom_hline(
    yintercept = pt_slice$host_m3_total,
    linetype   = "dashed",
    linewidth  = 0.8,
    color      = "tomato2"
  ) +
  annotate("point", 
           x= pt_slice$k,
           y = max_density_m3,
           shape = 24, 
           fill = "white", 
           stroke = 1, 
           size = 1.5, 
           color = "tomato2") +
  annotate("point", 
           x= pt_slice$k,
           y = min_density_m3,
           shape = 24, 
           fill = "white", 
           stroke = 1, 
           size = 1.5, 
           color = "tomato2") +
  annotate(
    "point",
    x = pt_slice$k,
    y = pt_slice$host_m3_total,
    shape  = 21,
    fill   = "white",
    stroke = 1.2,
    size   = 3,
    color  = "tomato2"
  ) +
  annotate(
    "text",
    x = 1.75,
    y = 3,  # slight vertical offset
    label = paste0("R[0] == ", round(R0_target, 2)),
    parse = TRUE,
    color = "tomato2",
    size = 5,
    fontface = "bold"
  )

print(bf_r0_plt)

ggsave("plots/bullfrog_multi_r0_plot.png",
       plot = bf_r0_plt,
       width = 6,
       height = 5,
       dpi = 300,
       bg = "white")


#--------------------------------------
#Community Green frog
#----------------------------------------

load("multi_gf_r0_df.RData")

pct_r01_gall <- mean(gf_df$R0 > 1, na.rm = TRUE) * 100
pct_r01_gall
#38.295
# observed counts
N_A<- ipm_dists[["R.clamitans"]]$df_density
N_T <- ipm_dists[["Racl.tadpole"]]$df_density

# pond volume
per_m3 <- calculate_pond_volume(diameter = 22, avg_depth = 0.5) / 1000

# observed total density (y-axis value)
obs_total_density_m3 <- (N_A + N_T) / per_m3
max_density_m3 <- 2550/per_m3
min_density_m3 <- 1300/per_m3

# observed movements (same units as x-axis)
m_A <- ipm_dists[["R.clamitans"]]$movement_distance   # adult
m_T <- ipm_dists[["Racl.tadpole"]]$movement_distance   # tadpole

# density-weighted average movement
obs_mv <- (N_A * m_A + N_T * m_T) / (N_A + N_T)

### Find R0 point

ipm_obs <- ipm_dists
N_A <- ipm_dists[["R.clamitans"]]$df_density
N_T <- ipm_dists[["Racl.tadpole"]]$df_density
m_A <- ipm_dists[["R.clamitans"]]$movement_distance
m_T <- ipm_dists[["Racl.tadpole"]]$movement_distance

# set observed densities by stage (total individuals in pond)
ipm_obs[["R.clamitans"]]$df_density <- N_A
ipm_obs[["Racl.tadpole"]]$df_density  <- N_T

# set observed movements by stage
ipm_obs[["R.clamitans"]]$movement_distance <- m_A
ipm_obs[["Racl.tadpole"]]$movement_distance  <- m_T

# compute exact combined R0 at that point
R0_obs_g <- compute_bullfrog_2stage_R0_current(
  ipm_d = ipm_obs,
  adult_name = "R.clamitans",
  tadpole_name = "Racl.tadpole",
  pond_diameter = 22,
  lower = -15, upper = 18, bins = 150
)

R0_obs_g

R0_target <- as.numeric(R0_obs_g)
y_obs <- obs_total_density_m3

# choose the closest y level in the grid
y0 <- gf_df %>%
  mutate(dy = abs(host_m3_total - y_obs)) %>%
  arrange(dy) %>%
  slice(1) %>%
  pull(host_m3_total)

# within that y0 slice, find meff that makes R0 closest to R0_obs
pt_slice <- gf_df %>%
  filter(host_m3_total == y0) %>%
  mutate(dR = abs(R0 - R0_target)) %>%
  arrange(dR) %>%
  slice(1)

pt_slice %>% select(k, host_m3_total, R0, dR)

### Graph for journal

gf_r0_plt <- ggplot(gf_df, aes(x = k, y = host_m3_total)) +
  geom_tile(aes(fill = R0)) +
  geom_contour(aes(z = R0), breaks = 1,
               colour = "white", linewidth = 0.6) +
  scale_fill_viridis_c(
    name   = expression(R[0]),
    limits = range(gf_df$R0, na.rm = TRUE)
  ) +
  guides(fill = guide_colorbar(
    barheight = unit(100, "pt"),
    barwidth  = unit(20,  "pt")
  )) +
  theme_minimal(base_size = 16) +
  labs(
    title = "<i>R. clamitans</i><br>(adult + tadpole)",
    x = NULL,
    y = NULL
  ) +
  theme(
    panel.grid   = element_blank(),
    axis.title   = element_blank(),
    axis.text    = element_text(size = 16),
    axis.ticks   = element_blank(),
    plot.title   = ggtext::element_markdown(
      hjust = 0.5,
      margin = margin(b = 0),
      lineheight = 1.1
    ),
    legend.title = element_text(size = 16)
  )

gf_r0_plt2 <- gf_r0_plt +
  geom_vline(
    xintercept = pt_slice$k,
    linetype   = "dashed",
    linewidth  = 0.8,
    color      = "tomato2"
  ) +
  geom_hline(
    yintercept = pt_slice$host_m3_total,
    linetype   = "dashed",
    linewidth  = 0.8,
    color      = "tomato2"
  ) +
  annotate("point", 
           x= pt_slice$k,
           y = max_density_m3,
           shape = 24, 
           fill = "white", 
           stroke = 1, 
           size = 1.5, 
           color = "tomato2") +
  annotate("point", 
           x= pt_slice$k,
           y = min_density_m3,
           shape = 24, 
           fill = "white", 
           stroke = 1, 
           size = 1.5, 
           color = "tomato2") +
  annotate(
    "point",
    x = pt_slice$k,
    y = pt_slice$host_m3_total,
    shape  = 21,
    fill   = "white",
    stroke = 1.2,
    size   = 3,
    color  = "tomato2"
  ) +
  annotate(
    "text",
    x = 1.7,
    y = 8,  # slight vertical offset
    label = paste0("R[0] == ", round(R0_target, 2)),
    parse = TRUE,
    color = "tomato2",
    size = 5,
    fontface = "bold"
  )

print(gf_r0_plt2)

ggsave("plots/greenfrog_multi_r0_plot.png",
       plot = gf_r0_plt2,
       width = 6,
       height = 5,
       dpi = 300,
       bg = "white")


#--------------------------------------
#Community bullfrog + green frog
#----------------------------------------
load("R0_all_df.RData")

N_A <- 261
N_T <- 2010
# pond volume
per_m3 <- calculate_pond_volume(diameter = 22, avg_depth = 0.5) / 1000

# observed total density (y-axis value)
obs_total_density_m3 <- (N_A + N_T) / per_m3
max_density_m3 <- (720 + 2550)/per_m3
min_density_m3 <- (310 + 1300)/per_m3


all_df <- R0_var_all %>%
  mutate(host_m3_total = Ntot / per_m3) %>%
  filter(is.finite(R0))

sp_comm <- c("R.clamitans", "R.catesbeianus", "Raca.tadpole", "Racl.tadpole")

R0_comm <- compute_R0_stage_community(
  species_list = sp_comm,
  ipm_dists = ipm_dists,
  pond_diameter = 22,
  lower = -15, upper = 18, bins = 150,
  include_emigration = TRUE,
  include_infect_then_meta = TRUE,
  z_decay = 0.3
)

R0_comm

R0_target <- as.numeric(R0_comm)
y_obs <- obs_total_density_m3

# choose the closest y level in the grid
y0 <- all_df %>%
  mutate(dy = abs(host_m3_total - y_obs)) %>%
  arrange(dy) %>%
  slice(1) %>%
  pull(host_m3_total)

# within that y0 slice, find meff that makes R0 closest to R0_obs
pt_slice <- all_df %>%
  filter(host_m3_total == y0) %>%
  mutate(dR = abs(R0 - R0_target)) %>%
  arrange(dR) %>%
  slice(1)

pt_slice %>% select(k, host_m3_total, R0, dR)

### Graph for journal

library(ggtext)


all_r0_plt <- ggplot(all_df, aes(x = k, y = host_m3_total)) +
  geom_tile(aes(fill = R0)) +
  geom_contour(aes(z = R0), breaks = 1,
               colour = "white", linewidth = 0.6) +
  scale_fill_viridis_c(
    name   = expression(R[0]),
    limits = range(all_df$R0, na.rm = TRUE)
  ) +
  guides(fill = guide_colorbar(
    barheight = unit(100, "pt"),
    barwidth  = unit(20,  "pt")
  )) +
  # coord_fixed(
  #   ratio = diff(range(gf_df$meff)) /
  #           diff(range(gf_df$host_m3_total))
  # ) +
  theme_minimal(base_size = 16) +
  labs(
    title = "<i>R. catesbeianus</i>, <i>R. clamitans</i><br>(adults + tadpoles)",
    x = NULL,
    y = NULL
  ) +
  theme(
    panel.grid   = element_blank(),
    axis.title   = element_blank(),
    axis.text    = element_text(size = 16),
    axis.ticks   = element_blank(),
    plot.title   = ggtext::element_markdown(
      hjust = 0.5,
      margin = margin(b = 0),
      lineheight = 1.1
    ),
    # plot.title   = element_text(
    #   hjust = 0.5,
    #   face  = "bold",
    #   margin = margin(b = 0)
    # ),
    legend.title = element_text(size = 16)
  )

all_r0_plt2 <- all_r0_plt +
  geom_vline(
    xintercept = pt_slice$k,
    linetype   = "dashed",
    linewidth  = 0.8,
    color      = "tomato2"
  ) +
  geom_hline(
    yintercept = pt_slice$host_m3_total,
    linetype   = "dashed",
    linewidth  = 0.8,
    color      = "tomato2"
  ) +
  annotate("point", 
           x= pt_slice$k,
           y = max_density_m3,
           shape = 24, 
           fill = "white", 
           stroke = 1, 
           size = 1.5, 
           color = "tomato2") +
  annotate("point", 
           x= pt_slice$k,
           y = min_density_m3,
           shape = 24, 
           fill = "white", 
           stroke = 1, 
           size = 1.5, 
           color = "tomato2") +
  annotate(
    "point",
    x = pt_slice$k,
    y = pt_slice$host_m3_total,
    shape  = 21,
    fill   = "white",
    stroke = 1.2,
    size   = 3,
    color  = "tomato2"
  ) +
  annotate(
    "text",
    x = 1.85,
    y = 10,  # slight vertical offset
    label = paste0("R[0] == ", round(R0_target, 2)),
    parse = TRUE,
    color = "tomato2",
    size = 5,
    fontface = "bold"
  )

print(all_r0_plt2)

ggsave("plots/all_r0_plot.png", plot = all_r0_plt2, width = 6, height = 5, 
       dpi = 300, 
       bg = "white")


#--------------------------------------
#Combine plots
#----------------------------------------


gf_multi_plot <- as.ggplot(
  rasterGrob(readPNG("plots/greenfrog_multi_r0_plot.png"),
             interpolate = TRUE)
)

bf_multi_plot <- as.ggplot(
  rasterGrob(readPNG("plots/bullfrog_multi_r0_plot.png"),
             interpolate = TRUE)
)

all_r0_plot <- as.ggplot(
  rasterGrob(readPNG("plots/all_r0_plot.png"),
             interpolate = TRUE)
)

# universal theme
pub_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    axis.ticks = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      size = 12,
      face = "plain",
      margin = margin(b = 4)
    ),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    plot.margin = margin(2, 2, 2, 2)
  )

strip_margins <- function(p) {
  p + theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2))
}


all_r0_plot    <- strip_margins(all_r0_plot)
gf_multi_plot <- strip_margins(gf_multi_plot)
bf_multi_plot <- strip_margins(bf_multi_plot)

pgrid2 <- cowplot::plot_grid(
  bf_multi_plot, gf_multi_plot, all_r0_plot,
  ncol = 3,
  labels = c("A", "B", "C"),
  label_size = 12,
  label_fontface = "bold",
  label_x = 0.03,
  label_y = 0.95,
  hjust = 0,
  vjust = 1
)

left_margin   <- 0.08
bottom_margin <- 0.10
right_margin  <- 0.03
top_margin    <- 0.03

final2 <- ggdraw() +
  draw_plot(
    pgrid2,
    x = left_margin,
    y = bottom_margin,
    width  = 1 - left_margin - right_margin,
    height = 1 - bottom_margin - top_margin
  ) +
  annotate(
    "text",
    x = left_margin + (1 - left_margin - right_margin)/2,
    y = bottom_margin * 0.72,
    size = 6,
    hjust = 0.5,
    label = "paste('Scaled movement distance (', k %.% frac(m^3, 6~plain('days')), ')')",
    parse = TRUE
  ) +
  annotate(
    "text",
    x = left_margin * 0.65,
    y = bottom_margin + (1 - bottom_margin - top_margin)/2,
    angle = 90,
    size = 6,
    hjust = 0.5,
    label = "paste('Total host density (per ', m^3, ')')",
    parse = TRUE
  )
print(final2)

ggsave(
  "plots/R0_all_multi_final.png",
  final2,
  width = 12,
  height = 5,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)




