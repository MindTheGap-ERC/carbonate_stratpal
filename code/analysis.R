library(admtools)
library(StratPal)
library(configr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tomledit)
library(ggpubr)
tag1 = "platform"
tag2 = "ramp"

adm_data_pl = read.csv(paste0("data/", tag1, "_adm.csv"))
wd_data_pl = read.csv(paste0("data/", tag1, "_wd.csv"))
adm_data_ra = read.csv(paste0("data/", tag2, "_adm.csv"))
wd_data_ra = read.csv(paste0("data/", tag2, "_wd.csv"))
t = adm_data_pl$time..Myr.
t_steps = adm_data_pl$timestep...

#### Metadata ####
meta_data_pl = tomledit::read_toml(paste0("data/", tag1, ".toml"))
l = tomledit::from_toml(meta_data_pl)
n_locations = length(l$locations)
loc_list = list()
distances = c()
for (loc in seq_len(n_locations)){ # convert to km
  loc_list[[loc]] = list(x = l$locations[[loc]]$x[[1]]/ 1000,
                         y = l$locations[[loc]]$y[[1]] / 1000)
  distances =  c(distances, l$locations[[loc]]$x[[1]]/ 1000)
}

#### Age-depth models ####
adm_list_pl = list()
adm_list_ra = list()
for (i in 1:(length(adm_data_pl)-2)){
  adm_list_pl[[i]] = tp_to_adm(t = adm_data_pl$time..Myr.,
                            h = adm_data_pl[, paste0("adm_", i, "..m.")])
  
  adm_list_ra[[i]] = tp_to_adm(t = adm_data_ra$time..Myr.,
                               h = adm_data_ra[, paste0("adm_", i, "..m.")])
}

#### Water depth ####
wd_pl = list()
wd_ra = list()
for (i in 1:(length(wd_data_pl)-2)){
  wd_pl[[i]] = list(t = wd_data_pl$time..Myr.,
                 wd = pmax(wd_data_pl[, paste0("wd_", i, "..m.")], 0))
  wd_ra[[i]] = list(t = wd_data_ra$time..Myr.,
                    wd = pmax(wd_data_ra[, paste0("wd_", i, "..m.")], 0))
}

#### sea level ####
period1 = 2.0 # Myr
amplitude1 = 15.0 # m
period2 = 0.2 # Myr
amplitude2 = 2.5 # m
sl = list(t = t,
          sl = amplitude1 * sin(2 * pi * t/ period1) + amplitude2 * sin(2 * pi * t / period2))

#### Base plots: SL, adms, wd ####
path_base_fig = "figs/base/"
if (!dir.exists(path_base_fig)) dir.create(path_base_fig, recursive = TRUE)
df = data.frame(t = sl$t, sl = sl$sl)
df |> ggplot(aes(x = t, y = sl)) +
  geom_line()
ggsave(filename = paste0(path_base_fig, "sl.png"))

wd_comb = list(ra = wd_ra, pl = wd_pl)
adm_comb = list(ra = adm_list_ra, pl = adm_list_pl)
cases = c("ra", "pl")
for (i in seq_len(n_locations)){
  for (case in cases){
  print(case)
  wd = wd_comb[[case]][[i]]
  adm = adm_comb[[case]][[i]]
  loc = distances[i]
  png(filename = paste0(path_base_fig, "wd_time_", loc, "km_", case,  ".png"))
  plot(wd$t, wd$wd, type = "l")
  title(main = paste0(case, " ", loc, "km"))
  dev.off()
  png(filename = paste0(path_base_fig, "adm_", loc, "km" , case, ".png"))
  plot(adm, lty_destr = 0)
  title(main = paste0(case, " ", loc, "km"))
  dev.off()
  png(filename = paste0(path_base_fig, "wd_strat_", loc, "km_", case,  ".png"))
  wd_loc = list(t = wd$t, y = wd$wd) |>
    time_to_strat(adm) |>
    plot(type = "l")
  title(main = paste0(case, " ", loc, "km"))
  dev.off()
  }
}



## water depth in the stratigraphic domain
col_names = c("wd_strat", "h", "dist", "tag")
df_wdstr = data.frame(matrix(ncol = 4, nrow = 0))
names(df_wdstr) = col_names

for (case in cases){
  for (i in seq_len(n_locations)){
    wd = wd_comb[[case]][[i]]
    adm = adm_comb[[case]][[i]]
    wd_loc = list(t = wd$t, y = wd$wd) |>
      time_to_strat(adm)
    h = wd_loc$h - min(wd_loc$h, na.rm = TRUE)
    h = h/max(h, na.rm = TRUE)
    df_temp = data.frame(wd_strat = wd_loc$y,
                         h = h,
                         dist = rep(distances[i], length(wd$t)),
                         tag = rep(case, length(wd$t)))
    df_wdstr = rbind(df_wdstr, df_temp)
  }
}

df_wdstr$dist = factor(df_wdstr$dist, levels = distances)  

df_wdstr |>
  filter(tag == "ra") |>
  ggplot( aes(x = h, y = wd_strat, color = dist)) + 
  geom_line()

ggsave(paste0(path_base_fig, "sl_strat_ramps_comp.png"))


df_wdstr |>
  filter(tag == "pl") |>
  ggplot( aes(x = h, y = wd_strat, color = dist)) + 
  geom_line()

ggsave(paste0(path_base_fig, "sl_strat_pl_comp.png"))


#### Plot extinction scenarios ####
LST_rate = approxfun(x = c(1.25, 1.5, 1.75), y = c(1, 25,1), rule = 2)
HST_rate = approxfun(x = c(2.25, 2.5, 2.75), y = c(1, 25,1), rule = 2)
TST_rate = approxfun(x = c(1.75, 2, 2.25), y = c(1, 25,1), rule = 2)
RST_rate = approxfun(x = c(0.75, 1, 1.25, 2.75, 3, 3.25), y = c(1, 25,1, 1, 25, 1), rule = 2)

get_lo_ab = function(case, adm, n = NULL){
  # last occurrences wit abundant fossils & no niche pref
  if (case == "const"){x = p3(from = min_time(adm), to = max_time(adm), rate = 1, n = n)}
  if (case == "LST"){x = p3_var_rate(x = LST_rate, from = min_time(adm), to = max_time(adm),  n = n, f_max = 100)}
  if (case == "HST"){x = p3_var_rate(x = HST_rate, from = min_time(adm), to = max_time(adm),  n = n, f_max = 100)}
  if (case == "TST"){x = p3_var_rate(x = TST_rate, from = min_time(adm), to = max_time(adm),  n = n, f_max = 100)}
  if (case == "RST"){x = p3_var_rate(x = LST_rate, from = min_time(adm), to = max_time(adm), n = n, f_max = 100)}
  y = x |> time_to_strat(adm, destructive = FALSE)
  return(y)
}

path_ex_fig = "figs/extinctions/"
if (!dir.exists(path_ex_fig)) dir.create(path_ex_fig, recursive = TRUE)
binwidth = 3 # width of sampling bins
n_LO = 1000 # no of last occ
sed_cases = c("ra", "pl")
cases = c("const", "HST", "TST", "RST", "LST")
for (i in seq_len(n_locations)){
  for (sed_case in sed_cases){
  adm = adm_comb[[sed_case]][[i]]
  pos = distances[i]
  plot(adm)
  for (case in cases){
    height = get_lo_ab(case = case,
                       adm = adm,
                       n = n_LO)
    data.frame(height) |> ggplot(aes(x = height)) +
      geom_histogram(binwidth = binwidth) +
      ylab("Height [m]") +
      xlab("# LO") +
      coord_flip() +
      ggtitle(paste0("# LO ", pos, " km from shore ", sed_case ))
    ggsave(filename = paste0(path_ex_fig, case, "_", pos, "km", sed_case, ".png"))
  }
  }
}


#### Sed Stats ####
path_sed_stats = "figs/sed_stats/"
if (!dir.exists(path_sed_stats)) dir.create(path_sed_stats, recursive = TRUE)

inc = c()
height = c()
no_hiat = c()
median_hiat = c()
min_hiat = c()
max_hiat = c()
case = c()
for (sed_case in sed_cases){
  adm_list = adm_comb[[sed_case]]
  for (i in seq_along(adm_list)){
    inc = c(inc, adm_list[[i]] |> get_incompleteness())
    height = c(height, adm_list[[i]] |> admtools::max_height())
    no_hiat = c(no_hiat, adm_list[[i]] |> get_hiat_no())
    median_hiat = c(median_hiat, adm_list[[i]] |> get_hiat_duration() |> median())
    min_hiat = c(min_hiat, adm_list[[i]] |> get_hiat_duration() |> min())
    max_hiat = c(max_hiat, adm_list[[i]] |> get_hiat_duration() |> max())
    case = c(case, sed_case)
  }
}
dist = rep(distances, 2)
median_hiat[is.na(median_hiat)] = 0
min_hiat[is.infinite(min_hiat)] = 0
max_hiat[is.infinite(max_hiat)] = 0
case = factor(case, levels = sed_cases)
#dist = factor(dist, levels = distances)

df_sed = data.frame(case = case, dist = dist, inc = inc, height = height, no_hiat = no_hiat,
                    median_hiat = median_hiat, min_hiat = min_hiat, max_hiat = max_hiat)

ggplot(df_sed, aes(x = dist, y = inc, color = case)) +
  geom_line()  +
  ylim(c(0,1)) + 
  ggtitle(label = "Incompleteness")
ggsave(filename = paste0(path_sed_stats, "completeness.png"))


df_sed |> ggplot(aes(x = dist, y = height, color = case)) +
  geom_line() +
  ggtitle("Accumulated sediment")
ggsave(filename = paste0(path_sed_stats, "height.png"))

df_sed |> ggplot(aes(x = dist, y = no_hiat, group = case, color = case)) +
  geom_line() +
  ggtitle("Number of hiatuses")
ggsave(filename = paste0(path_sed_stats, "no_hiatuses.png") )

df_sed |> ggplot(aes(x = dist, y = max_hiat, group = case, color = case)) +
  geom_line() +
  ggtitle("Max hiat duration")
ggsave(filename = paste0(path_sed_stats, "max_hiat.png") )

df_sed |> ggplot(aes(x = dist, y = min_hiat, group = case, color = case)) +
  geom_line() +
  ggtitle("Min hiat duration")
ggsave(filename = paste0(path_sed_stats, "min_hiat.png") )
  

df_sed |> ggplot(aes(x = dist, y = median_hiat, group = case, color = case)) +
  geom_line() +
  ggtitle("Median hiat duration")
ggsave(filename = paste0(path_sed_stats, "median_hiat.png") )


## ridgeline plots
names = c("distance", "system", "hiat_duration")
df = data.frame( matrix(ncol = 0, nrow = length(names)))
names(df) = df
for (i in seq_along(distances)){
  hiat_dur = adm_list_pl[[i]] |> get_hiat_duration()
  df_t = data.frame(distance = rep(distances[i], length(hiat_dur)),
                    system = rep("platform", length(hiat_dur)),
                    hiat_duration = hiat_dur)
  df = rbind(df, df_t)
  hiat_dur = adm_list_ra[[i]] |> get_hiat_duration()
  df_t = data.frame(distance = rep(distances[i], length(hiat_dur)),
                    system = rep("ramp", length(hiat_dur)),
                    hiat_duration = hiat_dur)
  df = rbind(df, df_t)
}
df$distance = factor(df$distance, levels = distances)

pos_interest = c(3,6,9, 12, 15, 18, 21)

p1 = df |> filter(system == "platform" & distance %in% pos_interest) |>
  ggplot(aes(x = hiat_duration, y = distance, fill = distance)) +
  geom_density_ridges() +
  theme_ridges() +
  scale_x_log10() +
  ggtitle("Platform")

p2 = df |> filter(system == "ramp" & distance %in% pos_interest) |>
  ggplot(aes(x = hiat_duration, y = distance, fill = distance)) +
  geom_density_ridges() +
  theme_ridges() +
  scale_x_log10() +
  ggtitle("Ramp")

p3 = ggpubr::ggarrange(p1, p2, nrow = 1, ncol = 2, common.legend = TRUE)
p3

ggsave("figs/hiatus_duration_comp.png", p3)


#### Niches ####
path_niches = "figs/niches/"
if (!dir.exists(path_niches)) dir.create(path_niches, recursive = TRUE)
n_niches = 100
min_depth = max(min(unlist(wd)), 0)
max_depth = max(unlist(wd))
depth_range <- seq(min_depth, max_depth, length.out = 100)
uniform_seq <- seq(0, 1, length.out = n_niches)
optima <- max_depth * (uniform_seq)^2
min_width <- 2   
max_width <- 10  
niche_widths <- min_width + (optima / max(optima)) * (max_width - min_width)

niche_list <- list()
for (i in 1:n_niches) {
  niche_list[[i]] <- StratPal::snd_niche(opt = optima[i],
                                     tol = niche_widths[i],
                                     cutoff_val = 0)
}
niche_val = list(niche_list)
for (i in 1:100) {
  niche = niche_list[[i]]
  niche_val[[i]]=niche(depth_range)
}
niche_mat <- do.call(cbind, niche_val)
png(filename = paste0(path_niches, "all_niches_time.png"))
matplot(depth_range, niche_mat, type = "l",
        lty = 1, col = rgb(0,0,1,0.2),
        xlab = "Water depth [m]",
        ylab = "Niche value",
        main = "Niche optima across different depths")
dev.off()
png(filename = paste0(path_niches, "niche_optima_time.png"))
hist(optima, breaks = 20, col = "blue", border = "white",
     main = "Distribution of Niche Optima Across Depth",
     xlab = "Water Depth (m)", ylab = "Number of Niches")
dev.off()


for (i in seq_len(n_locations)){
  adm = adm_list[[i]]
  pos = distances[i]
  gc = approxfun(x = wd[[i]]$t, y = wd[[i]]$wd)

  # all niches in the time domain
  niches_fossils_t <- do.call(rbind, lapply(seq_along(niche_list), function(j) {
    x <- p3(rate = 300, from = min(t), to = max(t)) |> 
      apply_niche(niche_def = niche_list[[j]], gc = gc)
    if (length(x) > 0 && any(!is.na(x))) {
      data.frame(t = x, niche = as.factor(j))
    } else {
      # Insert a row with NA to ensure this niche appears in the plot/legend
      data.frame(t = NA, niche = as.factor(j))
    }
  }))
  
  ggplot(niches_fossils_t, aes(x = t, fill = niche, color = niche)) +
    geom_density(alpha = 0.12, color = "black") +
    scale_fill_discrete(drop = FALSE) +
    scale_color_discrete(drop = FALSE) +
    labs(
      title = "Probability for fossil presence",
      x = "Time [Myr]",
      y = "Density probability"
    )
  ggsave(filename = paste0(path_niches, "niche_time_", pos, "km.png" ))
  
  # all niches in the stratigraphic domain
  niches_fossils_strat <- do.call(rbind, lapply(seq_along(niche_list), function(j) {
    x <- p3(rate = 300, from = min(t), to = max(t)) |> 
      apply_niche(niche_def = niche_list[[j]], gc = gc) |>                    
      time_to_strat(adm, destructive = TRUE, out_dom_val_h = "default")
    if (length(x) > 0 && any(!is.na(x))) {
      data.frame(t = x, niche = as.factor(j))
    } else {
      # Insert a row with NA to ensure this niche appears in the plot/legend
      data.frame(t = NA, niche = as.factor(j))
    }
  }))
  
  ggplot(niches_fossils_strat, aes(x = t, fill = niche, color = niche)) +
    geom_density(alpha = 0.12, color = "black") +
    scale_fill_discrete(drop = FALSE) +
    scale_color_discrete(drop = FALSE) +
    labs(
      x = "Stratigraphic position [m]",
      y = "Density probability"
    )
  ggsave(filename = paste0(path_niches, "niche_strat_", pos, "km.png" ))
}


#### Time content of section ####
path_time_cont = "figs/time_cont/"
if(!dir.exists(path_time_cont)){
  dir.create(path_time_cont, recursive = TRUE)
}

for (i in 1:n_locations){
  adm = adm_list[[i]]
  plot(adm)
  h_list = get_hiat_list(adm)
  durs = get_hiat_duration(adm)
  png(paste0(path_time_cont, "hiat_dur", i, ".png"))
  plot(NULL, 
       xlim  = c(0, max(durs)),
       ylim = c(admtools::min_height(adm), admtools::max_height(adm)),
       xlab = "Hiatus duration [Myr]",
       ylab = "Stratigraphic position [m]",
       main = paste0("Hiatus duration at pos ", i))
  for (j in seq_along(durs)){
    points(x = durs[j],
           y = h_list[[j]]["height"])
    lines(x = c(0, durs[j]),
          y = c( rep( h_list[[j]]["height"], 2)))
  }
  dev.off()
}


## sed rates ##
for (i in 1:n_locations){
  png(paste0(path_time_cont, "condensation", i, ".png"))
  adm = adm_list[[i]]
  f = condensation_fun(adm)
  h = seq(admtools::min_height(adm), admtools::max_height(adm), by = 0.01)
  plot(x = f(h),
       y = h,
       type = "l",
       xlab = "Condensation [Myr/m]",
       ylab = "Stratigraphic Height [m]",
       main = "Condensation")
  dev.off()
}

#### range offset ####
adm = adm_comb[["pl"]][[1]]
times_ext = seq(0, max_time(adm), by = 0.01)
range_offset_t = rep(NA, length(times_ext))
range_offset_h = rep(NA, length(times_ext))
for (i in seq_along(times_ext)){
  x = range_offset(t_ext = times_ext[i], rate =  100, adm)
  range_offset_h[i] = x["h"]
  range_offset_t[i] = x["t"]
}

plot(times_ext, range_offset_t, type = "l")
plot(times_ext, range_offset_h, type = "l")


#### condensation ratio plot ####
df = data.frame(matrix(nrow = 0, ncol = 4))
names = c("cond_ratio", "proportion_height", "distance", "system")
names(df) = names

bin_width_m = 2
distances_km  = paste(distances, "km")

for (i in seq_along(distances_km)){
  adm = adm_list_ra[[i]]
  h = seq(admtools::min_height(adm), admtools::max_height(adm), by = bin_width_m)  
  adm_2 = tp_to_adm(t = c(admtools::min_time(adm), admtools::max_time(adm)),
                    h = c(admtools::min_height(adm), admtools::max_height(adm)))
  t1 = strat_to_time(h, adm) |> diff()
  t2 = strat_to_time(h, adm_2) |> diff()
  cond_ratio = t1/t2
  proportion_height = ((h[1:(length(h)-1)] +  h[2:length(h)])/2) - admtools::min_height(adm)
  proportion_height = proportion_height/(admtools::max_height(adm)- admtools::min_height(adm))
  df2 = data.frame(cond_ratio = cond_ratio,
                   proportion_height = proportion_height,
                   distance = rep(distances_km[i], length(cond_ratio)),
                   system = rep("ramp", length(cond_ratio)))
  df = rbind(df, df2)
  adm = adm_list_pl[[i]]
  h = seq(admtools::min_height(adm), admtools::max_height(adm), by = bin_width_m)  
  adm_2 = tp_to_adm(t = c(admtools::min_time(adm), admtools::max_time(adm)),
                    h = c(admtools::min_height(adm), admtools::max_height(adm)))
  t1 = strat_to_time(h, adm) |> diff()
  t2 = strat_to_time(h, adm_2) |> diff()
  cond_ratio = t1/t2
  proportion_height = ((h[1:(length(h)-1)] +  h[2:length(h)])/2) - admtools::min_height(adm)
  proportion_height = proportion_height/(admtools::max_height(adm)- admtools::min_height(adm))
  df2 = data.frame(cond_ratio = cond_ratio,
                   proportion_height = proportion_height,
                   distance = rep(distances_km[i], length(cond_ratio)),
                   system = rep("platform", length(cond_ratio)))
  df = rbind(df, df2)
}

df$distance = factor(df$distance, levels = distances_km)
df$system = factor(df$system)

pos_interest = paste(c(3,6,9, 12, 15, 18, 21), "km")
y_lim = range(df$cond_ratio)
y_lim[1] = 0.95 * y_lim[1]
y_lim[2] = 1.1 * y_lim[2]
title_ramp = "Ramp Geometry"
title_platform = "Platform Geometry"
y_label = "Inflation of Rates [-]"
x_label = "Relative Height [-]"
legend_label = "Distance from Shore"
p1 = df |> 
  filter(system == "ramp" &
           distance %in% pos_interest) |>
  ggplot(aes(x = proportion_height, y = cond_ratio, color = distance)) +
  geom_line()  +
  coord_flip() +
  scale_y_log10(limits = y_lim) +
  labs(y = y_label,
       x = x_label,
       title = title_ramp,
       color = legend_label)
p1

p2 = df |> 
  filter(system == "platform" &
           distance %in% pos_interest) |>
  ggplot(aes(x = proportion_height, y = cond_ratio, color = distance)) +
  geom_line() +
  coord_flip() +
  scale_y_log10(limits = y_lim) +
  labs(y = y_label,
       x = x_label,
       title = title_platform,
       color = legend_label)
p2

p3 = ggpubr::ggarrange(p2, p1,  nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom", labels = LETTERS[1:2])
p3
ggsave("figs/inflation_of_rates.png", plot = p3)

#### AGe-depth model plots ####
names = c("time", "height", "distance", "system")
df = data.frame(matrix(nrow = 0, ncol = length(names)))
names(df) = names

title_platform = "Platform Geometry"
title_ramp = "Ramp Geometry"
x_lab_title = "Elapsed Model Time [Myr]"
y_lab_title = "Height [m]"
legend_label = "Distance from Shore"

distances_km = paste(distances, "km")

for (i in seq_along(distances_km)){
  adm = adm_list_ra[[i]]
  t = get_T_tp(adm)
  h = time_to_strat(t, adm, destructive = TRUE)
  df2 = data.frame(time = t,
                   height = h,
                   distance = rep(distances_km[i], length(time)),
                   system = rep("ramp", length(time)))
  df = rbind(df, df2)
  adm = adm_list_pl[[i]]
  t = get_T_tp(adm)
  h = time_to_strat(t, adm, destructive = TRUE)
  df2 = data.frame(time = t,
                   height = h,
                   distance = rep(distances_km[i], length(time)),
                   system = rep("platform", length(time)))
  df = rbind(df, df2)
}

df$distance = factor(df$distance, levels = distances_km)
df$system = factor(df$system)

pos_interest = paste(c(3,6,9, 12, 15, 18, 21), "km")

p1 = df |>
  filter(system == "ramp" &
           distance %in% pos_interest) |>
  ggplot(aes(x = time, y = height, color = distance)) +
  geom_line() +
  labs(title = title_ramp,
       x = x_lab_title,
       y = y_lab_title,
       col = legend_label)
p1
p2 = df |>
  filter(system == "platform" &
           distance %in% pos_interest) |>
  ggplot(aes(x = time, y = height, color = distance)) +
  geom_line() +
  labs(title = title_platform,
       x = x_lab_title,
       y = y_lab_title,
       col = legend_label)
p2
p3 = ggpubr::ggarrange(p1, p2, nrow = 1, ncol = 2,
                       common.legend = TRUE,
                       legend = "bottom",
                       labels = LETTERS[1:2])
p3

ggsave("figs/age_depth_overview.png", p3)

#### Water depth overview ####
## time domain ##
title_platform = "Platform Geometry"
title_ramp = "Ramp Geometry"
x_lab_title = "Elapsed Model Time [Myr]"
y_lab_title = "Water Depth [m]"
legend_label = "Distance from Shore"

col_names = c("wd", "t", "dist", "tag")
df_wd = data.frame(matrix(ncol = length(col_names), nrow = 0))
names(df_wd) = col_names
for (case in cases){
  for (i in seq_len(n_locations)){
    wd = wd_comb[[case]][[i]]
    df_temp = data.frame(wd = wd$wd,
                         t = wd$t,
                         dist = rep(distances_km[i], length(wd$t)),
                         tag = rep(case, length(wd$t)))
    df_wd = rbind(df_wd, df_temp)
  }
}

df_wd$dist = factor(df_wd$dist, levels = distances_km)

p1 = df_wd |>
  filter(tag == "ra" & dist %in% pos_interest) |>
  ggplot( aes(x = t, y = wd, color = dist)) + 
  geom_line() +
  labs(title = title_ramp,
       x = x_lab_title,
       y = y_lab_title,
       col = legend_label)
p1

p2 = df_wd |>
  filter(tag == "pl" & dist %in% pos_interest) |>
  ggplot( aes(x = t, y = wd, color = dist)) + 
  geom_line() +
  labs(title = title_platform,
       x = x_lab_title,
       y = y_lab_title,
       col = legend_label)
p2

p3 = ggpubr::ggarrange(p2, p1, ncol = 2, nrow = 1, labels = LETTERS[1:2], legend = "bottom", common.legend = TRUE)
p3

ggsave("figs/water_depth_comparison_time_domain.png")


## stratigraphic domain
names = c("wd", "proportion_height", "distance", "system")
df = data.frame(matrix(nrow = 0, ncol = length(names)))
names(df) = names

distances_km  = paste(distances, "km")
sampling_distance_m = 1 # m
for (case in names(wd_comb)){
  wd_case = wd_comb[[case]]
  adm_case = adm_comb[[case]]
    for (i in seq_along(distances_km)){
      adm = adm_case[[i]]
      wd = wd_case[[i]]
      max_height = admtools::max_height(adm)
      min_height = admtools::min_height(adm)
      h = seq(min_height, max_height, by = sampling_distance_m)
      t = strat_to_time(h, adm)
      wds = approx(x = wd$t, y = wd$wd, xout = t)$y
      h_prop = h/(max_height - min_height)
      df2 = data.frame(wd = wds, proportion_height = h_prop, distance = rep(distances_km[i], length(h_prop)), system = rep(case, length(h_prop)))
      df = rbind(df, df2)
  }
}

df$distance = factor(df$distance, levels = distances_km)
res_distances = c(3,6,9, 12, 15, 18) |> paste("km")
title_ramp = "Ramp Geometry"
title_platform = "Platform Geometry"
wd_label = "Water Depth [m]"
height_label = "Relative Height [m]"
legend_title = "Distance from Shore"

p1 = df |>
  filter(system == "pl" & distance %in% res_distances) |>
  ggplot(aes(y = wd, x = proportion_height, col = distance)) +
  geom_line() +
  coord_flip() +
  labs(title = title_platform,
       y = wd_label,
       x = height_label,
       col = legend_title)
p1

p2 = df |>
  filter(system == "ra" & distance %in% res_distances) |>
  ggplot(aes(y = wd, x = proportion_height, col = distance)) +
  geom_line() +
  coord_flip() +
  labs(title = title_platform,
       y = wd_label,
       x = height_label,
       col = legend_title)
p2

p3 = ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, labels = LETTERS[1:2], legend = "bottom", common.legend = TRUE)
p3
ggsave("figs/water_depth_strat_domain.png", p3)


#### Last occurrences ####
library(ggridges)
case = "pl"
adm_used = list("top" = adm_comb[[case]][[1]],
                "slope" = adm_comb[[case]][[10]])
rates = c(2, 5, 10, 30, 100, 1000)
df = data.frame(l_occ_h = NULL, rate = NULL, loc = NULL)
for (loc in names(adm_used)){
  adm = adm_used[[loc]]
  for (j in seq_along(rates)){
    rate = rates[j]
    t_ext = p3(1000, from = min_time(adm), to = max_time(adm), n = 1000)
    l_occ_t = rep(NA, length(t_ext))
    l_occ_h = rep(NA, length(t_ext))
    for (i in seq_along(t_ext)){
      x = last_occ(t_ext = t_ext[i], rate, adm = adm)
      l_occ_h[i] = x["h"]
      l_occ_t[i] = x["t"]
    }

    df = rbind(df, data.frame(l_occ_h = l_occ_h, rate = rep(rate, length(l_occ_h)), loc = rep(loc, length(l_occ_h))))
  }
}

df$rate = factor(df$rate, levels = rates)
df$loc = factor(df$loc, levels = c("top", "slope"))
p1 = df |>
  filter(loc == "top") |>
  ggplot(aes(x = l_occ_h, y = rate, fill = rate)) +
  geom_density_ridges(stat = "binline") +
  labs(title = "LOs platform top")
p1

p2 = df |> 
  filter(loc == "slope") |>
  ggplot(aes(x = l_occ_h, y = rate, fill = rate)) +
  geom_density_ridges(stat = "binline") +
  labs(title = "LOs platform slope")
p2

p3 = ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, labels = LETTERS[1:2], common.legend = TRUE, legend = "bottom") 

    df = rbind(df,
               data.frame(l_occ_h = l_occ_h,
                          rate = rep(rate, length(l_occ_h))))
  }
  df$rate = factor(df$rate, levels = rates)
  
  p1 = df |>
    filter( is.finite(l_occ_h)) |>
    ggplot(aes(x = l_occ_h, y = rate, fill = rate)) +
    geom_density_ridges(stat = "binline", bins = 30, scale = 1) 
  
  if (plot_st){
    st_sep_time_mod = c(0,st_sep_time, 4)
    st_sep_strat = time_to_strat(st_sep_time_mod, adm, destructive = FALSE)
    st_pres = !(diff(st_sep_strat) == 0)
    st_names_used = st_labels[st_pres]
    st_sep_strat_used = st_sep_strat[st_pres]
    df_text = data.frame(time = 0.5* (head(st_sep_strat_used, -1) + tail(st_sep_strat_used, -1)),
                         height = rep(2.5, length(st_sep_strat_used)-1),
                         label = st_names_used)
    
    p1 = p1 +
      geom_vline(xintercept = st_sep_strat_used[2:(length(st_sep_strat_used)-1)],
                 linetype = "dashed",
                 color = "black") +
      geom_text(data = df_text,
                aes(x = time,
                    y = height,
                    label = label),
                inherit.aes = FALSE)
  }
  p1 = p1 +
    labs(title = title,
         x = y_axis_label,
         y = x_axis_label,
         fill = legend_title) +
    coord_flip()
  return(p1)
}

p1 = plot_lo_by_rate(case = "ra",
                     pos = 12, 
                     rates = c(3,10,30,100),
                     title = 'Platform top')

plot_lo_comparison = function(){
  rates = c(3,10,30,100)
  p1 = plot_lo_by_rate(case = "pl",
                       pos = 3, 
                       rates = rates,
                       title = 'Platform top')
  p2 = plot_lo_by_rate(case = "pl",
                       pos = 12, 
                       rates = rates,
                       title = 'Platform ramp')
  p3 = ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, labels = LETTERS[1:2], common.legend = TRUE, legend = "bottom") 
  return(p3)
}

p = plot_lo_comparison()
p
ggsave("figs/ms_unknown_fig_no.png", p)

#### Plots: Incompleteness, gap statistics, and gap distribution ####

plot_completeness = function(pos = seq(1.5, 21, by = 1.5)){
  stopifnot(all(pos %in% distances))
  
  names = c("dist", "case", "com")
  df = data.frame( matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
  for (case in names(adm_comb)){
    for (i in seq_along(pos)){
      adm = adm_comb[[case]][[i]]
      com = admtools::get_completeness(adm)
      df3 = data.frame(dist = rep(pos[i], length(com)),
                       case = rep(case, length(com)),
                       com = com)
      df = rbind(df, df3)
    }
  }
  
  p = df |>
  ggplot(aes(x = dist, y = com, color = case)) + 
  geom_line(linewidth = 3) +
  labs(x = "distance from shore",
       y = "completeness",
       color = "geometry") +
  ylim(c(0,1)) +
    labs(x = "Distance from shore [km]",
         y = "Completeness [-]",
         title = "Stratigraphic Completeness",
         color = "Geometry") +
    scale_color_discrete(labels = c("Platform", "Ramp")) +
    theme(legend.position = c(0.1, 0.9))
return(p)
}
p = plot_completeness()
p
ggsave(filename = "figs/completeness.png",
       plot = p)

plot_gap_statistics = function(pos = seq(1.5, 21, by = 1.5)){
  stopifnot(all(pos %in% distances))
  names = c("max", "quant_1", "quant_3", "median", "dist", "case")
  df = data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
  for (case in names(adm_comb)){
    for (i in seq_along(distances_km)){
      adm = adm_comb[[case]][[i]]
      hiat = get_hiat_duration(adm)
      
      df2 = data.frame(
        max = max(hiat),
        quant_1 = quantile(hiat, probs = 0.25, names = FALSE),
        quant_3 = quantile(hiat, probs = 0.75, names = FALSE),
        median = median(hiat),
        n_hiat = length(hiat),
        acc_sed = get_total_thickness(adm),
        dist = distances[i],
        case  = case
      )
      df = rbind(df, df2)
    }
  }
  df$max[is.infinite(df$max)] = 0
  df$quant_1[is.na(df$quant_1)] = 0
  df$quant_3[is.na(df$quant_3)] = 0
  df$median[is.na(df$median)] = 0
  df$case = factor(df$case)
  
  p1 = df |> 
    pivot_longer(cols = c("max", "quant_1", "quant_3", "median"),
                 names_to = "measure",
                 values_to = "value") |>
    filter(case == "ra") |>
    ggplot(aes(x = dist, y = value, color = measure)) +
    geom_line(linewidth = 3) +
    labs(x = "Distance from shore [km]",
         y = "Value [Myr]",
         title = "Ramp",
         color = "Statistic") +
    scale_color_discrete(labels = c("maximum", "median", "1st quantile", "3rd quantile")) +
    theme(legend.position = c(0.9, 0.8)) +
    scale_y_log10()
  p2 = df  |>
    pivot_longer(cols = c("max", "quant_1", "quant_3", "median"),
                 names_to = "measure",
                 values_to = "value") |>
    filter(case == "pl") |>
    ggplot(aes(x = dist, y = value, color = measure)) +
    geom_line(linewidth = 3)  +
    labs(x = "Distance from shore [km]",
         y = "Value [Myr]",
         title = "Platform",
         color = "Statistic") +
    scale_color_discrete(labels = c("maximum", "median", "1st quantile", "3rd quantile")) +
    theme(legend.position = c(0.9, 0.8)) + 
    scale_y_log10()
  p3 = ggpubr::ggarrange(p2, p1, ncol = 2, nrow = 1, labels = LETTERS[1:2],common.legend = TRUE, legend = "bottom")
  return(p3)
}
p = plot_gap_statistics()
p
ggsave(filename = "figs/gap_statistics.png",
       plot = p)

plot_hiat_duration = function(pos = seq(3, 21, by = 3)){
  stopifnot(all(pos %in% distances))
  
  names = c("dist", "case", "hiat_duration")
  df = data.frame( matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
  for (case in names(adm_comb)){
    for (i in seq_along(pos)){
      adm = adm_comb[[case]][[i]]
      hiat = admtools::get_hiat_duration(adm)
      df3 = data.frame(dist = rep(pos[i], length(hiat)),
                       case = rep(case, length(hiat)),
                       hiat_duration = hiat)
      df = rbind(df, df3)
    }
  }
  df$dist = factor(df$dist)
  labels = paste(pos, "km")
  p1 = df |> 
    filter(case == "pl") |>
    ggplot(aes(x = hiat_duration, y = dist, fill = dist)) +
    geom_density_ridges(bandwidth = 0.1) +
    theme_ridges()+
    geom_vline(xintercept = c(0.001, 0.1, 1),
               linetype = "dashed") +
    scale_x_log10() +
    ggtitle("Platform") + 
    theme_classic()+
    labs(x = "Hiatus duration [Myr]",
         y = "Frequency",
         fill = "Distance from shore") +
    scale_fill_discrete(labels = labels) 
  p1
  
  p2 = df |> filter(case == "ra") |>
    ggplot(aes(x = hiat_duration, y = dist, fill = dist)) +
    geom_density_ridges(bandwidth = 0.1) +
    theme_ridges() +
    geom_vline(xintercept = c(0.001, 0.1, 1),
               linetype = "dashed") +
    scale_x_log10() +
    theme_classic() +
    ggtitle("Ramp") +
    labs(x = "Hiatus duration [Myr]",
         y = "Frequency",
         fill = "Distance from shore") +
    scale_fill_discrete(labels = labels)
  p2
  p3 = ggpubr::ggarrange(p1, p2, nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom", labels = LETTERS[1:2])
  return(p3)
}
p = plot_hiat_duration()
p
ggsave("figs/hiatus_duration_comp.png", p)

plot_no_of_hiat = function(pos = seq(1.5, 21, by = 1.5)){
  stopifnot(all(pos %in% distances))
  
  names = c("dist", "case", "n")
  df = data.frame( matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
  for (case in names(adm_comb)){
    for (i in seq_along(pos)){
      adm = adm_comb[[case]][[i]]
      th = admtools::get_hiat_no(adm)
      df3 = data.frame(dist = rep(pos[i], length(th)),
                       case = rep(case, length(th)),
                       n = th)
      df = rbind(df, df3)
    }
  }
  
  p = df |>
    ggplot(aes(x = dist, y = n, color = case)) +
    geom_line(linewidth = 3) +
    labs(title = "Number of Hiatuses",
         y = "Number of hiatuses",
         x = "Distance from shore [km]",
         color = "Geometry") +
    scale_color_discrete(labels = c("Platform", "Ramp")) +
    theme(legend.position = c(0.1,0.9)) +
    ylim(c(0, max(df$n)))
  return(p)
}
p = plot_no_of_hiat()
p
ggsave("figs/no_of_hiatuses.png", p)

plot_accumulated_sediment = function(pos = seq(3, 21, by = 1.5)){
  stopifnot(all(pos %in% distances))
  
  names = c("dist", "case", "thickness")
  df = data.frame( matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
  for (case in names(adm_comb)){
    for (i in seq_along(pos)){
      adm = adm_comb[[case]][[i]]
      th = admtools::get_total_thickness(adm)
      df3 = data.frame(dist = rep(pos[i], length(th)),
                       case = rep(case, length(th)),
                       thickness = th)
      df = rbind(df, df3)
    }
  }
  
  p = df |>
    ggplot(aes(x = dist, y = thickness, color = case)) +
    geom_line(linewidth = 3) +
    labs(title = "Section thickness",
         y = "Section thickness [m]",
         x = "Distance from shore [km]",
         color = "Geometry") +
    scale_color_discrete(labels = c("Platform", "Ramp")) +
    theme(legend.position = c(0.1, 0.9))
  return(p)
}
p = plot_accumulated_sediment()
p
ggsave("figs/section_thickness.png", p)

#### Extinctions by systems tract ####

plot_ext_scenario_comparison = function(rate, dist, case, title){
  stopifnot(dist %in% distances)
  stopifnot(case %in% cases)
  
  n_locc = 1000
  
  names = c("ext_scen", "case", "dist", "l_occ_h")
  df = data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  rate = rate

  for (scen in names(ext_scen)){
        f = ext_scen[[scen]]
        adm = adm_comb[[case]][[which(distances == dist)]]
        t_ext = p3_var_rate(x = f,
                            from = admtools::min_time(adm),
                            to = admtools::max_time(adm),
                            n = n_locc,
                            f_max = 25)
        l_occ_h = rep(NA, length(t_ext))
        for (j in 1:length(t_ext)){
          l_occ = last_occ(t_ext = t_ext[j],
                           rate = rate,
                           adm = adm)
          l_occ_h[j] = l_occ["h"]
        }
        df2 = data.frame(ext_scen = rep(scen, length(l_occ_h)),
                         case = rep(case, length(l_occ_h)),
                         dist = rep(dist, length(l_occ_h)),
                         l_occ_h = l_occ_h)
        df = rbind(df, df2)
  }
  df$case = factor(df$case)
  df$ext_scen = factor(df$ext_scen, levels = c("constant", "HST", "RST", "LST", "TST"),
                       ordered = TRUE)
  
  st_sep_time_mod = c(0,st_sep_time, 4)
  st_sep_strat = time_to_strat(st_sep_time_mod, adm, destructive = FALSE)
  st_pres = !(diff(st_sep_strat) == 0)
  st_names = st_labels
  st_names_used = st_names[st_pres]
  st_sep_strat_used = st_sep_strat[st_pres]
  df_text = data.frame(time = 0.5* (head(st_sep_strat_used, -1) + tail(st_sep_strat_used, -1)),
                       height = rep(2.5, length(st_sep_strat_used)-1),
                       label = st_names_used)
  
  p = df |>
    filter(!is.na(l_occ_h)) |>
    ggplot(aes(x = l_occ_h, y = ext_scen, fill = ext_scen)) +
    geom_density_ridges(stat = "binline",
                        breaks = seq(0,admtools::max_height(adm), length.out = 30),
                        scale = 0.9)+
    geom_vline(xintercept = tail(head(st_sep_strat_used, -1),-1),
               linetype = "dashed",
               color = "black") +
    geom_text(data = df_text,
              aes(x = time,
                  y = height,
                  label = label),
              inherit.aes = FALSE) +
    coord_flip() +
    labs(title = title,
         y = "Extinction scenario",
         x = "Stratigraphic height [m]",
         fill = "Extinction scenario")
  return(p)
}

p1 = plot_ext_scenario_comparison(rate = 10, dist = 15, case = "pl", title = "Platform proximal slope")
p2 = plot_ext_scenario_comparison(rate = 10, dist = 3, case = "pl", title = "Platform top")
p3 = ggpubr::ggarrange(p2,p1, ncol = 2, nrow =1, common.legend = TRUE, legend = "bottom", labels = LETTERS[1:2])

p3
ggsave("figs/last_occ.png", p3)

