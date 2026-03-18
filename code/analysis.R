library(admtools)
library(StratPal)
library(configr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tomledit)
library(ggpubr)
library(ggridges)

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

wd_comb = list(ra = wd_ra, pl = wd_pl)
adm_comb = list(ra = adm_list_ra, pl = adm_list_pl)
cases = c("ra", "pl")

distances_km = paste(distances, "km")

#### Base plots: SL, adms, wd ####
path_base_fig = "figs/base/"
if (!dir.exists(path_base_fig)) dir.create(path_base_fig, recursive = TRUE)
df = data.frame(t = sl$t, sl = sl$sl)
df |> ggplot(aes(x = t, y = sl)) +
  geom_line()
ggsave(filename = paste0(path_base_fig, "sl.png"))


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

#### Niches ####


# plot specific niche
plot_niche_strat = function(){
  
}
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

#### Plot: Age-depth models ####
plot_age_depth_models = function(){
  pos_interest = paste(c(3,6,9, 12, 15, 18, 21), "km")
  st_sep_time = c(-0.2,seq(0.25, 3.75, by = 0.5), 4.2)
  names = c("time", "height", "distance", "system")
  df = data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
  title_platform = "Platform"
  title_ramp = "Ramp"
  x_lab_title = "Elapsed Model Time [Myr]"
  y_lab_title = "Height [m]"
  legend_label = "Distance from Shore"
  distances_km = paste(distances, "km")
  
  for (case in cases){
    for (i in seq_along(distances_km)){
      adm = adm_comb[[case]][[i]]
      t = get_T_tp(adm)
      h = time_to_strat(t, adm, destructive = TRUE)
      df2 = data.frame(time = t,
                       height = h,
                       distance = rep(distances_km[i], length(time)),
                       system = rep(case, length(time)))
      df = rbind(df, df2)
    }
  }
  df$distance = factor(df$distance, levels = distances_km)
  df$system = factor(df$system)
  

  
  max_height = df |> filter(system == "ra") 
  max_height = max_height$height[!is.na(max_height$height)] |> max()
  df_text = data.frame(time = 0.5* (head(st_sep_time, -1) + tail(st_sep_time, -1)),
                      height = rep(max_height, length(st_sep_time)-1),
                      label = c("TST","HST", "RST","LST","TST","HST", "RST", "LST", "TST"))
  p1 = df |>
    filter(system == "ra" &
             distance %in% pos_interest) |>
    ggplot(aes(x = time, y = height, color = distance)) +
    geom_line() +
    labs(title = title_ramp,
         x = x_lab_title,
         y = y_lab_title,
         col = legend_label) +
    geom_vline(xintercept = head(tail(st_sep_time, -1), -1),
          linetype = "dashed",
          color = "black",
          lwd = 0.5) +
    geom_text(data = df_text, aes(x = time,y = height,label = label),inherit.aes = FALSE)
  
  max_height = df |> filter(system == "pl") 
  max_height = max_height$height[!is.na(max_height$height)] |> max()
  df_text = data.frame(time = 0.5* (head(st_sep_time, -1) + tail(st_sep_time, -1)),
                       height = rep(max_height, length(st_sep_time)-1),
                       label = c("TST","HST", "RST","LST","TST","HST", "RST", "LST", "TST"))
  
  p2 = df |>
    filter(system == "pl" &
             distance %in% pos_interest) |>
    ggplot(aes(x = time, y = height, color = distance)) +
    geom_line() +
    labs(title = title_platform,
         x = x_lab_title,
         y = y_lab_title,
         col = legend_label)  +
    geom_vline(xintercept = head(tail(st_sep_time, -1), -1),
               linetype = "dashed",
               color = "black",
               lwd = 0.5) +
    geom_text(data = df_text, aes(x = time,y = height,label = label),inherit.aes = FALSE)
  p3 = ggpubr::ggarrange(p2, p1, nrow = 1, ncol = 2,
                         common.legend = TRUE,
                         legend = "bottom",
                         labels = LETTERS[1:2])
  return(p3)
}
p = plot_age_depth_models()
p
ggsave("figs/age_depth_overview.png", p)


#### Plot: Water depths in time domain ####
plot_wd_time_domain = function(){
  st_sep_time = c(-0.2,seq(0.25, 3.75, by = 0.5), 4.2)
  title_platform = "Platform"
  title_ramp = "Ramp"
  x_lab_title = "Elapsed Model Time [Myr]"
  y_lab_title = "Water Depth [m]"
  legend_label = "Distance from Shore"
  pos_interest = paste(c(3,6,9, 12, 15, 18, 21), "km")
  
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
  df_text = data.frame(time = 0.5* (head(st_sep_time, -1) + tail(st_sep_time, -1)),
                       height = rep(125, length(st_sep_time)-1),
                       label = c("TST", "HST", "RST","LST","TST","HST", "RST", "LST", "TST"))
  
  p1 = df_wd |>
    filter(tag == "ra" & dist %in% pos_interest) |>
    ggplot( aes(x = t, y = wd, color = dist)) + 
    geom_line() +
    labs(title = title_ramp,
         x = x_lab_title,
         y = y_lab_title,
         col = legend_label)  +
    geom_vline(xintercept = head(tail(st_sep_time, -1),-1),
               linetype = "dashed",
               color = "black",
               lwd = 0.5) +
    geom_text(data = df_text,
              aes(x = time,
                  y = height,
                  label = label),
              inherit.aes = FALSE)
  p2 = df_wd |>
    filter(tag == "pl" & dist %in% pos_interest) |>
    ggplot( aes(x = t, y = wd, color = dist)) + 
    geom_line() +
    labs(title = title_platform,
         x = x_lab_title,
         y = y_lab_title,
         col = legend_label) +
    geom_vline(xintercept = head(tail(st_sep_time, -1),-1),
               linetype = "dashed",
               color = "black",
               lwd = 0.5) +
    geom_text(data = df_text,
              aes(x = time,
                  y = height,
                  label = label),
              inherit.aes = FALSE)
  p3 = ggpubr::ggarrange(p2, p1, ncol = 2, nrow = 1, labels = LETTERS[1:2], legend = "bottom", common.legend = TRUE)
  return(p3)
}
p = plot_wd_time_domain()
p
ggsave("figs/water_depth_comparison_time_domain.png", p)

#### Plot: Water depth in the stratigraphic domain ####
plot_wd_strat_domain = function(){
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
  res_distances = c(3,6,9, 12, 15, 18, 21) |> paste("km")
  title_ramp = "Ramp"
  title_platform = "Platform"
  wd_label = "Water Depth [m]"
  height_label = "Relative Height [-]"
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
  p2 = df |>
    filter(system == "ra" & distance %in% res_distances) |>
    ggplot(aes(y = wd, x = proportion_height, col = distance)) +
    geom_line() +
    coord_flip() +
    labs(title = title_ramp,
         y = wd_label,
         x = height_label,
         col = legend_title)
  p3 = ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, labels = LETTERS[1:2], legend = "bottom", common.legend = TRUE)
  return(p3)
}

p = plot_wd_strat_domain()
p
ggsave("figs/water_depth_strat_domain.png", p)

#### Plot: condensation ratio ####

plot_condensation_ratio = function(){
  df = data.frame(matrix(nrow = 0, ncol = 4))
  names = c("cond_ratio", "proportion_height", "distance", "system")
  names(df) = names
  bin_width_m = 4
  distances_km  = paste(distances, "km")
  for (case in names(adm_comb)){
    for (i in seq_along(distances_km)){
      adm = adm_comb[[case]][[i]]
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
                       system = rep(case, length(cond_ratio)))
      df = rbind(df, df2)
    }
  }
  df$distance = factor(df$distance, levels = distances_km)
  df$system = factor(df$system)
  
  pos_interest = paste(c(3,6,9,12,15,18,21), "km")
  y_lim = range(df$cond_ratio)
  y_lim[1] = 0.95 * y_lim[1]
  y_lim[2] = 1.1 * y_lim[2]
  title_ramp = "Ramp"
  title_platform = "Platform"
  y_label = "Condensation Ratio [-]"
  x_label = "Relative Height [-]"
  legend_label = "Distance from Shore"
  p1 = df |> 
    filter(system == "ra" &
             distance %in% pos_interest) |>
    ggplot(aes(x = proportion_height, y = cond_ratio, color = distance)) +
    geom_line() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    coord_flip() +
    scale_y_log10(limits = y_lim) +
    labs(y = y_label,
         x = x_label,
         title = title_ramp,
         color = legend_label)
  p2 = df |> 
    filter(system == "pl" &
             distance %in% pos_interest) |>
    ggplot(aes(x = proportion_height, y = cond_ratio, color = distance)) +
    geom_line() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    coord_flip() +
    scale_y_log10(limits = y_lim) +
    labs(y = y_label,
         x = x_label,
         title = title_platform,
         color = legend_label)
  p3 = ggpubr::ggarrange(p2, p1,  nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom", labels = LETTERS[1:2])
  return(p3)
}

p = plot_condensation_ratio()
p
ggsave("figs/condensation_ratio.png", plot = p)

#### Plot: Cumulative distribution of time vs height ####
plot_prop_time_vs_height = function(){
  names = c("prop_time", "prop_height", "dist", "case")
  df = data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
  x_label = "Proportion cumulative time [-]"
  y_label = "Proportion cumulative height [-]"
  title_platform = "Platform"
  title_ramp = "Ramp"
  legend_label = "Distance from shore"
  
  for (case in names(adm_comb)){
    for (i in seq_along(distances_km)){
      adm = adm_comb[[case]][[i]]
      r = diff(get_L_tp(adm))/diff(get_T_tp(adm)) 
      r = r |> sort() |> cumsum()
      t_rel = seq_along(r)/length(r)
      df2 = data.frame(prop_time = t_rel,
                       prop_height = r/max(r),
                       dist = rep(distances_km[i], length(t_rel)),
                       system = rep(case, length(t_rel)))
      df = rbind(df, df2)
    }
  }
  df$dist = factor(df$dist, levels = distances_km)
  df$system = factor(df$system)
  pos_interest = paste(c(3,6,9,12,15,18,21), "km")
  
  p1 = df |>
    filter(system == "pl" & dist %in% pos_interest) |>
    ggplot(aes(x = prop_time, y = prop_height, color = dist)) +
    geom_line() +
    labs(y = y_label,
         x = x_label,
         title = title_platform,
         color = legend_label) + 
    geom_segment(x = 0,
                 y = 0,
                 xend = 1,
                 yend = 1,
                 color = "black",
                 linetype = "dashed")
  p2 = df |>
    filter(system == "ra" & dist %in% pos_interest) |>
    ggplot(aes(x = prop_time, y = prop_height, color = dist)) +
    geom_line() +
    labs(y = y_label,
         x = x_label,
         title = title_ramp,
         color = legend_label) + 
    geom_segment(x = 0,
                 y = 0,
                 xend = 1,
                 yend = 1,
                 color = "black",
                 linetype = "dashed")
  p3 = ggpubr::ggarrange(p1, p2,  nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom", labels = LETTERS[1:2])
}
p = plot_prop_time_vs_height()
p
ggsave("figs/prop_time_vs_height.png", plot = p)

#### Plot: Abundance bias on last occurrences ####
plot_lo_comparison = function(case = "pl"){
  adm_used = list("top" = adm_comb[[case]][[1]],
                  "slope" = adm_comb[[case]][[10]])
  rates = c(2, 5, 10, 30, 100, 300)
  names = c("l_occ_h","rate", "loc")
  df = data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  if (case == "pl"){
    title_top = "Platform top"
    title_slope = "Platform slope"
  }
  if (case == "ra"){
    title_top = "Ramp top"
    title_slope = "Ramp slope"
  }

  y_axis_label = "Stratigraphic height [m]"
  legend_title = "Fossil abundance [#/Myr]"
  x_axis_label = "Last occurrences"
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
      }
      df = rbind(df, data.frame(l_occ_h = l_occ_h,
                                rate = rep(rate, length(l_occ_h)), 
                                loc = rep(loc, length(l_occ_h))))
    }
  }
  
  df$rate = factor(df$rate, levels = rates)
  df$loc = factor(df$loc, levels = c("top", "slope"))
  st_sep_time = c(0,seq(0.25, 3.75, by = 0.5), 4)
  st_sep_strat = time_to_strat(st_sep_time, adm_used[["top"]], destructive = FALSE)
  st_pres = !(diff(st_sep_strat) == 0)
  st_names = c("TST", "HST", "RST","LST","TST","HST", "RST", "LST", "TST")
  st_names_used = st_names[st_pres]
  st_sep_strat_used = st_sep_strat[st_pres]
  df_text = data.frame(time = 0.5* (head(st_sep_strat_used, -1) + tail(st_sep_strat_used, -1)),
                       height = rep(2.5, length(st_sep_strat_used)-1),
                       label = st_names_used)
  p1 = df |>
    filter(loc == "top" & is.finite(l_occ_h)) |>
    ggplot(aes(x = l_occ_h, y = rate, fill = rate)) +
    geom_density_ridges(stat = "binline", bins = 30, scale = 1) +
    labs(title = title_top,
         x = y_axis_label,
         y = x_axis_label,
         fill = legend_title) +
    geom_vline(xintercept = st_sep_strat_used[2:(length(st_sep_strat_used)-1)],
               linetype = "dashed",
               color = "black") +
    geom_text(data = df_text,
              aes(x = time,
                  y = height,
                  label = label),
              inherit.aes = FALSE) +
    coord_flip()
  
  st_sep_time = c(0,seq(0.25, 3.75, by = 0.5), 4)
  st_sep_strat = time_to_strat(st_sep_time, adm_used[["slope"]], destructive = FALSE)
  st_pres = !(diff(st_sep_strat) == 0)
  st_names = c("TST", "HST", "RST","LST","TST","HST", "RST", "LST", "TST")
  st_names_used = st_names[st_pres]
  st_sep_strat_used = st_sep_strat[st_pres]
  df_text = data.frame(time = 0.5* (head(st_sep_strat_used, -1) + tail(st_sep_strat_used, -1)),
                       height = rep(2.5, length(st_sep_strat_used)-1),
                       label = st_names_used)
  p2 = df |>
    filter(loc == "slope" & is.finite(l_occ_h)) |>
    ggplot(aes(x = l_occ_h, y = rate, fill = rate)) +
    geom_density_ridges(stat = "binline", bins = 30, scale = 1) +
    labs(title = title_slope,
         x = y_axis_label,
         y = x_axis_label,
         fill = legend_title) +
    geom_vline(xintercept = st_sep_strat_used[2:(length(st_sep_strat_used)-1)],
               linetype = "dashed",
               color = "black") +
    geom_text(data = df_text,
              aes(x = time,
                  y = height,
                  label = label),
              inherit.aes = FALSE) +
    coord_flip()
  p3 = ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, labels = LETTERS[1:2], common.legend = TRUE, legend = "bottom") 
  return(p3)  
}
p = plot_lo_comparison()
p
ggsave("figs/last_occ_pl.png", p)
p = plot_lo_comparison(case = "ra")
p
ggsave("figs/last_occ_ra.png", p)

#### Plots: Incompleteness, gap statistics, and gap distribution ####
names = c("comp", "max", "quant_1", "quant_3", "median", "dist", "case")
df = data.frame(matrix(nrow = 0, ncol = length(names)))
names(df) = names

names = c("dist", "case", "hiat_duration")
df1 = data.frame( matrix(nrow = 0, ncol = length(names)))
names(df1) = names

for (case in names(adm_comb)){
  for (i in seq_along(distances_km)){
    adm = adm_comb[[case]][[i]]
    comp = get_completeness(adm)
    hiat = get_hiat_duration(adm)
    df3 = data.frame(dist = rep(distances[i], length(hiat)),
                     case = rep(case, length(hiat)),
                     hiat_duration = hiat)
    df2 = data.frame(
      comp = comp,
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
    df1 = rbind(df1, df3)
  }
}
df$max[is.infinite(df$max)] = 0
df$quant_1[is.na(df$quant_1)] = 0
df$quant_3[is.na(df$quant_3)] = 0
df$median[is.na(df$median)] = 0

df$case = factor(df$case)
df1$case = factor(df1$case)
df1$dist = factor(df1$dist)

plot_completeness = function(){
  p = df |>
  ggplot(aes(x = dist, y = comp, color = case)) + 
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

plot_gap_statistics = function(){
  p1 = df |> 
    select(-comp) |>
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
    theme(legend.position = c(0.9, 0.8))
  p2 = df |> 
    select(-comp) |>
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
    theme(legend.position = c(0.9, 0.8))
  p3 = ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, labels = LETTERS[1:2],common.legend = TRUE, legend = "bottom")
  return(p3)
}
p = plot_gap_statistics()
p
ggsave(filename = "figs/gap_statistics.png",
       plot = p)

plot_hiat_duration = function(){
  pos_interest = c(3,6,9, 12, 15)
  
  p1 = df1 |> 
    filter(case == "pl" & dist %in% pos_interest) |>
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
    scale_fill_discrete(labels = c("3 km", "6 km", "9 km", "12 km", "15 km")) 
  p1
  
  p2 = df1 |> filter(case == "ra" & dist %in% pos_interest) |>
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
    scale_fill_discrete(labels = c("3 km", "6 km", "9 km", "12 km", "15 km"))
  p2
  p3 = ggpubr::ggarrange(p1, p2, nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom", labels = LETTERS[1:2])
  return(p3)
}
p = plot_hiat_duration()
p
ggsave("figs/hiatus_duration_comp.png", p)

plot_no_of_hiat = function(){
  p = df |>
    ggplot(aes(x = dist, y = n_hiat, color = case)) +
    geom_line(linewidth = 3) +
    labs(title = "Number of Hiatuses",
         y = "Number of hiatuses",
         x = "Distance from shore [km]",
         color = "Geometry") +
    scale_color_discrete(labels = c("Platform", "Ramp")) +
    theme(legend.position = c(0.1,0.9))
  return(p)
}
p = plot_no_of_hiat()
p
ggsave("figs/no_of_hiatuses.png", p)

plot_accumulated_sediment = function(){
  p = df |>
    ggplot(aes(x = dist, y = acc_sed, color = case)) +
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

#### Extinction patterns ####
ext_scen = list()
ext_scen[["HST"]] = approxfun(x = c(0.25,0.5, 0.75,2.25, 2.5, 2.75), y = c(1, 25,1, 1, 25, 1), rule = 2)
ext_scen[["RST"]] = approxfun(x = c(0.75, 1, 1.25, 2.75, 3, 3.25), y = c(1, 25,1, 1, 25, 1), rule = 2)
ext_scen[["LST"]] = approxfun(x = c(1.25, 1.5, 1.75, 3.25, 3.5, 3.75), y = c(1, 25,1, 1,25,1), rule = 2)
ext_scen[["TST"]] = approxfun(x = c(0, 0.25,1.75, 2, 2.25, 3.75,4), y = c(25,1,1, 25,1,1,25), rule = 2)
ext_scen[["constant"]] = approxfun(x = c(0, 4), y = c(1,1), rule = 2)
# for (scen in names(ext_scen)){
#   f = ext_scen[[scen]]
#   t = seq(0,4, by = 0.001)
#   plot(t, f(t), main = scen, ylim = c(0, 25))
#   lines(sl$t, 10 *sl$sl/max(sl$sl) + 10)
# }

plot_ext_scenario_comparison = function(rate, dist, case, title){
  stopifnot(dist %in% distances)
  stopifnot(case %in% cases)
  names = c("ext_scen", "case", "dist", "l_occ_h")
  df = data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  rate = rate
  n_locc = 1000
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
  df$ext_scen = factor(df$ext_scen)
  
  st_sep_time = c(0,seq(0.25, 3.75, by = 0.5), 4)
  st_sep_strat = time_to_strat(st_sep_time, adm, destructive = FALSE)
  st_pres = !(diff(st_sep_strat) == 0)
  st_names = c("TST", "HST", "RST","LST","TST","HST", "RST", "LST", "TST")
  st_names_used = st_names[st_pres]
  st_sep_strat_used = st_sep_strat[st_pres]
  df_text = data.frame(time = 0.5* (head(st_sep_strat_used, -1) + tail(st_sep_strat_used, -1)),
                       height = rep(2.5, length(st_sep_strat_used)-1),
                       label = st_names_used)
  
  p = df |>
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
ggsave(filename = "figs/lo_comp_different_ext_scenarios.png",
       plot = p3)

plot_ext_scen_by_rate = function(rate){
  p1 = plot_ext_scenario_comparison(rate = rate, dist = 15, case = "pl", title = "Platform proximal slope")
  p2 = plot_ext_scenario_comparison(rate = rate, dist = 3, case = "pl", title = "Platform top")
  p3 = ggpubr::ggarrange(p2,p1, ncol = 2, nrow =1, common.legend = TRUE, legend = "bottom", labels = LETTERS[1:2])
  name = paste0("figs/lo_comp_different_ext_scenarios_rate_", rate, ".png")
  ggsave(filename = name,
         plot = p3)
}
for (rate in c(1,3,10,30,100,300,1000,3000)){
  plot_ext_scen_by_rate(rate = rate)
}

#### Spatioal correlation of extinction scenarios ####

plot_spat_corr_ext = function(rate, case, ext_sce){
sel_loc = seq(3, 15, by = 3)
ext_s = ext_scen[[ext_sce]]
case = case
stopifnot(dist %in% distances)
stopifnot(case %in% cases)
names = c("case", "dist", "l_occ_h")
df = data.frame(matrix(nrow = 0, ncol = length(names)))
names(df) = names
rate = rate
n_locc = 1000
for (dist in sel_loc){
  f = ext_s
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
  df2 = data.frame(case = rep(case, length(l_occ_h)),
                   dist = rep(dist, length(l_occ_h)),
                   l_occ_h = l_occ_h)
  df = rbind(df, df2)
}
df$case = factor(df$case)
df$dist = factor(df$dist, levels = sel_loc, ordered = TRUE)

# correlation lines
st_sep_time = c(0,seq(0.25, 3.75, by = 0.5),4)
names = c("h", "dist", "trans")
df3 = data.frame(matrix(nrow = 0, ncol = length(names)))
names(df3) = names
for (dist in sel_loc){
  adm = adm_comb[[case]][[which(distances == dist)]]
  df4 = data.frame(h = time_to_strat(st_sep_time, adm, destructive = FALSE),
                   dist = rep(dist, length(st_sep_time)),
                   trans = letters[1:(length(st_sep_time))])
  df3 = rbind(df3, df4)
}
df3  = df3|> arrange(dist, h)
df3$dist = factor(df3$dist, levels = sel_loc, ordered = TRUE)
df3$trans = factor(df3$trans, ordered = TRUE)

names = c("pos", "labs", "name")
df_text = data.frame(matrix(nrow = 0, ncol = length(names)))
names(df_text) = names
st_sep = c(0,seq(0.25, 3.75, by = 0.5), 4)
labs = c("TST", "HST", "RST","LST","TST","HST", "RST", "LST", "TST")
for (plot_pos in sel_loc){
  pos1 = time_to_strat(st_sep, adm_comb[[case]][[which(distances == plot_pos)]], destructive = FALSE)
  is_pres = diff(pos1) != 0
  height =  0.5* (head(pos1, -1) + tail(pos1, -1))
  df_text1 = data.frame(pos = height[is_pres],
                       labs = labs[is_pres],
                       dist = factor(rep(plot_pos, sum(is_pres))))
  df_text = rbind(df_text, df_text1)
}



p = df   |>
  filter(!is.na(l_occ_h))|>
  ggplot(aes(x = l_occ_h, y = dist, group = dist, fill = dist)) +
  geom_density_ridges(stat = "binline",
                      bins = 40,
                      scale = 1)+
  geom_path(data = df3 ,
                aes(x = h, y = dist, group = trans)) +
  geom_text(data = df_text |> filter(dist %in% c(3, 15)),
            aes(x = pos, y = dist, label = labs),
            position = position_nudge(y = 0.2)) +
  coord_flip() +
  labs(y = "Distance",
       x = "Stratigraphic height [m]",
       fill = "Distance")
return(p)

}
for (case in c("ra", "pl")){
  for (ext_scenario in c("HST", "TST", "RST", "LST", "constant")){
    print(case)

    p = plot_spat_corr_ext(rate = 10, case = case, ext_sce = ext_scenario)
    p
    ggsave(filename = paste0("figs/ext_comp",case, ext_scenario,".png"))
  }
}


p
st_sep_time = c(0,seq(0.25, 3.75, by = 0.5), 4)
st_sep_strat = time_to_strat(st_sep_time, adm, destructive = FALSE)
st_pres = !(diff(st_sep_strat) == 0)
st_names = c("TST", "HST", "RST","LST","TST","HST", "RST", "LST", "TST")
st_names_used = st_names[st_pres]
st_sep_strat_used = st_sep_strat[st_pres]
df_text = data.frame(time = 0.5* (head(st_sep_strat_used, -1) + tail(st_sep_strat_used, -1)),
                     height = rep(2.5, length(st_sep_strat_used)-1),
                     label = st_names_used)

p = df |>
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



#### Plot range offset ####
min_depth = 0
max_depth = max(unlist(wd_comb))
min_width = 3
max_width = 50
n_subdiv  = 15
min_rate = 1
max_rate = 100
n_rates = 30
depth_optima = max_depth * seq(0,1, length.out = n_subdiv)^2
niche_width = min_width + (max_width - min_width)* seq(0,1,length.out = n_subdiv)^2
n_list = list()
k = 1
for (i in 1:n_subdiv){
  for (j in 1:n_subdiv){
    n_list[[k]] = list(opt = depth_optima[i],
                       width = niche_width[j])
    k = k + 1
  }
}
rate_range = min_rate + (max_rate - min_rate) * seq(0,1, length.out = n_rates)^2
rate_range
names = c("case", "distance", "opt", "width", "t_ext", "r_offset_h", "r_offset_t","h_ext", "rate")
df = data.frame(matrix(nrow = 0, ncol = length(names)))
names(df) = names

case = "ra"
dist = 3
adm = adm_comb[[case]][[which(distances == dist)]]
get_range_offset = function(){
  for (i in 1:1000){
    x = c("t" = NA, "h" = NA)
    while (is.na(x["t"])){
      t_ext = p3(1, 0,4, 1)
      rate = sample(rate_range, 1)
      ni = sample(n_list, 1)
      x = range_offset(t_ext = t_ext, 
                       rate = rate,
                       adm = adm,
                       niche = snd_niche(opt = ni[[1]]$opt,
                                         tol = ni[[1]]$width,
                                         cutoff_val = 0))
    }
    df2 = data.frame(case = case, 
                     distance = dist,
                     opt = ni[[1]]$opt,
                     width = ni[[1]]$width, 
                     t_ext = t_ext,
                     r_offset_h = x["h"],
                     r_offset_t = x["t"],
                     h_ext = time_to_strat(t_ext, adm),
                     rate = rate
                     )
    df = rbind(df, df2)
  }
  return(df)
}
hh = get_range_offset()
plot(hh$h_ext, hh$r_offset_t)
plot(hh$h_ext, hh$r_offset_h)
plot(hh$t_ext, hh$r_offset_t)
plot(hh$t_ext, hh$r_offset_h)
plot(hh$width, hh$r_offset_t)

opt = 1
width = 3
case = "ra"
dist = 3
r_foss = 0.1
adm = adm_comb[[case]][[which(distances == dist)]]
wd_loc = wd_comb[[case]][[which(distances == dist)]]
rate = 5
t_ext = p3(rate = 1, from  = 0, to = 4, n = 1)
niche = snd_niche(opt = opt, tol = width, cutoff_val = 0)
gc = approxfun(x = wd_loc$t,y = wd_loc$wd )
f_occ = p3(rate = r_foss, from = 0, to = t_ext) |>
  apply_niche(niche_def = niche,
              gc = gc) |>
  time_to_strat(adm)
f_occ = f_occ[!is.na(f_occ)]
highest_occ = max(f_occ)
h_ext = time_to_strat(t_ext, adm, destructive = FALSE)
t_last_occ = strat_to_time(highest_occ, adm)
offset_range_t = t_ext - t_last_occ
offset_range_h = h_ext - highest_occ

#### Some evol blobs ####
t = wd_comb$ra[[1]]$t
st = StratPal::random_walk(t = t, mu = 3)
plot(st)

names = c("h", "val", "offset", "group")
df = data.frame(matrix(nrow = 0, ncol = length(names)))
names(df) = names

sel_loc = seq(3, 21, by = 3)
case = "pl"
for (i in seq_along(sel_loc)){
  dist = sel_loc[i]
  adm = adm_comb[[case]][[which(dist == sel_loc)]]
  plot(adm)
  t_trans = time_to_strat(st, adm)
  t_trans$y = t_trans$y/max(t_trans$y, na.rm = TRUE) * 0.95
  df2 = data.frame(h = t_trans$h,
                   val = t_trans$y + i-1,
                   group = as.character(dist))  
  df = rbind(df, df2)
}
df$group = factor(df$group, ordered = TRUE)

st_sep_time = c(0,seq(0.25, 3.75, by = 0.5),4)
names = c("h", "dist", "trans")
df3 = data.frame(matrix(nrow = 0, ncol = length(names)))
names(df3) = names
for (i in seq_along(sel_loc)){
  dist = sel_loc[i]
  adm = adm_comb[[case]][[which(dist == sel_loc)]]
  df4 = data.frame(h = time_to_strat(st_sep_time, adm, destructive = FALSE),
                   dist = rep(i-1, length(st_sep_time)),
                   trans = letters[1:(length(st_sep_time))])
  df3 = rbind(df3, df4)
}
df3  = df3|> arrange(dist, h)
df3$trans = factor(df3$trans, ordered = TRUE)

df |> ggplot(aes(x = h, y = val, group = group, color = group)) +
  geom_line() +
  geom_hline(yintercept = 0:7) +
  geom_path(data = df3, aes(x = h, y = dist, group = trans), inherit.aes = FALSE) +
  coord_flip()
