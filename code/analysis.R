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

col_names = c("wd", "t", "dist", "tag")
df_wd = data.frame(matrix(ncol = 4, nrow = 0))
names(df_wd) = col_names
for (case in cases){
  for (i in seq_len(n_locations)){
    wd = wd_comb[[case]][[i]]
    df_temp = data.frame(wd = wd$wd,
                         t = wd$t,
                         dist = rep(distances[i], length(wd$t)),
                         tag = rep(case, length(wd$t)))
    df_wd = rbind(df_wd, df_temp)
  }
}

df_wd$dist = factor(df_wd$dist, levels = distances)

df_wd |>
  filter(tag == "ra") |>
ggplot( aes(x = t, y = wd, color = dist)) + 
  geom_line()

ggsave(paste0(path_base_fig, "sl_ramps_comp.png"))

df_wd |>
  filter(tag == "pl") |>
  ggplot( aes(x = t, y = wd, color = dist)) + 
  geom_line()

ggsave(paste0(path_base_fig, "sl_platform_comp.png"))


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
adm = adm_list[[1]]
times_ext = seq(0, max_time(adm), by = 0.01)
range_offset_t = rep(NA, length(times_ext))
range_offset_h = rep(NA, length(times_ext))
for (i in seq_along(times_ext)){
  x = range_offset(t_ext = times_ext[i], rate =  10, adm)
  range_offset_h[i] = x["h"]
  range_offset_t[i] = x["t"]
}

plot(times_ext, range_offset_t, type = "l")
plot(times_ext, range_offset_h, type = "l")

## last occurrences ##
adm = adm_list[[1]]
rates = c(2, 5, 10, 30, 100)
df = data.frame(l_occ_h = NULL, rate = NULL)
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
  df = rbind(df, data.frame(l_occ_h = l_occ_h, rate = rep(rate, length(l_occ_h))))
}
df$rate = factor(df$rate)

df |> ggplot(aes(x = l_occ_h, colour = rate, fill = rate)) +
  geom_histogram(position = "identity", alpha = 0.5) +
  coord_flip() + 
  labs(title = "Last occurrences as a function of fossil abundance",
       x = "Abundance",
       y = "Stratigraphic position")


#### condensation ratio plot ####
df = data.frame(matrix(nrow = 0, ncol = 4))
names = c("cond_ratio", "proportion_height", "distance", "system")
names(df) = names

bin_width_m = 2

for (i in seq_along(distances)){
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
                   distance = rep(distances[i], length(cond_ratio)),
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
                   distance = rep(distances[i], length(cond_ratio)),
                   system = rep("platform", length(cond_ratio)))
  df = rbind(df, df2)
}

df$distance = factor(df$distance, levels = distances)
df$system = factor(df$system)

pos_interest = c(3,6,9, 12, 15, 18, 21)
y_lim = range(df$cond_ratio)
y_lim[1] = 0.9 * y_lim[1]
y_lim[2] = 1.1 * y_lim[2]
p1 = df |> 
  filter(system == "ramp" &
           distance %in% pos_interest) |>
  ggplot(aes(x = proportion_height, y = cond_ratio, color = distance)) +
  geom_line() +
  ylim(y_lim) +
  coord_flip() +
  scale_y_log10() +
  labs(y = "Inflation of rates",
       x = "Relative height",
       title = "Inflation of rates on Ramp")
p1

p2 = df |> 
  filter(system == "platform" &
           distance %in% pos_interest) |>
  ggplot(aes(x = proportion_height, y = cond_ratio, color = distance)) +
  geom_line() +
  ylim(y_lim) +
  coord_flip() +
  scale_y_log10() +
  labs(y = "Inflation of rates",
       x = "Relative height",
       title = "Inflation of rates on Platform")
p2

p3 = ggpubr::ggarrange(p1, p2, nrow = 1, ncol = 2)
p3
ggsave("figs/inflation_of_rates.png", plot = p3)

#### AGe-depth model plots ####
names = c("time", "height", "distance", "system")
df = data.frame(matrix(nrow = 0, ncol = length(names)))
names(df) = names

for (i in seq_along(distances)){
  adm = adm_list_ra[[i]]
  t = get_T_tp(adm)
  h = time_to_strat(t, adm, destructive = TRUE)
  df2 = data.frame(time = t,
                   height = h,
                   distance = rep(distances[i], length(time)),
                   system = rep("ramp", length(time)))
  df = rbind(df, df2)
  adm = adm_list_pl[[i]]
  t = get_T_tp(adm)
  h = time_to_strat(t, adm, destructive = TRUE)
  df2 = data.frame(time = t,
                   height = h,
                   distance = rep(distances[i], length(time)),
                   system = rep("platform", length(time)))
  df = rbind(df, df2)
}

df$distance = factor(df$distance, levels = distances)
df$system = factor(df$system)

pos_interest = c(3,6,9, 12, 15, 18, 21)

p1 = df |>
  filter(system == "ramp" &
           distance %in% pos_interest) |>
  ggplot(aes(x = time, y = height, color = distance)) +
  geom_line() +
  labs(title = "ramp")
p1
p2 = df |>
  filter(system == "platform" &
           distance %in% pos_interest) |>
  ggplot(aes(x = time, y = height, color = distance)) +
  geom_line() +
  labs(title = "platform")
p2
p3 = ggpubr::ggarrange(p1, p2, nrow = 1, ncol = 2, common.legend = TRUE)
p3

ggsave("figs/age_depth_overview.png", p3)
