library(admtools)
library(StratPal)
library(configr)
library(ggplot2)
library(dplyr)
library(tidyr)

tag = "carbonate_stratpal_1"

adm_data = read.csv(paste0("data/", tag, "_adm.csv"))
wd_data = read.csv(paste0("data/", tag, "_wd.csv"))
t = adm_data$time..Myr.
t_steps = adm_data$timestep...

#### Metadata ####
meta_data = tomledit::read_toml(paste0("data/", tag, ".toml"))
l = tomledit::from_toml(meta_data)
n_locations = length(l$locations)
loc_list = list()
distances = c()
for (loc in seq_len(n_locations)){ # convert to km
  loc_list[[loc]] = list(x = l$locations[[loc]]$x[[1]]/ 1000,
                         y = l$locations[[loc]]$y[[1]] / 1000)
  distances =  c(distances, l$locations[[loc]]$x[[1]]/ 1000)
}

#### Age-depth models ####
adm_list = list()
for (i in 1:(length(adm_data)-2)){
  adm_list[[i]] = tp_to_adm(t = adm_data$time..Myr.,
                            h = adm_data[, paste0("adm_", i, "..m.")])
}

#### Water depth ####
wd = list()
for (i in 1:(length(wd_data)-2)){
  wd[[i]] = list(t = wd_data$time..Myr.,
                 wd = wd_data[, paste0("wd_", i, "..m.")])
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

for (i in seq_len(n_locations)){
  loc = distances[i]
  png(filename = paste0(path_base_fig, "wd_time_", loc, "km.png"))
  plot(wd[[i]]$t, wd[[i]]$wd, type = "l")
  dev.off()
  png(filename = paste0(path_base_fig, "adm_", loc, "km.png"))
  plot(adm_list[[i]], lty_destr = 0)
  dev.off()
  png(filename = paste0(path_base_fig, "wd_strat_", loc, "km.png"))
  wd_loc = list(t = wd[[i]]$t, y = wd[[i]]$wd) |>
    time_to_strat(adm_list[[i]]) |>
    plot(type = "l")
  dev.off()
}

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
cases = c("const", "HST", "TST", "RST", "LST")
for (i in seq_len(n_locations)){
  adm = adm_list[[i]]
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
      ggtitle(paste0("# LO ", pos, " km from shore"))
    ggsave(filename = paste0(path_ex_fig, case, "_", pos, "km.png"))
    
  }
}


#### Sed Stats ####
path_sed_stats = "figs/sed_stats/"
if (!dir.exists(path_sed_stats)) dir.create(path_sed_stats, recursive = TRUE)
comp = c()
for(i in seq_along(adm_list)){
  comp = c(comp, adm_list[[i]] |> get_completeness())
}
df = data.frame(dist = distances,
                comp = comp)

df |> ggplot(aes(x = dist, y = comp)) +
  geom_line() +
  ylim(c(0,1))
ggsave(filename = paste0(path_sed_stats, "completeness.png"))

## other parameters to track:
# 1. maximum/minimum/median/mean hiatus duration
# 2. thickness
# 3. sed rates?
# -> same parameters per systems tract

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
