library(admtools)
library(StratPal)
library(configr)
library(ggplot2)

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
  png(filename = paste0(path_base_fig, "wd_", loc, "km.png"))
  plot(wd[[i]]$t, wd[[i]]$wd, type = "l")
  dev.off()
  png(filename = paste0(path_base_fig, "adm_", loc, "km.png"))
  plot(adm_list[[i]], lty_destr = 0)
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

