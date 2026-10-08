case_labels = c(pl = "Platform", ra = "Ramp")

library(admtools)
library(StratPal)
library(configr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tomledit)
library(ggpubr)
library(ggridges)

#### Load data ####
source("code/load_data_and_constants.R")

niche1 = snd_niche(opt = 3, tol = 10, cutoff_val = 0)
niche2 = snd_niche(opt = 20, tol = 20, cutoff_val = 0)
niche3 = snd_niche(opt = 50, tol = 40, cutoff_val = 0)

niche_list = list("name1" = niche1,
                  "name2" = niche2,
                  "name3" = niche3)

# Row 1: sampling probability vs water depth
sampling_prob_vs_wd = function(niche_list, wd_min = -3, wd_res = 0.1){
  wd_max = max(c(sapply(wd_pl, function(x) max(x$wd)),
                 sapply(wd_ra, function(x) max(x$wd))))
  wds = seq(from = wd_min, to = wd_max, by = wd_res)
  
  df = bind_rows(lapply(names(niche_list), function(n)
    data.frame(wd = wds,
               sampling_prob = niche_list[[n]](wds),
               name = n)))
  
  ggplot(df, aes(x = wd, y = sampling_prob, colour = name)) +
    geom_line() +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Water depth [m]", y = "Sampling probability", colour = "Niche")
}

# Rows 2 and 3: sampling probability vs stratigraphic height
sampling_prob_vs_h = function(niche_list, case, pos, h_res = 1.2){
  idx    = which(distances == pos)
  adm    = adm_comb[[case]][[idx]]
  wd_loc = wd_comb[[case]][[idx]]
  wd_fun = approxfun(x = wd_loc$t, y = wd_loc$wd, rule = 2)
  
  h_sampled = seq(from = 0.5 * h_res, to = admtools::max_height(adm), by = h_res)
  wd_at_h   = h_sampled |> strat_to_time(adm) |> wd_fun()
  
  df = bind_rows(lapply(names(niche_list), function(n)
    data.frame(h = h_sampled,
               sampling_prob = niche_list[[n]](wd_at_h),
               name = n)))
  
  ggplot(df, aes(x = sampling_prob, y = h, colour = name)) +
    geom_path() +                          # path keeps the order along height
    scale_x_continuous(limits = c(0, 1)) +
    labs(x = "Sampling prob.", y = "Height [m]", colour = "Niche",
         title = paste0(case_labels[[case]], ", ", pos))
}

make_row = function(case){
  plots = lapply(positions_examined, function(pos)
    sampling_prob_vs_h(niche_list, case, pos))
  ggarrange(plotlist = plots, ncol = length(plots), nrow = 1, legend = "none")
}

p_wd = sampling_prob_vs_wd(niche_list)
leg  = get_legend(p_wd + theme(legend.position = "bottom"))

row1 = ggarrange(NULL, p_wd, NULL, nrow = 1, widths = c(1, 3, 1), legend = "none")
row2 = make_row("pl")
row3 = make_row("ra")

p_all = ggarrange(row1, row2, row3,
                  ncol = 1, nrow = 3,
                  heights = c(1, 1.2, 1.2),
                  labels = c("A", "B", "C"),
                  common.legend = TRUE, legend.grob = leg, legend = "bottom")
p_all
