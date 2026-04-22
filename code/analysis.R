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

#### Figure 2: Age-depth models ####
plot_age_depth_models = function(pos, plot_st = TRUE){
  #' pos: positions of adms
  #' plot_st: plot systems tracts?
  stopifnot(any( pos %in% distances))
  st_sep_time_mod = c(-0.2,st_sep_time, 4.2)
  names = c("time", "height", "distance", "system","lty")
  df = data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
  title_platform = "Platform"
  title_ramp = "Ramp"
  x_lab_title = "Elapsed Model Time [Myr]"
  y_lab_title = label_height_absolute
  legend_label = "Distance from Shore"
  distances_km = paste(pos, "km")
  
  for (case in cases){
    for (i in seq_along(pos)){
      adm = adm_comb[[case]][[which(distances == pos[i])]]
      sac = sac_comb[[case]][[which(distances == pos[i])]]
      t = get_T_tp(adm)
      h = time_to_strat(t, adm, destructive = TRUE)
      destr = is.na(h)
      h_sac = sac$h
      h_sac[!destr] = NA
      df2 = data.frame(time = rep(t, 2),
                       height = c(h, h_sac),
                       distance = rep(distances_km[i], 2 * length(t)),
                       system = rep(case, 2 * length(t)),
                       lty = c(rep("A", length(t)),rep("B", length(t))))
      df = rbind(df, df2)
    }
  }
  df$distance = factor(df$distance, levels = distances_km)
  df$system = factor(df$system)
  
  p1 = df |>
    filter(system == "ra") |>
    ggplot(aes(x = time,
               y = height, 
               color = distance, 
               linetype = lty)) +
    geom_line() +
    labs(title = title_ramp,
         x = x_lab_title,
         y = y_lab_title,
         col = legend_label) +
    guides(linetype = "none") +
    scale_color_hue()
  
  if (plot_st){
    max_height = df |> filter(system == "ra") 
    max_height = max_height$height[!is.na(max_height$height)] |> max()
    df_text = data.frame(time = 0.5* (head(st_sep_time_mod, -1) + tail(st_sep_time_mod, -1)),
                         height = rep(max_height, length(st_sep_time_mod)-1),
                         label = st_labels_en)
    
    p1 = p1 +  geom_vline(xintercept = head(tail(st_sep_time_mod, -1), -1),
                          linetype = "dashed",
                          color = "black",
                          lwd = 0.5) +
      geom_text(data = df_text,
                aes(x = time,y = height,label = label),
                inherit.aes = FALSE,
                size = st_text_size)
    
  }
  
  p2 = df |>
    filter(system == "pl") |>
    ggplot(aes(x = time, 
               y = height, 
               color = distance,
               linetype = lty)) +
    geom_line() +
    labs(title = title_platform,
         x = x_lab_title,
         y = y_lab_title,
         col = legend_label) +
    guides(linetype = "none") +
    scale_color_hue()
  if (plot_st){
    max_height = df |> filter(system == "pl") 
    max_height = max_height$height[!is.na(max_height$height)] |> max()
    df_text = data.frame(time = 0.5* (head(st_sep_time_mod, -1) + tail(st_sep_time_mod, -1)),
                         height = rep(max_height, length(st_sep_time_mod)-1),
                         label = st_labels_en)
    
    p2 = p2 +  geom_vline(xintercept = head(tail(st_sep_time_mod, -1), -1),
                          linetype = "dashed",
                          color = "black",
                          lwd = 0.5) +
      geom_text(data = df_text,
                aes(x = time,y = height,label = label),
                inherit.aes = FALSE,
                size = st_text_size)
    
  }
  p3 = ggpubr::ggarrange(p2, p1, nrow = 1, ncol = 2,
                         common.legend = TRUE,
                         legend = "bottom",
                         labels = LETTERS[1:2])
  return(p3)
}

plot_fig2 = function(){
  p = plot_age_depth_models(pos = positions_examined,
                            plot_st = TRUE)
  ggsave("figs/ms/fig2.png",
         plot = p,
         width = width_2_col_cm,
         height = 10,
         units = "cm",
         bg = "white")
}
plot_fig2()

#### Figure 3: Water depth in the stratigraphic domain ####
plot_wd_strat_domain = function(pos = seq(3, 21, by = 3)){
  # pos: distance from shore for plot
  stopifnot(all(pos %in% distances))
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
  res_distances = pos |> paste("km")
  title_ramp = "Ramp"
  title_platform = "Platform"
  wd_label = label_wd
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
         col = legend_title) +
    scale_color_hue()
  p2 = df |>
    filter(system == "ra" & distance %in% res_distances) |>
    ggplot(aes(y = wd, x = proportion_height, col = distance)) +
    geom_line() +
    coord_flip() +
    labs(title = title_ramp,
         y = wd_label,
         x = height_label,
         col = legend_title) +
    scale_color_hue()
  p3 = ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, labels = LETTERS[1:2], legend = "bottom", common.legend = TRUE)
  return(p3)
}
plot_fig3 = function(){
  p = plot_wd_strat_domain(pos = positions_examined)
  ggsave("figs/ms/fig3.png",
         plot = p,
         width = width_2_col_cm,
         height = 10,
         units = "cm",
         bg = "white")
}
plot_fig3()

#### Figure 4: condensation ratio ####
plot_condensation_ratio = function(pos){
  # pos: distance from shore
  stopifnot(all(pos %in% distances))
  # bin width over which condensation ratio is calculated
  bin_width_m = 4
  names = c("cond_ratio", "proportion_height", "distance", "system")
  df = data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  distances_km  = paste(pos, "km")
  for (case in names(adm_comb)){
    for (i in seq_along(distances_km)){
      adm = adm_comb[[case]][[which(pos[i] == distances)]]
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
  
  y_lim = range(df$cond_ratio)
  y_lim[1] = 0.95 * y_lim[1]
  y_lim[2] = 1.1 * y_lim[2]
  title_ramp = "Ramp"
  title_platform = "Platform"
  y_label = "Condensation Ratio [-]"
  x_label = "Relative Height [-]"
  legend_label = "Distance from Shore"
  p1 = df |> 
    filter(system == "pl") |>
    ggplot(aes(x = proportion_height, y = cond_ratio, color = distance)) +
    geom_step(direction = "mid") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    coord_flip() +
    scale_y_log10(limits = y_lim) +
    labs(y = y_label,
         x = x_label,
         title = title_platform,
         color = legend_label) +
    scale_color_hue()
  
  p2 = df |> 
    filter(system == "ra") |>
    ggplot(aes(x = proportion_height, y = cond_ratio, color = distance)) +
    geom_step(direction = "mid") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    coord_flip() +
    scale_y_log10(limits = y_lim) +
    labs(y = y_label,
         x = x_label,
         title = title_ramp,
         color = legend_label) +
    scale_color_hue()
  
  p3 = ggpubr::ggarrange(p1, p2,  nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom", labels = LETTERS[1:2])
  return(p3)
}

plot_fig4 = function(){
  p = plot_condensation_ratio(pos = positions_examined)
  ggsave("figs/ms/fig4.png",
         plot = p,
         width = width_2_col_cm,
         height = 10,
         units = "cm",
         bg = "white")
}
plot_fig4()

#### Figure 6: Abundance bias on last occurrences ####
plot_lo_by_rate = function(case, pos, rates, title, plot_st = TRUE){
  # case: "ra" or "pl", which simulation?
  # pos: distance from shore (in km)
  # rates: numeric vector, rate of foss. occ per Myr for each taxon
  # title: title for plot
  # plot_st: plot systems tracts?
  stopifnot(case %in% cases)
  stopifnot(pos %in% distances)
  # no of extinctions observed
  n_locc = n_ext
  
  y_axis_label = "Stratigraphic Height [m]"
  legend_title = "Sampling Rate [Fossils/Myr]"
  x_axis_label = "Sampling Rate [Fossils/Myr]"
  title = title
  
  names = c("l_occ_h","rate")
  df = data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
  adm = adm_comb[[case]][[which(pos == distances)]]
  
  for (j in seq_along(rates)){
    rate = rates[j]
    t_ext = p3(1, from = min_time(adm), to = max_time(adm), n = n_ext)
    l_occ_h = rep(NA, length(t_ext))
    if (is.finite(rate)){
      for (i in seq_along(t_ext)){
        x = last_occ(t_ext = t_ext[i], rate, adm = adm)
        l_occ_h[i] = x["h"]
      }
    }
    if (is.infinite(rate)){ # special case: infinate rate -> directly transform last occ 
      l_occ_h = time_to_strat(t_ext, adm, destructive = FALSE)
    }
    df = rbind(df,
               data.frame(l_occ_h = l_occ_h,
                          rate = rep(rate, length(l_occ_h))))
  }
  df$rate = factor(df$rate, levels = rates)
  
  p1 = df |>
    filter( is.finite(l_occ_h)) |>
    ggplot(aes(x = l_occ_h, y = rate, fill = rate)) +
    geom_density_ridges(stat = "binline", binwidth = 2, scale = 1) 
  
  if (plot_st){
    st_sep_time_mod = c(0,st_sep_time, 4)
    st_sep_strat = time_to_strat(st_sep_time_mod, adm, destructive = FALSE)
    st_pres = !(diff(st_sep_strat) == 0)
    st_names_used = st_labels_en[st_pres]
    st_sep_strat_used = st_sep_strat[st_pres]
    df_text = data.frame(time = 0.5* (head(st_sep_strat_used, -1) + tail(st_sep_strat_used, -1)),
                         height = rep(3.5, length(st_sep_strat_used)-1),
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
    coord_flip() + 
    theme(legend.position = "none")
  return(p1)
}

plot_fig6 = function(seed){
  set.seed(seed)
  p = plot_lo_by_rate(case = "pl",
                      pos = positions_examined[3], 
                      rates = rates_used,
                      title = '')
  ggsave(filename = "figs/ms/fig6.png",
         plot = p,
         width = width_2_col_cm,
         height = 10,
         units = "cm",
         bg = "white")
}
plot_fig6(seed = seed_main + 5)

#### Figure 7: Extinctions by systems tract ####
plot_ext_scenario_comparison = function(rate, dist, case, title){
  # rate: rate of fossil occurrences per Myr for each taxon
  # dist: distance from shore
  # case: "ra" or "pl", which simulation?
  # title: plot title
  stopifnot(dist %in% distances)
  stopifnot(case %in% cases)
  
  n_locc = n_ext
  
  names = c("ext_scen", "case", "dist", "l_occ_h")
  df = data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
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
  df$ext_scen = factor(df$ext_scen,
                       levels = c("constant", "HST", "RST", "LST", "TST"),
                       ordered = TRUE)
  
  st_sep_time_mod = c(0,st_sep_time, 4)
  st_sep_strat = time_to_strat(st_sep_time_mod, adm, destructive = FALSE)
  st_pres = !(diff(st_sep_strat) == 0)
  st_names = st_labels_en
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
         y = "Extinction Scenario",
         x = "Stratigraphic Height [m]",
         fill = "Extinction Scenario") +
    theme(legend.position = "none") +
    scale_y_discrete(labels = c("constant" = "Constant",
                                "HST" = "HST",
                                "RST" = "RST",
                                "LST" = "LST",
                                "TST" = "LST"))
  return(p)
}

plot_fig7 = function(seed){
  set.seed(seed)
  rate = systract_com_rate # sampling rate
  p1 = plot_ext_scenario_comparison(rate = rate,
                                    dist = positions_examined[1],
                                    case = "pl",
                                    title = "Platform Top")
  p2 = plot_ext_scenario_comparison(rate = rate,
                                    dist = positions_examined[3],
                                    case = "pl",
                                    title = "Platform Prograding Edge")
  p3 = ggpubr::ggarrange(p1,
                         p2,
                         ncol = 2, 
                         nrow =1,
                         labels = LETTERS[1:2])
  ggsave(filename = "figs/ms/fig7.png",
         plot = p3,
         width = width_2_col_cm,
         height = 10,
         units = "cm",
         bg = "white")
}
plot_fig7(seed = seed_main + 6)

#### Figure 8 and Figure 9: Spatial correlation of extinction scenarios ####
plot_spat_corr_ext = function(rate, 
                              case, 
                              ext_sce,
                              pos,
                              plot_st = TRUE,
                              plot_txt = TRUE){
  # rate: rate of foss. occ per Myr for each taxon
  # case: simulation case
  # ext_sce: extinction scenario
  # pos: distances from shore
  # plot_st: plot systems tracts?
  # plot_txt: plot systems tracts annotations?
  stopifnot(all(pos %in% distances))
  stopifnot(case %in% cases)
  stopifnot(ext_sce %in% names(ext_scen))
  ext_s = ext_scen[[ext_sce]]
  
  n_locc = n_ext
  
  names = c("dist", "l_occ_h")
  df = data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
  for (dist in pos){
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
    df2 = data.frame(dist = rep(dist, length(l_occ_h)),
                     l_occ_h = l_occ_h)
    df = rbind(df, df2)
  }
  df$dist = factor(df$dist, levels = pos, ordered = TRUE)
  
  p = df   |>
    filter(!is.na(l_occ_h))|>
    ggplot(aes(x = l_occ_h, y = dist, group = dist, fill = dist)) +
    geom_density_ridges(stat = "binline",
                        bins = 40,
                        scale = 1) +
    scale_color_hue() +
    theme(legend.position = "none")
  if (plot_st){
    # correlation lines
    st_sep_time_mod = c(0,st_sep_time,4)
    names = c("h", "dist", "trans")
    df3 = data.frame(matrix(nrow = 0, ncol = length(names)))
    names(df3) = names
    for (dist in pos){
      adm = adm_comb[[case]][[which(distances == dist)]]
      df4 = data.frame(h = time_to_strat(st_sep_time_mod, adm, destructive = FALSE),
                       dist = rep(dist, length(st_sep_time_mod)),
                       trans = letters[1:(length(st_sep_time_mod))])
      df3 = rbind(df3, df4)
    }
    df3$dist = factor(df3$dist, levels = pos, ordered = TRUE)
    p = p + geom_path(data = df3 ,
                      aes(x = h, y = dist, group = trans))
  }
  if (plot_txt){
    names = c("pos", "labs", "name")
    df_text = data.frame(matrix(nrow = 0, ncol = length(names)))
    names(df_text) = names
    st_sep = c(0,st_sep_time, 4)
    labs = st_labels_en
    for (plot_pos in pos){
      pos1 = time_to_strat(st_sep, adm_comb[[case]][[which(distances == plot_pos)]], destructive = FALSE)
      is_pres = diff(pos1) != 0
      height =  0.5* (head(pos1, -1) + tail(pos1, -1))
      df_text1 = data.frame(pos = height[is_pres],
                            labs = labs[is_pres],
                            dist = factor(rep(plot_pos, sum(is_pres))))
      df_text = rbind(df_text, df_text1)
      p = p + geom_text(data = df_text |> filter(dist %in% c(7.5, 12)),
                        aes(x = pos, y = dist, label = labs),
                        position = position_nudge(y = -  0.1),
                        size = st_text_size)
    }
    p = p + coord_flip() +
      labs(y = "Distance from Shore [km]",
           x = "Stratigraphic Height [m]",
           fill = "Distance from Shore [km]")
  }
  return(p)
}

plot_fig8 = function(seed){
  set.seed(seed)
  rate = 10
  p = plot_spat_corr_ext(rate = rate,
                         case = "pl",
                         ext_sce = "LST",
                         pos = positions_examined[2:5])
  p = p + labs(title = "Platform")
  ggsave(filename = "figs/ms/fig8.png",
         plot = p,
         width = width_2_col_cm,
         height = 10,
         units = "cm",
         bg = "white")
}
plot_fig8(seed = seed_main + 7)

plot_fig9 = function(seed){
  set.seed(seed)
  rate = 10
  p = plot_spat_corr_ext(rate = rate, 
                         case = "ra", 
                         ext_sce = "LST",
                         pos = positions_examined)
  p = p + 
    labs(title = "Ramp")
  ggsave(filename = "figs/ms/fig9.png",
         plot = p,
         width = width_2_col_cm,
         height = 10,
         units = "cm",
         bg = "white")
}
plot_fig9(seed = seed_main + 8)

#### Supp Figure 5: Sedimentation rates ####
plot_sed_rate = function(pos = positions_examined,
                         subset = 100,
                         plot_st = TRUE){
  # pos: distances from shore
  # subset: integer, how many steps of the original time increment are combined. Effectively temporal resolutions
  # plot_st: plot systems tracts?
  
  ext_factor = 1.1
  st_sep_time_mod = c(-0.2, st_sep_time, 4.2)
  
  names = c("t", "s", "pos", "case")
  df = data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
  for (case in cases){
    for (i in seq_along(pos)){
      adm = adm_comb[[case]][[which(pos[i] == distances)]]
      t = get_T_tp(adm)
      h = get_L_tp(adm)
      sel = seq(1, length(t), by = subset)
      h_sel = h[sel]
      t_sel = t[sel]
      t_midpoint = 0.5 * (tail(t_sel,-1) + head(t_sel, -1))
      sed = diff(h_sel)/diff(t_sel)
      df = rbind(df,
                 data.frame(t = t_midpoint,
                            s = sed,
                            pos = rep(paste(pos[i], "km"), length(sed)),
                            case = rep(case, length(sed))))
    }
  }
  df$case = factor(df$case)
  df$pos = factor(df$pos, levels = paste(pos, "km"))
  ylim = range(df$s)
  p1 = df |>
    filter(case == "pl") |>
    ggplot(aes(x = t, y = s, color = pos)) +
    geom_step(direction = "mid") +
    labs(x = label_time_emt,
         y = "Sedimentation Rate [m/Myr]",
         title = "Platform",
         color = "Distance from Shore") +
    ylim(ext_factor * ylim) + 
    scale_color_hue()
  
  if (plot_st){
    max_sed = max(ylim)
    df_text = data.frame(time = 0.5* (head(st_sep_time_mod, -1) + tail(st_sep_time_mod, -1)),
                         height = rep(1.05 * max_sed, length(st_sep_time_mod)-1),
                         label = st_labels_en)
    p1 = p1 +
      geom_vline(xintercept = head(tail(st_sep_time_mod, -1),-1),
                 linetype = "dashed",
                 color = "black",
                 lwd = 0.5) +
      geom_text(data = df_text,
                aes(x = time,
                    y = height,
                    label = label),
                inherit.aes = FALSE,
                size = st_text_size)
  }
  
  p2 = df |>
    filter(case == "ra") |>
    ggplot(aes(x = t, y = s, color = pos)) +
    geom_step(direction = "mid") +
    labs(x = label_time_emt,
         y = "Sedimentation rate [m/Myr]",
         title = "Ramp",
         color = "Distance") +
    ylim(ext_factor * ylim) +
    scale_color_hue()
  if (plot_st){
    max_sed = max(ylim)
    df_text = data.frame(time = 0.5* (head(st_sep_time_mod, -1) + tail(st_sep_time_mod, -1)),
                         height = rep(1.05 * max_sed, length(st_sep_time_mod)-1),
                         label = st_labels_en)
    p2 = p2 +
      geom_vline(xintercept = head(tail(st_sep_time_mod, -1),-1),
                 linetype = "dashed",
                 color = "black",
                 lwd = 0.5) +
      geom_text(data = df_text,
                aes(x = time,
                    y = height,
                    label = label),
                inherit.aes = FALSE,
                size = st_text_size)
  }
  
  p3 = ggarrange(p1,
                 p2, 
                 nrow = 1,
                 ncol = 2,
                 common.legend = TRUE,
                 legend = "bottom",
                 labels = LETTERS[1:2])
  return(p3)
  
}

plot_sfig5 = function(){
  p = plot_sed_rate()
  ggsave(filename = "figs/sm/sfig5_sedrate.png",
         plot = p,
         width = width_2_col_cm,
         height = 10,
         units = "cm",
         bg = "white")
}

plot_sfig5()

#### Supp Figure 6: Incompleteness, gap statistics, and gap distribution ####
plot_completeness = function(){
  
  names = c("dist", "case", "comp")
  df = data.frame( matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
  for (case in names(adm_comb)){
    for (i in seq_along(distances)){
      adm = adm_comb[[case]][[i]]
      comp = get_completeness(adm) * 100
      df1 = data.frame(dist = distances[i],
                       case = case,
                       comp = comp)
      df = rbind(df, df1)
    }
  }
  
  p = df |>
    ggplot(aes(x = dist, y = comp, color = case)) + 
    geom_line(linewidth = 3) +
    ylim(c(0,100)) +
    labs(x = "Distance from Shore [km]",
         y = "Completeness [%]",
         title = "Stratigraphic Completeness",
         color = "Geometry") +
    scale_color_discrete(labels = c("Platform", "Ramp")) +
    theme(legend.position = c(0.1, 0.8))
  return(p)
}

plot_sfig6 = function(){
  p = plot_completeness()
  ggsave(filename = "figs/sm/sfig6_completeness.png",
         plot = p,
         width = width_2_col_cm,
         height = 10,
         units = "cm",
         bg = "white")
}
plot_sfig6()

#### Supp Figure 7: hiatus durations ####
plot_hiat_duration = function(pos){
  # pos: distances from shore
  stopifnot(all(pos %in% distances))
  
  dur_highlighted = c(0.001, 0.1, 1) # highlighted hiatus durations
  
  names = c("dist", "case", "hiat_duration")
  df = data.frame( matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
  for (case in names(adm_comb)){
    for (i in seq_along(pos)){
      adm = adm_comb[[case]][[which(pos[i] == distances)]]
      hiat = get_hiat_duration(adm)
      df1 = data.frame(dist = rep(pos[i], length(hiat)),
                       case = rep(case, length(hiat)),
                       hiat_duration = hiat)
      df = rbind(df, df1)
    }
  }
  df$dist = factor(df$dist)
  
  p1 = df |> 
    filter(case == "pl") |>
    ggplot(aes(x = hiat_duration, y = dist, fill = dist)) +
    geom_density_ridges(stat = "density_ridges", bandwidth = 0.08, scale = 1) +
    theme_ridges()+
    geom_vline(xintercept = dur_highlighted,
               linetype = "dashed") +
    scale_x_log10() +
    ggtitle("Platform") + 
    theme_classic()+
    labs(x = "Hiatus Duration [Myr]",
         y = "Distance from Shore [km]")+
    scale_color_hue() +
    theme(legend.position = "none")
  
  p2 = df |> 
    filter(case == "ra") |>
    ggplot(aes(x = hiat_duration, y = dist, fill = dist)) +
    geom_density_ridges(stat = "density_ridges",
                        bandwidth = 0.08, 
                        scale = 1) +
    theme_ridges()+
    geom_vline(xintercept = dur_highlighted,
               linetype = "dashed") +
    scale_x_log10() +
    ggtitle("Ramp") + 
    theme_classic()+
    labs(x = "Hiatus Duration [Myr]",
         y = "Distance from Shore [km]") + 
    scale_color_hue() +
    theme(legend.position = "none")
  p3 = ggpubr::ggarrange(p1,
                         p2, 
                         nrow = 1,
                         ncol = 2,
                         labels = LETTERS[1:2])
  return(p3)
}

plot_sfig7 = function(){
  p = plot_hiat_duration(pos = positions_examined)
  ggsave("figs/sm/sfig7_hiatus_duration.png",
         plot = p,
         width = width_2_col_cm,
         height = 10,
         units = "cm",
         bg = "white")
}
plot_sfig7()

#### Supp Figure 8: Plot number of hiatuses ####
plot_no_of_hiat = function(){
  
  names = c("dist", "case", "n")
  df = data.frame( matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
  for (case in names(adm_comb)){
    for (i in seq_along(distances)){
      df1 = data.frame(dist = distances[i],
                       case = case,
                       n = get_hiat_no(adm_comb[[case]][[i]]))
      df = rbind(df, df1)
    }
  }
  
  p = df |>
    ggplot(aes(x = dist, y = n, color = case)) +
    geom_line(linewidth = 3) +
    labs(title = "Number of Hiatuses",
         y = "Number of Hiatuses",
         x = "Distance from Shore [km]",
         color = "Geometry") +
    scale_color_discrete(labels = c("Platform", "Ramp")) +
    theme(legend.position = c(0.1,0.8))
  return(p)
}

plot_sfig8 = function(){
  p = plot_no_of_hiat()
  ggsave("figs/sm/sfig8_no_of_hiatuses.png",
         plot = p,
         width = width_2_col_cm,
         height = 10,
         units = "cm",
         bg = "white")
}
plot_sfig8()

#### Supp Figure 9: Plot section thickness ####
plot_accumulated_sediment = function(){
  subsidence_total = 80 # total subsidence in m over the 4 Myr simulated
  names = c("dist", "case", "th")
  df = data.frame( matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
  for (case in names(adm_comb)){
    for (i in seq_along(distances)){
      adm = adm_comb[[case]][[i]]
      th = get_total_thickness(adm)
      df1 = data.frame(dist = distances[i],
                       case = case,
                       th = th)
      df = rbind(df, df1)
    }
  }
  
  p = df |>
    ggplot(aes(x = dist, y = th, color = case)) +
    geom_line(linewidth = 3) +
    labs(title = "Section Thickness",
         y = "Section Thickness [m]",
         x = "Distance from Shore [km]",
         color = "Geometry") +
    scale_color_discrete(labels = c("Platform", "Ramp")) +
    theme(legend.position = c(0.1, 0.85)) +
    ylim(c(0,max(df$th))) +
    geom_hline(yintercept = subsidence_total, linetype = "dashed") 
  return(p)
}

plot_sfig9 = function(){
  p = plot_accumulated_sediment()
  ggsave("figs/sm/sfig9_section_thickness.png",
         plot = p,
         width = width_2_col_cm,
         height = 10,
         units = "cm",
         bg = "white")
}
plot_sfig9()

#### Supp Figure 10: Plot gap statistics ####
plot_gap_statistics = function(){
  plot_index = seq(1, length(distances), by = 5) # which indices to plot
  
  names = c("max", "quant_1", "quant_3", "median", "dist", "case")
  df = data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
  for (case in names(adm_comb)){
    for (i in seq_along(distances[plot_index])){
      adm = adm_comb[[case]][[plot_index[i]]]
      hiat = get_hiat_duration(adm)
      df1 = data.frame(
        max = max(hiat),
        quant_1 = quantile(hiat, probs = 0.25, names = FALSE),
        quant_3 = quantile(hiat, probs = 0.75, names = FALSE),
        median = median(hiat),
        n_hiat = length(hiat),
        acc_sed = get_total_thickness(adm),
        dist = distances[plot_index[i]],
        case  = case
      )
      df = rbind(df, df1)
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
    geom_line(linewidth = 2) +
    labs(x = "Distance from Shore [km]",
         y = "Value [Myr]",
         title = "Ramp",
         color = "Statistic") +
    scale_color_discrete(labels = c("Maximum", "Median", "1st Quartile", "3rd Quartile")) +
    theme(legend.position = c(0.9, 0.8))
  p2 = df  |>
    pivot_longer(cols = c("max", "quant_1", "quant_3", "median"),
                 names_to = "measure",
                 values_to = "value") |>
    filter(case == "pl") |>
    ggplot(aes(x = dist, y = value, color = measure)) +
    geom_line(linewidth = 2)  +
    labs(x = "Distance from Shore [km]",
         y = "Value [Myr]",
         title = "Platform",
         color = "Statistic") +
    scale_color_discrete(labels = c("Maximum", "Median", "1st Quartile", "3rd Quartile")) +
    theme(legend.position = c(0.9, 0.8))
  p3 = ggpubr::ggarrange(p2,
                         p1, 
                         ncol = 2,
                         nrow = 1, 
                         labels = LETTERS[1:2],
                         common.legend = TRUE,
                         legend = "bottom")
  return(p3)
}

plot_sfig10 = function(){
  p = plot_gap_statistics()
  ggsave(filename = "figs/sm/sfig10_gap_statistics.png",
         plot = p,
         width = width_2_col_cm,
         height = 10,
         units = "cm",
         bg = "white")
}
plot_sfig10()

#### Supp Figure 11: Cumulative distribution of time vs height ####
plot_prop_time_vs_height = function(pos){
  stopifnot(all(pos %in% distances))
  
  x_label = "Proportion Cumulative Time [-]"
  y_label = "Proportion Cumulative Height [-]"
  title_platform = "Platform"
  title_ramp = "Ramp"
  legend_label = "Distance from Shore"
  
  names = c("prop_time", "prop_height", "dist", "case")
  df = data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
  for (case in names(adm_comb)){
    for (i in seq_along(pos)){
      adm = adm_comb[[case]][[which(pos[i] == distances)]]
      r = diff(get_L_tp(adm))/diff(get_T_tp(adm)) 
      r = r |> sort() |> cumsum()
      t_rel = seq_along(r)/length(r)
      df2 = data.frame(prop_time = t_rel,
                       prop_height = r/max(r),
                       dist = rep(paste(pos[i], "km"), length(t_rel)),
                       system = rep(case, length(t_rel)))
      df = rbind(df, df2)
    }
  }
  
  df$dist = factor(df$dist, ordered = TRUE, levels = paste(pos, "km"))
  df$system = factor(df$system)
  
  p1 = df |>
    filter(system == "pl") |>
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
                 linetype = "dashed") +
    scale_color_hue()
  p2 = df |>
    filter(system == "ra" ) |>
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
                 linetype = "dashed") +
    scale_color_hue()
  p3 = ggpubr::ggarrange(p1, p2,  nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom", labels = LETTERS[1:2])
}

plot_sfig11 = function(){
  p = plot_prop_time_vs_height(pos = positions_examined)
  ggsave("figs/sm/sfig11_prop_time_vs_height.png",
         plot = p,
         width = width_2_col_cm,
         height = 10,
         units = "cm",
         bg = "white")
}
plot_sfig11()

#### Supp Figure 12: Water depths in time domain ####
plot_wd_time_domain = function(pos, plot_st = TRUE){
  # pos: positions for water depth
  # plot_st: plot systems tracts?
  stopifnot(any( pos %in% distances))
  st_sep_time_mod = c(-0.2,st_sep_time, 4.2)
  title_platform = "Platform"
  title_ramp = "Ramp"
  x_lab_title = "Elapsed Model Time [Myr]"
  y_lab_title = label_wd
  legend_label = "Distance from Shore"
  pos_interest = paste(pos, "km")
  
  col_names = c("wd", "t", "dist", "tag")
  df_wd = data.frame(matrix(ncol = length(col_names), nrow = 0))
  names(df_wd) = col_names
  for (case in cases){
    for (i in seq_along(pos)){
      wd = wd_comb[[case]][[which(pos[i] == distances)]]
      df_temp = data.frame(wd = wd$wd,
                           t = wd$t,
                           dist = rep(pos[i], length(wd$t)),
                           tag = rep(case, length(wd$t)))
      df_wd = rbind(df_wd, df_temp)
    }
  }
  df_wd$dist = factor(paste(df_wd$dist, "km"),
                      levels = paste(pos, "km"),
                      ordered = TRUE)
  
  p1 = df_wd |>
    filter(tag == "ra") |>
    ggplot( aes(x = t, y = wd, color = dist)) + 
    geom_line() +
    labs(title = title_ramp,
         x = x_lab_title,
         y = y_lab_title,
         col = legend_label) +
    scale_color_hue()
  if (plot_st){
    max_wd = df_wd |> filter(tag == "ra") |> summarise(m = max(wd)) |> as.numeric()
    df_text = data.frame(time = 0.5* (head(st_sep_time_mod, -1) + tail(st_sep_time_mod, -1)),
                         height = rep(0.8 * max_wd, length(st_labels_en)),
                         label = st_labels_en)
    p1 = p1 +
      geom_vline(xintercept = head(tail(st_sep_time_mod, -1),-1),
                 linetype = "dashed",
                 color = "black",
                 lwd = 0.5) +
      geom_text(data = df_text,
                aes(x = time,
                    y = height,
                    label = label),
                inherit.aes = FALSE,
                size = st_text_size)
  }
  
  p2 = df_wd |>
    filter(tag == "pl") |>
    ggplot( aes(x = t, y = wd, color = dist)) + 
    geom_line() +
    labs(title = title_platform,
         x = x_lab_title,
         y = y_lab_title,
         col = legend_label) +
    scale_color_hue()
  if (plot_st){
    max_wd = df_wd |> filter(tag == "pl") |> summarise(m = max(wd)) |> as.numeric()
    
    df_text = data.frame(time = 0.5* (head(st_sep_time_mod, -1) + tail(st_sep_time_mod, -1)),
                         height = rep(0.8 * max_wd, length(st_sep_time_mod)-1),
                         label = st_labels_en)
    
    p2 = p2 +
      geom_vline(xintercept = head(tail(st_sep_time_mod, -1),-1),
                 linetype = "dashed",
                 color = "black",
                 lwd = 0.5) +
      geom_text(data = df_text,
                aes(x = time,
                    y = height,
                    label = label),
                inherit.aes = FALSE,
                size = st_text_size)
  }
  p3 = ggpubr::ggarrange(p2, p1, ncol = 2, nrow = 1, labels = LETTERS[1:2], legend = "bottom", common.legend = TRUE)
  return(p3)
}

plot_sfig12 = function(){
  p = plot_wd_time_domain(pos = positions_examined)
  ggsave(filename = "figs/sm/sfig12_wd_time_domain.png", 
         plot = p,
         width = width_2_col_cm,
         height = 10,
         units = "cm",
         bg = "white")
}
plot_sfig12()


#### Supp Figure 13: Extinction scenarios ####
plot_sfig13_extinction_scenarios_time = function(){
  y_limit = 1700
  names = c("t", "rate", "ext_scen")
  df = data.frame(matrix(nrow = 0, ncol = length(names)))
  for (ext_scenario in names(ext_scen)){
    f = ext_scen[[ext_scenario]]
    resc = ext_rate_rescaler_ext
    if (ext_scenario == "constant"){resc = ext_rate_rescaler_const}
    df2 = data.frame(t = t,
                     rate = f(t) * resc,
                     ext_scen = rep(ext_scenario, length(t)))
    df = rbind(df, df2)
  }
  df$ext_scen = factor(df$ext_scen)
  
  st_sep_time_mod = c(-0.2, st_sep_time, 4.2)
  df_text = data.frame(time = 0.5* (head(st_sep_time_mod, -1) + tail(st_sep_time_mod, -1)),
                       height = rep(y_limit, length(st_labels_en)),
                       label = st_labels_en)
  plot_list = list()
  titles = c("constant" = "Constant Extinction Rate",
             "TST" = "Extinction in TST",
             "HST" = "Extinction in HST",
             "RST" = "Extinction in RST",
             "LST" = "Extinction in LST")
  for (i in seq_along(ext_scen)){
    case = names(ext_scen)[i]
    p = df |>
      filter(ext_scen == case) |>
      ggplot(aes(x = t, y = rate, color = ext_scen)) +
      geom_line(linewidth = 2) +
      labs(x = label_time_emt,
           y = label_ext_rate,
           title = titles[case]) +
      ylim( c(0,y_limit)) +
      geom_vline(xintercept = st_sep_time,
                 linetype = "dashed") +
      theme(legend.position = "none") +
      geom_text(data = df_text,
                aes(x = time,
                    y = height,
                    label = label),
                inherit.aes = FALSE,
                size = st_text_size)
    plot_list[[i]] = p
  }
  p = ggarrange(plotlist = plot_list,
                nrow = 3,
                ncol = 2,
                labels = LETTERS[1:5])
  ggsave(filename = "figs/sm/sfig13_ext_rates_time.png",
         plot = p,
         width = width_2_col_cm,
         height = 20,
         units = "cm",
         bg = "white")
}
plot_sfig13_extinction_scenarios_time()



#### Supp Figure 14 and 15: Signor lipps effect in ramp and platform ####

plot_sfig14 = function(seed){
  # plot signor-lipps effect for all positions in platform
  set.seed(seed)
  rates = rates_used
  pos = positions_examined
  case = "pl"
  pl_list = list()
  for (i in seq_along(pos)){
    pl_list[[i]] = plot_lo_by_rate(case = case,
                                   pos = pos[i],
                                   rates = rates,
                                   title = paste("Platform,", pos[i], "km from Shore"))
  }
  p = ggarrange( plotlist = pl_list, nrow = 3, ncol = 2,
                 labels = LETTERS[1:5])
  ggsave(filename = "figs/sm/sfig14_lo_in_platform.png",
         plot = p,
         width = width_2_col_cm,
         height = 20,
         units = "cm",
         bg = "white")
}
plot_sfig14(seed = seed_supp + 13)

plot_sfig15= function(seed){
  # plot signor-lipps effect for all positions in ramp
  set.seed(seed)
  rates = rates_used
  pos = positions_examined
  case = "ra"
  pl_list = list()
  for (i in seq_along(pos)){
    pl_list[[i]] = plot_lo_by_rate(case = case,
                                   pos = pos[i],
                                   rates = rates,
                                   title = paste("Ramp,", pos[i],"km from Shore"))
  }
  p = ggarrange( plotlist = pl_list, nrow = 3, ncol = 2,
                 labels = LETTERS[1:5])
  ggsave(filename = "figs/sm/sfig15_lo_in_ramp.png", 
         plot = p,
         width = width_2_col_cm,
         height = 20,
         units = "cm",
         bg = "white")
}
plot_sfig15(seed = seed_supp + 14)

#### Supp Figures 16 and 17: Extinction scenarios ####
plot_sfig16 = function(seed){
  set.seed(seed)
  pos = positions_examined
  case = "pl"
  rate = systract_com_rate
  pl_list = list()
  for (i in seq_along(pos)){
    pl_list[[i]] = plot_ext_scenario_comparison(rate = rate,
                                                dist = pos[i],
                                                case = case,
                                                title = paste("Platform",
                                                              pos[i],
                                                              "km from Shore"))
  }
  
  p = ggarrange(plotlist = pl_list, nrow = 3, ncol = 2,
                labels = LETTERS[1:5])
  ggsave(filename = "figs/sm/sfig16_ext_scen_platform.png",
         plot = p,
         width = width_2_col_cm,
         height = 20,
         units = "cm",
         bg = "white")
}
plot_sfig16(seed = seed_supp + 15)

plot_sfig17 = function(seed){
  set.seed(seed)
  pos = positions_examined
  case = "ra"
  rate = systract_com_rate
  pl_list = list()
  for (i in seq_along(pos)){
    pl_list[[i]] = plot_ext_scenario_comparison(rate = rate,
                                                dist = pos[i],
                                                case = case,
                                                title = paste("Ramp",
                                                              pos[i],
                                                              "km from Shore"))
  }
  
  p = ggarrange(plotlist = pl_list, nrow = 3, ncol = 2,
                labels = LETTERS[1:5])
  ggsave(filename = "figs/sm/sfig17_ext_scen_ramp.png",
         plot = p,
         width = width_2_col_cm,
         height = 20,
         units = "cm",
         bg = "white")
}

plot_sfig17(seed = seed_supp + 16)

#### Supp Figures 18 and 19: Spatial correlations of extinction scenarios ####

plot_sfig18 = function(seed){
  set.seed(seed)
  rate = spat_comp_rate
  case = "pl"
  ext_scenario_names = names(ext_scen)
  pos = positions_examined
  titles = c("constant" = "Platform, Constant Extinction Rate",
             "HST" = "Platform, Extinction in HST",
             "RST" = "Platform, Extinction in RST",
             "LST" = "Platform, Extinction in LST",
             "TST" = "Platform, Extinction in TST")
  plot_list = list()
  for (i in seq_along(ext_scenario_names)){
    plot_list[[i]] = plot_spat_corr_ext(rate = rate,
                                        case = case,
                                        ext_sce = ext_scenario_names[i],
                                        pos = pos) +
      labs(title = titles[ext_scenario_names[i]])
  }
  p = ggarrange(plotlist = plot_list,
                ncol = 2,
                nrow = 3,
                labels = LETTERS[1:5])
  ggsave(filename = "figs/sm/sfig18_spat_corr_platform.png",
         plot = p,
         width = width_2_col_cm,
         height = 20,
         units = "cm",
         bg = "white")
  
}
plot_sfig18(seed = seed_supp + 17)

plot_sfig19 = function(seed){
  set.seed(seed)
  rate = spat_comp_rate
  case = "ra"
  ext_scenario_names = names(ext_scen)
  pos = positions_examined
  plot_list = list()
  titles = c("constant" = "Ramp, Constant Extinction Rate",
             "HST" = "Ramp, Extinction in HST",
             "RST" = "Ramp, Extinction in RST",
             "LST" = "Ramp, Extinction in LST",
             "TST" = "Ramp, Extinction in TST")
  for (i in seq_along(ext_scenario_names)){
    plot_list[[i]] = plot_spat_corr_ext(rate = rate,
                                        case = case,
                                        ext_sce = ext_scenario_names[i],
                                        pos = pos) +
      labs(title = titles[ext_scenario_names[i]])
  }
  p = ggarrange(plotlist = plot_list,
                ncol = 2,
                nrow = 3,
                labels = LETTERS[1:5])
  ggsave(filename = "figs/sm/sfig19_spat_corr_ramp.png",
         plot = p,
         width = width_2_col_cm,
         height = 20,
         units = "cm",
         bg = "white")
  
}

plot_sfig19(seed = seed_supp + 18)



  #### Parameter scan: extinction scenario by rate ####

set.seed(seed_supp + 19)
plot_ext_scen_by_rate = function(rate, pos = c(3,10.5)){
  # compare ext scenarios at different locations
  p1 = plot_ext_scenario_comparison(rate = rate, 
                                    dist = pos[1],
                                    case = "ra",
                                    title = "Platform top")
  p2 = plot_ext_scenario_comparison(rate = rate, 
                                    dist = pos[2],
                                    case = "ra", title = "Platform distal slope")
  p3 = ggpubr::ggarrange(p1,p2, ncol = 2, nrow =1, common.legend = TRUE, legend = "bottom", labels = LETTERS[1:2])
  name = paste0("figs/lo_comp_different_ext_scenarios_rate_", rate, ".png")
  ggsave(filename = name,
         plot = p3)
}

for (rate in c(1,3,10,30,100,300,1000,3000)){
  plot_ext_scen_by_rate(rate = rate)
}


#### UNUSED FIGURE: comparison of extinction scenarios ####

plot_lo_comparison = function(){
  rates = rates_shortened
  p1 = plot_lo_by_rate(case = "pl",
                       pos = 3, 
                       rates = rates,
                       title = 'Platform top')
  p2 = plot_lo_by_rate(case = "pl",
                       pos = 10.5, 
                       rates = rates,
                       title = 'Platform ramp')
  p3 = ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, labels = LETTERS[1:2], common.legend = TRUE, legend = "bottom") 
  return(p3)
}

plot_figUNKNOWN = function(seed){
  set.seed(seed)
  p = plot_lo_comparison()
  ggsave("figs/ms_unknown_fig_no.png", p)
}
plot_figUNKNOWN(seed = seed_supp + 20)

#### Figures: spatial comparison of extinction scenarios ####
set.seed(seed_supp + 21)
if (!dir.exists("figs/spat_comp/")){dir.create("figs/spat_comp/", recursive = TRUE)}
for (case in cases){
  for (ext_scenario in names(ext_scen)){
    for (rate in c(5,10,30)){
      p = plot_spat_corr_ext(rate = rate, 
                             case = case, 
                             ext_sce = ext_scenario,
                             pos = positions_examined)
      p
      ggsave(filename = paste0("figs/spat_comp/spat_comp_",case, "_", rate, "_", ext_scenario,".png"))
    }
  }
}

