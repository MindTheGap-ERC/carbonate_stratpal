tag1 = "platform"
tag2 = "ramp"

adm_data_pl = read.csv(paste0("data/", tag1, "_adm.csv"))
wd_data_pl = read.csv(paste0("data/", tag1, "_wd.csv"))
sac_data_pl = read.csv(paste0("data/", tag1, "_sac.csv"))
adm_data_ra = read.csv(paste0("data/", tag2, "_adm.csv"))
wd_data_ra = read.csv(paste0("data/", tag2, "_wd.csv"))
sac_data_ra = read.csv(paste0("data/", tag2, "_sac.csv"))
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
  adm_list_pl[[i]] = admtools::tp_to_adm(t = adm_data_pl$time..Myr.,
                               h = adm_data_pl[, paste0("adm_", i, "..m.")])
  
  adm_list_ra[[i]] = admtools::tp_to_adm(t = adm_data_ra$time..Myr.,
                               h = adm_data_ra[, paste0("adm_", i, "..m.")])
}
adm_comb = list(ra = adm_list_ra, pl = adm_list_pl)

#### Sediment accumulation curves ####
sac_list_pl = list()
sac_list_ra = list()
for (i in 1:(length(sac_data_pl)-2)){
  sac_list_pl[[i]] = admtools::tp_to_sac(t = sac_data_pl$time..Myr.,
                                         h = sac_data_pl[, paste0("sac_", i, "..m.")])
  sac_list_ra[[i]] = admtools::tp_to_sac(t = sac_data_ra$time..Myr.,
                                         h = sac_data_ra[, paste0("sac_", i, "..m.")])
}
sac_comb = list(ra = sac_list_ra,
                pl = sac_list_pl)
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

# shorthand forms for ramp and platform
cases = c("ra", "pl")

distances_km = paste(distances, "km")

# boundaries between systems tracts
st_sep_time = seq(0.25, 3.75, by = 0.5)
# (enumerated) labels of system tracts
st_labels = c("TST", "HST", "RST","LST","TST","HST", "RST", "LST", "TST")
st_labels_en = c("TST1", "HST1", "RST1","LST1","TST2","HST2", "RST2", "LST2", "TST3")

# Extinction scenarios: 25-fold in- and decrease in extinction rate in through specific systems tract
#+ one constant scenario as baseline
ext_scen = list()
ext_scen[["HST"]] = approxfun(x = c(0.25,0.5, 0.75,2.25, 2.5, 2.75), y = c(1, 25,1, 1, 25, 1), rule = 2)
ext_scen[["RST"]] = approxfun(x = c(0.75, 1, 1.25, 2.75, 3, 3.25), y = c(1, 25,1, 1, 25, 1), rule = 2)
ext_scen[["LST"]] = approxfun(x = c(1.25, 1.5, 1.75, 3.25, 3.5, 3.75), y = c(1, 25,1, 1,25,1), rule = 2)
ext_scen[["TST"]] = approxfun(x = c(0, 0.25,1.75, 2, 2.25, 3.75,4), y = c(25,1,1, 25,1,1,25), rule = 2)
ext_scen[["constant"]] = approxfun(x = c(0, 4), y = c(1,1), rule = 2)

#### Examined positions ####
positions_examined = c(3,7.5,10.5,12,18)

#### Rates of foss occ ####
rates_used = c(2,5,10,22,46,100)
rates_shortened = c(3,10,30,100)
spat_comp_rate = 30 # rate used for the spatial comparison
systract_com_rate = 10 # sampling rate used for comparison of systems tracts

# no of extinctions for conditioning
n_ext = 1000 
ext_rate_rescaler_ext = n_ext / 16 # rescaling factor to arrive at mean 1000 ext over 4 Myr
ext_rate_rescaler_const = 250
#### Random seeds ####
# taken from random.org
seed_main =  622834 # for figs in ms
seed_supp =  732759 # for figs in supp material

##### Plotting constants ####
width_2_col_cm = 18
width_1_col_cm = 9
st_text_size = 2.5 # text size of system tract labels (passed to geom_text)

## labels ##
label_height_absolute = "Stratigraphic Height [m]"
label_time_emt = "Elapsed Model Time [Myr]"
label_wd = "Water Depth [m]"
label_ext_rate = "Extinction Rate [#/Myr]"
q