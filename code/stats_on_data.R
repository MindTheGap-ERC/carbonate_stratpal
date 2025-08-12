library(admtools)

tag = "carbonate_stratpal_1"

adm_data = read.csv(paste0("data/", tag, "_adm.csv"))

adm_list = list()
for (i in 1:(length(adm_data)-2)){
  adm_list[[i]] = tp_to_adm(t = adm_data$time..Myr.,
                            h = adm_data[, paste0("adm_", i, "..m.")])
}

plot(adm_list[[1]])
plot(adm_list[[2]])
plot(adm_list[[3]])
plot(adm_list[[4]])

wd_data = read.csv(paste0("data/", tag, "_wd.csv"))
wd = list()
for (i in 1:(length(wd_data)-2)){
  wd[[i]] = list(t = wd_data$time..Myr.,
                 wd = wd_data[, paste0("wd_", i, "..m.")])
}
for (i in 1:4){
  plot(wd[[i]]$t, wd[[i]]$wd)
}

# sea level
period1 = 2.0 # Myr
amplitude1 = 15.0 # m
period2 = 0.2 # Myr
amplitude2 = 2.5 # m

sl = list(t = adm_data$time..Myr.,
          sl = amplitude1 * sin(2 * pi * t/ period1) + amplitude2 * sin(2 * pi * t / period2))
me = read_toml(paste0("data/", tag, ".toml"))
######################################### 4km         systems tract isolated
h4 = data_kitten$adm8..m. 
plot(h4,type='l')  
h4_mod = h4[t >= 1.25 & t <= 3.25]                          # modify to show 1.25-3.25 Myr
h4_mod_shifted = h4_mod - min(h4_mod)                       # shift to 0-2 Myr
adm4_mod = tp_to_adm(t = t_mod, h4_mod)                    
adm_4km = tp_to_adm(t = t_mod_shifted, h4_mod_shifted)
plot(adm_4km,lwd_acc = 2,lwd_destr = 0)

max_height_4km = max(get_height(adm_4km,destructive=FALSE))             # maximum height at 4 km
                                       
h4_LST = h4[ t >= 1.25 & t <= 1.75]                       # LST 4 km
h4_LST_shifted = h4_LST - min(h4_LST)
adm_LST_4km = tp_to_adm(t = LST_shifted, h4_LST_shifted)
plot(adm_LST_4km)
title('4km LST adm')
LST_height_4km = get_height(adm_4km,0.5)

h4_TST = h4[ t >= 1.75 & t <= 2.25]                       # TST 4 km
h4_TST_shifted = h4_TST - min(h4_LST)
adm_TST_4km = tp_to_adm(t = TST_shifted, h4_TST_shifted)
plot(adm_TST_4km)
title('4km TST adm')
TST_height_4km = get_height(adm_4km,1) - get_height(adm_4km,0.5)

h4_HST = h4[ t >= 2.25 & t <= 2.75]                       # HST 4 km
h4_HST_shifted = h4_HST - min(h4_LST)
adm_HST_4km = tp_to_adm(t = HST_shifted, h4_HST_shifted)
plot(adm_HST_4km)
title('4km HST adm')
HST_height_4km = get_height(adm_4km,1.5,destructive=FALSE) - get_height(adm_4km,1)

h4_FSST = h4[ t >= 2.75 & t <= 3.25]                      # FSST 4 km
h4_FSST_shifted = h4_FSST - min(h4_LST)
adm_FSST_4km = tp_to_adm(t = FSST_shifted, h4_FSST_shifted)
plot(adm_FSST_4km)
title('4km FSST adm')
FSST_height_4km = max_height_4km - get_height(adm_4km,1.5,destructive=FALSE)

adm_st_4km <- list(adm_LST_4km,adm_TST_4km,adm_HST_4km,adm_FSST_4km)

w_hiat_no_4km <- c()                          # Number of gaps per System Track at 4km.
for(hiat_no_4km in adm_st_4km){
  v_hiat_no_4km <- get_hiat_no(hiat_no_4km)
  w_hiat_no_4km <- c(w_hiat_no_4km,v_hiat_no_4km)
}
plot(w_hiat_no_4km,
     main = "number of hiatus per system tract in 4km adm",
     xlab = "LST,TST,HST,FSST",
     ylab = "number of hiatus",
     pch = c(16))
get_completeness(adm_4km)
get_hiat_no(adm_LST_4km)
get_hiat_no(adm_TST_4km)
get_hiat_no(adm_HST_4km)
get_hiat_no(adm_FSST_4km)
get_completeness(adm_LST_4km)
get_completeness(adm_TST_4km)
get_completeness(adm_HST_4km)
get_completeness(adm_FSST_4km)
tot_hiat_dur_LST_4km <- sum(get_hiat_duration(adm_LST_4km))
tot_hiat_dur_HST_4km <- sum(get_hiat_duration(adm_TST_4km))
tot_hiat_dur_TST_4km <- sum(get_hiat_duration(adm_HST_4km))
tot_hiat_dur_FSST_4km <- sum(get_hiat_duration(adm_FSST_4km))

###################################9km      system tracts & sed environments isolated 
h9 = data_kitten$adm18..m.                                  # Entire runtime (9 Myr)
plot(h9,type='l')  
h9_mod = h9[t >= 1.25 & t <= 3.25]                          # modify to show 1.25-3.25 Myr
h9_mod_shifted = h9_mod - min(h9_mod)                       # shift to 0-2 Myr
adm9_mod = tp_to_adm(t = t_mod, h9_mod)                    
adm_9km = tp_to_adm(t = t_mod_shifted, h9_mod_shifted)     
plot(adm_9km,lwd_acc = 2,lwd_destr = 0)
  
basin_boundary_9km = abline(v=(0.52))                                     # end of basin
slope_boundary_9km = abline(v=(0.66))                                     # end of slope
slope_dur_9km = 0.66-0.52                                                 # duration of slope          
interior_dur_9km = 2-0.66                                                 # duration of platform interior
basin_height_9km = get_height(adm_9km,0.52)                               # height of the basin at 9 km
slope_height_9km = get_height(adm_9km,0.66)-basin_height_9km              # height of the slope at 9 km
max_height_9km = max(get_height(adm_9km,destructive=FALSE))               # maximum height at 9 km
interior_height_9km = max_height_9km-get_height(adm_9km,0.66)             # height of platform interior

h9_LST = h9[ t >= 1.25 & t <= 1.75]                       #LST 9km
h9_LST_shifted = h9_LST - min(h9_LST)
adm_LST_9km = tp_to_adm(t = LST_shifted, h9_LST_shifted)
plot(adm_LST_9km)
title('9km LST adm')
LST_height_9km = get_height(adm_9km,0.5,destructive=FALSE)

h9_TST = h9[ t >= 1.75 & t <= 2.25]                       #TST 9km
h9_TST_shifted = h9_TST - min(h9_LST)
adm_TST_9km = tp_to_adm(t = TST_shifted, h9_TST_shifted)
plot(adm_TST_9km)
title('9km TST adm')  

h9_HST = h9[ t >= 2.25 & t <= 2.75]                       #HST 9km
h9_HST_shifted = h9_HST - min(h9_LST)
adm_HST_9km = tp_to_adm(t = HST_shifted, h9_HST_shifted)
plot(adm_HST_9km)
title('9km HST adm')

h9_FSST = h9[ t >= 2.75 & t <= 3.25]                      #FSST 9km
h9_FSST_shifted = h9_FSST - min(h9_LST)
adm_FSST_9km = tp_to_adm(t = FSST_shifted, h9_FSST_shifted)
plot(adm_FSST_9km)
title('9km FSST adm')

adm_st_9km <- list(adm_LST_9km,adm_TST_9km,adm_HST_9km,adm_FSST_9km)

w_hiat_no_9km <- c()                          # Number of hiatuses per System Track at 9km.
for(hiat_no_9km in adm_st_9km){
  v_hiat_no_9km <- get_hiat_no(hiat_no_9km)
  w_hiat_no_9km <- c(w_hiat_no_9km,v_hiat_no_9km)
}
plot(w_hiat_no_9km,
     main = "number of hiatus per system tract in 9km adm",
     xlab = "LST,TST,HST,FSST",
     ylab = "number of hiatus",
     pch = c(16))
# stats per systems tract
get_completeness(adm_LST_9km)
get_completeness(adm_TST_9km)
get_completeness(adm_HST_9km)
get_completeness(adm_FSST_9km)
get_completeness(adm_slope_9km)
sum(get_hiat_duration(adm_LST_9km))
sum(get_hiat_duration(adm_TST_9km))
sum(get_hiat_duration(adm_HST_9km))
sum(get_hiat_duration(adm_FSST_9km))

###################################11km      system tracts & sed environments isolated 
h11 = data_kitten$adm22..m.                                  # Entire runtime (11 Myr)
plot(h11, type = 'l')  
h11_mod = h11[t >= 1.25 & t <= 3.25]                          # modify to show 1.25–3.25 Myr
h11_mod_shifted = h11_mod - min(h11_mod)                     # shift to 0–2 Myr
adm11_mod = tp_to_adm(t = t_mod, h11_mod)                    
adm_11km = tp_to_adm(t = t_mod_shifted, h11_mod_shifted)     
plot(adm_11km, lwd_acc = 2, lwd_destr = 0)
  
basin_boundary_11km = abline(v=(1.6))                                     # end of basin
slope_boundary_11km = abline(v=(1.78))                                    # end of slope
slope_dur_11km = 1.78-1.6                                                 # duration of slope
time_to_strat(1.78,adm_11km); time_to_strat(1.6,adm_11km)                 # in the depth domain
interior_dur_11km = 2-1.78                                                # duration of platform interior
basin_height_11km = get_height(adm_11km,1.6)                              # height of the basin
slope_height_11km = get_height(adm_11km,1.78)-basin_height_11km           # height of the slope
max_height_11km = max(get_height(adm_11km,destructive=FALSE))             # maximum height at 11 km
interior_height_11km = max_height_11km-get_height(adm_11km,1.78)          # height of platform interior


h11_LST = h11[ t >= 1.25 & t <= 1.75]                       #LST 11km
h11_LST_shifted = h11_LST - min(h11_LST)
adm_LST_11km = tp_to_adm(t = LST_shifted, h11_LST_shifted)
plot(adm_LST_11km)
title('11km LST adm')
LST_height_11km = get_height(adm_11km,0.5,destructive=FALSE)

h11_TST = h11[ t >= 1.75 & t <= 2.25]                       #TST 11km
h11_TST_shifted = h11_TST - min(h11_LST)
adm_TST_11km = tp_to_adm(t = TST_shifted, h11_TST_shifted)
plot(adm_TST_11km)
title('11km TST adm')
TST_height_11km = get_height(adm_11km,1) - get_height(adm_11km,0.5)

h11_HST = h11[ t >= 2.25 & t <= 2.75]                       #HST 11km
h11_HST_shifted = h11_HST - min(h11_LST)
adm_HST_11km = tp_to_adm(t = HST_shifted, h11_HST_shifted)
plot(adm_HST_11km)
title('11km HST adm')
HST_height_11km = get_height(adm_11km,1.5) - get_height(adm_11km,1)

h11_FSST = h11[ t >= 2.75 & t <= 3.25]                      #FSST 11km
h11_FSST_shifted = h11_FSST - min(h11_FSST)
adm_FSST_11km = tp_to_adm(t = FSST_shifted, h11_FSST_shifted)
plot(adm_FSST_11km)
title('11km FSST adm')
FSST_height_11km = max_height_11km - get_height(adm_11km,1.5,destructive=FALSE)

w_hiat_no_11km <- c()                          # Number of hiatuses per Systems Track at 11km.
for(hiat_no_11km in adm_st_11km){
  v_hiat_no_11km <- get_hiat_no(hiat_no_11km)
  w_hiat_no_11km <- c(w_hiat_no_11km,v_hiat_no_11km)
}
plot(w_hiat_no_11km,
     main = "number of hiatus per system tract in 11km adm",
     xlab = "LST,TST,HST,FSST",
     ylab = "number of hiatus",
     pch = c(16))

get_completeness(adm_LST_11km)
get_completeness(adm_TST_11km)
get_completeness(adm_HST_11km)
get_completeness(adm_FSST_11km)
get_completeness(adm_slope_11km)
sum(get_hiat_duration(adm_LST_11km))
sum(get_hiat_duration(adm_TST_11km))
sum(get_hiat_duration(adm_HST_11km))
sum(get_hiat_duration(adm_FSST_11km))
get_hiat_no(adm_LST_11km)
get_hiat_no(adm_TST_11km)
get_hiat_no(adm_HST_11km)
get_hiat_no(adm_FSST_11km)

###################################15km      system tracts isolated 
###  15 km       adm                                       
h15 = data_kitten$adm30..m.                                   # Entire runtime (15 Myr)
plot(h15, type = 'l')  
h15_mod = h15[t >= 1.25 & t <= 3.25]                          # modify to show 1.25–3.25 Myr
h15_mod_shifted = h15_mod - min(h15_mod)                      # shift to 0–2 Myr
adm15_mod = tp_to_adm(t = t_mod, h15_mod)                    
adm_15km = tp_to_adm(t = t_mod_shifted, h15_mod_shifted)     
plot(adm_15km, lwd_acc = 2, lwd_destr = 0)

max_height_15km = max(get_height(adm_15km,destructive=FALSE))   # maximum height at 15 km

h15_LST = h15[ t >= 1.25 & t <= 1.75]                       #LST 15km
h15_LST_shifted = h15_LST - min(h15_LST)
adm_LST_15km = tp_to_adm(t = LST_shifted, h15_LST_shifted)
plot(adm_LST_15km);
title('15km LST adm')
LST_height_15km = get_height(adm_15km,0.5,destructive=FALSE)

h15_TST = h15[ t >= 1.75 & t <= 2.25]                       #TST 15km
h15_TST_shifted = h15_TST - min(h15_LST)
adm_TST_15km = tp_to_adm(t = TST_shifted, h15_TST_shifted)
plot(adm_TST_15km)
title('15km TST adm')

h15_HST = h15[ t >= 2.25 & t <= 2.75]                       #HST 15km
h15_HST_shifted = h15_HST - min(h15_LST)
adm_HST_15km = tp_to_adm(t = HST_shifted, h15_HST_shifted)
plot(adm_HST_15km)
title('15km HST adm')

h15_FSST = h15[ t >= 2.75 & t <= 3.25]                      #FSST 15km
h15_FSST_shifted = h15_FSST - min(h15_LST)
adm_FSST_15km = tp_to_adm(t = FSST_shifted, h15_FSST_shifted)
plot(adm_FSST_15km)
title('15km FSST adm')

adm_st_15km <- list(adm_LST_15km,adm_TST_15km,adm_HST_15km,adm_FSST_15km)

w_hiat_no_15km <- c()                          # Number of gaps per System Track at 15km.
for(hiat_no_15km in adm_st_15km){
  v_hiat_no_15km <- get_hiat_no(hiat_no_15km)
  w_hiat_no_15km <- c(w_hiat_no_15km,v_hiat_no_15km)
}
plot(w_hiat_no_15km,
     main = "number of hiatus per system tract in 15km adm",
     xlab = "LST,TST,HST,FSST",
     ylab = "number of hiatus",
     pch = c(16))

get_completeness(adm_LST_15km)
get_completeness(adm_TST_15km)
get_completeness(adm_HST_15km)
get_completeness(adm_FSST_15km)
sum(get_hiat_duration(adm_LST_15km))
sum(get_hiat_duration(adm_TST_15km))
sum(get_hiat_duration(adm_HST_15km))
sum(get_hiat_duration(adm_FSST_15km))

# Compare completeness per sedimentary environment
# First isolate 1.25 Myr - 3.25 Myr per ADM

###  1 km       adm
h1 = data_kitten$adm2..m. 
plot(h1,type='l')  
h1_mod = h1[t >= 1.25 & t <= 3.25]                          # modify to show 1.25-3.25 Myr
h1_mod_shifted = h1_mod - min(h1_mod)                       # shift to 0-2 Myr
adm1_mod = tp_to_adm(t = t_mod, h1_mod)                    
adm_1km = tp_to_adm(t = t_mod_shifted, h1_mod_shifted)
plot(adm_1km,lwd_acc = 2,lwd_destr = 0)

###  2 km       adm
h2 = data_kitten$adm4..m. 
plot(h2,type='l')  
h2_mod = h2[t >= 1.25 & t <= 3.25]                          # modify to show 1.25-3.25 Myr
h2_mod_shifted = h2_mod - min(h2_mod)                       # shift to 0-2 Myr
adm2_mod = tp_to_adm(t = t_mod, h2_mod)                    
adm_2km = tp_to_adm(t = t_mod_shifted, h2_mod_shifted)
plot(adm_2km,lwd_acc = 2,lwd_destr = 0)

###  3 km       adm
h3 = data_kitten$adm6..m. 
plot(h3,type='l')  
h3_mod = h3[t >= 1.25 & t <= 3.25]                          # modify to show 1.25-3.25 Myr
h3_mod_shifted = h3_mod - min(h3_mod)                       # shift to 0-2 Myr
adm3_mod = tp_to_adm(t = t_mod, h3_mod)                    
adm_3km = tp_to_adm(t = t_mod_shifted, h3_mod_shifted)
plot(adm_3km,lwd_acc = 2,lwd_destr = 0)

###  4 km       adm
h4 = data_kitten$adm8..m. 
plot(h4,type='l')  
h4_mod = h4[t >= 1.25 & t <= 3.25]                          # modify to show 1.25-3.25 Myr
h4_mod_shifted = h4_mod - min(h4_mod)                       # shift to 0-2 Myr
adm4_mod = tp_to_adm(t = t_mod, h4_mod)                    
adm_4km = tp_to_adm(t = t_mod_shifted, h4_mod_shifted)
plot(adm_4km,lwd_acc = 2,lwd_destr = 0)

###  5 km       adm
h5 = data_kitten$adm10..m. 
plot(h5,type='l')  
h5_mod = h5[t >= 1.25 & t <= 3.25]                          # modify to show 1.25-3.25 Myr
h5_mod_shifted = h5_mod - min(h5_mod)                       # shift to 0-2 Myr
adm5_mod = tp_to_adm(t = t_mod, h5_mod)                    
adm_5km = tp_to_adm(t = t_mod_shifted, h5_mod_shifted)
plot(adm_5km,lwd_acc = 2,lwd_destr = 0)
                              
###  6 km       adm
h6 = data_kitten$adm12..m. 
plot(h6,type='l')                                           # Entire runtime (6 Myr)
h6_mod = h6[t >= 1.25 & t <= 3.25]                          # modify to show 1.25-3.25 Myr
h6_mod_shifted = h6_mod - min(h6_mod)                       # shift to 0-2 Myr
adm6_mod = tp_to_adm(t = t_mod, h6_mod)                    
adm_6km = tp_to_adm(t = t_mod_shifted, h6_mod_shifted)
plot(adm_6km,lwd_acc = 2,lwd_destr = 0)

###  7 km       adm                                       
h7 = data_kitten$adm14..m.                                  # Entire runtime (7 Myr)
plot(h7,type='l')
h7_mod = h7[t >= 1.25 & t <= 3.25]                          # modify to show 1.25-3.25 Myr
h7_mod_shifted = h7_mod - min(h7_mod)                       # shift to 0-2 Myr
adm7_mod = tp_to_adm(t = t_mod, h7_mod)                    
adm_7km = tp_to_adm(t = t_mod_shifted, h7_mod_shifted)
plot(adm_7km,lwd_acc = 2,lwd_destr = 0)

###  8 km       adm                                       
h8 = data_kitten$adm16..m.                                  # Entire runtime (8 Myr)
plot(h8,type='l')  
h8_mod = h8[t >= 1.25 & t <= 3.25]                          # modify to show 1.25-3.25 Myr
h8_mod_shifted = h8_mod - min(h8_mod)                       # shift to 0-2 Myr
adm8_mod = tp_to_adm(t = t_mod, h8_mod)                    
adm_8km = tp_to_adm(t = t_mod_shifted, h8_mod_shifted)     
plot(adm_8km,lwd_acc = 2,lwd_destr = 0)

###  9 km       adm                                       
h9 = data_kitten$adm18..m.                                  # Entire runtime (9 Myr)
plot(h9,type='l')  
h9_mod = h9[t >= 1.25 & t <= 3.25]                          # modify to show 1.25-3.25 Myr
h9_mod_shifted = h9_mod - min(h9_mod)                       # shift to 0-2 Myr
adm9_mod = tp_to_adm(t = t_mod, h9_mod)                    
adm_9km = tp_to_adm(t = t_mod_shifted, h9_mod_shifted)     
plot(adm_9km,lwd_acc = 2,lwd_destr = 0)

###  10 km       adm                                       
h10 = data_kitten$adm20..m.                                  # Entire runtime (10 Myr)
plot(h10,type='l')  
h10_mod = h10[t >= 1.25 & t <= 3.25]                          # modify to show 1.25-3.25 Myr
h10_mod_shifted = h10_mod - min(h10_mod)                       # shift to 0-2 Myr
adm10_mod = tp_to_adm(t = t_mod, h10_mod)                    
adm_10km = tp_to_adm(t = t_mod_shifted, h10_mod_shifted)     
plot(adm_10km,lwd_acc = 2,lwd_destr = 0)

###  11 km       adm 
h11 = data_kitten$adm22..m.                                  # Entire runtime (11 Myr)
plot(h11, type = 'l')  
h11_mod = h11[t >= 1.25 & t <= 3.25]                          # modify to show 1.25–3.25 Myr
h11_mod_shifted = h11_mod - min(h11_mod)                     # shift to 0–2 Myr
adm11_mod = tp_to_adm(t = t_mod, h11_mod)                    
adm_11km = tp_to_adm(t = t_mod_shifted, h11_mod_shifted)     
plot(adm_11km, lwd_acc = 2, lwd_destr = 0)

###  12 km       adm                                       
h12 = data_kitten$adm24..m.                                  # Entire runtime (12 Myr)
plot(h12, type = 'l')  
h12_mod = h12[t >= 1.25 & t <= 3.25]                          # modify to show 1.25–3.25 Myr
h12_mod_shifted = h12_mod - min(h12_mod)                     # shift to 0–2 Myr
adm12_mod = tp_to_adm(t = t_mod, h12_mod)                    
adm_12km = tp_to_adm(t = t_mod_shifted, h12_mod_shifted)     
plot(adm_12km, lwd_acc = 2, lwd_destr = 0)

###  13 km       adm                                       
h13 = data_kitten$adm26..m.                                  # Entire runtime (13 Myr)
plot(h13, type = 'l')  
h13_mod = h13[t >= 1.25 & t <= 3.25]                          # modify to show 1.25–3.25 Myr
h13_mod_shifted = h13_mod - min(h13_mod)                     # shift to 0–2 Myr
adm13_mod = tp_to_adm(t = t_mod, h13_mod)                    
adm_13km = tp_to_adm(t = t_mod_shifted, h13_mod_shifted)     
plot(adm_13km, lwd_acc = 2, lwd_destr = 0)

###  14 km       adm                                       
h14 = data_kitten$adm28..m.                                  # Entire runtime (14 Myr)
plot(h14, type = 'l')  
h14_mod = h14[t >= 1.25 & t <= 3.25]                          # modify to show 1.25–3.25 Myr
h14_mod_shifted = h14_mod - min(h14_mod)                     # shift to 0–2 Myr
adm14_mod = tp_to_adm(t = t_mod, h14_mod)                    
adm_14km = tp_to_adm(t = t_mod_shifted, h14_mod_shifted)     
plot(adm_14km, lwd_acc = 2, lwd_destr = 0)

###  15 km       adm                                       
h15 = data_kitten$adm30..m.                                  # Entire runtime (15 Myr)
plot(h15, type = 'l')  
h15_mod = h15[t >= 1.25 & t <= 3.25]                          # modify to show 1.25–3.25 Myr
h15_mod_shifted = h15_mod - min(h15_mod)                     # shift to 0–2 Myr
adm15_mod = tp_to_adm(t = t_mod, h15_mod)                    
adm_15km = tp_to_adm(t = t_mod_shifted, h15_mod_shifted)     
plot(adm_15km, lwd_acc = 2, lwd_destr = 0)

###  16 km       adm 
h16 = data_kitten$adm32..m.                                  # Entire runtime (16 Myr)
plot(h16, type = 'l')  
h16_mod = h16[t >= 1.25 & t <= 3.25]                          # modify to show 1.25–3.25 Myr
h16_mod_shifted = h16_mod - min(h16_mod)                     # shift to 0–2 Myr
adm16_mod = tp_to_adm(t = t_mod, h16_mod)                    
adm_16km = tp_to_adm(t = t_mod_shifted, h16_mod_shifted)     
plot(adm_16km, lwd_acc = 2, lwd_destr = 0)

###  17 km       adm 
h17 = data_kitten$adm34..m.                                  # Entire runtime (17 Myr)
plot(h17, type = 'l')  
h17_mod = h17[t >= 1.25 & t <= 3.25]                          # modify to show 1.25–3.25 Myr
h17_mod_shifted = h17_mod - min(h17_mod)                     # shift to 0–2 Myr
adm17_mod = tp_to_adm(t = t_mod, h17_mod)                    
adm_17km = tp_to_adm(t = t_mod_shifted, h17_mod_shifted)     
plot(adm_17km, lwd_acc = 2, lwd_destr = 0)

###  18 km       adm
h18 = data_kitten$adm36..m.                                  # Entire runtime (18 Myr)
plot(h18, type = 'l')  
h18_mod = h18[t >= 1.25 & t <= 3.25]                          # modify to show 1.25–3.25 Myr
h18_mod_shifted = h18_mod - min(h18_mod)                     # shift to 0–2 Myr
adm18_mod = tp_to_adm(t = t_mod, h18_mod)                    
adm_18km = tp_to_adm(t = t_mod_shifted, h18_mod_shifted)     
plot(adm_18km, lwd_acc = 2, lwd_destr = 0)

###  19 km       adm
h19 = data_kitten$adm38..m.                                  # Entire runtime (19 Myr)
plot(h19, type = 'l')  
h19_mod = h19[t >= 1.25 & t <= 3.25]                          # modify to show 1.25–3.25 Myr
h19_mod_shifted = h19_mod - min(h19_mod)                     # shift to 0–2 Myr
adm19_mod = tp_to_adm(t = t_mod, h19_mod)                    
adm_19km = tp_to_adm(t = t_mod_shifted, h19_mod_shifted)     
plot(adm_19km, lwd_acc = 2, lwd_destr = 0)

###  20 km       adm
h20 = data_kitten$adm40..m.                                  # Entire runtime (20 Myr)
plot(h20, type = 'l')  
h20_mod = h20[t >= 1.25 & t <= 3.25]                          # modify to show 1.25–3.25 Myr
h20_mod_shifted = h20_mod - min(h20_mod)                     # shift to 0–2 Myr
adm20_mod = tp_to_adm(t = t_mod, h20_mod)                    
adm_20km = tp_to_adm(t = t_mod_shifted, h20_mod_shifted)     
plot(adm_20km, lwd_acc = 2, lwd_destr = 0)

#Incompleteness percentage per ADM #################
w_incompleteness <- c()
for(incompleteness in adm){
  v_incompleteness <- get_incompleteness(incompleteness)
  w_incompleteness <- c(w_incompleteness,v_incompleteness)
}
plot(w_incompleteness,
     main = "incompleteness per adm 1 to 20km",
     xlab = "km",
     ylab = "incompleteness percentage", 
     pch = c(16))

#Completeness percentage per adm #################
w_completeness <- c()
for(completeness in adm){
  v_completeness <- get_completeness(completeness)
  w_completeness <- c(w_completeness,v_completeness)
}
plot(w_completeness,
     main = "completeness per adm 1 to 20km",
     xlab = "km",
     ylab = "percentage",
     pch = c(16))

#completeness average _ interior
(get_completeness(adm_1km)+get_completeness(adm_2km)+get_completeness(adm_3km)+get_completeness(adm_4km)+get_completeness(adm_5km)+get_completeness(adm_6km)+get_completeness(adm_7km))/7
#completeness average _ all states
(get_completeness(adm_8km)+get_completeness(adm_9km)+get_completeness(adm_10km)+get_completeness(adm_11km))/4
#completeness average _ basin
(get_completeness(adm_12km)+get_completeness(adm_13km)+get_completeness(adm_14km)+get_completeness(adm_15km)+get_completeness(adm_16km)+get_completeness(adm_17km)+get_completeness(adm_18km)+get_completeness(adm_19km)+get_completeness(adm_20km))/9

#Number of hiatuses per ADM ##############################
w_hiat_no <- c()
for(hiat_no in adm){
  v_hiat_no <- get_hiat_no(hiat_no)
  w_hiat_no <- c(w_hiat_no,v_hiat_no)
}
plot(w_hiat_no,
     main = "number of gaps per adm 1 to 20km",
     xlab = "km",
     ylab = "number of gaps",
     pch = c(16))
#hiat_no average _ platform interior
(get_hiat_no(adm_1km)+get_hiat_no(adm_2km)+get_hiat_no(adm_3km)+get_hiat_no(adm_4km)+get_hiat_no(adm_5km)+get_hiat_no(adm_6km)+get_hiat_no(adm_7km))/7
#hiat_no average _ all states of platform
(get_hiat_no(adm_8km)+get_hiat_no(adm_9km)+get_hiat_no(adm_10km)+get_hiat_no(adm_11km))/4
#hiat_no average _ basin
(get_hiat_no(adm_13km)+get_hiat_no(adm_14km)+get_hiat_no(adm_15km)+get_hiat_no(adm_16km)+get_hiat_no(adm_17km)+get_hiat_no(adm_18km)+get_hiat_no(adm_19km)+get_hiat_no(adm_20km))/8

#Average duration of hiatus per ADM ##################
w_hiat_dur <- c()
for(hiat_dur in adm){
  v_hiat_dur <- mean(get_hiat_duration(hiat_dur))
  w_hiat_dur <- c(w_hiat_dur,v_hiat_dur)
}
plot(w_hiat_dur,
     main = "average duration of hiatus per adm 1 to 20km",
     xlab = "km",
     ylab = "Myr",
     pch = c(16))
get_hiat_duration(adm_4km)                # List of duration of all hiatuses at given distance
get_hiat_duration(adm_9km)
get_hiat_duration(adm_11km)
get_hiat_duration(adm_15km)
#hiat_duration average _ platform interior
(mean(get_hiat_duration(adm_1km))+mean(get_hiat_duration(adm_2km))+mean(get_hiat_duration(adm_3km))+mean(get_hiat_duration(adm_4km))+mean(get_hiat_duration(adm_5km))+mean(get_hiat_duration(adm_4km))+mean(get_hiat_duration(adm_5km)))/7
#hiat_duration average _ all stages of platform
(mean(get_hiat_duration(adm_8km))+mean(get_hiat_duration(adm_9km))+mean(get_hiat_duration(adm_10km))+mean(get_hiat_duration(adm_11km)))/4
#hiat_duration average basin
(mean(get_hiat_duration(adm_13km))+mean(get_hiat_duration(adm_14km))+mean(get_hiat_duration(adm_15km))+mean(get_hiat_duration(adm_16km))+mean(get_hiat_duration(adm_17km))+mean(get_hiat_duration(adm_18km))+mean(get_hiat_duration(adm_19km))+mean(get_hiat_duration(adm_20km)))/8


#Total hiatus duration per adm #######################
w_hiat_sum <- c()
for(hiat_sum in adm){
  v_hiat_sum <- sum(get_hiat_duration(hiat_sum))
  w_hiat_sum <- c(w_hiat_sum,v_hiat_sum)
}
plot(w_hiat_sum,
     main = "total duration of gaps per adm 1 to 20km",
     xlab = "km",
     ylab = "Myr",
     pch = c(16))

#Total thickness per ADM ######################################
w_thickness <- c()
for(thickness in adm){
  v_thickness <- get_total_thickness(thickness)
  w_thickness <- c(w_thickness,v_thickness)
}
plot(w_thickness,
     main = "total thickness",
     xlab = "ADM km",
     ylab = "meter",
     pch = c(16))   


