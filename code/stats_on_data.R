data_kitten = read.csv("data/alcap-example7_adm.csv")
water_depth = read.csv("data/alcap-example7_wd.csv")

library(admtools)
library(StratPal)

#ADM's: 500m:500m:15km
adm_1km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm2..m.)
adm_2km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm4..m.)
adm_3km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm6..m.)
adm_4km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm8..m.)
adm_5km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm10..m.)
adm_6km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm12..m.)
adm_7km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm14..m.)
adm_8km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm16..m.)
adm_9km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm18..m.)
adm_10km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm20..m.)
adm_11km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm22..m.)
adm_12km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm24..m.)
adm_13km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm26..m.)
adm_14km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm28..m.)
adm_15km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm30..m.)
adm_16km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm32..m.)
adm_17km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm34..m.)
adm_18km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm36..m.)
adm_19km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm38..m.)
adm_20km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm40..m.)

adm <- list(adm_1km,adm_2km,adm_3km,adm_4km,adm_5km,adm_6km,adm_7km,adm_8km,
            adm_9km,adm_10km,adm_11km,adm_12km,adm_13km,adm_14km,adm_15km,
            adm_16km,adm_17km,adm_18km,adm_19km,adm_20km)
df<-adm
t = data_kitten$time..Myr.
######################################### 4km         systems tract isolated
LST = t[ t >= 1.25 & t <= 1.75]                          #LST
LST_shifted = LST - min(LST)

h4 = data_kitten$adm8..m.                            
h4_LST = h4[ t >= 1.25 & t <= 1.75]
h4_LST_shifted = h4_LST - min(h4_LST)
adm_LST_4km = tp_to_adm(t = LST_shifted, h4_LST_shifted)
plot(adm_LST_4km)
title('4km LST adm')

TST = t[ t >= 1.75 & t <= 2.25]                          #TST
TST_shifted = TST - min(LST)

h4_TST = h4[ t >= 1.75 & t <= 2.25]
h4_TST_shifted = h4_TST - min(h4_LST)
adm_TST_4km = tp_to_adm(t = TST_shifted, h4_TST_shifted)
plot(adm_TST_4km)
title('4km TST adm')

HST = t[ t >= 2.25 & t <= 2.75]                          #HST
HST_shifted = HST - min(LST)

h4_HST = h4[ t >= 2.25 & t <= 2.75]
h4_HST_shifted = h4_HST - min(h4_LST)
adm_HST_4km = tp_to_adm(t = HST_shifted, h4_HST_shifted)
plot(adm_HST_4km)
title('4km HST adm')

FSST = t[ t >= 2.75 & t <= 3.25]                          #FSST
FSST_shifted = FSST - min(LST)

h4_FSST = h4[ t >= 2.75 & t <= 3.25]
h4_FSST_shifted = h4_FSST - min(h4_LST)
adm_FSST_4km = tp_to_adm(t = FSST_shifted, h4_FSST_shifted)
plot(adm_FSST_4km)
title('4km FSST adm')

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

###################################9km      system tracts isolated 
h9 = data_kitten$adm18..m.

h9_LST = h9[ t >= 1.25 & t <= 1.75]                       #LST 9km
h9_LST_shifted = h9_LST - min(h9_LST)
adm_LST_9km = tp_to_adm(t = LST_shifted, h9_LST_shifted)
plot(adm_LST_9km)
title('9km LST adm')
abline(v=(0.25),col="coral",lwd=3,lty='dashed')

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

slope9km = t[t >= 0.55 & t <= 0.66]                      #isolate position of slope
slope9km_shifted = slope9km - 1.25
h9_slope = h9[t >= 1.8 & t <= 1.91]
adm_slope_9km = tp_to_adm(t = slope9km, h9_slope)
plot(adm_slope_9km)

adm_st_9km <- list(adm_LST_9km,adm_TST_9km,adm_HST_9km,adm_FSST_9km)

w_hiat_no_9km <- c()                          # Number of gaps per System Track at 9km.
for(hiat_no_9km in adm_st_9km){
  v_hiat_no_9km <- get_hiat_no(hiat_no_9km)
  w_hiat_no_9km <- c(w_hiat_no_9km,v_hiat_no_9km)
}
plot(w_hiat_no_9km,
     main = "number of hiatus per system tract in 9km adm",
     xlab = "LST,TST,HST,FSST",
     ylab = "number of hiatus",
     pch = c(16))

get_completeness(adm_LST_9km)
get_completeness(adm_TST_9km)
get_completeness(adm_HST_9km)
get_completeness(adm_FSST_9km)
get_completeness(adm_slope_9km)
sum(get_hiat_duration(adm_LST_9km))
sum(get_hiat_duration(adm_TST_9km))
sum(get_hiat_duration(adm_HST_9km))
sum(get_hiat_duration(adm_FSST_9km))

###################################11km      system tracts isolated 
h11 = data_kitten$adm22..m.

h11_LST = h11[ t >= 1.25 & t <= 1.75]                       #LST 11km
h11_LST_shifted = h11_LST - min(h11_LST)
adm_LST_11km = tp_to_adm(t = LST_shifted, h11_LST_shifted)
plot(adm_LST_11km)
title('11km LST adm')

h11_TST = h11[ t >= 1.75 & t <= 2.25]                       #TST 11km
h11_TST_shifted = h11_TST - min(h11_LST)
adm_TST_11km = tp_to_adm(t = TST_shifted, h11_TST_shifted)
plot(adm_TST_11km)
title('11km TST adm')

h11_HST = h11[ t >= 2.25 & t <= 2.75]                       #HST 11km
h11_HST_shifted = h11_HST - min(h11_LST)
adm_HST_11km = tp_to_adm(t = HST_shifted, h11_HST_shifted)
plot(adm_HST_11km)
title('11km HST adm')

h11_FSST = h11[ t >= 2.75 & t <= 3.25]                      #FSST 11km
h11_FSST_shifted = h11_FSST - min(h11_FSST)
adm_FSST_11km = tp_to_adm(t = FSST_shifted, h11_FSST_shifted)
plot(adm_FSST_11km)
title('11km FSST adm')

slope11km = t[t >= 1.69 & t <= 1.78]                            #isolate position of slope
slope11km_shifted = slope11km - 1.25
h11_slope = h11[t >= 2.94 & t <= 3.03]
adm_slope_11km = tp_to_adm(t = slope11km, h11_slope)
plot(adm_slope_11km)                                            #plot just the slope

adm_st_11km <- list(adm_LST_11km,adm_TST_11km,adm_HST_11km,adm_FSST_11km)

w_hiat_no_11km <- c()                          # Number of gaps per System Track at 11km.
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
h15 = data_kitten$adm30..m.

h15_LST = h15[ t >= 1.25 & t <= 1.75]                       #LST 15km
h15_LST_shifted = h15_LST - min(h15_LST)
adm_LST_15km = tp_to_adm(t = LST_shifted, h15_LST_shifted)
plot(adm_LST_15km);
title('15km LST adm')

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

#Incompleteness percentage per adm #################
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
(get_completeness(adm_1km)+get_completeness(adm_2km)+get_completeness(adm_3km)+get_completeness(adm_4km)+get_completeness(adm_5km))/5
#completeness average _ all states
(get_completeness(adm_7km)+get_completeness(adm_8km)+get_completeness(adm_9km)+get_completeness(adm_10km)+get_completeness(adm_11km))/5
#completeness average _ basin
(get_completeness(adm_13km)+get_completeness(adm_14km)+get_completeness(adm_15km)+get_completeness(adm_16km)+get_completeness(adm_17km)+get_completeness(adm_18km)+get_completeness(adm_19km)+get_completeness(adm_20km))/8


#Number of gaps per adm ##############################
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
(get_hiat_no(adm_1km)+get_hiat_no(adm_2km)+get_hiat_no(adm_3km)+get_hiat_no(adm_4km)+get_hiat_no(adm_5km))/5
#hiat_no average _ all states of platform
(get_hiat_no(adm_7km)+get_hiat_no(adm_8km)+get_hiat_no(adm_9km)+get_hiat_no(adm_10km)+get_hiat_no(adm_11km))/5
#hiat_no average _ basin
(get_hiat_no(adm_13km)+get_hiat_no(adm_14km)+get_hiat_no(adm_15km)+get_hiat_no(adm_16km)+get_hiat_no(adm_17km)+get_hiat_no(adm_18km)+get_hiat_no(adm_19km)+get_hiat_no(adm_20km))/8

#Average duration of hiatus per adm ##################
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
#hiat_duration average _ platform interior
(mean(get_hiat_duration(adm_1km))+mean(get_hiat_duration(adm_2km))+mean(get_hiat_duration(adm_3km))+mean(get_hiat_duration(adm_4km))+mean(get_hiat_duration(adm_5km)))/5
#hiat_duration average _ all stages of platform
(mean(get_hiat_duration(adm_6km))+mean(get_hiat_duration(adm_7km))+mean(get_hiat_duration(adm_8km))+mean(get_hiat_duration(adm_9km))+mean(get_hiat_duration(adm_10km))+mean(get_hiat_duration(adm_11km)))/5
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

#Height per adm ######################################
w_height <- c()
for(height in adm){
  v_height <- get_height(height)
  w_height <- c(w_height,v_height)
}
plot(w_height,
     main = "height",
     xlab = "?",
     ylab = "meter",
     pch = c(16))    

#Total thickness per adm ######################################
w_thickness <- c()
for(thickness in adm){
  v_thickness <- get_total_thickness(thickness)
  w_thickness <- c(w_thickness,v_thickness)
}
plot(w_thickness,
     main = "total thickness",
     xlab = "adm km",
     ylab = "meter",
     pch = c(16))   

