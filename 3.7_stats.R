data_kitten = read.csv("C:/Users/sidne/OneDrive/Documenten/AA Utrecht/Guided research/Git repositories/Stratigraphic Paleobiology for Phenotypic Evolution/CarboKitten.jl/data/output/alcap-example2_adm.csv")
#grid: 20km x 50km
#time: 4 Myr
#Sea level from Holland & Patzkowsky 2015

library(admtools)
library(StratPal)

#ADM's: 500m:500m:15km
# Computer name "3.6"
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

t = data_kitten$time..Myr.
######################################### 4km         system tract isolated
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
abline(v=(0.75),col="coral",lwd=2,lty='dashed')
abline(v=(0.76),col="coral",lwd=2,lty='dashed')

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
sum(get_hiat_duration(adm_LST_9km))
sum(get_hiat_duration(adm_TST_9km))
sum(get_hiat_duration(adm_HST_9km))
sum(get_hiat_duration(adm_FSST_9km))

###################################12km      system tracts isolated 
h12 = data_kitten$adm24..m.

h12_LST = h12[ t >= 1.25 & t <= 1.75]                       #LST 12km
h12_LST_shifted = h12_LST - min(h12_LST)
adm_LST_12km = tp_to_adm(t = LST_shifted, h12_LST_shifted)
plot(adm_LST_12km)
title('12km LST adm')

h12_TST = h12[ t >= 1.75 & t <= 2.25]                       #TST 12km
h12_TST_shifted = h12_TST - min(h12_LST)
adm_TST_12km = tp_to_adm(t = TST_shifted, h12_TST_shifted)
plot(adm_TST_12km)
title('12km TST adm')

h12_HST = h12[ t >= 2.25 & t <= 2.75]                       #HST 12km
h12_HST_shifted = h12_HST - min(h12_LST)
adm_HST_12km = tp_to_adm(t = HST_shifted, h12_HST_shifted)
plot(adm_HST_12km)
title('12km HST adm')

h12_FSST = h12[ t >= 2.75 & t <= 3.25]                      #FSST 12km
h12_FSST_shifted = h12_FSST - min(h12_LST)
adm_FSST_12km = tp_to_adm(t = FSST_shifted, h12_FSST_shifted)
plot(adm_FSST_12km)
title('12km FSST adm')

adm_st_12km <- list(adm_LST_12km,adm_TST_12km,adm_HST_12km,adm_FSST_12km)

w_hiat_no_12km <- c()                          # Number of gaps per System Track at 12km.
for(hiat_no_12km in adm_st_12km){
  v_hiat_no_12km <- get_hiat_no(hiat_no_12km)
  w_hiat_no_12km <- c(w_hiat_no_12km,v_hiat_no_12km)
}
plot(w_hiat_no_12km,
     main = "number of hiatus per system tract in 12km adm",
     xlab = "LST,TST,HST,FSST",
     ylab = "number of hiatus",
     pch = c(16))

get_completeness(adm_LST_12km)
get_completeness(adm_TST_12km)
get_completeness(adm_HST_12km)
get_completeness(adm_FSST_12km)
sum(get_hiat_duration(adm_LST_12km))
sum(get_hiat_duration(adm_TST_12km))
sum(get_hiat_duration(adm_HST_12km))
sum(get_hiat_duration(adm_FSST_12km))
get_hiat_no(adm_LST_12km)
get_hiat_no(adm_TST_12km)
get_hiat_no(adm_HST_12km)
get_hiat_no(adm_FSST_12km)

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
#plot(get_height(adm_15km))

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

