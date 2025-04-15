data_kitten = read.csv("data/alcap-example7_adm.csv")
water_depth = read.csv("data/alcap-example7_wd.csv")

#grid: 20km x 50km
#time: 4 Myr
#Sea level from Holland & Patzkowsky 2015

library(admtools)
library(StratPal)
library(ggplot2)
library(dplyr)

#Changes of system tracts
LST_to_TST = 0.5      # Lowstand system tract to Transgressive system tract
TST_to_HST = 1        # Transgressive system tract to Highstand system tract
HST_to_FSST = 1.5     # Highstand system tract to Falling Stage system tract

### Eustatic Sea Level curve
AMP1 = 15                           # Amplitude and period of 3rd order fluctuation
PER1 = 2                            
AMP2 = 2.5                          # Amplitude and period of 4th order fluctuation
PER2 = 0.2
t = seq(1.25, 3.25, by=0.001)       # time frame
t_mod = t[t >= 1.25 & t <= 3.25]    # Subset data
t_mod_shifted = t_mod - min(t_mod)  # shifted from 0 to 2 Myr
sl = (AMP1*sin(2*pi*t/PER1)) + (AMP2*sin(2*pi*t/PER2))    # the curve
plot(t_mod_shifted,sl,type='l',     # Plot curve
     main="Eustatic Sea Level",
     xlab="time [Myr]",
     ylab="[m]")
abline(v=0.5,col="cyan4",lwd=3,lty='dashed')    # separate different system tracts
abline(v=1,col="cyan4",lwd=3,lty='dashed')
abline(v=1.5,col="cyan4",lwd=3,lty='dashed')
text(x=0.2, y=0, "LST",cex=0.9,col='darkblue',font=3)
text(x=0.75, y=12, "TST",cex=0.9,col='darkblue',font=3)
text(x=1.25, y=0, "HST",cex=0.9,col='darkblue',font=3)
text(x=1.75, y=12, "FSST",cex=0.9,col='darkblue',font=3)

# 40 age depth models. From 500m to 20km, every 500m
t = data_kitten$time..Myr.
t_mod = t[t >= 1.25 & t <= 3.25]                           # Subset data
t_mod_shifted = t_mod - min(t_mod)

###  4 km     adm
h4 = data_kitten$adm8..m.      
plot(h4,type='l')                                           # entire runtime (4 Myr)
h4_mod = h4[t >= 1.25 & t <= 3.25]                          # modify to show 1.25-3.25 Myr
h4_mod_shifted = h4_mod - min(h4_mod)                       # shift to 0-2 Myr
adm4_mod = tp_to_adm(t = t_mod, h4_mod)                    
adm_4km = tp_to_adm(t = t_mod_shifted, h4_mod_shifted)     
plot(adm_4km,lwd_acc = 2,lwd_destr = 0)                     # plot 1.25 - 3.25 Myr                              
title(main = "age depth model 4km",
      xlab = "time [Myr]", ylab = "height [m]",
      cex.main = 1.2)
abline(v=(LST_to_TST),col="cyan4",lwd=3,lty="dashed")
abline(v=(TST_to_HST),col="cyan4",lwd=3,lty="dashed")
abline(v=(HST_to_FSST),col="cyan4",lwd=3,lty="dashed")
text(x=0.21, y=80, "LST",cex=0.8,col='darkblue',font=3)
text(x=0.75, y=80, "TST",cex=0.8,col='darkblue',font=3)
text(x=1.25, y=20, "HST",cex=0.8,col='darkblue',font=3)
text(x=1.8, y=20, "FSST",cex=0.8,col='darkblue',font=3)

str(adm_4km)

adm_4km_h <- adm_4km$h        # Extract height from adm_4km
adm_4km_t <- adm_4km$t        # and time

df_adm_4km <- data.frame(T_Myr = adm_4km_t,H_4km = adm_4km_h)   # Create a data frame

df_adm_4km <- df_adm_4km %>%                 # Calculate height difference between data points and assign line thickness
  mutate(diff_H_4km = c(NA, diff(H_4km)),
    line_width = ifelse(diff_H_4km == 0, 0.1, 2))

df_segments <- df_adm_4km %>%                # Prepare data for plotting with geom_segment
  mutate(                                    # Each row represents a segment from point i to i+1
    T_Myr_next = lead(T_Myr),
    H_4km_next = lead(H_4km)
  ) %>% filter(!is.na(T_Myr_next))           # Remove the last row (no next point)

# Line plot
ggplot(df_segments) + geom_segment(aes(x=T_Myr, y=H_4km, xend=T_Myr_next,
    yend=H_4km_next, size=line_width)) + 
  scale_size_identity() + theme_minimal() + 
  labs(x="Time [Myr]", y="Stratigraphic heigh [m]", 
       title="Age Depth Model at 4 km") + 
  geom_vline(xintercept=LST_to_TST,color="cyan4", linetype="dashed",size=1) + 
  geom_vline(xintercept=TST_to_HST,color="cyan4", linetype="dashed",size=1) + 
  geom_vline(xintercept=HST_to_FSST,color="cyan4",linetype="dashed",size=1) +
  geom_text(aes(x=0.25,y=80,label="LST"),color="darkblue",size=4) + 
  geom_text(aes(x=0.75,y=80,label="TST"),color="darkblue",size=4) +
  geom_text(aes(x=1.25,y=20,label="HST"),color="darkblue",size=4) +
  geom_text(aes(x=1.75,y=20,label="FSST"),color="darkblue",size=4)

w4 = water_depth$wd8..m.                                    # Water Depth (wd) at 4 km
plot(w4,type='l',main="complete wd 4km",                    # Entire runtime (4 Myr)
     xlab="time [Myr]",ylab="[m]")
w4_mod = w4[t >= 1.25 & t <= 3.25]                          # modify and shift
w4_mod_shifted = w4_mod - min(w4_mod)
wd4_mod = tp_to_adm(t = t_mod, w4_mod)                    
wd_4km = tp_to_adm(t = t_mod_shifted, w4_mod_shifted) 
plot(wd_4km)                                                # 1.25 - 3.25 Myr
title(main="water depth at 4km",
      ylab="[m]",
      xlab="time [Myr]")

###  6 km       adm
h6 = data_kitten$adm12..m. 
plot(h6,type='l')                                           # Entire runtime (4 Myr)
###  7 km       adm                                       
h7 = data_kitten$adm14..m.                                 
plot(h7,type='l')                                           # Entire runtime (4 Myr)
###  8 km       adm
h8 = data_kitten$adm16..m.
plot(h8,type='l')                                           # Entire runtime (4 Myr)

###  9 km       adm
h9 = data_kitten$adm18..m.    
plot(h9,type='l')                                           # Entire runtime (4 Myr)
h9_mod = h9[ t >= 1.25 & t <= 3.25]                         # modify to show 1.25-3.25 Myr
h9_mod_shifted = h9_mod - min(h9_mod)                       # shift to 0-2 Myr
adm9_mod = tp_to_adm(t = t_mod, h9_mod)                     
adm_9km = tp_to_adm(t = t_mod_shifted, h9_mod_shifted)      
plot(adm_9km,lwd_acc = 2,lwd_destr = 0)                     # plot 1.25 - 3.25 Myr
title(main = "age depth model 9km",
      xlab = "time [Myr]", ylab = "height [m]",
      cex.main = 1.2)
abline(v=(LST_to_TST),col="cyan4",lwd=3,lty="dashed")
abline(v=(TST_to_HST),col="cyan4",lwd=3,lty="dashed")
abline(v=(HST_to_FSST),col="cyan4",lwd=3,lty="dashed")
text(x=0.22, y=50, "LST",cex=0.8,col='darkblue',font=3)
text(x=0.78, y=50, "TST",cex=0.8,col='darkblue',font=3)
text(x=1.28, y=50, "HST",cex=0.8,col='darkblue',font=3)
text(x=1.8, y=50, "FSST",cex=0.8,col='darkblue',font=3)

adm_9km_h <- adm_9km$h        # Extract height from adm_9km
adm_9km_t <- adm_9km$t        # and time

df_adm_9km <- data.frame(T_Myr = adm_9km_t,H_9km = adm_9km_h)   # Create a data frame

df_adm_9km <- df_adm_9km %>%                 # Calculate height difference between data points and assign line thickness
  mutate(diff_H_9km = c(NA, diff(H_9km)),
         line_width = ifelse(diff_H_9km == 0, 0.1, 2))

df_segments <- df_adm_9km %>%                # Prepare data for plotting with geom_segment
  mutate(                                    # Each row represents a segment from point i to i+1
    T_Myr_next = lead(T_Myr),
    H_9km_next = lead(H_9km)
  ) %>% filter(!is.na(T_Myr_next))           # Remove the last row (no next point)

# Line plot
ggplot(df_segments) + geom_segment(aes(x=T_Myr, y=H_9km, xend=T_Myr_next,
                                       yend=H_9km_next, size=line_width)) + 
  scale_size_identity() + theme_minimal() + 
  labs(x="Time [Myr]", y="Stratigraphic heigh [m]", 
       title="Age Depth Model at 9 km") + 
  geom_vline(xintercept=LST_to_TST,color="cyan4", linetype="dashed",size=1) + 
  geom_vline(xintercept=TST_to_HST,color="cyan4", linetype="dashed",size=1) + 
  geom_vline(xintercept=HST_to_FSST,color="cyan4",linetype="dashed",size=1) +
  geom_text(aes(x=0.25,y=80,label="LST"),color="darkblue",size=4) + 
  geom_text(aes(x=0.75,y=80,label="TST"),color="darkblue",size=4) +
  geom_text(aes(x=1.25,y=20,label="HST"),color="darkblue",size=4) +
  geom_text(aes(x=1.75,y=20,label="FSST"),color="darkblue",size=4)

w9 = water_depth$wd18..m.                                   # Water depth (wd) at 9 km
plot(w9,type='l',main="complete wd 9km",                    # Entire runtime (4 Myr)
     xlab="time [Myr]",ylab="[m]")
w9_mod = w9[ t >= 1.25 & t <= 3.25]                         # modify and shift
w9_mod_shifted = w9_mod - min(w9_mod)
wd9_mod = tp_to_adm(t = t_mod, w9_mod)                    
wd_9km = tp_to_adm(t = t_mod_shifted, w9_mod_shifted)
plot(wd_9km)                                                # Plot 1.25 - 3.25 Myr
title(main="water depth at 9km",
      ylab="[m]",
      xlab="time [Myr]")

### 10km adm
h10 = data_kitten$adm20..m.                           
plot(h10,type='l')                                          # Entire runtime (4 Myr)
h10_mod = h10[ t >= 1.25 & t <= 3.25]
h10_mod_shifted = h10_mod - min(h10_mod)
adm10_mod = tp_to_adm(t = t_mod, h10_mod)                  
adm_10km = tp_to_adm(t = t_mod_shifted, h10_mod_shifted)   
plot(adm_10km,lwd_acc = 2,lwd_destr = 0)                    # 1.25 - 3.25 Myr

### 11km adm
h11 = data_kitten$adm22..m.      
plot(h11,type='l')                                          # Entire runtime (4 Myr)
h11_mod = h11[ t >= 1.25 & t <= 3.25]                       # Modify to show 1.25-3.25 Myr
h11_mod_shifted = h11_mod - min(h11_mod)                    # Shift to 0-2 Myr
adm11_mod = tp_to_adm(t = t_mod, h11_mod)                  
adm_11km = tp_to_adm(t = t_mod_shifted, h11_mod_shifted)
plot(adm_11km,lwd_acc = 2,lwd_destr = 0)                    # Plot 1.25 - 3.25
title('11km adm')
abline(v=(LST_to_TST),col="cyan4",lwd=3,lty='dashed')
abline(v=(TST_to_HST),col="cyan4",lwd=3,lty='dashed')
abline(v=(HST_to_FSST),col="cyan4",lwd=3,lty='dashed')
text(x=0.22, y=80, "LST",cex=0.8,col='darkblue',font=3)
text(x=0.78, y=80, "TST",cex=0.8,col='darkblue',font=3)
text(x=1.28, y=80, "HST",cex=0.8,col='darkblue',font=3)
text(x=1.9, y=80, "FSST",cex=0.8,col='darkblue',font=3)

adm_11km_h <- adm_11km$h        # Extract height from adm_11km
adm_11km_t <- adm_11km$t        # and time

df_adm_11km <- data.frame(T_Myr = adm_11km_t,H_11km = adm_11km_h)   # Create a data frame

df_adm_11km <- df_adm_11km %>%                 # Calculate height difference between data points and assign line thickness
  mutate(diff_H_11km = c(NA, diff(H_11km)),
         line_width = ifelse(diff_H_11km == 0, 0.1, 2))

df_segments <- df_adm_11km %>%                # Prepare data for plotting with geom_segment
  mutate(                                    # Each row represents a segment from point i to i+1
    T_Myr_next = lead(T_Myr),
    H_11km_next = lead(H_11km)
  ) %>% filter(!is.na(T_Myr_next))           # Remove the last row (no next point)

# Line plot
ggplot(df_segments) + geom_segment(aes(x=T_Myr, y=H_11km, xend=T_Myr_next,
                                       yend=H_11km_next, size=line_width)) + 
  scale_size_identity() + theme_minimal() + 
  labs(x="Time [Myr]", y="Stratigraphic heigh [m]", 
       title="Age Depth Model at 11 km") + 
  geom_vline(xintercept=LST_to_TST,color="cyan4", linetype="dashed",size=1) + 
  geom_vline(xintercept=TST_to_HST,color="cyan4", linetype="dashed",size=1) + 
  geom_vline(xintercept=HST_to_FSST,color="cyan4",linetype="dashed",size=1) +
  geom_text(aes(x=0.25,y=100,label="LST"),color="darkblue",size=4) + 
  geom_text(aes(x=0.75,y=100,label="TST"),color="darkblue",size=4) +
  geom_text(aes(x=1.25,y=100,label="HST"),color="darkblue",size=4) +
  geom_text(aes(x=1.95,y=100,label="FSST"),color="darkblue",size=4)

w11 = water_depth$wd22..m.                                  # Water depth (wd) at 11 km                                
plot(w11,type='l',main="complete wd 11km",                  # Entire runtime (4 Myr)
     xlab="time [Myr]",ylab="[m]")
w11_mod = w11[ t >= 1.25 & t <= 3.25]                       # Modify and shift
w11_mod_shifted = w11_mod - min(w11_mod)
wd11_mod = tp_to_adm(t = t_mod, w11_mod)                    
wd_11km = tp_to_adm(t = t_mod_shifted, w11_mod_shifted)     
plot(wd_11km)                                               # plot 1.25-3.25 Myr
title(main="water depth at 11km",
      ylab="[m]",
      xlab="time [Myr]")

###  12 km    adm
h12 = data_kitten$adm24..m. 
plot(h12,type='l')                                          # Entire runtime (4 Myr)
h12_mod = h12[ t >= 1.25 & t <= 3.25]                       # Modify to show 1.25-3.25 Myr
h12_mod_shifted = h12_mod - min(h12_mod)                    # Sift to 0-2 Myr
adm12_mod = tp_to_adm(t = t_mod, h12_mod)           
adm_12km = tp_to_adm(t = t_mod_shifted, h12_mod_shifted) 
plot(adm_12km,lwd_acc = 2,lwd_destr = 0)                    # 1.25 - 3.25 Myr
abline(v=(LST_to_TST),col="cyan4",lwd=3,lty='dashed')
abline(v=(TST_to_HST),col="cyan4",lwd=3,lty='dashed')
abline(v=(HST_to_FSST),col="cyan4",lwd=3,lty='dashed')
title(main = "age depth model 12km",
      xlab = "time [Myr]", ylab = "height [m]",
      cex.main = 1.2)
text(x=0.21, y=40, "Lowstand 
     system tract",cex=0.6,col='darkblue',font=3)
text(x=0.7, y=50, "Transgressive 
     system tract",cex=0.6,col='darkblue',font=3)
text(x=1.2, y=60, "Highstand
     system tract",cex=0.6,col='darkblue',font=3)
text(x=1.7, y=70, "Falling stage 
     system tract",cex=0.6,col='darkblue',font=3)

### 15 km    adm
h15 = data_kitten$adm30..m.   
plot(h15,type='l')                                          # Entire runtime (4 Myr)
h15_mod = h15[ t >= 1.25 & t <= 3.25]                       # Modify to show 1.25-3.25 Myr
h15_mod_shifted = h15_mod - min(h15_mod)                    # Shift to 0-2 Myr
adm15_mod = tp_to_adm(t = t_mod, h15_mod)                  
adm_15km = tp_to_adm(t = t_mod_shifted, h15_mod_shifted)   
plot(adm_15km,lwd_acc = 2,lwd_destr = 0)                    # 1.25 - 3.25 Myr
abline(v=(LST_to_TST),col="cyan4",lwd=3,lty='dashed')
abline(v=(TST_to_HST),col="cyan4",lwd=3,lty='dashed')
abline(v=(HST_to_FSST),col="cyan4",lwd=3,lty='dashed')
title(main = "age depth model 15km",
      xlab = "time [Myr]", ylab = "height [m]",
      cex.main = 1.2)
text(x=0.25, y=20, "LST",cex=0.8,col='darkblue',font=3)
text(x=0.75, y=20, "TST",cex=0.8,col='darkblue',font=3)
text(x=1.25, y=10, "HST",cex=0.8,col='darkblue',font=3)
text(x=1.75, y=10, "FSST",cex=0.8,col='darkblue',font=3)

adm_15km_h <- adm_15km$h        # Extract height from adm_15km
adm_15km_t <- adm_15km$t        # and time

df_adm_15km <- data.frame(T_Myr = adm_15km_t,H_15km = adm_15km_h)   # Create a data frame

df_adm_15km <- df_adm_15km %>%                 # Calculate height difference between data points and assign line thickness
  mutate(diff_H_15km = c(NA, diff(H_15km)),
         line_width = ifelse(diff_H_15km == 0, 0.1, 2))

df_segments <- df_adm_15km %>%                # Prepare data for plotting with geom_segment
  mutate(                                     # Each row represents a segment from point i to i+1
    T_Myr_next = lead(T_Myr),
    H_15km_next = lead(H_15km)
  ) %>% filter(!is.na(T_Myr_next))            # Remove the last row (no next point)

# Line plot
ggplot(df_segments) + geom_segment(aes(x=T_Myr, y=H_15km, xend=T_Myr_next,
                                       yend=H_15km_next, size=line_width)) + 
  scale_size_identity() + theme_minimal() + 
  labs(x="Time [Myr]", y="Stratigraphic heigh [m]", 
       title="Age Depth Model at 15 km") + 
  geom_vline(xintercept=LST_to_TST,color="cyan4", linetype="dashed",size=1) + 
  geom_vline(xintercept=TST_to_HST,color="cyan4", linetype="dashed",size=1) + 
  geom_vline(xintercept=HST_to_FSST,color="cyan4",linetype="dashed",size=1) +
  geom_text(aes(x=0.25,y=40,label="LST"),color="darkblue",size=4) + 
  geom_text(aes(x=0.75,y=40,label="TST"),color="darkblue",size=4) +
  geom_text(aes(x=1.25,y=40,label="HST"),color="darkblue",size=4) +
  geom_text(aes(x=1.95,y=40,label="FSST"),color="darkblue",size=4)

w15 = water_depth$wd30..m.                                    # Water depth (wd) at 15 km                               
plot(w15,type='l',main="complete wd 15km",                    # Entire runtime (4 Myr)
     xlab="time [Myr]")
w15_mod = w15[ t >= 1.25 & t <= 3.25]                         # Modify and shift
w15_mod_shifted = w15_mod - min(w15_mod)
wd15_mod = tp_to_adm(t = t_mod, w15_mod)                    
wd_15km = tp_to_adm(t = t_mod_shifted, w15_mod_shifted)     
plot(wd_15km)                                                 # Plot 1.25 - 3.25
title(main="water depth at 15 km",
      ylab="[m]",
      xlab="time [Myr]")


###############################################################################
### 4km       constant extinction rate
ext0.4 <- p3(rate = 500, from = min_time(adm_4km), to = max_time(adm_4km)) |>    # constant rate in time domain
  time_to_strat(adm_4km, destructive = FALSE)                       # transform into depth domain
  #hist(main = "Last occurrence at 4 km with constant extinction rate",
  #    ylab = "Last occurrence per meter",
  #    xlab = "Stratigraphic height [m]",
  #    breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100),
  #    col="white")
#abline(v=(time_to_strat(LST_to_TST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
#abline(v=(time_to_strat(TST_to_HST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
#abline(v=(time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
#text(x=8, y=95, "LST",cex=0.9,col='darkblue',font=3)
#text(x=45, y=95, "TST",cex=0.9,col='darkblue',font=3)
#text(x=82, y=95, "HST",cex=0.9,col='darkblue',font=3)
#text(x=97, y=95, "FSST",cex=0.9,col='darkblue',font=3)
qplot(ext0.4, geom="histogram",binwidth=1) + coord_flip() + geom_vline(
  xintercept=time_to_strat(LST_to_TST,adm_4km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
    xintercept=time_to_strat(TST_to_HST,adm_4km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) +
  labs(title = "4 km with constant extinction rate (sc.1)",x = "stratigraphic height [m]",y = "Last occurrence") +
  geom_text(aes(x=10,y=110,label="LST"),color="deepskyblue4",size=4) + 
  geom_text(aes(x=50,y=110,label="TST"),color="deepskyblue4",size=4) +
  geom_text(aes(x=85,y=110,label="HST"),color="deepskyblue4",size=4) +
  geom_text(aes(x=102,y=110,label="FSST"),color="deepskyblue4",size=4) + 
  theme_minimal()


### 9km        constant extinction rate
ext0.9 <- p3(rate = 500, from = min_time(adm_9km), to = max_time(adm_9km)) |>    # constant rate in time domain
  time_to_strat(adm_9km, destructive = FALSE)                          # transform into depth domain
  #hist(main = "Last occurrence at 9 km with constant extinction rate",
  #   ylab = "Last occurrence per meter",
  #    xlab = "Stratigraphic height [m]",
  #    breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100),
  #    col="white")
#abline(v=(time_to_strat(LST_to_TST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
#abline(v=(time_to_strat(TST_to_HST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
#abline(v=(time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
#text(x=4, y=165, "LST",cex=0.9,col='darkblue',font=3)
#text(x=67, y=165, "TST",cex=0.9,col='darkblue',font=3)
#text(x=138, y=165, "HST",cex=0.9,col='darkblue',font=3)
#text(x=157, y=165, "FSST",cex=0.9,col='darkblue',font=3)
qplot(ext0.9, geom="histogram",binwidth=1) + coord_flip() + geom_vline(
  xintercept=time_to_strat(LST_to_TST,adm_9km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
    xintercept=time_to_strat(TST_to_HST,adm_9km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) +
  labs(title = "9 km with constant extinction rate (sc.1)",x = "stratigraphic height [m]",y = "Last occurrence") +
  geom_text(aes(x=0,y=100,label="LST"),color="deepskyblue4",size=4) + 
  geom_text(aes(x=70,y=100,label="TST"),color="deepskyblue4",size=4) +
  geom_text(aes(x=140,y=100,label="HST"),color="deepskyblue4",size=4) +
  geom_text(aes(x=160,y=100,label="FSST"),color="deepskyblue4",size=4) + 
  theme_minimal()

### 11km      constant extinction rate
ext0.11 <- p3(rate = 500, from = min_time(adm_11km), to = max_time(adm_11km)) |>  # constant rate in time domain
  time_to_strat(adm_11km, destructive = FALSE)                       # transform into depth domain
  #hist(main = "Last occurrence at 11 km with constant extinction rate ",
  #    ylab = "Last occurrence per meter",
  #    xlab = "Stratigraphic height [m]",
  #    breaks = seq(from = min_height(adm_11km), to = max_height(adm_11km), length.out = 100),
  #    col="white")
#abline(v=(time_to_strat(LST_to_TST,adm_11km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
#abline(v=(time_to_strat(TST_to_HST,adm_11km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
#abline(v=(time_to_strat(HST_to_FSST,adm_11km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
#text(x=3, y=80, "LST",cex=0.9,col='darkblue',font=3)
#text(x=13, y=80, "TST",cex=0.9,col='darkblue',font=3)
#text(x=23, y=80, "HST",cex=0.9,col='darkblue',font=3)
#text(x=50, y=80, "FSST",cex=0.9,col='darkblue',font=3)
qplot(ext0.11, geom="histogram",binwidth=1) + coord_flip() + geom_vline(
  xintercept=time_to_strat(LST_to_TST,adm_11km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
    xintercept=time_to_strat(TST_to_HST,adm_11km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(HST_to_FSST,adm_11km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) +
  labs(title = "11 km with constant extinction rate (sc.1)",x = "stratigraphic height [m]",y = "Last occurrence") +
  geom_text(aes(x=0,y=90,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=13,y=90,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=23,y=90,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=100,y=90,label="FSST"),color="deepskyblue3",size=4)

### 12km     constant extinction rate
p3(rate = 500, from = min_time(adm_12km), to = max_time(adm_12km)) |>  # constant rate in time domain
  time_to_strat(adm_12km, destructive = FALSE)                         # transform into depth domain
  #hist(main = "Last occurrence at 12 km with constant extinction rate ",
  #    ylab = "Last occurrence per meter",
  #    xlab = "Stratigraphic height [m]",
  #    breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100),
  #    col="white")
#abline(v=(time_to_strat(LST_to_TST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
#abline(v=(time_to_strat(TST_to_HST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
#abline(v=(time_to_strat(HST_to_FSST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
#text(x=5, y=30, "LST",cex=0.9,col='darkblue',font=3)
#text(x=15, y=30, "TST",cex=0.9,col='darkblue',font=3)
#text(x=24, y=30, "HST",cex=0.9,col='darkblue',font=3)
#text(x=32, y=30, "FSST",cex=0.9,col='darkblue',font=3)
#(time_to_strat(HST_to_FSST,adm_12km,destructive=FALSE))

### 15km      constant extinction rate
ext0.15 <- p3(rate = 500, from = min_time(adm_15km), to = max_time(adm_15km)) |>  # constant rate in time domain
  time_to_strat(adm_15km, destructive = FALSE)                       # transform into depth domain
  #hist(main = "Last occurrence at 15 km with constant extinction rate",
  #    ylab = "Last occurrence per meter",
  #    xlab = "Stratigraphic height [m]",
  #    breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100),
  #    col="white")
#abline(v=(time_to_strat(LST_to_TST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
#abline(v=(time_to_strat(TST_to_HST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
#abline(v=(time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
#text(x=4, y=45, "LST",cex=0.9,col='darkblue',font=3)
#text(x=13, y=45, "TST",cex=0.9,col='darkblue',font=3)
#text(x=22, y=45, "HST",cex=0.9,col='darkblue',font=3)
#text(x=30, y=45, "FSST",cex=0.9,col='darkblue',font=3)
qplot(ext0.15, geom="histogram",binwidth=0.5) + coord_flip() + geom_vline(
  xintercept=time_to_strat(LST_to_TST,adm_15km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
    xintercept=time_to_strat(TST_to_HST,adm_15km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) +
  labs(title = "15 km with constant extinction rate (sc.1)",x = "stratigraphic height [m]",y = "Last occurrence") +
  geom_text(aes(x=4,y=90,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=13,y=90,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=21,y=90,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=32,y=90,label="FSST"),color="deepskyblue3",size=4)

###
max_height(adm_4km)
max_height(adm_9km)
max_height(adm_11km)
max_height(adm_15km)

################################################################################
### Extinction event simulations
# Shape of extinction
tp = c(0,0.5,0.75,1,2)                            # points in time
ext_rate = c(1,1,25,1,1)                          # rate of extinction
plot(tp,ext_rate,type='l',lwd=3,                  # graph simply showing the extinction
     main="extinction event",
     xlab="time [Myr]",
     ylab="extinction rate [species/m]")

### Lowstand Systems Tract extincion - time: 1.25-1.75 Myr (shifted to 0-0.5 Myr)
#4km LST
h1.4 <- p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |> 
  time_to_strat(adm_4km, destructive = FALSE)                    # transform into depth domain
  #hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100))

ext1.4 <- cut(h1.4$breaks, c(-Inf,(time_to_strat(0.5,adm_4km,destructive=FALSE)),Inf))  # distinguishing the timing of extinction event
plot(h1.4,main = "Last occurrence of taxa at 4 km",                  # histogram with coloring
     sub = "Extinction during Lowstand Systems Tract ",
     font.sub=3,
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("red","white")[ext1.4])                          # extinction is red, rest is white 
  abline(v=(time_to_strat(0.25,adm_4km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=7, y=160, "LST",cex=0.9,col='darkblue',font=3)
text(x=45, y=160, "TST",cex=0.9,col='darkblue',font=3)
text(x=83, y=160, "HST",cex=0.9,col='darkblue',font=3)
text(x=98, y=160, "FSST",cex=0.9,col='darkblue',font=3)


#9km LST
h1.9 <- p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |> 
  time_to_strat(adm_9km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100))

ext1.9 <- cut(h1.9$breaks, c(-Inf,(time_to_strat(0.5,adm_9km,destructive=FALSE)),Inf))  # distinguishing the timing of extinction event
plot(h1.9,main = "Last occurrence of taxa at 9 km",                  # histogram with coloring
     sub = "Mass extinction during the Lowstand System Tract",
     font.sub=3,
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("red","white")[ext1.9])
  abline(v=(time_to_strat(0.25,adm_9km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=5, y=100, "LST",cex=0.9,col='darkblue',font=3)
text(x=75, y=100, "TST",cex=0.9,col='darkblue',font=3)
text(x=138, y=100, "HST",cex=0.9,col='darkblue',font=3)
text(x=155, y=100, "FSST",cex=0.9,col='darkblue',font=3)
time_to_strat(0.5,adm_4km,destructive=FALSE)

#11km LST
h1.11 <- p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_11km, destructive = FALSE) |>                    # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences  
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_11km), to = max_height(adm_11km), length.out = 100))

ext1.11 <- cut(h1.11$breaks, c(-Inf,(time_to_strat(0.5,adm_11km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h1.11,main = "Last occurrence of taxa at 11 km",                # histogram with coloring
     sub = "Extinction during the Lowstand Systems Tract",
     font.sub=3,
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("red","white")[ext1.11])
abline(v=(time_to_strat(0.25,adm_11km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
abline(v=(time_to_strat(LST_to_TST,adm_11km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')  
abline(v=(time_to_strat(TST_to_HST,adm_11km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
abline(v=(time_to_strat(HST_to_FSST,adm_11km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=3, y=150, "LST",cex=0.9,col='darkblue',font=3)
text(x=13.5, y=160, "TST",cex=0.9,col='darkblue',font=3)
text(x=21.5, y=150, "HST",cex=0.9,col='darkblue',font=3)
text(x=50, y=150, "FSST",cex=0.9,col='darkblue',font=3)

#12km LST
h1.12 <- p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_12km, destructive = FALSE) |>                    # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences  
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100))

ext1.12 <- cut(h1.12$breaks, c(-Inf,(time_to_strat(0.5,adm_12km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h1.12,main = "Last occurrence of taxa at 12 km",                # histogram with coloring
     sub = "Extinction during the Lowstand Systems Tract",
     font.sub=3,
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("red","white")[ext1.12])
  abline(v=(time_to_strat(0.25,adm_12km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')  
  abline(v=(time_to_strat(TST_to_HST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=3.8, y=70, "LST",cex=0.9,col='darkblue',font=3)
text(x=13.5, y=70, "TST",cex=0.9,col='darkblue',font=3)
text(x=21.5, y=70, "HST",cex=0.9,col='darkblue',font=3)
text(x=50, y=70, "FSST",cex=0.9,col='darkblue',font=3)

#15km LST
h1.15 <- p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |> 
  time_to_strat(adm_15km, destructive = FALSE) |>                    # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences  
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100))

ext1.15 <- cut(h1.15$breaks, c(-Inf,(time_to_strat(0.5,adm_15km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h1.15,main= "Last occurrence of taxa at 15 km",                 # histogram with coloring
     sub = "Extinction during the Lowstand Systems Tract",
     font.sub=3,
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("red","white")[ext1.15])
  abline(v=(time_to_strat(0.25,adm_15km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=3, y=60, "LST",cex=0.9,col='darkblue',font=3)
text(x=13, y=60, "TST",cex=0.9,col='darkblue',font=3)
text(x=21, y=60, "HST",cex=0.9,col='darkblue',font=3)
text(x=30, y=60, "FSST",cex=0.9,col='darkblue',font=3)

### Transgressive Systems Tract extinction - time: 1.75-2.25 Myr (shifted to 0.5-1 Myr)
#4km TST
h2.4 <- p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_4km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100))

ext2.4 <- cut(h2.4$breaks, c(-Inf,(time_to_strat(0.5,adm_4km,destructive=FALSE)),(time_to_strat(1,adm_4km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h2.4,main= "Last occurrence of taxa at 4 km",                   # histogram with coloring
     sub = "Extinction during the Transgressive Systems Tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext2.4])
  abline(v=(time_to_strat(0.75,adm_4km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=10, y=45, "LST",cex=0.9,col='darkblue',font=3)
text(x=45, y=45, "TST",cex=0.9,col='darkblue',font=3)
text(x=82, y=45, "HST",cex=0.9,col='darkblue',font=3)
text(x=100, y=45, "FSST",cex=0.9,col='darkblue',font=3)

#9km TST
h2.9 <- p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_9km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100))

ext2.9 <- cut(h2.9$breaks, c(-Inf,(time_to_strat(0.5,adm_9km,destructive=FALSE)),(time_to_strat(1,adm_9km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h2.9,main= "Last occurrence of taxa at 9 km",                   # histogram with coloring
     sub = "Extinction during the Transgressive Systems Tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext2.9])
  abline(v=(time_to_strat(0.75,adm_9km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate 
  abline(v=(time_to_strat(LST_to_TST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=5, y=58, "LST",cex=0.9,col='darkblue',font=3)
text(x=60, y=58, "TST",cex=0.9,col='darkblue',font=3)
text(x=140, y=58, "HST",cex=0.9,col='darkblue',font=3)
text(x=153, y=50, "FSST",cex=0.9,col='darkblue',font=3)
abline(v=(111))
strat_to_time(111,adm_9km)
#11km TST
h2.11 <- p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_11km, destructive = FALSE) |>                    # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_11km), to = max_height(adm_11km), length.out = 100))

ext2.11 <- cut(h2.11$breaks, c(-Inf,(time_to_strat(0.5,adm_11km,destructive=FALSE)),(time_to_strat(1,adm_11km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h2.11,main= "Last occurrence of taxa at 11 km",                 # histogram with coloring
     sub = "Extinction during the Transgressive Systems Tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext2.11])
abline(v=(time_to_strat(0.75,adm_11km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
abline(v=(time_to_strat(LST_to_TST,adm_11km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
abline(v=(time_to_strat(TST_to_HST,adm_11km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
abline(v=(time_to_strat(HST_to_FSST,adm_11km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=4, y=160, "LST",cex=0.9,col='darkblue',font=3)
text(x=13.5, y=140, "TST",cex=0.9,col='darkblue',font=3)
text(x=21.5, y=160, "HST",cex=0.9,col='darkblue',font=3)
text(x=50, y=160, "FSST",cex=0.9,col='darkblue',font=3)

#12km TST
h2.12 <- p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_12km, destructive = FALSE) |>                    # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100))

ext2.12 <- cut(h2.12$breaks, c(-Inf,(time_to_strat(0.5,adm_12km,destructive=FALSE)),(time_to_strat(1,adm_12km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h2.12,main= "Last occurrence of taxa at 12 km",                 # histogram with coloring
     sub = "Extinction during the Transgressive Systems Tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext2.12])
  abline(v=(time_to_strat(0.75,adm_12km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=4, y=190, "LST",cex=0.9,col='darkblue',font=3)
text(x=13.5, y=190, "TST",cex=0.9,col='darkblue',font=3)
text(x=21.5, y=190, "HST",cex=0.9,col='darkblue',font=3)
text(x=50, y=190, "FSST",cex=0.9,col='darkblue',font=3)

#15km TST
h2.15 <- p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_15km, destructive = FALSE) |>                    # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100))

ext2.15 <- cut(h2.15$breaks, c(-Inf,(time_to_strat(0.5,adm_15km,destructive=FALSE)),(time_to_strat(1,adm_15km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h2.15,main= "Last occurrence of taxa at 15 km",                 # histogram with coloring
     sub = "Extinction during the Transgressive Systems Tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext2.15])
  abline(v=(time_to_strat(0.75,adm_15km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=4, y=120, "LST",cex=0.9,col='darkblue',font=3)
text(x=14, y=120, "TST",cex=0.9,col='darkblue',font=3)
text(x=21.5, y=120, "HST",cex=0.9,col='darkblue',font=3)
text(x=30, y=120, "FSST",cex=0.9,col='darkblue',font=3)

#### Highstand System Tract extinction - time: 2.25-2.75 Myr (shifted to 1-1.5 Myr)
#4km HST
h3.4 <- p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_4km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100))
time_to_strat(1,adm_4km)
time_to_strat(1.25,adm_4km)
time_to_strat(1.5,adm_4km,destructive=FALSE)
ext3.4 <- cut(h3.4$breaks, c(-Inf,(time_to_strat(1,adm_4km,destructive=FALSE)),(time_to_strat(1.5,adm_4km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h3.4,main= "Last occurrence of taxa at 4 km",                   # histogram with coloring
     sub = "Extinction during the Highstand Systems Tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext3.4])
  abline(v=(time_to_strat(1.25,adm_4km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=10, y=145, "LST",cex=0.9,col='darkblue',font=3)
text(x=45, y=145, "TST",cex=0.9,col='darkblue',font=3)
text(x=81, y=145, "HST",cex=0.9,col='darkblue',font=3)
text(x=100, y=145, "FSST",cex=0.9,col='darkblue',font=3)

#9km HST
h3.9 <- p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_9km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100))

ext3.9 <- cut(h3.9$breaks, c(-Inf,(time_to_strat(1,adm_9km,destructive=FALSE)),(time_to_strat(1.5,adm_9km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h3.9,main= "Last occurrence of taxa at 9 km",                  # histogram with coloring
     sub = "Extinction during the Highstand Systems Tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext3.9])
  abline(v=(time_to_strat(1.25,adm_9km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=5, y=165, "LST",cex=0.9,col='darkblue',font=3)
text(x=75, y=165, "TST",cex=0.9,col='darkblue',font=3)
text(x=135, y=165, "HST",cex=0.9,col='darkblue',font=3)
text(x=154, y=165, "FSST",cex=0.9,col='darkblue',font=3)

#11km HST
h3.11 <- p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_11km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_11km), to = max_height(adm_11km), length.out = 100))

ext3.11 <- cut(h3.11$breaks, c(-Inf,(time_to_strat(1,adm_11km,destructive=FALSE)),(time_to_strat(1.5,adm_11km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h3.11,main= "Last occurrence of taxa at 11 km",                  # histogram with coloring
     sub = "Extinction during the Highstand Systems Tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext3.11])
abline(v=(time_to_strat(1.25,adm_11km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
abline(v=(time_to_strat(LST_to_TST,adm_11km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
abline(v=(time_to_strat(TST_to_HST,adm_11km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
abline(v=(time_to_strat(HST_to_FSST,adm_11km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=3, y=150, "LST",cex=0.9,col='darkblue',font=3)
text(x=12, y=130, "TST",cex=0.9,col='darkblue',font=3)
text(x=23, y=150, "HST",cex=0.9,col='darkblue',font=3)
text(x=70, y=150, "FSST",cex=0.9,col='darkblue',font=3)

#12km HST
h3.12 <- p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_12km, destructive = FALSE) |>                    # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100))

ext3.12 <- cut(h3.12$breaks, c(-Inf,(time_to_strat(1,adm_12km,destructive=FALSE)),(time_to_strat(1.5,adm_12km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h3.12,main= "Last occurrence of taxa at 12 km",                 # histogram with coloring
     sub = "Mass extinction during the Highstand system tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext3.12])
  abline(v=(time_to_strat(1.25,adm_12km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=4, y=61, "LST",cex=0.9,col='darkblue',font=3)
text(x=14, y=61, "TST",cex=0.9,col='darkblue',font=3)
text(x=22, y=61, "HST",cex=0.9,col='darkblue',font=3)
text(x=55, y=61, "FSST",cex=0.9,col='darkblue',font=3)

#15km HST
h3.15 <- p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_15km, destructive = FALSE) |>                    # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100))

ext3.15 <- cut(h3.15$breaks, c(-Inf,(time_to_strat(1,adm_15km,destructive=FALSE)),(time_to_strat(1.5,adm_15km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h3.15,main= "Last occurrence of taxa at 15 km",                 # histogram with coloring
     sub = "Extinction during the Highstand Systems Tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext3.15])
  abline(v=(time_to_strat(1.25,adm_15km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=4, y=60, "LST",cex=0.9,col='darkblue',font=3)
text(x=13, y=60, "TST",cex=0.9,col='darkblue',font=3)
text(x=22.5, y=60, "HST",cex=0.9,col='darkblue',font=3)
text(x=30, y=60, "FSST",cex=0.9,col='darkblue',font=3)

### Falling Stage System Tract extinction - time: 2.75-3.25 Myr (shifted to 1.5-2 Myr)
#4km FSST
h4.4 <- p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_4km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100))

ext4.4 <- cut(h4.4$breaks, c(-Inf,(time_to_strat(1.5,adm_4km,destructive=FALSE)),(time_to_strat(2,adm_4km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h4.4,main= "Last occurrence of taxa at 4 km",                   # histogram with coloring
     sub = "Extinction during the Falling Stage Systems Tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext4.4])
  abline(v=(time_to_strat(1.75,adm_4km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=10, y=220, "LST",cex=0.9,col='darkblue',font=3)
text(x=45, y=220, "TST",cex=0.9,col='darkblue',font=3)
text(x=82, y=220, "HST",cex=0.9,col='darkblue',font=3)
text(x=98, y=220, "FSST",cex=0.9,col='darkblue',font=3)

#9km FSST
h4.9 <- p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_9km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100))

ext4.9 <- cut(h4.9$breaks, c(-Inf,(time_to_strat(1.5,adm_9km,destructive=FALSE)),(time_to_strat(2,adm_9km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h4.9,main= "Last occurrence of taxa at 9 km",                   # histogram with coloring
     sub = "Extinction during the Falling Stage Systems Tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext4.9])
  abline(v=(time_to_strat(1.75,adm_9km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=5, y=350, "LST",cex=0.9,col='darkblue',font=3)
text(x=70, y=350, "TST",cex=0.9,col='darkblue',font=3)
text(x=138, y=350, "HST",cex=0.9,col='darkblue',font=3)
text(x=155, y=350, "FSST",cex=0.9,col='darkblue',font=3)

#11km FSST
h4.11 <- p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_11km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_11km), to = max_height(adm_11km), length.out = 100))

ext4.11 <- cut(h4.11$breaks, c(-Inf,(time_to_strat(1.5,adm_11km,destructive=FALSE)),(time_to_strat(2,adm_11km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h4.11,main= "Last occurrence of taxa at 11 km",                   # histogram with coloring
     sub = "Extinction during the Falling Stage Systems Tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext4.11])
abline(v=(time_to_strat(1.75,adm_11km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
abline(v=(time_to_strat(LST_to_TST,adm_11km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
abline(v=(time_to_strat(TST_to_HST,adm_11km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
abline(v=(time_to_strat(HST_to_FSST,adm_11km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=3, y=100, "LST",cex=0.9,col='darkblue',font=3)
text(x=12, y=90, "TST",cex=0.9,col='darkblue',font=3)
text(x=22, y=100, "HST",cex=0.9,col='darkblue',font=3)
text(x=80, y=100, "FSST",cex=0.9,col='darkblue',font=3)

#12km FSST
h4.12 <- p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_12km, destructive = FALSE) |>                    # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100))

ext4.12 <- cut(h4.12$breaks, c(-Inf,(time_to_strat(1.5,adm_12km,destructive=FALSE)),(time_to_strat(2,adm_12km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h4.12,main= "Last occurrence of taxa at 12 km",                 # histogram with coloring
     sub = "Mass extinction during the Falling Stage system tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext4.12])
  abline(v=(time_to_strat(1.75,adm_12km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=4, y=52, "LST",cex=0.9,col='darkblue',font=3)
text(x=13, y=52, "TST",cex=0.9,col='darkblue',font=3)
text(x=22, y=52, "HST",cex=0.9,col='darkblue',font=3)
text(x=50, y=52, "FSST",cex=0.9,col='darkblue',font=3)

#15km FSST
h4.15 <- p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_15km, destructive = FALSE) |>                    # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100))

ext4.15 <- cut(h4.15$breaks, c(-Inf,(time_to_strat(1.5,adm_15km,destructive=FALSE)),(time_to_strat(2,adm_15km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h4.15,main= "Last occurrence of taxa at 15 km",                 # histogram with coloring
     sub = "Extinction during the Falling Stage Systems Tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext4.15])
  abline(v=(time_to_strat(1.75,adm_15km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=4, y=70, "LST",cex=0.9,col='darkblue',font=3)
text(x=13, y=70, "TST",cex=0.9,col='darkblue',font=3)
text(x=21.5, y=70, "HST",cex=0.9,col='darkblue',font=3)
text(x=29, y=70, "FSST",cex=0.9,col='darkblue',font=3)
abline(v=31.85)
strat_to_time(31.85,adm_15km)
time_to_strat(1.75,adm_15km,destructive=FALSE)
######################################################################

