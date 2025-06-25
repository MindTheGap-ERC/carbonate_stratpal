library(StratPal) #bio and ecology
library(admtools) #stratigraphy

#task is to create 100 niches described by water depth from CarboKitten model

#Created using AI sources and the help of dr. Emilia Jarochowska, Niklas Hohmann, and dr. Charlotte Summers 

Anna_wd <- read.csv("data/example_Anna_wd.csv") 
Anna_adm <- read.csv("data/example_Anna_adm.csv")

#parameters
set.seed(124)              
n_niches <- 100              # Number of niches
depth_range <- seq(min(Anna_wd), max(Anna_wd), length.out = 100)  # Water depths from -4.10 to 42.44 meters, min and max across all 4 locations 

# Generate optima more densely toward shallow depth
# Step 1: Create a uniform sequence between 0 and 1
uniform_seq <- seq(0, 1, length.out = n_niches)
# Step 2: Apply a transformation to bias toward 0 (shallow)
# There are optima at wd 0 since those would represent intertidal species, who are quite important when it comes to marine coastal environments (Andrades et al. 2019) and (Bradley et al. 2020)
optima <- max(Anna_wd) * (uniform_seq)^2  # Quadratic: more values near minimum depth (bias towards 0) # -4.10 + (42.44 + 4.10)
#for bias toward shallow: 'log1p(uniform_seq) / log1p(1)'

# Define niche width as a function of depth (e.g., linearly increasing)
min_width <- 2   # Narrowest at shallowest depth
max_width <- 10  # Widest at deepest depth
niche_widths <- min_width + (optima / max(optima)) * (max_width - min_width) #skewed towards 0 (see figure 2 from Reijmer, 2021)

# You can change this to a nonlinear function if you'd like (e.g., log, exp, etc.):
#niche_widths <- min_width + log1p(optima) / log1p(max(optima)) * (max_width - min_width)

# Create list
niches <- list()
for (i in 1:n_niches) {
  niches[[i]] <- StratPal::snd_niche(opt = optima[i],
                                     tol = niche_widths[i],
                                     cutoff_val = 0)
}

niche_val = list(niches)
for (i in 1:100) {
  niche = niches[[i]]
  niche_val[[i]]=niche(depth_range)
} 

niche_mat <- do.call(cbind, niche_val) #converting to matrix (used AI for this line)

#plot niches from list
#used AI to find out how to plot 
matplot(depth_range, niche_mat, type = "l",
        lty = 1, col = rgb(0,0,1,0.2),
        xlab = "Water depth [m]",
        ylab = "Niche value",
        main = "Niche optima across different depths")

#(using AI) plotting them as depth ranges to show more clearly the range and abundance differences across the depths
#note: the ranges continue after 42 for the highest numbered niches, because their optima fall close to or around 42, which means although it does not show they still have the broadest suitable depth range 
threshold <- 0.005  # Set your threshold; chose 0.005 since that means there is only a 0.5% likelihood of finding the organism at that depth

n_niches <- length(niches)
min_depths <- numeric(n_niches)
max_depths <- numeric(n_niches)

for (i in 1:n_niches) {
  values <- niches[[i]](depth_range)
  # Use a threshold to define "present" (e.g., >0 or >0.01)
  present <- which(values > threshold)
  if (length(present) > 0) {
    min_depths[i] <- depth_range[min(present)]
    max_depths[i] <- depth_range[max(present)]
  } else {
    min_depths[i] <- NA
    max_depths[i] <- NA
  }
}

plot(NA, xlim = c(min(Anna_wd), max(Anna_wd)), ylim = c(1, n_niches),
     xlab = "Water depth [m]", ylab = "Niche",
     main = "Depth range of each niche", yaxt = "n")
axis(2, at = 1:n_niches, labels = 1:n_niches, las = 2, cex.axis = 0.5)

for (i in 1:n_niches) {
  if (!is.na(min_depths[i]) && !is.na(max_depths[i])) {
    segments(min_depths[i], i, max_depths[i], i, lwd = 2, col = "blue")
  }
}

# Plot histogram of niche optima
hist(optima, breaks = 20, col = "blue", border = "white",
     main = "Distribution of Niche Optima Across Depth",
     xlab = "Water Depth (m)", ylab = "Number of Niches")

#adapted from:
vignette("event_data")

#change numbers to 1, 2, 3, or 4 depending on which location
#Niche modeling
#how gradient changes with time (in this case water depth)
t = Anna_wd$time..Myr.           # time steps of the model
wd = Anna_wd$wd1..m.   # water depth 2 km offshore at model time steps
gc = approxfun(t, wd)         # define function that defines how the gradient changes with time (gc = *G*radient *C*hange)

plot(gc(t), t, 
     type = "o", 
     xlab = "Time [Myr]",
     ylab = "Water depth [m]",
     main = "Water depth in section 1.5 km from shore")

#niche model in time dimension

set.seed(124)

niches_applied_t <- list()
for (i in 1:100) {
  niches_applied_t[[i]] <- p3(rate = 300, from = min(t), to = max(t)) |> 
    apply_niche(niche_def = niches[[i]], gc = gc)
}

# used AI to combine all fossil times into a data frame, including empty niches
niches_fossils_t <- do.call(rbind, lapply(seq_along(niches_applied_t), function(i) {
  x <- as.numeric(niches_applied_t[[i]])
  if (length(x) > 0 && any(!is.na(x))) {
    data.frame(t = x, niche = as.factor(i))
  } else {
    # Insert a row with NA to ensure this niche appears in the plot/legend
    data.frame(t = NA, niche = as.factor(i))
  }
}))

#new plot ensuring all niches are represented in the legend so that their numbering is meaningful 

library(ggplot2)

ggplot(niches_fossils_t, aes(x = t, fill = niche, color = niche)) +
  geom_density(alpha = 0.12, color = "black") +
  scale_fill_discrete(drop = FALSE) +
  scale_color_discrete(drop = FALSE) +
  labs(
    title = "Probability for fossil presence 1.5 km from shore",
    x = "Time [Myr]",
    y = "Density probability"
  )

#niche model in terms of stratigraphy 

strat_adm = tp_to_adm(Anna_adm$time..Myr., Anna_adm$adm1..m.)

set.seed(124)

niches_applied_strat <- list()
for (i in 1:100) {
  niches_applied_strat[[i]] <- p3(rate = 300, from = min(t), to = max(t)) |> 
    apply_niche(niche_def = niches[[i]], gc = gc) |>                    
    time_to_strat(strat_adm, destructive = TRUE, out_dom_val_h = "default")# transform into strat. domain, destroy fossils that coincide with hiatuses 
}

# used AI to combine all fossil times into a data frame, including empty niches
niches_fossils_strat <- do.call(rbind, lapply(seq_along(niches_applied_strat), function(i) {
  x <- as.numeric(niches_applied_strat[[i]])
  if (length(x) > 0 && any(!is.na(x))) {
    data.frame(t = x, niche = as.factor(i))
  } else {
    # Insert a row with NA to ensure this niche appears in the plot/legend
    data.frame(t = NA, niche = as.factor(i))
  }
}))

#new plot ensuring all niches are represented in the legend so that their numbering is meaningful 

library(ggplot2)

ggplot(niches_fossils_strat, aes(x = t, fill = niche, color = niche)) +
  geom_density(alpha = 0.12, color = "black") +
  scale_fill_discrete(drop = FALSE) +
  scale_color_discrete(drop = FALSE) +
  labs(
    title = "Probability for fossil presence 1.5 km from shore",
    x = "Stratigraphic position [m]",
    y = "Density probability"
  )

#water depth as presented by the stratigraphy 

wds <- list("t" = t, "y" = wd) |>    # create list with time - water depth information
  time_to_strat(strat_adm)    # transform into the strat. domain
  
plot(wds, orientation = "du",    # plot water depth information in the stratigraphic domain
       type = "l",
       ylab = "Stratigraphic position [m]",
       xlab = "Water depth [m]",
       main = "Water depth in section 1.5 km from shore") 

#wds <- list("t" = t, "y" = wd) |> # create list with time - water depth information
#  time_to_strat(strat_adm) # transform into the strat. domain

wds <- as.data.frame(wds)

#to be p2
ggplot(wds, aes(x = h, y = y)) +
  geom_line() +
  labs(
    x = "Stratigraphic position [m]",
    y = "Water depth [m]",
    title = "Water depth in section 1.5 km from shore"
  ) +
  theme_minimal()

#niche model in terms of last occurrences (strat domain)

set.seed(124)

niches_applied_strat_occ <- list()
for (i in 1:100) {
  niches_applied_strat_occ[[i]] <- p3(rate = 300, from = min(t), to = max(t)) |> 
    apply_niche(niche_def = niches[[i]], gc = gc) |>                    
    time_to_strat(strat_adm, destructive = FALSE, out_dom_val_h = "default") # transform into strat. domain, do not destroy fossils that coincide with hiatuses for last occurrences
}

niches_strat_occ <- do.call(rbind, lapply(seq_along(niches_applied_strat_occ), function(i) {
  x <- as.numeric(niches_applied_strat_occ[[i]])
  if (length(x) > 0 && any(!is.na(x))) {
    data.frame(t = x, niche = as.factor(i))
  } else {
    # Insert a row with NA to ensure this niche appears in the plot/legend
   data.frame(t = NA, niche = as.factor(i))
  }
}))

#AI for how to find the min for each grouped niche and how to plot them based on their niche group
library(dplyr)

niches_strat_occ <- niches_strat_occ |>
  mutate(
    niche_num = as.integer(as.character(niche)),
    niche_group = cut(niche_num,
                      breaks = seq(0, 100, by = 10),
                      labels = paste(seq(1, 91, by = 10), seq(10, 100, by = 10), sep = "-"),
                      right = TRUE,
                      include.lowest = TRUE)
  )

niches_last_occ <- niches_strat_occ |>
  filter(!is.na(t)) |>
  group_by(niche, niche_group) |>
  summarise(min_t = min(t), .groups = "drop")

library(ggplot2)

#p1 <- 
ggplot(niches_last_occ, aes(x = min_t, fill = niche_group)) +
  geom_histogram(binwidth = 1, color = "black", position = "identity") +
  scale_fill_viridis_d() +
  labs(
    x = "Stratigraphic height [m]",
    y = "# Last occurrences",
    fill = "Niche group",
    title = "Last occurrences 1.5 km offshore"
  ) +
  theme_minimal() +
  coord_cartesian(xlim = c(min(Anna_adm$adm1..m.), max(Anna_adm$adm1..m.)))

#together with water depth
install.packages("remotes")
install.packages("gridGraphics")
remotes::install_version("cowplot", version = "0.9.4", repos = "http://cran.us.r-project.org")
library(cowplot)
plot_grid(p1, p2, ncol = 1)


#"normal" histograms 

#niches_strat_occ <- do.call(rbind, lapply(seq_along(niches_applied_strat_occ), function(i) {
#  x <- as.numeric(niches_applied_strat_occ[[i]])
#  if (length(x) > 0 && any(!is.na(x))) {
#    data.frame(t = x, niche = as.factor(i))
#  } else {
#    # Insert a row with NA to ensure this niche appears in the plot/legend
#   data.frame(t = NA, niche = as.factor(i))
#  }
#}))

#AI for how to find the min for each grouped niche
#library(dplyr)
#niches_last_occ <- niches_strat_occ |>
#  filter(!is.na(t)) |>
#  group_by(niche) |>
#  summarise(min_t = min(t)) |>
#  pull(min_t)

#hist(niches_last_occ,
#     xlab = "Stratigraphic height [m]",
#     main = "Last occurrences 1.5 km offshore",
#     ylab = "# Last occurrences",
#     breaks = seq(from = min(Anna_adm$adm1..m.), to = max(Anna_adm$adm1..m.), length.out = 100))
