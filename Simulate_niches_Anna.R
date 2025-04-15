library(StratPal) #bio and ecology
library(admtools) #stratigraphy

#task is to create 100 niches described by water depth from CarboKitten model

#Created using chatgpt 

#parameters
set.seed(123)              
n_niches <- 100              # Number of niches
depth_range <- seq(0, 80, length.out = 400)  # Water depths from 0 to 80 meters: to be changed when I get the CarboKitten range

# Generate optima more densely toward shallow depth
# Step 1: Create a uniform sequence between 0 and 1
uniform_seq <- seq(0, 1, length.out = n_niches)
# Step 2: Apply a transformation to bias toward 0 (shallow)
# Try squaring or taking sqrt, log, etc.
optima <- 80 * (uniform_seq)^2  # Quadratic: more values near 0
#other options for bias toward deep: 'sqrt(uniform_seq)' or 'exp(x) / exp(1)'
#for bias toward shallow: 'log1p(uniform_seq) / log1p(1)'

# Define niche width as a function of depth (e.g., linearly increasing)
min_width <- 2   # Narrowest at shallowest depth
max_width <- 20  # Widest at deepest depth
niche_widths <- min_width + (optima / max(optima)) * (max_width - min_width)

# You can change this to a nonlinear function if you'd like (e.g., log, exp, etc.):
#niche_widths <- min_width + log1p(optima) / log1p(max(optima)) * (max_width - min_width)

# Create response matrix
niche_matrix <- matrix(0, nrow = n_niches, ncol = length(depth_range))
for (i in 1:n_niches) {
  niche_matrix[i, ] <- exp(-((depth_range - optima[i])^2) / (2 * niche_widths[i]^2))
}

# Plot to see how optima cluster near shallow depths
matplot(depth_range, t(niche_matrix[seq(1, 80, by = 8), ]), type = "l",
        lty = 1, col = rainbow(10), lwd = 2,
        xlab = "Water Depth (m)", ylab = "Niche Suitability",
        main = "Niches Clustered at Shallow Depths")
legend("topright", legend = paste("Niche", seq(1, 80, by = 8)),
       col = rainbow(10), lty = 1, lwd = 2)

# Plot histogram of niche optima
hist(optima, breaks = 20, col = "skyblue", border = "white",
     main = "Distribution of Niche Optima Across Depth",
     xlab = "Water Depth (m)", ylab = "Number of Niches")