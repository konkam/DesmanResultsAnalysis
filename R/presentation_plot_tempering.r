plot_tempering<-function(){# Load necessary libraries
library(ggplot2)
library(dplyr)

# Set seed for reproducibility
set.seed(123)

# Define the range for the x-axis
x <- seq(-5, 5, length.out = 1000)

# Define unimodal and bimodal distributions
unimodal <- function(x) {
  dnorm(x, mean = 0, sd = 1)
}

bimodal <- function(x) {
  0.5 * dnorm(x, mean = -2, sd = 0.5) + 0.5 * dnorm(x, mean = 2, sd = 0.5)
}

# Define a function to create intermediate distributions
intermediate <- function(x, alpha) {
  (1 - alpha) * unimodal(x) + alpha * bimodal(x)
}

# Generate data for the intermediate distributions
alphas <- seq(0, 1, length.out = 11)
data <- data.frame(
  x = rep(x, times = length(alphas)),
  density = unlist(lapply(alphas, function(a) intermediate(x, a))),
  alpha = factor(rep(alphas, each = length(x)))
)

# Plot the transition from unimodal to bimodal distribution
tempering_plot=ggplot(data, aes(x = x, y = density, color = alpha)) +
  geom_line() +
  scale_color_viridis_d(name = "Tempering\n(alpha)") +
  labs(
    title = "Transition progressive (Tempering)",
    x = "x",
    y = "Density"
  ) +
  theme_minimal()
library(tikzDevice)
tikz("tempering.tex", standAlone = TRUE, width = 8, height = 4)
print(tempering_plot)
dev.off()
}