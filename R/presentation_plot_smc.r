presentation_plot_smc<-function(){# Load necessary libraries
library(ggplot2)
library(dplyr)
library(cowplot)

# Set seed for reproducibility
set.seed(123)

# Parameters
N <- 20  # Number of particles
time_t <- 0  # Initial time step

# Generate initial particles
theta_t <- rnorm(N, mean = 0, sd = 1)
weights <- dnorm(theta_t)

# Define a model of transition for SMC
transition_model <- function(theta) {
  theta + rnorm(length(theta), mean = 0, sd = 0.5)
}

# Define a likelihood function for weights calculation
likelihood <- function(theta) {
  dnorm(theta, mean = 0, sd = 1)
}

# Resampling step based on weights
resample <- function(particles, weights) {
  sample(particles, size = length(particles), replace = TRUE, prob = weights)
}

# SMC step
resampled_particles <- resample(theta_t, weights)
new_particles <- transition_model(resampled_particles)

weights2 <-rep(mean(weights),N)
weights3 <- dnorm(new_particles)
# Prepare data for plotting
data <- data.frame(
  x = ordered(rep(c("t", "Resampled", "t+1"), each = N),c("t", "Resampled", "t+1")),
  xx = rep(1:3,each=N),
  weight=c(weights,weights2,weights3),
  y = c(theta_t, resampled_particles, new_particles),
  group = rep(1:N, 3)
)|>
  dplyr::mutate(xx=xx+rnorm(length(xx),sd=.15*exp(-y^4)))

# Target distribution for density plot
theta_values <- seq(min(new_particles) - 1, max(new_particles) + 1, length.out = 100)
target_density <- dnorm(theta_values, mean = 0, sd = 1)

# Plotting the SMC steps
p1 <- ggplot(data, aes(x = xx, 
                       y = y, group = group)) +
  geom_point(aes(size=weight)) +
  geom_line(linetype = "dashed") +
  labs(
    title = "Sequential Monte Carlo Step",
    x = "Step",
    y = expression(theta)
  ) +
  theme_minimal()+ylim(-2,2)+
  scale_x_continuous(breaks = 1:3, labels = c("t", "Resampled", "t+1")) 

# Plotting the target density
p2 <- ggplot(data.frame(x = theta_values, y = target_density), aes(y, x)) +
  geom_path() +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ylim(-2,2)

# Combine plots using cowplot
combined_plot <- plot_grid(
  p1, p2, 
  ncol = 2, 
  align = "h", 
  rel_widths = c(3, 1)#,
#  labels = c("", "Density"),
#  label_x = c(0.2, 0.5), 
#  label_y = c(1, 0.5)
)

# Display combined plot
print(combined_plot)
library(tikzDevice)
tikz("smc_plot.tex", standAlone = TRUE, width = 8, height = 4)
print(combined_plot)
dev.off()
}