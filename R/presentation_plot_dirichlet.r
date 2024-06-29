# Load necessary libraries
plot_dirichlet<-function(){
library(ggplot2)
library(ggtern)

# Set seed for reproducibility
set.seed(123)

# Number of samples
N <- 100000

# Shape parameters (small values)
alpha <- c(0.05, 0.05, 0.05)

# Generate samples from Dirichlet distribution
rrdirichlet <- function(N, alpha) {
  samples <- matrix(rgamma(N * length(alpha), shape = alpha, scale = 1), ncol = length(alpha), byrow = TRUE)
  samples <- samples / rowSums(samples)
  return(samples)
}

samples <- rrdirichlet(N, alpha)
samples_df <- data.frame(samples)
colnames(samples_df) <- c("A", "C", "G")

samples2 <- rrdirichlet(N, rep(.5,3))
samples_df2 <- data.frame(samples2)
colnames(samples_df2) <- c("A", "C", "G")
# Plotting the Dirichlet distribution on a 3-simplex

ggtern(data = samples_df, aes(x = A, y = C, z = G)) +
  geom_point(alpha = 0.1) +
  labs(
    #title = "Dirichlet Distribution on a 3-Simplex",
    x = "(1,0,0)",
    y = "(0,1,0)",
    z = "(0,0,1)"
  ) +
  theme_minimal() +
  theme_showarrows()


ggtern(data = samples_df2, aes(x = A, y = C, z = G)) +
  geom_point(alpha = 0.01) +
  labs(
   # title = "Dirichlet Distribution on a 3-Simplex",
    x = "(1,0,0)",
    y = "(0,1,0)",
    z = "(0,0,1)"
  ) +
  theme_minimal() +
  theme_showarrows()




project_simplex_3d <- function(point=NULL) {
  if(is.null(point)){point=rdirichlet_smallalpha(5000,alpha = rep(.1,4))}
  A=matrix(
    c(0,0,0,
      1,0,0,
      1/2,sqrt(3)/2,0,
      1/2,sqrt(3)/6,sqrt(6)/3),4,3,byrow=TRUE)
  proj_x=point%*%A
  plyr::aaply(proj_x,1,
                function(x){paste0(
                  "\\fill[rouge_ensae,  opacity=\\alphaa] (",paste(round(x,3),collapse=","),") circle (.2pt);"
)
})|>paste(collapse="
         ")}


project_simplex_3d(diag(4))

project_simplex_3d()|>clipr::write_clip()}
  
