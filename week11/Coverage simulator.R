library(ggplot2)
setwd("/Users/cmdb/qb25-answers/week11")

plot_coverage_histogram <- function(coverage_level) {
  # Read coverage data
  coverage_file <- paste0('coverage_', coverage_level, 'x.txt')
  coverage_data <- scan(coverage_file, quiet = TRUE)
  # Read distribution data
  dist_file <- paste0('distributions_', coverage_level, 'x.txt')
  distributions <- read.table(dist_file, header = TRUE)
  # Histogram
  max_cov <- max(coverage_data)
  hist_data <- table(factor(coverage_data, levels = 0:max_cov))
  hist_df <- data.frame(
    coverage = 0:max_cov,
    count = as.numeric(hist_data))
  genome_size <- length(coverage_data)   # Scale distribution
  distributions$poisson_scaled <- distributions$poisson * genome_size
  distributions$normal_scaled <- distributions$normal * genome_size
  
  p <- ggplot() +
    geom_bar(data = hist_df, 
             aes(x = coverage, y = count, fill = "Simulation"), 
             stat = "identity", alpha = 0.7) +
    geom_line(data = distributions, 
              aes(x = coverage, y = poisson_scaled, color = "Poisson"), 
              size = 1.2) +
    geom_point(data = distributions, 
               aes(x = coverage, y = poisson_scaled, color = "Poisson"), 
               size = 2) +
    geom_line(data = distributions, 
              aes(x = coverage, y = normal_scaled, color = "Normal"), 
              size = 1.2) +
    geom_point(data = distributions, 
               aes(x = coverage, y = normal_scaled, color = "Normal"), 
               size = 2) +
    scale_fill_manual(name = "", values = c("Simulation" = "steelblue")) +
    scale_color_manual(name = "", 
                       values = c("Poisson" = "red", "Normal" = "green4")) +
    labs(title = paste0(coverage_level, "x Coverage Simulation"),
         subtitle = paste0("1 Mbp genome, 100bp reads (lambda = ", coverage_level, ")"),
         x = "Coverage Depth",
         y = "Frequency") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10))
  
  # Save plot
  output_file <- paste0('ex1_', coverage_level, 'x_cov.png')
  ggsave(output_file, plot = p, width = 10, height = 6, dpi = 300)
  }

# Generate all plots
plot_coverage_histogram(3)
plot_coverage_histogram(10)
plot_coverage_histogram(30)