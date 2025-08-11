# July 28 2025
# Visualization of vcfeval results for Independent Study

setwd("~/uleth/independent-study/project")

# Import results
library(tidyverse)

data <- read.csv("./vcfeval-data-P1.tsv", header = TRUE, sep = "\t")

# Add sample-ID column to df
library(dplyr)

data <- data %>%
  mutate(`Sample-ID` = paste(`Variant_Caller`, `Variant_Level`, Threshold, sep = "_")) %>%
  select(`Sample-ID`, everything())

# Plot sensitivity vs precision
library(ggplot2)
  
  # Filter data for Q0 threshold only
  q0_data <- data[data$Threshold == "Q0", ]
  # Create plot
  p <- ggplot(q0_data, aes(x = Precision, y = Sensitivity, color = Variant_Caller)) +
    geom_point(size = 4, alpha = 0.8, 
               position = position_dodge(width = 0.02)) +
    # Add manual colours
    scale_color_manual(values = c("clair3" = "#49043c", 
                                  "freebayes" = "#FD5E53", 
                                  "lofreq" = "cyan"),
                       name = " ") +
    # Create facet with Variant_Level
    facet_wrap(~ Variant_Level, ncol = 2, 
               labeller = labeller(Variant_Level = function(x) paste("Variant Level:", x, "%"))) +
    # Axis labels
    labs(
         x = "Precision",
         y = "Sensitivty") +
    
    xlim(0, 1) + ylim(0, 1) +
    # Apply serif font and clean theme
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 12, face = "bold", family = "serif"),
          strip.text = element_text(size = 12, face = "bold", family = "serif"),
          axis.title.x = element_text(size = 12, face = "bold", family = "serif"),
          axis.title.y = element_text(size = 12, face = "bold", family = "serif"),
          axis.text.x = element_text(family = "serif"),
          axis.text.y = element_text(family = "serif"),
          panel.border = element_rect(color = "darkgray", fill = NA, linewidth = 0.6),
          panel.background = element_rect(fill = "white", color = NA),
          strip.background = element_rect(fill = "#f0f0f0", color = "darkgrey"))
  
  print(p)

  
  
  
  # Plot F1 score, sensitivity, and precision against Q for each variant caller
  
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  # Reshape data from wide to long format for plotting multiple metrics
  data_long <- data %>%
    pivot_longer(cols = c(Precision, Sensitivity, F1_score), 
                 names_to = "Metric", 
                 values_to = "Value")
  
  # Convert Threshold to numeric for proper ordering
  # Extract numeric part from Q0, Q5, Q10, etc.
  data_long$Threshold_num <- as.numeric(gsub("Q", "", data_long$Threshold))
  
  # Create the plot
  p2 <- ggplot(data_long, aes(x = Threshold_num, y = Value, color = Metric)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(size = 2.5, alpha = 0.9) +
    
    # Custom colors for each metric
    scale_color_manual(values = c("Precision" = "#E31A1C",    # Red
                                  "Sensitivity" = "#1F78B4",   # Blue  
                                  "F1_score" = "#33A02C"),     # Green
                       name = "Performance Metric",
                       labels = c("F1_score" = "F1 Score", 
                                  "Precision" = "Precision",
                                  "Sensitivity" = "Sensitivity")) +
    
    # Create facets with rows by Variant_Caller, columns by Variant_Level
    facet_grid(Variant_Caller ~ paste("Variant Level", Variant_Level),
               scales = "free_x") +
    
    # Customize x-axis
    scale_x_continuous(name = "Quality Threshold", 
                       breaks = c(0, 5, 10, 15, 20),
                       labels = c("Q0", "Q5", "Q10", "Q15", "Q20")) +
    
    # Set y-axis from 0 to 1
    scale_y_continuous(name = "Performance Score", 
                       limits = c(0, 1),
                       breaks = seq(0, 1, 0.25)) +
    
    # Main title
    labs(title = "Performance Metrics Across Quality Thresholds",
         subtitle = "Precision, Sensitivity, and F1 Score by Variant Caller and Variant Level") +
    
    # Apply serif font and clean theme
    theme_minimal(base_family = "serif") +
    theme(
      # Title formatting
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 20)),
      
      # Axis formatting
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      
      # Facet labels
      strip.text = element_text(size = 11, face = "bold", margin = margin(5, 5, 5, 5)),
      strip.background = element_rect(fill = "grey90", color = "grey70"),
      
      # Legend
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.position = "bottom",
      legend.box = "horizontal",
      
      # Panel formatting
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey95", size = 0.5),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
      
      # Overall margins
      plot.margin = margin(20, 20, 20, 20)
    )
  
  # Display the plot
  print(p2)
  
  
  
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  # Reshape data from wide to long format for plotting multiple metrics
  data_long <- data %>%
    pivot_longer(cols = c(Precision, Sensitivity, F1_score), 
                 names_to = "Metric", 
                 values_to = "Value")
  
  # Convert Threshold to numeric for proper ordering
  # Extract numeric part from Q0, Q5, Q10, etc.
  data_long$Threshold_num <- as.numeric(gsub("Q", "", data_long$Threshold))
  
  # Set factor levels for Variant_Level to control column order (20, 15, 10, 5)
  data_long$Variant_Level <- factor(data_long$Variant_Level, levels = c("20", "15", "10", "5"))
  
  # Create the plot
  p3 <- ggplot(data_long, aes(x = Threshold_num, y = Value, color = Metric)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(size = 2.5, alpha = 0.9) +
    
    # Custom colors for each metric
    scale_color_manual(values = c("Precision" = "orange",
                                  "Sensitivity" = "#1F78B4", 
                                  "F1_score" = "#33A02C"),
                       name = "Performance Metric",
                       labels = c("F1_score" = "F1 Score", 
                                  "Precision" = "Precision",
                                  "Sensitivity" = "Sensitivity")) +
    
    # Create facets with rows by Variant_Caller, columns by Variant_Level
    facet_grid(Variant_Caller ~ factor(paste("Variant Level", Variant_Level, "%"), 
                                       levels = c("Variant Level 20 %", "Variant Level 15 %", 
                                                  "Variant Level 10 %", "Variant Level 5 %")),
               scales = "free_x") +
    
    # Customize x-axis
    scale_x_continuous(name = "Quality Threshold", 
                       breaks = c(0, 5, 10, 15, 20),
                       labels = c("Q0", "Q5", "Q10", "Q15", "Q20")) +
    
    # Set y-axis from 0 to 1
    scale_y_continuous(name = "Performance Score", 
                       limits = c(0, 1),
                       breaks = seq(0, 1, 0.25)) +
    
    # Apply serif font and clean theme
    theme_minimal(base_family = "serif") +
    theme(
      # Axis formatting
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      
      # Facet labels
      strip.text = element_text(size = 11, face = "bold", margin = margin(5, 5, 5, 5)),
      strip.background = element_rect(fill = "grey90", color = "grey70"),
      
      # Legend
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.position = "bottom",
      legend.box = "horizontal",
      
      # Panel formatting
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey95", size = 0.5),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
      
      # Overall margins
      plot.margin = margin(20, 20, 20, 20)
    )
  
  # Display the plot
  print(p3)
  
