# Raman Spectral Analysis for Alpha-Synuclein Fibrillation Assay
# Author: Dr Nathan Coles
# Date: 08/01/25
# Description:
# This script performs data preprocessing, statistical analysis, 
# PCA, UMAP visualization, and random forest modeling on Raman spectral data.

# ------------------------------
# Step 1: Install Required Libraries
# ------------------------------
# The following line ensures required libraries are installed.
install.packages(c("dplyr", "openxlsx", "reshape2", "readxl", "ggplot2", 
                   "factoextra", "cowplot", "gridExtra", "randomForest", 
                   "caret", "umap", "RamanMP"))

# ------------------------------
# Step 2: Load Libraries
# ------------------------------
# Data Manipulation
library(dplyr)          # Provides tools for data wrangling and manipulation

# File Handling
library(openxlsx)       # For reading and writing Excel files
library(readxl)         # For reading Excel files

# Data Reshaping
library(reshape2)       # Functions for reshaping data (e.g., wide to long format)

# Raman-Specific Functions
library(RamanMP)        # Specialized package for Raman spectral data analysis

# Visualization
library(ggplot2)        # Core package for creating plots
library(cowplot)        # For combining multiple ggplot2 plots into a single figure
library(gridExtra)      # Advanced grid-based plotting
library(factoextra)     # Visualization tools for PCA and other multivariate analyses

# Statistical Analysis
library(dunn.test)      # Post hoc analysis using Dunn's test

# Dimensionality Reduction
library(umap)           # For UMAP visualization of complex datasets

# Machine Learning
library(randomForest)   # Random forest for classification and regression
library(caret)          # Framework for cross-validation and model training

# Miscellaneous
library(GGally)         # Creates pairwise plots for exploratory data analysis
# ---------------------------------





# Data Preprocessing and Reshaping for Raman Spectral Analysis

# ---------------------------------
# Set Working Directory
# ---------------------------------
# Users should set their working directory to the folder containing their data files.
setwd("/path/to/your/data/")

# ---------------------------------
# Load and Normalize Data
# ---------------------------------
# Load Raman spectral data from an Excel file
Protein <- read_excel("Protein.xlsx")

# Normalize the data using Standard Normal Variate (SNV) normalization
Protein_result <- norm.SNV(as.matrix(Protein))  # Applies SNV normalization
Protein_result_2 <- as.data.frame(Protein_result)  # Convert matrix to data frame for further processing

# Save normalized data to an Excel file for review
write.xlsx(Protein_result_2, "Protein_result.xlsx")  

# Reload the normalized data for analysis
Protein_result_Final <- read_excel("Protein_result.xlsx")

# ---------------------------------
# Define Data Reshaping Function
# ---------------------------------
# This function reshapes the data to create a tidy format with replicates and frequency values
reshape_data <- function(data) {
  # Number of replicates per timepoint
  n_replicates <- 10
  
  # Extract the frequency values to use as column names
  freq_values <- data$freq[1:1024]
  
  # Initialize an empty data frame for reshaped data
  reshaped_data <- data.frame()
  
  # Counter for naming samples
  sample_counter <- 1
  
  # Loop through each timepoint column in the dataset
  for (day in names(data)[-1]) {
    # Select the data for the current day
    day_data <- data %>% select(freq, !!sym(day))
    
    # Split the data into replicates
    for (i in 1:n_replicates) {
      start_row <- (i - 1) * 1024 + 1
      end_row <- i * 1024
      replicate_data <- day_data[start_row:end_row, ]
      replicate_data <- t(replicate_data[, -1])  # Transpose data for correct format
      replicate_df <- data.frame(replicate_data)
      colnames(replicate_df) <- freq_values  # Assign frequency values as column names
      replicate_df$Sample <- paste0("Sample ", sample_counter)  # Add sample identifier
      sample_counter <- sample_counter + 1
      reshaped_data <- bind_rows(reshaped_data, replicate_df)
    }
  }
  
  return(reshaped_data)
}

# ---------------------------------
# Reshape the Normalized Data
# ---------------------------------
Protein_grouping <- reshape_data(Protein_result_Final)  # Apply reshaping function

# Reorder columns to place 'Sample' as the first column
Protein_grouping <- Protein_grouping %>% select(Sample, everything())

# Add a row with frequency labels
freq_values <- Protein_result_Final$freq[1:1024]
wave_row <- c("#Wave", as.character(freq_values))
Protein_grouping <- rbind(wave_row, Protein_grouping)

# Remove the duplicated frequency row (first row)
Protein_grouping <- Protein_grouping[-1, ]

# Save the reshaped data to a CSV file for further analysis
write.csv(Protein_grouping, "Protein_grouping.csv", row.names = FALSE)

# ---------------------------------
# Verify Reshaped Data
# ---------------------------------
# Display the first few rows of the reshaped data to verify correctness
head(Protein_grouping)

# ---------------------------------
# Melt the Normalized Data
# ---------------------------------
# Convert the data to long format for plotting or further statistical analysis
Protein_melted_result <- melt(Protein_result_Final, id.vars = "freq")  # 'freq' is retained as an identifier
# ---------------------------------





# ---------------------------------
# Raman Spectral Data Visualization
# Description: Functions for plotting Raman spectral data, including line plots, stacked plots, and heatmaps.
# ---------------------------------

# Define color mapping for treatment days
variable_colors <- c("day_0" = "blue", "day_1" = "purple", "day_4" = "orange", "day_7" = "darkkhaki")

# ---------------------------------
# Line Plot of SNV-Normalised Raman Data
# ---------------------------------
# Function to plot average Raman spectral intensities with customizable legend labels
plot_SNVRaman <- function(data_subset, title_suffix, colors, legend_title = "Treatment", legend_labels = NULL) {
  # Calculate average intensities for each frequency and treatment group
  averages <- data_subset %>% 
    group_by(freq, variable) %>% 
    summarize(mean_value = mean(value, na.rm = TRUE)) %>%
    ungroup()  # Ensure data is ungrouped for plotting
  
  # Generate the line plot
  ggplot(averages, aes(x = freq, y = mean_value, color = variable)) +
    geom_line(size = 1.2) +
    labs(
      x = "Wavenumber (cm^-1)",
      y = "Average Intensity (SNV-normalised)",
      color = legend_title
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 30),
      axis.title = element_text(size = 30),
      legend.text = element_text(size = 30),
      axis.line = element_line(size = 1.5)
    ) +
    scale_color_manual(values = colors, labels = legend_labels) +  # Customize legend labels
    scale_x_continuous(limits = c(600, 1750), breaks = seq(600, 1750, by = 100)) +
    scale_y_continuous(labels = NULL)  # Suppress y-axis tick labels for simplicity
}

# Define custom legend labels
legend_labels <- c("Day 0", "Day 1", "Day 4", "Day 7")

# Generate and display the line plot
print(plot_SNVRaman(Protein_melted_result, "Protein", variable_colors, "Day", legend_labels))

# ---------------------------------
# Vertically Stacked Line Plot
# ---------------------------------
# Function to create stacked line plots with vertical offsets for treatment groups
plot_Stacked_SNVRaman <- function(data_subset, title_suffix, colors, legend_title = "Treatment", legend_labels = NULL) {
  # Calculate average Raman intensity for each frequency and treatment group
  averages <- data_subset %>% 
    group_by(freq, variable) %>% 
    summarize(mean_value = mean(value, na.rm = TRUE)) %>%
    ungroup()
  
  # Add vertical offsets for stacking
  vertical_offsets <- c("day_0" = 0, "day_1" = 5, "day_4" = 10, "day_7" = 15)
  averages <- averages %>%
    mutate(mean_value_stacked = mean_value + vertical_offsets[variable])
  
  # Generate the stacked line plot
  ggplot(averages, aes(x = freq, y = mean_value_stacked, color = variable)) +
    geom_line(size = 1.2) +
    labs(
      title = paste("Vertically Stacked Plot of SNV Normalised Raman Spectral Intensities -", title_suffix),
      x = "Wavenumber (cm^-1)",
      y = "Stacked Intensity (SNV-normalised)",
      color = legend_title
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 20),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      axis.line = element_line(size = 1.5)
    ) +
    scale_color_manual(values = colors, labels = legend_labels) +  # Customize legend labels
    scale_x_continuous(limits = c(600, 1750), breaks = seq(600, 1750, by = 100)) +
    scale_y_continuous(labels = NULL)  # Suppress y-axis tick labels for aesthetics
}

# Generate and display the stacked plot
print(plot_Stacked_SNVRaman(Protein_melted_result, "Protein", variable_colors, "Day", legend_labels))

# ---------------------------------
# Heatmap of Raman Intensities
# ---------------------------------
# Function to create a heatmap of Raman spectral intensities with custom sample titles
create_heatmap <- function(data_subset, title, sample_titles) {
  ggplot(data_subset, aes(x = freq, y = variable, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-2, 9)) +
    labs(
      title = title,
      x = "Wavenumber (cm^-1)",
      y = "Treatment"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
      axis.text.y = element_text(size = 15),
      legend.text = element_text(size = 15),
      axis.line = element_line(size = 1.5)
    ) +
    scale_x_continuous(breaks = seq(600, 1750, by = 100)) +
    guides(fill = guide_colorbar(title = "Intensity")) +  # Add legend title
    scale_y_discrete(labels = sample_titles)  # Apply custom sample titles
}

# Define sample titles
sample_titles <- c("Day 0", "Day 1", "Day 4", "Day 7")

# Generate and display the heatmap
heatmap_plot_Protein <- create_heatmap(Protein_melted_result, "Heatmap of Raman Intensities - Protein", sample_titles)
print(heatmap_plot_Protein)
# ---------------------------------





# ---------------------------------
# PCA and Scree Plot Analysis for Raman Spectral Data
# Description: This script performs PCA analysis on Raman spectral data, generates scree plots, 
#              calculates cumulative variance explained, and visualizes PCA results.
# ---------------------------------
# Define File Paths
# Update these paths to match your local directory structure
raman_file_path <- "/path/to/your/data/"
grouping_file_path <- "/path/to/your/data/"

# Load Raman Spectral Data and Grouping Information
raman_data <- read.csv(raman_file_path, header = TRUE, stringsAsFactors = FALSE)

# Remove the first column (sample numbers) and ensure all data columns are numeric
raman_data <- apply(raman_data[, -1], 2, as.numeric)

# Load grouping information
grouping <- read.csv(grouping_file_path)

# ---------------------------------
# Principal Component Analysis (PCA)
# Perform PCA on the Raman spectral data
pca <- prcomp(raman_data, scale. = TRUE)

# ---------------------------------
# Scree Plot for PCA
# Create a data frame of variance explained by each principal component
pca_data <- data.frame(
  Principal_Component = 1:length(pca$sdev),
  Variance_Explained = pca$sdev^2
)

# Subset the data to include only the first 20 principal components
pca_data_subset <- pca_data[1:20, ]

# Generate a scree plot for the first 20 principal components
ggplot(pca_data_subset, aes(x = Principal_Component, y = Variance_Explained)) +
  geom_point(shape = 21, size = 4, fill = "black", color = "black") +  # Filled circles
  geom_path(color = "black") +  # Connecting lines
  labs(
    x = "Number of Principal Components",
    y = "Eigenvalue",
    title = "Scree Plot Demonstrating Eigenvalues for Protein"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ---------------------------------
# Cumulative Variance Explained
# Calculate cumulative variance explained by principal components
cumulative_variance <- cumsum(summary(pca)$importance[2, ])

# Create a data frame for the cumulative variance plot
cumulative_data <- data.frame(
  Principal_Component = 1:length(cumulative_variance),
  Cumulative_Variance = cumulative_variance * 100  # Convert to percentage
)

# Subset the data to include only the first 20 principal components
cumulative_data_2 <- cumulative_data[1:20, ]

# Generate a cumulative scree plot
ggplot(cumulative_data_2, aes(x = Principal_Component, y = Cumulative_Variance)) +
  geom_line(color = "black") +
  geom_point(shape = 21, size = 4, fill = "black", color = "black") +
  labs(
    x = "Number of Principal Components",
    y = "Variance Explained (%)",
    title = "Cumulative Variance Explained by Principal Components"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ---------------------------------
# PCA Result Visualization
# Check the number of principal components
num_components <- ncol(pca$x)

# Create a data frame with PCA results and grouping information
pc_data <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = if (num_components >= 2) pca$x[, 2] else rep(0, nrow(pca$x)),
  group = as.factor(grouping$Group)
)

# Get the proportion of variance explained for PC1 and PC2
pca_var_raw <- summary(pca)$importance[2, ]

# Visualize PCA results with ellipses for group clustering
ggplot(pc_data, aes(x = PC1, y = PC2, color = group, fill = group)) +
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", level = 0.6, alpha = 0, linewidth = 1.5) +
  scale_color_manual(values = c("blue", "purple", "orange", "darkkhaki")) +
  labs(
    x = paste0("PC1 (", round(pca_var_raw[1] * 100, 2), "% Variance)"),
    y = paste0("PC2 (", round(pca_var_raw[2] * 100, 2), "% Variance)"),
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(size = 1.5),
    text = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.text = element_text(size = 20)
  )
# ---------------------------------





# ---------------------------------
# PCA Pairwise Plot Visualization
# Description: This script generates pairwise plots of principal components (PC1, PC2, PC3, etc.)
#              using grouping information to visualize clustering in PCA results.
# ---------------------------------
# Merge PCA Results with Grouping Information
if (num_components >= 2) {
  # Create a data frame for PC1 and PC2
  pc_data <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
} else {
  # Handle cases with fewer than 2 components
  pc_data <- data.frame(PC1 = pca$x[, 1], PC2 = rep(0, length(pca$x[, 1])))
}

# Add grouping information
pc_data$group <- as.factor(grouping$Group)

# Add PC3 if available
if (num_components >= 3) {
  pc_data$PC3 <- pca$x[, 3]
}

# ---------------------------------
# Inspect the Data
# Display the structure of the data frame
str(pc_data)

# View the first few rows of the data frame
head(pc_data)

# ---------------------------------
# Pairwise Plot of Principal Components
# Create pairwise plots for the first few principal components
# Color-code points by their group for better visual clustering
ggpairs(
  pc_data, 
  columns = 1:4,  # Include PC1, PC2, PC3, and group if available
  aes(color = group)  # Color by group
) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )




# --------------------------------------------
# Loadings Analysis and Visualization for PCA
# --------------------------------------------
# This script performs the following tasks:
# 1. Identifies influential wavenumbers in PCA loadings.
# 2. Visualizes the loadings with highlighted influential wavenumbers.
# 3. Performs statistical analysis (ANOVA, Kruskal-Wallis, Dunn's test).
# --------------------------------------------
# Convert row names of PCA loadings to numeric wavenumbers
loadings <- as.data.frame(pca$rotation)
loadings$Wavenumber <- as.numeric(sub("X", "", rownames(loadings)))

# Ensure no NA values in Wavenumber column
if (any(is.na(loadings$Wavenumber))) {
  stop("NA values found in Wavenumber. Check row names.")
}

# Set threshold for influential loadings
threshold <- 0.05

# Identify influential wavenumbers
influential_wavenumbers <- lapply(loadings[, -ncol(loadings)], function(x) abs(x) > threshold)

# Plot loadings and highlight influential wavenumbers
for (i in 1:min(4, ncol(loadings) - 1)) {
  # Prepare data for plotting
  loading_data <- data.frame(
    Wavenumber = loadings$Wavenumber,
    Loading = loadings[, i],
    Influential = influential_wavenumbers[[i]]
  )
  
  # Add sorting and labeling
  loading_data <- loading_data[order(loading_data$Wavenumber), ]
  loading_data$Influential <- factor(loading_data$Influential, levels = c(TRUE, FALSE))
  
  # Plot loadings
  ggplot(loading_data, aes(x = Wavenumber, y = Loading, color = Influential)) +
    geom_hline(yintercept = c(-threshold, threshold), linetype = "dashed") +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
    labs(
      title = paste("Loadings for PC", i, "- Protein"),
      x = "Wavenumber (cm^-1)",
      y = "Loading"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 20),
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      axis.line = element_line(size = 1.5)
    ) +
    scale_x_continuous(breaks = seq(600, 1750, by = 100)) +
    scale_y_continuous()
}

# --------------------------------------------
# Raman Spectral Intensities Analysis
# --------------------------------------------
# Analyze a specific frequency (e.g., 725.901367 cm^-1)
specific_freq <- Protein_melted_result %>%
  filter(freq == 725.901367)

# Define day colors and calculate mean intensities
days_to_plot <- unique(specific_freq$variable)
day_colors <- c("blue", "purple", "orange", "darkkhaki")

# Prepare and plot jitter plots
plots_list <- list()
for (i in seq_along(days_to_plot)) {
  day <- days_to_plot[i]
  protein_plot <- ggplot(specific_freq %>% filter(variable == day), aes(x = variable, y = value)) +
    geom_jitter(width = 0.1, alpha = 0.7, size = 3, color = day_colors[i]) +
    geom_hline(yintercept = mean(specific_freq$value[specific_freq$variable == day]), color = "black") +
    theme_minimal() +
    labs(x = "", y = "") +
    theme(axis.text.x = element_text(size = 15))
  plots_list[[i]] <- protein_plot
}

# Combine plots into a grid
combined_plots <- plot_grid(plotlist = plots_list, ncol = 2)
print(combined_plots)

# --------------------------------------------
# Statistical Tests
# --------------------------------------------
# Perform one-way ANOVA
anova_results <- aov(value ~ variable, data = specific_freq)
anova_summary <- summary(anova_results)

# Normality test
shapiro_test <- shapiro.test(residuals(anova_results))

# Perform Kruskal-Wallis test for non-parametric data
kruskal_test <- kruskal.test(value ~ variable, data = specific_freq)

# Perform Dunn's test for post hoc comparisons
dunn_test_result <- dunn.test(specific_freq$value, specific_freq$variable, method = "bonferroni")

# Print results
cat("ANOVA Summary:\n")
print(anova_summary)
cat("\nShapiro Test for Normality:\n")
print(shapiro_test)
cat("\nKruskal-Wallis Test:\n")
print(kruskal_test)
cat("\nDunn's Test:\n")
print(dunn_test_result)

# --------------------------------------------





# ---------------------------------
# UMAP Analysis and Random Forest Classification for Protein Data
# Description: This script performs UMAP dimensionality reduction and trains a Random Forest model for supervised learning.
# --------------------------------------------


# --- Initialization ---
num_seeds <- 10  # Number of UMAP seeds for averaging
all_umap_results <- list()


# --- UMAP Dimensionality Reduction ---
for (i in 1:num_seeds) {
  set.seed(123)  # Ensure reproducibility
  umap_result <- umap(raman_data)  # Perform UMAP
  all_umap_results[[i]] <- data.frame(
    UMAP1 = umap_result$layout[, 1],
    UMAP2 = umap_result$layout[, 2]
  )
}

# Average UMAP results over multiple seeds
average_umap <- Reduce("+", all_umap_results) / num_seeds
average_umap$group <- as.factor(grouping$Group)  # Add group labels

# --- Define Custom Colors ---
custom_colors <- c("Day 0" = "blue", "Day 1" = "purple", "Day 4" = "orange", "Day 7" = "darkkhaki")

# --- UMAP Visualization ---
umap_plot <- ggplot(average_umap, aes(x = UMAP1, y = UMAP2, color = group, fill = group)) +
  geom_point(size = 3) +
  stat_ellipse(aes(color = group, fill = group),
               geom = "polygon", level = 0.6, alpha = 0, linewidth = 1.5) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.line = element_line(size = 1.5)
  ) +
  labs(
    title = "UMAP Plot for Protein",
    x = "UMAP1", y = "UMAP2", color = "Group", fill = "Group"
  )
print(umap_plot)

# --- Supervised Learning with Random Forest ---
# Define Random Forest parameters
num_trees <- 100  # Number of trees
num_folds <- 10   # Number of folds for cross-validation

# Define control parameters for cross-validation
train_control <- trainControl(method = "cv", number = num_folds)

# Train Random Forest model with cross-validation
set.seed(123)  # Ensure reproducibility
rf_model <- train(
  group ~ UMAP1 + UMAP2,
  data = average_umap,
  method = "rf",
  ntree = num_trees,
  trControl = train_control
)

# Print model summary
print(rf_model)

# Extract cross-validated metrics
accuracy <- rf_model$results$Accuracy
kappa <- rf_model$results$Kappa

# Print model performance
cat("Random Forest Cross-Validated Accuracy:\n")
print(accuracy)
cat("\nKappa Statistics:\n")
print(kappa)
# --------------------------------------------











