# ==============================================================================
# ===============4==common method bias==========================================
# ==============================================================================

library(openxlsx)      # For reading Excel files
library(psych)         # For Harman's single factor test
library(dplyr)         # For data manipulation

# 1. Clear current environment
rm(list = ls())
# 2. Set working directory
setwd("Your_fold/code")
# 3. Load data
data_origin <- read.xlsx("CMB_test.xlsx", sheet = "revere_collected")   # Load origin data with reversed
# 4. Select variables for analysis
cmb_data <- data_origin[, c('CRP.code', 'TM_1', 'TM_2', 'SM_1', 'SM_2', 'SM_3', 'PO_1', 'PO_2', 
                            'PO_3', 'PO_4', 'PO_5', 'PO_6', 'PO_7', 'PO_8', 'TBS_1', 'TBS_2', 
                            'TBS_3', 'TBS_4', 'TBS_5', 'TBS_6', 'TBS_7', 'DS_1', 'DS_2', 'DS_3',
                            'MC_1', 'MC_2', 'MC_3', 'PS_1', 'PS_2', 'DT_1', 'DT_2', 'DT_3', 'eff',
                            'robot_eff_adding')]
# Remove unnecessary variables
cmb_data <- cmb_data %>% select(-CRP.code, -eff, -robot_eff_adding)

# Check for missing values
if (any(is.na(cmb_data))) stop("Data contains missing values. Please handle them before analysis.")

# 5. Harman's single factor test the reverse coding data
harman_result <- principal(cmb_data, nfactors = 1, rotate = "none")
variance_explained <- harman_result$Vaccounted[2, 1]  # Extract variance explained by single factor

cat("Original data: Variance explained by single factor: ", round(variance_explained * 100, 2), "%\n")
if (variance_explained > 0.40) {
  cat("Warning: Original data may have common method bias.\n")
} else {
  cat("Original data shows no significant common method bias.\n")
}

head(cmb_data$PO_2,30)  # View values before reverse coding
# 6. Reverse coding  , converting  the reversed variables back to their original values
reverse_vars <- c("PO_2", "PO_4", "PO_5","PO_6", "PO_8","PO_7",
                  "TBS_1", "TBS_2", "TBS_4", "TBS_5","TBS_6","TBS_7",
                  "TM_1", "TM_2",
                  "DS_1","DS_3",
                  "PS_1","PS_2",
                  "SM_1",
                  "DT_1","DT_3",
                  "MC_1")

# Only transform variables that need reverse coding, keep others unchanged
cmb_data_2 <- cmb_data %>%
  mutate(across(all_of(reverse_vars), ~ 8 - .))
head(cmb_data_2$PO_2,30)  # View values after reverse coding

# 7. Perform Harman's single factor test the reverse-collected data
harman_result_reversed <- principal(cmb_data_2, nfactors = 1, rotate = "none")
variance_explained_reversed <- harman_result_reversed$Vaccounted[2, 1]

cat("After reverse coding: Variance explained by single factor: ", round(variance_explained_reversed * 100, 2), "%\n")
if (variance_explained_reversed > 0.40) {
  cat("Warning: Data after reverse coding may still have common method bias.\n")
} else {
  cat("Data after reverse coding shows no significant common method bias.\n")
}

#======== Alternative CMB method implementing Lindell & Whitney (2001) approach ===============
# Function to check for missing values
check_missing <- function(data) {
  if(any(is.na(data))) {
    stop("Data contains missing values. Please handle them before analysis.")
  }
}

# Calculate mean correlation coefficient for a single variable
# In Lindell & Whitney (2001) method, we should select the variable with the lowest mean_correlation as marker variable.
# This is because: A marker variable should theoretically have little or no relation to the main constructs in the study
# Lower mean correlation indicates lower association with other variables
# Such variables are more suitable for estimating common method bias impact

calculate_mean_correlation <- function(cor_matrix, var_index) {
  correlations <- cor_matrix[var_index, -var_index]
  mean(abs(correlations))
}

# Find potential marker variable
find_marker_variable <- function(data) {
  # Check data
  check_missing(data)
  
  # Calculate correlation matrix (only once)
  cor_matrix <- cor(data)
  n_vars <- ncol(data)
  
  # Calculate statistics for each variable
  results <- data.frame(
    variable = names(data),
    mean_correlation = sapply(1:n_vars, function(i) calculate_mean_correlation(cor_matrix, i)),
    stringsAsFactors = FALSE
  )
  
  # Add maximum correlation
  results$max_correlation <- apply(abs(cor_matrix), 1, function(x) max(x[x != 1]))
  
  # Sort by mean correlation
  results <- results[order(results$mean_correlation), ]
  return(results)
}

# Perform CMB analysis
perform_cmb_analysis <- function(data, marker_var, threshold = 0.05) {
  # Check if marker variable exists in dataset
  if(!marker_var %in% names(data)) {
    stop("Specified marker variable not found in dataset")
  }
  
  # Get other variables
  other_vars <- setdiff(names(data), marker_var)
  
  # Calculate correlations and p-values
  cor_results <- corr.test(data[c(marker_var, other_vars)])
  cor_matrix <- cor_results$r
  p_matrix <- cor_results$p
  
  # Extract marker variable correlations and p-values
  marker_correlations <- cor_matrix[marker_var, other_vars]
  marker_pvalues <- p_matrix[marker_var, other_vars]
  
  # Create results dataframe
  marker_results <- data.frame(
    variable = other_vars,
    correlation = marker_correlations,
    p_value = marker_pvalues,
    stringsAsFactors = FALSE
  )
  
  # Calculate mean absolute correlation
  mean_correlation <- mean(abs(marker_correlations))
  
  # Return results
  return(list(
    marker_results = marker_results,
    mean_correlation = mean_correlation,
    significant_vars = marker_results[marker_results$p_value < threshold, ]
  ))
}

# Main analysis workflow
main_analysis <- function(data, top_n = 5) {
  # 1. Find potential marker variables
  cat("Starting to find potential marker variables...\n")
  potential_markers <- find_marker_variable(data)
  
  # Print top N potential marker variables
  cat("\n=== Top", top_n, "Potential Marker Variables ===\n")
  print(head(potential_markers, top_n))
  
  # 2. Select variable with lowest correlation as marker variable
  marker_var <- potential_markers$variable[1]
  cat("\nSelected", marker_var, "as marker variable\n")
  
  # 3. Perform CMB analysis
  cat("\nPerforming CMB analysis...\n")
  results <- perform_cmb_analysis(data, marker_var)
  
  # 4. Output results
  cat("\n=== CMB Analysis Results ===\n")
  cat("Marker Variable:", marker_var, "\n")
  cat("Mean Absolute Correlation:", round(results$mean_correlation, 3), "\n")
  
  if(nrow(results$significant_vars) > 0) {
    cat("\nSignificantly Correlated Variables (p < 0.05):\n")
    print(results$significant_vars)
  }
  
  # 5. Output conclusion
  cat("\n=== Analysis Conclusion ===\n")
  if(results$mean_correlation > 0.2) {
    cat("Warning: Data may have common method bias (mean correlation >0.2)\n")
  } else {
    cat("Common method bias level in data is acceptable (mean correlation ≤0.2)\n")
  }
  
  return(results)
}

# Run analysis
results <- main_analysis(cmb_data_2)  # Data after reverse coding for PLS procedure
results <- main_analysis(cmb_data)   # Original reverse-collected data