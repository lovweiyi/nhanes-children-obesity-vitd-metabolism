# ============================================================
# NHANES Data Network Analysis - Complete Version with Statistical Inference
# Collinearity Screening + Centrality Analysis + Statistical Inference + Stability Analysis + Multi-Center Analysis
# ============================================================

# 1. Load necessary packages ------------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(bootnet)
library(qgraph)
library(igraph)
library(ggplot2)
library(corrplot)
library(glue)
library(Matrix)
library(RColorBrewer)
library(patchwork)
library(reshape2)
library(car)
library(corrplot)
library(cluster)
library(ppcor)        # For partial correlation and p-values
library(bootnet)      # For network stability analysis
library(boot)         # For bootstrap confidence intervals

# Set working directory
setwd("~/1.NHANES/1.Public-database")

# 2. Create result directory structure ----------------------------------------------------------------
cat("Step 1: Creating result directory structure...\n")

# Main directory
if (!dir.exists("results")) dir.create("results")

# Subdirectories
subdirs <- c(
  "results/collinearity_analysis",
  "results/partial_correlations",
  "results/network_stability", 
  "results/centrality",
  "results/Multi_Center_Analysis",
  "results/network_analysis_before_selection",
  "results/network_analysis_before_selection/plots",
  "results/network_analysis_after_selection",
  "results/network_analysis_after_selection/plots"
)

for (dir in subdirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

cat("✓ Directory structure created\n")

# 3. Data loading and preprocessing ------------------------------------------------------
cat("\nStep 2: Data loading and preprocessing...\n")

# Read data
data_orignal <- read_excel("data.xlsx")
data <- data_orignal

# Remove specified column
data <- data[, -41]

# Define factor and numeric columns
factor_cols <- c(1, 2, 4, 5, 40)
all_cols <- seq_along(data)

# Convert data types
data <- data %>%
  mutate(
    across(all_of(factor_cols), as.factor),
    across(all_of(setdiff(all_cols, factor_cols)), as.numeric)
  ) %>%
  as.data.frame()

cat("✓ Data preprocessing completed\n")

# 4. Data filtering and cleaning ------------------------------------------------------
cat("Step 3: Data filtering and cleaning...\n")

# Filter 7-18 year old samples
data <- data %>%
  filter(RIDAGEYR >= 7 & RIDAGEYR <= 18)

# Remove unnecessary columns
cols_to_remove <- c("IN", "INSI", "FERSI", "RBF", "RBFSI", "FOLSI", "FOL",
                    "LBXIN", "LBDINSI", "LBDRBF", "LBXRBFSI", "LBDFOL", 
                    "LBXFOLSI", "LBDFERSI")

data <- data[, !names(data) %in% cols_to_remove]

# Rename columns
colnames(data)[12:32] <- c(
  "TSH", "HDL", "TC", "GLU", "Vit_D", "Vit_D2", "Vit_D3", 
  "VFA", "VFM", "VFV", "BMC_Head", "BMD_Head", 
  "BMD_Lower_Limb", "BMD_L_Spine", "BMD_Pelvis", 
  "Total_BMC", "Total_BMD", "Testo", "E2", "SHBG", "hs_CRP"
)

# 5. Extract continuous variables ----------------------------------------------------------------
cat("Step 4: Extracting continuous variables...\n")

continuous_vars <- data %>%
  dplyr::select(
    TSH, HDL, TC, GLU, Vit_D, Vit_D3,
    VFA, VFM, VFV, Total_BMC, Total_BMD, Testo,
    E2, SHBG, hs_CRP
  ) %>%
  mutate(across(everything(), as.numeric))

cat("✓ Continuous variables extracted\n")

# 6. Data quality check and cleaning ------------------------------------------------------------
cat("Step 5: Data quality check and cleaning...\n")

# Remove rows with too many missing values (more than 30% of variables)
threshold <- ncol(continuous_vars) * 0.3
continuous_vars_clean <- continuous_vars[rowSums(is.na(continuous_vars)) < threshold, ]

cat("Original sample size:", nrow(continuous_vars), "\n")
cat("After removing rows with >30% missing values:", nrow(continuous_vars_clean), "\n")

# Impute missing values with median
continuous_vars_imputed <- continuous_vars_clean %>%
  mutate(across(everything(), ~ifelse(is.na(.), median(., na.rm = TRUE), .)))

cat("✓ Data quality check completed\n")
cat("Final sample size:", nrow(continuous_vars_imputed), "\n")
cat("Initial number of variables:", ncol(continuous_vars_imputed), "\n")

# 7. Collinearity analysis and variable selection ---------------------------------------------------------
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("STEP 6: Collinearity Analysis and Variable Selection\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# 7.1 Calculate correlation matrix
cat("\n6.1 Calculating correlation matrix...\n")
cor_matrix <- cor(continuous_vars_imputed, use = "pairwise.complete.obs")

# Save correlation matrix
write.csv(cor_matrix, "results/collinearity_analysis/correlation_matrix.csv")

# 7.2 Identify highly correlated pairs
cat("\n6.2 Identifying highly correlated variable pairs...\n")
high_cor_threshold <- 0.8

high_cor_pairs <- which(abs(cor_matrix) > high_cor_threshold & 
                          upper.tri(cor_matrix), arr.ind = TRUE)

if (nrow(high_cor_pairs) > 0) {
  high_cor_df <- data.frame(
    Var1 = colnames(cor_matrix)[high_cor_pairs[, 1]],
    Var2 = colnames(cor_matrix)[high_cor_pairs[, 2]],
    Correlation = cor_matrix[high_cor_pairs]
  )
  write.csv(high_cor_df, "results/collinearity_analysis/high_correlation_pairs.csv", 
            row.names = FALSE)
  cat("Found", nrow(high_cor_pairs), "pairs of highly correlated variables\n")
}

# 7.3 Variable selection function
select_representative_variables <- function(data, cor_matrix) {
  
  # Define problematic variable groups
  problematic_groups <- list()
  
  # Group 1: Adiposity indicators
  adiposity_vars <- c("VFA", "VFM", "VFV")
  adiposity_in_data <- intersect(adiposity_vars, colnames(data))
  if (length(adiposity_in_data) > 1) {
    problematic_groups[["Adiposity_Indicators"]] <- adiposity_in_data
  }
  
  # Group 2: Bone indicators
  bone_vars <- c("Total_BMC", "Total_BMD")
  bone_in_data <- intersect(bone_vars, colnames(data))
  if (length(bone_in_data) > 1) {
    problematic_groups[["Bone_Indicators"]] <- bone_in_data
  }
  
  # Group 3: Vitamin D indicators
  vitd_vars <- c("Vit_D", "Vit_D3")
  vitd_in_data <- intersect(vitd_vars, colnames(data))
  if (length(vitd_in_data) > 1) {
    problematic_groups[["Vitamin_D_Indicators"]] <- vitd_in_data
  }
  
  # Selection rules
  selection_rules <- list(
    Adiposity_Indicators = list(
      representative = "VFA",
      reason = "VFA is the most clinically used visceral fat indicator"
    ),
    Bone_Indicators = list(
      representative = "Total_BMD",
      reason = "Total bone mineral density is the clinical gold standard for bone health"
    ),
    Vitamin_D_Indicators = list(
      representative = "Vit_D",
      reason = "Total vitamin D includes D2 and D3, better reflecting overall vitamin D status"
    )
  )
  
  # Apply selection rules
  all_vars <- colnames(data)
  selected_vars <- c()
  removed_vars <- c()
  removal_log <- data.frame(
    Variable = character(),
    Reason = character(),
    Group = character(),
    stringsAsFactors = FALSE
  )
  
  # First add non-problematic group variables
  non_problematic_vars <- setdiff(all_vars, unlist(problematic_groups))
  selected_vars <- c(selected_vars, non_problematic_vars)
  
  # Process each problematic group
  for (group_name in names(problematic_groups)) {
    group_vars <- problematic_groups[[group_name]]
    
    if (group_name %in% names(selection_rules)) {
      representative <- selection_rules[[group_name]]$representative
      reason <- selection_rules[[group_name]]$reason
      
      selected_vars <- c(selected_vars, representative)
      
      other_vars <- setdiff(group_vars, representative)
      removed_vars <- c(removed_vars, other_vars)
      
      for (var in other_vars) {
        removal_log <- rbind(removal_log, data.frame(
          Variable = var,
          Reason = reason,
          Group = group_name,
          stringsAsFactors = FALSE
        ))
      }
    } else {
      selected_vars <- c(selected_vars, group_vars)
    }
  }
  
  selected_vars <- unique(selected_vars)
  selected_vars <- sort(selected_vars)
  
  return(list(
    selected_variables = selected_vars,
    removed_variables = removed_vars,
    removal_log = removal_log,
    problematic_groups = problematic_groups,
    selection_rules = selection_rules
  ))
}

# 7.4 Apply variable selection
cat("\n6.3 Applying variable selection...\n")
selection_results <- select_representative_variables(continuous_vars_imputed, cor_matrix)

# Save selection results
write.csv(selection_results$removal_log, 
          "results/collinearity_analysis/variable_selection_log.csv", 
          row.names = FALSE)

cat("Variable selection results:\n")
cat("Original variables:", ncol(continuous_vars_imputed), "\n")
cat("Selected variables:", length(selection_results$selected_variables), "\n")
cat("Removed variables:", length(selection_results$removed_variables), "\n")
cat("\nSelected variables:", paste(selection_results$selected_variables, collapse = ", "), "\n")

# 7.5 Create final dataset
cat("\n6.4 Creating final dataset...\n")
continuous_vars_final <- continuous_vars_imputed[, selection_results$selected_variables]

cat("Final dataset dimensions:", dim(continuous_vars_final), "\n")

# 8. Calculate partial correlations and statistical inference ----------------------------------------------------------
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("STEP 7: Partial Correlation Analysis and Statistical Inference\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# 8.1 Calculate partial correlation matrix and p-values
cat("\n7.1 Calculating partial correlation matrix and p-values...\n")

# Use ppcor package to calculate partial correlations and p-values
tryCatch({
  pcor_result <- pcor(continuous_vars_final, method = "pearson")
  
  # Extract partial correlation matrix
  pcor_matrix <- pcor_result$estimate
  pvalue_matrix <- pcor_result$p.value
  
  cat("✓ Partial correlation calculation successful\n")
  
  # 8.2 Save results
  cat("\n7.2 Saving partial correlation results...\n")
  
  # Save partial correlation matrix
  write.csv(pcor_matrix, 
            "results/partial_correlations/partial_correlation_matrix.csv")
  
  # Save p-value matrix
  write.csv(pvalue_matrix, 
            "results/partial_correlations/pvalue_matrix.csv")
  
  # Create combined matrix with partial correlations and p-values
  combined_matrix <- matrix("", nrow = nrow(pcor_matrix), ncol = ncol(pcor_matrix))
  rownames(combined_matrix) <- rownames(pcor_matrix)
  colnames(combined_matrix) <- colnames(pcor_matrix)
  
  for (i in 1:nrow(pcor_matrix)) {
    for (j in 1:ncol(pcor_matrix)) {
      if (i != j) {
        combined_matrix[i, j] <- sprintf("%.3f (p=%.4f)", 
                                         pcor_matrix[i, j], 
                                         pvalue_matrix[i, j])
      } else {
        combined_matrix[i, j] <- "1.000"
      }
    }
  }
  
  write.csv(combined_matrix, 
            "results/partial_correlations/partial_correlation_matrix_with_pvalues.csv")
  
  # 8.3 Identify significant partial correlation edges
  cat("\n7.3 Identifying significant partial correlation edges...\n")
  
  # Set significance level
  alpha <- 0.05
  
  # Extract significant edges
  significant_edges <- which(pvalue_matrix < alpha & upper.tri(pvalue_matrix), arr.ind = TRUE)
  
  if (length(significant_edges) > 0) {
    sig_edges_df <- data.frame(
      Var1 = rownames(pcor_matrix)[significant_edges[, 1]],
      Var2 = colnames(pcor_matrix)[significant_edges[, 2]],
      Partial_Correlation = pcor_matrix[significant_edges],
      P_Value = pvalue_matrix[significant_edges],
      stringsAsFactors = FALSE
    ) %>%
      arrange(P_Value)
    
    # Add significance markers
    sig_edges_df$Significance <- ifelse(sig_edges_df$P_Value < 0.001, "***",
                                        ifelse(sig_edges_df$P_Value < 0.01, "**",
                                               ifelse(sig_edges_df$P_Value < 0.05, "*", "")))
    
    write.csv(sig_edges_df, 
              "results/partial_correlations/significant_partial_correlations.csv", 
              row.names = FALSE)
    
    cat("Found", nrow(sig_edges_df), "significant partial correlation edges (p < 0.05)\n")
    
    # 8.4 Visualize significant partial correlations
    cat("\n7.4 Creating partial correlation visualizations...\n")
    
    # Create partial correlation heatmap
    png("results/partial_correlations/partial_correlation_heatmap.png", 
        width = 1400, height = 1200, res = 150)
    
    corrplot(pcor_matrix, 
             type = "upper", 
             method = "color",
             tl.col = "black", 
             tl.srt = 45,
             title = "Partial Correlation Matrix",
             mar = c(0, 0, 3, 0))
    
    dev.off()
    
    # Create p-value heatmap
    png("results/partial_correlations/pvalue_heatmap.png", 
        width = 1400, height = 1200, res = 150)
    
    corrplot(pvalue_matrix, 
             type = "upper", 
             method = "color",
             tl.col = "black", 
             tl.srt = 45,
             title = "P-Value Matrix for Partial Correlations",
             mar = c(0, 0, 3, 0))
    
    dev.off()
    
    cat("✓ Partial correlation visualization completed\n")
  }
  
}, error = function(e) {
  cat("Partial correlation calculation failed:", e$message, "\n")
  cat("Trying alternative method...\n")
  
  # Use correlation matrix as alternative
  cor_matrix_final <- cor(continuous_vars_final, use = "pairwise.complete.obs")
  write.csv(cor_matrix_final, 
            "results/partial_correlations/correlation_matrix_final.csv")
  
  # Calculate p-values for correlations
  cor_test_results <- corr.test(continuous_vars_final, use = "pairwise")
  write.csv(cor_test_results$p, 
            "results/partial_correlations/correlation_pvalues.csv")
})

# ============================================================================
# Function: Run categorical network analysis (set line type based on significance) - Complete centrality version
# ============================================================================
run_categorical_network_analysis_complete <- function(
    analysis_data, 
    output_dir,
    analysis_name = "Network Analysis",
    gamma_value = 0.2
) {
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("Starting ", analysis_name, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(file.path(output_dir, "plots"))) {
    dir.create(file.path(output_dir, "plots"), recursive = TRUE)
  }
  
  # 1. Check correlation matrix positive definiteness
  cat("\nChecking correlation matrix positive definiteness...\n")
  
  is_positive_definite <- function(mat) {
    eigenvalues <- eigen(mat, only.values = TRUE)$values
    all(eigenvalues > 1e-10)
  }
  
  cor_matrix <- cor(analysis_data)
  
  # Check for high collinearity
  high_cor_pairs <- which(abs(cor_matrix) > 0.8 & upper.tri(cor_matrix), arr.ind = TRUE)
  if (nrow(high_cor_pairs) > 5) {
    cat("Warning: High collinearity detected, adjusting network estimation parameters\n")
    gamma_value <- max(0.05, gamma_value - 0.1)  # Reduce gamma for denser network
    cat("Adjusted gamma value:", gamma_value, "\n")
  }
  
  if (!is_positive_definite(cor_matrix)) {
    cat("Warning: Correlation matrix is not positive definite\n")
    cat("Adjusting correlation matrix to be positive definite...\n")
    
    # Use nearPD function to find nearest positive definite matrix
    near_pd <- Matrix::nearPD(cor_matrix, corr = TRUE, keepDiag = TRUE)
    cor_matrix_adj <- as.matrix(near_pd$mat)
    
    cat("Correlation matrix adjusted to be positive definite\n")
  } else {
    cor_matrix_adj <- cor_matrix
    cat("Correlation matrix is positive definite\n")
  }
  
  # 2. Network estimation - Use more flexible settings for high collinearity data
  cat("\nEstimating network...\n")
  set.seed(123)
  
  # Check data dimensions
  cat("Data dimensions:", dim(analysis_data), "\n")
  cat("Using gamma value:", gamma_value, "\n")
  
  tryCatch({
    # Use EBICglasso to estimate network
    net_all <- estimateNetwork(
      analysis_data,
      default = "EBICglasso",
      tuning = gamma_value,
      corMethod = "none",
      covMat = cor_matrix_adj,
      verbose = FALSE
    )
    cat("Network estimation successful\n")
  }, error = function(e) {
    cat("Standard EBICglasso failed:", e$message, "\n")
    
    # Try using more relaxed settings
    cat("Trying more relaxed settings for estimation...\n")
    
    tryCatch({
      # Use smaller gamma value
      net_all <- estimateNetwork(
        analysis_data,
        default = "EBICglasso",
        tuning = 0.05,  # Very small gamma value
        corMethod = "none",
        covMat = cor_matrix_adj,
        verbose = FALSE
      )
      cat("Estimation successful with relaxed settings\n")
    }, error = function(e2) {
      cat("Relaxed settings also failed:", e2$message, "\n")
      
      # Try using huge package
      cat("Trying estimation with huge package...\n")
      
      if (!require("huge")) install.packages("huge")
      library(huge)
      
      huge_result <- huge(as.matrix(analysis_data), method = "glasso")
      huge_select <- huge.select(huge_result, criterion = "stars")
      
      net_matrix <- as.matrix(huge_select$opt.icov)
      colnames(net_matrix) <- rownames(net_matrix) <- colnames(analysis_data)
      
      net_all <- list(graph = net_matrix, results = list(opt.icov = net_matrix))
      cat("Estimation successful with huge package\n")
    })
  })
  
  # 3. Safely extract and process network matrix
  cat("\nExtracting and processing network matrix...\n")
  
  # Safely extract network matrix
  safe_extract_network <- function(net_obj) {
    if (is.list(net_obj)) {
      if ("graph" %in% names(net_obj)) {
        return(net_obj$graph)
      } else if ("results" %in% names(net_obj)) {
        if ("opt.icov" %in% names(net_obj$results)) {
          return(net_obj$results$opt.icov)
        }
      }
    }
    return(net_obj)
  }
  
  network_matrix_raw <- safe_extract_network(net_all)
  
  # Ensure it's a matrix
  if (!is.matrix(network_matrix_raw)) {
    network_matrix_raw <- as.matrix(network_matrix_raw)
  }
  
  cat("Raw network matrix dimensions:", dim(network_matrix_raw), "\n")
  
  # Get analysis data variables
  data_vars <- colnames(analysis_data)
  n_data_vars <- length(data_vars)
  
  # Check and handle dimensions
  if (ncol(network_matrix_raw) == n_data_vars) {
    # Dimensions match, use directly
    network_matrix <- network_matrix_raw
    rownames(network_matrix) <- colnames(network_matrix) <- data_vars
    cat("Dimensions match, variable names set\n")
  } else {
    cat("Warning: Dimension mismatch! Raw dimensions:", dim(network_matrix_raw), ", Expected:", n_data_vars, "\n")
    
    # Create new empty matrix
    network_matrix <- matrix(0, nrow = n_data_vars, ncol = n_data_vars)
    rownames(network_matrix) <- colnames(network_matrix) <- data_vars
    
    # Try to match existing names
    if (!is.null(colnames(network_matrix_raw))) {
      common_vars <- intersect(colnames(network_matrix_raw), data_vars)
      if (length(common_vars) > 0) {
        cat("Found", length(common_vars), "common variables\n")
        
        # Get indices
        raw_indices <- match(common_vars, colnames(network_matrix_raw))
        new_indices <- match(common_vars, data_vars)
        
        # Copy values
        for (i in 1:length(common_vars)) {
          for (j in 1:length(common_vars)) {
            network_matrix[new_indices[i], new_indices[j]] <- 
              network_matrix_raw[raw_indices[i], raw_indices[j]]
          }
        }
      }
    } else {
      # No names, assume first n variables match
      n_common <- min(ncol(network_matrix_raw), n_data_vars)
      cat("No variable names, using first", n_common, "variables\n")
      
      network_matrix[1:n_common, 1:n_common] <- 
        network_matrix_raw[1:n_common, 1:n_common]
    }
  }
  
  cat("Final network matrix dimensions:", dim(network_matrix), "\n")
  cat("Network nodes:", ncol(network_matrix), "\n")
  
  # 4. Calculate p-values for network partial correlations
  cat("\nCalculating p-values for network partial correlations...\n")
  
  calculate_partial_correlation_pvalues <- function(partial_corr_matrix, n_samples) {
    p <- ncol(partial_corr_matrix)
    p_value_matrix <- matrix(NA, nrow = p, ncol = p)
    
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        r <- partial_corr_matrix[i, j]
        if (!is.na(r) && r != 0) {
          # Calculate t-statistic (simplified method)
          t_stat <- r * sqrt((n_samples - p) / (1 - r^2))
          df <- n_samples - p
          
          # Calculate two-sided p-value
          p_value <- 2 * pt(-abs(t_stat), df)
          p_value_matrix[i, j] <- p_value
          p_value_matrix[j, i] <- p_value
        }
      }
    }
    
    diag(p_value_matrix) <- 0
    colnames(p_value_matrix) <- rownames(p_value_matrix) <- colnames(partial_corr_matrix)
    
    return(p_value_matrix)
  }
  
  n_samples <- nrow(analysis_data)
  network_p_values <- calculate_partial_correlation_pvalues(network_matrix, n_samples)
  
  # Apply FDR correction
  network_p_values_fdr <- matrix(p.adjust(as.vector(network_p_values), method = "fdr"),
                                 nrow = nrow(network_p_values),
                                 ncol = ncol(network_p_values))
  colnames(network_p_values_fdr) <- rownames(network_p_values_fdr) <- colnames(network_matrix)
  
  # 5. Define variable categories and colors
  cat("\nDefining variable categories and colors...\n")
  
  # Create variable category mapping
  variable_categories <- list(
    "Metabolic Indicators" = c("TSH", "HDL", "TC", "GLU", "hs_CRP"),
    "Vitamin D" = c("Vit_D", "Vit_D3"),
    "Adiposity Indicators" = c("VFA", "VFM", "VFV"),
    "Bone Metabolism" = c("Total_BMC", "Total_BMD"),
    "Sex Hormones" = c("Testo", "E2", "SHBG")
  )
  
  # Define color scheme (Set2 palette, suitable for categorical variables)
  category_colors <- c(
    "Metabolic Indicators" = "#66C2A5",    # Green
    "Vitamin D" = "#FC8D62",               # Orange
    "Adiposity Indicators" = "#8DA0CB",    # Blue
    "Bone Metabolism" = "#E78AC3",         # Pink
    "Sex Hormones" = "#A6D854",            # Light green
    "Other" = "#B3B3B3"                    # Gray
  )
  
  # Assign category and color to each variable
  assign_variable_category <- function(var_name) {
    for(category in names(variable_categories)) {
      if(var_name %in% variable_categories[[category]]) {
        return(category)
      }
    }
    return("Other")
  }
  
  # Get all variable names
  all_variables <- colnames(network_matrix)
  
  # Create category and color vectors
  variable_categories_vec <- sapply(all_variables, assign_variable_category)
  node_colors_categorical <- sapply(variable_categories_vec, function(cat) category_colors[cat])
  
  # Save category information
  category_info <- data.frame(
    Variable = all_variables,
    Category = variable_categories_vec,
    Color = node_colors_categorical,
    stringsAsFactors = FALSE
  )
  
  write.csv(category_info, file.path(output_dir, "variable_categories.csv"), row.names = FALSE)
  
  cat("Variable category information:\n")
  print(category_info)
  
  # 6. Create line type matrix (based on significance)
  cat("\nCreating line type matrix (based on significance)...\n")
  
  # Create line type matrix: significant connections (p<0.05) use solid line(1), non-significant connections use dashed line(2)
  lty_matrix <- matrix(1, nrow = nrow(network_matrix), ncol = ncol(network_matrix))
  
  # Set non-significant connections (p>=0.05) to dashed line
  lty_matrix[network_p_values >= 0.05 & network_p_values > 0] <- 2
  
  # Set diagonal and zero-weight edges to 0 (not displayed)
  diag(lty_matrix) <- 0
  lty_matrix[network_matrix == 0] <- 0
  
  cat("Line type matrix created:\n")
  cat("  Solid lines (p<0.05):", sum(lty_matrix == 1, na.rm = TRUE)/2, "edges\n")
  cat("  Dashed lines (p>=0.05):", sum(lty_matrix == 2, na.rm = TRUE)/2, "edges\n")
  
  # 7. Calculate network statistics (complete centrality metrics)
  cat("\nCalculating network statistics (complete centrality metrics)...\n")
  
  # Extract network matrix (remove diagonal)
  network_matrix_no_diag <- network_matrix
  diag(network_matrix_no_diag) <- 0
  
  # Calculate network density - using non-zero edges
  total_possible_edges <- (ncol(network_matrix_no_diag)^2 - ncol(network_matrix_no_diag)) / 2
  existing_edges <- sum(abs(network_matrix_no_diag) > 0) / 2  # Only count non-zero connections
  network_density <- existing_edges / total_possible_edges
  
  # Count edge strength distribution
  edge_strengths <- abs(network_matrix_no_diag[upper.tri(network_matrix_no_diag)])
  edge_strengths <- edge_strengths[edge_strengths > 0]  # Only consider non-zero edges
  
  # Calculate descriptive statistics
  edge_summary <- if(length(edge_strengths) > 0) {
    data.frame(
      Min = min(edge_strengths),
      Q1 = quantile(edge_strengths, 0.25),
      Median = median(edge_strengths),
      Mean = mean(edge_strengths),
      Q3 = quantile(edge_strengths, 0.75),
      Max = max(edge_strengths),
      SD = sd(edge_strengths)
    )
  } else {
    data.frame(
      Min = 0, Q1 = 0, Median = 0, Mean = 0, Q3 = 0, Max = 0, SD = 0
    )
  }
  
  # Calculate complete centrality metrics
  calculate_complete_centrality <- function(net_matrix) {
    # Ensure diagonal is 0
    diag(net_matrix) <- 0
    
    # 1. Strength centrality (sum of absolute weights)
    strength <- colSums(abs(net_matrix))
    
    # 2. Degree centrality (number of non-zero connections)
    degree <- colSums(abs(net_matrix) > 0)
    
    # 3. Expected influence (weighted sum, considering positive and negative)
    expected_influence <- colSums(net_matrix)
    
    # Initialize other centrality metrics
    betweenness <- rep(NA, ncol(net_matrix))
    closeness <- rep(NA, ncol(net_matrix))
    eigenvector <- rep(NA, ncol(net_matrix))
    
    # Try to calculate other centrality metrics using igraph
    if (requireNamespace("igraph", quietly = TRUE)) {
      tryCatch({
        # Create undirected weighted graph
        g <- igraph::graph_from_adjacency_matrix(
          net_matrix,
          mode = "undirected",
          weighted = TRUE,
          diag = FALSE
        )
        
        # Calculate betweenness centrality
        betweenness <- igraph::betweenness(g, weights = NA, normalized = FALSE)
        
        # Calculate closeness centrality
        closeness <- igraph::closeness(g, weights = NA, normalized = FALSE)
        
        # Calculate eigenvector centrality (handling negative weights)
        eigen_result <- tryCatch({
          igraph::eigen_centrality(g, weights = NULL)
        }, error = function(e) {
          # If fails, use absolute value weights
          abs_g <- igraph::graph_from_adjacency_matrix(
            abs(net_matrix),
            mode = "undirected",
            weighted = TRUE,
            diag = FALSE
          )
          igraph::eigen_centrality(abs_g, weights = NULL)
        })
        eigenvector <- eigen_result$vector
        
        cat("igraph centrality calculation successful\n")
      }, error = function(e) {
        cat("igraph centrality calculation warning:", e$message, "\n")
      })
    }
    
    # Return results
    return(list(
      strength = strength,
      degree = degree,
      expected_influence = expected_influence,
      betweenness = betweenness,
      closeness = closeness,
      eigenvector = eigenvector
    ))
  }
  
  # Calculate centrality
  centrality_results <- calculate_complete_centrality(network_matrix)
  
  # Create complete centrality data frame
  centrality_df <- data.frame(
    Variable = colnames(network_matrix),
    Strength = round(centrality_results$strength, 4),
    Degree = centrality_results$degree,
    ExpectedInfluence = round(centrality_results$expected_influence, 4),
    Betweenness = round(centrality_results$betweenness, 4),
    Closeness = round(centrality_results$closeness, 4),
    Eigenvector = round(centrality_results$eigenvector, 4),
    stringsAsFactors = FALSE
  )
  
  # Add category information
  centrality_df$Category <- variable_categories_vec
  
  # Sort by strength centrality
  centrality_df <- centrality_df[order(centrality_df$Strength, decreasing = TRUE), ]
  
  # Print top 10 nodes
  cat("\nCentrality ranking (top 10):\n")
  print(head(centrality_df, 10))
  
  # 8. Draw categorical colored network graph (set line type based on significance)
  cat("\nDrawing categorical colored network graph (set line type based on significance)...\n")
  
  # Check if network matrix contains too many zero values
  zero_count <- sum(network_matrix == 0) - ncol(network_matrix)  # Subtract diagonal
  total_possible <- ncol(network_matrix) * (ncol(network_matrix) - 1)
  zero_ratio <- zero_count / total_possible
  
  # 8.1 Main network graph - Show all non-zero connections, set line type based on significance
  png(file.path(output_dir, "plots", "network_categorical_all_connections.png"), 
      width = 1600, height = 1400, res = 150)
  
  if (zero_ratio > 0.9 && nrow(high_cor_pairs) > 5) {  
    # If over 90% of connections are zero and high collinearity exists, create information graph
    cat("Warning: Network is very sparse (zero connection ratio:", round(zero_ratio*100, 1), "%), creating information graph\n")
    
    plot.new()
    text(0.5, 0.8, paste(analysis_name, "\nNetwork Structure Analysis", sep = ""), cex = 1.8, font = 2)
    text(0.5, 0.6, paste("Variables:", ncol(network_matrix), 
                         "\nSample size:", n_samples,
                         "\nNetwork density:", round(network_density, 4),
                         "\nZero connection ratio:", round(zero_ratio*100, 1), "%",
                         "\nHigh collinearity detected:", nrow(high_cor_pairs), "pairs"), 
         cex = 1.4)
    text(0.5, 0.3, "Network is very sparse, likely due to high collinearity\nPlease check network analysis after collinearity screening", 
         cex = 1.4, col = "red")
    
    # Add variable category legend
    legend("bottom", 
           legend = unique(variable_categories_vec),
           fill = category_colors[match(unique(variable_categories_vec), names(category_colors))],
           title = "Variable Categories",
           horiz = FALSE,
           cex = 1.0,
           ncol = 2)
  } else {
    tryCatch({
      # Draw network graph - Show all non-zero connections, set line type based on significance
      edge_width <- ifelse(abs(network_matrix) > 0, 2.5, 0)
      diag(edge_width) <- 0
      
      # Adjust margins for better display
      par(mar = c(1, 1, 3, 1))
      
      qgraph(network_matrix,
             layout = "spring",
             title = paste(analysis_name, "\nSolid: p<0.05 | Dashed: p≥0.05 | Sample size:", n_samples, sep = ""),
             labels = colnames(network_matrix),
             color = node_colors_categorical,
             vsize = 9,
             edge.width = edge_width,
             lty = lty_matrix,
             posCol = "#0D47A1",
             negCol = "#B71C1C",
             fade = FALSE,
             label.cex = 1.1,
             title.cex = 1.4,
             edge.labels = FALSE,
             minimum = 0,
             repulsion = 0.8,
             mar = c(5, 5, 5, 5))
      
      # Add legend explanation
      legend("bottomleft", 
             legend = c("Positive correlation (p<0.05)", "Negative correlation (p<0.05)", 
                        "Positive correlation (p≥0.05)", "Negative correlation (p≥0.05)"),
             lty = c(1, 1, 2, 2),
             lwd = 2.5,
             col = c("#0D47A1", "#B71C1C", "#0D47A1", "#B71C1C"),
             bty = "n",
             cex = 0.9,
             title = "Edge Types")
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Network graph drawing failed:", e$message), cex = 1.5)
    })
  }
  
  dev.off()
  cat("Categorical colored network graph saved: ", 
      file.path(output_dir, "plots", "network_categorical_all_connections.png"), "\n")
  
  # 8.2 Draw weighted network graph (edge width adjusted based on weight, line type based on significance)
  cat("\nDrawing weighted network graph (edge width adjusted based on weight, line type based on significance)...\n")
  
  png(file.path(output_dir, "plots", "network_weighted.png"), 
      width = 1600, height = 1400, res = 150)
  
  if (zero_ratio > 0.9 && nrow(high_cor_pairs) > 5) {
    # Create correlation matrix heatmap as alternative
    cor_matrix_plot <- cor(analysis_data)
    
    # Adjust margins for better display
    par(mar = c(0, 0, 4, 0))
    
    corrplot(cor_matrix_plot,
             method = "color",
             type = "upper",
             tl.col = "black",
             tl.srt = 45,
             title = paste(analysis_name, "\nCorrelation Matrix Heatmap (Alternative to Sparse Network)"),
             mar = c(0, 0, 3, 0))
  } else {
    # Calculate weighted edge width: base width 1 + absolute weight * 5
    edge_weighted_width <- 1 + abs(network_matrix) * 5
    diag(edge_weighted_width) <- 0
    edge_weighted_width[network_matrix == 0] <- 0
    
    # Adjust margins for better display
    par(mar = c(1, 1, 3, 1))
    
    qgraph(network_matrix,
           layout = "spring",
           title = paste(analysis_name, "- Weighted Network Graph\n",
                         "Edge width adjusted by connection strength | Solid: p<0.05 | Dashed: p≥0.05\nSample size: ", n_samples, sep = ""),
           labels = colnames(network_matrix),
           color = node_colors_categorical,
           vsize = 9,
           edge.width = edge_weighted_width,  # Weighted edge width
           lty = lty_matrix,  # Line type based on significance
           posCol = "#0D47A1",
           negCol = "#B71C1C",
           fade = FALSE,
           label.cex = 1.1,
           title.cex = 1.4,
           edge.labels = FALSE,
           minimum = 0,  # Show all non-zero connections
           repulsion = 0.8,
           mar = c(5, 5, 5, 5))
    
    # Add legend
    legend("bottomright", 
           legend = c("Solid: p<0.05", "Dashed: p≥0.05"),
           lty = c(1, 2),
           lwd = 2.5,
           col = "gray50",
           bty = "n",
           cex = 0.9,
           title = "Significance")
  }
  
  dev.off()
  cat("Weighted network graph saved\n")
  
  # 8.3 Draw simplified network graph showing only significant connections
  cat("\nDrawing simplified network graph showing only significant connections...\n")
  
  png(file.path(output_dir, "plots", "network_significant_only.png"), 
      width = 1600, height = 1400, res = 150)
  
  # Create matrix containing only significant connections
  network_significant <- network_matrix
  network_significant[network_p_values >= 0.05] <- 0
  diag(network_significant) <- 0
  
  # Calculate number of significant connections, ensure no NA values
  significant_edges_count <- sum(network_p_values < 0.05 & upper.tri(network_p_values), na.rm = TRUE)
  
  # Check if significant_edges_count is valid
  if (!is.na(significant_edges_count) && significant_edges_count > 0) {
    # Only draw graph with significant connections
    # Adjust margins for better display
    par(mar = c(1, 1, 3, 1))
    
    qgraph(network_significant,
           layout = "spring",
           title = paste(analysis_name, "- Significant Connections Network Graph\n",
                         "Only showing connections with p<0.05 | Significant connections: ", significant_edges_count, 
                         " | Sample size: ", n_samples, sep = ""),
           labels = colnames(network_matrix),
           color = node_colors_categorical,
           vsize = 9,
           edge.width = 2.5,
           posCol = "#0D47A1",
           negCol = "#B71C1C",
           fade = FALSE,
           label.cex = 1.1,
           title.cex = 1.4,
           edge.labels = FALSE,
           minimum = 0,  # Show all non-zero connections
           repulsion = 0.8,
           mar = c(5, 5, 5, 5))
  } else {
    plot.new()
    text(0.5, 0.5, "No significant connections found (p<0.05)", cex = 1.8)
  }
  
  dev.off()
  cat("Significant connections network graph saved\n")
  
  # 8.4 Draw variable category legend (saved separately)
  cat("\nDrawing variable category legend...\n")
  
  png(file.path(output_dir, "plots", "variable_categories_legend.png"), 
      width = 1000, height = 800, res = 150)
  
  # Create legend data
  legend_data <- data.frame(
    Category = unique(variable_categories_vec),
    Color = sapply(unique(variable_categories_vec), function(cat) category_colors[cat]),
    Count = sapply(unique(variable_categories_vec), function(cat) sum(variable_categories_vec == cat))
  )
  
  # Sort by category name
  legend_data <- legend_data[order(legend_data$Category), ]
  
  # Draw simple legend with adjusted layout
  par(mar = c(2, 2, 4, 2))
  plot(0, 0, type = "n", 
       xlim = c(0, 1), ylim = c(0, nrow(legend_data) + 1.5),
       axes = FALSE, xlab = "", ylab = "", main = "Variable Categories")
  
  # Add legend items
  for (i in 1:nrow(legend_data)) {
    y_pos <- nrow(legend_data) - i + 1
    
    # Draw color square
    rect(0.15, y_pos - 0.4, 0.25, y_pos + 0.4, 
         col = legend_data$Color[i], border = "black", lwd = 1)
    
    # Draw category name
    text(0.3, y_pos, legend_data$Category[i], pos = 4, cex = 1.3)
    
    # Draw variable count
    text(0.85, y_pos, paste("(", legend_data$Count[i], " variables)", sep = ""), pos = 4, cex = 1.1)
  }
  
  # Add network statistics information
  text(0.5, nrow(legend_data) + 1.2, 
       paste("Network Statistics: ", ncol(network_matrix), " variables, ", 
             existing_edges, " edges, Density=", round(network_density, 3), sep = ""), 
       cex = 1.2, font = 2)
  
  dev.off()
  cat("Variable category legend saved\n")
  
  # 9. Save network analysis results
  cat("\nSaving network analysis results...\n")
  
  # Save network matrix
  write.csv(network_matrix, file.path(output_dir, "network_matrix.csv"), row.names = TRUE)
  write.csv(network_p_values, file.path(output_dir, "network_p_values.csv"), row.names = TRUE)
  write.csv(network_p_values_fdr, file.path(output_dir, "network_p_values_fdr.csv"), row.names = TRUE)
  write.csv(lty_matrix, file.path(output_dir, "line_type_matrix.csv"), row.names = TRUE)
  write.csv(centrality_df, file.path(output_dir, "centrality_measures.csv"), row.names = FALSE)
  
  # Create network edge list
  edge_list <- data.frame()
  
  for (i in 1:(ncol(network_matrix)-1)) {
    for (j in (i+1):ncol(network_matrix)) {
      if (network_matrix[i, j] != 0) {  # Save all non-zero connections
        edge_list <- rbind(edge_list, data.frame(
          Node1 = colnames(network_matrix)[i],
          Category1 = variable_categories_vec[i],
          Node2 = colnames(network_matrix)[j],
          Category2 = variable_categories_vec[j],
          Partial_Correlation = round(network_matrix[i, j], 4),
          Abs_Correlation = round(abs(network_matrix[i, j]), 4),
          P_value = round(network_p_values[i, j], 6),
          P_value_FDR = round(network_p_values_fdr[i, j], 6),
          Significant = ifelse(network_p_values[i, j] < 0.05, "Yes", "No"),
          Line_Type = ifelse(network_p_values[i, j] < 0.05, "Solid", "Dashed"),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Sort by absolute partial correlation coefficient
  edge_list <- edge_list[order(edge_list$Abs_Correlation, decreasing = TRUE), ]
  write.csv(edge_list, file.path(output_dir, "edge_list.csv"), row.names = FALSE)
  
  # Save edge strength statistics
  write.csv(edge_summary, file.path(output_dir, "edge_strength_summary.csv"), row.names = FALSE)
  
  # Save line type statistics
  line_type_summary <- data.frame(
    Line_Type = c("Solid(p<0.05)", "Dashed(p≥0.05)"),
    Count = c(sum(lty_matrix == 1, na.rm = TRUE)/2, sum(lty_matrix == 2, na.rm = TRUE)/2)
  )
  write.csv(line_type_summary, file.path(output_dir, "line_type_summary.csv"), row.names = FALSE)
  
  # Save network statistics summary
  sink(file.path(output_dir, "network_statistics.txt"))
  cat("=== Network Analysis Statistical Summary (Set Line Type Based on Significance) ===\n\n")
  cat("Analysis Name:", analysis_name, "\n")
  cat("Analysis Time:", Sys.time(), "\n")
  cat("Sample Size:", n_samples, "\n")
  cat("Number of Variables:", ncol(network_matrix), "\n\n")
  
  cat("Network Settings:\n")
  cat("Network Estimation Method: EBICglasso\n")
  cat("Regularization Parameter (gamma):", gamma_value, "\n")
  cat("Display Settings: Show all non-zero connections\n")
  cat("Line Type Settings: Significant connections (p<0.05) use solid lines, non-significant connections (p≥0.05) use dashed lines\n")
  cat("Edge Width Settings: Edge width adjusted based on connection strength\n\n")
  
  cat("Variable Category Statistics:\n")
  for(category in unique(variable_categories_vec)) {
    var_count <- sum(variable_categories_vec == category)
    cat(category, ": ", var_count, " variables (", 
        round(var_count/length(variable_categories_vec)*100, 1), "%)\n", sep = "")
  }
  cat("\n")
  
  cat("Network Connection Statistics:\n")
  cat("Total Possible Edges:", total_possible_edges, "\n")
  cat("Non-Zero Edges:", existing_edges, "\n")
  cat("Network Density:", round(network_density, 3), "\n")
  cat("Sparsity:", round(1 - network_density, 3), "\n\n")
  
  cat("Edge Significance Statistics:\n")
  significant_edges <- sum(network_p_values < 0.05 & upper.tri(network_p_values), na.rm = TRUE)
  non_significant_edges <- existing_edges - significant_edges
  cat("Significant Edges (p < 0.05):", significant_edges, "\n")
  cat("Non-Significant Edges (p ≥ 0.05):", non_significant_edges, "\n")
  cat("Significant Edge Proportion:", round(significant_edges/existing_edges*100, 1), "%\n\n")
  
  cat("Line Type Statistics:\n")
  cat("Solid Edges (p<0.05):", sum(lty_matrix == 1, na.rm = TRUE)/2, "\n")
  cat("Dashed Edges (p≥0.05):", sum(lty_matrix == 2, na.rm = TRUE)/2, "\n\n")
  
  cat("Edge Strength Statistics (Non-Zero Connections):\n")
  if (length(edge_strengths) > 0) {
    cat("Number of Edges:", length(edge_strengths), "\n")
    cat("Minimum:", round(min(edge_strengths), 4), "\n")
    cat("Maximum:", round(max(edge_strengths), 4), "\n")
    cat("Mean:", round(mean(edge_strengths), 4), "\n")
    cat("Median:", round(median(edge_strengths), 4), "\n")
    cat("Standard Deviation:", round(sd(edge_strengths), 4), "\n\n")
  } else {
    cat("No non-zero connections\n\n")
  }
  
  cat("Centrality Ranking (Top 5):\n")
  top5_strength <- head(centrality_df, 5)
  for (k in 1:nrow(top5_strength)) {
    cat(k, ". ", top5_strength$Variable[k], " (", top5_strength$Category[k], ")", 
        ": Strength=", top5_strength$Strength[k],
        ", Degree=", top5_strength$Degree[k],
        ", Expected Influence=", top5_strength$ExpectedInfluence[k], "\n", sep = "")
  }
  
  cat("\nStrongest Connections (Top 10):\n")
  if (nrow(edge_list) > 0) {
    top10_edges <- head(edge_list, 10)
    for (k in 1:nrow(top10_edges)) {
      cat(k, ". ", top10_edges$Node1[k], " (", top10_edges$Category1[k], ") - ", 
          top10_edges$Node2[k], " (", top10_edges$Category2[k], ")", 
          ": r = ", round(top10_edges$Partial_Correlation[k], 3), 
          ", p = ", format(top10_edges$P_value[k], scientific = TRUE, digits = 3),
          " (", ifelse(top10_edges$Significant[k] == "Yes", "Significant", "Non-significant"), 
          ", ", top10_edges$Line_Type[k], ")\n", sep = "")
    }
  } else {
    cat("No connections found\n")
  }
  
  # Count connections by category
  cat("\nConnection Statistics by Category:\n")
  if (nrow(edge_list) > 0) {
    category_pairs <- data.frame()
    for (i in 1:length(unique(variable_categories_vec))) {
      for (j in i:length(unique(variable_categories_vec))) {
        cat1 <- unique(variable_categories_vec)[i]
        cat2 <- unique(variable_categories_vec)[j]
        
        # Count intra-category connections
        if (i == j) {
          # Same category connections
          cat_vars <- names(variable_categories_vec[variable_categories_vec == cat1])
          if (length(cat_vars) > 1) {
            possible_edges <- (length(cat_vars)^2 - length(cat_vars)) / 2
            actual_edges <- sum(edge_list$Category1 == cat1 & edge_list$Category2 == cat1)
            density <- round(actual_edges / possible_edges, 3)
            
            category_pairs <- rbind(category_pairs, data.frame(
              Category1 = cat1,
              Category2 = cat1,
              Possible_Edges = possible_edges,
              Actual_Edges = actual_edges,
              Density = density,
              stringsAsFactors = FALSE
            ))
          }
        } else {
          # Different category connections
          cat1_vars <- names(variable_categories_vec[variable_categories_vec == cat1])
          cat2_vars <- names(variable_categories_vec[variable_categories_vec == cat2])
          
          possible_edges <- length(cat1_vars) * length(cat2_vars)
          actual_edges <- sum((edge_list$Category1 == cat1 & edge_list$Category2 == cat2) |
                                (edge_list$Category1 == cat2 & edge_list$Category2 == cat1))
          density <- round(actual_edges / possible_edges, 3)
          
          category_pairs <- rbind(category_pairs, data.frame(
            Category1 = cat1,
            Category2 = cat2,
            Possible_Edges = possible_edges,
            Actual_Edges = actual_edges,
            Density = density,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    # Sort by density
    category_pairs <- category_pairs[order(category_pairs$Density, decreasing = TRUE), ]
    
    for (k in 1:nrow(category_pairs)) {
      if (category_pairs$Category1[k] == category_pairs$Category2[k]) {
        cat(category_pairs$Category1[k], " intra-category connection density: ", category_pairs$Density[k], 
            " (", category_pairs$Actual_Edges[k], "/", category_pairs$Possible_Edges[k], ")\n", sep = "")
      } else {
        cat(category_pairs$Category1[k], " - ", category_pairs$Category2[k], 
            " connection density: ", category_pairs$Density[k], 
            " (", category_pairs$Actual_Edges[k], "/", category_pairs$Possible_Edges[k], ")\n", sep = "")
      }
    }
  } else {
    cat("No connections to count\n")
  }
  
  sink()
  
  # Save inter-category connection statistics
  if (exists("category_pairs") && nrow(category_pairs) > 0) {
    write.csv(category_pairs, file.path(output_dir, "category_connection_density.csv"), row.names = FALSE)
  }
  
  cat("Network analysis completed! Results saved to ", output_dir, "\n")
  
  # Return network matrix and centrality results for subsequent use
  return(list(
    network_matrix = network_matrix,
    centrality_df = centrality_df,
    edge_list = edge_list,
    variable_categories = variable_categories_vec,
    n_samples = n_samples,
    network_p_values = network_p_values,
    high_cor_pairs = high_cor_pairs
  ))
}

# ============================================================================
# 9.5 Supplementary Analysis: Pre-collinearity screening correlation matrix visualization
# ============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Supplementary Analysis: Pre-Collinearity Screening Correlation Matrix Visualization\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

# Create collinearity visualization directory
collinearity_vis_dir <- "results/network_analysis_before_selection/collinearity_visualization"
if (!dir.exists(collinearity_vis_dir)) {
  dir.create(collinearity_vis_dir, recursive = TRUE)
}

# Calculate correlation matrix
cor_matrix_before <- cor(continuous_vars_imputed)

# Save correlation matrix
write.csv(cor_matrix_before, 
          file.path(collinearity_vis_dir, "correlation_matrix_before_selection.csv"))

# Create correlation matrix heatmap
png(file.path(collinearity_vis_dir, "correlation_heatmap_before_selection.png"), 
    width = 1400, height = 1200, res = 150)

# Adjust margins for better display
par(mar = c(0, 0, 4, 0))

corrplot(cor_matrix_before,
         method = "color",
         type = "upper",
         tl.col = "black",
         tl.srt = 45,
         title = "Correlation Matrix Before Collinearity Screening",
         mar = c(0, 0, 3, 0),
         diag = FALSE)

dev.off()

# Create collinearity network graph (based on high correlation threshold)
cat("\nCreating high correlation threshold-based network graph...\n")

# Set high correlation threshold
high_cor_threshold_vis <- 0.7

# Create high correlation adjacency matrix
high_cor_matrix <- ifelse(abs(cor_matrix_before) > high_cor_threshold_vis, 
                          cor_matrix_before, 0)
diag(high_cor_matrix) <- 0

# Count high correlation connections
high_cor_edges <- sum(high_cor_matrix != 0) / 2
cat("High correlation connections (r >", high_cor_threshold_vis, "):", high_cor_edges, "\n")

# Draw high correlation network graph
png(file.path(collinearity_vis_dir, "high_correlation_network.png"), 
    width = 1600, height = 1400, res = 150)

# Define variable categories and colors (consistent with main analysis)
variable_categories_before <- list(
  "Metabolic Indicators" = c("TSH", "HDL", "TC", "GLU", "hs_CRP"),
  "Vitamin D" = c("Vit_D", "Vit_D3"),
  "Adiposity Indicators" = c("VFA", "VFM", "VFV"),
  "Bone Metabolism" = c("Total_BMC", "Total_BMD"),
  "Sex Hormones" = c("Testo", "E2", "SHBG")
)

# Assign category and color to each variable
assign_category_before <- function(var_name) {
  for(category in names(variable_categories_before)) {
    if(var_name %in% variable_categories_before[[category]]) {
      return(category)
    }
  }
  return("Other")
}

all_vars_before <- colnames(continuous_vars_imputed)
variable_categories_vec_before <- sapply(all_vars_before, assign_category_before)

# Use same color scheme
category_colors <- c(
  "Metabolic Indicators" = "#66C2A5",    # Green
  "Vitamin D" = "#FC8D62",               # Orange
  "Adiposity Indicators" = "#8DA0CB",    # Blue
  "Bone Metabolism" = "#E78AC3",         # Pink
  "Sex Hormones" = "#A6D854",            # Light green
  "Other" = "#B3B3B3"                    # Gray
)

node_colors_before <- sapply(variable_categories_vec_before, function(cat) {
  if (cat %in% names(category_colors)) {
    return(category_colors[cat])
  } else {
    return(category_colors["Other"])
  }
})

# Adjust margins for better display
par(mar = c(1, 1, 3, 1))

# Draw network graph
qgraph(high_cor_matrix,
       layout = "spring",
       title = paste("High Correlation Network Before Collinearity Screening\n",
                     "Threshold: |r| >", high_cor_threshold_vis, 
                     " | High Correlation Connections:", high_cor_edges, sep = ""),
       labels = colnames(high_cor_matrix),
       color = node_colors_before,
       vsize = 9,
       posCol = "#0D47A1",
       negCol = "#B71C1C",
       fade = FALSE,
       label.cex = 1.1,
       title.cex = 1.4,
       edge.labels = FALSE,
       minimum = high_cor_threshold_vis,
       maximum = 1,
       mar = c(5, 5, 5, 5))

dev.off()

# Create intra-group correlation analysis
cat("\nAnalyzing intra-group correlations...\n")

group_cor_analysis <- data.frame()

# Analyze each variable group
for(group_name in names(variable_categories_before)) {
  group_vars <- variable_categories_before[[group_name]]
  
  if(length(group_vars) > 1) {
    # Extract intra-group correlation matrix
    group_cor <- cor_matrix_before[group_vars, group_vars]
    
    # Calculate intra-group average correlation
    group_mean_cor <- mean(abs(group_cor[upper.tri(group_cor)]))
    group_max_cor <- max(abs(group_cor[upper.tri(group_cor)]))
    
    # Count high correlation pairs
    high_cor_pairs_group <- sum(abs(group_cor[upper.tri(group_cor)]) > high_cor_threshold_vis)
    total_pairs_group <- length(group_cor[upper.tri(group_cor)])
    
    group_cor_analysis <- rbind(group_cor_analysis, data.frame(
      Group = group_name,
      Variables = paste(group_vars, collapse = ", "),
      N_Variables = length(group_vars),
      Mean_Correlation = round(group_mean_cor, 3),
      Max_Correlation = round(group_max_cor, 3),
      High_Cor_Pairs = high_cor_pairs_group,
      Total_Pairs = total_pairs_group,
      High_Cor_Ratio = round(high_cor_pairs_group/total_pairs_group, 3),
      stringsAsFactors = FALSE
    ))
  }
}

# Save intra-group correlation analysis
write.csv(group_cor_analysis, 
          file.path(collinearity_vis_dir, "within_group_correlation_analysis.csv"),
          row.names = FALSE)

cat("\nIntra-group correlation analysis:\n")
print(group_cor_analysis)

# Create intra-group correlation visualization
if(nrow(group_cor_analysis) > 0) {
  png(file.path(collinearity_vis_dir, "within_group_correlation_plot.png"), 
      width = 1400, height = 900, res = 150)
  
  # Prepare plot data
  plot_data <- group_cor_analysis
  plot_data$Group <- factor(plot_data$Group, levels = plot_data$Group[order(plot_data$Mean_Correlation, decreasing = TRUE)])
  
  # Create bar plot
  bar_plot <- ggplot(plot_data, aes(x = Group, y = Mean_Correlation, fill = Group)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = paste0(Mean_Correlation, "\n(n=", N_Variables, ")")), 
              vjust = -0.5, size = 5) +
    labs(
      title = "Intra-Group Mean Correlation Before Collinearity Screening",
      x = "Variable Group",
      y = "Mean Absolute Correlation Coefficient",
      caption = paste("High correlation proportion (based on", high_cor_threshold_vis, "threshold) shown in labels")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.caption = element_text(size = 10, color = "gray40"),
      legend.position = "none"
    ) +
    ylim(0, max(plot_data$Mean_Correlation) * 1.2)
  
  print(bar_plot)
  dev.off()
}

cat("\n✓ Pre-collinearity screening supplementary analysis completed\n")
cat("✓ Results saved to:", collinearity_vis_dir, "\n")

# ============================================================================
# 9. Part 1: Network analysis before collinearity screening (complete centrality version)
# ============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Part 1: Network Analysis Before Collinearity Screening (Complete Centrality Version)\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

# Before collinearity screening, use smaller gamma value (0.1) for denser network
network_results_before <- run_categorical_network_analysis_complete(
  analysis_data = continuous_vars_imputed,
  output_dir = "results/network_analysis_before_selection",
  analysis_name = "Network Analysis Before Collinearity Screening",
  gamma_value = 0.1  # Use smaller gamma for denser network
)

# ============================================================================
# 10. Part 2: Network analysis after collinearity screening (complete centrality version)
# ============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Part 2: Network Analysis After Collinearity Screening (Complete Centrality Version)\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

network_results_after <- run_categorical_network_analysis_complete(
  analysis_data = continuous_vars_final,
  output_dir = "results/network_analysis_after_selection",
  analysis_name = "Network Analysis After Collinearity Screening",
  gamma_value = 0.2  # Keep original gamma value
)

# ============================================================================
# 11. Centrality analysis and subsequent processing
# ============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("STEP 11: Centrality Analysis and Subsequent Processing\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# 11.1 Extract post-screening network results for subsequent analysis
cat("\nExtracting post-screening network results for subsequent analysis...\n")
network_matrix_after <- network_results_after$network_matrix
centrality_df_after <- network_results_after$centrality_df
network_p_values_after <- network_results_after$network_p_values

# 11.2 Extract pre-screening network results for comparison
cat("\nExtracting pre-screening network results for comparison...\n")
network_matrix_before <- network_results_before$network_matrix
centrality_df_before <- network_results_before$centrality_df
network_p_values_before <- network_results_before$network_p_values

# 11.3 Calculate comprehensive centrality metrics for post-screening data
cat("\nCalculating comprehensive centrality metrics for post-screening data...\n")

# Add ranking for each centrality metric
for (col in c("Strength", "Closeness", "Betweenness", "Eigenvector", "Degree", "ExpectedInfluence")) {
  centrality_df_after[[paste0(col, "_Rank")]] <- rank(-centrality_df_after[[col]], ties.method = "min")
}

# Calculate comprehensive ranking (mean)
rank_cols <- c("Strength_Rank", "Closeness_Rank", "Betweenness_Rank", 
               "Eigenvector_Rank", "Degree_Rank", "ExpectedInfluence_Rank")
centrality_df_after$Overall_Rank <- rowMeans(centrality_df_after[, rank_cols])

# Sort by comprehensive ranking
centrality_df_after <- centrality_df_after[order(centrality_df_after$Overall_Rank), ]

# Save centrality results
write.csv(centrality_df_after, 
          "results/centrality/centrality_results_after_selection.csv", 
          row.names = FALSE)

# 11.4 Calculate comprehensive centrality metrics for pre-screening data (for comparison)
cat("\nCalculating comprehensive centrality metrics for pre-screening data (for comparison)...\n")

for (col in c("Strength", "Closeness", "Betweenness", "Eigenvector", "Degree", "ExpectedInfluence")) {
  centrality_df_before[[paste0(col, "_Rank")]] <- rank(-centrality_df_before[[col]], ties.method = "min")
}

centrality_df_before$Overall_Rank <- rowMeans(centrality_df_before[, rank_cols])
centrality_df_before <- centrality_df_before[order(centrality_df_before$Overall_Rank), ]

write.csv(centrality_df_before, 
          "results/centrality/centrality_results_before_selection.csv", 
          row.names = FALSE)

# Identify top 6 central nodes (post-screening)
top6_center_vars_after <- centrality_df_after$Variable[1:6]
top6_center_names_after <- top6_center_vars_after

cat("\nTop 6 Central Nodes After Screening:\n")
for (i in 1:6) {
  cat(i, ". ", top6_center_vars_after[i], 
      " (Overall Rank:", round(centrality_df_after$Overall_Rank[i], 2), 
      ", Strength Centrality:", round(centrality_df_after$Strength[i], 3), ")\n", sep = "")
}

# Identify top 6 central nodes (pre-screening)
top6_center_vars_before <- centrality_df_before$Variable[1:6]
cat("\nTop 6 Central Nodes Before Screening:\n")
for (i in 1:6) {
  cat(i, ". ", top6_center_vars_before[i], 
      " (Overall Rank:", round(centrality_df_before$Overall_Rank[i], 2), 
      ", Strength Centrality:", round(centrality_df_before$Strength[i], 3), ")\n", sep = "")
}

# 11.5 Centrality analysis visualization
cat("\n11.5 Creating centrality visualizations...\n")

# Create centrality directory
if (!dir.exists("results/centrality/plots")) {
  dir.create("results/centrality/plots", recursive = TRUE)
}

# 11.5.1 Radar chart for top 6 central nodes after screening
top6_nodes <- centrality_df_after$Variable[1:6]
top6_data <- centrality_df_after[centrality_df_after$Variable %in% top6_nodes, ]

# Ensure sorted by overall rank (from 1 to 6)
top6_data <- top6_data[order(top6_data$Overall_Rank), ]

# Reorder factor levels to ensure facets display in rank order
top6_data$Variable <- factor(top6_data$Variable, levels = top6_data$Variable)

# Standardization
normalize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

top6_norm <- top6_data
cent_measures <- c("Strength", "Closeness", "Betweenness", "Eigenvector", "Degree", "ExpectedInfluence")

for (measure in cent_measures) {
  top6_norm[[paste0(measure, "_norm")]] <- normalize(top6_norm[[measure]])
}

# Convert to long format
top6_long <- top6_norm %>%
  dplyr::select(Variable, ends_with("_norm")) %>%
  pivot_longer(cols = -Variable, names_to = "Centrality", values_to = "Value") %>%
  mutate(Centrality = gsub("_norm", "", Centrality))

# Set display order
top6_long$Centrality <- factor(top6_long$Centrality,
                               levels = cent_measures)

# Ensure Variable factor in rank order
top6_long$Variable <- factor(top6_long$Variable, levels = top6_data$Variable)

# Add rank labels
rank_labels <- paste0("Rank ", 1:6, ": ", top6_data$Variable)
names(rank_labels) <- top6_data$Variable

# Create radar chart
node_colors <- c("#FF6B6B", "#4ECDC4", "#FFD166", "#95E1D3", "#A8DADC", "#DDA0DD")

p_radar <- ggplot(top6_long, aes(x = Centrality, y = Value, group = Variable, color = Variable)) +
  geom_polygon(fill = NA, size = 1.2, alpha = 0.8) +
  geom_point(size = 4) +
  geom_line(size = 1.2, alpha = 0.7) +
  facet_wrap(~ Variable, ncol = 3, labeller = labeller(Variable = rank_labels)) +
  coord_polar() +
  ylim(0, 1.1) +
  scale_color_manual(values = node_colors) +
  labs(
    title = "Top 6 Central Nodes After Screening - Centrality Distribution",
    subtitle = "Standardized Centrality Metrics (0-1 Range) - Sorted by Overall Rank",
    x = "",
    y = "Standardized Value"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 14),
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "none",
    panel.spacing = unit(25, "pt")
  )

ggsave("results/centrality/plots/centrality_radar_top6_after_selection.png", 
       p_radar, width = 22, height = 14, dpi = 300)

cat("✓ Post-screening radar chart created\n")

# 11.5.2 Pre- and post-screening centrality comparison
cat("\nCreating pre- and post-screening centrality comparison graph...\n")

# Mark pre- and post-screening data
centrality_df_before$Selection <- "Before Screening"
centrality_df_after$Selection <- "After Screening"

# Merge data
common_vars <- intersect(centrality_df_before$Variable, centrality_df_after$Variable)
centrality_combined <- rbind(
  centrality_df_before[centrality_df_before$Variable %in% common_vars, ],
  centrality_df_after[centrality_df_after$Variable %in% common_vars, ]
)

# Draw strength centrality comparison
p_strength_compare <- ggplot(centrality_combined, aes(x = Variable, y = Strength, fill = Selection)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(
    title = "Strength Centrality Comparison: Before vs After Screening",
    x = "Variable",
    y = "Strength Centrality"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  scale_fill_brewer(palette = "Set2")

ggsave("results/centrality/plots/strength_comparison_before_after.png", 
       p_strength_compare, width = 16, height = 10, dpi = 300)

# 11.6 Save key results
cat("\n11.6 Saving key results...\n")

# Save top 6 node data after screening
top6_summary <- centrality_df_after[1:6, c("Variable", "Strength", "Closeness", "Betweenness", 
                                           "Eigenvector", "Degree", "ExpectedInfluence", "Overall_Rank")]
write.csv(top6_summary, 
          "results/centrality/top6_central_nodes_after_selection.csv", 
          row.names = FALSE)

# Save pre- and post-screening comparison summary
comparison_summary <- data.frame(
  Variable = common_vars,
  Strength_Before = centrality_df_before$Strength[match(common_vars, centrality_df_before$Variable)],
  Strength_After = centrality_df_after$Strength[match(common_vars, centrality_df_after$Variable)],
  Rank_Before = centrality_df_before$Overall_Rank[match(common_vars, centrality_df_before$Variable)],
  Rank_After = centrality_df_after$Overall_Rank[match(common_vars, centrality_df_after$Variable)]
)

# Calculate rank changes
comparison_summary$Rank_Change <- comparison_summary$Rank_Before - comparison_summary$Rank_After
comparison_summary$Strength_Change <- comparison_summary$Strength_After - comparison_summary$Strength_Before

write.csv(comparison_summary, 
          "results/centrality/centrality_comparison_before_after.csv", 
          row.names = FALSE)

cat("✓ Centrality analysis completed\n")

# ============================================================================
# 12. Multi-Center In-Depth Analysis (Based on Post-Screening Results)
# ============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PART 3: Multi-Center In-Depth Analysis (Based on Post-Screening Results)\n")
cat("Analyzing Top 6 Central Nodes\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Create multi-center analysis directory
multi_center_dir <- "results/Multi_Center_Analysis"
if (!dir.exists(multi_center_dir)) {
  dir.create(multi_center_dir, recursive = TRUE)
}

# 12.1 Center variable analysis function
analyze_center_variable_advanced <- function(center_var, network_matrix, center_name, centrality_df, network_p_values = NULL) {
  
  cat("In-depth analysis of central variable:", center_name, "(", center_var, ")\n")
  
  if (!(center_var %in% colnames(network_matrix))) {
    cat("Warning: Variable", center_var, "not in network matrix\n")
    return(NULL)
  }
  
  # Extract center variable relationships
  center_index <- which(colnames(network_matrix) == center_var)
  center_relationships <- network_matrix[, center_index]
  
  # Create analysis data frame
  center_df <- data.frame(
    Variable = names(center_relationships),
    Partial_Correlation = as.numeric(center_relationships),
    Abs_Correlation = abs(as.numeric(center_relationships)),
    stringsAsFactors = FALSE
  )
  
  # Remove self-correlation
  center_df <- center_df[center_df$Variable != center_var, ]
  
  # If p-value information exists, add p-values
  if (!is.null(network_p_values)) {
    p_values <- network_p_values[center_index, ]
    center_df$P_Value <- p_values[match(center_df$Variable, names(p_values))]
    center_df$Significant <- ifelse(center_df$P_Value < 0.05, "Yes", "No")
    center_df$Line_Type <- ifelse(center_df$P_Value < 0.05, "Solid", "Dashed")
  }
  
  # Sort by absolute correlation
  center_df <- center_df[order(abs(center_df$Partial_Correlation), decreasing = TRUE), ]
  
  # Get top and bottom 5 associations
  top5 <- head(center_df, 5)
  bottom5 <- tail(center_df[order(abs(center_df$Partial_Correlation), decreasing = FALSE), ], 5)
  
  # Category analysis
  variable_categories <- list(
    "Adiposity_Indicators" = c("VFA", "VFM", "VFV"),
    "Vitamin_D" = c("Vit_D", "Vit_D3"),
    "Sex_Hormones" = c("Testo", "E2", "SHBG"),
    "Bone_Metabolism" = c("Total_BMC", "Total_BMD"),
    "Metabolic_Indicators" = c("TSH", "HDL", "TC", "GLU", "hs_CRP")
  )
  
  # Adjust categories
  if (center_var %in% c("E2", "SHBG", "Testo")) {
    variable_categories[["Sex_Hormones"]] <- setdiff(variable_categories[["Sex_Hormones"]], center_var)
  }
  
  # Calculate average association for each category
  category_analysis <- list()
  for (category_name in names(variable_categories)) {
    vars <- variable_categories[[category_name]]
    vars <- vars[vars %in% center_df$Variable]
    
    if (length(vars) > 0) {
      cat_data <- center_df[center_df$Variable %in% vars, ]
      avg_corr <- mean(abs(cat_data$Partial_Correlation))
      direction_ratio <- sum(cat_data$Partial_Correlation > 0) / nrow(cat_data)
      
      category_analysis[[category_name]] <- list(
        data = cat_data,
        avg_corr = avg_corr,
        direction_ratio = direction_ratio,
        n_vars = nrow(cat_data)
      )
    }
  }
  
  # Get centrality information
  center_info <- centrality_df[centrality_df$Variable == center_var, ]
  
  return(list(
    center_var = center_var,
    center_name = center_name,
    relationships = center_df,
    top5 = top5,
    bottom5 = bottom5,
    category_analysis = category_analysis,
    centrality_info = center_info
  ))
}

# 12.2 Execute multi-center analysis
cat("\nStarting multi-center in-depth analysis (top 6 central nodes)...\n")

# Analyze top 6 central nodes (after screening)
top6_center_vars <- top6_center_vars_after

for (i in 1:length(top6_center_vars)) {
  center_var <- top6_center_vars[i]
  center_name <- center_var
  
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("Analyzing center node", i, ":", center_name, "(", center_var, ")\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")
  
  # Create analysis directory for this node
  center_dir <- file.path(multi_center_dir, center_var)
  if (!dir.exists(center_dir)) {
    dir.create(center_dir, recursive = TRUE)
    dir.create(file.path(center_dir, "plots"), recursive = TRUE)
  }
  
  # Analyze center variable
  results <- analyze_center_variable_advanced(center_var, network_matrix_after, center_name, 
                                              centrality_df_after, network_p_values_after)
  
  if (!is.null(results)) {
    # Save analysis results
    write.csv(results$relationships, 
              file.path(center_dir, paste0(center_var, "_partial_correlations.csv")), 
              row.names = FALSE)
    write.csv(results$top5, 
              file.path(center_dir, paste0(center_var, "_top5_associations.csv")), 
              row.names = FALSE)
    
    # Generate brief report
    report_file <- file.path(center_dir, paste0(center_var, "_analysis_summary.txt"))
    sink(report_file)
    cat("=== ", center_name, " Analysis Summary ===\n\n")
    cat("Centrality Information:\n")
    cat("Strength Centrality:", results$centrality_info$Strength, "\n")
    cat("Degree Centrality:", results$centrality_info$Degree, "\n")
    cat("Expected Influence:", results$centrality_info$ExpectedInfluence, "\n")
    cat("Overall Rank:", round(results$centrality_info$Overall_Rank, 2), "\n\n")
    
    cat("Strongest Associations (Top 5):\n")
    for (j in 1:nrow(results$top5)) {
      cat(j, ". ", results$top5$Variable[j], ": ", 
          round(results$top5$Partial_Correlation[j], 4), 
          ifelse("P_Value" %in% colnames(results$top5), 
                 paste(" (p=", round(results$top5$P_Value[j], 4), ")", sep = ""), ""),
          "\n", sep = "")
    }
    sink()
    
    cat("✓ Single-center analysis completed\n")
    cat("✓ Results saved to:", center_dir, "\n")
    cat("✓ Report file:", report_file, "\n\n")
  }
}

cat("✓ Multi-center analysis completed\n")

# ============================================================================
# 13. Single-Center In-Depth Analysis Function
# ============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PART 4: Single-Center In-Depth Analysis\n")
cat("Advanced visualization for center nodes\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Create single-center analysis directory
single_center_dir <- "results/Single_Center_Analysis"
if (!dir.exists(single_center_dir)) {
  dir.create(single_center_dir, recursive = TRUE)
}

# 13.1 Enhanced visualization function for center variables
visualize_center_variable_enhanced <- function(center_var, network_matrix, centrality_df, 
                                               variable_categories_vec, output_dir, 
                                               n_samples, network_p_values = NULL) {
  
  cat("\nCreating enhanced visualizations for center variable:", center_var, "\n")
  
  # Create directory for this center variable
  center_vis_dir <- file.path(output_dir, center_var)
  if (!dir.exists(center_vis_dir)) {
    dir.create(center_vis_dir, recursive = TRUE)
    dir.create(file.path(center_vis_dir, "plots"), recursive = TRUE)
  }
  
  if (!(center_var %in% colnames(network_matrix))) {
    cat("Warning: Variable", center_var, "not in network matrix\n")
    return(NULL)
  }
  
  # Extract center variable relationships
  center_index <- which(colnames(network_matrix) == center_var)
  center_relationships <- network_matrix[, center_index]
  
  # Create analysis data frame
  center_df <- data.frame(
    Variable = names(center_relationships),
    Partial_Correlation = as.numeric(center_relationships),
    Abs_Correlation = abs(as.numeric(center_relationships)),
    stringsAsFactors = FALSE
  )
  
  # Remove self-correlation
  center_df <- center_df[center_df$Variable != center_var, ]
  
  # If p-value information exists, add p-values
  if (!is.null(network_p_values)) {
    p_values <- network_p_values[center_index, ]
    center_df$P_Value <- p_values[match(center_df$Variable, names(p_values))]
    center_df$Significant <- ifelse(center_df$P_Value < 0.05, "Yes", "No")
  }
  
  # Sort by absolute correlation
  center_df <- center_df[order(abs(center_df$Partial_Correlation), decreasing = TRUE), ]
  
  # Define color scheme for categories
  category_colors_enhanced <- c(
    "Metabolic Indicators" = "#66C2A5",    # Green
    "Vitamin D" = "#FC8D62",               # Orange
    "Adiposity Indicators" = "#8DA0CB",    # Blue
    "Bone Metabolism" = "#E78AC3",         # Pink
    "Sex Hormones" = "#A6D854",            # Light green
    "Other" = "#B3B3B3"                    # Gray
  )
  
  # Assign categories
  center_df$Category <- variable_categories_vec[match(center_df$Variable, names(variable_categories_vec))]
  center_df$Category[is.na(center_df$Category)] <- "Other"
  
  # Ensure category factor is properly ordered
  center_df$Category <- factor(center_df$Category, 
                               levels = names(category_colors_enhanced))
  
  # Add significance markers if available
  if ("P_Value" %in% colnames(center_df)) {
    center_df$Significance <- ifelse(center_df$P_Value < 0.001, "***",
                                     ifelse(center_df$P_Value < 0.01, "**",
                                            ifelse(center_df$P_Value < 0.05, "*", "")))
  }
  
  # 13.1.1 Enhanced relationship bar chart
  cat("Creating enhanced relationship bar chart...\n")
  
  p_bar_enhanced <- ggplot(center_df, aes(x = reorder(Variable, Partial_Correlation), 
                                          y = Partial_Correlation, 
                                          fill = Category)) +
    geom_bar(stat = "identity", width = 0.7, alpha = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
    geom_hline(yintercept = c(-0.1, 0.1), linetype = "dotted", color = "gray70", alpha = 0.7, size = 0.5) +
    
    # Add significance markers if available
    {if ("Significance" %in% colnames(center_df) && any(center_df$Significance != "")) 
      geom_text(aes(label = Significance), vjust = -0.8, size = 6, fontface = "bold")} +
    
    scale_fill_manual(values = category_colors_enhanced) +
    coord_flip() +
    labs(
      title = paste(center_var, "Partial Correlation Relationships"),
      subtitle = paste("Unique associations after controlling for other variables | Sample size:", n_samples),
      x = "Variables",
      y = "Partial Correlation Coefficient",
      caption = "Positive = positive correlation, Negative = negative correlation\n* p < 0.05, ** p < 0.01, *** p < 0.001"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20, margin = margin(b = 10)),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 14, margin = margin(b = 15)),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      plot.caption = element_text(size = 10, color = "gray50", hjust = 0),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    ylim(min(center_df$Partial_Correlation) * 1.1, max(center_df$Partial_Correlation) * 1.1)
  
  ggsave(file.path(center_vis_dir, "plots", paste0(center_var, "_relationships_bar_enhanced.png")), 
         p_bar_enhanced, width = 16, height = 14, dpi = 300)
  
  # 13.1.2 Subnetwork visualization
  cat("Creating subnetwork visualization...\n")
  
  # Select variables with strongest connections (top 8 + center variable)
  strong_connections <- center_df[abs(center_df$Partial_Correlation) > 0.05, "Variable"]
  if (length(strong_connections) > 8) {
    strong_connections <- head(strong_connections, 8)
  }
  
  if (length(strong_connections) > 0) {
    focus_vars <- c(center_var, strong_connections)
    focus_vars <- focus_vars[focus_vars %in% colnames(network_matrix)]
    
    if (length(focus_vars) > 1) {
      subnetwork_matrix <- network_matrix[focus_vars, focus_vars]
      
      # Get colors for subnetwork
      subnetwork_colors <- category_colors_enhanced[variable_categories_vec[focus_vars]]
      subnetwork_colors[is.na(subnetwork_colors)] <- category_colors_enhanced["Other"]
      
      # Highlight center node
      node_sizes <- ifelse(focus_vars == center_var, 15, 10)
      
      png(file.path(center_vis_dir, "plots", paste0(center_var, "_subnetwork_enhanced.png")), 
          width = 1600, height = 1400, res = 150)
      
      # Adjust margins for better display
      par(mar = c(1, 1, 3, 1))
      
      qgraph(subnetwork_matrix,
             layout = "spring",
             title = paste(center_var, "Association Subnetwork"),
             labels = focus_vars,
             color = subnetwork_colors,
             vsize = node_sizes,
             esize = 12,
             edge.width = 3,
             posCol = "#0D47A1",
             negCol = "#B71C1C",
             label.cex = 1.4,
             title.cex = 1.8,
             edge.labels = TRUE,
             edge.label.cex = 1.0,
             mar = c(5, 5, 5, 5))
      dev.off()
    }
  }
  
  # 13.1.3 Correlation direction pie chart
  cat("Creating correlation direction pie chart...\n")
  
  direction_summary <- data.frame(
    Direction = c("Strong Positive\n(r > 0.1)", "Weak Positive\n(0 < r ≤ 0.1)", 
                  "Weak Negative\n(-0.1 ≤ r < 0)", "Strong Negative\n(r < -0.1)", 
                  "No Correlation\n(|r| ≤ 0.01)"),
    Count = c(
      sum(center_df$Partial_Correlation > 0.1),
      sum(center_df$Partial_Correlation > 0.01 & center_df$Partial_Correlation <= 0.1),
      sum(center_df$Partial_Correlation >= -0.1 & center_df$Partial_Correlation < -0.01),
      sum(center_df$Partial_Correlation < -0.1),
      sum(abs(center_df$Partial_Correlation) <= 0.01)
    )
  )
  
  # Filter out zero counts
  direction_summary <- direction_summary[direction_summary$Count > 0, ]
  
  p_pie_enhanced <- ggplot(direction_summary, aes(x = "", y = Count, fill = Direction)) +
    geom_bar(stat = "identity", width = 1, color = "white", size = 0.5) +
    coord_polar("y", start = 0) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = paste(center_var, "Association Direction Distribution"),
      fill = "Association Direction",
      caption = paste("Total variables:", nrow(center_df))
    ) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 10)),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      plot.caption = element_text(size = 11, color = "gray40", hjust = 0.5, margin = margin(t = 10))
    )
  
  ggsave(file.path(center_vis_dir, "plots", paste0(center_var, "_direction_pie_enhanced.png")), 
         p_pie_enhanced, width = 12, height = 10, dpi = 300)
  
  # 13.1.4 Strength comparison with other variables
  cat("Creating strength comparison visualization...\n")
  
  # Get centrality information for all variables
  all_strength <- centrality_df$Strength
  names(all_strength) <- centrality_df$Variable
  
  # Create comparison data
  strength_comparison <- data.frame(
    Variable = names(all_strength),
    Strength = all_strength,
    Is_Center = ifelse(names(all_strength) == center_var, "Center", "Other"),
    stringsAsFactors = FALSE
  )
  
  # Sort by strength
  strength_comparison <- strength_comparison[order(strength_comparison$Strength, decreasing = TRUE), ]
  
  # Highlight top 10
  strength_comparison$Rank <- 1:nrow(strength_comparison)
  strength_comparison$Highlight <- ifelse(strength_comparison$Rank <= 10, "Top 10", "Other")
  strength_comparison$Highlight[strength_comparison$Variable == center_var] <- "Center"
  
  p_strength_compare <- ggplot(strength_comparison, aes(x = reorder(Variable, Strength), y = Strength, 
                                                        fill = Highlight)) +
    geom_bar(stat = "identity", width = 0.7, alpha = 0.9) +
    scale_fill_manual(values = c("Center" = "#FF6B6B", "Top 10" = "#4ECDC4", "Other" = "#B3B3B3")) +
    coord_flip() +
    labs(
      title = paste("Strength Centrality Comparison:", center_var, "vs Other Variables"),
      subtitle = "Centrality ranking based on network analysis",
      x = "Variables",
      y = "Strength Centrality",
      caption = paste("Center variable rank:", which(strength_comparison$Variable == center_var), 
                      "of", nrow(strength_comparison))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 10)),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 14, margin = margin(b = 15)),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 11),
      axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 10)),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      plot.caption = element_text(size = 10, color = "gray50", hjust = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(center_vis_dir, "plots", paste0(center_var, "_strength_comparison.png")), 
         p_strength_compare, width = 16, height = 14, dpi = 300)
  
  # Save data
  write.csv(center_df, file.path(center_vis_dir, paste0(center_var, "_detailed_relationships.csv")), 
            row.names = FALSE)
  
  cat("✓ Enhanced visualizations created for", center_var, "\n")
  cat("✓ Results saved to:", center_vis_dir, "\n")
  
  return(center_vis_dir)
}

# 13.2 Execute enhanced single-center visualization for top 6 central nodes
cat("\nStarting enhanced single-center visualization for top 6 central nodes...\n")

for (i in 1:length(top6_center_vars)) {
  center_var <- top6_center_vars[i]
  
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("Creating enhanced visualizations for center node", i, ":", center_var, "\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")
  
  # Get variable categories vector
  variable_categories_vec_after <- network_results_after$variable_categories
  
  # Create enhanced visualizations
  center_vis_dir <- visualize_center_variable_enhanced(
    center_var = center_var,
    network_matrix = network_matrix_after,
    centrality_df = centrality_df_after,
    variable_categories_vec = variable_categories_vec_after,
    output_dir = single_center_dir,
    n_samples = network_results_after$n_samples,
    network_p_values = network_p_values_after
  )
  
  if (!is.null(center_vis_dir)) {
    cat("✓ Enhanced single-center visualization completed for", center_var, "\n")
    cat("✓ Results saved to:", center_vis_dir, "\n\n")
  }
}

cat("✓ Enhanced single-center visualization completed\n")

# 13.3 Generate comprehensive report for each center node
generate_comprehensive_center_report <- function(center_var, results, output_dir, selection_results, 
                                                 n_samples, pcor_result = NULL) {
  
  if (is.null(results)) return(NULL)
  
  center_name <- results$center_name
  
  # Generate detailed report
  report_file <- file.path(output_dir, paste0(center_var, "_comprehensive_analysis_report.txt"))
  
  sink(report_file)
  
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat(center_name, "Comprehensive Network Analysis Report\n")
  cat("Generated:", Sys.time(), "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  cat("1. Basic Information\n")
  cat("Central Variable:", center_name, "(", center_var, ")\n")
  cat("Total Variables Analyzed:", nrow(results$relationships) + 1, "\n")
  cat("Sample Size:", n_samples, "\n")
  cat("Age Range: 7-18 years\n\n")
  
  cat("2. Variable Selection Summary\n")
  cat("Original Variables:", ncol(continuous_vars_imputed), "\n")
  cat("Selected Variables:", length(selection_results$selected_variables), "\n")
  cat("Removed Variables:", length(selection_results$removed_variables), "\n\n")
  
  cat("3. Centrality Information\n")
  cat("Strength Centrality:", round(results$centrality_info$Strength, 3), "\n")
  cat("Closeness Centrality:", round(results$centrality_info$Closeness, 3), "\n")
  cat("Betweenness Centrality:", round(results$centrality_info$Betweenness, 3), "\n")
  cat("Eigenvector Centrality:", round(results$centrality_info$Eigenvector, 3), "\n")
  cat("Degree:", results$centrality_info$Degree, "\n")
  cat("Expected Influence:", round(results$centrality_info$ExpectedInfluence, 3), "\n")
  cat("Overall Rank:", round(results$centrality_info$Overall_Rank, 2), "\n\n")
  
  cat("4. Key Association Analysis\n")
  cat("4.1 Top 5 Strongest Associations:\n")
  for (i in 1:nrow(results$top5)) {
    var <- results$top5$Variable[i]
    corr <- results$top5$Partial_Correlation[i]
    p_value <- ifelse("P_Value" %in% colnames(results$top5), 
                      results$top5$P_Value[i], NA)
    
    cat("   ", i, ". ", var, ": ", round(corr, 4), sep = "")
    
    if (!is.na(p_value)) {
      cat(" (p = ", round(p_value, 4), ")", sep = "")
    }
    
    cat(" - ", ifelse(corr > 0, "Positive", "Negative"), " correlation\n", sep = "")
  }
  
  cat("\n4.2 Top 5 Weakest Associations:\n")
  for (i in 1:nrow(results$bottom5)) {
    var <- results$bottom5$Variable[i]
    corr <- results$bottom5$Partial_Correlation[i]
    p_value <- ifelse("P_Value" %in% colnames(results$bottom5), 
                      results$bottom5$P_Value[i], NA)
    
    cat("   ", i, ". ", var, ": ", round(corr, 4), sep = "")
    
    if (!is.na(p_value)) {
      cat(" (p = ", round(p_value, 4), ")", sep = "")
    }
    
    cat("\n")
  }
  
  cat("\n5. Association Statistics Summary\n")
  total_vars <- nrow(results$relationships)
  strong_pos <- sum(results$relationships$Partial_Correlation > 0.1)
  weak_pos <- sum(results$relationships$Partial_Correlation > 0.01 & 
                    results$relationships$Partial_Correlation <= 0.1)
  weak_neg <- sum(results$relationships$Partial_Correlation >= -0.1 & 
                    results$relationships$Partial_Correlation < -0.01)
  strong_neg <- sum(results$relationships$Partial_Correlation < -0.1)
  no_corr <- sum(abs(results$relationships$Partial_Correlation) <= 0.01)
  avg_abs_corr <- mean(abs(results$relationships$Partial_Correlation))
  
  cat("Total Associations:", total_vars, "\n")
  cat("Strong Positive (r > 0.1):", strong_pos, " (", round(strong_pos/total_vars*100, 1), "%)\n", sep = "")
  cat("Weak Positive (0.01 < r ≤ 0.1):", weak_pos, " (", round(weak_pos/total_vars*100, 1), "%)\n", sep = "")
  cat("Weak Negative (-0.1 ≤ r < -0.01):", weak_neg, " (", round(weak_neg/total_vars*100, 1), "%)\n", sep = "")
  cat("Strong Negative (r < -0.1):", strong_neg, " (", round(strong_neg/total_vars*100, 1), "%)\n", sep = "")
  cat("No Correlation (|r| ≤ 0.01):", no_corr, " (", round(no_corr/total_vars*100, 1), "%)\n", sep = "")
  cat("Mean Absolute Association Strength:", round(avg_abs_corr, 4), "\n")
  
  # If p-value information exists
  if ("P_Value" %in% colnames(results$relationships)) {
    sig_count <- sum(results$relationships$P_Value < 0.05, na.rm = TRUE)
    cat("Significant Associations (p < 0.05):", sig_count, " (", round(sig_count/total_vars*100, 1), "%)\n", sep = "")
  }
  
  cat("\n6. Biological and Clinical Significance\n")
  
  # Provide specific biological interpretation based on center variable
  if (center_var == "SHBG") {
    cat("- Sex Hormone Binding Globulin (SHBG) is a key protein regulating sex hormone bioavailability.\n")
    cat("- Associations with sex hormones (Testo, E2) reflect SHBG's role in modulating free hormone levels.\n")
    cat("- Relationships with metabolic indicators (TC, HDL, GLU) suggest SHBG's involvement in glucose and lipid metabolism.\n")
    cat("- Association with inflammatory marker (hs_CRP) may indicate chronic inflammation's impact on SHBG levels.\n")
    cat("- Relationships with adiposity indicators (VFA, VFM) support the concept of adipose tissue as an endocrine organ.\n")
    cat("- Clinical significance: SHBG levels can serve as a biomarker for metabolic health; low SHBG is often associated with insulin resistance and metabolic syndrome.\n")
    cat("- In adolescents, changes in SHBG levels may reflect dynamic changes in pubertal development and metabolic status.\n")
  } else if (center_var %in% c("VFA", "VFM", "VFV")) {
    cat("- Visceral fat indicators reflect central adiposity.\n")
    cat("- Closely related to metabolic syndrome and insulin resistance.\n")
    cat("- May serve as a hub connecting adipose tissue endocrine function and systemic metabolism.\n")
    cat("- In adolescents, visceral fat accumulation is a key predictor of future metabolic disease risk.\n")
  } else if (center_var %in% c("Testo", "E2")) {
    cat("- Sex hormones play crucial roles in pubertal development and metabolic regulation.\n")
    cat("- Associated with body composition, bone metabolism, lipid metabolism, and multiple systems.\n")
    cat("- May reflect central position in hormone-metabolic axis.\n")
    cat("- In adolescents, sex hormone levels are crucial for normal growth and development.\n")
  } else if (center_var %in% c("Total_BMD", "Total_BMC")) {
    cat("- Bone mineral density/content reflects skeletal health and development.\n")
    cat("- Associated with growth, nutrition, hormonal status, and physical activity.\n")
    cat("- In adolescents, bone accumulation is crucial for achieving peak bone mass.\n")
  } else {
    cat("- This variable acts as a central hub in the biomarker network.\n")
    cat("- Its central position suggests it may regulate or be regulated by multiple physiological processes.\n")
    cat("- Further research is needed to understand its specific biological mechanisms in adolescents.\n")
  }
  
  cat("\n7. Methodology\n")
  cat("7.1 Data Processing\n")
  cat("- Age Range: 7-18 years\n")
  cat("- Missing Data: Median imputation\n")
  cat("- Variable Selection: Based on collinearity analysis\n")
  cat("- Selected Variables:", paste(colnames(continuous_vars_final), collapse = ", "), "\n\n")
  
  cat("7.2 Network Analysis\n")
  cat("- Estimation Method: EBICglasso\n")
  cat("- Regularization Parameter: 0.2\n")
  cat("- Centrality Metrics: Strength, Closeness, Betweenness, Eigenvector, Degree, Expected Influence\n")
  cat("- Single-Center Analysis: Detailed examination of central node's network relationships\n\n")
  
  cat("8. Limitations and Future Directions\n")
  cat("8.1 Limitations\n")
  cat("- Cross-sectional design limits causal inference\n")
  cat("- Sample limited to adolescents (7-18 years)\n")
  cat("- Network analysis is exploratory; mechanisms require experimental validation\n")
  cat("- Variable selection may exclude some biologically relevant variables\n\n")
  
  cat("8.2 Future Directions\n")
  cat("- Conduct longitudinal studies to examine temporal relationships\n")
  cat("- Stratified analysis by sex to explore gender differences\n")
  cat("- Integrate genomic and metabolomic data for multi-omics analysis\n")
  cat("- Validate findings in independent cohorts\n")
  cat("- Examine interactions with lifestyle factors (diet, physical activity)\n")
  
  cat("\n9. Visualization Summary\n")
  cat("The following visualizations have been generated:\n")
  cat("1. Enhanced relationship bar chart showing all associations\n")
  cat("2. Subnetwork visualization focusing on strongest connections\n")
  cat("3. Association direction distribution pie chart\n")
  cat("4. Strength centrality comparison with other variables\n")
  cat("5. Comprehensive statistical summary tables\n")
  
  sink()
  
  cat("✓ Comprehensive report generated for", center_var, "\n")
  
  return(report_file)
}

# 13.4 Generate comprehensive reports for top 6 central nodes
cat("\nGenerating comprehensive reports for top 6 central nodes...\n")

for (i in 1:length(top6_center_vars)) {
  center_var <- top6_center_vars[i]
  
  # Get analysis results
  results <- analyze_center_variable_advanced(center_var, network_matrix_after, center_var, 
                                              centrality_df_after, network_p_values_after)
  
  if (!is.null(results)) {
    # Generate comprehensive report
    report_file <- generate_comprehensive_center_report(
      center_var = center_var,
      results = results,
      output_dir = file.path(single_center_dir, center_var),
      selection_results = selection_results,
      n_samples = network_results_after$n_samples,
      pcor_result = if(exists("pcor_result")) pcor_result else NULL
    )
    
    cat("✓ Comprehensive report generated for", center_var, "\n")
    cat("  Report file:", report_file, "\n\n")
  }
}

cat("✓ Comprehensive reporting completed\n")

# 14. Create final summary report -------------------------------------------------------
cat("\nCreating final summary report...\n")

final_report_file <- "results/final_summary_report.txt"
sink(final_report_file)

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("NHANES Network Analysis - Final Summary Report\n")
cat("Generated:", Sys.time(), "\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("1. Executive Summary\n")
cat("This analysis applied network methods to study relationships between biomarkers in NHANES dataset for adolescents aged 7-18.\n")
cat("Analysis included:\n")
cat("1. Collinearity-based variable selection\n")
cat("2. Network estimation using EBICglasso (line types set based on significance)\n")
cat("3. Partial correlation analysis and statistical inference\n")
cat("4. Network analysis before and after collinearity screening\n")
cat("5. Identification of top 6 central nodes and in-depth analysis\n\n")

cat("2. Data Characteristics\n")
cat("Sample Size:", nrow(continuous_vars_final), "\n")
cat("Variables Before Screening:", ncol(continuous_vars_imputed), "\n")
cat("Variables After Screening:", ncol(continuous_vars_final), "\n")
cat("Selected Variables:", paste(colnames(continuous_vars_final), collapse = ", "), "\n\n")

cat("3. Variable Selection Summary\n")
cat("Original Variables:", ncol(continuous_vars_imputed), "\n")
cat("Selected Variables:", length(selection_results$selected_variables), "\n")
cat("Removed Variables:", length(selection_results$removed_variables), "\n")
cat("Removal Reason: High intra-group collinearity\n\n")

cat("4. Network Analysis Settings\n")
cat("Network Estimation Method: EBICglasso\n")
cat("Regularization Parameter (gamma): 0.2\n")
cat("Display Settings: Show all non-zero connections\n")
cat("Line Type Settings: Significant connections (p<0.05) use solid lines, non-significant connections (p≥0.05) use dashed lines\n\n")

cat("5. Centrality Analysis Results (After Screening)\n")
cat("Top 6 Central Nodes (by Overall Rank):\n")
for (i in 1:min(6, nrow(centrality_df_after))) {
  cat(i, ". ", centrality_df_after$Variable[i], 
      " (Overall Rank:", round(centrality_df_after$Overall_Rank[i], 2), 
      ", Strength:", round(centrality_df_after$Strength[i], 3), 
      ", Degree:", centrality_df_after$Degree[i], ")\n", sep = "")
}
cat("\n")

cat("6. Partial Correlation Analysis Results\n")
if (exists("pcor_result")) {
  cat("Partial correlation matrix calculated and saved\n")
  cat("Significant partial correlation edges:", ifelse(exists("sig_edges_df"), nrow(sig_edges_df), "Not calculated"), "\n")
} else {
  cat("Partial correlation analysis not completed\n")
}
cat("\n")

cat("7. Multi-Center Analysis\n")
cat("In-depth analysis conducted for each of the top 6 central nodes:\n")
for (i in 1:6) {
  cat(i, ". ", top6_center_vars[i], " (Detailed analysis saved in: results/Multi_Center_Analysis/", 
      top6_center_vars[i], "/)\n", sep = "")
}
cat("\n")

cat("8. Enhanced Single-Center Visualization\n")
cat("Enhanced visualizations created for each central node with:\n")
cat("- Detailed relationship bar charts\n")
cat("- Subnetwork visualizations\n")
cat("- Association direction distributions\n")
cat("- Strength centrality comparisons\n")
cat("- Comprehensive statistical reports\n")
cat("Results saved in: results/Single_Center_Analysis/\n\n")

cat("9. Output File Structure\n")
cat("All analysis outputs saved in the following directories:\n")
cat("1. results/collinearity_analysis/ - Collinearity analysis results\n")
cat("2. results/partial_correlations/ - Partial correlation and statistical inference results\n")
cat("3. results/centrality/ - Centrality analysis results\n")
cat("4. results/Multi_Center_Analysis/ - Multi-center detailed analysis results\n")
cat("5. results/Single_Center_Analysis/ - Enhanced single-center visualization results\n")
cat("6. results/network_analysis_before_selection/ - Pre-collinearity screening network analysis\n")
cat("7. results/network_analysis_after_selection/ - Post-collinearity screening network analysis\n\n")

cat("10. Key Findings\n")
cat("Network analysis identified 6 central hub nodes, indicating these variables play key integrative roles\n")
cat("in adolescent biomarker networks. Discovery of these central nodes suggests:\n")
cat("- Complex network structure of multi-system physiological regulation\n")
cat("- Tight interactions between pubertal development, metabolic, and endocrine systems\n")
cat("- Potential of key biomarkers as targets for health assessment and intervention\n\n")

cat("11. Recommendations and Next Steps\n")
cat("Recommendations:\n")
cat("1. Review detailed analysis reports for each central node\n")
cat("2. Examine visualizations to understand network relationships\n")
cat("3. Consider subgroup analysis by sex and age groups\n")
cat("4. Validate findings in independent datasets\n")
cat("5. Explore biological mechanisms suggested by network structure\n\n")

sink()

cat("✓ Final summary report saved:", final_report_file, "\n")

# 15. Final Summary ---------------------------------------------------------
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("Analysis Completed!\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

cat("\nCompleted Analysis Steps:\n")
cat("1. ✓ Data preprocessing and cleaning\n")
cat("2. ✓ Collinearity analysis and variable selection\n")
cat("3. ✓ Partial correlation analysis and statistical inference\n")
cat("4. ✓ Network analysis before collinearity screening (complete centrality version)\n")
cat("5. ✓ Network analysis after collinearity screening (complete centrality version)\n")
cat("6. ✓ Centrality analysis (pre- and post-screening comparison)\n")
cat("7. ✓ Multi-center in-depth analysis (top 6 central nodes)\n")
cat("8. ✓ Enhanced single-center visualization (top 6 central nodes)\n")
cat("9. ✓ Comprehensive reporting\n")
cat("10. ✓ Generate comprehensive reports\n\n")

cat("Key Outputs:\n")
cat("1. Variable Selection: Reduced from", ncol(continuous_vars_imputed), "to", 
    ncol(continuous_vars_final), "variables\n")
cat("2. Central Node Identification (After Screening): Top 6 central nodes identified\n")
cat("3. Detailed Analysis Saved to: results/Multi_Center_Analysis/\n")
cat("4. Enhanced Visualizations Saved to: results/Single_Center_Analysis/\n")
cat("5. Partial Correlation and p-value matrices saved\n")
cat("6. Network analysis with line types set based on significance executed:\n")
cat("   - Significant connections (p<0.05) use solid lines\n")
cat("   - Non-significant connections (p≥0.05) use dashed lines\n")
cat("   - Before collinearity screening: results/network_analysis_before_selection/\n")
cat("   - After collinearity screening: results/network_analysis_after_selection/\n\n")

cat("Output Directory Structure:\n")
cat("results/\n")
cat("├── collinearity_analysis/          # Collinearity analysis results\n")
cat("├── partial_correlations/           # Partial correlation statistical inference\n")
cat("├── centrality/                     # Centrality analysis\n")
cat("│   ├── centrality_results_before_selection.csv\n")
cat("│   ├── centrality_results_after_selection.csv\n")
cat("│   ├── top6_central_nodes_after_selection.csv\n")
cat("│   └── plots/                      # Centrality visualizations\n")
cat("├── Multi_Center_Analysis/          # Multi-center detailed analysis\n")
cat("│   ├── [Center Node 1]/            # 1st central node\n")
cat("│   ├── [Center Node 2]/            # 2nd central node\n")
cat("│   ├── [Center Node 3]/            # 3rd central node\n")
cat("│   ├── [Center Node 4]/            # 4th central node\n")
cat("│   ├── [Center Node 5]/            # 5th central node\n")
cat("│   └── [Center Node 6]/            # 6th central node\n")
cat("├── Single_Center_Analysis/         # Enhanced single-center visualization\n")
cat("│   ├── [Center Node 1]/            # 1st central node visualizations\n")
cat("│   ├── [Center Node 2]/            # 2nd central node visualizations\n")
cat("│   ├── [Center Node 3]/            # 3rd central node visualizations\n")
cat("│   ├── [Center Node 4]/            # 4th central node visualizations\n")
cat("│   ├── [Center Node 5]/            # 5th central node visualizations\n")
cat("│   └── [Center Node 6]/            # 6th central node visualizations\n")
cat("├── network_analysis_before_selection/ # Pre-collinearity screening network analysis\n")
cat("│   ├── collinearity_visualization/    # Collinearity visualizations\n")
cat("│   └── plots/                         # Network graphs\n")
cat("└── network_analysis_after_selection/  # Post-collinearity screening network analysis\n")
cat("    └── plots/                         # Network graphs\n\n")

cat("=== All Analysis Completed ===\n")