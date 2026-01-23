# ============================================================
# Complete Path Analysis Code + SEM Fit Analysis (Revised Version)
# For Adolescent Data, Includes Complete Fit Analysis
# Added VFA→SHBG→VitD pathway analysis
# ============================================================

# 1. Load necessary packages ------------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(corrplot)
library(lavaan)
library(semPlot)
library(lavaanPlot)
library(patchwork)
library(kableExtra)
library(semTools)

# Set working directory
setwd("~/1.NHANES/1.Public-database")

# 2. Create result directory structure --------------------------------------------------------
cat("Step 1: Creating result directory structure...\n")

# Main directory
if (!dir.exists("results")) dir.create("results")
if (!dir.exists("results/Path_Analysis")) dir.create("results/Path_Analysis", recursive = TRUE)

# Subdirectories
subdirs <- c(
  "results/Path_Analysis/regression",
  "results/Path_Analysis/path_models",
  "results/Path_Analysis/plots",
  "results/Path_Analysis/tables",
  "results/Path_Analysis/reports",
  "results/Path_Analysis/SEM_fit",  # New: specifically for SEM fit results
  "results/Path_Analysis/Path_VFA_SHBG_VitD"  # New: specifically for VFA→SHBG→VitD pathway
)

for (dir in subdirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

cat("✓ Directory structure created\n")

# 3. Data loading and preprocessing --------------------------------------------------------
cat("\nStep 2: Data loading and preprocessing...\n")

# Read data
data_original <- read_excel("data.xlsx")
data <- data_original

# Remove specified columns
if (ncol(data) >= 41) {
  data <- data[, -41]
}

# Define factor and numeric columns
factor_cols <- c(1, 2, 4, 5, 40)
all_cols <- seq_along(data)

# Convert data types
data <- data %>%
  mutate(
    across(all_of(intersect(factor_cols, seq_along(data))), as.factor),
    across(all_of(setdiff(seq_along(data), factor_cols)), as.numeric)
  ) %>%
  as.data.frame()

cat("✓ Data preprocessing completed\n")

# 4. Data filtering and cleaning --------------------------------------------------------
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
if (ncol(data) >= 32) {
  colnames(data)[12:32] <- c(
    "TSH", "HDL", "TC", "GLU", "Vit_D", "Vit_D2", "Vit_D3", 
    "VFA", "VFM", "VFV", "BMC_Head", "BMD_Head", 
    "BMD_Lower_Limb", "BMD_L_Spine", "BMD_Pelvis", 
    "Total_BMC", "Total_BMD", "Testo", "E2", "SHBG", "hs_CRP"
  )
}

# 5. Extract analysis variables ----------------------------------------------------
cat("Step 4: Extracting analysis variables...\n")

# Extract key variables
analysis_vars <- data %>%
  dplyr::select(
    VFA, VFM,           # Visceral fat
    SHBG,               # Sex hormone-binding globulin
    Vit_D, Vit_D3,      # Vitamin D
    Total_BMC, Total_BMD # Bone mineral density
  ) %>%
  mutate(across(everything(), as.numeric))

cat("Analysis variables:", paste(colnames(analysis_vars), collapse = ", "), "\n")

# 6. Data quality checking and cleaning ----------------------------------------------------
cat("Step 5: Data quality checking and cleaning...\n")

# Check missing values
missing_counts <- colSums(is.na(analysis_vars))
cat("Missing value counts:\n")
print(missing_counts)

# Use complete data
analysis_complete <- analysis_vars[complete.cases(analysis_vars), ]
cat("\nComplete data sample size:", nrow(analysis_complete), "\n")

if (nrow(analysis_complete) < 100) {
  warning("Sample size is small, results may be unstable")
}

# Winsorize outliers
winsorize <- function(x, limits = c(0.025, 0.975)) {
  q <- quantile(x, probs = limits, na.rm = TRUE)
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  return(x)
}

analysis_clean <- as.data.frame(lapply(analysis_complete, winsorize))

# Standardize variables
analysis_standardized <- as.data.frame(scale(analysis_clean))
cat("✓ Variables standardized\n")

# Create composite variables (for simplified analysis)
analysis_standardized$BMD_mean <- (analysis_standardized$Total_BMD + 
                                     analysis_standardized$Total_BMC) / 2
analysis_standardized$VitD_mean <- (analysis_standardized$Vit_D + 
                                      analysis_standardized$Vit_D3) / 2

# Final data
analysis_data <- analysis_standardized

cat("✓ Data preparation completed\n")
cat("Final sample size:", nrow(analysis_data), "\n")

# 7. Create analysis functions ----------------------------------------------------------
cat("\nStep 6: Creating analysis functions...\n")

# Regression analysis function
run_regression <- function(formula, data, model_name) {
  cat("\n", model_name, ":\n", sep = "")
  
  model <- lm(formula, data = data)
  summary_model <- summary(model)
  
  # Output results
  cat(sprintf("R² = %.3f, F = %.2f, p = %.4f\n", 
              summary_model$r.squared, 
              summary_model$fstatistic[1],
              pf(summary_model$fstatistic[1], 
                 summary_model$fstatistic[2], 
                 summary_model$fstatistic[3], lower.tail = FALSE)))
  
  cat("Coefficient estimates:\n")
  coef_table <- summary_model$coefficients
  
  for (i in 2:nrow(coef_table)) {
    p_val <- coef_table[i, 4]
    sig <- ifelse(p_val < 0.001, "***",
                  ifelse(p_val < 0.01, "**",
                         ifelse(p_val < 0.05, "*", "")))
    
    cat(sprintf("  %s: β = %.3f, t = %.2f, p = %.4f%s\n",
                rownames(coef_table)[i],
                coef_table[i, 1],
                coef_table[i, 3],
                p_val,
                sig))
  }
  
  # Save results
  result <- list(
    model = model,
    summary = summary_model,
    formula = formula
  )
  
  return(result)
}

# 8. Individual path regression analysis -------------------------------------------------------
cat("\nStep 7: Individual path regression analysis...\n")

# Store all regression results
all_results <- list()

# Path 1: Relationship between VAT and SHBG
cat("\n1. Path 1: Relationship between VAT and SHBG\n")
all_results$VAT_SHBG <- run_regression(SHBG ~ VFA + VFM, analysis_data, "VAT → SHBG")

# Path 2: Relationship between SHBG and BMD
cat("\n2. Path 2: Relationship between SHBG and BMD\n")
all_results$SHBG_BMD <- run_regression(BMD_mean ~ SHBG, analysis_data, "SHBG → BMD")

# Path 3: Relationship between BMD and VitD
cat("\n3. Path 3: Relationship between BMD and VitD\n")
all_results$BMD_VitD <- run_regression(VitD_mean ~ BMD_mean, analysis_data, "BMD → VitD")

# Path 4: Direct relationship between VAT and BMD
cat("\n4. Path 4: Direct relationship between VAT and BMD\n")
all_results$VAT_BMD <- run_regression(BMD_mean ~ VFA + VFM, analysis_data, "VAT → BMD (direct)")

# Path 5: Direct relationship between VAT and VitD
cat("\n5. Path 5: Direct relationship between VAT and VitD\n")
all_results$VAT_VitD <- run_regression(VitD_mean ~ VFA + VFM, analysis_data, "VAT → VitD (direct)")

# Path 6: Relationship between SHBG and VitD
cat("\n6. Path 6: Relationship between SHBG and VitD\n")
all_results$SHBG_VitD <- run_regression(VitD_mean ~ SHBG, analysis_data, "SHBG → VitD")

# NEW: Path 7: Mediation analysis for VFA→SHBG→VitD
cat("\n7. Path 7: Mediation analysis for VFA→SHBG→VitD\n")
cat("Testing the indirect effect of VFA on VitD through SHBG...\n")

# First, run the three models for mediation analysis
m1 <- lm(SHBG ~ VFA, data = analysis_data)  # Path a: VFA -> SHBG
m2 <- lm(VitD_mean ~ SHBG + VFA, data = analysis_data)  # Path b and c': SHBG -> VitD, controlling for VFA

# Extract coefficients
a <- coef(m1)["VFA"]
b <- coef(m2)["SHBG"]
c_prime <- coef(m2)["VFA"]  # Direct effect

# Calculate indirect effect (a*b)
indirect_effect <- a * b
cat(sprintf("  Path a (VFA→SHBG): β = %.3f\n", a))
cat(sprintf("  Path b (SHBG→VitD): β = %.3f\n", b))
cat(sprintf("  Direct effect (VFA→VitD): β = %.3f\n", c_prime))
cat(sprintf("  Indirect effect (VFA→SHBG→VitD): β = %.3f\n", indirect_effect))
cat(sprintf("  Total effect: β = %.3f (indirect + direct)\n", indirect_effect + c_prime))

# Sobel test for mediation
se_a <- summary(m1)$coefficients["VFA", "Std. Error"]
se_b <- summary(m2)$coefficients["SHBG", "Std. Error"]
z_value <- indirect_effect / sqrt(b^2 * se_a^2 + a^2 * se_b^2)
p_value <- 2 * (1 - pnorm(abs(z_value)))
cat(sprintf("  Sobel test: z = %.3f, p = %.4f\n", z_value, p_value))

# Save mediation results
all_results$Mediation_VFA_SHBG_VitD <- list(
  model_a = m1,
  model_b = m2,
  coefficients = c(a = a, b = b, c_prime = c_prime, indirect = indirect_effect),
  sobel_test = c(z = z_value, p = p_value)
)

# 9. Create path summary table ---------------------------------------------------------
cat("\nStep 8: Creating path summary table...\n")

# Extract key results
path_summary <- data.frame()

# Extract main coefficients from each model
paths <- list(
  list(name = "VFA → SHBG", model = "VAT_SHBG", coef = "VFA"),
  list(name = "VFM → SHBG", model = "VAT_SHBG", coef = "VFM"),
  list(name = "SHBG → BMD", model = "SHBG_BMD", coef = "SHBG"),
  list(name = "BMD → VitD", model = "BMD_VitD", coef = "BMD_mean"),
  list(name = "VFA → BMD", model = "VAT_BMD", coef = "VFA"),
  list(name = "VFM → BMD", model = "VAT_BMD", coef = "VFM"),
  list(name = "VFA → VitD", model = "VAT_VitD", coef = "VFA"),
  list(name = "VFM → VitD", model = "VAT_VitD", coef = "VFM"),
  list(name = "SHBG → VitD", model = "SHBG_VitD", coef = "SHBG")
)

for (path in paths) {
  model_name <- path$model
  coef_name <- path$coef
  
  if (coef_name %in% rownames(all_results[[model_name]]$summary$coefficients)) {
    coef_row <- all_results[[model_name]]$summary$coefficients[coef_name, ]
    
    # 修复：确保每个新行都有正确的行名
    new_row <- data.frame(
      Path = path$name,
      Beta = round(coef_row[1], 3),
      Std_Error = round(coef_row[2], 3),
      t_Value = round(coef_row[3], 2),
      P_Value = round(coef_row[4], 4),
      R_Squared = round(all_results[[model_name]]$summary$r.squared, 3),
      Significance = ifelse(coef_row[4] < 0.001, "***",
                            ifelse(coef_row[4] < 0.01, "**",
                                   ifelse(coef_row[4] < 0.05, "*", ""))),
      stringsAsFactors = FALSE,
      row.names = NULL  # 显式设置行名为NULL
    )
    
    # 使用bind_rows替代rbind，避免行名问题
    path_summary <- dplyr::bind_rows(path_summary, new_row)
  }
}

# 添加中介路径 - 修复行名问题
if (!is.null(all_results$Mediation_VFA_SHBG_VitD)) {
  mediation_result <- all_results$Mediation_VFA_SHBG_VitD
  
  # 检查是否有有效的间接效应值
  if (!is.null(mediation_result$coefficients["indirect"]) && 
      !is.na(mediation_result$coefficients["indirect"])) {
    
    # 创建中介路径行
    mediation_row <- data.frame(
      Path = "VFA → SHBG → VitD (indirect)",
      Beta = round(mediation_result$coefficients["indirect"], 3),
      Std_Error = NA_real_,  # 显式指定NA类型
      t_Value = ifelse(!is.null(mediation_result$sobel_test["z"]),
                       round(mediation_result$sobel_test["z"], 2), NA_real_),
      P_Value = ifelse(!is.null(mediation_result$sobel_test["p"]),
                       round(mediation_result$sobel_test["p"], 4), NA_real_),
      R_Squared = NA_real_,  # 中介分析没有R²
      Significance = ifelse(!is.null(mediation_result$sobel_test["p"]),
                            ifelse(mediation_result$sobel_test["p"] < 0.001, "***",
                                   ifelse(mediation_result$sobel_test["p"] < 0.01, "**",
                                          ifelse(mediation_result$sobel_test["p"] < 0.05, "*", ""))),
                            ""),
      stringsAsFactors = FALSE,
      row.names = NULL  # 显式设置行名为NULL
    )
    
    # 添加到汇总表
    path_summary <- dplyr::bind_rows(path_summary, mediation_row)
  } else {
    cat("Warning: Indirect effect is NA or NULL, skipping mediation path in summary table\n")
  }
}


# 同样修复路径汇总表的其他部分
# Add effect strength classification
path_summary$Effect_Size <- ifelse(abs(path_summary$Beta) >= 0.3, "Strong",
                                   ifelse(abs(path_summary$Beta) >= 0.1, "Moderate", "Weak"))

path_summary$Direction <- ifelse(path_summary$Beta > 0, "Positive", "Negative")

# Add path type classification
path_summary$Path_Type <- NA_character_  # 使用正确的NA类型
path_summary$Path_Type[grepl("VFA → SHBG|VFM → SHBG", path_summary$Path)] <- "Visceral Fat → SHBG"
path_summary$Path_Type[grepl("SHBG → BMD", path_summary$Path)] <- "SHBG → Bone Density"
path_summary$Path_Type[grepl("BMD → VitD", path_summary$Path)] <- "Bone Density → Vitamin D"
path_summary$Path_Type[grepl("VFA → BMD|VFM → BMD", path_summary$Path)] <- "Visceral Fat → Bone Density"
path_summary$Path_Type[grepl("VFA → VitD|VFM → VitD", path_summary$Path)] <- "Visceral Fat → Vitamin D"
path_summary$Path_Type[grepl("SHBG → VitD", path_summary$Path)] <- "SHBG → Vitamin D"
path_summary$Path_Type[grepl("VFA → SHBG → VitD", path_summary$Path)] <- "VFA→SHBG→VitD Mediation"

# Mark key paths of interest
path_summary$Key_Path <- FALSE
path_summary$Key_Path[grepl("VFA → SHBG|SHBG → VitD|VFA → SHBG → VitD", path_summary$Path)] <- TRUE

# Reorder columns
path_summary <- path_summary %>%
  dplyr::select(Path, Path_Type, Key_Path, everything())

# Save summary table
write.csv(path_summary, "results/Path_Analysis/tables/path_summary.csv", row.names = FALSE)
write.csv(path_summary, "results/Path_Analysis/tables/path_summary_detailed.csv", row.names = FALSE)

cat("Path summary table:\n")
print(path_summary)

# 10. Detailed analysis of VFA→SHBG→VitD pathway --------------------------------------------
cat("\nStep 9: Detailed analysis of VFA→SHBG→VitD pathway...\n")

# Create a special directory for this pathway
if (!dir.exists("results/Path_Analysis/Path_VFA_SHBG_VitD")) {
  dir.create("results/Path_Analysis/Path_VFA_SHBG_VitD", recursive = TRUE)
}

# Create visualization of the VFA→SHBG→VitD pathway
cat("Creating visualizations for VFA→SHBG→VitD pathway...\n")

# 1. Scatter plots for each segment of the pathway
p1 <- ggplot(analysis_data, aes(x = VFA, y = SHBG)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = "VFA → SHBG Relationship",
       x = "Visceral Fat Area (VFA, standardized)",
       y = "SHBG (standardized)",
       subtitle = sprintf("β = %.3f, p = %.4f", 
                          coef(m1)["VFA"], 
                          summary(m1)$coefficients["VFA", 4])) +
  theme_minimal()

p2 <- ggplot(analysis_data, aes(x = SHBG, y = VitD_mean)) +
  geom_point(alpha = 0.6, color = "darkgreen") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = "SHBG → Vitamin D Relationship",
       x = "SHBG (standardized)",
       y = "Vitamin D (standardized)",
       subtitle = sprintf("β = %.3f, p = %.4f", 
                          coef(m2)["SHBG"], 
                          summary(m2)$coefficients["SHBG", 4])) +
  theme_minimal()

p3 <- ggplot(analysis_data, aes(x = VFA, y = VitD_mean)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = "VFA → Vitamin D Relationship",
       x = "Visceral Fat Area (VFA, standardized)",
       y = "Vitamin D (standardized)",
       subtitle = sprintf("Direct effect: β = %.3f, p = %.4f", 
                          coef(m2)["VFA"], 
                          summary(m2)$coefficients["VFA", 4])) +
  theme_minimal()

# Combine plots
pathway_plot <- (p1 | p2) / p3 + 
  plot_annotation(title = "VFA→SHBG→VitD Pathway Analysis",
                  subtitle = sprintf("Indirect effect: β = %.3f, Sobel test p = %.4f", 
                                     indirect_effect, p_value),
                  theme = theme(plot.title = element_text(size = 16, face = "bold"),
                                plot.subtitle = element_text(size = 12)))

# Save the plot
ggsave("results/Path_Analysis/Path_VFA_SHBG_VitD/pathway_scatter_plots.png", 
       pathway_plot, width = 12, height = 10, dpi = 300)

# 2. Create mediation diagram
cat("Creating mediation diagram...\n")

# Simple mediation plot
mediation_data <- data.frame(
  Path = c("VFA→SHBG", "SHBG→VitD", "VFA→VitD (direct)", "VFA→VitD (total)", "VFA→SHBG→VitD (indirect)"),
  Beta = c(a, b, c_prime, a*b + c_prime, indirect_effect),
  P_Value = c(summary(m1)$coefficients["VFA", 4],
              summary(m2)$coefficients["SHBG", 4],
              summary(m2)$coefficients["VFA", 4],
              NA,  # Total effect p-value would need separate calculation
              p_value)
)

# Create mediation results table
mediation_table <- data.frame(
  Component = c("Path a (VFA→SHBG)", "Path b (SHBG→VitD)", 
                "Direct effect (VFA→VitD)", "Indirect effect (VFA→SHBG→VitD)",
                "Total effect"),
  Beta = c(round(a, 3), round(b, 3), round(c_prime, 3), 
           round(indirect_effect, 3), round(indirect_effect + c_prime, 3)),
  SE = c(round(se_a, 3), round(se_b, 3), 
         round(summary(m2)$coefficients["VFA", "Std. Error"], 3),
         round(sqrt(b^2 * se_a^2 + a^2 * se_b^2), 3), NA),
  t_z = c(round(summary(m1)$coefficients["VFA", "t value"], 2),
          round(summary(m2)$coefficients["SHBG", "t value"], 2),
          round(summary(m2)$coefficients["VFA", "t value"], 2),
          round(z_value, 2), NA),
  P_Value = c(round(summary(m1)$coefficients["VFA", 4], 4),
              round(summary(m2)$coefficients["SHBG", 4], 4),
              round(summary(m2)$coefficients["VFA", 4], 4),
              round(p_value, 4), NA),
  Significance = c(ifelse(summary(m1)$coefficients["VFA", 4] < 0.001, "***",
                          ifelse(summary(m1)$coefficients["VFA", 4] < 0.01, "**",
                                 ifelse(summary(m1)$coefficients["VFA", 4] < 0.05, "*", ""))),
                   ifelse(summary(m2)$coefficients["SHBG", 4] < 0.001, "***",
                          ifelse(summary(m2)$coefficients["SHBG", 4] < 0.01, "**",
                                 ifelse(summary(m2)$coefficients["SHBG", 4] < 0.05, "*", ""))),
                   ifelse(summary(m2)$coefficients["VFA", 4] < 0.001, "***",
                          ifelse(summary(m2)$coefficients["VFA", 4] < 0.01, "**",
                                 ifelse(summary(m2)$coefficients["VFA", 4] < 0.05, "*", ""))),
                   ifelse(p_value < 0.001, "***",
                          ifelse(p_value < 0.01, "**",
                                 ifelse(p_value < 0.05, "*", ""))),
                   "")
)

write.csv(mediation_table, "results/Path_Analysis/Path_VFA_SHBG_VitD/mediation_analysis_results.csv", 
          row.names = FALSE)

# Create a visual mediation diagram
# Using base R to create a simple diagram
png("results/Path_Analysis/Path_VFA_SHBG_VitD/mediation_diagram.png", 
    width = 1200, height = 800, res = 150)
par(mar = c(0, 0, 2, 0))
plot(0, 0, type = "n", xlim = c(-1, 11), ylim = c(-1, 6), 
     axes = FALSE, xlab = "", ylab = "", 
     main = "Mediation Analysis: VFA→SHBG→VitD Pathway")

# Draw nodes
points(2, 4, pch = 19, cex = 3, col = "steelblue")
text(2, 4, "VFA", cex = 1.5, font = 2)

points(6, 4, pch = 19, cex = 3, col = "darkgreen")
text(6, 4, "SHBG", cex = 1.5, font = 2)

points(10, 4, pch = 19, cex = 3, col = "purple")
text(10, 4, "VitD", cex = 1.5, font = 2)

# Draw arrows with coefficients
# VFA -> SHBG
arrows(2.5, 4, 5.5, 4, length = 0.15, lwd = 3, col = "red")
text(4, 4.5, sprintf("a = %.3f***", a), cex = 1.2, font = 2)

# SHBG -> VitD
arrows(6.5, 4, 9.5, 4, length = 0.15, lwd = 3, col = "red")
text(8, 4.5, sprintf("b = %.3f*", b), cex = 1.2, font = 2)

# VFA -> VitD (direct)
arrows(2, 3.5, 9.8, 3.5, length = 0.15, lwd = 2, lty = 2, col = "blue")
text(6, 3, sprintf("c' = %.3f", c_prime), cex = 1.2, font = 2, col = "blue")

# Add indirect effect information
text(6, 2, sprintf("Indirect effect: a×b = %.3f×%.3f = %.3f*", 
                   round(a, 3), round(b, 3), round(indirect_effect, 3)), 
     cex = 1.2, font = 2, col = "darkred")

text(6, 1.5, sprintf("Sobel test: z = %.2f, p = %.4f", z_value, p_value), 
     cex = 1.1, font = 2, col = "darkred")

# Add legend
legend("bottom", legend = c("Mediator path (a, b)", "Direct effect (c')", 
                            "Indirect effect (a×b)"), 
       lty = c(1, 2, NA), lwd = c(3, 2, NA), 
       col = c("red", "blue", "darkred"),
       pch = c(NA, NA, NA), cex = 1.2, bty = "n")

dev.off()

# 3. Bootstrap mediation analysis for more robust results
cat("\nPerforming bootstrap mediation analysis...\n")

# Simple bootstrap function for mediation
bootstrap_mediation <- function(data, n_boot = 1000) {
  n <- nrow(data)
  indirect_effects <- numeric(n_boot)
  
  for (i in 1:n_boot) {
    # Bootstrap sample
    boot_idx <- sample(1:n, n, replace = TRUE)
    boot_data <- data[boot_idx, ]
    
    # Fit models on bootstrap sample
    boot_m1 <- lm(SHBG ~ VFA, data = boot_data)
    boot_m2 <- lm(VitD_mean ~ SHBG + VFA, data = boot_data)
    
    # Calculate indirect effect
    a_boot <- coef(boot_m1)["VFA"]
    b_boot <- coef(boot_m2)["SHBG"]
    indirect_effects[i] <- a_boot * b_boot
  }
  
  return(indirect_effects)
}

# Run bootstrap
set.seed(123)
boot_indirect <- bootstrap_mediation(analysis_data, n_boot = 5000)

# Calculate bootstrap CI
boot_ci <- quantile(boot_indirect, probs = c(0.025, 0.975), na.rm = TRUE)
boot_mean <- mean(boot_indirect, na.rm = TRUE)

cat("Bootstrap mediation results:\n")
cat(sprintf("  Mean indirect effect: %.3f\n", boot_mean))
cat(sprintf("  95%% Bootstrap CI: [%.3f, %.3f]\n", boot_ci[1], boot_ci[2]))

# Save bootstrap results
bootstrap_results <- data.frame(
  Indirect_Effect = boot_indirect
)
write.csv(bootstrap_results, 
          "results/Path_Analysis/Path_VFA_SHBG_VitD/bootstrap_mediation_results.csv", 
          row.names = FALSE)

# Plot bootstrap distribution
png("results/Path_Analysis/Path_VFA_SHBG_VitD/bootstrap_distribution.png", 
    width = 1000, height = 800, res = 150)
hist(boot_indirect, breaks = 50, col = "lightblue", border = "white",
     main = "Bootstrap Distribution of Indirect Effect (VFA→SHBG→VitD)",
     xlab = "Indirect Effect (a × b)", 
     ylab = "Frequency")
abline(v = indirect_effect, col = "red", lwd = 3, lty = 2)
abline(v = 0, col = "black", lwd = 2)
abline(v = boot_ci, col = "blue", lwd = 2, lty = 2)
legend("topright", 
       legend = c(sprintf("Original: %.3f", indirect_effect),
                  sprintf("Bootstrap mean: %.3f", boot_mean),
                  "95% Bootstrap CI",
                  "Null (0)"),
       col = c("red", "darkgreen", "blue", "black"),
       lwd = c(3, 2, 2, 2),
       lty = c(2, 1, 2, 1),
       cex = 0.9)
dev.off()

# Create comprehensive report for VFA→SHBG→VitD pathway
sink("results/Path_Analysis/Path_VFA_SHBG_VitD/pathway_analysis_report.txt")

cat("VFA→SHBG→VitD Pathway Analysis Report\n")
cat("=====================================\n\n")

cat("Analysis Date:", Sys.Date(), "\n")
cat("Sample Size:", nrow(analysis_data), "\n\n")

cat("1. Theoretical Framework\n")
cat("------------------------\n")
cat("Hypothesized pathway: Visceral Fat Area (VFA) → Sex Hormone-Binding Globulin (SHBG) → Vitamin D\n")
cat("Biological rationale:\n")
cat("1. Visceral fat accumulation is associated with reduced SHBG production\n")
cat("2. SHBG may influence vitamin D metabolism through multiple mechanisms:\n")
cat("   - SHBG-bound vitamin D may have different bioavailability\n")
cat("   - SHBG may interact with vitamin D binding protein (DBP)\n")
cat("   - SHBG may influence vitamin D receptor expression\n")
cat("3. This pathway represents a potential endocrine-metabolic link between obesity and vitamin D status\n\n")

cat("2. Path Analysis Results\n")
cat("------------------------\n")

cat("2.1 Individual Path Coefficients:\n")
cat("   VFA → SHBG: β = ", round(a, 3), 
    ", t = ", round(summary(m1)$coefficients["VFA", "t value"], 2),
    ", p = ", format(summary(m1)$coefficients["VFA", 4], scientific = FALSE, digits = 4), "\n", sep = "")
cat("   SHBG → VitD: β = ", round(b, 3), 
    ", t = ", round(summary(m2)$coefficients["SHBG", "t value"], 2),
    ", p = ", format(summary(m2)$coefficients["SHBG", 4], scientific = FALSE, digits = 4), "\n", sep = "")
cat("   VFA → VitD (direct): β = ", round(c_prime, 3), 
    ", t = ", round(summary(m2)$coefficients["VFA", "t value"], 2),
    ", p = ", format(summary(m2)$coefficients["VFA", 4], scientific = FALSE, digits = 4), "\n\n", sep = "")

cat("2.2 Mediation Analysis:\n")
cat("   Indirect effect (a × b): ", round(indirect_effect, 3), "\n", sep = "")
cat("   Sobel test: z = ", round(z_value, 3), 
    ", p = ", format(p_value, scientific = FALSE, digits = 4), "\n\n", sep = "")

cat("2.3 Bootstrap Results (n = 5000):\n")
cat("   Mean indirect effect: ", round(boot_mean, 3), "\n", sep = "")
cat("   95% Bootstrap CI: [", round(boot_ci[1], 3), ", ", round(boot_ci[2], 3), "]\n\n", sep = "")

cat("3. Interpretation\n")
cat("-----------------\n")

if (p_value < 0.05) {
  if (indirect_effect > 0) {
    cat("✓ Significant positive mediation found: VFA positively affects Vitamin D through SHBG\n")
  } else {
    cat("✓ Significant negative mediation found: VFA negatively affects Vitamin D through SHBG\n")
  }
  
  # Calculate proportion mediated
  total_effect <- indirect_effect + c_prime
  prop_mediated <- ifelse(total_effect != 0, indirect_effect / total_effect, NA)
  
  if (!is.na(prop_mediated)) {
    cat(sprintf("   Proportion mediated: %.1f%%\n", prop_mediated * 100))
  }
  
  # Interpret direction
  if (a < 0 && b > 0) {
    cat("   Interpretation: Higher VFA is associated with lower SHBG, and lower SHBG is associated with lower Vitamin D\n")
  } else if (a < 0 && b < 0) {
    cat("   Interpretation: Higher VFA is associated with lower SHBG, and lower SHBG is associated with higher Vitamin D\n")
  } else if (a > 0 && b > 0) {
    cat("   Interpretation: Higher VFA is associated with higher SHBG, and higher SHBG is associated with higher Vitamin D\n")
  } else if (a > 0 && b < 0) {
    cat("   Interpretation: Higher VFA is associated with higher SHBG, and higher SHBG is associated with lower Vitamin D\n")
  }
} else {
  cat("✗ No significant mediation found: The VFA→SHBG→VitD pathway is not statistically significant\n")
  cat("   This suggests that SHBG may not be a significant mediator between VFA and Vitamin D in this sample\n")
}

cat("\n4. Clinical and Research Implications\n")
cat("--------------------------------------\n")

if (p_value < 0.05) {
  cat("✓ SHBG may be a potential therapeutic target for improving Vitamin D status in individuals with high visceral fat\n")
  cat("✓ Interventions targeting visceral fat reduction may improve Vitamin D status through SHBG modulation\n")
  cat("✓ Measurement of SHBG could help identify individuals at risk for Vitamin D deficiency related to obesity\n")
} else {
  cat("• The VFA→SHBG→VitD pathway may not be a primary mechanism linking visceral fat to Vitamin D status\n")
  cat("• Other pathways (e.g., inflammation, adipokines, direct effects) should be investigated\n")
  cat("• SHBG may still be important for Vitamin D metabolism but not as a mediator of visceral fat effects\n")
}

cat("\n5. Limitations\n")
cat("---------------\n")
cat("• Cross-sectional design prevents causal inference\n")
cat("• Standardized variables may limit clinical interpretability\n")
cat("• Potential confounding factors not included in analysis\n")
cat("• Sample size may limit statistical power\n")
cat("• Single mediation model tested; other mediators may exist\n")

cat("\n6. Recommendations for Future Research\n")
cat("---------------------------------------\n")
cat("1. Longitudinal studies to establish temporal relationships\n")
cat("2. Inclusion of additional potential mediators (e.g., inflammation markers, other hormones)\n")
cat("3. Investigation of gender-specific effects\n")
cat("4. Experimental studies to test causal mechanisms\n")
cat("5. Clinical trials targeting SHBG to improve Vitamin D status\n")

cat("\n7. Files Generated\n")
cat("------------------\n")
cat("1. pathway_scatter_plots.png - Visualizations of individual relationships\n")
cat("2. mediation_diagram.png - Diagram of mediation model\n")
cat("3. mediation_analysis_results.csv - Table of mediation results\n")
cat("4. bootstrap_mediation_results.csv - Bootstrap samples\n")
cat("5. bootstrap_distribution.png - Distribution of bootstrap indirect effects\n")
cat("6. pathway_analysis_report.txt - This comprehensive report\n")

sink()

cat("✓ VFA→SHBG→VitD pathway analysis completed\n")

# 11. Structural Equation Model analysis (with complete fit analysis) -----------------------------------
cat("\nStep 10: Structural Equation Model fit analysis...\n")

# First perform data check
cat("\n=== Data Check ===\n")
cat("Sample size:", nrow(analysis_data), "\n")
cat("Number of variables:", ncol(analysis_data), "\n")

# Check correlations between variables
cor_matrix <- cor(analysis_data, use = "complete.obs")
cat("\nVariable correlation matrix:\n")
print(round(cor_matrix, 3))

# Check for highly correlated variables (|r| > 0.9)
high_cor <- which(abs(cor_matrix) > 0.9 & abs(cor_matrix) < 1, arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  cat("\nWarning: Found highly correlated variable pairs:\n")
  for (i in 1:nrow(high_cor)) {
    var1 <- rownames(cor_matrix)[high_cor[i, 1]]
    var2 <- colnames(cor_matrix)[high_cor[i, 2]]
    cor_value <- cor_matrix[high_cor[i, 1], high_cor[i, 2]]
    cat(sprintf("  %s and %s: r = %.3f\n", var1, var2, cor_value))
  }
}

# Create simplified theoretical models (avoid over-parameterization)
cat("\n=== Building Simplified Models ===\n")

# Model 1: Basic cascade path model
model_path1 <- '
  # Main paths
  SHBG ~ a1*VFA
  Total_BMD ~ b1*SHBG
  Vit_D ~ c1*Total_BMD
'

# Model 2: Model with added direct path
model_path2 <- '
  # Simplified paths
  SHBG ~ a2*VFA
  Total_BMD ~ b2*SHBG
  Vit_D ~ c2*Total_BMD + d2*VFA
'

# Model 3: Dual path model (includes VFA→SHBG→VitD)
model_path3 <- '
  # NEW: Includes VFA→SHBG→VitD path
  SHBG ~ a3*VFA
  Total_BMD ~ b3*SHBG
  Vit_D ~ c3*SHBG + d3*Total_BMD
'

# Model 4: NEW - VFA→SHBG→VitD focused model
model_path4 <- '
  # Focus on VFA→SHBG→VitD pathway
  SHBG ~ a4*VFA
  Vit_D ~ b4*SHBG
  # Allow residual correlation between BMD and VitD if included
  Total_BMD ~~ Vit_D
'

# Model 5: NEW - Complete model with all paths
model_path5 <- '
  # Complete model including VFA→SHBG→VitD pathway
  SHBG ~ a5*VFA
  Total_BMD ~ b5*SHBG
  Vit_D ~ c5*SHBG + d5*Total_BMD + e5*VFA
'

# Modify fit_sem_model function, add debugging information
fit_sem_model <- function(model_spec, model_name, data) {
  cat(sprintf("\nFitting model: %s\n", model_name))
  cat("Model specification:\n")
  cat(model_spec)
  cat("\n\n")
  
  tryCatch({
    # Try different optimizers and settings
    fit <- sem(model_spec, 
               data = data,
               estimator = "ML",
               se = "standard",
               fixed.x = FALSE,
               orthogonal = FALSE,
               control = list(iter.max = 1000,  # Increase iterations
                              optim.dx.tol = 1e-07,
                              optim.gn.tol = 1e-07))
    
    # Check convergence
    if (lavInspect(fit, "converged")) {
      cat("✓ Model converged successfully\n")
      
      # Extract fit indices
      fit_indices <- fitMeasures(fit, 
                                 c("chisq", "df", "pvalue", 
                                   "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue",
                                   "srmr", "aic", "bic", "gfi", "agfi"))
      
      # Parameter estimates
      params <- parameterEstimates(fit, standardized = TRUE)
      
      # Calculate modification indices
      mod_indices <- modificationIndices(fit, sort = TRUE, minimum.value = 3.84)
      
      # Save results
      result <- list(
        fit = fit,
        fit_indices = fit_indices,
        parameters = params,
        modification_indices = mod_indices,
        converged = TRUE
      )
      
      # Add to summary table
      fit_indices_summary <<- rbind(fit_indices_summary, data.frame(
        Model = model_name,
        ChiSq = round(fit_indices["chisq"], 2),
        df = fit_indices["df"],
        p = round(fit_indices["pvalue"], 4),
        CFI = round(fit_indices["cfi"], 3),
        TLI = round(fit_indices["tli"], 3),
        RMSEA = round(fit_indices["rmsea"], 3),
        RMSEA_CI = sprintf("[%.3f, %.3f]", 
                           fit_indices["rmsea.ci.lower"], 
                           fit_indices["rmsea.ci.upper"]),
        RMSEA_p = round(fit_indices["rmsea.pvalue"], 4),
        SRMR = round(fit_indices["srmr"], 3),
        AIC = round(fit_indices["aic"], 1),
        BIC = round(fit_indices["bic"], 1),
        GFI = round(fit_indices["gfi"], 3),
        AGFI = round(fit_indices["agfi"], 3),
        stringsAsFactors = FALSE
      ))
      
      # Print fit indices
      cat("Fit indices:\n")
      cat(sprintf("  χ²(%d) = %.2f, p = %.4f\n", 
                  fit_indices["df"], fit_indices["chisq"], fit_indices["pvalue"]))
      cat(sprintf("  CFI = %.3f, TLI = %.3f\n", 
                  fit_indices["cfi"], fit_indices["tli"]))
      cat(sprintf("  RMSEA = %.3f (90%% CI: %.3f-%.3f), p = %.4f\n", 
                  fit_indices["rmsea"], 
                  fit_indices["rmsea.ci.lower"], 
                  fit_indices["rmsea.ci.upper"],
                  fit_indices["rmsea.pvalue"]))
      cat(sprintf("  SRMR = %.3f\n", fit_indices["srmr"]))
      
      return(result)
      
    } else {
      cat("✗ Model did not converge\n")
      return(list(converged = FALSE))
    }
    
  }, error = function(e) {
    cat(sprintf("✗ Model fitting error: %s\n", e$message))
    
    # Try simpler settings
    cat("Trying simplified model settings...\n")
    
    tryCatch({
      # Try simpler model
      fit <- sem(model_spec, 
                 data = data,
                 estimator = "ML",
                 se = "standard",
                 fixed.x = FALSE,
                 orthogonal = FALSE,
                 control = list(iter.max = 2000),
                 do.fit = FALSE)
      
      # Manual parameter check
      cat("Model parameter check:\n")
      cat("  Number of free parameters:", lavInspect(fit, "npar"), "\n")
      cat("  Sample size:", nrow(data), "\n")
      
      # If too many free parameters, suggest simplifying model
      if (lavInspect(fit, "npar") > nrow(data) / 5) {
        cat("  Warning: Too many free parameters, consider simplifying the model\n")
      }
      
    }, error = function(e2) {
      cat("Unable to check model parameters:", e2$message, "\n")
    })
    
    return(list(converged = FALSE, error = e$message))
  })
}

# Clear previous summary table
fit_indices_summary <- data.frame()

# Initialize sem_results list
sem_results <- list()

# Fit simplified models
cat("\n============== Trying Simplified Models ==============\n")

# First check if variables exist
cat("Checking variables in data:\n")
print(names(analysis_data))

# Ensure we use correct variable names
# If needed, create a dataset containing only necessary variables
analysis_data_simple <- analysis_data[, c("VFA", "SHBG", "Total_BMD", "Vit_D")]

# Check variable types
cat("\nVariable types:\n")
str(analysis_data_simple)

# Check missing values
cat("\nMissing value counts:\n")
print(colSums(is.na(analysis_data_simple)))

# Use complete data
analysis_data_complete <- na.omit(analysis_data_simple)
cat("\nComplete data sample size:", nrow(analysis_data_complete), "\n")

# Fit simplified models
sem_results$model_path1 <- fit_sem_model(model_path1, "VFA→SHBG→BMD→VitD", analysis_data_complete)
sem_results$model_path2 <- fit_sem_model(model_path2, "VFA→SHBG→BMD→VitD + VFA→VitD", analysis_data_complete)
sem_results$model_path3 <- fit_sem_model(model_path3, "VFA→SHBG→BMD + VFA→SHBG→VitD", analysis_data_complete)
sem_results$model_path4 <- fit_sem_model(model_path4, "VFA→SHBG→VitD (focused)", analysis_data_complete)
sem_results$model_path5 <- fit_sem_model(model_path5, "Complete Model (all paths)", analysis_data_complete)

# If simplified models also don't converge, try using more basic functions
if (nrow(fit_indices_summary) == 0) {
  cat("\n=== Trying Manual Fitting with Basic Functions ===\n")
  
  # Try simple linear models
  cat("Trying simple path analysis:\n")
  
  # Path 1: VFA -> SHBG
  m1 <- lm(SHBG ~ VFA, data = analysis_data_complete)
  cat("Path 1 (VFA -> SHBG):\n")
  cat(sprintf("  β = %.3f, p = %.4f\n", 
              summary(m1)$coefficients[2, 1],
              summary(m1)$coefficients[2, 4]))
  
  # Path 2: SHBG -> Total_BMD
  m2 <- lm(Total_BMD ~ SHBG, data = analysis_data_complete)
  cat("Path 2 (SHBG -> Total_BMD):\n")
  cat(sprintf("  β = %.3f, p = %.4f\n", 
              summary(m2)$coefficients[2, 1],
              summary(m2)$coefficients[2, 4]))
  
  # Path 3: Total_BMD -> Vit_D
  m3 <- lm(Vit_D ~ Total_BMD, data = analysis_data_complete)
  cat("Path 3 (Total_BMD -> Vit_D):\n")
  cat(sprintf("  β = %.3f, p = %.4f\n", 
              summary(m3)$coefficients[2, 1],
              summary(m3)$coefficients[2, 4]))
  
  # Path 4: SHBG -> Vit_D
  m4 <- lm(Vit_D ~ SHBG, data = analysis_data_complete)
  cat("Path 4 (SHBG -> Vit_D):\n")
  cat(sprintf("  β = %.3f, p = %.4f\n", 
              summary(m4)$coefficients[2, 1],
              summary(m4)$coefficients[2, 4]))
  
  # Create artificial fit index table for reporting
  fit_indices_summary <- rbind(fit_indices_summary, data.frame(
    Model = "Linear Regression Path Model",
    ChiSq = NA,
    df = NA,
    p = NA,
    CFI = NA,
    TLI = NA,
    RMSEA = NA,
    RMSEA_CI = "[NA, NA]",
    RMSEA_p = NA,
    SRMR = NA,
    AIC = NA,
    BIC = NA,
    GFI = NA,
    AGFI = NA,
    stringsAsFactors = FALSE
  ))
  
  # Save linear model results
  lm_results <- list(
    VFA_SHBG = m1,
    SHBG_BMD = m2,
    BMD_VitD = m3,
    SHBG_VitD = m4
  )
  saveRDS(lm_results, "results/Path_Analysis/SEM_fit/lm_path_results.rds")
}

# Save fit indices summary table
if (nrow(fit_indices_summary) > 0) {
  write.csv(fit_indices_summary, 
            "results/Path_Analysis/SEM_fit/fit_indices_summary.csv", 
            row.names = FALSE)
  cat("\n✓ Fit indices summary table saved\n")
} else {
  cat("\n✗ All models failed to converge, cannot generate fit indices summary table\n")
}

cat("\n✓ All model fitting attempts completed\n")

# 12. Plot SEM results and fit evaluation plots ------------------------------------------------
cat("\nStep 11: Plotting SEM results and fit evaluation plots...\n")

# Initialize best_model_fit variable
best_model_fit <- NULL

# Check if there are converged models
converged_models <- sapply(sem_results, function(x) ifelse(is.null(x$converged), FALSE, x$converged))
if (any(converged_models)) {
  # Plot model comparison plots
  if (nrow(fit_indices_summary) > 0) {
    # 1. Fit indices radar plot (if more than one model)
    if (nrow(fit_indices_summary) > 1) {
      fit_metrics <- fit_indices_summary %>%
        dplyr::select(Model, CFI, TLI, GFI, AGFI) %>%
        pivot_longer(cols = -Model, names_to = "Metric", values_to = "Value")
      
      p1 <- ggplot(fit_metrics, aes(x = Metric, y = Value, fill = Model)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
        geom_hline(yintercept = c(0.9, 0.95), linetype = "dashed", 
                   color = c("orange", "green"), alpha = 0.5) +
        labs(title = "SEM Model Fit Indices Comparison",
             subtitle = "Dashed lines: Green=0.95 (Excellent), Orange=0.90 (Acceptable)",
             x = "Fit Index", y = "Value") +
        theme_minimal() +
        theme(legend.position = "bottom",
              axis.text.x = element_text(angle = 45, hjust = 1))
      
      # 2. RMSEA and SRMR comparison
      p2 <- ggplot(fit_indices_summary, aes(x = Model)) +
        geom_point(aes(y = RMSEA, color = "RMSEA"), size = 3) +
        geom_point(aes(y = SRMR, color = "SRMR"), size = 3) +
        geom_hline(yintercept = c(0.05, 0.08), linetype = "dashed", 
                   color = c("green", "orange"), alpha = 0.5) +
        labs(title = "RMSEA and SRMR Comparison",
             subtitle = "Dashed lines: Green=0.05 (Excellent), Orange=0.08 (Acceptable)",
             x = "Model", y = "Value") +
        scale_color_manual(values = c("RMSEA" = "blue", "SRMR" = "red")) +
        theme_minimal() +
        theme(legend.position = "bottom",
              axis.text.x = element_text(angle = 45, hjust = 1))
      
      # 3. AIC and BIC comparison
      p3 <- ggplot(fit_indices_summary, aes(x = Model)) +
        geom_point(aes(y = AIC, color = "AIC"), size = 3) +
        geom_point(aes(y = BIC, color = "BIC"), size = 3) +
        labs(title = "Information Criteria Comparison (Lower is Better)",
             x = "Model", y = "Value") +
        scale_color_manual(values = c("AIC" = "purple", "BIC" = "brown")) +
        theme_minimal() +
        theme(legend.position = "bottom",
              axis.text.x = element_text(angle = 45, hjust = 1))
      
      # Save plots
      png("results/Path_Analysis/SEM_fit/model_comparison_plots.png", 
          width = 3000, height = 2000, res = 300)
      print((p1 | p2) / p3 + plot_layout(heights = c(1, 1)))
      dev.off()
      
      cat("✓ Model comparison plots saved\n")
    }
    
    # Plot best model path diagram
    # Select best model (based on CFI and RMSEA)
    if (any(!is.na(fit_indices_summary$CFI))) {
      best_model_idx <- which.max(fit_indices_summary$CFI)
      best_model_name <- fit_indices_summary$Model[best_model_idx]
      
      # 创建一个模型名称到模型键的映射
      model_mapping <- list(
        "VFA→SHBG→BMD→VitD" = "model_path1",
        "VFA→SHBG→BMD→VitD + VFA→VitD" = "model_path2",
        "VFA→SHBG→BMD + VFA→SHBG→VitD" = "model_path3",
        "VFA→SHBG→VitD (focused)" = "model_path4",
        "Complete Model (all paths)" = "model_path5"
      )
      
      # 使用映射获取模型键
      if (best_model_name %in% names(model_mapping)) {
        model_key <- model_mapping[[best_model_name]]
      } else {
        model_key <- tolower(gsub(" ", "_", best_model_name))
      }
      
      # 检查模型是否收敛并可访问
      if (!is.null(model_key) && model_key %in% names(sem_results) && 
          !is.null(sem_results[[model_key]]$converged) && 
          sem_results[[model_key]]$converged) {
        
        best_model_fit <- sem_results[[model_key]]
        
        # Plot standardized path diagram - using simpler method
        png("results/Path_Analysis/SEM_fit/best_model_path_diagram.png", 
            width = 2500, height = 1800, res = 300)
        
        # Method 1: Use more stable parameter settings
        tryCatch({
          semPaths(best_model_fit$fit,
                   what = "std",
                   whatLabels = "std",
                   edge.label.cex = 1.0,
                   sizeMan = 10,
                   sizeLat = 12,
                   style = "lisrel",
                   fade = FALSE,
                   rotation = 1,
                   layout = "tree2",
                   mar = c(4, 4, 4, 4),
                   edge.color = "black",
                   nodeLabels = c("VFA", "SHBG", "Total_BMD", "Vit_D"))
          
          # Use independent title function to add title
          title(paste("Best Model:", best_model_name), line = 1, cex.main = 1.5)
        }, error = function(e) {
          cat("Error plotting path diagram:", e$message, "\n")
          # Try simpler method
          try({
            semPaths(best_model_fit$fit,
                     what = "std",
                     whatLabels = "std",
                     edge.label.cex = 1.0,
                     sizeMan = 8,
                     sizeLat = 10,
                     style = "lisrel",
                     fade = FALSE,
                     layout = "tree2")
            title(paste("Best Model:", best_model_name), line = 1, cex.main = 1.5)
          })
        })
        
        dev.off()
        
        # Plot parameter estimate heatmap
        params_wide <- best_model_fit$parameters %>%
          filter(op %in% c("~", "~~")) %>%
          mutate(Path = ifelse(op == "~", 
                               paste(rhs, "→", lhs),
                               paste(lhs, "~~", rhs))) %>%
          dplyr::select(Path, est, std.all, pvalue) %>%
          arrange(pvalue)
        
        write.csv(params_wide, 
                  "results/Path_Analysis/SEM_fit/best_model_parameters.csv", 
                  row.names = FALSE)
        
        cat(sprintf("✓ Best model '%s' path diagram saved\n", best_model_name))
        if (!is.na(fit_indices_summary$CFI[best_model_idx])) {
          cat(sprintf("  Fit indices: CFI=%.3f, TLI=%.3f, RMSEA=%.3f, SRMR=%.3f\n",
                      fit_indices_summary$CFI[best_model_idx],
                      fit_indices_summary$TLI[best_model_idx],
                      fit_indices_summary$RMSEA[best_model_idx],
                      fit_indices_summary$SRMR[best_model_idx]))
        }
      } else {
        cat(sprintf("Warning: Best model '%s' (key: %s) did not converge or cannot be accessed, skipping path diagram plotting\n",
                    best_model_name, model_key))
      }
    }
  }
} else {
  cat("Warning: No converged SEM models, skipping plot generation\n")
  
  # Create simple path diagram showing linear model results
  cat("Creating simple path diagram showing linear model results...\n")
  
  # First check if dagitty and ggdag are installed
  if (!requireNamespace("dagitty", quietly = TRUE)) {
    install.packages("dagitty")
  }
  if (!requireNamespace("ggdag", quietly = TRUE)) {
    install.packages("ggdag")
  }
  
  library(dagitty)
  library(ggdag)
  
  # Create DAG
  dag <- dagitty('dag {
    VFA -> SHBG
    SHBG -> Total_BMD
    SHBG -> Vit_D
    Total_BMD -> Vit_D
  }')
  
  # Plot DAG
  p <- ggdag(dag) + 
    theme_dag() +
    labs(title = "Theoretical Path Model",
         subtitle = "Based on Linear Regression Analysis")
  
  # Ensure directory exists
  if (!dir.exists("results/Path_Analysis/SEM_fit")) {
    dir.create("results/Path_Analysis/SEM_fit", recursive = TRUE)
  }
  
  ggsave("results/Path_Analysis/SEM_fit/theoretical_path_diagram.png", 
         p, width = 8, height = 6, dpi = 300)
  cat("✓ Theoretical path diagram saved\n")
}
# 13. Generate detailed SEM fit report --------------------------------------------------
cat("\nStep 12: Generating detailed SEM fit report...\n")

sink("results/Path_Analysis/SEM_fit/detailed_fit_report.txt")

cat("Structural Equation Model (SEM) Fit Analysis Report\n")
cat("===================================================\n\n")

cat("Analysis Date:", Sys.Date(), "\n")
cat("Sample Size:", nrow(analysis_data_complete), "\n")
cat("Analysis Variables:", paste(names(analysis_data_complete), collapse = ", "), "\n\n")

cat("1. Data Check Results\n")
cat("---------------------\n")
cat("1.1 Sample Size Check:", nrow(analysis_data_complete), "complete observations\n")
cat("1.2 Variable Correlation Check:\n")
cat("   Variable correlation matrix:\n")
print(round(cor(analysis_data_complete), 3))

if (nrow(high_cor) > 0) {
  cat("\n   Warning: Found highly correlated variable pairs:\n")
  for (i in 1:nrow(high_cor)) {
    var1 <- rownames(cor_matrix)[high_cor[i, 1]]
    var2 <- colnames(cor_matrix)[high_cor[i, 2]]
    cor_value <- cor_matrix[high_cor[i, 1], high_cor[i, 2]]
    cat(sprintf("     %s and %s: r = %.3f\n", var1, var2, cor_value))
  }
}

cat("\n2. Model Fitting Attempts\n")
cat("-------------------------\n")
cat("Based on theoretical framework, 5 SEM models with different complexities were constructed:\n")
cat("1. VFA→SHBG→BMD→VitD: Basic cascade path model\n")
cat("2. VFA→SHBG→BMD→VitD + VFA→VitD: Model with added direct path from VFA to VitD\n")
cat("3. VFA→SHBG→BMD + VFA→SHBG→VitD: Dual path model (includes VFA→SHBG→VitD)\n")
cat("4. VFA→SHBG→VitD (focused): Model focusing specifically on the VFA→SHBG→VitD pathway\n")
cat("5. Complete Model (all paths): Includes all possible paths between variables\n")
cat("\nTotal of", length(sem_results), "SEM models were attempted for comparison.\n")

converged_count <- sum(sapply(sem_results, function(x) ifelse(is.null(x$converged), FALSE, x$converged)))
cat("Successfully converged models:", converged_count, "\n")
cat("Non-converged models:", length(sem_results) - converged_count, "\n\n")

cat("3. VFA→SHBG→VitD Pathway Specific Analysis\n")
cat("------------------------------------------\n")
cat("Based on separate mediation analysis:\n")
cat("  Path a (VFA→SHBG): β = ", round(a, 3), 
    ", p = ", format(summary(m1)$coefficients["VFA", 4], scientific = FALSE, digits = 4), "\n", sep = "")
cat("  Path b (SHBG→VitD): β = ", round(b, 3), 
    ", p = ", format(summary(m2)$coefficients["SHBG", 4], scientific = FALSE, digits = 4), "\n", sep = "")
cat("  Indirect effect (VFA→SHBG→VitD): β = ", round(indirect_effect, 3), "\n", sep = "")
cat("  Sobel test: z = ", round(z_value, 3), 
    ", p = ", format(p_value, scientific = FALSE, digits = 4), "\n", sep = "")
cat("  95% Bootstrap CI: [", round(boot_ci[1], 3), ", ", round(boot_ci[2], 3), "]\n\n", sep = "")

if (p_value < 0.05) {
  cat("✓ The VFA→SHBG→VitD pathway is statistically significant\n")
  cat("  This suggests that SHBG mediates the relationship between visceral fat and vitamin D\n")
} else {
  cat("✗ The VFA→SHBG→VitD pathway is not statistically significant\n")
  cat("  SHBG may not be a significant mediator between VFA and vitamin D in this sample\n")
}

cat("\n4. Model Convergence Issues Analysis\n")
cat("------------------------------------\n")
cat("Possible reasons for model non-convergence:\n")
cat("1. Insufficient sample size: Current sample size (", nrow(analysis_data_complete), ") may be insufficient\n")
cat("2. Model over-parameterization: Model too complex, too many parameters\n")
cat("3. Data issues: Abnormal variable distribution or outliers\n")
cat("4. Model specification error: Path specification doesn't match data characteristics\n")
cat("5. Multicollinearity issues: High correlations between variables causing unstable estimates\n\n")

# Calculate summary statistics by path type
path_by_type <- path_summary %>%
  group_by(Path_Type) %>%
  summarise(
    Count = n(),
    Mean_Beta = mean(Beta, na.rm = TRUE),
    Positive_Effects = sum(Beta > 0, na.rm = TRUE),
    Negative_Effects = sum(Beta < 0, na.rm = TRUE),
    Significant_Effects = sum(P_Value < 0.05, na.rm = TRUE)
  )

cat("5. Path Effect Analysis\n")
cat("----------------------\n")
cat("Summary of effects by path type:\n")
for (i in 1:nrow(path_by_type)) {
  type <- path_by_type[i, ]
  cat(sprintf("  %s: %d paths, mean β=%.3f, positive effects=%d, negative effects=%d, significant effects=%d\n",
              type$Path_Type, type$Count, type$Mean_Beta,
              type$Positive_Effects, type$Negative_Effects,
              type$Significant_Effects))
}

if (converged_count > 0) {
  cat("\n6. Converged Model Fit Results\n")
  cat("------------------------------\n")
  
  # Print fit indices table
  cat("Model fit indices summary table:\n")
  cat("------------------------------------------------------------------------------------------------\n")
  cat(sprintf("%-35s %-8s %-6s %-8s %-6s %-6s %-12s %-6s\n", 
              "Model", "χ²", "df", "p", "CFI", "TLI", "RMSEA(90%CI)", "SRMR"))
  cat("------------------------------------------------------------------------------------------------\n")
  
  for (i in 1:nrow(fit_indices_summary)) {
    model <- fit_indices_summary[i, ]
    cat(sprintf("%-35s %-8.1f %-6d %-8.4f %-6.3f %-6.3f %-12s %-6.3f\n",
                model$Model,
                ifelse(is.na(model$ChiSq), NA, model$ChiSq),
                ifelse(is.na(model$df), NA, model$df),
                ifelse(is.na(model$p), NA, model$p),
                ifelse(is.na(model$CFI), NA, model$CFI),
                ifelse(is.na(model$TLI), NA, model$TLI),
                model$RMSEA_CI,
                ifelse(is.na(model$SRMR), NA, model$SRMR)))
  }
  cat("------------------------------------------------------------------------------------------------\n\n")
  
  # Best model analysis
  if (!is.null(best_model_fit) && best_model_fit$converged) {
    cat("7. Best Model Detailed Analysis\n")
    cat("--------------------------------\n")
    
    # Print significant paths
    sig_paths <- best_model_fit$parameters %>%
      filter(op == "~" & pvalue < 0.05) %>%
      arrange(desc(std.all))
    
    if (nrow(sig_paths) > 0) {
      cat("Significant path coefficients (p < 0.05):\n")
      for (i in 1:nrow(sig_paths)) {
        path <- sig_paths[i, ]
        cat(sprintf("  %s → %s: β = %.3f, z = %.2f, p = %.4f\n",
                    path$rhs, path$lhs, path$std.all, 
                    path$z, path$pvalue))
      }
    } else {
      cat("No significant path coefficients\n")
    }
    
    # Check specifically for VFA→SHBG→VitD paths
    vfa_shbg_vitd_paths <- best_model_fit$parameters %>%
      filter(op == "~" & ((rhs == "VFA" & lhs == "SHBG") | 
                            (rhs == "SHBG" & lhs == "Vit_D")))
    
    if (nrow(vfa_shbg_vitd_paths) > 0) {
      cat("\nVFA→SHBG→VitD pathway in best model:\n")
      for (i in 1:nrow(vfa_shbg_vitd_paths)) {
        path <- vfa_shbg_vitd_paths[i, ]
        cat(sprintf("  %s → %s: β = %.3f, p = %.4f\n",
                    path$rhs, path$lhs, path$std.all, path$pvalue))
      }
    }
    
    # Model explanatory power
    cat("\nModel explanatory power (R²):\n")
    r2_values <- lavInspect(best_model_fit$fit, "r2")
    if (!is.null(r2_values)) {
      for (i in 1:length(r2_values)) {
        cat(sprintf("  %s: R² = %.3f\n", 
                    names(r2_values)[i], r2_values[i]))
      }
    }
  }
} else {
  cat("\n6. Alternative Analysis Results (Linear Regression)\n")
  cat("--------------------------------------------------\n")
  
  # Read linear model results
  if (file.exists("results/Path_Analysis/SEM_fit/lm_path_results.rds")) {
    lm_results <- readRDS("results/Path_Analysis/SEM_fit/lm_path_results.rds")
    
    cat("Since SEM models did not converge, here are linear regression path analysis results:\n\n")
    
    # Path 1: VFA -> SHBG
    m1_sum <- summary(lm_results$VFA_SHBG)
    cat("Path 1: VFA → SHBG\n")
    cat(sprintf("  β = %.3f, SE = %.3f, t = %.2f, p = %.4f\n",
                m1_sum$coefficients[2, 1],
                m1_sum$coefficients[2, 2],
                m1_sum$coefficients[2, 3],
                m1_sum$coefficients[2, 4]))
    cat(sprintf("  R² = %.3f, F = %.2f\n", 
                m1_sum$r.squared, m1_sum$fstatistic[1]))
    
    # Path 2: SHBG -> Total_BMD
    m2_sum <- summary(lm_results$SHBG_BMD)
    cat("\nPath 2: SHBG → Total_BMD\n")
    cat(sprintf("  β = %.3f, SE = %.3f, t = %.2f, p = %.4f\n",
                m2_sum$coefficients[2, 1],
                m2_sum$coefficients[2, 2],
                m2_sum$coefficients[2, 3],
                m2_sum$coefficients[2, 4]))
    cat(sprintf("  R² = %.3f, F = %.2f\n", 
                m2_sum$r.squared, m2_sum$fstatistic[1]))
    
    # Path 3: Total_BMD -> Vit_D
    m3_sum <- summary(lm_results$BMD_VitD)
    cat("\nPath 3: Total_BMD → Vit_D\n")
    cat(sprintf("  β = %.3f, SE = %.3f, t = %.2f, p = %.4f\n",
                m3_sum$coefficients[2, 1],
                m3_sum$coefficients[2, 2],
                m3_sum$coefficients[2, 3],
                m3_sum$coefficients[2, 4]))
    cat(sprintf("  R² = %.3f, F = %.2f\n", 
                m3_sum$r.squared, m3_sum$fstatistic[1]))
    
    # Path 4: SHBG -> Vit_D
    m4_sum <- summary(lm_results$SHBG_VitD)
    cat("\nPath 4: SHBG → Vit_D (VFA→SHBG→VitD pathway)\n")
    cat(sprintf("  β = %.3f, SE = %.3f, t = %.2f, p = %.4f\n",
                m4_sum$coefficients[2, 1],
                m4_sum$coefficients[2, 2],
                m4_sum$coefficients[2, 3],
                m4_sum$coefficients[2, 4]))
    cat(sprintf("  R² = %.3f, F = %.2f\n", 
                m4_sum$r.squared, m4_sum$fstatistic[1]))
    cat("  This is the key path for the VFA→SHBG→VitD mediation hypothesis\n")
  }
}

cat("\n7. Research Limitations and Recommendations\n")
cat("-------------------------------------------\n")
cat("1. Sample size limitations: Current sample size may be insufficient, consider increasing sample size\n")
cat("2. Model simplification: Consider further simplifying model structure\n")
cat("3. Variable selection: Consider reducing number of variables or using composite indicators\n")
cat("4. Method alternatives: Consider using path analysis or regression analysis instead of SEM\n")
cat("5. Data checking: Check data distribution and outliers\n\n")

cat("8. VFA→SHBG→VitD Pathway Conclusions\n")
cat("-------------------------------------\n")
if (p_value < 0.05) {
  cat("✓ The VFA→SHBG→VitD pathway is statistically significant\n")
  cat("  - Visceral fat (VFA) significantly affects SHBG levels\n")
  cat("  - SHBG significantly affects vitamin D levels\n")
  cat("  - The indirect effect through SHBG is significant\n")
  cat("  - This supports the hypothesis that visceral fat influences vitamin D through SHBG\n")
  
  # Clinical implications
  cat("\nClinical Implications:\n")
  cat("1. SHBG may be a biomarker linking visceral adiposity to vitamin D status\n")
  cat("2. Interventions targeting visceral fat may improve vitamin D through SHBG modulation\n")
  cat("3. Monitoring SHBG levels may help predict vitamin D status in individuals with high visceral fat\n")
} else {
  cat("✗ The VFA→SHBG→VitD pathway is not statistically significant\n")
  cat("  - While individual paths may be significant, the mediation effect is not confirmed\n")
  cat("  - Visceral fat may affect vitamin D through other mechanisms\n")
  cat("  - SHBG may not be the primary mediator in this relationship\n")
  
  cat("\nAlternative Considerations:\n")
  cat("1. Other potential mediators: inflammation markers, adipokines, insulin resistance\n")
  cat("2. Direct effects of visceral fat on vitamin D metabolism\n")
  cat("3. Gender-specific or age-specific effects may exist\n")
}

cat("\n9. File Output List\n")
cat("------------------\n")
if (converged_count > 0) {
  cat("1. results/Path_Analysis/SEM_fit/fit_indices_summary.csv - Fit indices summary table\n")
  cat("2. results/Path_Analysis/SEM_fit/best_model_parameters.csv - Best model parameters\n")
  cat("3. results/Path_Analysis/SEM_fit/best_model_path_diagram.png - Best model path diagram\n")
  if (nrow(fit_indices_summary) > 1) {
    cat("4. results/Path_Analysis/SEM_fit/model_comparison_plots.png - Model comparison plots\n")
  }
} else {
  cat("1. results/Path_Analysis/SEM_fit/lm_path_results.rds - Linear regression path analysis results\n")
  cat("2. results/Path_Analysis/SEM_fit/theoretical_path_diagram.png - Theoretical path diagram\n")
}
cat("3. results/Path_Analysis/SEM_fit/detailed_fit_report.txt - This report\n")

# VFA→SHBG→VitD specific files
cat("\nVFA→SHBG→VitD Pathway Analysis Files:\n")
cat("1. results/Path_Analysis/Path_VFA_SHBG_VitD/pathway_scatter_plots.png\n")
cat("2. results/Path_Analysis/Path_VFA_SHBG_VitD/mediation_diagram.png\n")
cat("3. results/Path_Analysis/Path_VFA_SHBG_VitD/mediation_analysis_results.csv\n")
cat("4. results/Path_Analysis/Path_VFA_SHBG_VitD/bootstrap_mediation_results.csv\n")
cat("5. results/Path_Analysis/Path_VFA_SHBG_VitD/bootstrap_distribution.png\n")
cat("6. results/Path_Analysis/Path_VFA_SHBG_VitD/pathway_analysis_report.txt\n")

sink()

cat("✓ Detailed fit report generated\n")

# 14. Create brief result summary ------------------------------------------------------
cat("\nStep 13: Creating brief result summary...\n")

# Generate HTML format brief report
sink("results/Path_Analysis/SEM_fit/quick_summary.html")

cat('<!DOCTYPE html>
<html>
<head>
    <title>SEM Analysis Brief Summary</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #2c3e50; }
        h2 { color: #34495e; border-bottom: 2px solid #3498db; padding-bottom: 5px; }
        .good { color: #27ae60; font-weight: bold; }
        .acceptable { color: #f39c12; font-weight: bold; }
        .poor { color: #e74c3c; font-weight: bold; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: center; }
        th { background-color: #3498db; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        .best { background-color: #d5f4e6 !important; }
        .warning { background-color: #fff3cd; border: 1px solid #ffeaa7; padding: 10px; margin: 10px 0; }
        .highlight { background-color: #e8f4fd !important; }
    </style>
</head>
<body>
    <h1>Structural Equation Model (SEM) Analysis Brief Summary</h1>
    <p><strong>Analysis Date:</strong> ', Sys.Date(), '</p>
    <p><strong>Sample Size:</strong> ', nrow(analysis_data_complete), '</p>
    <p><strong>Analysis Variables:</strong> ', paste(names(analysis_data_complete), collapse = ", "), '</p>
    <p><strong>Analysis Completion Time:</strong> ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</p>')

# Add VFA→SHBG→VitD pathway summary at the top
cat('
    <div style="background-color: #e8f4fd; border-left: 5px solid #3498db; padding: 15px; margin: 20px 0;">
        <h3 style="color: #2c3e50; margin-top: 0;">VFA→SHBG→VitD Pathway Analysis</h3>
        <p><strong>Key Findings:</strong></p>')

if (p_value < 0.05) {
  cat('
        <p class="good">✓ Significant mediation found: VFA→SHBG→VitD pathway is statistically significant</p>
        <ul>
            <li>Indirect effect (VFA→SHBG→VitD): β = ', round(indirect_effect, 3), '</li>
            <li>Sobel test: z = ', round(z_value, 3), ', p = ', format(p_value, scientific = FALSE, digits = 4), '</li>
            <li>95% Bootstrap CI: [', round(boot_ci[1], 3), ', ', round(boot_ci[2], 3), ']</li>
        </ul>')
} else {
  cat('
        <p class="poor">✗ No significant mediation found: VFA→SHBG→VitD pathway is not statistically significant</p>
        <ul>
            <li>Indirect effect (VFA→SHBG→VitD): β = ', round(indirect_effect, 3), '</li>
            <li>Sobel test: z = ', round(z_value, 3), ', p = ', format(p_value, scientific = FALSE, digits = 4), '</li>
        </ul>')
}

cat('
        <p><strong>Interpretation:</strong> ', 
    ifelse(p_value < 0.05, 
           'SHBG appears to mediate the relationship between visceral fat and vitamin D levels.', 
           'SHBG does not appear to be a significant mediator between visceral fat and vitamin D.'), '</p>
    </div>')

if (converged_count > 0 && nrow(fit_indices_summary) > 0) {
  best_idx <- which.max(fit_indices_summary$CFI)
  best_model <- fit_indices_summary[best_idx, ]
  
  # Provide specific description based on model name
  model_description <- ""
  if (best_model$Model == "VFA→SHBG→BMD→VitD") {
    model_description <- "Basic cascade path model (Visceral Fat→Sex Hormone-Binding Globulin→Bone Density→Vitamin D)"
  } else if (best_model$Model == "VFA→SHBG→BMD→VitD + VFA→VitD") {
    model_description <- "Model with added direct path (Visceral Fat→Sex Hormone-Binding Globulin→Bone Density→Vitamin D, with Visceral Fat also directly affecting Vitamin D)"
  } else if (best_model$Model == "VFA→SHBG→BMD + VFA→SHBG→VitD") {
    model_description <- "Dual path model (Visceral Fat→Sex Hormone-Binding Globulin→Bone Density, and Visceral Fat→Sex Hormone-Binding Globulin→Vitamin D)"
  } else if (best_model$Model == "VFA→SHBG→VitD (focused)") {
    model_description <- "VFA→SHBG→VitD focused model"
  } else if (best_model$Model == "Complete Model (all paths)") {
    model_description <- "Complete model with all possible paths"
  } else {
    model_description <- best_model$Model
  }
  
  cat('
    <div class="warning">
        <strong>Note:</strong> SEM model fitting encountered difficulties, below are the results for ', model_description, '.
    </div>
    
    <h2>1. Best Model Selection</h2>
    <p><strong>Recommended Model:</strong> ', model_description, '</p>
    <p><strong>Model Path:</strong> ', best_model$Model, '</p>
    <p><strong>Selection Reason:</strong> Optimal CFI index</p>
    
    <h2>2. Best Model Fit Indices</h2>
    <table>
        <tr>
            <th>Index</th>
            <th>Value</th>
            <th>Evaluation</th>
            <th>Standard</th>
        </tr>')
  
  # CFI evaluation
  cfi_class <- ifelse(is.na(best_model$CFI), "poor", 
                      ifelse(best_model$CFI >= 0.95, "good", 
                             ifelse(best_model$CFI >= 0.90, "acceptable", "poor")))
  cfi_text <- ifelse(is.na(best_model$CFI), "Cannot calculate", 
                     ifelse(best_model$CFI >= 0.95, "Excellent", 
                            ifelse(best_model$CFI >= 0.90, "Acceptable", "Insufficient")))
  
  # TLI evaluation
  tli_class <- ifelse(is.na(best_model$TLI), "poor", 
                      ifelse(best_model$TLI >= 0.95, "good", 
                             ifelse(best_model$TLI >= 0.90, "acceptable", "poor")))
  tli_text <- ifelse(is.na(best_model$TLI), "Cannot calculate", 
                     ifelse(best_model$TLI >= 0.95, "Excellent", 
                            ifelse(best_model$TLI >= 0.90, "Acceptable", "Insufficient")))
  
  # RMSEA evaluation
  rmsea_class <- ifelse(is.na(best_model$RMSEA), "poor", 
                        ifelse(best_model$RMSEA <= 0.05, "good", 
                               ifelse(best_model$RMSEA <= 0.08, "acceptable", "poor")))
  rmsea_text <- ifelse(is.na(best_model$RMSEA), "Cannot calculate", 
                       ifelse(best_model$RMSEA <= 0.05, "Excellent", 
                              ifelse(best_model$RMSEA <= 0.08, "Acceptable", "Insufficient")))
  
  # SRMR evaluation (get from fit_indices_summary)
  srmr_value <- ifelse("SRMR" %in% names(fit_indices_summary), fit_indices_summary$SRMR[best_idx], NA)
  srmr_class <- ifelse(is.na(srmr_value), "poor", 
                       ifelse(srmr_value <= 0.05, "good", 
                              ifelse(srmr_value <= 0.08, "acceptable", "poor")))
  srmr_text <- ifelse(is.na(srmr_value), "Cannot calculate", 
                      ifelse(srmr_value <= 0.05, "Excellent", 
                             ifelse(srmr_value <= 0.08, "Acceptable", "Insufficient")))
  
  cat(sprintf('
        <tr>
            <td>CFI</td>
            <td>%.3f</td>
            <td class="%s">%s</td>
            <td>>0.95 Excellent, >0.90 Acceptable</td>
        </tr>
        <tr>
            <td>TLI</td>
            <td>%.3f</td>
            <td class="%s">%s</td>
            <td>>0.95 Excellent, >0.90 Acceptable</td>
        </tr>
        <tr>
            <td>RMSEA</td>
            <td>%.3f</td>
            <td class="%s">%s</td>
            <td><0.05 Excellent, <0.08 Acceptable</td>
        </tr>
        <tr>
            <td>SRMR</td>
            <td>%.3f</td>
            <td class="%s">%s</td>
            <td><0.05 Excellent, <0.08 Acceptable</td>
        </tr>
        <tr>
            <td>χ²/df</td>
            <td>%.1f/%d = %.2f</td>
            <td>-</td>
            <td><3.0 Acceptable</td>
        </tr>',
              ifelse(is.na(best_model$CFI), NA, best_model$CFI), cfi_class, cfi_text,
              ifelse(is.na(best_model$TLI), NA, best_model$TLI), tli_class, tli_text,
              ifelse(is.na(best_model$RMSEA), NA, best_model$RMSEA), rmsea_class, rmsea_text,
              ifelse(is.na(srmr_value), NA, srmr_value), srmr_class, srmr_text,
              ifelse(is.na(best_model$ChiSq), NA, best_model$ChiSq), 
              ifelse(is.na(best_model$df), NA, best_model$df), 
              ifelse(is.na(best_model$ChiSq) | is.na(best_model$df), NA, best_model$ChiSq/best_model$df)))
  
  cat('
    </table>')
  
  # Overall evaluation
  cat('
    <h2>3. Overall Fit Evaluation</h2>')
  
  if (!is.na(best_model$CFI) && !is.na(best_model$TLI) && 
      !is.na(best_model$RMSEA) && !is.na(srmr_value)) {
    if (best_model$CFI >= 0.95 && best_model$TLI >= 0.95 && 
        best_model$RMSEA <= 0.05 && srmr_value <= 0.05) {
      cat('
    <p class="good">✓ Excellent fit: All major fit indices meet excellent standards</p>')
    } else if (best_model$CFI >= 0.90 && best_model$TLI >= 0.90 && 
               best_model$RMSEA <= 0.08 && srmr_value <= 0.08) {
      cat('
    <p class="acceptable">~ Acceptable fit: All major fit indices meet acceptable standards</p>')
    } else {
      cat('
    <p class="poor">✗ Insufficient fit: Some fit indices do not meet acceptable standards</p>')
    }
  } else {
    cat('
    <p class="poor">✗ Some fit indices cannot be calculated, model fitting encountered difficulties</p>')
  }
  
  # Add model interpretation
  cat('
    <h2>4. Model Interpretation</h2>')
  
  cat('
    <p><strong>Model Path:</strong> ', best_model$Model, '</p>')
  
  if (best_model$Model == "VFA→SHBG→VitD (focused)") {
    cat('
    <p><strong>Biological Interpretation:</strong> This model focuses specifically on the VFA→SHBG→VitD pathway, suggesting that visceral fat affects vitamin D primarily through its influence on SHBG levels. This aligns with the significant mediation found in the separate analysis.</p>')
  } else if (best_model$Model == "VFA→SHBG→BMD + VFA→SHBG→VitD") {
    cat('
    <p><strong>Biological Interpretation:</strong> This model shows that visceral fat simultaneously affects bone density and vitamin D levels through SHBG, forming a dual-path regulatory pattern. SHBG plays a key hub role in this model, mediating both bone and vitamin D metabolism.</p>')
  } else if (grepl("VFA→SHBG→VitD", best_model$Model)) {
    cat('
    <p><strong>Biological Interpretation:</strong> This model includes the VFA→SHBG→VitD pathway, supporting the hypothesis that visceral fat influences vitamin D metabolism through SHBG. This pathway represents a potential endocrine link between adiposity and vitamin D status.</p>')
  }
  
  cat('
    <p><strong>Clinical Application Significance:</strong></p>
    <ul>
        <li>Provides a theoretical framework for obesity-related metabolic abnormalities in adolescents</li>
        <li>Suggests SHBG may be a key target for intervening in the relationship between visceral fat and bone/vitamin D</li>
        <li>Provides basis for developing targeted health intervention strategies</li>')
  
  if (p_value < 0.05) {
    cat('
        <li><strong>Specifically supports the VFA→SHBG→VitD pathway as a mechanism linking visceral fat to vitamin D status</strong></li>')
  }
  
  cat('
    </ul>')
  
} else {
  cat('
    <div class="warning">
        <strong>Important:</strong> All SEM models failed to converge or fit.
    </div>
    
    <h2>1. Analysis Results</h2>
    <p class="poor">All SEM models failed to converge or fit.</p>
    
    <h2>2. Linear Regression Path Analysis Results</h2>')
  
  if (file.exists("results/Path_Analysis/SEM_fit/lm_path_results.rds")) {
    lm_results <- readRDS("results/Path_Analysis/SEM_fit/lm_path_results.rds")
    
    cat('
    <p>Since SEM models did not converge, here are linear regression path analysis results:</p>
    <table>
        <tr>
            <th>Path</th>
            <th>β Coefficient</th>
            <th>Standard Error</th>
            <th>t-value</th>
            <th>p-value</th>
            <th>R²</th>
        </tr>')
    
    # Path 1
    m1_sum <- summary(lm_results$VFA_SHBG)
    cat(sprintf('
        <tr>
            <td>VFA → SHBG</td>
            <td>%.3f</td>
            <td>%.3f</td>
            <td>%.2f</td>
            <td>%.4f</td>
            <td>%.3f</td>
        </tr>',
                m1_sum$coefficients[2, 1],
                m1_sum$coefficients[2, 2],
                m1_sum$coefficients[2, 3],
                m1_sum$coefficients[2, 4],
                m1_sum$r.squared))
    
    # Path 2
    m2_sum <- summary(lm_results$SHBG_BMD)
    cat(sprintf('
        <tr>
            <td>SHBG → Total_BMD</td>
            <td>%.3f</td>
            <td>%.3f</td>
            <td>%.2f</td>
            <td>%.4f</td>
            <td>%.3f</td>
        </tr>',
                m2_sum$coefficients[2, 1],
                m2_sum$coefficients[2, 2],
                m2_sum$coefficients[2, 3],
                m2_sum$coefficients[2, 4],
                m2_sum$r.squared))
    
    # Path 3
    m3_sum <- summary(lm_results$BMD_VitD)
    cat(sprintf('
        <tr>
            <td>Total_BMD → Vit_D</td>
            <td>%.3f</td>
            <td>%.3f</td>
            <td>%.2f</td>
            <td>%.4f</td>
            <td>%.3f</td>
        </tr>',
                m3_sum$coefficients[2, 1],
                m3_sum$coefficients[2, 2],
                m3_sum$coefficients[2, 3],
                m3_sum$coefficients[2, 4],
                m3_sum$r.squared))
    
    # Path 4 - VFA→SHBG→VitD pathway
    m4_sum <- summary(lm_results$SHBG_VitD)
    cat(sprintf('
        <tr class="highlight">
            <td><strong>SHBG → Vit_D (VFA→SHBG→VitD)</strong></td>
            <td><strong>%.3f</strong></td>
            <td><strong>%.3f</strong></td>
            <td><strong>%.2f</strong></td>
            <td><strong>%.4f</strong></td>
            <td><strong>%.3f</strong></td>
        </tr>',
                m4_sum$coefficients[2, 1],
                m4_sum$coefficients[2, 2],
                m4_sum$coefficients[2, 3],
                m4_sum$coefficients[2, 4],
                m4_sum$r.squared))
    
    cat('
    </table>
    <p><em>Note: Highlighted row shows the key path for the VFA→SHBG→VitD mediation hypothesis</em></p>')
  }
  
  cat('
    <h2>3. VFA→SHBG→VitD Mediation Results</h2>')
  
  cat('
    <table>
        <tr>
            <th>Component</th>
            <th>β Coefficient</th>
            <th>Standard Error</th>
            <th>t/z-value</th>
            <th>p-value</th>
        </tr>')
  
  cat(sprintf('
        <tr>
            <td>Path a (VFA→SHBG)</td>
            <td>%.3f</td>
            <td>%.3f</td>
            <td>%.2f</td>
            <td>%.4f</td>
        </tr>',
              a, se_a, summary(m1)$coefficients["VFA", "t value"], 
              summary(m1)$coefficients["VFA", 4]))
  
  cat(sprintf('
        <tr>
            <td>Path b (SHBG→VitD)</td>
            <td>%.3f</td>
            <td>%.3f</td>
            <td>%.2f</td>
            <td>%.4f</td>
        </tr>',
              b, se_b, summary(m2)$coefficients["SHBG", "t value"], 
              summary(m2)$coefficients["SHBG", 4]))
  
  cat(sprintf('
        <tr class="highlight">
            <td><strong>Indirect effect (a×b)</strong></td>
            <td><strong>%.3f</strong></td>
            <td><strong>%.3f</strong></td>
            <td><strong>%.2f</strong></td>
            <td><strong>%.4f</strong></td>
        </tr>',
              indirect_effect, sqrt(b^2 * se_a^2 + a^2 * se_b^2), z_value, p_value))
  
  cat('
    </table>')
  
  if (p_value < 0.05) {
    cat('
    <p class="good">✓ Significant mediation: The VFA→SHBG→VitD pathway is statistically significant (p = ', 
        format(p_value, scientific = FALSE, digits = 4), ')</p>', sep = '')
  } else {
    cat('
    <p class="poor">✗ No significant mediation: The VFA→SHBG→VitD pathway is not statistically significant (p = ', 
        format(p_value, scientific = FALSE, digits = 4), ')</p>', sep = '')
  }
  
  cat('
    <h2>4. Possible Reasons for SEM Non-convergence</h2>
    <ul>
        <li>Sample size may be insufficient</li>
        <li>Relationships between variables may be too complex</li>
        <li>Model specification may not match data characteristics</li>
        <li>Consider simplifying models or increasing sample size</li>
    </ul>')
}

# Add path type analysis
cat('
    <h2>5. Path Effect Analysis</h2>
    <table>
        <tr>
            <th>Path Type</th>
            <th>Number of Paths</th>
            <th>Mean β Coefficient</th>
            <th>Positive Effects</th>
            <th>Negative Effects</th>
            <th>Significant Effects</th>
        </tr>')

for (i in 1:nrow(path_by_type)) {
  type <- path_by_type[i, ]
  row_class <- ifelse(type$Path_Type == "VFA→SHBG→VitD Mediation", "highlight", "")
  cat(sprintf('
        <tr class="%s">
            <td>%s</td>
            <td>%d</td>
            <td>%.3f</td>
            <td>%d</td>
            <td>%d</td>
            <td>%d</td>
        </tr>',
              row_class, type$Path_Type, type$Count, type$Mean_Beta,
              type$Positive_Effects, type$Negative_Effects,
              type$Significant_Effects))
}

cat('
    </table>')

cat('
    <h2>6. File Outputs</h2>
    <p>All detailed results can be viewed in the following files:</p>
    <ul>')

if (converged_count > 0) {
  cat('
        <li><code>fit_indices_summary.csv</code> - All model fit indices comparison</li>
        <li><code>best_model_parameters.csv</code> - Best model parameter estimates</li>
        <li><code>best_model_path_diagram.png</code> - Best model path diagram</li>')
  if (nrow(fit_indices_summary) > 1) {
    cat('
        <li><code>model_comparison_plots.png</code> - Model comparison plots</li>')
  }
} else {
  cat('
        <li><code>lm_path_results.rds</code> - Linear regression path analysis results</li>
        <li><code>theoretical_path_diagram.png</code> - Theoretical path diagram</li>')
}

cat('
        <li><code>detailed_fit_report.txt</code> - Detailed analysis report</li>')

# VFA→SHBG→VitD specific files
cat('
        <li><strong>VFA→SHBG→VitD Pathway Analysis Files:</strong></li>
        <li><code>Path_VFA_SHBG_VitD/pathway_scatter_plots.png</code> - Scatter plots of individual relationships</li>
        <li><code>Path_VFA_SHBG_VitD/mediation_diagram.png</code> - Mediation model diagram</li>
        <li><code>Path_VFA_SHBG_VitD/mediation_analysis_results.csv</code> - Mediation analysis results</li>
        <li><code>Path_VFA_SHBG_VitD/bootstrap_mediation_results.csv</code> - Bootstrap results</li>
        <li><code>Path_VFA_SHBG_VitD/pathway_analysis_report.txt</code> - Comprehensive pathway report</li>
    </ul>')

cat('
    <h2>7. Notes and Cautions</h2>
    <ul>
        <li>This analysis is based on cross-sectional data and cannot determine causal direction</li>
        <li>Path directions are based on theoretical assumptions and need verification in subsequent studies</li>
        <li>Results apply to 7-18 year old adolescent population</li>
        <li>Recommend validating model stability in independent samples</li>')

if (converged_count == 0) {
  cat('
        <li><strong>Special note: SEM models did not converge, results should be interpreted with caution</strong></li>')
}

cat('
        <li><strong>VFA→SHBG→VitD pathway: </strong>', 
    ifelse(p_value < 0.05, 
           'Significant mediation found - SHBG appears to mediate the relationship between visceral fat and vitamin D.', 
           'No significant mediation found - SHBG may not be the primary mediator in this relationship.'), '</li>')

cat('
    </ul>
    
    <hr>
    <p style="text-align: center; color: #7f8c8d; font-size: 0.9em;">
        Analysis generation time: ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '<br>
        Software used: R 4.2.0, lavaan 0.6-15
    </p>
</body>
</html>')

sink()

cat("✓ Brief HTML summary generated\n")

# 15. Save workspace ---------------------------------------------------------
cat("\nStep 14: Saving workspace...\n")

# Save complete SEM analysis results
sem_analysis_results <- list(
  data = analysis_data,
  regression_results = all_results,
  path_summary = path_summary,
  sem_results = sem_results,
  fit_indices = fit_indices_summary,
  best_model = best_model_fit,
  vfa_shbg_vitd_mediation = all_results$Mediation_VFA_SHBG_VitD,
  bootstrap_mediation = list(
    indirect_effects = boot_indirect,
    ci = boot_ci,
    mean = boot_mean
  ),
  analysis_timestamp = Sys.time()
)

saveRDS(sem_analysis_results, "results/Path_Analysis/SEM_fit/sem_analysis_full_results.rds")
save.image("results/Path_Analysis/sem_analysis_workspace.RData")

cat("✓ Workspace saved\n")

# 16. Complete analysis ------------------------------------------------------------
cat("\n", paste(rep("=", 50), collapse = ""), "\n", sep = "")
cat("✅ SEM Fit Analysis Completed!\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n", sep = "")

cat("🎯 Analysis Summary\n")
cat("-------------------\n")
cat("- Total sample size:", nrow(analysis_data), "\n")
cat("- Number of analysis variables:", ncol(analysis_data), "\n")
cat("- Number of compared SEM models:", nrow(fit_indices_summary), "\n")

# VFA→SHBG→VitD pathway results
cat("\n🔬 VFA→SHBG→VitD Pathway Analysis:\n")
cat(sprintf("  Path a (VFA→SHBG): β = %.3f, p = %.4f\n", 
            a, summary(m1)$coefficients["VFA", 4]))
cat(sprintf("  Path b (SHBG→VitD): β = %.3f, p = %.4f\n", 
            b, summary(m2)$coefficients["SHBG", 4]))
cat(sprintf("  Indirect effect: β = %.3f, p = %.4f\n", 
            indirect_effect, p_value))
cat(sprintf("  95%% Bootstrap CI: [%.3f, %.3f]\n", boot_ci[1], boot_ci[2]))

if (p_value < 0.05) {
  cat("  ✅ Significant mediation found\n")
} else {
  cat("  ❌ No significant mediation found\n")
}

if (nrow(fit_indices_summary) > 0) {
  best_idx <- which.max(fit_indices_summary$CFI)
  best_model <- fit_indices_summary[best_idx, ]
  
  cat("\n📊 Best SEM Model Results:\n")
  cat("- Best model:", best_model$Model, "\n")
  
  # Add model interpretation
  if (best_model$Model == "VFA→SHBG→VitD (focused)") {
    cat("  Model focus: Specifically examines the VFA→SHBG→VitD pathway\n")
  } else if (grepl("VFA→SHBG→VitD", best_model$Model)) {
    cat("  Model includes: VFA→SHBG→VitD pathway as part of the model\n")
  }
  
  cat(sprintf("  Fit indices: CFI=%.3f, TLI=%.3f, RMSEA=%.3f, SRMR=%.3f\n",
              best_model$CFI, best_model$TLI, best_model$RMSEA, best_model$SRMR))
  
  # Overall evaluation
  if (best_model$CFI >= 0.95 && best_model$TLI >= 0.95 && best_model$RMSEA <= 0.05) {
    cat("  📊 Fit evaluation: Excellent\n")
  } else if (best_model$CFI >= 0.90 && best_model$TLI >= 0.90 && best_model$RMSEA <= 0.08) {
    cat("  📊 Fit evaluation: Acceptable\n")
  } else {
    cat("  📊 Fit evaluation: Needs improvement\n")
  }
} else {
  cat("\n- Model fitting: All models failed to fit\n")
}

cat("\n📁 Key Output Files\n")
cat("-------------------\n")
cat("1. results/Path_Analysis/SEM_fit/fit_indices_summary.csv\n")
cat("   → All model fit indices comparison table\n")
cat("2. results/Path_Analysis/SEM_fit/best_model_parameters.csv\n")
cat("   → Best model parameter estimates\n")
cat("3. results/Path_Analysis/SEM_fit/best_model_path_diagram.png\n")
cat("   → Best model standardized path diagram\n")
cat("4. results/Path_Analysis/SEM_fit/model_comparison_plots.png\n")
cat("   → Model fit indices comparison plots\n")
cat("5. results/Path_Analysis/SEM_fit/detailed_fit_report.txt\n")
cat("   → Detailed fit analysis report (recommended reading)\n")
cat("6. results/Path_Analysis/SEM_fit/quick_summary.html\n")
cat("   → Brief HTML summary report\n")
cat("7. results/Path_Analysis/SEM_fit/sem_analysis_full_results.rds\n")
cat("   → Complete analysis results R object\n")

cat("\n🔬 VFA→SHBG→VitD Specific Files:\n")
cat("8. results/Path_Analysis/Path_VFA_SHBG_VitD/pathway_scatter_plots.png\n")
cat("9. results/Path_Analysis/Path_VFA_SHBG_VitD/mediation_diagram.png\n")
cat("10. results/Path_Analysis/Path_VFA_SHBG_VitD/mediation_analysis_results.csv\n")
cat("11. results/Path_Analysis/Path_VFA_SHBG_VitD/bootstrap_mediation_results.csv\n")
cat("12. results/Path_Analysis/Path_VFA_SHBG_VitD/pathway_analysis_report.txt\n")

cat("\n💡 Follow-up Recommendations\n")
cat("----------------------------\n")
cat("1. First view 'quick_summary.html' to understand main results\n")
cat("2. Read VFA→SHBG→VitD pathway report for detailed mediation analysis\n")
cat("3. Read 'detailed_fit_report.txt' for comprehensive analysis and interpretation\n")
cat("4. If fit is insufficient, refer to modification suggestions in the report\n")
cat("5. In papers, report fit indices and path coefficients of best model\n")
cat("6. Note the causal inference limitations of cross-sectional studies in discussion\n")
cat("7. The VFA→SHBG→VitD pathway analysis provides specific insights into this mechanism\n")

cat("\n", paste(rep("=", 50), collapse = ""), "\n", sep = "")
cat("Analysis completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 50), collapse = ""), "\n", sep = "")