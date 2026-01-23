# NHANES Nutritional Status and Health Outcomes Analysis Pipeline

This repository contains the complete R analysis pipeline for a study investigating the associations between nutritional status (malnutrition and obesity) and various health outcomes in children and adolescents using data from the National Health and Nutrition Examination Survey (NHANES).

## Project Structure

The analysis is divided into eight sequential steps, each implemented in a dedicated R script:

*   `01_data_download_nhanesA.R`: Downloads and merges raw NHANES data from multiple survey cycles (1999-2018).
*   `02_participant_selection.R`: Filters participants based on age and selects relevant variables.
*   `03_definition_nutritional_groups.R`: Defines nutritional groups (e.g., malnutrition, normal, obesity, metabolic syndrome) using established clinical criteria and BMI z-scores.
*   `04_descriptive_statistics.R`: Generates descriptive statistics (Table 1) and performs multinomial logistic regression for subgroup analysis (Table 2).
*   `05_restricted_cubic_splines.R`: Explores potential non-linear relationships between Vitamin D and key health indicators using Restricted Cubic Splines (RCS).
*   `06_subgroup_analysis_forest_plot.R`: Creates forest plots to visualize the results of the subgroup analyses from Step 4.
*   `07_network_graph_analysis.R`: Constructs and analyzes a partial correlation network to understand the interrelationships among key continuous variables.
*   `08_SME.R`: Performs Structural Equation Modeling (SEM) to test hypothesized causal pathways linking nutritional status, vitamin D levels, adiposity, and health outcomes.

## Prerequisites

Before running the scripts, ensure you have the following R packages installed: `nhanesA`, `dplyr`, `readxl`, `openxlsx`, `nnet`, `car`, `rms`, `Hmisc`, `ggplot2`, `forestploter`, `gridExtra`, `bootnet`, `qgraph`, `igraph`, `ppcor`, `zscorer`, `anthro`, `lavaan`, `semPlot`, `tidyverse`, `knitr`, `kableExtra`, and other dependencies as listed in each script.

You can install them with:
```R
# Example for a few key packages; refer to each script for a full list
install.packages(c("nhanesA", "dplyr", "rms", "qgraph", "forestploter", "lavaan", "semPlot"))

How to Run£º
1.Set your working directory: Open the project folder in RStudio and set the working directory to this location (Session > Set Working Directory > To Source File Location).
2.Run the scripts in order: Execute the R scripts sequentially from 01 to 07. Each script depends on the output of the previous one(s).
3.Important Note: Script 02 will generate the core dataset data.csv in the working directory, which is used by all subsequent scripts.
4.Output: The final outputs (tables, figures, and network analysis results) will be saved in the working directory or its subdirectories (e.g., the results folder created by 07_network_graph_analysis.R).

Data Source:
This study utilizes publicly available data from the NHANES database (https://www.cdc.gov/nchs/nhanes/).

Author
Wei Yi

```
