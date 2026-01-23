library(forestploter)
library(grid)
library(readxl)
library(cowplot)  # for plot_grid and save_plot

# 自定义森林图主题
tm <- forest_theme(
  base_size = 12,
  ci_lwd = 1.5,
  
  # 参考线样式（OR = 1 的竖线）
  refline_gp = gpar(
    lwd = 1.5,
    lty = "dashed",
    col = "grey20"
  ),
  
  # 汇总效应（summary）的填充色和边框色
  summary_fill = "#0000FF",
  summary_col = "#0000FF",
  
  # 脚注样式
  footnote_gp = gpar(
    cex = 0.8,
    fontface = "italic",
    col = "grey30"
  ),
  
  ci_pch = 15,
  legend_name = "Group",
  legend_value = c("Malnutrition", "Obesity"),
  core = list(bg_params = list(fill = c("#FFFFFF", "#f5f7f6"), col = NA))
)

# 设置工作目录（请确保路径中不含中文或特殊字符，或使用英文路径）
setwd("C:/Users/loveweiyi/Desktop/data/Nutrition/table2")

# 要处理的文件编号
file_numbers <- 8:16
plots_list <- list()

# 变量名映射（原始列名 → 显示名称）
variable_mapping <- c(
  "Genger"     = "Malnutrition",
  "Vidms"      = "Vit D",
  "VD3MS"      = "Vit D3",
  "TST"        = "Testo",
  "EST"        = "E2",
  "SHBG"       = "SHBG",
  "HSCRP"      = "hs-CRP",
  "FATA"       = "VFA",
  "FATM"       = "VFM",
  "FATV"       = "VFV",
  "DXXHEBMD"   = "BMD (Head)",
  "DXXLLBMD"   = "BMD (Lower Limb)",
  "DXXLSBMD"   = "BMD (L Spine)",
  "DXXPEBMD"   = "BMD (Pelvis)",
  "DXDTOBMC"   = "Total BMC",
  "DXDTOBMD"   = "Total BMD"
)

# 分类映射（原始类别 → 显示名称）
category_mapping <- c(
  "营养不良" = "Malnutrition",
  "肥胖"     = "Obesity"
)

# 循环处理每个文件
for (i in file_numbers) {
  
  # 读取第一个 Excel 文件（例如：table_2_8.xlsx）
  filename1 <- paste0("table_2_", i, ".xlsx")
  data1 <- read_excel(filename1)
  
  # 应用变量和分类映射
  data1$Variable <- ifelse(data1$Variable %in% names(variable_mapping), 
                           variable_mapping[data1$Variable], 
                           data1$Variable)
  data1$Category <- ifelse(data1$Category %in% names(category_mapping), 
                           category_mapping[data1$Category], 
                           data1$Category)
  
  # 提取有效行（第3到20行）
  df1 <- data1[3:20, ]
  df1$`OR(95%CI)` <- paste0(df1$Estimate, " (", df1$X2.5.., ", ", df1$X97.5.., ")")
  
  # 处理缺失值
  df1[, 1:6][is.na(df1[, 1:6])] <- " "
  
  # 重命名列
  colnames(df1)[2] <- "Group"
  colnames(df1)[6] <- "P"
  
  # 添加空白列用于对齐（可选）
  df1$` ` <- paste(rep("   ", nrow(df1)), collapse = "  ")
  
  # 设置颜色：Malnutrition 用红色，Obesity 用青色
  colors1 <- ifelse(df1$Group == "Malnutrition", "#FF6B6B", "#4ECDC4")
  
  # 绘制第一个森林图（OR 尺度）
  p1 <- forest(
    data = df1[, c(1, 2, 7, 8, 6)],
    lower = df1$X2.5..,
    upper = df1$X97.5..,
    est = df1$Estimate,
    ci_column = 4,
    sizes = 1,
    ref_line = 1,
    ticks_at = c(0, 1, 2),
    xlim = c(0, 2),
    color = colors1
  )
  
  # 读取第二个 Excel 文件（例如：肥胖table_2_8.xlsx）
  filename2 <- paste0("肥胖table_2_", i, ".xlsx")
  data2 <- read_excel(filename2)
  
  # 同样应用映射
  data2$Variable <- ifelse(data2$Variable %in% names(variable_mapping), 
                           variable_mapping[data2$Variable], 
                           data2$Variable)
  data2$Category <- ifelse(data2$Category %in% names(category_mapping), 
                           category_mapping[data2$Category], 
                           data2$Category)
  
  # 提取有效行（第1到12行）
  df2 <- data2[1:12, ]
  df2$`log(OR) (95%CI)` <- paste0(df2$Estimate, " (", df2$X2.5.., ", ", df2$X97.5.., ")")
  
  df2[, 1:6][is.na(df2[, 1:6])] <- " "
  colnames(df2)[2] <- "Group"
  colnames(df2)[6] <- "P"
  df2$` ` <- paste(rep("   ", nrow(df2)), collapse = "  ")
  
  colors2 <- ifelse(df2$Group == "Malnutrition", "#FF6B6B", "#4ECDC4")
  
  # 绘制第二个森林图（log(OR) 尺度）
  p2 <- forest(
    data = df2[, c(1, 2, 7, 8, 6)],
    lower = df2$X2.5..,
    upper = df2$X97.5..,
    est = df2$Estimate,
    ci_column = 4,
    sizes = 1,
    ref_line = 1,
    ticks_at = c(-30, 0, 50),
    xlim = c(-30, 50),
    color = colors2
  )
  
  # 垂直拼接两个图
  combined_plot <- plot_grid(p1, p2, 
                             ncol = 1, 
                             align = 'v',
                             rel_heights = c(0.9, 0.5))
  
  # 保存图像
  filename3 <- paste0("combined_plot_", i, ".png")
  save_plot(
    filename = filename3,
    plot = combined_plot,
    dpi = 300,
    base_width = 10,
    base_height = 9.5
  )
}
