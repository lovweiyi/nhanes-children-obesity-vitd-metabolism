library(Hmisc)
library(corrplot)
library(ggplot2)
library(reshape2)


data_orignal <- read_csv("C:/Users/loveweiyi/Desktop/data.csv")

data <- data_orignal
data <- data %>% mutate(across(c(1,2,4,5,40), as.factor))
data <- as.data.frame(data)
str(data)
# data <- data[,-c(12:15)]

data <- data[data$RIDAGEYR > 7 & data$RIDAGEYR < 17, ]
cols_to_remove <- c("IN", "INSI", "FERSI", "RBF", "RBFSI", "FOLSI","FOL","LBXIN","LBDINSI","LBDRBF","LBXRBFSI", "LBDFOL","LBXFOLSI","LBDFERSI")
data <- data[, !names(data) %in% cols_to_remove]
colnames(data)
colnames(data)[12:32] <- c(
  "TSH", "HDL", "TC", "GLU", "Vit D", "Vit D2", "Vit D3", 
  "VFA", "VFM", "VFV", "BMC (Head)", "BMD (Head)", 
  "BMD (Lower Limb)", "BMD (L Spine)", "BMD (Pelvis)", 
  "Total BMC", "Total BMD", "Testo", "E2", "SHBG", "hs-CRP"
)
#----------------------------数据集生成------------------------------#
data$group2 <- factor(data$Group, levels = c("消瘦", "正常", "代谢综合征", "肥胖"))



#----------------------------差异性柱状图------------------------------#


#----------------------------相关性分析------------------------------#
df <- data[,c(12:33)]

# 计算相关性矩阵和P值矩阵
cor_results <- rcorr(as.matrix(df[, sapply(df, is.numeric)]), type = "pearson") # 默认皮尔逊相关系数
cor_matrix <- cor_results$r  # 提取相关性矩阵
p_values <- cor_results$P    # 提取P值矩阵

# 定义一个函数用于根据P值生成显著性标记
significance_markers <- function(p) {
  ifelse(p < 0.001, "***",  # P < 0.001
         ifelse(p < 0.01, "**",  # 0.001 <= P < 0.01
                ifelse(p < 0.05, "*", "")))  # 0.01 <= P < 0.05
}

# 根据P值生成显著性标记矩阵
sig_markers <- significance_markers(p_values)

# 可视化相关性矩阵并标注显著性标记
corrplot(cor_matrix, 
         method = "circle", 
         type = "upper", 
         tl.cex = 0.8,
         p.mat = p_values, 
         sig.level = c(0.001, 0.01, 0.05), 
         insig = "label_sig", 
         pch.cex = 0.8, 
         pch.col = "black")


#----------------------------RCS分析------------------------------#
library(rms)
# 使用 ols 函数拟合 RCS 模型
# 首先需要定义一个 datadist 对象
colnames(data)
dd <- datadist(data)
options(datadist = "dd")

colnames(data)
#"FATA", "FATM", "FATV", 
#"TOBMC","TOBMD","TST","EST",
#"TST","EST","SHBG","HSCRP"
dependent_vars <- c("VFV", "`Total BMD`", "Testo", "E2", "SHBG", "`hs-CRP`")
independent_var <- "`Vit D`"

# 定义节点数
knots <- 4

library(gridExtra)  # 用于拼接多个 ggplot 图形

# 初始化一个列表存储每个子图
plot_list <- list()

for (dep_var in dependent_vars) {
  # 构建公式
  formula <- as.formula(paste(independent_var, "~ rcs(", dep_var, ",", knots, ")"))
  
  # 拟合 RCS 模型
  fit <- ols(formula, data = data)
  
  # 获取 ANOVA 表
  anova_table <- as.data.frame(anova(fit))
  
  # 提取总 P 值和非线性 P 值
  total_p_value <- anova_table[3, 5]
  nonlinear_p_value <- anova_table[2, 5]
  
  # 预测数据
  pred_data <- Predict(fit)
  
  # 动态计算 P 值标签位置
  x_pos <- max(pred_data[, 1]) - max(pred_data[, 1]) / 8  # x 轴最大值
  y_pos <- max(pred_data$yhat) + max(pred_data$yhat) / 8  # y 轴最大值
  
  # 生成 P 值标签
  p_label <- paste0("P (Total) = ", round(total_p_value, 4), 
                    "\nP (Nonlinear) = ", round(nonlinear_p_value, 3))
  
  # 绘制单个子图
  p <- ggplot(pred_data) +
    annotate("text", x = x_pos*0.6, y = y_pos, label = p_label, hjust = 0, vjust = 1, size = 5) +  # 加大标注文本字体
    theme_classic() +
    theme(
      text = element_text(size = 25),  # 设置全局字体大小为 25
      axis.text = element_text(size = 20, face = "bold"),  # 坐标轴刻度字体加大并加粗
      axis.title = element_text(size = 25, face = "bold"),  # 坐标轴标题字体加大并加粗
      panel.background = element_rect(fill = "#E6F7FF", color = NA)  # 浅蓝色背景，无边框
    )
  
  # 将当前子图添加到列表中
  plot_list[[dep_var]] <- p
}

# 使用 grid.arrange 拼接子图
do.call(grid.arrange, c(plot_list, ncol = 3, nrow = 2))

# 或者保存拼接后的图为 PNG 文件
ggsave("combined_rcs_plot.png", arrangeGrob(grobs = plot_list, ncol = 3, nrow = 2), width = 15, height = 10)

