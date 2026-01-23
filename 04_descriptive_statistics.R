library(readr)
library(MatchIt)
library(nnet)
library(car)
library(mvtnorm)
library(broom)
library(dplyr)
library(ggpubr)
library(gtable)
library(grid)
library(stats)
library(reportROC)
library(pROC)
library(reshape2)
library(forestplot)
library(descriptr)
library(survival)
library(pastecs)
library(compareGroups)
library(stringr)
library(CBCgrps)
library(openxlsx)
library(boot)
library(devEMF)
library(tableone)
library(mgcv)
library(dplyr)
library(labelled)
library(tableone)
library(survey)
library(dplyr)

data_orignal <- read_csv("C:/Users/loveweiyi/Desktop/data.csv")

data <- data_orignal
data <- data %>% mutate(across(c(1,2,4,5,40), as.factor))
data <- as.data.frame(data)
str(data)
# data <- data[,-c(12:15)]

data <- data[data$RIDAGEYR > 7 & data$RIDAGEYR <= 18, ]
#   维生素D三等分
# quantiles <- quantile(data$LBXVIDMS, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
# print(quantiles)


#----#代谢综合征8-12岁没有人

#------------------------------表1表2绘制---------------------------------------#
## 1. 构造 2 岁一组的新变量 age_grp2
data <- data %>% 
  mutate(
    RIDAGEYR = case_when(
      RIDAGEYR %in% 8:9  ~ "8-9",
      RIDAGEYR %in% 10:11 ~ "10-11",
      RIDAGEYR %in% 12:13 ~ "12-13",
      RIDAGEYR %in% 14:15 ~ "14-15",
      RIDAGEYR %in% 16:17 ~ "16-17"
    )
  ) %>% 
  mutate(RIDAGEYR = factor(RIDAGEYR,
                           levels = c("8-9","10-11","12-13","14-15","16-17")))

table(data$RIDAGEYR)

data <- data %>% mutate(across(c(1:5), as.factor))
data <- as.data.frame(data)
data$group2 <- factor(data$Group, levels = c("消瘦", "正常", "代谢综合征", "肥胖"))
data$group3 <- factor(data$Group, levels = c("正常","消瘦", "肥胖"))

colnames(data)
data <- data %>%
  mutate(
    VD_Tertiles = ntile(LBXVIDMS, 3),
    VD_Tertiles = factor(VD_Tertiles, labels = c("T1", "T2", "T3"))
  )

table1 <- CreateTableOne(vars = c("RIAGENDR", "RIDAGEYR", "RIDRETH1",   "LBDTRSI", "LBDHDDSI", "LBDTCSI",
                                  "LBDGLUSI", "LBDFERSI", "LBXVIDMS",  "DXXVFATV","DXDTOBMC", "DXDTOBMD",  
                                  "LBXTST", "LBXEST", "LBXSHBG", "LBXHSCRP", "PAD680"),
                         strata = "group3",
                         data = data,
                         test = TRUE)

table1 <- print(table1, showAllLevels = TRUE)
#write.xlsx(table1,file = "./table11.xlsx",colNames = T,rowNames= T)

#"RIDAGEYR","RIDRETH1","income",,"DXXHEA"
#"LBDRBF", "LBXRBFSI", "LBDFOL", "LBXFOLSI", "LBXIN", "LBDINSI"
colnames(data)
varlist= c( "RIAGENDR","RIDAGEYR","BMXBMI","RIDRETH1","income",
            "LBDTRSI", "LBDHDDSI", "LBDTCSI", "LBDGLUSI", "LBDFERSI", "LBXVIDMS", "LBXVD2MS", "LBXVD3MS", 
            "DXXVFATA", "DXXVFATM", "DXXVFATV", "DXXHEBMC", "DXXHEBMD", "DXXLLBMD", 
            "DXXLSBMD", "DXXPEBMD", "DXDTOBMC", "DXDTOBMD", 
            "LBXTST", "LBXEST", "LBXSHBG", "LBXHSCRP")


table(data$RIDAGEYR)
#8-9 10-11 12-13 14-15 16-17 
#334   273   289   264   255
# 定义变量范围
v <- c("8-9","10-11","12-13","14-15","16-17")

# 循环处理每个 var 值
for (i in c(1:5)) {
  # 筛选数据框中 RIDAGEYR 等于当前 v 的行
  var <- v[i]
  df <- data[data$RIDAGEYR %in% var, ]
  
  # 调用 multigrps 函数进行分组分析
  b <- multigrps(
    df, 
    "group3", 
    varlist = c(
      "RIAGENDR", "RIDAGEYR","BMXBMI", "LBXVIDMS", "LBXVD2MS","LBXVD3MS", "DXXVFATA", "DXXVFATM", 
      "DXXVFATV", "DXXHEBMC", "DXXHEBMD", "DXXLLBMD", "DXXLSBMD", "DXXPEBMD", 
      "DXDTOBMC", "DXDTOBMD", "LBXTST", "LBXEST", "LBXSHBG", "LBXHSCRP"
    ), 
    ShowStatistic = TRUE
  )
  
  # 打印结果
  print(b, quote = TRUE)
  
  # 将结果写入 Excel 文件
  write.xlsx(
    data.frame(b), 
    file = paste0("table_", v, ".xlsx"), 
    colNames = FALSE, 
    rowNames = FALSE
  )
}



#-----------------------亚组分析------------------------------------------#
table(data$group3)
data$group2 <- factor(data$group2, levels = c("正常","消瘦",  "代谢综合征", "肥胖"))
data$group3 <- factor(data$Group, levels = c("正常","消瘦", "肥胖"))
varlist <- c("RIAGENDR","BMXBMI","LBXVIDMS", "LBXVD3MS", "DXXVFATA", "DXXVFATM", "DXXVFATV", 
             "DXXHEBMC", "DXXHEBMD", "DXXLLBMD", "DXXLSBMD", "DXXPEBMD","DXDTOBMC","DXDTOBMD", 
             "LBXTST", "LBXEST", "LBXSHBG", "LBXHSCRP")


# 初始化存储结果的列表
vad <- c(8:16)

for (v in vad) {
  df <- data[data$RIDAGEYR %in% c(v), ]
  
  results_list <- list()
  # 遍历每个变量
  for (var in varlist) {
    # 构建多分类逻辑回归模型
    formula <- as.formula(paste("group3 ~", var))
    model_multinom <- nnet::multinom(formula, data = df, na.action = na.omit)
    
    # 提取模型摘要
    model_summary <- summary(model_multinom)
    
    # 计算OR值及95%置信区间
    or_values <- exp(coef(model_multinom))
    ci_values <- exp(confint(model_multinom))
    
    # 显著性检验：提取p值
    anova_test <- Anova(model_multinom, type = "III")
    p_value <- anova_test$`Pr(>Chisq)`[1]  # 提取变量的p值
    
    # 初始化置信区间的下限和上限
    ci_lower <- c()
    ci_upper <- c()
    
    # 遍历每个类别（假设类别为：消瘦、代谢综合征、肥胖）
    categories <- c("消瘦",  "肥胖")
    for (category in categories) {
      ci_lower <- c(ci_lower, ci_values[2, "2.5 %", category])
      ci_upper <- c(ci_upper, ci_values[2, "97.5 %", category])
    }
    
    # 展平 OR 值
    or_values_flat <- as.vector(or_values[, -1])  # 排除截距项
    
    # 创建结果数据框
    result <- data.frame(
      Variable = rep(var, length(categories)),
      Category = categories,
      Estimate = round(or_values_flat, digits = 2),
      `2.5 %` = round(ci_lower, digits = 2),
      `97.5 %` = round(ci_upper, digits = 2),
      P_Value = round(rep(p_value, length(categories)), digits = 4)  # 添加p值
    )
    
    # 将结果添加到列表中
    results_list[[var]] <- result
  }
  
  # 合并所有结果
  final_result <- do.call(rbind, results_list)
  print(final_result, quote = TRUE)
  # 写入Excel文件
  #write.xlsx(data.frame(final_result),  file = paste0("table_2_", v, ".xlsx"), , colNames = TRUE, rowNames = FALSE)
}












#----------------------------倾向性匹配-----------------------------------------@
# 创建一个新的二分类变量
data$Obese_vs_Normal <- ifelse(data$Group == "肥胖", 1, 
                               ifelse(data$Group == "正常", 0, NA))

# 去掉 NA 值
data1 <- data[!is.na(data$Obese_vs_Normal),]

# 执行匹配
match_model <- matchit(Obese_vs_Normal ~ RIDAGEYR, 
                       data = data1, 
                       method = "nearest", 
                       ratio = 5)  # 1:1 匹配

# 查看匹配结果
summary(match_model)

# 获取匹配后的数据集
matched_data <- match.data(match_model)

# 2. 检查匹配质量
# 可以查看标准化均值差(SMD)，确保匹配后SMD < 0.1
summary(match_model)$sum.all # 匹配前的平衡情况
summary(match_model)$sum.matched # 匹配后的平衡情况

matched_data <- matched_data[,colnames(matched_data)%in%colnames(data)]
data2 <- data[is.na(data$Obese_vs_Normal),]
data <- rbind(matched_data,data2)
