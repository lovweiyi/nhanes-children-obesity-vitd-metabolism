library(dplyr)
# 定义计算分组百分位数的函数
# 计算各RIDAGEYRRIAGENDR的BMXWAIST、血压分位数
colnames(dat)
result1 <- dat %>%
  group_by(RIDAGEYR, RIAGENDR) %>%
  summarise(
    waist_q25 = quantile(BMXWAIST, 0.25, na.rm = TRUE),
    waist_q50 = quantile(BMXWAIST, 0.50, na.rm = TRUE),
    waist_q75 = quantile(BMXWAIST, 0.75, na.rm = TRUE),
    waist_q90 = quantile(BMXWAIST, 0.90, na.rm = TRUE)
  ) %>%
  ungroup()
result1
result2 <- dat %>%
  group_by(RIDAGEYR, RIAGENDR) %>%
  summarise(
    sbp_q25 = quantile(BPXSY1, 0.25, na.rm = TRUE),
    sbp_q50 = quantile(BPXSY1, 0.50, na.rm = TRUE),
    sbp_q75 = quantile(BPXSY1, 0.75, na.rm = TRUE),
    sbp_q90 = quantile(BPXSY1, 0.90, na.rm = TRUE)
  ) %>%
  ungroup()
result2
result3 <- dat %>%
  group_by(RIDAGEYR, RIAGENDR) %>%
  summarise(
    dbp_q25 = quantile(BPXDI1, 0.25, na.rm = TRUE),
    dbp_q50 = quantile(BPXDI1, 0.50, na.rm = TRUE),
    dbp_q75 = quantile(BPXDI1, 0.75, na.rm = TRUE),
    dbp_q90 = quantile(BPXDI1, 0.90, na.rm = TRUE)
  ) %>%
  ungroup()
result3

# ---------------------------- 阈值判断函数优化 ---------------------------- #
# 中心性肥胖判断（保持不变）
check_obesity <- function(age, sex, waist) {
  case_when(
    age == 6 & sex == "Male"   ~ waist >= 66,
    age == 6 & sex == "Female"   ~ waist >= 66.4,
    age == 7 & sex == "Male"   ~ waist >= 72.5,
    age == 7 & sex == "Female" ~ waist >= 73.2,
    age == 8 & sex == "Male"   ~ waist >= 78.4,
    age == 8 & sex == "Female" ~ waist >= 78.4,
    age == 9 & sex == "Male"   ~ waist >= 83.4,
    age == 9 & sex == "Female" ~ waist >= 83.4,
    age == 10 & sex == "Male"   ~ waist >= 87.1,
    age == 10 & sex == "Female"   ~ waist >= 86.1,
    age == 11 & sex == "Male"   ~ waist >= 92.5,
    age == 11 & sex == "Female" ~ waist >= 92,
    age == 12 & sex == "Male"   ~ waist >= 95.5,
    age == 12 & sex == "Female" ~ waist >= 95.9,
    age == 13 & sex == "Male"   ~ waist >= 99.2,
    age == 13 & sex == "Female" ~ waist >= 98.6,
    age == 14 & sex == "Male"   ~ waist >= 102,
    age == 14 & sex == "Female" ~ waist >= 99,
    age == 15 & sex == "Male"   ~ waist >= 103.1,
    age == 15 & sex == "Female" ~ waist >= 100.4,
    age >= 16 & sex == "Male"   ~ waist >= 106.2,
    age >= 16 & sex == "Female" ~ waist >= 102.1,
    age >= 17 & sex == "Male"   ~ waist >= 104.9,
    age >= 17 & sex == "Female" ~ waist >= 104.1,
    age >= 18 & sex == "Male"   ~ waist >= 107.8,
    age >= 18 & sex == "Female" ~ waist >= 107.5,
    TRUE                        ~ FALSE
  )
}

# 血压判断函数优化
check_bp1 <- function(age, sex, sbp) {
  case_when(
    age == 6 & sex == "Male"   ~ sbp >= 105,
    age == 6 & sex == "Female"   ~ sbp >= 106,
    age == 7 & sex == "Male"   ~ sbp >= 108,
    age == 7 & sex == "Female" ~ sbp >= 106,
    age == 8 & sex == "Male"   ~ sbp >= 112,
    age == 8 & sex == "Female" ~ sbp >= 110,
    age == 9 & sex == "Male"   ~ sbp >= 114,
    age == 9 & sex == "Female" ~ sbp >= 114,
    age == 10 & sex == "Male"   ~ sbp >= 114,
    age == 10 & sex == "Female"   ~ sbp >= 116,
    age == 11 & sex == "Male"   ~ sbp >= 116,
    age == 11 & sex == "Female" ~ sbp >= 116,
    age == 12 & sex == "Male"   ~ sbp >= 118,
    age == 12 & sex == "Female" ~ sbp >= 120,
    age == 13 & sex == "Male"   ~ sbp >= 122,
    age == 13 & sex == "Female" ~ sbp >= 120,
    age == 14 & sex == "Male"   ~ sbp >= 124,
    age == 14 & sex == "Female" ~ sbp >= 118,
    age == 15 & sex == "Male"   ~ sbp >= 126,
    age == 15 & sex == "Female" ~ sbp >= 120,
    age >= 16 & sex == "Male"   ~ sbp >= 128,
    age >= 16 & sex == "Female" ~ sbp >= 120,
    age == 17 & sex == "Male"   ~ sbp >= 128,
    age == 17 & sex == "Female" ~ sbp >= 120,
    age >= 18 & sex == "Male"   ~ sbp >= 128,
    age >= 18 & sex == "Female" ~ sbp >= 122,
    TRUE                        ~ FALSE
  )
}

check_bp2 <- function(age, sex, dbp) {
  case_when(
    age == 6 & sex == "Male"   ~ dbp >= 60,
    age == 6 & sex == "Female"   ~ dbp >= 60,
    age == 7 & sex == "Male"   ~ dbp >= 62,
    age == 7 & sex == "Female" ~ dbp >= 64,
    age == 8 & sex == "Male"   ~ dbp >= 66,
    age == 8 & sex == "Female" ~ dbp >= 66,
    age == 9 & sex == "Male"   ~ dbp >= 68,
    age == 9 & sex == "Female" ~ dbp >= 68,
    age == 10 & sex == "Male"   ~ dbp >= 70,
    age == 10 & sex == "Female"   ~ dbp >= 68,
    age == 11 & sex == "Male"   ~ dbp >= 70,
    age == 11 & sex == "Female" ~ dbp >= 70,
    age == 12 & sex == "Male"   ~ dbp >= 70,
    age == 12 & sex == "Female" ~ dbp >= 70,
    age == 13 & sex == "Male"   ~ dbp >= 72,
    age == 13 & sex == "Female" ~ dbp >= 72,
    age == 14 & sex == "Male"   ~ dbp >= 72,
    age == 14 & sex == "Female" ~ dbp >= 72,
    age == 15 & sex == "Male"   ~ dbp >= 74,
    age == 15 & sex == "Female" ~ dbp >= 74,
    age >= 16 & sex == "Male"   ~ dbp >= 76,
    age >= 16 & sex == "Female" ~ dbp >= 74,
    age == 17 & sex == "Male"   ~ dbp >= 76,
    age == 17 & sex == "Female" ~ dbp >= 74,
    age >= 18 & sex == "Male"   ~ dbp >= 76,
    age >= 18 & sex == "Female" ~ dbp >= 74,
    TRUE                        ~ FALSE
  )
}


# ---------------------------- 更新IDF诊断函数 ---------------------------- #
diagnose_IDF <- function(age, sex, waist, 
                         sbp, dbp, 
                         fpg, diabetes, 
                         tg, hdl_c) {
  # 中心性肥胖判断
  obesity <- check_obesity(age, sex, waist)
  # 血压判断（向量化处理）
  hypertension <- check_bp1(age, sex, sbp) | check_bp2(age, sex, dbp)
  # 其他代谢指标判断
  hyperglycemia <- fpg >= 5.6 | diabetes >= 1
  hyperTG <- tg >= 1.7
  lipid <- case_when(
    age >= 10 & age < 16 ~ hdl_c < 1.03,
    sex == "male"        ~ hdl_c < 1.03,
    sex == "female"     ~ hdl_c < 1.29,
    TRUE                ~ FALSE
  )
  # 综合诊断
  obesity & (hypertension + hyperglycemia + hyperTG + lipid >= 2)
}

# ---------------------------- 执行诊断 ---------------------------- #
diagnosis_result <- diagnose_IDF(
  dat$RIDAGEYR,
  dat$RIAGENDR,
  dat$BMXWAIST,
  dat$BPXSY1,
  dat$BPXDI1,
  dat$LBXGLU,
  dat$胰岛素,  # 确保这是糖尿病诊断状态
  dat$甘油三酯,
  dat$LBDHDD
)

# 将结果添加到数据框
dat$IDF_diagnosis <- diagnosis_result


# ---------------------------- Cook标准判断函数 ---------------------------- #
diagnose_Cook <- function(age, sex, waist, 
                          sbp, dbp, 
                          fpg, tg, hdl_c) {
  # 中心性肥胖判断
  obesity <- check_obesity(age, sex, waist)
  
  # 血压判断（向量化处理）
  hypertension <- check_bp1(age, sex, sbp) | check_bp2(age, sex, dbp)
  
  # 其他代谢指标判断
  hyperglycemia <- fpg >= 5.1
  hyperTG <- tg >= 1.24
  lipid <- hdl_c < 1.03  # 注意这里修正为小于号，因为LBDHDD低于阈值才是异常
  
  # 综合诊断（使用逐元素逻辑运算符 &）
  result <-   rowSums(cbind(obesity,hypertension, hyperglycemia, hyperTG, lipid)) >= 2
  
  return(result)
}

check_bp2(dat$RIDAGEYR, dat$RIAGENDR, dat$BPXDI1)
dat$Cook_result <- diagnose_Cook(
  dat$RIDAGEYR,
  dat$RIAGENDR,
  dat$BMXWAIST,
  dat$BPXSY1,
  dat$BPXDI1,
  dat$LBDGLUSI,
  dat$LBDTRSI,
  dat$LBDHDDSI
)

table(dat$Cook_result)

# ---------------------------- 中国标准判断函数 ---------------------------- #
diagnose_China <- function(age, sex, waist, sbp, dbp, fpg, diabetes, tg, hdl_c, non_hdl_c) {
  # 获取BMXWAIST的P50值（根据RIDAGEYR和RIAGENDR）
  p50_waist <- get_p50_waist(age, sex)
  
  # 获取血压的P90值（根据RIDAGEYR和RIAGENDR）
  p90_sbp <- get_p90_sbp(age, sex)
  p90_dbp <- get_p90_dbp(age, sex)
  
  # 五项判断
  obesity <- waist >= p50_waist
  hypertension <- sbp >= p90_sbp || dbp >= p90_dbp
  hyperglycemia <- fpg >= 5.6 || diabetes
  hyperTG <- tg >= 1.47
  lipid <- hdl_c < 1.03 || non_hdl_c >= 3.76
  
  # 综合诊断
  sum(obesity, hypertension, hyperglycemia, hyperTG, lipid) >= 3
}

write.csv(dat, "dat.csv", row.names = FALSE)



#---------------------------营养不良判断------------------------------------ #
result4 <- dat %>%
  group_by(RIDAGEYR, RIAGENDR) %>%
  summarise(
    BMI_q5 = quantile(BMXBMI, 0.05, na.rm = TRUE),
    BMI_q25 = quantile(BMXBMI, 0.25, na.rm = TRUE),
    BMI_q50 = quantile(BMXBMI, 0.50, na.rm = TRUE),
    BMI_q75 = quantile(BMXBMI, 0.75, na.rm = TRUE),
    BMI_q90 = quantile(BMXBMI, 0.90, na.rm = TRUE),
    BMI_q95 = quantile(BMXBMI, 0.95, na.rm = TRUE),
    BMI_sd = sd(BMXBMI, na.rm = TRUE)  # 添加标准差的计算
  ) %>%
  ungroup()


result4$BMIZ <- (result4$BMI_q5+result4$BMI_q95)/2

# 定义函数：获取RIDAGEYRRIAGENDR对应的参考均值和标准差
get_bmi_ref <- function(RIDAGEYR, RIAGENDR) {
  match_row <- result4[result4$RIDAGEYR == RIDAGEYR & result4$RIAGENDR == RIAGENDR, ]
  if (nrow(match_row) == 1) {
    return(list(BMIZ = match_row$BMIZ, BMI_sd = match_row$BMI_sd,
                BMI_q95 = match_row$BMI_q95,BMI_q5 = match_row$BMI_q5))
  } else {
    return(list(BMIZ = NA, BMI_sd = NA, BMI_q95 = NA, BMI_q5 = NA))  # 返回 NA 用于后续处理
  }
}

# 动态生成营养状态的规则
df <- dat %>%
  rowwise() %>%  # 按行处理
  mutate(
    # 获取参考均值和标准差
    ref = list(get_bmi_ref(RIDAGEYR, RIAGENDR)),
    BMIZ = ref$BMIZ,  # 提取参考均值
    BMI_sd = ref$BMI_sd,  # 提取参考标准差
    BMI_q95 = ref$BMI_q95,
    BMI_q5 = ref$BMI_q5,
    
    # 根据参考值计算营养状态
    nutritional_status = case_when(
      is.na(BMIZ) | is.na(BMI_sd) ~ "无法计算",  # 如果参考值缺失，则无法计算
      # BMXBMI < (BMIZ - 3 * BMI_sd) ~ "重度营养不良",
      # BMXBMI >= (BMIZ - 3 * BMI_sd) & BMXBMI < (BMIZ - 2 * BMI_sd) ~ "中度营养不良",
      # BMXBMI >= (BMIZ - 2 * BMI_sd) & BMXBMI < (BMIZ - 1 * BMI_sd) ~ "轻度营养不良",
      BMXBMI >= BMI_q95 ~ "肥胖",
      BMXBMI <= BMI_q5 ~ "消瘦",
      TRUE ~ "正常"
    )
  ) %>%
  ungroup()         # 取消分组

df <- df[,-c(34)]
table(df$nutritional_status)
#肥胖、消瘦判定，
#install.packages("anthro")  # 安装 anthro 包
library(anthro)             # 加载 anthro 包

# 将性别转换为数值：1 = male, 2 = female
df$RIAGENDR <- ifelse(df$RIAGENDR == "Male", 1, 
                       ifelse(df$RIAGENDR == "Female", 2, NA))

df$month <- df$RIDAGEYR*(365.25/12)
# df$z_scores <- anthro_zscores(
#   sex = dat$RIAGENDR,  # sex：1=male，2=female
#   age = dat$month,  # 将年龄转换为天数
#   weight = dat$BMXWT,
#   lenhei = dat$BMXHT
# )
# dat <-dat[,-c(40:47)]
library(zscorer)  
df <- addWGSR(data = df, sex = "RIAGENDR", firstPart = "BMXWT", 
        secondPart = "BMXHT", thirdPart = "month", index = "bfa", 
        output = "bmiAgeZ", digits = 2)

df <- df[!is.na(df$RIAGENDR),]
df <- df[!is.na(df$RIDAGEYR),]
df <- df[!is.na(df$BMXWT),]
df <- df[!is.na(df$BMXHT),]

df$groups <- ifelse(df$bmiAgeZ > 2, "Overweight/Obese",
             ifelse(df$bmiAgeZ >= -2 & df$bmiAgeZ <= 2, "Normal weight",
             ifelse(df$bmiAgeZ >= -3 & df$bmiAgeZ < -2, "Thinness", "Severe Thinness")))

library(dplyr)
df <- df %>%
    mutate(
    Group = case_when(
      nutritional_status == "肥胖" & Cook_result == TRUE ~ "肥胖合并代谢综合征",
      nutritional_status == "消瘦" & Cook_result == TRUE ~ "消瘦合并代谢综合征",
      nutritional_status == "肥胖" ~ "肥胖",
      nutritional_status == "消瘦" ~ "消瘦",
      Cook_result == TRUE ~ "代谢综合征",
      TRUE ~ "正常"
    )
  )
table(df$Group)
table(df$groups,df$nutritional_status)
table(df$nutritional_status,df$Cook_result)

extdat<- df %>%
  filter(
    (groups == "Normal weight" & nutritional_status == "正常") |
      (groups == "Normal weight" & nutritional_status == "消瘦") |
      (groups == "Overweight/Obese" & nutritional_status == "肥胖") |
      (groups %in% c("Severe Thinness", "Thinness"))
  )
extdat <- extdat[!is.na(extdat$DXXHEBMD),]
extdat <-extdat[,-c(6,34:37)]
table(extdat$nutritional_status)
table(extdat$Group)
table(extdat$groups,extdat$nutritional_status)
table(extdat$nutritional_status,extdat$Cook_result)




# 计算缺失值比例  #
missing_proportion_mergedata <- colSums(is.na(extdat)) / nrow(extdat)
print(missing_proportion_mergedata)
