#SEQN、SDDSRVYR、RIAGENDR、RIDAGEYR、RIDRETH1、BMXWT、BMXHT、BMXBMI、BMXWAIST、BPXSY1、BPXDI1、RIDEXPRG（怀孕）
#LBXGLU、LBXCRP、LBXTR、LBDLDL、LBDHDD、INDHHIN2(家庭收入)、DXXVATV: DXA 内脏脂肪体积、LBXVIDMS、LBXVD2MS、LBXVD3MS
#DXAHEBV：全身骨密度（Bone Density）相关变量。DXXHEA：全身区域的某种面积（Area）。
#DXXHEBMC：全身骨矿物质含量（Bone Mineral Content）。DXXHEBMD：全身骨密度（Bone Mineral Density）。
#DIQ220:糖尿病、LBDFERSI	Ferritin (ug/L)、LBXFER	Ferritin (ng/mL)、LBDFOL	Serum folate (ng/mL)
#LBDRBF	RBC folate (ng/mL)、LBXFOLSI	Serum folate ( nmol/L)、LBXRBFSI	RBC folate (nmol/L)
#LBDSF1LC	5-Methyl-tetrahydrofolate comment code、LBDSF2LC	Folic acid comment code、
#LBXSF1SI	5-Methyl-tetrahydrofolate (nmol/L)、LBXSF2SI	Folic acid (nmol/L)


# 定义需要的列名
selected_columns <- c(
  "SEQN", "SDDSRVYR", "RIAGENDR", "RIDAGEYR", "RIDRETH1", 
  "BMXWT", "BMXHT", "BMXBMI", "BMXWAIST", "BPXSY1", "BPXDI1", 
  "RIDEXPRG", "LBDGLUSI", "LBXVD2MS","LBXVIDMS", "LBDTRSI", "LBDHDDSI", "LBDTCSI", "LBXHSCRP",
  "INDHHIN2", "DXXVFATA", "LBXVIDMS", "LBXVD2MS", "LBXVD3MS", "DXXVFATM", "DXXVFATV","DXXHEA", "DXXHEBMC", "DXXHEBMD",
  "LBDFOL","LBDRBF","LBXFOLSI","LBXRBFSI","LBDFERSI"
)

# 筛选出这些列
df <- final_merged_df[colnames(final_merged_df) %in% selected_columns]

# 检查是否所有列都成功提取
missing_columns <- setdiff(selected_columns, colnames(df))
if (length(missing_columns) > 0) {
  warning("以下列未找到: ", paste(missing_columns, collapse = ", "))
}


df <- df[!is.na(df$LBXVIDMS),]
df$SDDSRVYR <- gsub(".*?(\\d{4}-\\d{4}).*", "\\1", df$SDDSRVYR)

dat <- df[df$RIDAGEYR > 9 & df$RIDAGEYR < 19, ]
dat <- dat[!is.na(dat$RIDAGEYR),]

# 计算缺失值比例  #
missing_proportion_mergedata <- colSums(is.na(dat)) / nrow(dat)
print(missing_proportion_mergedata)





