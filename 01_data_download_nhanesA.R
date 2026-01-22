library(nhanesA)
nhanesSearch("WTMEC2YR")   # 查看所有周期 MEC 权重
nhanesSearch("weight")     # 关键词检索

#Demographics (DEMO) - Dietary (DIET) - Examination (EXAM) - Laboratory (LAB) - Questionnaire (Q).
nhanesTables('LAB', 2013)
nhanesTableVars('LAB', 'BMX_J')
EXAM.DXX <- nhanesTableVars('EXAM', 'DXX_J')
write.csv(EXAM.DXX, "EXAM.DXX.字典.csv", row.names = FALSE)


# A: 1999
demo <- nhanes('DEMO')
bpx <- nhanes('BPX')
bmx <- nhanes('BMX')
dxxag <- nhanes('DXXAG')
VID <- nhanes('VID_A')
VITAEC <- nhanes('VITAEC_A')
TRIGLY <- nhanes('TRIGLY_A')
HDL <- nhanes('HDL_A')
TCHOL <- nhanes('TCHOL_A')
GLU <- nhanes('GLU_A')
CRP <- nhanes('CRP_A')
PAQ <- nhanes("PAQ_A")
ALB <- nhanes("BIOPRO_A")

data_list <- list(demo, bpx, bmx, dxxag, VID, VITAEC, TRIGLY, HDL, TCHOL, GLU, CRP,PAQ,ALB)
data_list <- Filter(Negate(is.null), data_list)
demo99 <- Reduce(function(x, y) merge(x, y, by = "SEQN", all = TRUE), data_list)

# B: 2001
demo_b <- nhanes('DEMO_B')
bpx_b <- nhanes('BPX_B')
bmx_b <- nhanes('BMX_B')
dxxag_b <- nhanes('DXXAG_B')
VID_B <- nhanes('VID_B')
VITAEC_B <- nhanes('VITAEC_B')
TRIGLY_B <- nhanes('TRIGLY_B')
HDL_B <- nhanes('HDL_B')
TCHOL_B <- nhanes('TCHOL_B')
GLU_B <- nhanes('GLU_B')
CRP_B <- nhanes('CRP_B')
PAQ_B <- nhanes("PAQ_B")
ALB_B <- nhanes("BIOPRO_B")
data_list_b <- list(demo_b, bpx_b, bmx_b, dxxag_b, VID_B, VITAEC_B, TRIGLY_B, HDL_B, TCHOL_B, GLU_B, CRP_B,PAQ_B,ALB_B)
data_list_b <- Filter(Negate(is.null), data_list_b)
demo01 <- Reduce(function(x, y) merge(x, y, by = "SEQN", all = TRUE), data_list_b)

# C: 2003
demo_c <- nhanes('DEMO_C')
bpx_c <- nhanes('BPX_C')
bmx_c <- nhanes('BMX_C')
dxxag_c <- nhanes('DXXAG_C')
VID_C <- nhanes('VID_C')
VITAEC_C <- nhanes('VITAEC_C')
TRIGLY_C <- nhanes('TRIGLY_C')
HDL_C <- nhanes('HDL_C')
TCHOL_C <- nhanes('TCHOL_C')
GLU_C <- nhanes('GLU_C')
CRP_C <- nhanes('CRP_C')
PAQ_C <- nhanes("PAQ_C")
ALB_C <- nhanes("BIOPRO_C")
data_list_c <- list(demo_c, bpx_c, bmx_c, dxxag_c, VID_C, VITAEC_C, TRIGLY_C, HDL_C, TCHOL_C, GLU_C, CRP_C,PAQ_C,ALB_C)
data_list_c <- Filter(Negate(is.null), data_list_c)
demo03 <- Reduce(function(x, y) merge(x, y, by = "SEQN", all = TRUE), data_list_c)

# D: 2005
demo_d <- nhanes('DEMO_D')
bpx_d <- nhanes('BPX_D')
bmx_d <- nhanes('BMX_D')
dxxag_d <- nhanes('DXXAG_D')
VID_D <- nhanes('VID_D')
VITAEC_D <- nhanes('VITAEC_D')
TRIGLY_D <- nhanes('TRIGLY_D')
HDL_D <- nhanes('HDL_D')
TCHOL_D <- nhanes('TCHOL_D')
GLU_D <- nhanes('GLU_D')
CRP_D <- nhanes('CRP_D')
DIQ_D <- nhanes('DIQ_D')
DLQ_D <- nhanes('DLQ_D')
FERTIN_D <- nhanes('FERTIN_D')
FETIB_D <- nhanes('FETIB_D')
FOLATE_D <- nhanes('FOLATE_D')
FOLFMS_D <- nhanes('FOLFMS_D')
PAQ_D <- nhanes("PAQ_D")
ALB_D <- nhanes("BIOPRO_D")
data_list_d <- list(demo_d, bpx_d, bmx_d, dxxag_d, VID_D, VITAEC_D, TRIGLY_D, HDL_D, TCHOL_D, GLU_D, CRP_D,
                    DIQ_D,DLQ_D,FERTIN_D,FETIB_D,FOLATE_D,FOLFMS_D,FOLFMS_D,PAQ_D,ALB_D)
data_list_d <- Filter(Negate(is.null), data_list_d)
demo05 <- Reduce(function(x, y) merge(x, y, by = "SEQN", all = TRUE), data_list_d)

# E: 2007
demo_e <- nhanes('DEMO_E')
bpx_e <- nhanes('BPX_E')
bmx_e <- nhanes('BMX_E')
dxxag_e <- nhanes('DXXAG_E')
VID_E <- nhanes('VID_E')
VITAEC_E <- nhanes('VITAEC_E')
TRIGLY_E <- nhanes('TRIGLY_E')
HDL_E <- nhanes('HDL_E')
TCHOL_E <- nhanes('TCHOL_E')
GLU_E <- nhanes('GLU_E')
CRP_E <- nhanes('CRP_E')
DIQ_E <- nhanes('DIQ_E')
DLQ_E <- nhanes('DLQ_E')
FERTIN_E <- nhanes('FERTIN_E')
FETIB_E <- nhanes('FETIB_E')
FOLATE_E <- nhanes('FOLATE_E')
FOLFMS_E <- nhanes('FOLFMS_E')
PAQ_E <- nhanes("PAQ_E")
ALB_E <- nhanes("BIOPRO_E")

data_list_e <- list(demo_e, bpx_e, bmx_e, dxxag_e, VID_E, VITAEC_E, TRIGLY_E, HDL_E, TCHOL_E, GLU_E, CRP_E,
                    DIQ_E,DLQ_E,FERTIN_E,FETIB_E,FOLATE_E,FOLFMS_E,FOLFMS_E,PAQ_E,ALB_E)
data_list_e <- Filter(Negate(is.null), data_list_e)
demo07 <- Reduce(function(x, y) merge(x, y, by = "SEQN", all = TRUE), data_list_e)

# F: 2009
demo_f <- nhanes('DEMO_F')
bpx_f <- nhanes('BPX_F')
bmx_f <- nhanes('BMX_F')
dxxag_f <- nhanes('DXXAG_F')
VID_F <- nhanes('VID_F')
VITAEC_F <- nhanes('VITAEC_F')
TRIGLY_F <- nhanes('TRIGLY_F')
HDL_F <- nhanes('HDL_F')
TCHOL_F <- nhanes('TCHOL_F')
GLU_F <- nhanes('GLU_F')
CRP_F <- nhanes('CRP_F')
DIQ_F <- nhanes('DIQ_F')
DLQ_F <- nhanes('DLQ_F')
FERTIN_F <- nhanes('FERTIN_F')
FETIB_F <- nhanes('FETIB_F')
FOLATE_F <- nhanes('FOLATE_F')
FOLFMS_F <- nhanes('FOLFMS_F')
PAQ_F <- nhanes("PAQ_F")
ALB_F <- nhanes("BIOPRO_F")

data_list_f <- list(demo_f, bpx_f, bmx_f, dxxag_f, VID_F, VITAEC_F, TRIGLY_F, HDL_F, TCHOL_F, GLU_F, CRP_F,
                    DIQ_F,DLQ_F,FERTIN_F,FETIB_F,FOLATE_F,FOLFMS_F,FOLFMS_F,PAQ_F,ALB_F)
data_list_f <- Filter(Negate(is.null), data_list_f)
demo09 <- Reduce(function(x, y) merge(x, y, by = "SEQN", all = TRUE), data_list_f)

# G: 2011
demo_g <- nhanes('DEMO_G')
bpx_g <- nhanes('BPX_G')
bmx_g <- nhanes('BMX_G')
dxxag_g <- nhanes('DXXAG_G')
VID_G <- nhanes('VID_G')
VITAEC_G <- nhanes('VITAEC_G')
TRIGLY_G <- nhanes('TRIGLY_G')
HDL_G <- nhanes('HDL_G')
TCHOL_G <- nhanes('TCHOL_G')
GLU_G <- nhanes('GLU_G')
CRP_G <- nhanes('CRP_G')
DXX_G <- nhanes('DXX_G')
DIQ_G <- nhanes('DIQ_G')
DLQ_G <- nhanes('DLQ_G')
FERTIN_G <- nhanes('FERTIN_G')
FETIB_G <- nhanes('FETIB_G')
FOLATE_G <- nhanes('FOLATE_G')
FOLFMS_G <- nhanes('FOLFMS_G')
TST_G <- nhanes('TST_G')
PAQ_G <- nhanes("PAQ_G")
ALB_G <- nhanes("BIOPRO_G")
data_list_g <- list(demo_g, bpx_g, bmx_g, dxxag_g, VID_G, VITAEC_G, TRIGLY_G, HDL_G, TCHOL_G, GLU_G, CRP_G,DXX_G,
                    DIQ_G,DLQ_G,FERTIN_G,FETIB_G,FOLATE_G,FOLFMS_G,FOLFMS_G,TST_G,PAQ_G,ALB_G)
data_list_g <- Filter(Negate(is.null), data_list_g)
demo11 <- Reduce(function(x, y) merge(x, y, by = "SEQN", all = TRUE), data_list_g)



# H: 2013
demo_h <- nhanes('DEMO_H')
bpx_h <- nhanes('BPX_H')
bmx_h <- nhanes('BMX_H')
dxxag_h <- nhanes('DXXAG_H')
VID_H <- nhanes('VID_H')
VITAEC_H <- nhanes('VITAEC_H')
TRIGLY_H <- nhanes('TRIGLY_H')
HDL_H <- nhanes('HDL_H')
TCHOL_H <- nhanes('TCHOL_H')
GLU_H <- nhanes('GLU_H')
CRP_H <- nhanes('CRP_H')
DXX_H <- nhanes('DXX_H')
DIQ_H <- nhanes('DIQ_H')
DLQ_H <- nhanes('DLQ_H')
FERTIN_H <- nhanes('FERTIN_H')
FETIB_H <- nhanes('FETIB_H')
FOLATE_H <- nhanes('FOLATE_H')
FOLFMS_H <- nhanes('FOLFMS_H')
TST_H <- nhanes('TST_H')
PAQ_H <- nhanes("PAQ_H")
ALB_H <- nhanes("BIOPRO_H")
data_list_h <- list(demo_h, bpx_h, bmx_h, dxxag_h, VID_H, VITAEC_H, TRIGLY_H, HDL_H, TCHOL_H, GLU_H, CRP_H,DXX_H,
                    DIQ_H,DLQ_H,FERTIN_H,FETIB_H,FOLATE_H,FOLFMS_H,FOLFMS_H,TST_H,PAQ_H,ALB_H)
data_list_h <- Filter(Negate(is.null), data_list_h)
demo13 <- Reduce(function(x, y) merge(x, y, by = "SEQN", all = TRUE), data_list_h)

# I: 2015
demo_i <- nhanes('DEMO_I')
bpx_i <- nhanes('BPX_I')
bmx_i <- nhanes('BMX_I')
dxxag_i <- nhanes('DXXAG_I')
VID_I <- nhanes('VID_I')
VITAEC_I <- nhanes('VITAEC_I')
TRIGLY_I <- nhanes('TRIGLY_I')
HDL_I <- nhanes('HDL_I')
TCHOL_I <- nhanes('TCHOL_I')
GLU_I <- nhanes('GLU_I')
CRP_I <- nhanes('HSCRP_I')
DXX_I  <- nhanes('DXX_I')
DIQ_I <- nhanes('DIQ_I')
DLQ_I <- nhanes('DLQ_I')
FERTIN_I <- nhanes('FERTIN_I')
FETIB_I <- nhanes('FETIB_I')
FOLATE_I <- nhanes('FOLATE_I')
FOLFMS_I <- nhanes('FOLFMS_I')
TST_I <- nhanes('TST_I')
PAQ_I <- nhanes("PAQ_I")
ALB_I <- nhanes("BIOPRO_I")
data_list_i <- list(demo_i, bpx_i, bmx_i, dxxag_i, VID_I, VITAEC_I, TRIGLY_I, HDL_I, TCHOL_I, GLU_I, CRP_I,DXX_I,
                    DIQ_I,DLQ_I,FERTIN_I,FETIB_I,FOLATE_I,FOLFMS_I,FOLFMS_I,TST_I,PAQ_I,ALB_I)
data_list_i <- Filter(Negate(is.null), data_list_i)
demo15 <- Reduce(function(x, y) merge(x, y, by = "SEQN", all = TRUE), data_list_i)

# J: 2017
demo_j <- nhanes('DEMO_J')
bpx_j <- nhanes('BPX_J')
bmx_j <- nhanes('BMX_J')
dxxag_j <- nhanes('DXXAG_J')
VID_J <- nhanes('VID_J')
VITAEC_J <- nhanes('VITAEC_J')
TRIGLY_J <- nhanes('TRIGLY_J')
HDL_J <- nhanes('HDL_J')
TCHOL_J <- nhanes('TCHOL_J')
GLU_J <- nhanes('GLU_J')
CRP_J <- nhanes('HSCRP_J')
DXX_J   <- nhanes('DXX_J')
DIQ_J <- nhanes('DIQ_J')
DLQ_J <- nhanes('DLQ_J')
FERTIN_J <- nhanes('FERTIN_J')
FETIB_J <- nhanes('FETIB_J')
FOLATE_J <- nhanes('FOLATE_J')
FOLFMS_J <- nhanes('FOLFMS_J')
PAQ_J <- nhanes("PAQ_J")
ALB_J <- nhanes("BIOPRO_J")
data_list_j <- list(demo_j, bpx_j, bmx_j, dxxag_j, VID_J, VITAEC_J, TRIGLY_J, HDL_J, TCHOL_J, GLU_J, CRP_J,DXX_J,
                    DIQ_J,DLQ_J,FERTIN_J,FETIB_J,FOLATE_J,FOLFMS_J,FOLFMS_J,PAQ_J,ALB_J)
data_list_j <- Filter(Negate(is.null), data_list_j)
demo17 <- Reduce(function(x, y) merge(x, y, by = "SEQN", all = TRUE), data_list_j)

merged_data_list <- list(demo01,demo03,demo05,demo07,demo09,demo11,demo13,demo15,demo17,demo99)


# 定义一个函数来合并两个 data.frame 并处理列名不一致的情况
merge_two_dfs <- function(df1, df2) {
  # 获取两个 data.frame 的列名
  cols_df1 <- colnames(df1)
  cols_df2 <- colnames(df2)
  
  # 找到两个 data.frame 的所有列名（去重后的并集）
  all_cols <- union(cols_df1, cols_df2)
  
  # 检查 df1 是否缺少某些列，如果缺少则添加并填充 "NA"
  missing_cols_df1 <- setdiff(all_cols, cols_df1)
  if (length(missing_cols_df1) > 0) {
    for (col in missing_cols_df1) {
      df1[[col]] <- NA
    }
  }
  
  # 检查 df2 是否缺少某些列，如果缺少则添加并填充 "NA"
  missing_cols_df2 <- setdiff(all_cols, cols_df2)
  if (length(missing_cols_df2) > 0) {
    for (col in missing_cols_df2) {
      df2[[col]] <- NA
    }
  }
  
  # 确保两个 data.frame 的列顺序一致
  df1 <- df1[, all_cols]
  df2 <- df2[, all_cols]
  
  # 合并两个 data.frame
  merged_df <- rbind(df1, df2)
  
  return(merged_df)
}

# 使用 Reduce 函数逐步合并 list 中的所有 data.frame
final_merged_df <- Reduce(merge_two_dfs, merged_data_list)
dim(final_merged_df)



#-----------------------单位查询---------------------------------#
# 运行 nhanesTableVars 并存储结果
results <- list(
  GLU_J = nhanesTableVars('LAB', 'GLU_J'),
  VID_J = nhanesTableVars('LAB', 'VID_J'),
  VITAEC_J = nhanesTableVars('LAB', 'VITAEC_J'),
  TRIGLY_J = nhanesTableVars('LAB', 'TRIGLY_J'),
  HDL_J = nhanesTableVars('LAB', 'HDL_J'),
  TCHOL_J = nhanesTableVars('LAB', 'TCHOL_J'),
  HSCRP_J = nhanesTableVars('LAB', 'HSCRP_J'),
  HSCRP_J = nhanesTableVars('LAB', 'HSCRP_J'),
  DIQ_E <- nhanesTableVars('Q', 'DIQ_E'),
  FERTIN_E <- nhanesTableVars('LAB', 'FERTIN_E'),
  FOLATE_E <- nhanesTableVars('LAB', 'FOLATE_E'),
  FOLFMS_E <- nhanesTableVars('LAB', 'FOLFMS_E'),
  nhanesTableVars('EXAM', 'DXXAG_J')
)


combined_results <- bind_rows(results, .id = "Variable")
write.csv(combined_results, "nhanes_results_summary_单位.csv", row.names = FALSE)

