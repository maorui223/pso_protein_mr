setwd('E:/UKB/MR/micro/r_0.001/pso_protein')
exposure_finn = read.csv('./finn/exposure_finn.csv')

pso_ukb = fread('PheCode_696.4_SAIGE_MACge20.txt.vcf.gz')
rm(pso_train,pso_ukb)
gc()


3225*7.1968
#################################normal analysis
chd_out_dat_finn_pso_ukb1 <- format_data(pso_ukb,
                                        type="outcome",
                                        snp_col = "ID",
                                        snps = exposure_finn$SNP,
                                        beta_col = "beta",
                                        se_col = "sebeta",
                                        eaf_col = "af",
                                        effect_allele_col = "ALT",
                                        other_allele_col = "REF",
                                        pval_col = "pval",
                                        ncase_col = 'num_cases',
                                        ncontrol_col = 'num_controls',
                                        chr_col = '#CHROM',
                                        pos_col = 'POS')

dat_finn_pso_ukb <- harmonise_data(exposure_finn, chd_out_dat_finn_pso_ukb1)
dat_finn_pso_ukb <- harmonise_data_modifed(exposure_finn, chd_out_dat_finn_pso_ukb,r2_thershold = 0.9, build = "37", pop = "EUR")

dat_finn_pso_ukb = dat_finn_pso_ukb %>%
  mutate(R2 = 2*(1-eaf.exposure)*eaf.exposure*(beta.exposure^2))
dat_finn_pso_ukb = dat_finn_pso_ukb %>%
  distinct(exposure,SNP,beta.exposure,eaf.exposure, .keep_all = T)
heterogeneity_finn_pso_ukb = mr_heterogeneity(dat_finn_pso_ukb)
pleiotropy_finn_pso_ukb = mr_pleiotropy_test(dat_finn_pso_ukb)
write.csv(heterogeneity_finn_pso_ukb, './finn/heterogeneity_finn_pso_ukb.csv')
write.csv(pleiotropy_finn_pso_ukb, './finn/pleiotropy_finn_pso_ukb.csv')

list = mr_method_list()
list1 = list[c(1,8),]$obj
dat_finn_pso_ukb$samplesize.exposure = 10708
res_finn_pso_ukb <- mr_modified(dat_finn_pso_ukb, prop_var_explained = T)
res_finn_pso_ukb1 = res_finn_pso_ukb %>%
  filter(pval > 0)
length(unique(res_finn_pso_ukb1$exposure))
write.csv(res_finn_pso_ukb1, './finn/res_finn_pso_ukb.csv')
res_finn_pso_ukb2 = res_finn_pso_ukb1 %>%
  filter(pval < 0.05/1320)
res_pso_ukb_merge = merge(res_finn_pso_ukb2, res_finn_pso_ukb2,by.x = 'exposure', by.y  = 'id.exposure')
res_finn_pso_ukb2$id.exposure = res_finn_pso_ukb2$exposure
final_result = generate_odds_ratios(final_result)
APOF_finn_pso_ukb = subset(final_result, id.exposure == 'rs1274496')
APOF_finn_pso_ukb$pve = 0.11
APOF_finn_pso_ukb <- APOF_finn_pso_ukb[, c(1:9, 15, 10:14)]
res_finn_pso_ukb1 = rbind(res_finn_pso_ukb1, APOF_finn_pso_ukb)
dat_finn_pso_ukb$samplesize.outcome = 2237+398199
dat_finn_pso_ukb$ncase.outcome = 2237
dat_finn_pso_ukb$ncontrol.outcome = 398199

steiger_pso_ukb = dat_finn_pso_ukb %>%
  steiger_filtering()
dat_finn_pso_ukb$steiger_pval = steiger_pso_ukb$steiger_pval
write.csv(dat_finn_pso_ukb, '../../../finn/dat_finn_pso_ukb.csv')

###################################APOF
a = find_proxy(snp = 'rs2643623')
chd_out_dat_finn_pso_ukb_apof <- format_data(pso_ukb,
                                             type="outcome",
                                             snp_col = "ID",
                                             snps = a,
                                             beta_col = "beta",
                                             se_col = "sebeta",
                                             eaf_col = "af",
                                             effect_allele_col = "ALT",
                                             other_allele_col = "REF",
                                             pval_col = "pval",
                                             ncase_col = 'num_cases',
                                             ncontrol_col = 'num_controls',
                                             chr_col = '#CHROM',
                                             pos_col = 'POS')

APOF = subset(exposure_finn, exposure == 'APOF')
a = as.data.frame(a)
n <- 134 # 重复的次数
# 选择要重复的行
row_to_repeat <- APOF[1, ]
# 重复选择的行
APOF <- do.call("rbind", lapply(1:n, function(x) row_to_repeat))
APOF$SNP = a$a
APOF = subset(APOF, SNP %in% chd_out_dat_finn_pso_ukb_apof$SNP) 
APOF$samplesize.exposure = 10708
APOF$eaf.exposure = 0.9324
APOF = cbind(APOF, chd_out_dat_finn_pso_ukb_apof[,1:6])
colnames(APOF) <- gsub("outcome", "exposure", colnames(APOF))
APOF$id.exposure = APOF$SNP
# 初始化一个空的list来存储每行的结果
list_of_results <- list()
# 对APOF数据框的每一行进行循环
for(i in 1:nrow(APOF)){
  # 使用try()函数包裹可能导致错误的代码
  result <- try({
    # 提取当前行
    current_row <- APOF[i, ]
    # 和其他数据框进行整合
    harmonised_data <- harmonise_data(current_row, chd_out_dat_finn_pso_ukb_apof)
    # 进行MR分析
    TwoSampleMR::mr(harmonised_data, method_list = list1)
  }, silent = TRUE)  # silent = TRUE 可以阻止错误消息的打印
  # 检查结果是否为"try-error"类，如果不是，将其存储到list中
  if(!inherits(result, "try-error")) {
    list_of_results[[i]] <- result
  } else {
    # 如果你想保存出错的信息，你可以这样做
    list_of_results[[i]] <- list(error = TRUE, message = as.character(result))
  }
}
# 如果需要的话，将所有结果合并为一个数据框
final_result <- do.call("rbind", list_of_results)
write.csv(final_result, './finn/APOF_mr.csv')




# chd_out_dat_finn_pso_train <- format_data(pso_train,
#                                           type="outcome", 
#                                           phenotype_col = 'pso',
#                                           snp_col = "rsid",
#                                           snps = exposure_finn$SNP,
#                                           beta_col = "beta",
#                                           se_col = "se",
#                                           #eaf_col = "af_alt",
#                                           effect_allele_col = "effect_allele",
#                                           other_allele_col = "other_allele",
#                                           pval_col = "p",
#                                           #gene_col = 'nearest_genes',
#                                           chr_col = 'chrom',
#                                           pos_col = 'pos')
# 
# head(chd_out_dat_finn_pso_train)
# chd_out_dat_finn_pso_train = get_eaf_from_1000G(chd_out_dat_finn_pso_train, path = 'E:/UKB/MR/micro/r_0.001/drug/1kg.v3/EUR/', type = "outcome")
# dat_finn_pso_train <- harmonise_data(exposure_finn, chd_out_dat_finn_pso_train)
# dat_finn_pso_train = dat_finn_pso_train %>%
#   mutate(R2 = 2*(1-eaf.exposure)*eaf.exposure*(beta.exposure^2))
# dat_finn_pso_train = dat_finn_pso_train %>%
#   distinct(exposure,SNP,beta.exposure,eaf.exposure, .keep_all = T)
# heterogeneity_finn_pso_train = mr_heterogeneity(dat_finn_pso_train)
# pleiotropy_finn_pso_train = mr_pleiotropy_test(dat_finn_pso_train)
# write.csv(heterogeneity_finn_pso_train, './finn/heterogeneity_finn_pso_train.csv')
# write.csv(pleiotropy_finn_pso_train, './finn/pleiotropy_finn_pso_train.csv')
# 
# list = mr_method_list()
# list1 = list[c(1,8),]$obj
# dat_finn_pso_train$samplesize.exposure = 10708
# res_finn_pso_train <- mr_modified(dat_finn_pso_train, prop_var_explained = T)
# res_finn_pso_train1 = res_finn_pso_train %>%
#   filter(pval > 0)
# length(unique(res_finn_pso_train1$exposure))
# write.csv(res_finn_pso_train1, './finn/res_finn_pso_train.csv')
# res_finn_pso_train2 = res_finn_pso_train1 %>%
#   filter(pval < 0.05/100)
# res_pso_train_merge = merge(res_finn_pso_train2, res_finn_pso_train2,by.x = 'exposure', by.y  = 'id.exposure')
# res_finn_pso_train2$id.exposure = res_finn_pso_train2$exposure
