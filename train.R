library(data.table)
library(TwoSampleMR)
library(tidyverse)
library(coloc)
library(rco)
source('E:/UKB/MR/micro/r_0.001/Get_MR-main/1.0/Get_MR1.0.r')
exposure_train = read.csv('train/plasma_protein.csv')
head(exposure_train)
exposure_train$id.exposure = exposure_train$exposure
pso_train <- fread('tsoi_2012_23143594_pso_efo0000676_1_ichip.sumstats.tsv.gz')
head(pso_train)
chd_out_dat_train_pso_train <- format_data(pso_train, 
                                           type="outcome", 
                                           phenotype_col = 'pso',
                                           snp_col = "rsid",
                                           snps = exposure_train$SNP,
                                           beta_col = "beta",
                                           se_col = "se",
                                           #eaf_col = "af_alt",
                                           effect_allele_col = "effect_allele",
                                           other_allele_col = "other_allele",
                                           pval_col = "p",
                                           #gene_col = 'nearest_genes',
                                           chr_col = 'chrom',
                                           pos_col = 'pos')

chd_out_dat_train_pso_train = get_eaf_from_1000G(chd_out_dat_train_pso_train, path = 'E:/UKB/MR/micro/r_0.001/drug/1kg.v3/EUR/', type = "outcome")
head(chd_out_dat_train_pso_train)
dat_train_pso_train <- harmonise_data(exposure_train, chd_out_dat_train_pso_train)
dat_train_pso_train = dat_train_pso_train %>%
  mutate(R2 = 2*(1-eaf.exposure)*eaf.exposure*(beta.exposure^2))
dat_train_pso_train = dat_train_pso_train %>%
  distinct(exposure,SNP,beta.exposure,eaf.exposure, .keep_all = T)
heterogeneity_train_pso_train = mr_heterogeneity(dat_train_pso_train)
pleiotropy_train_pso_train = mr_pleiotropy_test(dat_train_pso_train)
write.csv(heterogeneity_train_pso_train, './train/heterogeneity_train_pso_train.csv')
write.csv(pleiotropy_train_pso_train, './train/pleiotropy_train_pso_train.csv')

list = mr_method_list()
list1 = list[c(1,8),]$obj
res_train_pso_train <- mr_modified(dat_train_pso_train, prop_var_explained = T)
res_train_pso_train1 = res_train_pso_train %>%
  filter(pval > 0)
length(unique(res_train_pso_train1$exposure))
write.csv(res_train_pso_train1, './train/res_train_pso_train.csv')
res_train_pso_train2 = res_train_pso_train1 %>%
  filter(pval < 0.05/52)
format(0.05/52,scientific = TRUE)
0.05/1399
dat_train_pso_train$samplesize.outcome = 10588 + 22806
dat_train_pso_train$ncase.outcome = 10588
dat_train_pso_train$ncontrol.outcome = 22806

steiger_pso_train = dat_train_pso_train %>%
  steiger_filtering()
dat_train_pso_train$steiger_pval = steiger_pso_train$steiger_pval
write.csv(dat_train_pso_train, '../../../train/dat_train_pso_train.csv')



dat_apof_train = dat_train_pso_train %>%
  filter(exposure == 'APOF')
dat_apof_train$effect_allele.exposure = 'C'
dat_apof_train$other_allele.exposure = 'G'
dat_apof_train$effect_allele.outcome = 'G'
dat_apof_train$other_allele.outcome = 'C'
res_apof_train <- mr_modified(dat_apof_train, prop_var_explained = T)
