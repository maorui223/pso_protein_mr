library(PheWAS)
phewas_decode = exposure_decode %>%
  filter(id.exposure %in% c(res_fa_merge_all$id.exposure, res_tl_merge_all$id.exposure, res_fi_merge$id.exposure, res_decode_long992$id.exposure, res_finn_long902$id.exposure))
phewas_decode = exposure_decode %>%
  filter(id.exposure %in% c('APOF', 'IL12B'))

write.csv(phewas_decode, 'phewas_decode2.csv')
colnames(decode_phewas_result2)[14] <- 'OR'

APOF_phewas_result = decode_phewas_result2 %>%
  filter(id.exposure == 'APOF')
IL12B_phewas_result = decode_phewas_result2 %>%
  filter(id.exposure == 'IL12B')

write.csv(APOF_phewas_result, './decode/phewas/APOF_phewas_result.csv')
write.csv(IL12B_phewas_result, './decode/phewas/IL12B_phewas_result.csv')

APOF_phewas_result = read.csv('./decode/phewas/APOF_phewas_result.csv')
APOF_phewas_result$FDR =p.adjust(APOF_phewas_result$p, method = 'fdr')
library(dplyr)
library(tidyverse)
APOF_phewas_result_filtered <- APOF_phewas_result %>%
  filter(!str_detect(`Phenotype Description`, "Psoriasis"))

phewasManhattan(APOF_phewas_result_filtered, title="APOF",y.axis.interval = 3, OR.size = T,OR.direction = T, annotate.level = 0.01)






rm(met953, met954, met956)
gc()
##################################################mediate MR
met_id = read.csv('./decode/medi/met.csv')
met_id1 = met_id$id
target_gene_exposure = exposure_decode %>%
  filter(id.exposure == 'APOF')
chd_out_dat <- extract_outcome_data(snps = target_gene_exposure$SNP, outcomes = met_id1)
dat <- harmonise_data(target_gene_exposure, chd_out_dat) 
res_expo_medi = TwoSampleMR::mr(dat, method_list = list1)
res_expo_medi = generate_odds_ratios(res_expo_medi)

res_expo_medi_met_a = res_expo_medi %>%
  filter(grepl("met-a", id.outcome))
res_expo_medi_met_c = res_expo_medi %>%
  filter(grepl("met-c", id.outcome))
res_expo_medi_met_d = res_expo_medi %>%
  filter(grepl("met-d", id.outcome))
write.csv(res_expo_medi_met_a, './decode/medi/res_expo_medi_met_a.csv')
write.csv(res_expo_medi_met_c, './decode/medi/res_expo_medi_met_c.csv')
write.csv(res_expo_medi_met_d, './decode/medi/res_expo_medi_met_d.csv')
write.csv(res_expo_medi_met_a1, './decode/medi/res_expo_medi_met_a_filter.csv')
write.csv(res_expo_medi_met_c1, './decode/medi/res_expo_medi_met_c_filter.csv')
write.csv(res_expo_medi_met_d1, './decode/medi/res_expo_medi_met_d_filter.csv')

length(unique(res_expo_medi_met_a$id.outcome))
res_expo_medi_met_a1 = res_expo_medi_met_a %>%
  filter(pval < 0.05)

length(unique(res_expo_medi_met_c$id.outcome))
res_expo_medi_met_c1 = res_expo_medi_met_c %>%
  filter(pval < 0.05)

length(unique(res_expo_medi_met_d$id.outcome))
res_expo_medi_met_d1 = res_expo_medi_met_d %>%
  filter(pval < 0.05)

medi_id = unique(c(res_expo_medi_met_a1$id.outcome, res_expo_medi_met_c1$id.outcome))

medi_snp <- extract_instruments(outcomes = medi_id)
saveRDS(medi_snp, 'medi_snp.RDS')
write.csv(medi_snp, './decode/medi/medi_snp.csv')
chd_out_dat_medi_pso <- format_data(pso_finn, 
                                    type="outcome", 
                                    phenotype_col = 'pso',
                                    snp_col = "rsids",
                                    snps = medi_snp$SNP,
                                    beta_col = "beta",
                                    se_col = "sebeta",
                                    eaf_col = "af_alt",
                                    effect_allele_col = "alt",
                                    other_allele_col = "ref",
                                    pval_col = "pval",
                                    #gene_col = 'nearest_genes',
                                    chr_col = 'chrom',
                                    pos_col = 'pos')


dat_medi_pso <- harmonise_data(medi_snp, chd_out_dat_medi_pso)
res_medi_pso = TwoSampleMR::mr(dat_medi_pso)
res_medi_pso = generate_odds_ratios(res_medi_pso)

write.csv(res_medi_pso, './decode/medi/res_medi_pso.csv')



####################################delta
med_asym <- RMediation::medci(mu.x = 0.146561290322581, mu.y = 0.0965325547780991, type = 'all',
                              se.x = 0.069658064516129, se.y = 0.0378284934846459)
med_asym[[1]]
pval=2*pnorm(-abs(0.01414794/0.005104803))
pval


install.packages("E:/UKB/MR/micro/r_0.001/Get_MR-main/EA_brain_AD_MR-main/anchors_3.0-8.tar.gz", repos = NULL, type="source")
install.packages("vcfR")
library(vcfR)
met953 <- read.vcfR("./decode/medi/MVMR/met-c-953.vcf")
meta <- data.frame(met953@meta)
fix <- data.frame(met953@fix)
gt <- data.frame(met953@gt)
head(gt)
#整理数据
fix <- data.frame(met953@fix[,c(1:5)])
fix <- fix %>% 
  dplyr::select(ID,ALT,REF,everything())
gt <- data.frame(met953@gt[,2])

beta <- as.numeric(unlist(strsplit(as.character(gt$met953.gt...2.),split = ":"))[seq(1,nrow(gt)*6,6)])
head(beta)
se <- as.numeric(unlist(strsplit(as.character(gt$met953.gt...2.),split = ":"))[seq(2,nrow(gt)*6,6)])
p <- as.numeric(unlist(strsplit(as.character(gt$met953.gt...2.),split = ":"))[seq(3,nrow(gt)*6,6)])
eaf<- as.numeric(unlist(strsplit(as.character(gt$met953.gt...2.),split = ":"))[seq(4,nrow(gt)*6,6)])
MR_data <- data.frame(beta = beta,
                      se = se,
                      adjpvalue = p,
                      eaf = eaf)
MR_data$pvalue <- 10^(MR_data$adjpvalue*-1)
MR <- cbind(fix,MR_data)
head(MR)
rownames(MR) <- NULL
colnames(MR)[1] <- "SNP"
library(data.table)
fwrite(MR, file="met953.txt")



met954 = read.vcfR("./decode/medi/MVMR/met-c-954.vcf")
meta <- data.frame(met954@meta)
fix <- data.frame(met954@fix)
gt <- data.frame(met954@gt)
head(gt)
#整理数据
fix <- data.frame(met954@fix[,c(1:5)])
fix <- fix %>% 
  dplyr::select(ID,ALT,REF,everything())
gt <- data.frame(met954@gt[,2])

beta <- as.numeric(unlist(strsplit(as.character(gt$met954.gt...2.),split = ":"))[seq(1,nrow(gt)*6,6)])
head(beta)
se <- as.numeric(unlist(strsplit(as.character(gt$met954.gt...2.),split = ":"))[seq(2,nrow(gt)*6,6)])
p <- as.numeric(unlist(strsplit(as.character(gt$met954.gt...2.),split = ":"))[seq(3,nrow(gt)*6,6)])
eaf<- as.numeric(unlist(strsplit(as.character(gt$met954.gt...2.),split = ":"))[seq(4,nrow(gt)*6,6)])
MR_data <- data.frame(beta = beta,
                      se = se,
                      adjpvalue = p,
                      eaf = eaf)
MR_data$pvalue <- 10^(MR_data$adjpvalue*-1)
MR <- cbind(fix,MR_data)
head(MR)
rownames(MR) <- NULL
colnames(MR)[1] <- "SNP"
library(data.table)
fwrite(MR, file="met954.txt")

met956 = read.vcfR("./decode/medi/MVMR/met-c-956.vcf")
meta <- data.frame(met956@meta)
fix <- data.frame(met956@fix)
gt <- data.frame(met956@gt)
head(gt)
#整理数据
fix <- data.frame(met956@fix[,c(1:5)])
fix <- fix %>% 
  dplyr::select(ID,ALT,REF,everything())
gt <- data.frame(met956@gt[,2])

beta <- as.numeric(unlist(strsplit(as.character(gt$met956.gt...2.),split = ":"))[seq(1,nrow(gt)*6,6)])
head(beta)
se <- as.numeric(unlist(strsplit(as.character(gt$met956.gt...2.),split = ":"))[seq(2,nrow(gt)*6,6)])
p <- as.numeric(unlist(strsplit(as.character(gt$met956.gt...2.),split = ":"))[seq(3,nrow(gt)*6,6)])
eaf<- as.numeric(unlist(strsplit(as.character(gt$met956.gt...2.),split = ":"))[seq(4,nrow(gt)*6,6)])
MR_data <- data.frame(beta = beta,
                      se = se,
                      adjpvalue = p,
                      eaf = eaf)
MR_data$pvalue <- 10^(MR_data$adjpvalue*-1)
MR <- cbind(fix,MR_data)
head(MR)
rownames(MR) <- NULL
colnames(MR)[1] <- "SNP"
library(data.table)
fwrite(pso_finn, "./decode/medi/MVMR/pso_finn.txt")

head(APOF)
colnames(APOF)[5] = 'REF'
colnames(APOF)[6] = 'ALT'
colnames(APOF)[4] = 'SNP'
colnames(APOF)[8] = 'pvalue'
colnames(APOF)[10] = 'se'
colnames(APOF)[12] = 'eaf'
colnames(APOF)[7] = 'beta'
colnames(APOF)[1] = 'CHROM'
colnames(APOF)[2] = 'POS'
fwrite(APOF, "./APOF_exposure.txt")

setwd('./decode/medi/MVMR/')
gc()



