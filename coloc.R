###############coloc#################
setwd('./decode/coloc/')
rm(HCC_finn,HCC_UKB)
gc()
NCAN = fread('./coloc/15573_110_NCAN_CSPG3.txt.gz')
CILP2 = fread('./coloc/8841_65_CILP2_CILP2.txt.gz')
PSG5 = fread('./coloc/9314_9_PSG5_PSG5.txt.gz')
PSG5 = fread('./coloc/GCST90090624_buildGRCh37.tsv')

head(NCAN)
head(CILP2)
str(CILP2)
NCAN <- NCAN[Chrom == 'chr19' &
               Pos > 19211958  - 100*1000 &
               Pos < 19211958  + 100*1000]
CILP2 <- CILP2[Chrom == 'chr19' &
                 Pos > 19538265  - 100*1000 &
                 Pos < 19538265  + 100*1000]
PSG5 <- PSG5[chromosome == '19' &
               base_pair_location > 43167743  - 100*1000 &
               base_pair_location < 43167743  + 100*1000]



NCAN$pvalue = 10^(-NCAN$LOG10P)
chr19 = fread('./coloc/olink_rsid_map_mac5_info03_b0_7_chr19_patched_v2.tsv')
head(chr19)
NCAN = merge(NCAN, chr19, by = 'ID')
head(NCAN)
NCAN = drop_na(NCAN, ImpMAF)
NCAN = drop_na(NCAN, Pval)
NCAN = as.data.frame(NCAN)
write.csv(NCAN, './coloc/NCAN_coloc.csv')

NCAN = drop_na(NCAN, Pval)
NCAN = drop_na(NCAN, ImpMAF)
NCAN = as.data.frame(NCAN)
NCAN$N = 5363

CILP2 = drop_na(CILP2, p_value)
CILP2 = drop_na(CILP2, effect_allele_frequency)
CILP2 = as.data.frame(CILP2)
CILP2$N = 5363

PSG5 = drop_na(PSG5, p_value)
PSG5 = drop_na(PSG5, effect_allele_frequency)
PSG5 = as.data.frame(PSG5)

head(HCC_UKB)
HCC_UKB = drop_na(HCC_UKB, p_value)
HCC_UKB = drop_na(HCC_UKB, effect_allele_frequency)
HCC_UKB = as.data.frame(HCC_UKB)


NCAN_coloc = coloc_test_hcc(exposure_dat = NCAN, outcome_dat = HCC_UKB)
CILP2_coloc = coloc_test_hcc(exposure_dat = CILP2, outcome_dat = HCC_UKB)
PSG5_coloc = coloc_test_hcc(exposure_dat = PSG5, outcome_dat = HCC_UKB)


options(ggrepel.max.overlaps = Inf)
coloc_plot_hcc(exposure_dat = NCAN, exposure = "NCAN",combine = T)



coloc_result = rbind(NCAN_coloc, CILP2_coloc)
coloc_result = as.data.frame(coloc_result)
write.csv(coloc_result, './decode/coloc/coloc_result_all.csv')


################SMR NCAN##################
head(NCAN_smr)
NCAN_smr$Chr = gsub('chr', '',NCAN_smr$Chr)
NCAN_smr =NCAN[Chrom == 'chr19' &
                 Pos > 19211958  - 1000*1000 &
                 Pos < 19211958  + 1000*1000 &
                 Pval < 1.57e-3]
NCAN_smr = NCAN_smr %>%
  select(4,1,2,5,6,12,7,10,8)

NCAN_smr$Probe = 'ENSG00000130287'
NCAN_smr$Probe_Chr = 19
NCAN_smr$Probe_bp = 19211958
NCAN_smr$Gene = 'NCAN'
NCAN_smr$Orientation = 'N'
NCAN_smr$Probe_Chr = as.character(NCAN_smr$Probe_Chr)


NCAN_smr<- NCAN_smr[!grepl("[a-zA-Z].*[a-zA-Z]", NCAN_smr$effectAllele), ]
NCAN_smr<- NCAN_smr[!grepl("[a-zA-Z].*[a-zA-Z]", NCAN_smr$otherAllele), ]
head(NCAN_smr)

NCAN_smr = as.data.frame(NCAN_smr)
NCAN_smr <- NCAN_smr[!duplicated(NCAN_smr[, 1:3]), ]


NCAN_smr <-NCAN_smr %>%  drop_na()
str(NCAN_smr)
NCAN_smr <-NCAN_smr %>% 
  select(1:6,10:14,7:9)

colnames(NCAN_smr)[1:6] = c("SNP","Chr","BP","A1","A2","Freq")
colnames(NCAN_smr)[12:14] = c("b","se","p")
NCAN_smr$SNP = sub(",.*", "", NCAN_smr$SNP)
write.table(NCAN_smr,file = "./SMR/NCAN_smr.txt", quote =FALSE,row.names = F)

####################HCC_finn#################
head(HCC_finn)
HCC_finn$N = 314193 + 500
hcc_smr = HCC_finn %>%
  select(5,4,3,11,9,10,7,14)
hcc_smr <- format_data(hcc_smr, 
                       type="exposure",
                       snp_col = "rsids",
                       beta_col = "beta",
                       se_col = "sebeta",
                       eaf_col = "af_alt",
                       effect_allele_col = "alt",
                       other_allele_col = "ref",
                       pval_col = "pval",
                       samplesize_col = 'N')

hcc_smr = hcc_smr %>%
  select(1:8)
colnames(hcc_smr) = c("SNP","A1","A2","freq","b","se","p","n")
str(hcc_smr)
hcc_smr<- hcc_smr[!grepl("[a-zA-Z].*[a-zA-Z]", hcc_smr$A1), ]
hcc_smr<- hcc_smr[!grepl("[a-zA-Z].*[a-zA-Z]", hcc_smr$A2), ]
# 假设data是您的数据框

# 仅基于某列（例如第一列）删除空白行
hcc_smr <- hcc_smr[hcc_smr[[1]] != "", ]
hcc_smr$SNP = sub(",.*", "", hcc_smr$SNP)

fwrite(hcc_smr,file = "./SMR/hcc_smr.txt", quote =FALSE,row.names = F)


################SMR CILP2##################
head(CILP2)
CILP2_smr =CILP2[Chrom == 'chr19' &
                   Pos > 19538265  - 1000*1000 &
                   Pos < 19538265  + 1000*1000 &
                 Pval < 1.57e-3]
CILP2_smr = CILP2_smr %>%
  select(4,1,2,5,6,12,7,10,8)

CILP2_smr$Probe = 'ENSG00000160161'
CILP2_smr$Probe_Chr = 19
CILP2_smr$Probe_bp = 19538265
CILP2_smr$Gene = 'CILP2'
CILP2_smr$Orientation = 'N'
CILP2_smr$Probe_Chr = as.character(CILP2_smr$Probe_Chr)


CILP2_smr<- CILP2_smr[!grepl("[a-zA-Z].*[a-zA-Z]", CILP2_smr$effectAllele), ]
CILP2_smr<- CILP2_smr[!grepl("[a-zA-Z].*[a-zA-Z]", CILP2_smr$otherAllele), ]
head(CILP2_smr)

CILP2_smr = as.data.frame(CILP2_smr)
CILP2_smr <- CILP2_smr[!duplicated(CILP2_smr[, 1:3]), ]


CILP2_smr <-CILP2_smr %>%  drop_na()
str(CILP2_smr)
CILP2_smr <-CILP2_smr %>% 
  select(1:6,10:14,7:9)

colnames(CILP2_smr)[1:6] = c("SNP","Chr","BP","A1","A2","Freq")
colnames(CILP2_smr)[12:14] = c("b","se","p")
CILP2_smr$SNP = sub(",.*", "", CILP2_smr$SNP)
CILP2_smr$Chr = gsub('chr', '',CILP2_smr$Chr)
write.table(CILP2_smr,file = "./SMR/CILP2_smr.txt", quote =FALSE,row.names = F)

################SMR PSG5##################
head(PSG5)
PSG5_smr =PSG5[Chrom == 'chr19' &
                 Pos > 43167743  - 1000*1000 &
                 Pos < 43167743  + 1000*1000 &
                   Pval < 1.57e-3]
PSG5_smr = PSG5_smr %>%
  select(4,1,2,5,6,12,7,10,8)

PSG5_smr$Probe = 'ENSG00000204941'
PSG5_smr$Probe_Chr = 19
PSG5_smr$Probe_bp = 43167743
PSG5_smr$Gene = 'PSG5'
PSG5_smr$Orientation = 'N'
PSG5_smr$Probe_Chr = as.character(PSG5_smr$Probe_Chr)


PSG5_smr<- PSG5_smr[!grepl("[a-zA-Z].*[a-zA-Z]", PSG5_smr$effectAllele), ]
PSG5_smr<- PSG5_smr[!grepl("[a-zA-Z].*[a-zA-Z]", PSG5_smr$otherAllele), ]
head(PSG5_smr)

PSG5_smr = as.data.frame(PSG5_smr)
PSG5_smr <- PSG5_smr[!duplicated(PSG5_smr[, 1:3]), ]


PSG5_smr <-PSG5_smr %>%  drop_na()
str(PSG5_smr)
PSG5_smr <-PSG5_smr %>% 
  select(1:6,10:14,7:9)

colnames(PSG5_smr)[1:6] = c("SNP","Chr","BP","A1","A2","Freq")
colnames(PSG5_smr)[12:14] = c("b","se","p")
PSG5_smr$SNP = sub(",.*", "", PSG5_smr$SNP)
PSG5_smr$Chr = gsub('chr', '',PSG5_smr$Chr)
write.table(PSG5_smr,file = "./SMR/PSG5_smr.txt", quote =FALSE,row.names = F)
