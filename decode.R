exposure_decode = read.csv('decode/pQpso_finn_deCODE_clumped.csv')
head(exposure_decode)
########################################################pso_finn########
pso_finn <- fread('./finngen_R9_L12_PSORIASIS.gz')
rm(pso_finn,pso_train)
gc()
head(pso_finn)
chd_out_dat_decode_pso_finn <- format_data(pso_finn, 
                                           type="outcome", 
                                           phenotype_col = 'pso',
                                           snp_col = "rsids",
                                           snps = exposure_decode$SNP,
                                           beta_col = "beta",
                                           se_col = "sebeta",
                                           eaf_col = "af_alt",
                                           effect_allele_col = "alt",
                                           other_allele_col = "ref",
                                           pval_col = "pval",
                                           #gene_col = 'nearest_genes',
                                           chr_col = 'chrom',
                                           pos_col = 'pos')
head(chd_out_dat_decode_pso_finn)
dat_decode_pso_finn <- harmonise_data(exposure_decode, chd_out_dat_decode_pso_finn)
dat_decode_pso_finn = dat_decode_pso_finn %>%
  mutate(R2 = 2*(1-eaf.exposure)*eaf.exposure*(beta.exposure^2))
dat_decode_pso_finn = dat_decode_pso_finn %>%
  distinct(exposure,SNP,beta.exposure,eaf.exposure, .keep_all = T)
heterogeneity_decode_pso_finn = mr_heterogeneity(dat_decode_pso_finn)
pleiotropy_decode_pso_finn = mr_pleiotropy_test(dat_decode_pso_finn)
write.csv(heterogeneity_decode_pso_finn, './decode/heterogeneity_decode_pso_finn.csv')
write.csv(pleiotropy_decode_pso_finn, './decode/pleiotropy_decode_pso_finn.csv')

list = mr_method_list()
list1 = list[c(1,8),]$obj
dat_decode_pso_finn$samplesize.exposure = 35559
res_decode_pso_finn <- mr_modified(dat_decode_pso_finn, prop_var_explained = T)
res_decode_pso_finn1 = res_decode_pso_finn %>%
  filter(pval > 0)
length(unique(res_decode_pso_finn1$exposure))
setwd('../..')
write.csv(res_decode_pso_finn1, './decode/res_decode_pso_finn.csv')
res_decode_pso_finn2 = res_decode_pso_finn1 %>%
  filter(pval < 0.05/1410)
res_pso_finn_merge = merge(res_finn_pso_finn2, res_decode_pso_finn2,by.x = 'exposure', by.y  = 'id.exposure')
dat_decode_pso_finn$samplesize.outcome = 364071+9267
dat_decode_pso_finn$ncase.outcome = 9267
dat_decode_pso_finn$ncontrol.outcome = 364071
dat_decode_pso_finn$exposure = dat_decode_pso_finn$id.exposure
steiger_pso_finn = dat_decode_pso_finn %>%
  steiger_filtering()
dat_decode_pso_finn$steiger_pval = steiger_pso_finn$steiger_pval
write.csv(dat_decode_pso_finn, '../../dat_decode_pso_finn.csv')

0.05/1410
9267/(364071+9267)

###################meta#########
library(meta)

data1 = read.csv('./decode/meta.csv')
lnor<- log(data1[,"or"])
lnuci<- log(data1[,"or_uci95"])
lnlci<- log(data1[,"or_lci95"])
selnor<- (lnuci-lnlci)/(2*1.96)

?metagen
pfs=metagen(lnor,selnor,sm="OR",data=data1, studlab = data1$exposure_to_outcome,subset = c(Protein == 'IL12B'),  comb.fixed = FALSE)

settings.meta('revman5')
forest(pfs)
# chd_out_dat_decode_pso_train <- format_data(
#   pso_train,
#   type = "outcome",
#   phenotype_col = 'pso',
#   snp_col = "rsid",
#   snps = exposure_decode$SNP,
#   beta_col = "beta",
#   se_col = "se",
#   #eaf_col = "af_alt",
#   effect_allele_col = "effect_allele",
#   other_allele_col = "other_allele",
#   pval_col = "p",
#   #gene_col = 'nearest_genes',
#   chr_col = 'chrom',
#   pos_col = 'pos'
# )
# chd_out_dat_decode_pso_train = get_eaf_from_1000G(chd_out_dat_decode_pso_train,
#                                                   path = 'E:/UKB/MR/micro/r_0.001/drug/1kg.v3/EUR/',
#                                                   type = "outcome")
# head(chd_out_dat_decode_pso_train)
# dat_decode_pso_train <-
#   harmonise_data(exposure_decode, chd_out_dat_decode_pso_train)
# dat_decode_pso_train = dat_decode_pso_train %>%
#   mutate(R2 = 2 * (1 - eaf.exposure) * eaf.exposure * (beta.exposure ^ 2))
# dat_decode_pso_train = dat_decode_pso_train %>%
#   distinct(exposure, SNP, beta.exposure, eaf.exposure, .keep_all = T)
# heterogeneity_decode_pso_train = mr_heterogeneity(dat_decode_pso_train)
# pleiotropy_decode_pso_train = mr_pleiotropy_test(dat_decode_pso_train)
# write.csv(heterogeneity_decode_pso_train,
#           './decode/heterogeneity_decode_pso_train.csv')
# write.csv(pleiotropy_decode_pso_train,
#           './decode/pleiotropy_decode_pso_train.csv')
# 
# list = mr_method_list()
# list1 = list[c(1, 8), ]$obj
# dat_decode_pso_train$samplesize.exposure = 35559
# res_decode_pso_train <-
#   mr_modified(dat_decode_pso_train, prop_var_explained = T)
# res_decode_pso_train1 = res_decode_pso_train %>%
#   filter(pval > 0)
# length(unique(res_decode_pso_train1$exposure))
# setwd('../..')
# write.csv(res_decode_pso_train1, './decode/res_decode_pso_train.csv')
# res_decode_pso_train2 = res_decode_pso_train1 %>%
#   filter(pval < 0.05 / 121)
# res_pso_train_merge = merge(res_train_pso_train2,
#                             res_decode_pso_train2,
#                             by.x = 'exposure',
#                             by.y  = 'id.exposure')
# dat_decode_pso_train$samplesize.outcome = 364071 + 9267
# steiger_pso_train = dat_decode_pso_train %>%
#   filter(id.exposure %in% res_pso_train_merge$id.exposure) %>%
#   steiger_filtering()
# res_pso_train_merge$steiger_pval = steiger_pso_train$steiger_pval
# write.csv(res_pso_train_merge, './merge/res_pso_train_merge.csv')

################################################## volcano_plot
library(ggvenn)

volcano_plot <- function(.data,
                         number_comparasion = 1,
                         title = "(A)",
                         col_beta = "b",
                         col_size = "pve",
                         col_label = "exposure",
                         legend.position = "none",
                         colors = c(negative = "#008B8B", positive = "#B80000"))
{
  p_thershold <- 0.05 / number_comparasion
  
  message("The function is written by Jianfeng Lin")
  .data <- .data %>%
    mutate(regulate = ifelse(b > 0, "positive", "negative"))
  
  p <- .data %>%
    rename(beta := !!col_beta,
           size := !!col_size,
           label := !!col_label) %>%
    mutate(
      x = beta,
      y = -log10(pval),
      label = ifelse(pval < p_thershold, label, NA)
    ) %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(size = size, color = regulate), alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = -log10(p_thershold),
               linetype = 2) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6.5),
      legend.position = legend.position
    ) +
    labs(
      x = "ln(OR)",
      y = parse(text = "-log[10]*(italic(P)-value)"),
      title = title
    ) +
    scale_size(name = "PVE",
               breaks = c(0.2 * 1:3)) +
    ggrepel::geom_label_repel(aes(label = label), size = 3) +
    scale_color_manual(values = colors) # 使用自定义颜色
  
  plot(p)
}



res_decode_pso_finn1$exposure = res_decode_pso_finn1$id.exposure


list_of_pso_finn <- list(
  Decode = res_decode_pso_finn2$id.exposure,
  FinnGen = res_finn_pso_ukb2$id.exposure,
  validation = res_train_pso_train2$id.exposure
)


library(venn)
library(VennDiagram)
library(eulerr)

? ggvenn
create_venn <- function(list_of_sets, title, colors) {
  p <- ggvenn(list_of_sets,
              show_percentage = F,
              show_elements = T) +
    scale_fill_manual(values = colors) # 使用自定义颜色
  
  p +
    labs(title = title) +
    theme_minimal()
}




gridExtra::grid.arrange(
  #开始用数据画图
  volcano_plot(
    res_decode_pso_finn1,
    number_comparasion = 1410,
    colors = c(negative = "#66A61E", positive = "#E7298A"),
    legend.position = c(0.9, 0.9),
    title = "A  Decode FinnGen pso"
  ),
  volcano_plot(
    res_finn_pso_ukb1,
    colors = c(negative = "#008B8B", positive = "#B80000"),
    legend.position = c(0.9, 0.9),
    number_comparasion = 1320,
    title = "B  Finn UKB pso"
  ),
  volcano_plot(
    res_train_pso_train1,
    number_comparasion = 52,
    colors = c(negative = "#1B9E77", positive = "#FF7F00"),
    legend.position = c(0.9, 0.9),
    title = "C  validation pso"
  ),
  create_venn(
    list_of_sets = list_of_pso_finn,
    title = "D mr results of pso",
    colors = c("#E7298A", "#B80000", "#FF7F00")
  ),
  ncol = 2,
  nrow = 2
)



