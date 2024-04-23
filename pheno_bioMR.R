snp_decode = subset(exposure_decode, id.exposure == 'APOF')
snp_decode = snp_decode[,2]
snp_finn = subset(exposure_finn, exposure == 'APOF')
snp_finn =snp_finn[,2]
snp_train = subset(exposure_train, exposure == 'APOF')
snp_train = snp_train[,4]
snp_fa = rbind(snp_decode, snp_finn, snp_train)
library(phenoscanner)
snp_fa = as.data.frame(snp_fa)
colnames(snp_fa)[1] = 'SNP'
fa_pheno <- mr_phenoscanner(snp_fa)
fa_pheno_summary_raw <- fa_pheno$GWAS %>% 
mutate(catalog = "GWAS") %>% 
  bind_rows(fa_pheno$pQTL%>% 
              mutate(catalog = "pQTL")) %>% 
  bind_rows(fa_pheno$mQTL%>% 
              mutate(catalog = "mQTL")) %>% 
  filter(p != "NA",
         ancestry %in% c("Mixed", "European")) %>% 
  mutate(across(p, as.numeric)) %>% 
  group_by(trait, pmid) %>% 
  filter(p == min(p))
write.csv(fa_pheno_summary_raw, '../../pheno_summary.csv')



##################################################bio-MR
head(pso_finn)
pso_bio = pso_finn %>%
  filter(pval < 5e-8)
pso_bio = format_data(pso_bio, 
                      type="exposure",
                      snp_col = "rsids",
                      beta_col = "beta",
                      se_col = "sebeta",
                      eaf_col = "af_alt",
                      effect_allele_col = "alt",
                      other_allele_col = "ref",
                      pval_col = "pval",
                      #gene_col = 'nearest_genes',
                      chr_col = '#chrom',
                      pos_col = 'pos')
pso_bio = clean_expo(pso_bio, pval = 5e-8, clump = TRUE, kb = 10000, r2 = 0.001, LD_file='E:/UKB/MR/micro/r_0.001/drug/1kg.v3/EUR/EUR')

head(APOF)
APOF_bio<- format_data(APOF, 
               type="outcome",
               snps = pso_bio$SNP,
               snp_col = "rsids",
               beta_col = "Beta",
               se_col = "SE",
               eaf_col = "ImpMAF",
               effect_allele_col = "effectAllele",
               other_allele_col = "otherAllele",
               pval_col = "Pval",
               #gene_col = 'nearest_genes',
               chr_col = 'Chrom',
               pos_col = 'Pos')

dat_bio <- harmonise_data(pso_bio, APOF_bio)
heterogeneity_bio = mr_heterogeneity(dat_bio)
pleiotropy_bio = mr_pleiotropy_test(dat_bio)

bio_mr_results = mr(dat_bio)
bio_mr_results = generate_odds_ratios(bio_mr_results)
write.csv(bio_mr_results, '../../bio_mr_results.csv')
