mr_modified <- function(dat , prop_var_explained = T)#mr_modified：作者自己定义的函数
{
  mr_res <- mr(dat, method_list = list1)
  
  pve <- dat %>% 
    dplyr::select(id.exposure, beta.exposure, se.exposure, samplesize.exposure) %>% 
    dplyr::group_by(id.exposure) %>% 
    dplyr::summarise(pve = sum((beta.exposure^2)/(beta.exposure^2 + samplesize.exposure*se.exposure^2)))
  
  if(prop_var_explained)
  {
    mr_res <- mr_res %>% 
      dplyr::left_join(pve, by = "id.exposure")
  }
  mr_res = generate_odds_ratios(mr_res)
  return(mr_res)
}


# Table 2: Phenoscanner ---------------------------------------------------

mr_phenoscanner <- function(dat, 
                            catalog = c("GWAS", "eQTL", "pQTL", "mQTL", "methQTL"),
                            pvalue = 5e-08,
                            proxies = "EUR", 
                            r2 = 0.8,
                            build = 37)
{
  stopifnot("SNP" %in% names(dat))
  
  snpquery <- unique(dat$SNP)
  query_limit <- ceiling(0.01*length(snpquery))
  
  query_res <- list()
  for (i in 1:query_limit) {
    lo_i <- 100*i - 99
    up_i <- min(100*i,length(snpquery))
    snp_i <- snpquery[lo_i:up_i]
    print(paste0("Searching for ", snp_i))
    
    for (catalogue in catalog) {
      query_i <- try(
        phenoscanner::phenoscanner(
          snpquery = snp_i, catalogue = catalogue,
          proxies = proxies, pvalue = pvalue, r2 = r2, build = build)
      )
      if(!"try-error" %in% class(query_i))
      {
        if("results" %in% names(query_i))
        {
          query_i <- list(query_i$results)
          names(query_i) <- catalogue
          query_res <- append(query_res, query_i)
        }
      }
    }
  }
  
  return(query_res)
}


pleio_pheno <- mr_phenoscanner(table1) # the key information is the SNP column

# saveRDS(pleio_pheno, file = "phenoscanner_res.rds")
# pleio_pheno <- readRDS("phenoscanner_res.rds")

pleio_pheno_summary_raw <- pleio_pheno$GWAS %>% 
  mutate(catalog = "GWAS") %>% 
  bind_rows(pleio_pheno$pQTL%>% 
              mutate(catalog = "pQTL")) %>% 
  bind_rows(pleio_pheno$mQTL%>% 
              mutate(catalog = "mQTL")) %>% 
  filter(p != "NA",
         ancestry %in% c("Mixed", "European")) %>% 
  left_join(table1_plasma %>% 
              bind_rows(table1_csf) %>% 
              dplyr::select(exposure, Tissue:eaf.exposure) %>% 
              rename(snp = SNP),
            by = "snp") %>% 
  mutate(across(p, as.numeric)) %>% 
  group_by(trait, pmid) %>% 
  filter(p == min(p)) %>% 
  ungroup() %>% 
  dplyr::select(trait, Tissue,Protein = exposure, UniProt, snp, ref_a1, ref_a2, effect_allele.exposure:eaf.exposure, proxy, r2, rsid, a1, a2, pmid:p,catalog) %>% 
  mutate(ref_a1f = ifelse(ref_a1 == effect_allele.exposure, eaf.exposure, 1- eaf.exposure),
         ref_a1f_reliability = ifelse((ref_a1 == effect_allele.exposure&ref_a2 == other_allele.exposure)|
                                        (ref_a2 == effect_allele.exposure&ref_a1 == other_allele.exposure),
                                      "high", "low")) %>% 
  dplyr::select(trait, Tissue,Protein, UniProt, snp, effect_allele.exposure,ref_a1, ref_a2, ref_a1f,ref_a1f_reliability, proxy, r2, rsid, a1, a2, pmid:p,catalog) %>% 
  rename(effector = effect_allele.exposure,
         SNP = rsid,
         effect_allele = a1,
         other_allele = a2) %>% 
  unique() %>% 
  rename(proxy_SNP = SNP,
         SNP = snp,
         pval = p) %>% 
  mutate(effect_allele = effector,
         proxy_effect_allele = ifelse(effector == ref_a1, effect_allele, other_allele),
         beta = as.numeric(beta),
         se = as.numeric(se),
         beta = ifelse(effector == ref_a1, beta, -beta),
         proxy_SNP = ifelse(proxy==0,NA,proxy_SNP),
         r2 = ifelse(proxy==0,NA,r2),
         proxy_effect_allele = ifelse(proxy==0,NA,proxy_effect_allele)) %>% 
  dplyr::select(Tissue:SNP,trait,catalog,effect_allele,proxy_SNP,r2,proxy_effect_allele,beta,se,pval,ancestry,pmid) %>% 
  filter(!trait %in% c("Blood protein levels",
                       "Blood metabolite levels",
                       "Fc receptor-like protein 3",
                       "Alpha-2-HS-glycoprotein",
                       "CD5 antigen-like",
                       "SLAM family member 7")) %>% 
  mutate(Protein = factor(Protein, 
                          levels = c("FCRL3","TYMP","AHSG",
                                     "MMEL1","SLAMF7","CD5L"))) %>% 
  arrange(desc(Tissue), Protein)




# Table 2: Steiger filtering ----------------------------------------------
?steiger_filtering
harmonised_csf_MS %>%
  filter(SNP %in% table1_csf$SNP) %>% 
  steiger_filtering()

harmonised_plasma_MS %>%
  filter(SNP %in% table1_plasma$SNP) %>% 
  steiger_filtering()



# Table 2/ Figure S1: Bidirectional MR ---------------------------------------------

#' @source pubmed.ncbi.nlm.nih.gov/31604244/; for detail, contact linjf15
bimr_snp_MS <- read_excel("genetic_instrument_exposure_MS.xlsx")

#' @source pubmed.ncbi.nlm.nih.gov/34239129/; pubmed.ncbi.nlm.nih.gov/29875488/; pubmed.ncbi.nlm.nih.gov/34857953/
bimr_snp_protein <- read_excel("outcome_protein.xlsx")

bimr_harmonised_MS_protein <- harmonise_data(bimr_snp_MS, bimr_snp_protein)
bimr_mr_MS_protein <- mr(bimr_harmonised_MS_protein)

#' @description some functions are loaded from tutorial_main_analysis.R
bimr_table <- bimr_mr_MS_protein %>% 
  mutate(method = factor(method, 
                         levels = c("Inverse variance weighted","MR Egger","Weighted median",
                                    "Simple mode","Weighted mode")),
         outcome = factor(outcome,
                          levels = c("FCRL3","TYMP","AHSG","MMEL1","SLAMF7","CD5L"))) %>% 
  generate_odds_ratios() %>% 
  rename(Odds = or,Low = or_lci95,High = or_uci95) %>% 
  mutate(ci = sprintf("%.2f (%.2f, %.2f)",Odds,Low,High),
         pval = as.numeric(pval),
         ci = ifelse(is.na(pval),NA,ci),
         pval = ifelse(pval<0.001, scales::scientific(pval,digits = 3),sprintf("%.3f",pval))) %>% 
  dplyr::select(-exposure) %>% 
  rename(exposure = outcome) %>% 
  arrange(exposure, method) %>% 
  dplyr::select(exposure,method,ci,pval,Odds:High) %>% 
  mutate(exposure = c(NA,NA,"FCRL3",NA,NA,
                      NA,NA,"TYMP",NA,NA,
                      NA,NA,"AHSG",NA,NA,
                      NA,NA,"MMEL1",NA,NA,
                      NA,NA,"SLAMF7",NA,NA,
                      NA,NA,"CD5L",NA,NA)) %>% 
  add_title_for_table()


forestplot(labeltext = as.matrix(bimr_table[,c("exposure","method","ci","pval")]),
           mean=bimr_table$Odds,lower = bimr_table$Low,upper=bimr_table$High, 
           align = "l",
           is.summary = c(T,rep(c(F,F,F),12)),
           graph.pos = 3,
           hrzl_lines = list("2" = gpar(lty = 1,col = "black"),
                             "7" = gpar(lty = 1,col = "gray"),
                             "12" = gpar(lty = 1,col = "gray"),
                             "17" = gpar(lty = 1,col = "gray"),
                             "22" = gpar(lty = 1,col = "gray"),
                             "27" = gpar(lty = 1,col = "gray")),
           xlab = paste0(c("Lower ???  ??? Higher",paste0(rep("", 0),collapse = " ")),collapse = ""),
           zero = 1, #zero line
           graphwidth = unit(5,'cm'),
           colgap = unit(10,'mm'),
           lineheight = unit(6,'mm'), 
           col=fpColors(box='black',lines = 'black',zero = 'lightgray'),
           txt_gp=fpTxtGp(label=gpar(cex=0.85),ticks=gpar(cex=0.85),xlab=gpar(cex = 0.85)),
           xlog = F,
           xticks = c(0.8,1,1.25),
           clip = c(0.8, 1.25),
           xticks.digits = 2,
           lwd.xaxis = 1, 
           lwd.zero = 1,
           lwd.ci = 1,
           lty.ci = 1,
           ci.vertices = F,
           boxsize = 0.15,
           mar = unit(rep(0, times = 4), "mm"),
           new_page = T,
           fn.ci_norm = fpDrawNormalCI)


find_proxy <- function(snp, r2_threshold = 0.8, build = "37", pop = "EUR")
{
  stopifnot(build %in% c("37","38"))
  
  # Create and get a url
  server <- ifelse(build == "37","http://grch37.rest.ensembl.org","http://rest.ensembl.org")
  ext <- paste0("/ld/human/",snp,"/1000GENOMES:phase_3:", pop)
  url <- paste(server, ext, sep = "")
  print(paste0("searching proxies for ",snp," ......"))
  res <- httr::GET(url, httr::content_type("application/json"))
  
  # Converts http errors to R errors or warnings
  httr::stop_for_status(res)
  
  # Convert R objects from JSON
  res <- httr::content(res)
  proxy_df <- jsonlite::fromJSON(jsonlite::toJSON(res))
  proxy_df <- lapply(proxy_df, unlist)
  
  # Extract SNP associated with targeted SNP with r2 > r2_threshold
  snp_r2 <- proxy_df$r2[proxy_df$r2>r2_threshold]
  snp_v <- proxy_df$variation2[proxy_df$r2>r2_threshold]
  
  # If the query returns nothing, the output would be NULL
  # Else, returning a ordered vectors of proxy SNPs
  if(is.null(snp_v)) return(NULL)
  else return(snp_v[order(snp_r2, decreasing=T)])
}
