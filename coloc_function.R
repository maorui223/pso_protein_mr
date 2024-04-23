

coloc_test_fa <- function(exposure_dat,
                       outcome_dat,
                       type_exposure = "quant",
                       col_snp_exposure = "rsids",
                       col_beta_exposure = "Beta",
                       col_se_exposure = "SE",
                       col_pvalues_exposure = 'Pval',
                       sd_exposure = NA,
                       col_N_exposure = "N",
                       col_MAF_exposure = "ImpMAF",
                       type_outcome = "cc",
                       col_snp_outcome = "SNP",
                       col_beta_outcome = "LOG_OR",
                       col_se_outcome = "LOG_SE",
                       col_pvalues_outcome = NA,
                       col_N_outcome = NA,
                       col_MAF_outcome = NA,
                       prevalence_outcome = NA)
{
  cols_exposure <- c(col_pvalues_exposure,
                     col_N_exposure,
                     col_MAF_exposure,
                     col_beta_exposure,
                     col_se_exposure,
                     col_snp_exposure)
  cols_exposure <- cols_exposure[!is.na(cols_exposure)]
  
  cols_outcome <- c(col_pvalues_outcome,
                    col_N_outcome,
                    col_MAF_outcome,
                    col_beta_outcome,
                    col_se_outcome,
                    col_snp_outcome)
  cols_outcome <- cols_outcome[!is.na(cols_outcome)]
  
  stopifnot(all(cols_exposure %in% names(exposure_dat)))
  stopifnot(all(cols_outcome %in% names(outcome_dat)))
  
  snp_overlap <- intersect(exposure_dat[[col_snp_exposure]],
                           outcome_dat[[col_snp_outcome]])
  snp_overlap <- unique(snp_overlap)
  
  slice_1 <- function(dat, col_snp, list_snp)
  {
    dat <- dat[dat[[col_snp]] %in% list_snp,]
    names(dat) <- stringr::str_replace(names(dat), col_snp, "SNP")
    library(tidyverse)
    dat <- dat %>% 
      group_by(SNP) %>% 
      slice(1) %>% 
      ungroup() %>% 
      arrange(SNP)
    names(dat) <- stringr::str_replace(names(dat), "SNP", col_snp)
    dat
  }
  
  exposure_dat <- slice_1(exposure_dat, col_snp_exposure, snp_overlap)
  outcome_dat <- slice_1(outcome_dat, col_snp_outcome, snp_overlap)
  
  exposure_list <- list()
  outcome_list <- list()
  
  for (i in 1:8) {
    list_element <- c("pvalues","N","MAF",
                      "beta","varbeta","snp",
                      "type","sdY")[i]
    col_element <- c(col_pvalues_exposure, col_N_exposure, col_MAF_exposure,
                     col_beta_exposure, col_se_exposure, col_snp_exposure,
                     type_exposure, sd_exposure)[i]
    if(!is.na(col_element)){
      if(list_element %in% c("type","sdY"))
      {
        if(list_element=="sdY") col_element <- as.numeric(col_element)
        exposure_list[[list_element]] <- col_element
      }
      else
      {
        if(list_element == "varbeta")
        {
          exposure_list[[list_element]] <- exposure_dat[[col_element]]*exposure_dat[[col_element]]
        }
        else
        {
          exposure_list[[list_element]] <- exposure_dat[[col_element]]
        }
      }
    }
  }
  
  for (j in 1:8) {
    list_element <- c("pvalues","N","MAF",
                      "beta","varbeta","snp",
                      "type","s")[j]
    col_element <- c(col_pvalues_outcome, col_N_outcome, col_MAF_outcome,
                     col_beta_outcome, col_se_outcome, col_snp_outcome,
                     type_outcome, prevalence_outcome)[j]
    if(!is.na(col_element)){
      if(list_element %in% c("type","s"))
      {
        if(list_element=="s") col_element <- as.numeric(col_element)
        outcome_list[[list_element]] <- col_element
      }
      else
      {
        if(list_element == "varbeta")
        {
          outcome_list[[list_element]] <- outcome_dat[[col_element]]*outcome_dat[[col_element]]
        }
        else
        {
          outcome_list[[list_element]] <- outcome_dat[[col_element]]
        }
      }
    }
  }
  
  coloc_res <- suppressWarnings(
    coloc::coloc.abf(exposure_list,outcome_list, p1 = 1e-04, p2 = 1e-04,p12 = 1e-05)
  )
  
  coloc_res$summary
}



head(TL)
coloc_test_tl <- function(exposure_dat,
                          outcome_dat,
                          type_exposure = "quant",
                          col_snp_exposure = "rsids",
                          col_beta_exposure = "Beta",
                          col_se_exposure = "SE",
                          col_pvalues_exposure = 'Pval',
                          sd_exposure = NA,
                          col_N_exposure = "N",
                          col_MAF_exposure = "ImpMAF",
                          type_outcome = "cc",
                          col_snp_outcome = "variant_id",
                          col_beta_outcome = "beta",
                          col_se_outcome = "standard_error",
                          col_pvalues_outcome = NA,
                          col_N_outcome = NA,
                          col_MAF_outcome = NA,
                          prevalence_outcome = NA)
{
  cols_exposure <- c(col_pvalues_exposure,
                     col_N_exposure,
                     col_MAF_exposure,
                     col_beta_exposure,
                     col_se_exposure,
                     col_snp_exposure)
  cols_exposure <- cols_exposure[!is.na(cols_exposure)]
  
  cols_outcome <- c(col_pvalues_outcome,
                    col_N_outcome,
                    col_MAF_outcome,
                    col_beta_outcome,
                    col_se_outcome,
                    col_snp_outcome)
  cols_outcome <- cols_outcome[!is.na(cols_outcome)]
  
  stopifnot(all(cols_exposure %in% names(exposure_dat)))
  stopifnot(all(cols_outcome %in% names(outcome_dat)))
  
  snp_overlap <- intersect(exposure_dat[[col_snp_exposure]],
                           outcome_dat[[col_snp_outcome]])
  snp_overlap <- unique(snp_overlap)
  
  slice_1 <- function(dat, col_snp, list_snp)
  {
    dat <- dat[dat[[col_snp]] %in% list_snp,]
    names(dat) <- stringr::str_replace(names(dat), col_snp, "SNP")
    library(tidyverse)
    dat <- dat %>% 
      group_by(SNP) %>% 
      slice(1) %>% 
      ungroup() %>% 
      arrange(SNP)
    names(dat) <- stringr::str_replace(names(dat), "SNP", col_snp)
    dat
  }
  
  exposure_dat <- slice_1(exposure_dat, col_snp_exposure, snp_overlap)
  outcome_dat <- slice_1(outcome_dat, col_snp_outcome, snp_overlap)
  
  exposure_list <- list()
  outcome_list <- list()
  
  for (i in 1:8) {
    list_element <- c("pvalues","N","MAF",
                      "beta","varbeta","snp",
                      "type","sdY")[i]
    col_element <- c(col_pvalues_exposure, col_N_exposure, col_MAF_exposure,
                     col_beta_exposure, col_se_exposure, col_snp_exposure,
                     type_exposure, sd_exposure)[i]
    if(!is.na(col_element)){
      if(list_element %in% c("type","sdY"))
      {
        if(list_element=="sdY") col_element <- as.numeric(col_element)
        exposure_list[[list_element]] <- col_element
      }
      else
      {
        if(list_element == "varbeta")
        {
          exposure_list[[list_element]] <- exposure_dat[[col_element]]*exposure_dat[[col_element]]
        }
        else
        {
          exposure_list[[list_element]] <- exposure_dat[[col_element]]
        }
      }
    }
  }
  
  for (j in 1:8) {
    list_element <- c("pvalues","N","MAF",
                      "beta","varbeta","snp",
                      "type","s")[j]
    col_element <- c(col_pvalues_outcome, col_N_outcome, col_MAF_outcome,
                     col_beta_outcome, col_se_outcome, col_snp_outcome,
                     type_outcome, prevalence_outcome)[j]
    if(!is.na(col_element)){
      if(list_element %in% c("type","s"))
      {
        if(list_element=="s") col_element <- as.numeric(col_element)
        outcome_list[[list_element]] <- col_element
      }
      else
      {
        if(list_element == "varbeta")
        {
          outcome_list[[list_element]] <- outcome_dat[[col_element]]*outcome_dat[[col_element]]
        }
        else
        {
          outcome_list[[list_element]] <- outcome_dat[[col_element]]
        }
      }
    }
  }
  
  coloc_res <- suppressWarnings(
    coloc::coloc.abf(exposure_list,outcome_list, p1 = 1e-04, p2 = 1e-04,p12 = 1e-05)
  )
  
  coloc_res$summary
}





head(fi)
coloc_test_fi <- function(exposure_dat,
                          outcome_dat,
                          type_exposure = "quant",
                          col_snp_exposure = "rsids",
                          col_beta_exposure = "Beta",
                          col_se_exposure = "SE",
                          col_pvalues_exposure = 'Pval',
                          sd_exposure = NA,
                          col_N_exposure = "N",
                          col_MAF_exposure = "ImpMAF",
                          type_outcome = "cc",
                          col_snp_outcome = "variant_id",
                          col_beta_outcome = "beta",
                          col_se_outcome = "standard_error",
                          col_pvalues_outcome = NA,
                          col_N_outcome = NA,
                          col_MAF_outcome = NA,
                          prevalence_outcome = NA)
{
  cols_exposure <- c(col_pvalues_exposure,
                     col_N_exposure,
                     col_MAF_exposure,
                     col_beta_exposure,
                     col_se_exposure,
                     col_snp_exposure)
  cols_exposure <- cols_exposure[!is.na(cols_exposure)]
  
  cols_outcome <- c(col_pvalues_outcome,
                    col_N_outcome,
                    col_MAF_outcome,
                    col_beta_outcome,
                    col_se_outcome,
                    col_snp_outcome)
  cols_outcome <- cols_outcome[!is.na(cols_outcome)]
  
  stopifnot(all(cols_exposure %in% names(exposure_dat)))
  stopifnot(all(cols_outcome %in% names(outcome_dat)))
  
  snp_overlap <- intersect(exposure_dat[[col_snp_exposure]],
                           outcome_dat[[col_snp_outcome]])
  snp_overlap <- unique(snp_overlap)
  
  slice_1 <- function(dat, col_snp, list_snp)
  {
    dat <- dat[dat[[col_snp]] %in% list_snp,]
    names(dat) <- stringr::str_replace(names(dat), col_snp, "SNP")
    library(tidyverse)
    dat <- dat %>% 
      group_by(SNP) %>% 
      slice(1) %>% 
      ungroup() %>% 
      arrange(SNP)
    names(dat) <- stringr::str_replace(names(dat), "SNP", col_snp)
    dat
  }
  
  exposure_dat <- slice_1(exposure_dat, col_snp_exposure, snp_overlap)
  outcome_dat <- slice_1(outcome_dat, col_snp_outcome, snp_overlap)
  
  exposure_list <- list()
  outcome_list <- list()
  
  for (i in 1:8) {
    list_element <- c("pvalues","N","MAF",
                      "beta","varbeta","snp",
                      "type","sdY")[i]
    col_element <- c(col_pvalues_exposure, col_N_exposure, col_MAF_exposure,
                     col_beta_exposure, col_se_exposure, col_snp_exposure,
                     type_exposure, sd_exposure)[i]
    if(!is.na(col_element)){
      if(list_element %in% c("type","sdY"))
      {
        if(list_element=="sdY") col_element <- as.numeric(col_element)
        exposure_list[[list_element]] <- col_element
      }
      else
      {
        if(list_element == "varbeta")
        {
          exposure_list[[list_element]] <- exposure_dat[[col_element]]*exposure_dat[[col_element]]
        }
        else
        {
          exposure_list[[list_element]] <- exposure_dat[[col_element]]
        }
      }
    }
  }
  
  for (j in 1:8) {
    list_element <- c("pvalues","N","MAF",
                      "beta","varbeta","snp",
                      "type","s")[j]
    col_element <- c(col_pvalues_outcome, col_N_outcome, col_MAF_outcome,
                     col_beta_outcome, col_se_outcome, col_snp_outcome,
                     type_outcome, prevalence_outcome)[j]
    if(!is.na(col_element)){
      if(list_element %in% c("type","s"))
      {
        if(list_element=="s") col_element <- as.numeric(col_element)
        outcome_list[[list_element]] <- col_element
      }
      else
      {
        if(list_element == "varbeta")
        {
          outcome_list[[list_element]] <- outcome_dat[[col_element]]*outcome_dat[[col_element]]
        }
        else
        {
          outcome_list[[list_element]] <- outcome_dat[[col_element]]
        }
      }
    }
  }
  
  coloc_res <- suppressWarnings(
    coloc::coloc.abf(exposure_list,outcome_list, p1 = 1e-04, p2 = 1e-04,p12 = 1e-05)
  )
  
  coloc_res$summary
}




head(CD14)
coloc_test_long <- function(exposure_dat,
                          outcome_dat,
                          type_exposure = "quant",
                          col_snp_exposure = "rsids",
                          col_beta_exposure = "Beta",
                          col_se_exposure = "SE",
                          col_pvalues_exposure = 'Pval',
                          sd_exposure = NA,
                          col_N_exposure = "N",
                          col_MAF_exposure = "ImpMAF",
                          type_outcome = "cc",
                          col_snp_outcome = "SNP",
                          col_beta_outcome = "beta.outcome",
                          col_se_outcome = "se.outcome",
                          col_pvalues_outcome = NA,
                          col_N_outcome = NA,
                          col_MAF_outcome = NA,
                          prevalence_outcome = NA)
{
  cols_exposure <- c(col_pvalues_exposure,
                     col_N_exposure,
                     col_MAF_exposure,
                     col_beta_exposure,
                     col_se_exposure,
                     col_snp_exposure)
  cols_exposure <- cols_exposure[!is.na(cols_exposure)]
  
  cols_outcome <- c(col_pvalues_outcome,
                    col_N_outcome,
                    col_MAF_outcome,
                    col_beta_outcome,
                    col_se_outcome,
                    col_snp_outcome)
  cols_outcome <- cols_outcome[!is.na(cols_outcome)]
  
  stopifnot(all(cols_exposure %in% names(exposure_dat)))
  stopifnot(all(cols_outcome %in% names(outcome_dat)))
  
  snp_overlap <- intersect(exposure_dat[[col_snp_exposure]],
                           outcome_dat[[col_snp_outcome]])
  snp_overlap <- unique(snp_overlap)
  
  slice_1 <- function(dat, col_snp, list_snp)
  {
    dat <- dat[dat[[col_snp]] %in% list_snp,]
    names(dat) <- stringr::str_replace(names(dat), col_snp, "SNP")
    library(tidyverse)
    dat <- dat %>% 
      group_by(SNP) %>% 
      slice(1) %>% 
      ungroup() %>% 
      arrange(SNP)
    names(dat) <- stringr::str_replace(names(dat), "SNP", col_snp)
    dat
  }
  
  exposure_dat <- slice_1(exposure_dat, col_snp_exposure, snp_overlap)
  outcome_dat <- slice_1(outcome_dat, col_snp_outcome, snp_overlap)
  
  exposure_list <- list()
  outcome_list <- list()
  
  for (i in 1:8) {
    list_element <- c("pvalues","N","MAF",
                      "beta","varbeta","snp",
                      "type","sdY")[i]
    col_element <- c(col_pvalues_exposure, col_N_exposure, col_MAF_exposure,
                     col_beta_exposure, col_se_exposure, col_snp_exposure,
                     type_exposure, sd_exposure)[i]
    if(!is.na(col_element)){
      if(list_element %in% c("type","sdY"))
      {
        if(list_element=="sdY") col_element <- as.numeric(col_element)
        exposure_list[[list_element]] <- col_element
      }
      else
      {
        if(list_element == "varbeta")
        {
          exposure_list[[list_element]] <- exposure_dat[[col_element]]*exposure_dat[[col_element]]
        }
        else
        {
          exposure_list[[list_element]] <- exposure_dat[[col_element]]
        }
      }
    }
  }
  
  for (j in 1:8) {
    list_element <- c("pvalues","N","MAF",
                      "beta","varbeta","snp",
                      "type","s")[j]
    col_element <- c(col_pvalues_outcome, col_N_outcome, col_MAF_outcome,
                     col_beta_outcome, col_se_outcome, col_snp_outcome,
                     type_outcome, prevalence_outcome)[j]
    if(!is.na(col_element)){
      if(list_element %in% c("type","s"))
      {
        if(list_element=="s") col_element <- as.numeric(col_element)
        outcome_list[[list_element]] <- col_element
      }
      else
      {
        if(list_element == "varbeta")
        {
          outcome_list[[list_element]] <- outcome_dat[[col_element]]*outcome_dat[[col_element]]
        }
        else
        {
          outcome_list[[list_element]] <- outcome_dat[[col_element]]
        }
      }
    }
  }
  
  coloc_res <- suppressWarnings(
    coloc::coloc.abf(exposure_list,outcome_list, p1 = 1e-04, p2 = 1e-04,p12 = 1e-05)
  )
  
  coloc_res$summary
}





head(TL)
coloc_plot_tl <- function(exposure_dat,
                       exposure,
                       col_snp_exposure = "rsids",
                       col_pvalues_exposure = "Pval",
                       outcome_dat = TL,
                       outcome = "TL",
                       col_snp_outcome = "variant_id",
                       col_pvalues_outcome = "p_value",
                       path = tempdir(),
                       combine = FALSE,
                       legend = FALSE)
{
  library(locuscomparer)
  
  snp_overlap <- intersect(exposure_dat[[col_snp_exposure]],
                           outcome_dat[[col_snp_outcome]])
  snp_overlap <- unique(snp_overlap)
  
  exposure_dat <- exposure_dat[exposure_dat[[col_snp_exposure]] %in% snp_overlap,]
  outcome_dat <- outcome_dat[outcome_dat[[col_snp_outcome]] %in% snp_overlap,]
  
  exposure_dat <- exposure_dat[order(exposure_dat[[col_snp_exposure]]),]
  outcome_dat <- outcome_dat[order(outcome_dat[[col_snp_outcome]]),]
  
  exposure_dat <- data.frame(
    rsid = exposure_dat[[col_snp_exposure]],
    pval = exposure_dat[[col_pvalues_exposure]])
  
  outcome_dat <- data.frame(
    rsid = outcome_dat[[col_snp_outcome]],
    pval = outcome_dat[[col_pvalues_outcome]])
  
  write.table(exposure_dat,paste0(path,"/exposure_test.tsv"),sep = "\t",row.names = F,quote = F)
  write.table(outcome_dat,paste0(path,"/outcome_test.tsv"),sep = "\t",row.names = F,quote = F)
  p <- locuscompare(in_fn1 = paste0(path,"/outcome_test.tsv"), 
                    in_fn2 = paste0(path,"/exposure_test.tsv"), 
                    title1 = paste0(outcome," GWAS"),
                    title2 = paste0(exposure," pQTL"),
                    combine = combine,
                    legend = legend)
  p
}

head(fi)
coloc_plot_fi <- function(exposure_dat,
                          exposure,
                          col_snp_exposure = "rsids",
                          col_pvalues_exposure = "Pval",
                          outcome_dat = fi,
                          outcome = "fi",
                          col_snp_outcome = "variant_id",
                          col_pvalues_outcome = "p_value",
                          path = tempdir(),
                          combine = FALSE,
                          legend = FALSE)
{
  library(locuscomparer)
  
  snp_overlap <- intersect(exposure_dat[[col_snp_exposure]],
                           outcome_dat[[col_snp_outcome]])
  snp_overlap <- unique(snp_overlap)
  
  exposure_dat <- exposure_dat[exposure_dat[[col_snp_exposure]] %in% snp_overlap,]
  outcome_dat <- outcome_dat[outcome_dat[[col_snp_outcome]] %in% snp_overlap,]
  
  exposure_dat <- exposure_dat[order(exposure_dat[[col_snp_exposure]]),]
  outcome_dat <- outcome_dat[order(outcome_dat[[col_snp_outcome]]),]
  
  exposure_dat <- data.frame(
    rsid = exposure_dat[[col_snp_exposure]],
    pval = exposure_dat[[col_pvalues_exposure]])
  
  outcome_dat <- data.frame(
    rsid = outcome_dat[[col_snp_outcome]],
    pval = outcome_dat[[col_pvalues_outcome]])
  
  write.table(exposure_dat,paste0(path,"/exposure_test.tsv"),sep = "\t",row.names = F,quote = F)
  write.table(outcome_dat,paste0(path,"/outcome_test.tsv"),sep = "\t",row.names = F,quote = F)
  p <- locuscompare(in_fn1 = paste0(path,"/outcome_test.tsv"), 
                    in_fn2 = paste0(path,"/exposure_test.tsv"), 
                    title1 = paste0(outcome," GWAS"),
                    title2 = paste0(exposure," pQTL"),
                    combine = combine,
                    legend = legend)
  p
}


head(fa)
coloc_plot_fa <- function(exposure_dat,
                          exposure,
                          col_snp_exposure = "rsids",
                          col_pvalues_exposure = "Pval",
                          outcome_dat = fa,
                          outcome = "fa",
                          col_snp_outcome = "SNP",
                          col_pvalues_outcome = "P_BOLT_LMM_INF",
                          path = tempdir(),
                          combine = FALSE,
                          legend = FALSE)
{
  library(locuscomparer)
  
  snp_overlap <- intersect(exposure_dat[[col_snp_exposure]],
                           outcome_dat[[col_snp_outcome]])
  snp_overlap <- unique(snp_overlap)
  
  exposure_dat <- exposure_dat[exposure_dat[[col_snp_exposure]] %in% snp_overlap,]
  outcome_dat <- outcome_dat[outcome_dat[[col_snp_outcome]] %in% snp_overlap,]
  
  exposure_dat <- exposure_dat[order(exposure_dat[[col_snp_exposure]]),]
  outcome_dat <- outcome_dat[order(outcome_dat[[col_snp_outcome]]),]
  
  exposure_dat <- data.frame(
    rsid = exposure_dat[[col_snp_exposure]],
    pval = exposure_dat[[col_pvalues_exposure]])
  
  outcome_dat <- data.frame(
    rsid = outcome_dat[[col_snp_outcome]],
    pval = outcome_dat[[col_pvalues_outcome]])
  
  write.table(exposure_dat,paste0(path,"/exposure_test.tsv"),sep = "\t",row.names = F,quote = F)
  write.table(outcome_dat,paste0(path,"/outcome_test.tsv"),sep = "\t",row.names = F,quote = F)
  p <- locuscompare(in_fn1 = paste0(path,"/outcome_test.tsv"), 
                    in_fn2 = paste0(path,"/exposure_test.tsv"), 
                    title1 = paste0(outcome," GWAS"),
                    title2 = paste0(exposure," pQTL"),
                    combine = combine,
                    legend = legend)
  p
}

head(long)
coloc_plot_long99 <- function(exposure_dat,
                          exposure,
                          col_snp_exposure = "rsids",
                          col_pvalues_exposure = "Pval",
                          outcome_dat = long99_hg19,
                          outcome = "long99",
                          col_snp_outcome = "SNP",
                          col_pvalues_outcome = "P-value",
                          path = tempdir(),
                          combine = FALSE,
                          legend = FALSE)
{
  library(locuscomparer)
  
  snp_overlap <- intersect(exposure_dat[[col_snp_exposure]],
                           outcome_dat[[col_snp_outcome]])
  snp_overlap <- unique(snp_overlap)
  
  exposure_dat <- exposure_dat[exposure_dat[[col_snp_exposure]] %in% snp_overlap,]
  outcome_dat <- outcome_dat[outcome_dat[[col_snp_outcome]] %in% snp_overlap,]
  
  exposure_dat <- exposure_dat[order(exposure_dat[[col_snp_exposure]]),]
  outcome_dat <- outcome_dat[order(outcome_dat[[col_snp_outcome]]),]
  
  exposure_dat <- data.frame(
    rsid = exposure_dat[[col_snp_exposure]],
    pval = exposure_dat[[col_pvalues_exposure]])
  
  outcome_dat <- data.frame(
    rsid = outcome_dat[[col_snp_outcome]],
    pval = outcome_dat[[col_pvalues_outcome]])
  
  write.table(exposure_dat,paste0(path,"/exposure_test.tsv"),sep = "\t",row.names = F,quote = F)
  write.table(outcome_dat,paste0(path,"/outcome_test.tsv"),sep = "\t",row.names = F,quote = F)
  p <- locuscompare(in_fn1 = paste0(path,"/outcome_test.tsv"), 
                    in_fn2 = paste0(path,"/exposure_test.tsv"), 
                    title1 = paste0(outcome," GWAS"),
                    title2 = paste0(exposure," pQTL"),
                    combine = combine,
                    legend = legend)
  p
}



gridExtra::grid.arrange(
  coloc_plot(exposure_dat = coloc_FCRL3_dat, exposure = "FCRL3",combine = F)$locuscompare+ggtitle("(A)"),
  coloc_plot(exposure_dat = coloc_AHSG_dat, exposure = "AHSG",combine = F)$locuscompare+ggtitle("(B)"),
  coloc_plot(exposure_dat = coloc_TYMP_dat, exposure = "TYMP",combine = F)$locuscompare+ggtitle("(C)"),
  coloc_plot(exposure_dat = coloc_MMEL1_dat, exposure = "MMEL1",combine = F)$locuscompare+ggtitle("(D)"),
  coloc_plot(exposure_dat = coloc_SLAMF7_dat, exposure = "SLAMF7",combine = F)$locuscompare+ggtitle("(E)"),
  coloc_plot(exposure_dat = coloc_CD5L_dat, exposure = "CD5L",combine = F)$locuscompare+ggtitle("(F)"),
  nrow = 2, ncol = 3)