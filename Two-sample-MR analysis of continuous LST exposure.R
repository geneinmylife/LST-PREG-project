######### Two-sample-MR analysis of continuous LST exposure

library(devtools)
library(ggplot2)
library(TwoSampleMR)
library(readxl)
library(writexl)
library(dbplyr)
library(tidyverse)
library(plyr) #
library(data.table)
library(dplyr)
library(readr)
library(mr.raps)
library(MendelianRandomization)
setwd('/Users/Downloads/folder') #setwd to your own file folder

#exposure
exposure_dat <- read_exposure_data(
  filename = as.character("female_gcta_combined.tsv"),
  sep = "\t",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  effect_allele_col = "refA",
  #other_allele_col = "ALLELE0",
  eaf_col = "freq",
  pval_col = "p"
)
exposure_dat$exposure <- "LST_Female"

#outcome1-Ectopic Pregancy
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,  
  filename = "/Users/Downloads/Twosample/LST/-other/finngen_R11_O15_PREG_ECTOP.gz",  
  sep = "\t",  
  snp_col = "rsids",  
  beta_col = "beta",  
  se_col = "sebeta",  
  effect_allele_col = "alt",  
  other_allele_col = "ref",  
  eaf_col = "af_alt",  
  pval_col = "pval" 
)
outcome_dat$outcome <- "finngen_R11_O15_PREG_ECTOP"

#outcome2-Endometriosis
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,  
  filename = "/Users/Downloads/Twosample/LST/-other/finngen_R11_N14_ENDOMETRIOSIS.gz", 
  sep = "\t", 
  snp_col = "rsids",  
  beta_col = "beta",  
  se_col = "sebeta",  
  effect_allele_col = "alt", 
  other_allele_col = "ref",  
  eaf_col = "af_alt",  
  pval_col = "pval"  
)
outcome_dat$outcome <- "finngen_R11_N14_ENDOMETRIOSIS"

#outcome3-PCOS
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,  
  filename = "/Users/Downloads/Twosample/LST/-other/PCOS_PMC6300389.tsv", 
  sep = "\t", 
  snp_col = "SNP",  
  beta_col = "Effect	",  
  se_col = "StdErr",  
  effect_allele_col = "Effect_allele", 
  other_allele_col = "Other_allele",  
  eaf_col = "EAF",  
  pval_col = "Pvalue"  
)
outcome_dat$outcome <- "PCOS_PMC6300389"

#Conduct MR analysis
dat <- NULL 
try(dat <- harmonise_data(exposure_dat, outcome_dat,action=2))

 
mr_results <- NULL
mr_hetero <- NULL
mr_pleio <- NULL
mr_single <- NULL
try(mr_results <- mr(dat, method_list=c("mr_wald_ratio", "mr_ivw"))) 
mr_hetero <- mr_heterogeneity(dat)
mr_pleio <- mr_pleiotropy_test(dat) 
try(mr_single <- mr_singlesnp(dat))


# MR-Robust analysis
mr_robust_results <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
dat1 <- dat[ !(dat$ambiguous %in% c("TRUE")), ]
dat1 <- dat_to_MRInput(dat1)

mr_robust <- NULL
mr_robust <- MendelianRandomization::mr_ivw(dat1[[1]], robust = TRUE)

mr_robust_results <- data.frame(
  Method = "IVW",
  Estimate = mr_robust@Estimate,
  StdError = mr_robust@StdError,
  CILower = mr_robust@CILower,
  CIUpper = mr_robust@CIUpper,
  Pvalue = mr_robust@Pvalue 
)

# estimate the effects
res<-mr(dat,method_list=c("mr_ivw","mr_wald_ratio","mr_egger_regression","mr_weighted_median"))

res 
#estimate odds ratio and 95% confidence interval
exp(res$b[1])
exp(res$b[1]-1.96*res$se[1])
exp(res$b[1]+1.96*res$se[1])

# MR.raps analysis
b_exp <- dat$beta.exposure
b_out <- dat$beta.outcome
se_exp <- dat$se.exposure
se_out <- dat$se.outcome

mr_raps_results <- mr.raps(b_exp, b_out, se_exp, se_out, over.dispersion = FALSE, loss.function = "huber")

# De-biased IVW analysis
mr_debiased_ivw_results <- mr(dat, method_list = "mr_ivw_mre")

# Show MR results
mr_results
mr_robust
mr_robust_results
exp(res$b[1])
exp(res$b[1]-1.96*res$se[1])
exp(res$b[1]+1.96*res$se[1])
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
mr_debiased_ivw_results
mr_raps_results


##Save MR files
exposure <- NULL   
exposure <- "LST_female_cojo_5e8_PCOS_PMC6300389"#or other exposure-outcome names
result_file0 <- paste0("./results/",exposure,".harmonise.txt")
result_file <- paste0("./results/",exposure,".mr.txt")
result_file2 <- paste0("./results/",exposure,".mr_hetero.txt")
result_file3 <- paste0("./results/",exposure,".mr_pleio.txt")
result_file4 <- paste0("./results/",exposure,".mr_single.txt")
result_file5 <- paste0("./results/",exposure,".mr_robust_results.txt")
result_file6 <- paste0("./results/",exposure,".mr_raps_results.txt")
result_file7 <- paste0("./results/",exposure,".mr_debiased_ivw_results.txt")
if (exists("dat")==TRUE){ write.table(dat,file=result_file0,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_results")==TRUE){ write.table(mr_results,file=result_file,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_hetero")==TRUE){ write.table(mr_hetero,file=result_file2,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_pleio")==TRUE){write.table(mr_pleio,file=result_file3,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_single")==TRUE){write.table(mr_single,file=result_file4,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_robust_results")==TRUE){write.table(mr_robust_results,file=result_file5,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_raps_results")==TRUE){write.table(mr_raps_results,file=result_file6,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_debiased_ivw_results")==TRUE){write.table(mr_debiased_ivw_results,file=result_file7,sep="\t",col.names=T,row.names=F,quote=F)}
#################################################################


# 1. Create a scatter plot 
p1 <- mr_scatter_plot(res, dat)
p1 #have a look at the scatter plot
#save your plot using the png() function
png("scatter.png", width=6, height=6, units="in", res=300)
p1
dev.off()

# 2.	Create a funnel plot of the results. 
res_single <- mr_singlesnp(dat)
p2 <- mr_funnel_plot(res_single)
p2
#save your plot using the png() function
png("funnel.png", width=6, height=6, units="in", res=300)
p2
dev.off()

# 3. Create a forest plot of the results 
res_single <- mr_singlesnp(dat, all_method=c("mr_ivw", "mr_egger_regression","mr_weighted_median"))
p3 <- mr_forest_plot(res_single)
p3
#save your plot using the png() function
png("forest.png", width=6, height=6, units="in", res=300)
p3
dev.off()

# 3. Create a forest plot of the leave one out analysis results 
res_loo <- mr_leaveoneout(dat)
p4 <- mr_leaveoneout_plot(res_loo)
p4
png("loo.png", width=6, height=6, units="in", res=300)
p4
dev.off()


