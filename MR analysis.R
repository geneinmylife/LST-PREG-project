## Content Checklist ##
# 1. Two-sample-MR analysis
# 2. MVMR analysis
# 3. MR-PRESSO analysis
# 4. NLMR analysis


#########################################
####### 1. Two-sample-MR analysis #######

library(TwoSampleMR)
library(mr.raps)
library(MendelianRandomization)

exposure_dat <- read_exposure_data(
  filename = as.character("exposure.tsv"),
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "FRQ",
  pval_col = "P"
)
exposure_dat$exposure <- "LST"

outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,  
  filename = "outcome.tsv.gz", 
  sep = "\t",  
  snp_col = "rs_id",  
  beta_col = "beta",  
  se_col = "standard_error", 
  effect_allele_col = "effect_allele", 
  other_allele_col = "other_allele",  
  eaf_col = "effect_allele_frequency", 
  pval_col = "p_value" 
)
outcome_dat$outcome <- "ENDOM"

dat <- NULL 
try(dat <- harmonise_data(exposure_dat, outcome_dat,action=2))

mr_results <- NULL
mr_hetero <- NULL
mr_pleio <- NULL
mr_single <- NULL
try(mr_results <- mr(dat, method_list=c("mr_wald_ratio", "mr_ivw","mr_egger_regression","mr_weighted_median"))) 
mr_hetero <- mr_heterogeneity(dat)
mr_pleio <- mr_pleiotropy_test(dat) 
try(mr_single <- mr_singlesnp(dat))

mr_raps <- mr.raps(dat, over.dispersion = TRUE, loss.function = "huber")

dat1 <- dat[ !(dat$ambiguous %in% c("TRUE")), ]
dat1 <- dat_to_MRInput(dat1)
mr_robust <- MendelianRandomization::mr_ivw(dat1[[1]], robust = TRUE)

ham_dat_2 <- dat_to_MRInput(dat1)
MR_DIVW<-mr_divw(ham_dat_2[[1]])
MR_DIVW_results<-cbind('LST',MR_DIVW$SNPs,MR_DIVW$Estimate,MR_DIVW$StdError,MR_DIVW$CILower,MR_DIVW$CIUpper,MR_DIVW$Pvalue)


dat$rsq.exposure <- (get_r_from_pn(p=dat$pval.exposure, 
                                   n=dat$samplesize.exposure))^2
dat$rsq.outcome <- (get_r_from_lor(
  lor=dat$beta.outcome,af=dat$eaf.outcome,
  ncase=dat$ncase,ncontrol=dat$ncontrol,
  prevalence=dat$ncase/dat$samplesize.outcome))^2
harmdat <- dat
st <- psych::r.test( 
  n = harmdat$samplesize.exposure, 
  n2 = harmdat$samplesize.outcome, 
  r12 = sqrt(harmdat$rsq.exposure), 
  r34 = sqrt(harmdat$rsq.outcome))
harmdat$steiger_dir <- harmdat$rsq.exposure > harmdat$rsq.outcome
harmdat$steiger_pval <- pnorm(-abs(st$z)) * 2



####################################
######### 2. MVMR analysis #########
library(TwoSampleMR)
library(MendelianRandomization)
library(MVMR)

mvmr_input_mr <- mr_mvinput(
  bx = cbind( mvmr_data$beta_LST_Female, mvmr_data$beta_BMI]), 
  by = mvmr_data$beta.outcome,  
  bxse  = cbind( mvmr_data$se_LST_Female, mvmr_data$se_BMI),
  byse = mvmr_data$se.outcome,  
  snps = mvmr_data$SNP  
)

MRMVObject <- mr_mvivw(mvmr_input_mr, 
                       model = "default",
                       correl = FALSE,
                       distribution = "normal",
                       alpha = 0.05)

MRMVObject_egger<-mr_mvegger(mvmr_input_mr,
                             orientate = 1,
                             correl = FALSE,
                             distribution = "normal",
                             alpha = 0.05)


sres <- strength_mvmr(r_input = mvmr_input_mr, gencov = 0)

Robust MV-IVW <- mr_mvivw(mvmr_input_mr, model = "random", robust = TRUE)

wtd_median <- mr_mvmedian(mvmr_input_mr)

pres <- pleiotropy_mvmr(r_input = mvmr_input_mr, gencov  = 0)




#########################################
######### 3. MR-PRESSO analysis #########

library(MRPRESSO)

mr_presso_result <- mr_presso(
  BetaOutcome = "beta.outcome",  
  BetaExposure = c("beta_LST_Female", beta_col), 
  SdOutcome = "se.outcome",    
  SdExposure = c("se_LST_Female", se_col),  
  # 检验参数（文献通用设置）
  OUTLIERtest = TRUE,      
  DISTORTIONtest = TRUE,   
  data = as.data.frame(mvmr_data_clean),   
  NbDistribution = 5000,    
  SignifThreshold = 0.05   
)
outlier_test <- mr_presso_result$`MR-PRESSO results`$`Outlier Test`



################################
####### 4. NLMR analysis #######

library(SUMnlmr)

#double ranking method
summ_covar <- create_nlmr_summary(y = nonlinear_use$Endometriosis,x = nonlinear_use$LST,g = nonlinear_use$grsW,
                                 covar = matrix(data = c(nonlinear_use$age,
                                                         nonlinear_use$array,
                                                         nonlinear_use$PC1,
                                                         nonlinear_use$PC2,
                                                         nonlinear_use$PC3,
                                                         nonlinear_use$PC4,
                                                         nonlinear_use$PC5,
                                                         nonlinear_use$PC6,
                                                         nonlinear_use$PC7,
                                                         nonlinear_use$PC8,
                                                         nonlinear_use$PC9,
                                                         nonlinear_use$PC10), 
                                                ncol=12),
                                 family = "binomial",q = 4,strata_method = "ranked",extra_statistics = FALSE,report_GR = TRUE)

model <- with(summ_covar$summary, piecewise_summ_mr(by,bx,byse,bxse,xmean,xmin,xmax,
                                                    ci="bootstrap_se",
                                                    nboot=1000,
                                                    fig= TRUE,
                                                    family="gaussian",
                                                    seed=1234))
model_bxby<-summ_covar$summary
model_results<-data.frame(model$lace)
model_results$OR<-round(exp(model_results$beta),2)
model_results$OR_lci<-round(exp(model_results$lci),2)
model_results$OR_uci<-round(exp(model_results$uci),2)
model_results<-cbind(model_bxby,model_results)


