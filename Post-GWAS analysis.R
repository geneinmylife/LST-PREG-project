## Content Checklist ##
# 1.Heritability
# 2. cross-trait genetic correlations


library(ldscr)

######## 1.Heritability

female_lst$Z <- female_lst$BETA / female_lst$SE
male_lst$Z <- male_lst$BETA / male_lst$SE
cc_female <- ldscr::ldsc_h2(munged_sumstats = female_lst, ancestry = "EUR")
cc_male <- ldscr::ldsc_h2(munged_sumstats = male_lst, ancestry = "EUR")


######## 2. cross-trait genetic correlations

rg_res <- ldscr::ldsc_rg(
munged_sumstats = list(
  "LST" = female_lst,
  "PA" =  female_pa,
  "SDW" =  female_sdw,
  "SDC" = female_sdc,
  "TSD" = female_tsd
),
ancestry = "EUR"
)

