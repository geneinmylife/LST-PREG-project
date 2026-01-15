
########## REGENIE-Step1 ##########
#!/bin/bash
source /share/apps/anaconda3/etc/profile.d/conda.sh
conda activate regenie_env

GENOdir="/data/UKBB/qc_genotype/output"
SNP_DIR="/share/home/SNP/output"
PHENOdir="/share/home/PHENO/input"
OUTdir="/share/home/OUT/output"
regenie \
  --step 1 \
  --pgen ${GENOdir}/ukb_allchr_merged \
  --extract ${SNP_DIR}/combined_ld_pruned_snplist.txt \
  --phenoFile ${PHENOdir}/pheno.txt \
  --covarFile ${PHENOdir}/pheno.txt \
  --phenoCol cluster \
  --covarColList age,sex,array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --keep ${PHENOdir}/sample_id.txt \
  --bt \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix tmp_rg \
  --out ${OUTdir}/step1_bin

########## REGENIE-Step2 ##########
#!/bin/bash
source /share/apps/anaconda3/etc/profile.d/conda.sh
conda activate regenie_env

GENOdir="/data/UKBB/qc_genotype/output"
PHENOdir="/share/home/PHENO/input"
Step1dir="/share/home/OUT/output"
Step2dir="/share/home/OUT2/output"

regenie \
  --step 2 \
  --pgen ${GENOdir}/ukb_allchr_merged \
  --phenoFile ${PHENOdir}/pheno.txt \
  --covarFile ${PHENOdir}/pheno.txt \
  --phenoCol cluster \
  --covarColList age,sex,array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --keep ${PHENOdir}/sample_id.txt \
  --bsize 200 \
  --bt \
  --firth \
  --approx \
  --pThresh 0.999999 \
  --pred ${Step1dir}/step1_bin_pred.list \
  --out ${Step2dir}/step2_bin
done
