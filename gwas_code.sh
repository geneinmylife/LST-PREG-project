
########## REGENIE-Step1 ##########
conda activate regenie_env

GENOdir="/pathway/to/genotype/data"
SNP_DIR="/pathway/to/snplist"
PHENOdir="/pathway/to/phenotype/data"
OUTdir="/pathway/to/output"

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
conda activate regenie_env

GENOdir="/pathway/to/genotype/data"
SNP_DIR="/pathway/to/snplist"
PHENOdir="/pathway/to/phenotype/data"
OUTdir="/pathway/to/output"

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
