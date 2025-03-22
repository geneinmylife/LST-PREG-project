library(devtools)
library(ggplot2)
library(TwoSampleMR)
library(readxl)
library(writexl)
library(dbplyr)
library(tidyverse)
library(plyr)
library(dplyr)
library(MungeSumstats)
library(data.table)

######## combine the sex-specific-LST GWAS outcome 

folder_path <- "/Volumes/T7 SS/female24h/bgen"  
files <- list.files(folder_path, full.names = TRUE)
data_list <- lapply(files, read.table, header = TRUE)
merged_data_female <- bind_rows(data_list)
head(merged_data_female)
merged_data_female_dup <- merged_data_female[!duplicated(merged_data_female$SNP),]
write_tsv(merged_data_female_dup,"LST_female_24.txt")


folder_path <- "/Volumes/T7 SS/male24h/bgen" 
files <- list.files(folder_path, full.names = TRUE)
data_list <- lapply(files, read.table, header = TRUE)
merged_data_male <- bind_rows(data_list)
head(merged_data_male)
merged_data_male_dup <- merged_data_male[!duplicated(merged_data_male$SNP),]
write_tsv(merged_data_male_dup,"LST_male_24.txt")

#Combining binary LST GWAS outcome is similar.


