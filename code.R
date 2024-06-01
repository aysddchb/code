library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(MRPRESSO)

setwd('C:\\Users\\86182\\AppData\\Local\\R\\win-library\\4.3\\TwoSampleMR')
uc_IIBDGC <- fread("E:\\刘国威\\data\\uc\\convert_ieu-a-32.vcf.gz.csv",select = c(1,2,3,4,5,6,7,8,9),col.names = c('chr','pos','SNP','other_allele','effect_allele','eaf','beta','se','p'))
uc_ebi <- fread('') # get from ieu
uc_finn <- fread("E:\\刘国威\\data\\uc\\finngen_R10_K11_UC_STRICT2.gz",header = T, sep = '\t', fill = T, select = c(3,4,5,7,9,10,11),col.names = c('other_allele','effect_allele','SNP','p','beta','se','eaf'))
ich_ISGC <- fread("E:\\刘国威\\data\\ich\\ICH_GWAS_phase1_finalresults_allICH",select = c(1:4,6:8),col.names = c('SNP','effect_allele','other_allele','eaf','beta','se','p'))
ich_finn <- fread("E:\\刘国威\\data\\ich\\finngen_R10_I9_ICH.gz",header = T, sep = '\t', fill = T, select = c(3,4,5,7,9,10,11),col.names = c('other_allele','effect_allele','SNP','p','beta','se','eaf'))
mdd_pgc <- fread("E:\\刘国威\\data\\depressiom\\PGC_UKB_depression_genome-wide.txt",col.names = c('SNP','effect_allele','other_allele','eaf','beta','se','p'))
htn_ukb <- fread("E:\\刘国威\\data\\hypertension\\convert_ukb-b-12493.vcf.gz.csv",select = c(3,4,5,6,7,8,9),col.names = c('SNP','other_allele','effect_allele','eaf','beta','se','p'))

uc_IIBDGC <- fread("E:\\刘国威\\data\\uc\\convert_ieu-a-32.vcf.gz.csv")

write.csv(mdd_pgc,'E:\\刘国威\\data\\depressiom\\MDD.csv')

#### uc as exposure:
###### uc_IIBDGC as exposure:
b <- subset(uc_IIBDGC,p<5e-08)
write.csv(b, file="exposure.csv")
bmi<-system.file("exposure.csv",package = "TwoSampleMR")
x1<-read_exposure_data(
  filename = bmi,
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "p",
  clump = T
)
########## to ich_ISGC:
d<-merge(x1,ich_ISGC,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'UC(IIBDGC)'
dat$outcome <- 'ICH(ISGC)'
mr <- mr(dat)
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_uc(IIBDGC)_ich(ISGC).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\uc(IIBDGC)_ich(ISGC).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_uc(IIBDGC)_ich(ISGC).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_uc(IIBDGC)_ich(ISGC).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_uc(IIBDGC)_ich(ISGC).pdf'),width = 5,height = 4)
#########  to ich_finn:
d<-merge(x1,ich_finn,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'UC(IIBDGC)'
dat$outcome <- 'ICH(finn)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_uc(IIBDGC)_ich(finn).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\uc(IIBDGC)_ich(finn).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_uc(IIBDGC)_ich(finn).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_uc(IIBDGC)_ich(finn).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_uc(IIBDGC)_ich(finn).pdf'),width = 5,height = 4)
#########  to mdd_pgc:
d<-merge(x1,mdd_pgc,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'UC(IIBDGC)'
dat$outcome <- 'MDD(PGC)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_uc(IIBDGC)_ICH(Finn).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\uc(IIBDGC)_ICH(Finn).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_uc(IIBDGC)_mdd(pgc).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_uc(IIBDGC)_mdd(pgc).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_uc(IIBDGC)_mdd(pgc).pdf'),width = 5,height = 4)
#########  to htn_ukb:
d<-merge(x1,htn_ukb,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'UC(IIBDGC)'
dat$outcome <- 'HTN(UKB)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_uc(IIBDGC)_htn(ukb).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\uc(IIBDGC)_htn(ukb).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_uc(IIBDGC)_htn(ukb).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_uc(IIBDGC)_htn(ukb).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_uc(IIBDGC)_htn(ukb).pdf'),width = 5,height = 4)
###### uc_ebi as exposure:
x1 <- extract_instruments(outcomes = 'ebi-a-GCST90018933', p1=1, clump = F)
########## to ich_ISGC:
d<-merge(x1,ich_ISGC,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'UC(EBI)'
dat$outcome <- 'dat$outcome <- 'ICH(Finn)'
mr <- mr(dat)'
mr <- mr(dat)
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_UC(EBI)_ich(ISGC).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\UC(EBI)_ich(ISGC).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_UC(EBI)_ich(ISGC).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_UC(EBI)_ich(ISGC).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_UC(EBI)_ich(ISGC).pdf'),width = 5,height = 4)
#########  to ich_finn:
d<-merge(x1,ich_finn,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'UC(EBI)'
dat$outcome <- 'ICH(finn)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_UC(EBI)_ich(finn).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\UC(EBI)_ich(finn).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_UC(EBI)_ich(finn).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_UC(EBI)_ich(finn).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_UC(EBI)_ich(finn).pdf'),width = 5,height = 4)
#########  to mdd_pgc:
d<-merge(x1,mdd_pgc,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'UC(EBI)'
dat$outcome <- 'MDD(PGC)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_UC(EBI)_mdd(pgc).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\UC(EBI)_mdd(pgc).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_UC(EBI)_mdd(pgc).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_UC(EBI)_mdd(pgc).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_UC(EBI)_mdd(pgc).pdf'),width = 5,height = 4)
#########  to htn_ukb:
d<-merge(x1,htn_ukb,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'UC(EBI)'
dat$outcome <- 'HTN(UKB)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_UC(EBI)_htn(ukb).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\UC(EBI)_htn(ukb).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_UC(EBI)_htn(ukb).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_UC(EBI)_htn(ukb).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_UC(EBI)_htn(ukb).pdf'),width = 5,height = 4)
##### uc_finn as exposure:
b <- subset(uc_finn,p<5e-08)
write.csv(b, file="exposure.csv")
bmi<-system.file("exposure.csv",package = "TwoSampleMR")
x1<-read_exposure_data(
  filename = bmi,
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "p",
  clump = T
)
#########  to ich_ISGC:
d<-merge(x1,ich_ISGC,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'UC(Finn)'
dat$outcome <- 'dat$outcome <- 'ICH(Finn)'
mr <- mr(dat)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_UC(Finn)_ich(ISGC).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\UC(Finn)_ich(ISGC).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_UC(Finn)_ich(ISGC).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_UC(Finn)_ich(ISGC).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_UC(Finn)_ich(ISGC).pdf'),width = 5,height = 4)
#########  to mdd_pgc:
d<-merge(x1,mdd_pgc,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'UC(Finn)'
dat$outcome <- 'MDD(PGC)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_UC(Finn)_mdd(pgc).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\UC(Finn)_mdd(pgc).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_UC(Finn)_mdd(pgc).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_UC(Finn)_mdd(pgc).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_UC(Finn)_mdd(pgc).pdf'),width = 5,height = 4)
#########  to htn_ukb:
d<-merge(x1,htn_ukb,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'UC(Finn)'
dat$outcome <- 'HTN(UKB)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_UC(Finn)_htn(ukb).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\UC(Finn)_htn(ukb).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_UC(Finn)_htn(ukb).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_UC(Finn)_htn(ukb).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_UC(Finn)_htn(ukb).pdf'),width = 5,height = 4)

##### ich_ISGC as exposure:
b <- subset(ich_ISGC,p<5e-06)
write.csv(b, file="exposure.csv")
bmi<-system.file("exposure.csv",package = "TwoSampleMR")
x1<-read_exposure_data(
  filename = bmi,
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "p",
  clump = T
)
########## to uc_IIBDGC:
d<-merge(x1,uc_IIBDGC,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'dat$outcome <- 'ICH(Finn)'
mr <- mr(dat)'
dat$outcome <- 'UC(IIBDGC)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_dat$outcome <- 'ICH(Finn)'
mr <- mr(dat)_UC(IIBDGC).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\dat$outcome <- 'ICH(Finn)'
mr <- mr(dat)_UC(IIBDGC).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_ICH(ISGC)_UC(IIBDGC).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_ICH(ISGC)_UC(IIBDGC).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_ICH(ISGC)_UC(IIBDGC).pdf'),width = 5,height = 4)
#########  to uc_ebi:
outcome_dat <- extract_outcome_data(snps = x1$SNP,outcomes = 'ebi-a-GCST90018933')
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'ICH(ISGC)'
dat$outcome <- 'UC(EBI)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_ICH(ISGC)_UC(EBI).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\ICH(ISGC)_UC(EBI).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_ICH(ISGC)_UC(EBI).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_ICH(ISGC)_UC(EBI).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_ICH(ISGC)_UC(EBI).pdf'),width = 5,height = 4)
#########  to uc_finn:
d<-merge(x1,uc_finn,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'ICH(ISGC)'
dat$outcome <- 'UC(Finn)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_ICH(ISGC)_UC(Finn).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\ICH(ISGC)_UC(Finn).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_ICH(ISGC)_UC(Finn).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_ICH(ISGC)_UC(Finn).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_ICH(ISGC)_UC(Finn).pdf'),width = 5,height = 4)
#########  to mdd_pgc:
d<-merge(x1,mdd_pgc,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'ICH(ISGC)'
dat$outcome <- 'MDD(PGC)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_ICH(ISGC)_mdd(pgc).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\ICH(ISGC)_mdd(pgc).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_ICH(ISGC)_mdd(pgc).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_ICH(ISGC)_mdd(pgc).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_ICH(ISGC)_mdd(pgc).pdf'),width = 5,height = 4)
#########  to htn_ukb:
d<-merge(x1,htn_ukb,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'ICH(ISGC)'
dat$outcome <- 'HTN(UKB)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_ICH(ISGC)_htn(ukb).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\ICH(ISGC)_htn(ukb).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_ICH(ISGC)_htn(ukb).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_ICH(ISGC)_htn(ukb).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_ICH(ISGC)_htn(ukb).pdf'),width = 5,height = 4)
##### ich_finn as exposure:
b <- subset(ich_finn,p<5e-06)
write.csv(b, file="exposure.csv")
bmi<-system.file("exposure.csv",package = "TwoSampleMR")
x1<-read_exposure_data(
  filename = bmi,
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "p",
  clump = T
)
########## to uc_IIBDGC:
d<-merge(x1,uc_IIBDGC,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'ICH(Finn)'
dat$outcome <- 'UC(IIBDGC)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_ICH(Finn)_UC(IIBDGC).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\ICH(Finn)_UC(IIBDGC).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_ICH(Finn)_UC(IIBDGC).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_ICH(Finn)_UC(IIBDGC).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_ICH(Finn)_UC(IIBDGC).pdf'),width = 5,height = 4)
#########  to uc_ebi:
outcome_dat <- extract_outcome_data(snps = x1$SNP,outcomes = 'ebi-a-GCST90018933')
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'ICH(Finn)'
dat$outcome <- 'UC(EBI)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_ICH(Finn)_UC(EBI).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\ICH(Finn)_UC(EBI).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_ICH(Finn)_UC(EBI).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_ICH(Finn)_UC(EBI).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_ICH(Finn)_UC(EBI).pdf'),width = 5,height = 4)
#########  to mdd_pgc:
d<-merge(x1,mdd_pgc,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'ICH(Finn)'
dat$outcome <- 'MDD(PGC)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_ICH(Finn)_mdd(pgc).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\ICH(Finn)_mdd(pgc).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_ICH(Finn)_mdd(pgc).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_ICH(Finn)_mdd(pgc).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_ICH(Finn)_mdd(pgc).pdf'),width = 5,height = 4)
#########  to htn_ukb:
d<-merge(x1,htn_ukb,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'ICH(Finn)'
dat$outcome <- 'HTN(UKB)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_ICH(Finn)_htn(ukb).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\ICH(Finn)_htn(ukb).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_ICH(Finn)_htn(ukb).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_ICH(Finn)_htn(ukb).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_ICH(Finn)_htn(ukb).pdf'),width = 5,height = 4)

##### mdd_pgc as exposure:
b <- subset(mdd_pgc,p<5e-08)
write.csv(b, file="exposure.csv")
bmi<-system.file("exposure.csv",package = "TwoSampleMR")
x1<-read_exposure_data(
  filename = bmi,
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "p",
  clump = T
)
########## to uc_IIBDGC:
d<-merge(x1,uc_IIBDGC,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'MDD(PGC)'
dat$outcome <- 'UC(IIBDGC)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_MDD(PGC)_UC(IIBDGC).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\MDD(PGC)_UC(IIBDGC).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_MDD(PGC)_UC(IIBDGC).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_MDD(PGC)_UC(IIBDGC).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_MDD(PGC)_UC(IIBDGC).pdf'),width = 5,height = 4)
#########  to uc_ebi:
outcome_dat <- extract_outcome_data(snps = x1$SNP,outcomes = 'ebi-a-GCST90018933')
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'MDD(PGC)'
dat$outcome <- 'UC(EBI)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_MDD(PGC)_UC(EBI).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\MDD(PGC)_UC(EBI).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_MDD(PGC)_UC(EBI).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_MDD(PGC)_UC(EBI).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_MDD(PGC)_UC(EBI).pdf'),width = 5,height = 4)
#########  to uc_finn:
d<-merge(x1,uc_finn,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'MDD(PGC)'
dat$outcome <- 'UC(Finn)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_MDD(PGC)_UC(Finn).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\MDD(PGC)_UC(Finn).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_MDD(PGC)_UC(Finn).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_MDD(PGC)_UC(Finn).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_MDD(PGC)_UC(Finn).pdf'),width = 5,height = 4)
#########  to ich_ISGC:
d<-merge(x1,ich_ISGC,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'MDD(PGC)'
dat$outcome <- 'ICH(ISGC)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_MDD(PGC)_ICH(ISGC).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\MDD(PGC)_ICH(ISGC).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_MDD(PGC)_ICH(ISGC).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_MDD(PGC)_ICH(ISGC).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_MDD(PGC)_ICH(ISGC).pdf'),width = 5,height = 4)
#########  to ich_finn:
d<-merge(x1,ich_finn,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'MDD(PGC)'
dat$outcome <- 'ICH(Finn)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_MDD(PGC)_ICH(Finn).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\MDD(PGC)_ICH(Finn).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_MDD(PGC)_ICH(Finn).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_MDD(PGC)_ICH(Finn).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_MDD(PGC)_ICH(Finn).pdf'),width = 5,height = 4)
#########  to htn_ukb:
d<-merge(x1,htn_ukb,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'MDD(PGC)'
dat$outcome <- 'HTN(UKB)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_MDD(PGC)_htn(ukb).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\MDD(PGC)_htn(ukb).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_MDD(PGC)_htn(ukb).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_MDD(PGC)_htn(ukb).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_MDD(PGC)_htn(ukb).pdf'),width = 5,height = 4)
##### htn_ukb as exposure:
b <- subset(htn_ukb,p<5e-08)
write.csv(b, file="exposure.csv")
bmi<-system.file("exposure.csv",package = "TwoSampleMR")
x1<-read_exposure_data(
  filename = bmi,
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "p",
  clump = T
)
########## to uc_IIBDGC:
d<-merge(x1,uc_IIBDGC,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'HTN(UKB)'
dat$outcome <- 'UC(IIBDGC)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_HTN(UKB)_UC(IIBDGC).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\HTN(UKB)_UC(IIBDGC).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_HTN(UKB)_UC(IIBDGC).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_HTN(UKB)_UC(IIBDGC).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_HTN(UKB)_UC(IIBDGC).pdf'),width = 5,height = 4)
#########  to uc_ebi:
outcome_dat <- extract_outcome_data(snps = x1$SNP,outcomes = 'ebi-a-GCST90018933')
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'HTN(UKB)'
dat$outcome <- 'UC(EBI)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_HTN(UKB)_UC(EBI).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\HTN(UKB)_UC(EBI).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_HTN(UKB)_UC(EBI).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_HTN(UKB)_UC(EBI).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_HTN(UKB)_UC(EBI).pdf'),width = 5,height = 4)
#########  to uc_finn:
d<-merge(x1,uc_finn,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'HTN(UKB)'
dat$outcome <- 'UC(Finn)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_HTN(UKB)_UC(Finn).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\HTN(UKB)_UC(Finn).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_HTN(UKB)_UC(Finn).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_HTN(UKB)_UC(Finn).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_HTN(UKB)_UC(Finn).pdf'),width = 5,height = 4)
#########  to ich_ISGC:
d<-merge(x1,ich_ISGC,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'HTN(UKB)'
dat$outcome <- 'ICH(ISGC)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_HTN(UKB)_ICH(ISGC).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\HTN(UKB)_ICH(ISGC).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_HTN(UKB)_ICH(ISGC).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_HTN(UKB)_ICH(ISGC).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_HTN(UKB)_ICH(ISGC).pdf'),width = 5,height = 4)
#########  to ich_finn:
d<-merge(x1,ich_finn,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'HTN(UKB)'
dat$outcome <- 'ICH(Finn)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_HTN(UKB)_ICH(Finn).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\HTN(UKB)_ICH(Finn).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_HTN(UKB)_ICH(Finn).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_HTN(UKB)_ICH(Finn).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_HTN(UKB)_ICH(Finn).pdf'),width = 5,height = 4)
#########  to mdd_pgc:
d<-merge(x1,mdd_pgc,by="SNP")
write.csv(d,file = "outcome.csv")
out<-system.file("outcome.csv",package = "TwoSampleMR")
outcome_dat<-read_outcome_data(snps = x1$SNP,
                               filename = out,
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "eaf",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p")
dat<-harmonise_data(exposure_dat = x1,outcome_dat = outcome_dat)
dat$exposure <- 'HTN(UKB)'
dat$outcome <- 'MDD(PGC)'
mr <- mr(dat)
mr
mr_odd <- generate_odds_ratios(mr_res = mr)
writexl::write_xlsx(dat,'E:\\刘国威\\结果\\dat_HTN(UKB)_MDD(PGC).xlsx')
writexl::write_xlsx(mr_odd,'E:\\刘国威\\结果\\HTN(UKB)_MDD(PGC).xlsx')
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
mr_scatter_plot(mr,dat)
ggsave(paste0('E:\\刘国威\\结果\\scatter_HTN(UKB)_MDD(PGC).pdf'),width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave(paste0('E:\\刘国威\\结果\\funnel_HTN(UKB)_MDD(PGC).pdf'),width = 5,height = 4)
mr_pleiotropy_test(dat) 
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste0('E:\\刘国威\\结果\\leaveoneout_HTN(UKB)_MDD(PGC).pdf'),width = 5,height = 4)

######### meta-analyse
library(meta)
data1 <- readxl::read_excel("E:\\刘国威\\结果\\Meta_htn_ich.xlsx")
lnor<- log(data1[,'or'][[1]])
lnuci<- log(data1[,'or_uci95'][[1]])
lnlci<- log(data1[,'or_lci95'][[1]])
selnor<- (lnuci-lnlci)/(2*1.96)
pfs <- metagen(lnor,selnor,sm="OR",data=data1,studlab = data1$outcome)
summary(pfs)
forest(pfs,family='sans',col.square = 'skyblue',col.diamond.common = 'maroon',col.diamond.lines.random = 'maroon')


# 绘制森林图
library(forestplot)
d <- readxl::read_excel("E:\\刘国威\\结果\\HTN(UKB)_ICH(Finn).xlsx")
e <- readxl::read_excel("E:\\刘国威\\结果\\HTN(UKB)_ICH(ISGC).xlsx")
d <- d[3,]
e <- e[3,]
f <- rbind(d,e)

files <- list.files("E:\\刘国威\\结果",pattern = '.xlsx')[c(34:54,57:67)]
setwd("E:\\刘国威\\结果")
for (v in files) {
  data <- readxl::read_excel(v)[3,]
  d <- rbind(d,data)
}
library(tidyr)
d$aaa <- "("; d$bbb <- "-"; d$ccc <- ")"
d$or <- round(d$or,3)
d$or_lci95 <- round(d$or_lci95,3)
d$or_uci95 <- round(d$or_uci95,3)
d <-  unite(d, "or (95%CI)", c(or, aaa, or_lci95, bbb, or_uci95, ccc),sep = '',remove = F)
d[,c(1,2)] <- d[,c(4,3)]
names(d)[c(1,2)] <- names(d)[c(4,3)]
d <- d[-c(1:)]
d <- d[,-c(3,4,5,17:19)]
d_htn_ich <- d[c(1:3),]
d <- d[-c(1:3,5:7),]

writexl::write_xlsx(d,'forest.xlsx')
d <- readxl::read_excel('forest.xlsx')
forestplot(
  labeltext = as.matrix(d[,c(1,2,3,9,6)]),
  graph.pos = 4,
  mean = d$or, 
  lower = d$or_lci95,
  upper = d$or_uci95,
  zero = 1,
  boxsize = 0.1,
  line.margin = 1,
  col = fpColors(box = "black", lines = "gray",
                 summary = "black", zero = "lightgray", text = "black",
                 axes = "black", hrz_lines = "black")
)

