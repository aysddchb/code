library(TwoSampleMR)
library(data.table)
library(MRInstruments)
library(MRPRESSO)
library(ggplot2)
library(tidyverse)

uc_IIBDGC <- fread("E:\\刘国威\\data\\uc\\convert_ieu-a-32.vcf.gz.csv",select = c(1,2,3,4,5,6,7,8,9),col.names = c('chr','pos','SNP','other_allele','effect_allele','eaf','beta','se','p'))
uc_ebi <- fread('') # get from ieu
uc_finn <- fread("E:\\刘国威\\data\\uc\\finngen_R10_K11_UC_STRICT2.gz",header = T, sep = '\t', fill = T, select = c(3,4,5,7,9,10,11),col.names = c('other_allele','effect_allele','SNP','p','beta','se','eaf'))
ich_ISGC <- fread("E:\\刘国威\\data\\ich\\ICH_GWAS_phase1_finalresults_allICH",select = c(1:4,6:8),col.names = c('SNP','effect_allele','other_allele','eaf','beta','se','p'))
ich_finn <- fread("E:\\刘国威\\data\\ich\\finngen_R10_I9_ICH.gz",header = T, sep = '\t', fill = T, select = c(3,4,5,7,9,10,11),col.names = c('other_allele','effect_allele','SNP','p','beta','se','eaf'))
mdd_pgc <- fread("E:\\刘国威\\data\\depressiom\\PGC_UKB_depression_genome-wide.txt",col.names = c('SNP','effect_allele','other_allele','eaf','beta','se','p'))
htn_ukb <- fread("E:\\刘国威\\data\\hypertension\\convert_ukb-b-12493.vcf.gz.csv",select = c(3,4,5,6,7,8,9),col.names = c('SNP','other_allele','effect_allele','eaf','beta','se','p'))
l <- data.frame()

uc_ebi <- data.frame()

setwd('E:\\刘国威\\药靶\\sbp')
sbp<-extract_instruments(outcome="ieu-b-38",clump = T,kb=100,r2=0.3)

ADRB1<-subset(sbp,chr.exposure==10 & pos.exposure>115803625-100000 & pos.exposure<115806663+100000)
CACNA1D<-subset(sbp,chr.exposure==3 & pos.exposure>53528638-100000 & pos.exposure<53847760+100000)
CACNB2<-subset(sbp,chr.exposure==10 & pos.exposure>18429353-100000 & pos.exposure<18832486+100000)
SLC12A2<-subset(sbp,chr.exposure==5 & pos.exposure>127419458-100000 & pos.exposure<127525369+100000)

# ADRB1_gwas<-subset(ADRB1_gwas,eaf.exposure>0.01)
# CACNA1D_gwas<-subset(CACNA1D_gwas,eaf.exposure>0.01)
# CACNB2_gwas<-subset(CACNB2_gwas,eaf.exposure>0.01)
# SLC12A2_gwas<-subset(SLC12A2_gwas,eaf.exposure>0.01)
disease <- c('uc_ebi','uc_IIBDGC','uc_finn','ich_ISGC','ich_finn','mdd_pgc','htn_ukb')
genes <- c('ADRB1','CACNA1D','CACNB2','SLC12A2')
for (dis in disease) {
  for (gene in genes) {
    g <- get(gene)
    k <- data.frame(id=1:10)
    if(dis!='uc_ebi'){
      di <- get(dis)
      d<-merge(g,di,by = "SNP")
      write.csv(d,file = "C:\\Users\\86182\\AppData\\Local\\R\\win-library\\4.3\\TwoSampleMR\\outcome.csv")
      out<-system.file("outcome.csv",package = "TwoSampleMR")
      outcome_dat<-read_outcome_data(snps = g$SNP,
                                     filename = out,
                                     sep = ",",
                                     snp_col = "SNP",
                                     beta_col = "beta",
                                     se_col = "se",
                                     eaf_col = "eaf",
                                     effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",
                                     pval_col = "p")
    }
    else if(dis=='uc_ebi'){
      outcome_dat <- extract_outcome_data(snps = g$SNP,
                                          outcomes = 'ebi-a-GCST90018933')
    }
    dat <- harmonise_data(g,outcome_dat)
    dat$exposure <- gene
    dat$outcome <- dis
    writexl::write_xlsx(dat,paste0('dat_',gene,'_',dis,'.xlsx'))
    res3 <- mr(dat)
    res3
    mrTab=generate_odds_ratios(res3)
    mrTab$orDrug=1/mrTab$or
    beta=log(mrTab$orDrug)
    mrTab$or_lci95Drug=exp(beta-1.96*mrTab$se)
    mrTab$or_uci95Drug=exp(beta+1.96*mrTab$se)
    writexl::write_xlsx(mrTab,paste0('mrdrug_',gene,'_',dis,'.xlsx'))
    k$or <- mrTab$orDrug[mrTab$method=='Inverse variance weighted']
    k$or_l <- mrTab$or_lci95Drug[mrTab$method=='Inverse variance weighted']
    k$or_u <- mrTab$or_uci95Drug[mrTab$method=='Inverse variance weighted']
    k$p <- mrTab$pval[mrTab$method=='Inverse variance weighted']
    k$pleiotropy <- mr_pleiotropy_test(dat)$pval
    k$heterogeneity <- mr_heterogeneity(dat)$Q_pval[2]
    presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
    k$presso <- presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
    k$exposure <- gene
    k$outcome <- dis
    k$nsnp <- mrTab[mrTab$method=='Inverse variance weighted',6]
    l <- rbind(l,k[1,])
    mr_scatter_plot(res3,dat)
    ggsave(paste0('scatter_',gene,'_',dis,'.png'),width = 5,height = 4)
    mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
    ggsave(paste0('funnel_',gene,'_',dis,'.png'),width = 5,height = 4)
    mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
    ggsave(paste0('leaveoneout_',gene,'_',dis,'.png'),width = 5,height = 4)
  }
}
writexl::write_xlsx(l,'forest.xlsx')

dat <- harmonise_data(ADRB1,ECM1_UC_gwas)
writexl::write_xlsx(dat,'dat_ADRB1_uc_ebi.xlsx')
res3 <- mr(dat)
res3
mrTab=generate_odds_ratios(res3)
mrTab$orDrug=1/mrTab$or
beta=log(mrTab$orDrug)
mrTab$or_lci95Drug=exp(beta-1.96*mrTab$se)
mrTab$or_uci95Drug=exp(beta+1.96*mrTab$se)
writexl::write_xlsx(mrTab,'drugmr_ADRB1_uc_ebi.xlsx')
k$or <- mrTab$orDrug[mrTab$method=='Inverse variance weighted']
k$or_l <- mrTab$or_lci95Drug[mrTab$method=='Inverse variance weighted']
k$or_u <- mrTab$or_uci95Drug[mrTab$method=='Inverse variance weighted']
k$p <- mrTab$pval[mrTab$method=='Inverse variance weighted']
k$pleiotropy <- mr_pleiotropy_test(dat)$pval
k$heterogeneity <- mr_heterogeneity(dat)$Q_pval[2]
presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
k$presso <- presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
k$exposure <- 'ADRB1'
k$outcome <- 'uc_ebi'
k$nsnp <- mrTab[mrTab$method=='Inverse variance weighted',6]
l <- rbind(l,k[1,])
mr_scatter_plot(res3,dat)
ggsave('scatter_ADRB1_uc_ebi.png',width = 5,height = 4)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
ggsave('funnel_ADRB1_uc_ebi.png',width = 5,height = 4)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
ggsave('leaveoneout_ADRB1_uc_ebi.png',width = 5,height = 4)




setwd('E:\\刘国威\\药靶\\dbp')
dbp_gwas<-extract_instruments(outcome="ieu-b-39",clump = T,kb=100,r2=0.3)

ADRB1<-subset(dbp_gwas,chr.exposure==10 & pos.exposure>115803625-100000 & pos.exposure<115806663+100000)
CACNA1D<-subset(dbp_gwas,chr.exposure==3 & pos.exposure>53528638-100000 & pos.exposure<53847760+100000)
CACNB2<-subset(dbp_gwas,chr.exposure==10 & pos.exposure>18429353-100000 & pos.exposure<18832486+100000)
ll <- data.frame()
# ADRB1_gwas<-subset(ADRB1_gwas,eaf.exposure>0.01)
# CACNA1D_gwas<-subset(CACNA1D_gwas,eaf.exposure>0.01)
# CACNB2_gwas<-subset(CACNB2_gwas,eaf.exposure>0.01)

disease <- c('uc_ebi','uc_IIBDGC','uc_finn','ich_ISGC','ich_finn','mdd_pgc','htn_ukb')
genes <- c('ADRB1','CACNA1D','CACNB2')
for (dis in disease) {
  for (gene in genes) {
    g <- get(gene)
    k <- data.frame(id=1:10)
    if(dis!='uc_ebi'){
      di <- get(dis)
      d<-merge(g,di,by = "SNP")
      write.csv(d,file = "C:\\Users\\86182\\AppData\\Local\\R\\win-library\\4.3\\TwoSampleMR\\outcome.csv")
      out<-system.file("outcome.csv",package = "TwoSampleMR")
      outcome_dat<-read_outcome_data(snps = g$SNP,
                                     filename = out,
                                     sep = ",",
                                     snp_col = "SNP",
                                     beta_col = "beta",
                                     se_col = "se",
                                     eaf_col = "eaf",
                                     effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",
                                     pval_col = "p")
    }
    else if(dis=='uc_ebi'){
      outcome_dat <- extract_outcome_data(snps = g$SNP,
                                          outcomes = 'ebi-a-GCST90018933')
    }
    dat <- harmonise_data(g,outcome_dat)
    dat$exposure <- gene
    dat$outcome <- dis
    writexl::write_xlsx(dat,paste0('dat_',gene,'_',dis,'.xlsx'))
    res3 <- mr(dat)
    res3
    mrTab=generate_odds_ratios(res3)
    mrTab$orDrug=1/mrTab$or
    beta=log(mrTab$orDrug)
    mrTab$or_lci95Drug=exp(beta-1.96*mrTab$se)
    mrTab$or_uci95Drug=exp(beta+1.96*mrTab$se)
    writexl::write_xlsx(mrTab,paste0('mrdrug_',gene,'_',dis,'.xlsx'))
    k$or <- mrTab$orDrug[mrTab$method=='Inverse variance weighted']
    k$or_l <- mrTab$or_lci95Drug[mrTab$method=='Inverse variance weighted']
    k$or_u <- mrTab$or_uci95Drug[mrTab$method=='Inverse variance weighted']
    k$p <- mrTab$pval[mrTab$method=='Inverse variance weighted']
    k$pleiotropy <- mr_pleiotropy_test(dat)$pval
    k$heterogeneity <- mr_heterogeneity(dat)$Q_pval[2]
    presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
    k$presso <- presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
    k$exposure <- gene
    k$outcome <- dis
    k$nsnp <- mrTab[mrTab$method=='Inverse variance weighted',6]
    ll <- rbind(ll,k[1,])
    mr_scatter_plot(res3,dat)
    ggsave(paste0('scatter_',gene,'_',dis,'.png'),width = 5,height = 4)
    mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
    ggsave(paste0('funnel_',gene,'_',dis,'.png'),width = 5,height = 4)
    mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
    ggsave(paste0('leaveoneout_',gene,'_',dis,'.png'),width = 5,height = 4)
  }
}
writexl::write_xlsx(ll,'forest.xlsx')



ECM1_UC_gwas <- extract_outcome_data(snps = CACNB2_gwas$SNP,
                                     outcomes = 'ebi-a-GCST90018806')
dat <- harmonise_data(dbp_gwas,ECM1_UC_gwas)
res3 <- mr(dat)
res3
mrTab=generate_odds_ratios(res3)
mrTab$orDrug=1/mrTab$or
beta=log(mrTab$orDrug)
mrTab$or_lci95Drug=exp(beta-1.96*mrTab$se)
mrTab$or_uci95Drug=exp(beta+1.96*mrTab$se)


# forest
ll$or <- round(ll$or,3)
ll$or_l <- round(ll$or_l,3)
ll$or_u <- round(ll$or_u,3)
ll$aaa <- "("; ll$bbb <- "-"; ll$ccc <- ")"
ll <-  unite(ll, "or (95%CI)", c(or, aaa, or_l, bbb, or_u, ccc),sep = '',remove = F)
ll <- rbind(names(ll),ll)
ll[1,][c(3,4,5)] <- NA
ll$p[c(-1)] <- round(as.numeric(ll$p[c(-1)]),3)
library(forestplot)
forestplot(
  labeltext = as.matrix(ll[,c(10,11,12,2,6)]),
  graph.pos = 3,
  mean = as.numeric(ll$or),
  lower = as.numeric(ll$or_l),
  upper = as.numeric(ll$or_u),
  zero = 1,
  boxsize = 0.1,
  line.margin = 1,
  col = fpColors(box = "black", lines = "gray",
                  summary = "black", zero = "lightgray", text = "black",
                 axes = "black", hrz_lines = "black")
  )



il_10 <- read.csv("C:\\Users\\86182\\Downloads\\convert_prot-c-3722_49_2.vcf.gz.csv",header = T,fill = T)
b<-subset(il_10,pval<5e-05)
write.csv(b, file="C:\\Users\\86182\\AppData\\Local\\R\\win-library\\4.3\\TwoSampleMR\\exposure.csv")
#读取exposure数据
bmi<-system.file("exposure.csv",package = "TwoSampleMR")
#数据标准化
x1<-read_exposure_data(
  filename = bmi,
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect.allele",
  other_allele_col = "other.allele",
  eaf_col = "eaf",
  pval_col = "pval",
  clump = F
)

cccc <- subset(x1,chr.exposure==4 & pos.exposure>42313324-100000 & pos.exposure<42388442+100000)
ECM1_UC_gwas <- extract_outcome_data(snps = x1$SNP,
                                     outcomes = 'finn-b-K11_FIBROCHIRLIV')
dat <- harmonise_data(x1,ECM1_UC_gwas)
res3 <- mr(dat)
res3



