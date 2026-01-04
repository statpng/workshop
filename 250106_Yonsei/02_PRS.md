# Microbiome/PRS Analysis Workshop: 

## From QC/GWAS to PRS

이 파이프라인은 Raw Genotype 데이터의 전처리(QC, Imputation)부터 GWAS 분석을 수행하고, 이를 바탕으로 Polygenic Risk Score(PRS)를 계산하여 검증하는 전체 과정을 포함함.

## Part 1. Data Pre-processing & GWAS

여기서는 Phenotype 데이터 정리, 성별 확인, Imputation, PCA 산출, 그리고 최종 GWAS 분석까지 수행함.

### 1.1 Phenotype & Fam File Processing

데이터 인코딩 문제 해결하고 `.pheno`랑 `.fam` 파일 정리하는 단계임.

```r
library(tidyverse)
library(magrittr)

# Fixing broken "Proj Name" encoding: "MR \xc1ߵ\xb5\xc6\xf7\xb1\xe2" >> "MR"


# Pheno -------------------------------------------------------------------
df.pheno <- read.table("./data/data_org/SKKU_dataset_WB_FixEncoding.txt", header=TRUE) %>%
  mutate(FID = gsub("_", "", Array_Name),
         IID = FID,
         SEX = as.numeric(FEMALE == 0)) %>%
  select(FID, IID, everything(), -Array_Name)

# save it to ".pheno" file
write.table(df.pheno, file="./data/merge0422.pheno", quote = FALSE, row.names = FALSE)


# Fam ---------------------------------------------------------------------
df.fam <- read.table("./data/data_org/merge0422.fam", header=FALSE) %>%
  set_colnames(c("FID", "IID", "PID", "MID", "SEX", "PHENO")) %>%
  mutate(FID = gsub("_", "", FID),
         IID = FID
         ) %>%
  filter(IID %in% df.pheno$IID)



df.fam %>%
  write.table(file="./data/merge0422.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)

df.fam %>%
  select(FID, IID) %>%
  write.table(file="./data/merge0422.ID", quote = FALSE, row.names = FALSE, col.names = FALSE)

```

### 1.2 Format Conversion (PLINK)

기본 포맷 변환 작업임.

```r
# (Ped, Map) ----------------------------------------------------------------
system("cp ./data/data_org/merge0422.bim ./data/merge0422.bim")
system("cp ./data/data_org/merge0422.bed ./data/merge0422.bed")


# (bim, bed, fam) >> (ped, map), denoted by BBF to PM
system("~/plink --bfile ./data/merge0422 --keep ./data/merge0422.ID --recode --out ./data/merge0422")

# # (ped, map) >> (bim, bed, fam), denoted by PM to BBF
# system("~/plink --file ./data/merge0422 --out ./data/merge0422")
#
# # (ped, map) >> (vcf)
# system("~/plink --ped ./data/merge0422.ped --map ./data/merge0422.map --recode vcf --out ./data/merge0422_vcf")

```

### 1.3 Gender Check

성별 불일치 확인해서 이상한 샘플 걸러내는 작업.

```r
# Check Gender ------------------------------------------------------------
system("mkdir check-sex")
system("~/plink --ped ./data/merge0422.ped --check-sex --map ./data/merge0422.map --keep ./data/merge0422.ID --snps-only just-acgt --out ./check-sex/merge0422_check-sex")


df.checksex <- read.table("./check-sex/merge0422_check-sex.sexcheck", header=TRUE) %>%
  mutate(FSEX = case_when(abs(F)>0.8 ~ 1,
                          abs(F)<0.2 ~ 0,
                          .default = 9999))
df.nosex <- read.table("./check-sex/merge0422_check-sex.nosex")

df.sex_ambiguous <-
  df.checksex %>%
  filter(PEDSEX==0) # Ambiguous (no gender in .fam file)

# Appropriate
df.sex_ambiguous %>% {
  table(.$FSEX, .$SNPSEX)
}


df.SexDiscrepancy <-
  df.pheno %>%
  left_join(df.checksex, by="IID", suffix=c("", ".new")) %>%
  # filter(SEX != SEX)
  filter(SEX != FSEX)




read.table("./data/merge0422.ID", col.names = c("FID", "IID")) %>%
  filter(!IID %in% df.SexDiscrepancy$IID) %>%
  write.table(file="./data/merge0422.ID_sex", quote = FALSE, row.names = FALSE, col.names = FALSE)




table(df.fam$SEX, df.pheno$SEX)

```

### 1.4 Pre-Imputation QC

Imputation 하기 전에 빡세게 QC 돌리는 구간.

```r
# Pre-Imputation QC ------------------------------------------------------

# 1: pre-screening
system("mkdir Pre_QC")
system("~/plink --ped ./data/merge0422.ped --map ./data/merge0422.map --keep ./data/merge0422.ID_sex --snps-only just-acgt --geno 0.2 --mind 0.2 --recode --out ./Pre_QC/merge0422_geno0.2_mind0.2")

# 2: main SNP-QC
system("~/plink --ped ./Pre_QC/merge0422_geno0.2_mind0.2.ped --map ./Pre_QC/merge0422_geno0.2_mind0.2.map --chr 1-22 --snps-only just-acgt --geno 0.05 --maf 0.005 --hwe 0.000001 --recode --out ./Pre_QC/merge0422_geno0.2_mind0.2_maf0.005_geno0.05_hwe0.000001")

# 3: main Sample-QC
system("~/plink --ped ./Pre_QC/merge0422_geno0.2_mind0.2_maf0.005_geno0.05_hwe0.000001.ped --map ./Pre_QC/merge0422_geno0.2_mind0.2_maf0.005_geno0.05_hwe0.000001.map --genome --max 0.2 --mind 0.03 --recode --out ./Pre_QC/merge0422_geno0.2_mind0.2_maf0.005_geno0.05_hwe0.000001_genome0.2_mind0.03")

```

### 1.5 Imputation Preparation

서버에 올릴 VCF 파일 만드는 과정.

```r
# Imputation --------------------------------------------------------------


# Convert to vcf.gz -------------------------------------------------------

system("mkdir Imputation_vcf")
for( i in 1:22 ){
  system(paste0("~/plink --ped ./Pre_QC/merge0422_geno0.2_mind0.2_maf0.005_geno0.05_hwe0.000001_genome0.2_mind0.03.ped --map ./Pre_QC/merge0422_geno0.2_mind0.2_maf0.005_geno0.05_hwe0.000001_genome0.2_mind0.03.map --recode vcf --chr ", i, " --snps-only just-acgt --out ./Imputation_vcf/merge0422_geno0.2_mind0.2_maf0.005_geno0.05_hwe0.000001_genome0.2_mind0.03_chr",i))
  system(paste0("bcftools sort ./Imputation_vcf/merge0422_geno0.2_mind0.2_maf0.005_geno0.05_hwe0.000001_genome0.2_mind0.03_chr",i,".vcf -Oz -o ./Imputation_vcf/merge0422_geno0.2_mind0.2_maf0.005_geno0.05_hwe0.000001_genome0.2_mind0.03_chr",i,".vcf.gz"))
}




# Imputation server -------------------------------------------------------

# (1) https://coda.nih.go.kr/frt/index.do
# (2) 활용 - 한국인 임퓨테이션 서비스 - 서비스로 이동
#     Utilization - Korean Impression Service - Go to Service
# (3) Login: ID/PWD
# rsq Filter: off
# Phasing: Eagle v2.4
# Imputation: Minimac4
# Mode: QC & Imputation
# Check - "Generate Meta-imputation file"

```

### 1.6 Post-Imputation Processing

Imputation 끝나고 받은 파일 압축 풀고 정리하는 과정.

```r
# Post-QC after imputation ------------------------------------------------
PASSWORD = "i2Lgl5*wiBUpQY"

system("mkdir Imputatino_${PASSWORD}")

system('
PASSWORD="i2Lgl5*wiBUpQY"

for i in {1..22}
do
unzip -P ${PASSWORD} "./Imputation_${PASSWORD}/chr_${i}.zip"
done
')



system('
PASSWORD="i2Lgl5*wiBUpQY"

for i in {1..22}
do
gunzip -d ./Imputation_${PASSWORD}/chr${i}.info.gz
done
')
# gunzip -d ./Imputation_${PASSWORD}/chr${i}.dose.vcf.gz
# gunzip -d ./Imputation_${PASSWORD}/chr${i}.empiricalDose.vcf.gz

```

### 1.7 Post-Imputation QC & Merge

결과 파일 합치고(, MAF 등 기준) 다시 PLINK 파일로 변환.

```r
# Post-Imputation QC --------------------------------------------------

# # Imp > 0.8
# system('
# PASSWORD="kgjkNKJ?r7DoM4"
#
# for i in {1..22}
# do
# filename=chr${i}.dose
#
# bcftools view -i "R2>.8" -Oz ./Imputation_${PASSWORD}/$filename.vcf.gz > ./Imputation_${PASSWORD}/$filename.filtered.vcf.gz; tabix -p vcf ./Imputation_${PASSWORD}/$filename.filtered.vcf.gz;
#
# done
# ')

system("mkdir Post_QC")

system('
PATH="./Imputation_i2Lgl5*wiBUpQY"
cd $PATH

bcftools concat chr1.dose.vcf.gz chr2.dose.vcf.gz chr3.dose.vcf.gz chr4.dose.vcf.gz chr5.dose.vcf.gz chr6.dose.vcf.gz chr7.dose.vcf.gz chr8.dose.vcf.gz chr9.dose.vcf.gz chr10.dose.vcf.gz chr11.dose.vcf.gz chr12.dose.vcf.gz chr13.dose.vcf.gz chr14.dose.vcf.gz chr15.dose.vcf.gz chr16.dose.vcf.gz chr17.dose.vcf.gz chr18.dose.vcf.gz chr19.dose.vcf.gz chr20.dose.vcf.gz chr21.dose.vcf.gz chr22.dose.vcf.gz -Ou |
bcftools view -Ou -i "R2>0.8" |
bcftools annotate -Oz -x ID -I +"%CHROM:%POS:%REF:%ALT" -o ../Post_QC/Total.dose.R20.8.vcf.gz')
# This does not work. Please run the shell script above in terminal.




# MAF > 0.5%
# SNP call rate > 99% (geno < 1%)

system('
~/plink --vcf ./Post_QC/Total.dose.R20.8.vcf.gz --biallelic-only --allow-extra-chr 0 --autosome --maf 0.005 --geno 0.05 --hwe 1e-6 --make-bed --out ./Post_QC/Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001
')

```

### 1.8 PCA & Covariate Generation

LD Pruning 후 PCA 돌리고, 공변량 파일(Covariate file) 만드는 과정.

```r
# PCA ---------------------------------------------------------------------

system("mkdir PCA")

# Prune snps
system('
FILENAME="Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001"

~/plink --bfile ./Post_QC/${FILENAME} --indep-pairwise 1000 10 0.02 --autosome --exclude ../exclusion_regions_hg19.txt --out ./PCA/${FILENAME}_PrunedData
')


# Extract pruned SNPs and only these variants will be used for PC calculation
system('
FILENAME="Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001"

~/plink --bfile ./Post_QC/${FILENAME} --extract ./PCA/${FILENAME}_PrunedData.prune.in --make-bed --out ./PCA/${FILENAME}_PrunedData_forPCA
')


# Calculate/generate PCs based on pruned data set
system('
FILENAME="Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001"

~/plink --bfile ./PCA/${FILENAME}_PrunedData_forPCA --pca 30 --out ./PCA/${FILENAME}_PrunedData_forPCA_outPCA
')

# then use the .eigenvec file





# R
library(tidyverse)
df.pc = read.table("./PCA/Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_PrunedData_forPCA_outPCA.eigenvec", header=FALSE)

df.eigenvalue <- scan("./PCA/Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_PrunedData_forPCA_outPCA.eigenval")

df.pc <- df.pc[,-1]
# set names
names(df.pc)[1] <- "IID"
names(df.pc)[2:ncol(df.pc)] <- paste0("PC", 1:(ncol(df.pc)-1))
head(df.pc)

write.table(df.pc, file="./data/merge0422.PC", row.names=FALSE, quote=FALSE)


df.pheno = read.table("./data/merge0422.pheno", header=TRUE)
# df.pheno$FID <- paste0(df.pheno$FID, "_", df.pheno$FID)
df.pheno.pc <- left_join(df.pc, df.pheno, "IID") %>% mutate(AGE_SEX = AGE * SEX, AGE2 = AGE^2, AGE2_SEX = AGE2 * SEX) %>% select(FID, IID, AGE, SEX, AGE_SEX, AGE2, AGE2_SEX, everything())


library(magrittr)
df.fam <- read.table("./Post_QC/Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001.fam") %>%
  set_colnames(c("FID", "IID", "PID", "MID", "SEX", "PHENO"))

df.id <- read.table("./data/merge0422.ID")[,2]



df.fam <- left_join( df.fam %>% select(-SEX), df.pheno.pc %>% select(IID, SEX, SWB), "IID" ) %>% mutate(SEX = case_when(SEX==0~2, SEX==1~1, .default=0)) %>% select(FID, IID, PID, MID, SEX, PHENO)




write.table(df.pheno.pc, file="./data/merge0422.pheno_pc", row.names=FALSE, quote=FALSE)



write.table(df.fam, file="./Post_QC/Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender.fam", col.names=FALSE, row.names=FALSE, quote=FALSE)

system("cp ./Post_QC/Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001.bim ./Post_QC/Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender.bim")
system("cp ./Post_QC/Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001.bed ./Post_QC/Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender.bed")

```

### 1.9 Final GWAS Execution

최종적으로 GWAS (Linear Regression) 실행.

```r
# GWAS --------------------------------------------------------------------

library(tidyverse)
df.pheno <- read.table("./data/merge0422.pheno_pc", header=TRUE)

df.pheno %>% head(3)


system("mkdir GWAS")


# FINAL
system('
FILENAME="Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender"

~/plink --bfile ./Post_QC/${FILENAME} --autosome --snps-only just-acgt --covar ./data/merge0422.pheno_pc --covar-name AGE,SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --pheno ./data/merge0422.pheno_pc --pheno-name SWB --linear --adjust --out ./GWAS/${FILENAME}_SWB
')

# AGE,SEX,AGE_SEX,AGE2,AGE2_SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10

```

---

## Part 2. PRS Analysis

Part 1에서 만들어진 Clean Data와 Summary Statistics를 이용해서 PRS 점수 계산하고 검증하는 단계임. (이전 요청하셨던 코드 부분)

### 2.1 Setup & Helper Functions

```r
## https://github.com/FINNGEN/CS-PRS-pipeline
library(tidyverse)
library(data.table)
library(png.utils)


# 사용자 정의 함수
library(stringi)
striHelper <- function(x) stri_c(x[stri_order(x)], collapse = "")
png.str.sort <- function(string) vapply(stri_split_boundaries(string, type = "character"), striHelper, "")
"asjdkflasdf" %>% png.str.sort()

```

### 2.2 Data Loading & Harmonization

```r
# 작업 디렉토리 및 파일 설정
setwd("/Volumes/png2/LSH/merge0422-PRS")

bfile_path <- "/Volumes/png2/LSH/merge0422/FinalData"
summary_path <- "/Volumes/png2/LSH/merge0422-PRS/SummaryStatistics"
summary_header <- c("CHR", "POS", "SNP", "A1_EFFECT", "A2_NONEFFECT", "NMISS", "BETA", "SE", "PVALUE")


# 데이터 불러오기: bfile=.bim file,  sfile=SummaryStat file
bfile <- list.files(bfile_path, pattern="\\.(bim)", full.names=TRUE)
sfile <- list.files(summary_path, pattern="\\.tsv$", full.names=TRUE)

df_bim <- data.table::fread(bfile)
colnames(df_bim) <- c("bim.CHR", "bim.SNP", "bim.GD", "bim.POS", "bim.REF", "bim.ALT")

df_summary <- data.table::fread(sfile)
colnames(df_summary) <- summary_header
dim(df_summary)
# [1] 6825851        9

# .bim과 SummaryStat join하기
# > by: ("CHR", "POS") in summary  vs  ("bim.CHR", "bim.POS") in .bim
# >> 이를 위해 relationship = "many-to-many"로 설정함.
# > suffix: 중복되는 칼럼명이 있으면 SummaryStat 칼럼명은 그대로 ("") 두고 .bim의 칼럼명에 .new 추가하기 (e.g. pvalue >> pvalue.new).
df_summary2 <- df_summary %>% 
  left_join(df_bim, 
            by=c("CHR"="bim.CHR", "POS"="bim.POS"), 
            suffix=c("", ".new"), 
            relationship="many-to-many")




# SummaryStat과 .bim 간의 shared SNPs 비교
## the number of rows
df_summary2 %>% dim
# [1] 6837756      13
df_summary2 %>% filter(!is.na(bim.SNP)) %>% dim
# [1] 6294054      13
df_summary2 %>% filter(is.na(bim.SNP)) %>% dim
# [1] 543702     13

## the number of SNPs
df_summary2 %>% filter(!duplicated(SNP)) %>% dim
# [1] 6823957      13
df_summary2 %>% filter(!is.na(bim.SNP)) %>% filter(!duplicated(SNP)) %>% dim
# [1] 6280287      13
df_summary2 %>% filter(is.na(bim.SNP)) %>% filter(!duplicated(SNP)) %>% dim
# [1] 543670     13

# Shared SNPs 리스트 따로 저장해두기
df_summary2 %>% filter(!is.na(bim.SNP)) %>% arrange(PVALUE) %>% write.table(file="sumstat+bim_common.txt", quote=FALSE, row.names=F)
df_summary2 %>% filter(is.na(bim.SNP)) %>% arrange(PVALUE) %>% write.table(file="sumstat+bim_onlysumstat.txt", quote=FALSE, row.names=F)



# Multiallelic SNP의 경우에는 (CHR, POS)로 join되었더라도 다른 SNP을 가질 수 있음.
df_summary3 <- df_summary2 %>% filter(!is.na(bim.SNP))

df_summary3[,c("chr","pos","bim.allele1","bim.allele2"):=tstrsplit(bim.SNP, ":", fixed = TRUE, keep = 1:4)]
df_summary3[,c("A12","bim.allele12"):=list( paste(A1_EFFECT, A2_NONEFFECT, sep="") %>% png.str.sort(), paste(bim.allele1, bim.allele2, sep="") %>% png.str.sort() )]


df_summary3 %>% filter(A12 != bim.allele12) %>% filter(!duplicated(SNP)) %>% dim()
# [1] 11083    19
df_summary3 %>% filter(A12 == bim.allele12) %>% filter(!duplicated(SNP)) %>% dim()
# [1] 6278736       19

df_summary3 %>% filter(A12 != bim.allele12) %>% dim()
# [1] 13620    19
df_summary3 %>% filter(A12 == bim.allele12) %>% dim()
# [1] 6280434       19



df_summary3 %>% filter(A12 != bim.allele12) %>% arrange(PVALUE) %>%
  write.table("./sumstat+bim_DiffAlleles.txt", quote=FALSE, row.names=F, col.names=T)

df_summary3 %>% filter(A12 == bim.allele12) %>%
  write.table("./SummaryStatistics/GCST90104897_buildGRCh37_new.txt", quote=FALSE, row.names=F, col.names=T)


# c("CHR", "POS", "SNP", "A1_EFFECT", "A2_NONEFFECT", "NMISS", "BETA", "SE", "PVALUE")
df_summary3 %>% filter(A12 == bim.allele12) %>% filter(!duplicated(bim.SNP)) %>% 
  select(bim.SNP, SNP, CHR, POS, A1_EFFECT, A2_NONEFFECT) %>%
  write.table(file="Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender.rsID_conversion", quote=FALSE, row.names=F, col.names=F)

df_summary3 %>% filter(A12 == bim.allele12) %>% filter(!duplicated(bim.SNP)) %>% .$bim.SNP %>%
  write.table(file="Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS.snpList", quote=FALSE, row.names=FALSE, col.names=FALSE)

```

### 2.3 PRS Calculation (Clumping & Scoring)

```r
#
system("~/plink2 --bfile ./FinalData/Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender --extract ./Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS.snpList --make-bed --out Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered")



system("~/plink2 --bfile ./Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered --update-name Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender.rsID_conversion 2 1 --make-bed --out Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID")


system("~/plink2 --bfile Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID --freq --out Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID")





# --allow-extra-chr

# Remove multi-allelic SNPs
system("cut -f 2 Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID.bim | sort | uniq -d > Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID.dups")

system("~/plink --bfile Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID --exclude Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID.dups --clump ./SummaryStatistics/GCST90104897_buildGRCh37_new.txt --clump-field PVALUE --clump-kb 250 --clump-p1 1 --clump-r2 0.1 --clump-snp-field SNP --out Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump")

system("awk 'NR!=1{print $3}' Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump.clumped > Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump.snp")

system("awk '{print $3,$5}' Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump.clumped > Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump.snp_pvalue")


system('echo "5e-8 0 5e-8" > range_list
echo "5e-7 0 5e-7" >> range_list
echo "5e-6 0 5e-6" >> range_list
echo "5e-5 0 5e-5" >> range_list
echo "5e-4 0 5e-4" >> range_list
echo "5e-3 0 5e-3" >> range_list
echo "0.01 0 0.01" >> range_list
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list
echo "0.2 0 0.2" >> range_list
echo "0.3 0 0.3" >> range_list
echo "0.4 0 0.4" >> range_list
echo "0.5 0 0.5" >> range_list
echo "1 0 1" >> range_list')




system("~/plink2 --bfile Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID --score ./SummaryStatistics/GCST90104897_buildGRCh37_new.txt 3 4 7 header --q-score-range range_list Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump.snp_pvalue --extract Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump.snp --out Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump_PRS")

```

### 2.4 Threshold Optimization

```r
library(tidyverse)

p.threshold = c("5e-8","5e-7","5e-6","5e-5","5e-4","5e-3","0.01","0.05","0.1","0.2","0.3","0.4","0.5","1")[-1]


# Read in the phenotype file
pheno <- data.table::fread("./FinalData/merge0422.pheno_pc") %>%
  select(FID, IID, AGE, SEX, SWB, PC1:PC4) %>% as.data.frame

# We can then calculate the null model (model with PRS) using a linear regression
null.model <- lm(SWB~., data=pheno %>% {.[,!colnames(.)%in%c("FID","IID")]})
null.r2 <- summary(null.model)$r.squared

prs.result <- NULL
for(i in p.threshold){
  # Go through each p-value threshold
  prs <- read.table(paste0("./Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump_PRS.",i,".sscore"), header=FALSE)
  colnames(prs) <- c("FID", "IID", "CNT1", "CNT2", "SCORE")
  # Merge the prs with the phenotype matrix
  # We only want the FID, IID and PRS from the PRS file, therefore we only select the
  # relevant columns
  pheno.prs <- merge(pheno, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
  # Now perform a linear regression on Height with PRS and the covariates
  # ignoring the FID and IID from our model
  pheno.prs$SCORE <- scale(pheno.prs$SCORE)
  model <- lm((SWB)~., data=pheno.prs %>% {.[,!colnames(.)%in%c("FID","IID")]})
  # model R2 is obtained as
  model.r2 <- summary(model)$r.squared
  # R2 of PRS is simply calculated as the model R2 minus the null R2
  prs.r2 <- model.r2-null.r2
  # We can also obtain the coeffcient and p-value of association of PRS as follow
  prs.coef <- summary(model)$coeff["SCORE",]
  prs.beta <- as.numeric(prs.coef[1])
  prs.se <- as.numeric(prs.coef[2])
  prs.p <- as.numeric(prs.coef[4])
  # We can then store the results
  prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
}


prs.result
# Threshold            R2         P        BETA          SE
# 1       5e-7 0.0038947048 0.4259553 -0.06314598 0.07911224
# 2       5e-6 0.0002015533 0.8564159 -0.01473470 0.08130273
# 3       5e-5 0.0102604553 0.1955525 -0.10273400 0.07903862
# 4       5e-4 0.0044460682 0.3949131 -0.06767693 0.07933497
# 5       5e-3 0.0091189674 0.2225228  0.09682509 0.07906424
# 6       0.01 0.0103551130 0.1934928  0.10233056 0.07836372
# 7       0.05 0.0043088217 0.4023240  0.06599004 0.07858539
# 8        0.1 0.0007470564 0.7275426  0.02747658 0.07872692
# 9        0.2 0.0001091313 0.8940846  0.01048044 0.07859307
# 10       0.3 0.0027545846 0.5032768  0.05275520 0.07863701
# 11       0.4 0.0035742110 0.4456964  0.06012240 0.07864172
# 12       0.5 0.0019251601 0.5758545  0.04413411 0.07872556
# 13         1 0.0014657285 0.6254832  0.03866840 0.07906907

# Best result is:
prs.result[which.max(prs.result$R2),]
# Threshold         R2         P      BETA       SE
# 6       0.01 0.01035511 0.1934928 0.1023306 0.07836372

```

### 2.5 Visualization & Advanced Analysis

```r
# Visualization -----------------------------------------------------------

df_prs <- read.table(paste0("./Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump_PRS.",p.threshold[which.max(prs.result$R2)],".sscore"))
df_pheno <- data.table::fread("./FinalData/merge0422.pheno_pc")

df_prs %>% head
df_prs[,5] %>% summary
df_prs[,5] %>% hist


pdf(file="Figure-Scatter-PRS0.01_vs_SWB.pdf", width=7, height=5)
cbind.data.frame(PRS=df_prs[,5] %>% scale, SWB=df_pheno$SWB) %>%
  ggpubr::ggscatter(x="PRS", y="SWB", shape=18,
                    xlab="Standardized PRS", ylab="SWB",
                    add="reg.line", conf.int = TRUE,
                    add.params = list(color = "blue", fill = "gray50"), # Customize reg. line
                    cor.coef = TRUE,
                    cor.coeff.args = list(method = "pearson",
                                          label.x = max(df_prs[,5] %>% scale)*0.75,
                                          label.y = max(df_pheno$SWB)*0.95,
                                          label.sep = "\n")
  )
dev.off()


# pdf(file="Figure-Scatter-Decile_max-PRS0.01_vs_SWB.pdf", width=7, height=5)
# cbind.data.frame(PRS=df_prs[,5] %>% scale, SWB=df_pheno$SWB) %>%
#   mutate(Decile_Group_PRS=dplyr::ntile(PRS, 5)) %>% 
#   group_by(Decile_Group_PRS) %>%
#   summarise(maxSWB = max(SWB)) %>% 
#   # ggpubr::ggboxplot(x="Decile_Group_PRS", y="SWB", fill="Decile_Group_PRS")
#   ggpubr::ggscatter(x="Decile_Group_PRS", y="maxSWB", shape=18,
#                     xlab="PRS decile groups", ylab="SWB",
#                     add="reg.line", conf.int = TRUE,
#                     add.params = list(color = "blue", fill = "gray50"), # Customize reg. line
#                     cor.coef = TRUE,
#                     cor.coeff.args = list(method = "pearson",
#                                           # label.x = max(df_prs[,5] %>% scale)*0.75,
#                                           label.y = 2.7,
#                                           label.sep = "\n")
#   )
# dev.off()
# 
# 
# pdf(file="Figure-Scatter-Decile_mean-PRS0.01_vs_SWB.pdf", width=7, height=5)
# cbind.data.frame(PRS=df_prs[,5] %>% scale, SWB=df_pheno$SWB) %>%
#   mutate(Decile_Group_PRS=dplyr::ntile(PRS, 10)) %>% 
#   group_by(Decile_Group_PRS) %>%
#   summarise(meanSWB = mean(SWB)) %>% 
#   # ggpubr::ggboxplot(x="Decile_Group_PRS", y="SWB", fill="Decile_Group_PRS")
#   ggpubr::ggscatter(x="Decile_Group_PRS", y="meanSWB", shape=18,
#                     xlab="PRS decile groups", ylab="SWB",
#                     add="reg.line", conf.int = TRUE,
#                     add.params = list(color = "blue", fill = "gray50"), # Customize reg. line
#                     cor.coef = TRUE,
#                     cor.coeff.args = list(method = "pearson",
#                                           # label.x = max(df_prs[,5] %>% scale)*0.75,
#                                           # label.y = max(df_pheno$SWB)*0.95,
#                                           label.sep = "\n")
#   )
# dev.off()

cor(df_prs[,5], df_pheno$SWB)





# ggplot2 is a handy package for plotting
library(ggplot2)
# generate a pretty format for p-value output
prs.result$print.p <- round(prs.result$P, digits = 3)
prs.result$print.p[!is.na(prs.result$print.p) & prs.result$print.p == 0] <-
  format(prs.result$P[!is.na(prs.result$print.p) & prs.result$print.p == 0], digits = 2)
prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)
# Initialize ggplot, requiring the threshold as the x-axis
# (use factor so that it is uniformly distributed)
p = ggplot(data = prs.result, aes(x = factor(Threshold, levels = p.threshold), y = R2)) +
  # Specify that we want to print p-value on top of the bars
  # geom_text(
  #   aes(label = paste(print.p)),
  #   # vjust = -1.5,
  #   hjust = -0.2,
  #   angle = 90,
  #   cex = 5,
  #   parse = T
  # )  +
  # Specify the range of the plot, *1.25 to provide enough space for the p-values
  scale_y_continuous(limits = c(0, max(prs.result$R2) * 1.25)) +
  # Specify the axis labels
  xlab(expression(italic(P) - value ~ threshold)) +
  ylab(expression(paste("PRS model fit:  ", R ^ 2))) +
  # Draw a bar plot
  geom_bar(aes(fill = -log10(P)), stat = "identity") +
  # Specify the colors
  scale_fill_gradient2(
    low = "dodgerblue",
    high = "firebrick",
    mid = "dodgerblue",
    midpoint = 1e-4,
    name = bquote(atop(-log[10] ~ model ~ italic(P) - value))
  ) +
  # Some beautification of the plot
  theme_classic() + theme(
    legend.position="none",
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(caption = "Adjusted for age, sex, and the top four PCs.")
# save the plot
ggsave("Figure-Pred_R2_barplot.pdf", p, height = 5, width = 10)









# Regression --------------------------------------------------------------

df_prs <- read.table("/Volumes/png2/LSH/merge0422-PRS/Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump_PRS.0.01.sscore") %>% {scale(.[,5]) %>% as.vector}

df_pheno <- data.table::fread("./FinalData/merge0422.pheno_pc")


fit <- lm(df_pheno$SWB ~ df_pheno$SEX + df_pheno$AGE+df_prs, 
          subset=(df_pheno$AGE<30))
fit %>% summary

library(png.utils)
cbind(PRS=df_prs, AGE=df_pheno$AGE, SWB=df_pheno$SWB) %>% 
  pairs(upper.panel = NULL,
        cex = 1.5, pch = 18, # adjustcolor(4, .4),
        cex.labels = 2, font.labels = 2)



fit %>% plot
library(olsrr)
fit %>% ols_plot_cooksd_bar()

fit <- lm(df_pheno$SWB ~ df_pheno$SEX + df_pheno$AGE+df_prs, 
          subset=!(1:nrow(df_pheno)) %in% c(5,11,34,43,47,65,106,107,130,146))

fit %>% summary


library(robustbase) 
rfit <- lmrob(df_pheno$SWB ~ df_pheno$SEX + df_pheno$AGE+df_prs)
rfit %>% summary


library(broom)
df_tmp <- df_pheno %>% as.data.frame() %>% dplyr::select(SWB:PANAS_pa)


out <- NULL
out_age20 <- NULL
out_age20to35 <- NULL
for( k in 1:ncol(df_tmp) ){
  tmp <- cbind.data.frame( variable=colnames(df_tmp)[k], lm(scale(df_tmp[,k]) ~ df_pheno$SEX + df_pheno$AGE+scale(df_prs)) %>% broom::tidy() %>% filter(term=="scale(df_prs)") )
  out <- rbind.data.frame(out, tmp)
   
  tmp <- cbind.data.frame( variable=colnames(df_tmp)[k], lm(scale(df_tmp[,k]) ~ df_pheno$SEX + df_pheno$AGE+scale(df_prs), subset=df_pheno$AGE<30) %>% broom::tidy() %>% filter(term=="scale(df_prs)") )
  out_age20 <- rbind.data.frame(out_age20, tmp)
   
  tmp <- cbind.data.frame( variable=colnames(df_tmp)[k], lm(scale(df_tmp[,k]) ~ df_pheno$SEX + df_pheno$AGE+scale(df_prs), subset=df_pheno$AGE<35) %>% broom::tidy() %>% filter(term=="scale(df_prs)") )
  out_age20to35 <- rbind.data.frame(out_age20to35, tmp)
}


out
# > out
# variable          term    estimate  std.error  statistic    p.value
# 1        SWB scale(df_prs)  0.14407643 0.07730388  1.8637670 0.06415313
# 2  PWB_total scale(df_prs)  0.15189150 0.07727960  1.9654798 0.05105791
# 3     PWB_AU scale(df_prs)  0.18033386 0.07718020  2.3365301 0.02068027
# 4     PWB_EM scale(df_prs)  0.09260991 0.07778871  1.1905315 0.23556851
# 5     PWB_PG scale(df_prs)  0.08543715 0.07775540  1.0987937 0.27347838
# 6     PWB_PL scale(df_prs)  0.16866620 0.07720297  2.1847112 0.03033788
# 7     PWB_SA scale(df_prs)  0.10430445 0.07753952  1.3451779 0.18043628
# 8     PWB_PR scale(df_prs) -0.05487384 0.07822462 -0.7014907 0.48399719
# 9       SWLS scale(df_prs)  0.14407643 0.07730388  1.8637670 0.06415313
# 10  PANAS_pa scale(df_prs)  0.09571047 0.07625593  1.2551217 0.21123105

out_age20
# > out_age20
# variable          term    estimate  std.error  statistic    p.value
# 1        SWB scale(df_prs)  0.17761998 0.07819010  2.2716430 0.02444066
# 2  PWB_total scale(df_prs)  0.18389854 0.07909526  2.3250259 0.02132676
# 3     PWB_AU scale(df_prs)  0.19904436 0.07945169  2.5052252 0.01323787
# 4     PWB_EM scale(df_prs)  0.12902525 0.07917040  1.6297157 0.10512876
# 5     PWB_PG scale(df_prs)  0.11528803 0.07966300  1.4471968 0.14979823
# 6     PWB_PL scale(df_prs)  0.14685988 0.07886509  1.8621659 0.06441358
# 7     PWB_SA scale(df_prs)  0.13241067 0.07904748  1.6750775 0.09587236
# 8     PWB_PR scale(df_prs) -0.02277162 0.07976167 -0.2854958 0.77563359
# 9       SWLS scale(df_prs)  0.17761998 0.07819010  2.2716430 0.02444066
# 10  PANAS_pa scale(df_prs)  0.12240064 0.07811809  1.5668669 0.11912143

out_age20to35
# > out_age20to35
# variable          term    estimate  std.error  statistic    p.value
# 1        SWB scale(df_prs)  0.15784336 0.07869954  2.0056452 0.04655878
# 2  PWB_total scale(df_prs)  0.17589607 0.07823356  2.2483454 0.02590287
# 3     PWB_AU scale(df_prs)  0.19123022 0.07865353  2.4312986 0.01613519
# 4     PWB_EM scale(df_prs)  0.12576140 0.07816126  1.6089992 0.10956352
# 5     PWB_PG scale(df_prs)  0.10509437 0.07893882  1.3313395 0.18494776
# 6     PWB_PL scale(df_prs)  0.16222770 0.07876467  2.0596505 0.04103204
# 7     PWB_SA scale(df_prs)  0.12034853 0.07886374  1.5260312 0.12895206
# 8     PWB_PR scale(df_prs) -0.03549617 0.07943275 -0.4468708 0.65556521
# 9       SWLS scale(df_prs)  0.15784336 0.07869954  2.0056452 0.04655878
# 10  PANAS_pa scale(df_prs)  0.11024346 0.07760178  1.4206306 0.15734602

```
