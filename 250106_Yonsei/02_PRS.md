# Polygenic Risk Score (PRS) Analysis Pipeline

본 튜토리얼은 GWAS Summary Statistic과 Target Genotype 데이터를 사용하여 Polygenic Risk Score(PRS)를 계산하고, 표현형(Phenotype)과의 연관성을 분석하는 전체 파이프라인을 다룹니다.

## Prerequisites

이 분석을 수행하기 위해서는 다음의 소프트웨어와 R 패키지가 필요합니다.

* **Software:** [PLINK 1.9](https://www.cog-genomics.org/plink/) & [PLINK 2.0](https://www.cog-genomics.org/plink/2.0/)
* **R Packages:** `tidyverse`, `data.table`, `stringi`, `ggplot2`, `ggpubr`, `robustbase`, `broom`, `olsrr`
* **Custom Package:** `png.utils` (User defined)


## 1. Setup & Helper Functions

필요한 라이브러리를 로드하고 문자열 처리를 위한 사용자 정의 함수를 설정합니다.

```r
library(tidyverse)
library(data.table)
library(stringi)

# 사용자 정의 함수 (Custom Helper Functions)
# png.utils 패키지가 없는 경우 아래 함수들을 직접 정의하여 사용
striHelper <- function(x) stri_c(x[stri_order(x)], collapse = "")
png.str.sort <- function(string) vapply(stri_split_boundaries(string, type = "character"), striHelper, "")

# Test function
"asjdkflasdf" %>% png.str.sort()

```

## 2. Data Loading & Harmonization

작업 경로를 설정하고, Target 데이터(.bim)와 Base 데이터(Summary Statistics)를 불러옵니다. 이후 두 데이터 간의 SNP를 매칭하고 Allele를 정렬하여 Harmonization을 수행합니다.

> **Note:** 실습 환경에 맞춰 `setwd` 및 파일 경로(`bfile_path`, `summary_path`)를 수정해 주세요.

```r
# 작업 디렉토리 및 파일 설정
setwd("/Volumes/png2/LSH/merge0422-PRS") # 본인의 경로로 수정 필요

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

# .bim과 SummaryStat join하기 (Based on CHR, POS)
df_summary2 <- df_summary %>% 
  left_join(df_bim, 
            by=c("CHR"="bim.CHR", "POS"="bim.POS"), 
            suffix=c("", ".new"), 
            relationship="many-to-many")

# Shared SNPs 확인 및 저장
df_summary2 %>% filter(!is.na(bim.SNP)) %>% arrange(PVALUE) %>% 
  write.table(file="sumstat+bim_common.txt", quote=FALSE, row.names=F)

# Multiallelic SNP 처리 및 Allele 정렬 확인
df_summary3 <- df_summary2 %>% filter(!is.na(bim.SNP))
df_summary3[,c("chr","pos","bim.allele1","bim.allele2"):=tstrsplit(bim.SNP, ":", fixed = TRUE, keep = 1:4)]
df_summary3[,c("A12","bim.allele12"):=list( paste(A1_EFFECT, A2_NONEFFECT, sep="") %>% png.str.sort(), paste(bim.allele1, bim.allele2, sep="") %>% png.str.sort() )]

# Mismatched Alleles 확인 및 저장
df_summary3 %>% filter(A12 != bim.allele12) %>% arrange(PVALUE) %>%
  write.table("./sumstat+bim_DiffAlleles.txt", quote=FALSE, row.names=F, col.names=T)

# Matched Alleles 저장 (QC 통과한 SNP)
df_summary3 %>% filter(A12 == bim.allele12) %>%
  write.table("./SummaryStatistics/GCST90104897_buildGRCh37_new.txt", quote=FALSE, row.names=F, col.names=T)

# rsID Conversion 파일 생성
df_summary3 %>% filter(A12 == bim.allele12) %>% filter(!duplicated(bim.SNP)) %>% 
  select(bim.SNP, SNP, CHR, POS, A1_EFFECT, A2_NONEFFECT) %>%
  write.table(file="Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender.rsID_conversion", quote=FALSE, row.names=F, col.names=F)

# PRS용 SNP List 생성
df_summary3 %>% filter(A12 == bim.allele12) %>% filter(!duplicated(bim.SNP)) %>% .$bim.SNP %>%
  write.table(file="Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS.snpList", quote=FALSE, row.names=FALSE, col.names=FALSE)

```

## 3. PLINK Processing (Clumping & Scoring)

R의 `system()` 함수를 이용하여 PLINK 명령어를 실행합니다. 이 단계에서는 SNP Filtering, Clumping, 그리고 최종적인 Score 계산이 수행됩니다.

```r
# 1. Extract SNPs and Make BED
system("~/plink2 --bfile ./FinalData/Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender --extract ./Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS.snpList --make-bed --out Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered")

# 2. Update SNP IDs to rsIDs
system("~/plink2 --bfile ./Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered --update-name Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender.rsID_conversion 2 1 --make-bed --out Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID")

# 3. Frequency Check
system("~/plink2 --bfile Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID --freq --out Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID")

# 4. Remove Multi-allelic Duplicates
system("cut -f 2 Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID.bim | sort | uniq -d > Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID.dups")

# 5. Clumping (LD Pruning)
# Parameters: kb=250, p1=1, r2=0.1
system("~/plink --bfile Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID --exclude Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID.dups --clump ./SummaryStatistics/GCST90104897_buildGRCh37_new.txt --clump-field PVALUE --clump-kb 250 --clump-p1 1 --clump-r2 0.1 --clump-snp-field SNP --out Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump")

# 6. Extract Clumped SNPs
system("awk 'NR!=1{print $3}' Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump.clumped > Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump.snp")
system("awk '{print $3,$5}' Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump.clumped > Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump.snp_pvalue")

# 7. Create Range List for Thresholding
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

# 8. Calculate PRS Score
system("~/plink2 --bfile Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID --score ./SummaryStatistics/GCST90104897_buildGRCh37_new.txt 3 4 7 header --q-score-range range_list Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump.snp_pvalue --extract Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump.snp --out Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump_PRS")

```

## 4. Threshold Selection & Validation

다양한 P-value threshold에 대해 계산된 PRS를 로드하고, Phenotype(Subjective Well-Being, SWB)에 대한 설명력()을 비교하여 최적의 Threshold를 선정합니다.

```r
library(tidyverse)

p.threshold = c("5e-8","5e-7","5e-6","5e-5","5e-4","5e-3","0.01","0.05","0.1","0.2","0.3","0.4","0.5","1")[-1]

# Phenotype Load (Covariates included)
pheno <- data.table::fread("./FinalData/merge0422.pheno_pc") %>% 
  select(FID, IID, AGE, SEX, SWB, PC1:PC4) %>% as.data.frame

# Null Model (Only Covariates)
null.model <- lm(SWB~., data=pheno %>% {.[,!colnames(.)%in%c("FID","IID")]})
null.r2 <- summary(null.model)$r.squared

prs.result <- NULL
for(i in p.threshold){
  # Read PRS file for specific threshold
  prs <- read.table(paste0("./Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump_PRS.",i,".sscore"), header=FALSE)
  colnames(prs) <- c("FID", "IID", "CNT1", "CNT2", "SCORE")
  
  # Merge with Phenotype
  pheno.prs <- merge(pheno, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
  
  # Scale PRS
  pheno.prs$SCORE <- scale(pheno.prs$SCORE)
  
  # Full Model (Covariates + PRS)
  model <- lm((SWB)~., data=pheno.prs %>% {.[,!colnames(.)%in%c("FID","IID")]})
  
  # Calculate Partial R2 for PRS
  model.r2 <- summary(model)$r.squared
  prs.r2 <- model.r2 - null.r2
  
  # Store Statistics
  prs.coef <- summary(model)$coeff["SCORE",]
  prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs.r2, P=as.numeric(prs.coef[4]), BETA=as.numeric(prs.coef[1]), SE=as.numeric(prs.coef[2])))
}

# Best Threshold 확인
print(prs.result[which.max(prs.result$R2),])

```

## 5. Visualization

결과를 시각화합니다.

1. **Scatter Plot:** PRS(0.01 threshold)와 SWB 간의 관계
2. **Bar Plot:** P-value threshold 별 Model Fit () 비교

```r
library(ggplot2)
library(ggpubr)

# Best Threshold에 해당하는 PRS 데이터 로드
best_thresh <- p.threshold[which.max(prs.result$R2)]
df_prs <- read.table(paste0("./Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump_PRS.", best_thresh, ".sscore"))
df_pheno <- data.table::fread("./FinalData/merge0422.pheno_pc")

# 1. Scatter Plot
pdf(file="Figure-Scatter-PRS0.01_vs_SWB.pdf", width=7, height=5)
cbind.data.frame(PRS=df_prs[,5] %>% scale, SWB=df_pheno$SWB) %>%
  ggscatter(x="PRS", y="SWB", shape=18,
            xlab="Standardized PRS", ylab="SWB",
            add="reg.line", conf.int = TRUE,
            add.params = list(color = "blue", fill = "gray50"),
            cor.coef = TRUE,
            cor.coeff.args = list(method = "pearson", label.sep = "\n"))
dev.off()

# 2. R2 Bar Plot
p = ggplot(data = prs.result, aes(x = factor(Threshold, levels = p.threshold), y = R2)) +
  scale_y_continuous(limits = c(0, max(prs.result$R2) * 1.25)) +
  xlab(expression(italic(P) - value ~ threshold)) +
  ylab(expression(paste("PRS model fit:  ", R ^ 2))) +
  geom_bar(aes(fill = -log10(P)), stat = "identity") +
  scale_fill_gradient2(
    low = "dodgerblue", high = "firebrick", mid = "dodgerblue",
    midpoint = 1e-4,
    name = bquote(atop(-log[10] ~ model ~ italic(P) - value))
  ) +
  theme_classic() + 
  theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(caption = "Adjusted for age, sex, and the top four PCs.")

ggsave("Figure-Pred_R2_barplot.pdf", p, height = 5, width = 10)

```

## 6. Detailed Regression Analysis

연령별 Subgroup 분석 및 Robust Regression 등을 통해 결과를 검증합니다.

```r
library(broom)
library(robustbase)
library(olsrr)

# 데이터 준비
df_prs_vec <- read.table(paste0("./Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS_snpFiltered_rsID_clump_PRS.", best_thresh, ".sscore")) %>% 
  {scale(.[,5]) %>% as.vector}

# Subgroup Analysis Function
run_subgroup_analysis <- function(data, prs_vec, age_limit=NULL) {
  subset_condition <- if(is.null(age_limit)) rep(TRUE, nrow(data)) else data$AGE < age_limit
  
  df_tmp <- data %>% as.data.frame() %>% dplyr::select(SWB:PANAS_pa)
  out_df <- NULL
  
  for(k in 1:ncol(df_tmp)){
    res <- lm(scale(df_tmp[,k]) ~ data$SEX + data$AGE + scale(prs_vec), subset=subset_condition) %>% 
      broom::tidy() %>% 
      filter(term=="scale(prs_vec)")
    out_df <- rbind.data.frame(out_df, cbind.data.frame(variable=colnames(df_tmp)[k], res))
  }
  return(out_df)
}

# Run Analyses
out_all <- run_subgroup_analysis(df_pheno, df_prs_vec)
out_age20 <- run_subgroup_analysis(df_pheno, df_prs_vec, age_limit=30)
out_age35 <- run_subgroup_analysis(df_pheno, df_prs_vec, age_limit=35)

# 결과 확인
print(out_age20)

# Robust Regression Example
rfit <- lmrob(df_pheno$SWB ~ df_pheno$SEX + df_pheno$AGE + df_prs_vec)
summary(rfit)

```

---

### 💡 교수님, 추가로 확인하실 사항입니다.

1. **`png.utils` 패키지 의존성:** 코드 초반부에 `library(png.utils)`가 있는데, 이는 교수님께서 직접 만드신 패키지나 로컬 함수 모음으로 보입니다. GitHub에 올리실 때는 해당 패키지가 같이 업로드되어 있거나, 제가 작성해 드린 코드의 **Step 1**처럼 `striHelper`, `png.str.sort` 함수를 스크립트 내에 직접 정의해 주는 것이 실습생들에게 오류를 줄이는 방법일 것 같습니다. (위 마크다운에는 직접 정의하는 방식으로 넣어두었습니다.)
2. **경로(Path) 수정:** `/Volumes/png2/...` 와 같은 절대 경로는 워크숍 참가자들의 환경과 다를 수 있습니다. 실습용 데이터를 GitHub 레포지토리 내 `data/` 폴더 등에 넣고 상대 경로(`.` 또는 `./data`)를 사용하도록 안내하시면 더 좋을 것 같습니다.
3. **PLINK 실행 권한:** Mac/Linux 환경에 따라 `~/plink2` 경로가 다를 수 있으니, 워크숍 전에 환경 변수 설정($PATH)에 대해 간단히 언급해주시면 진행이 매끄러울 것입니다.

성공적인 워크숍 되시길 응원합니다!
