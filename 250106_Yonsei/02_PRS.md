# Microbiome/PRS Analysis Workshop Code

## 1. 라이브러리 및 사용자 정의 함수 설정

기본 패키지 로드하고, 문자열 정렬을 위한 `png.str.sort` 함수 정의하는 부분임.

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

## 2. 작업 경로 및 데이터 로드

경로 설정하고 `.bim` 파일이랑 Summary Statistics 파일 불러오는 과정임. 경로(`/Volumes/png2/...`)는 실습 환경에 맞춰서 수정해야 함.

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

```

## 3. 데이터 병합 (Harmonization)

Summary Statistics랑 Target Data(.bim)를 Chromosome이랑 Position 기준으로 합치는 부분임.

```r
# .bim과 SummaryStat join하기
# > by: ("CHR", "POS") in summary  vs  ("bim.CHR", "bim.POS") in .bim
# >> 이를 위해 relationship = "many-to-many"로 설정함.
# > suffix: 중복되는 칼럼명이 있으면 SummaryStat 칼럼명은 그대로 ("") 두고 .bim의 칼럼명에 .new 추가하기 (e.g. pvalue >> pvalue.new).
df_summary2 <- df_summary %>% 
  left_join(df_bim, 
            by=c("CHR"="bim.CHR", "POS"="bim.POS"), 
            suffix=c("", ".new"), 
            relationship="many-to-many")

```

## 4. SNP 매칭 확인 및 분류

제대로 합쳐졌는지 개수 확인(`dim`)하고, 공통된 SNP랑 Summary Stat에만 있는 SNP 따로 저장하는 코드임.

```r
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

```

## 5. Multiallelic SNP 처리 및 Allele 정렬

위치(CHR, POS)가 같아도 Allele가 다를 수 있어서 확인하는 과정임. `png.str.sort` 써서 순서 상관없이 매칭되는지 체크함.

```r
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

```

## 6. PLINK 입력 파일 생성

PLINK 돌리기 전에 rsID 변환 파일이랑 SNP 리스트 파일 만드는 작업임.

```r
# c("CHR", "POS", "SNP", "A1_EFFECT", "A2_NONEFFECT", "NMISS", "BETA", "SE", "PVALUE")
df_summary3 %>% filter(A12 == bim.allele12) %>% filter(!duplicated(bim.SNP)) %>% 
  select(bim.SNP, SNP, CHR, POS, A1_EFFECT, A2_NONEFFECT) %>%
  write.table(file="Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender.rsID_conversion", quote=FALSE, row.names=F, col.names=F)

df_summary3 %>% filter(A12 == bim.allele12) %>% filter(!duplicated(bim.SNP)) %>% .$bim.SNP %>%
  write.table(file="Total.dose_R20.8_MAF0.005_geno0.05_hwe0.000001_Gender_forPRS.snpList", quote=FALSE, row.names=FALSE, col.names=FALSE)

```

## 7. PLINK 실행 (System Calls)

R 내부에서 `system()` 명령어로 PLINK 1.9/2.0 실행하는 부분임. QC, Clumping, Scoring 다 여기서 함.

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

## 8. Threshold 검증 (R Loop)

여러 P-value threshold에 대해 PRS 계산해서  비교하는 루프임.

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

## 9. 시각화 (Visualization)

산점도 그리는 코드랑, 주석 처리된 Decile Plot 코드까지 다 포함되어 있음.

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

```

## 10. 상세 회귀분석 (Subgroup Analysis & Robust Regression)

나이대별로 쪼개서 분석하거나, 이상치 제거하고 Robust Regression 돌리는 심화 분석 코드임.

```r
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
