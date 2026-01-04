# Microbiome/PRS Analysis Workshop:

## CoDA Pipeline

본 워크숍에서는 마이크로바이옴 데이터의 특성(Sparsity, Compositionality)을 고려한 전처리, 다양성 분석, CoDA 회귀분석(Ternary Plot) 및 Differential Abundance Analysis를 실습함.

  
## [Section 0] 라이브러리 로드 및 환경 설정
  
  ```r
rm(list=ls())

# ==============================================================================
# Title: Compositional Data Analysis (CoDA) Pipeline for Microbiome
# Author: Prof.  (Statistics & CoDA Specialist)
# Date: 2026-01-04
# Description: 
#   Microbiome 데이터의 특성(Sparsity, Compositionality)을 고려한 
#   전처리, 다양성 분석, CoDA 회귀분석(Ternary Plot) 및 Differential Abundance Analysis
# ==============================================================================

# ------------------------------------------------------------------------------
# [Section 0] 라이브러리 로드 및 환경 설정
# ------------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")


if(FALSE){
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(c("phyloseq", "microbiome", "ANCOMBC", "ComplexHeatmap")) # 필요시 설치
  
  install.packages(
    "microViz",
    repos = c(davidbarnett = "[https://david-barnett.r-universe.dev](https://david-barnett.r-universe.dev)", getOption("repos"))
  )
  
  
  
  
}

# CRAN 패키지 로드
pacman::p_load(
  tidyverse, readxl, tibble, ggpubr, vegan, MASS, # 기본 통계 및 데이터 처리
  ggtern, ggrepel,                                # 시각화 확장 (Ternary plot 등)
  remotes, composition, Rfast, Compositional      # CoDA 관련 패키지
)

# Bioconductor 패키지 로드
library(phyloseq)
library(microbiome)
library(ComplexHeatmap)
library(ANCOMBC)
# library(circlize)


# GitHub 패키지 로드 (설치되어 있지 않다면 주석 해제 후 설치)
# remotes::install_github("david-barnett/microViz")
library(microViz)

# 시각화 테마 설정
theme_set(theme_bw())

```

  ## [Section 1] 데이터 로드 및 전처리 (Data Import & Wrangling)
  
  ```r
# ------------------------------------------------------------------------------
# [Section 1] 데이터 로드 및 전처리 (Data Import & Wrangling)
# ------------------------------------------------------------------------------
message(">>> [Step 1] 데이터 로드 및 Phyloseq 객체 생성...")

# 1-1. Raw Data 불러오기
raw_counts    <- read_excel("../data/Gutmicrobiome_genus_2019.xlsx")
raw_metadata <- read_excel("../data/Gutmicrobiome_survey2019.xlsx") %>% as.data.frame()
id_mapping    <- read_excel("../data/Gutmicrobiome_IDmatch.xlsx")

# 1-2. 샘플 ID 매핑 및 정렬
# 매핑 딕셔너리 생성 (Source ID -> Target ID)
id_dictionary <- setNames(unlist(id_mapping[,1]), unlist(id_mapping[,2]))
colnames(raw_counts)[-1] <- id_dictionary[colnames(raw_counts)[-1]]

# 샘플 ID 순서 정렬 (Numeric sorting)
sorted_indices <- order(as.numeric(colnames(raw_counts)[-1])) + 1
proc_counts     <- raw_counts[, c(1, sorted_indices)]

# 1-3. 메타데이터 필터링 및 정리
clean_metadata <- raw_metadata %>% 
  filter(ID %in% colnames(proc_counts)) %>% 
  column_to_rownames(var = "ID")

# 1-4. Taxonomy 분리 및 Matrix 변환
proc_counts_split <- proc_counts %>%
  tidyr::separate(col = Taxonomy, 
                  into = c("kingdom", "phylum", "class", "order", "family", "genus"), 
                  sep = ";", fill = "right") %>% 
  mutate(across(c(kingdom, phylum, class, order, family, genus),
                ~replace_na(., "Unclassified")))

# Taxonomy & OTU Matrix 분리
# (1~6열: Taxonomy, 7열~끝: OTU Counts)
tax_mat <- as.matrix(proc_counts_split[, c("kingdom", "phylum", "class", "order", "family","genus")])
otu_mat <- as.matrix(proc_counts_split[, 7:ncol(proc_counts_split)])
class(otu_mat) <- "numeric"

# 1-5. Phyloseq 객체 생성
ps_gut <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(clean_metadata)
)

print(ps_gut)

```

  
  ## [Tutorial] Phyloseq 기본 조작 실습 (Practice)
  
  GlobalPatterns 예제 데이터를 이용한 기본기 다지기.

```r
# -------------------------------------------------------------------------
# [준비 단계] 패키지 및 데이터 로드
# -------------------------------------------------------------------------
library(phyloseq)
library(dplyr)        # 파이프 연산자(%>%) 및 데이터 조작용
library(compositions) # CLR 변환(Center Log Ratio)용

# phyloseq 내장 예제 데이터인 'GlobalPatterns' 로드
data("GlobalPatterns")
ps <- GlobalPatterns  # 실습을 위해 객체명을 ps로 지정

# 데이터 기본 정보 확인 (Taxa, Sample 수 등)
print(ps) 

# -------------------------------------------------------------------------
# [실습 1] 데이터 필터링 및 전처리 (Pruning & Filtering)
# -------------------------------------------------------------------------

# 1. Total Count가 0인 균(Taxa) 제거
# 설명: 모든 샘플을 통틀어 한 번도 발견되지 않은 불필요한 Taxa를 제거하여 분석 속도 향상
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# (참고) filter_taxa 함수를 이용한 동일한 로직
# filter_taxa(ps, function(x) sum(x > 0) > 0, TRUE)


# 2. Read depth가 1,000 미만인 불량 샘플 제거
# 설명: 시퀀싱 깊이(Read 수)가 너무 낮은 샘플은 통계적 신뢰도가 낮으므로 제외
ps <- prune_samples(sample_sums(ps) >= 1000, ps)

# 확인용: 샘플별 Read sum 요약 통계 및 히스토그램 시각화
summary(sample_sums(ps))
hist(sample_sums(ps), main = "Histogram of Sample Read Depths", xlab = "Read Depth")


# 3. (고급) 노이즈 제거 필터링
# 설명: 적어도 3개 이상의 샘플에서, 10번(count) 이상 발견된 균만 남김
# 희소(Sparse)한 데이터가 많은 마이크로바이옴 특성상 분석 신뢰도를 높이기 위해 사용
ps_filtered <- filter_taxa(ps, function(x) sum(x > 10) > 3, TRUE)
# print(ps_filtered) # 필터링 후 남은 Taxa 수 확인 가능


# -------------------------------------------------------------------------
# [실습 2] 부분 데이터 추출 (Subsetting)
# -------------------------------------------------------------------------

# 1. 특정 Phylum만 추출하여 분석
# 설명: 'Firmicutes' 문(Phylum)에 해당하는 균들만 따로 떼어내어 분석
# (GlobalPatterns 데이터의 Rank 이름은 대문자 'Phylum'임에 유의)
ps_firmicutes <- subset_taxa(ps, Phylum == "Firmicutes")


# 2. 오염원(Contaminant) 제거
# 설명: 미토콘드리아(Mitochondria) 등 숙주 DNA나 오염으로 간주되는 분류군 제거
# (!= 연산자를 사용하여 해당되지 않는 것만 남김)
ps_clean <- subset_taxa(ps, Family != "Mitochondria")


# 3. 특정 샘플 그룹만 추출
# 설명: 전체 데이터 중 SampleType이 'Feces'(분변)인 샘플만 추출
ps_feces <- subset_samples(ps, SampleType == "Feces")


# -------------------------------------------------------------------------
# [실습 3] 데이터 변환 (Transformation)
# -------------------------------------------------------------------------

# 1. 상대 풍부도(Relative Abundance)로 변환
# 설명: 각 샘플의 Count를 전체 합으로 나누어 0~1 사이의 비율로 변환
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))


# 2. 백분율(%)로 변환
# 설명: 가독성을 위해 비율에 100을 곱해 퍼센트로 변환
ps_percent <- transform_sample_counts(ps, function(x) 100 * x / sum(x))


# 3. Log 변환 및 CLR 변환
# 설명: 데이터의 치우침(Skewness)을 줄이기 위해 Log 변환 수행. 
# 0값(Zero-count) 오류 방지를 위해 1 또는 0.5를 더해줌 (Pseudo-count)

# 3-1. Log(x+1) 변환
ps_log <- transform_sample_counts(ps, function(x) log(1 + x))

# 3-2. CLR (Centered Log Ratio) 변환 (조성 데이터 분석용)
# 설명: compositions 패키지의 clr 함수 사용 (0.5를 더해 0값 처리)
ps_clr <- transform_sample_counts(ps, function(x) compositions::clr(0.5 + x))

# 변환 결과 확인 (첫 5개 Taxa, 2개 샘플)
ps_clr %>% otu_table %>% .[1:5, 1:2]

# (검증) 수동 계산 로직과 비교 (두 번째 샘플 기준)
# 설명: ps 객체의 OTU table에서 2번째 컬럼(샘플)을 가져와 수동으로 CLR 계산
ps@otu_table[, 2] %>% 
  as.data.frame %>% 
  unlist %>% 
  {log(. + 0.5)} %>%    # Log 변환 (Pseudo-count 0.5)
  { . - mean(.)} %>%    # 기하평균(Log scale에서는 산술평균) 빼기
  head                  # 앞부분 출력


# -------------------------------------------------------------------------
# [실습 4] 분류 레벨 병합 (Agglomeration)
# -------------------------------------------------------------------------

# 1. Genus 레벨로 데이터 집계
# 설명: 같은 Genus에 속하는 하위 OTU/ASV들의 Count를 합침 (NA가 있으면 해당 Taxa 제거됨)
ps_genus <- tax_glom(ps, taxrank = "Genus")


# 2. Phylum 레벨로 집계
# 설명: 문(Phylum) 수준에서 전체적인 구성을 볼 때(예: Barplot) 주로 사용
ps_phylum <- tax_glom(ps, taxrank = "Phylum")

```

  
  ## [Section 2] 데이터 필터링 (Quality Control)
  
  ```r
# ------------------------------------------------------------------------------
# [Section 2] 데이터 필터링 (Quality Control)
# ------------------------------------------------------------------------------
message(">>> [Step 2] Sparsity 필터링 및 QC 수행...")

# 필터링 파라미터 설정
min_depth <- 1000        # 샘플 당 최소 Read 수
min_count <- 10          # Taxa 당 최소 Count
min_sample_prop <- 0.05  # 최소 출현 샘플 비율 (5%)
min_prevalence <- nsamples(ps_gut) * min_sample_prop

# 2-1. 오염원 및 미분류 제거, Read Depth 필터링
ps_clean <- ps_gut %>%
  subset_taxa(family != "Mitochondria" & order != "Chloroplast") %>%
  subset_taxa(phylum != "Unclassified") %>% 
  subset_taxa(!is.na(phylum)) %>%
  prune_samples(sample_sums(.) >= min_depth, .)

# 2-2. Prevalence Filtering (Low abundant taxa removal)
ps_clean <- prune_taxa(taxa_sums(ps_clean) >= min_count, ps_clean)
ps_clean <- prune_taxa(taxa_sums(ps_clean) >= min_prevalence, ps_clean)

# (선택) Phylum/Genus 레벨 병합이 필요한 경우 사용 가능한 추가 파이프라인
ps_clean2 <- ps_gut %>%
  subset_taxa(phylum != "Unclassified") %>% 
  subset_taxa(family != "Mitochondria" & order != "Chloroplast") %>%
  prune_samples(sample_sums(.) >= min_depth, .) %>% 
  prune_taxa(taxa_sums(.) >= min_count, .) %>% 
  prune_taxa(taxa_sums(.) >= min_prevalence, .) %>% 
  tax_glom(taxrank = "genus")

# 최종 데이터 확인
print(ps_clean)
summary(sample_sums(ps_clean))

```


  
  ## [Section 3] 다양성 분석 (Diversity Analysis)
  
  ```r
# ------------------------------------------------------------------------------
# [Section 3] 다양성 분석 (Diversity Analysis)
# ------------------------------------------------------------------------------
message(">>> [Step 3] Alpha & Beta Diversity 분석...")

# 3-1. Alpha Diversity (Shannon Index)
alpha_meta <- estimate_richness(ps_clean, measures = c("Shannon")) %>%
  cbind.data.frame(sample_data(ps_clean))

p_alpha <- ggviolin(alpha_meta, 
                    x = "PER_GENDER_M2", 
                    y = "Shannon",
                    fill = "PER_GENDER_M2", 
                    palette = "jco", 
                    add = "boxplot",
                    add.params = list(fill = "white")) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  stat_compare_means(comparisons = list(c("1", "2"))) + # 1:Male, 2:Female
  labs(title = "Alpha Diversity (Shannon)", x = "Group")

print(p_alpha)

# 3-2. Composition Barplot (Phylum Level)
p_bar <- ps_clean %>%
  comp_barplot(
    tax_level = "phylum",
    n_taxa = 10,
    tax_order = sum,
    other_name = "Other Phyla",
    palette = distinct_palette(10, pal = "kelly"),
    bar_width = 0.8
  ) +
  facet_grid(~PER_GENDER_M2, scales = "free_x", space = "free") +
  labs(title = "Microbial Composition (Phylum)", x = "Sample ID", y = "Relative Abundance") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
        strip.background = element_rect(fill = "#2C3E50"),
        strip.text = element_text(color = "white", face = "bold"))

print(p_bar)

# 3-3. Beta Diversity (PCA Biplot with CLR)
# CoDA 접근법: Aitchison Distance (CLR transformation + Euclidean)
p_beta <- ps_clean %>%
  tax_transform("clr", rank = "genus") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(
    color = "PER_GENDER_M2", shape = "PER_GENDER_M2",
    size = 3, alpha = 0.8,
    plot_taxa = 1:5,
    tax_lab_style = list(size = 4, fontface = "bold.italic")
  ) +
  stat_ellipse(aes(color = PER_GENDER_M2), type = "t") +
  labs(title = "PCA Biplot (Aitchison Distance)") + 
  theme(legend.position="top")

print(p_beta)

```


  
  ## [Section 4] Visualization: CoDA-based Heatmap
  
  ```r
# ------------------------------------------------------------------------------
# [Section 4] Visualization: CoDA-based Heatmap
# ------------------------------------------------------------------------------
message(">>> [Step 4] CLR 변환 기반 Heatmap 시각화...")

# 4-1. 상위 20개 Genus 추출 및 전처리
top20_taxa <- names(sort(taxa_sums(ps_clean), decreasing = TRUE))[1:20]
ps_heatmap <- prune_taxa(top20_taxa, ps_clean)

# OTU Matrix 추출 및 Transpose (Row: Sample, Col: Taxa)
mat_count <- as(otu_table(ps_heatmap), "matrix")
if(!taxa_are_rows(ps_heatmap)) mat_count <- t(mat_count)

# 4-2. CoDA Transformation (Pseudo-count + CLR + Z-score)
# CLR 변환 (Geometric Mean 보정)
mat_clr <- apply(mat_count + 0.5, 2, function(x) log(x) - mean(log(x)))
# Visualization을 위한 Z-score Scaling (Taxa 기준)
mat_plot <- t(scale(t(mat_clr)))

# 4-3. Heatmap Annotation & Drawing
meta_df <- data.frame(sample_data(ps_heatmap))
group_colors <- c("F" = "#E41A1C", "M" = "#377EB8") # 라벨에 맞게 조정 필요

ha <- ComplexHeatmap::HeatmapAnnotation(
  Group = factor(meta_df$PER_GENDER_M2, levels=1:2, labels=c("F","M")),
  col = list(Group = group_colors),
  annotation_name_gp = gpar(fontface = "bold")
)

ht <- ComplexHeatmap::Heatmap(
  mat_plot,
  name = "Rel. Abundance\n(Z-score of CLR)",
  top_annotation = ha,
  cluster_rows = TRUE, cluster_columns = TRUE,
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  show_row_names = TRUE, show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_title = "Clustering Analysis for Top 20 Taxa (CLR Transformed)"
)

draw(ht, merge_legend = TRUE)

```


  
  ## [Section 5] Advanced: Compositional Regression & Ternary Plot
  
  ```r
# ------------------------------------------------------------------------------
# [Section 5] Advanced: Compositional Regression & Ternary Plot
# ------------------------------------------------------------------------------
message(">>> [Step 5] CoDA 회귀분석 및 Ternary Plot...")

# 5-0. 사용자 정의 함수 정의 (회귀 및 시각화용)

# (1) CLR 기반 회귀 분석 엔진 (Customized Function)
comp.reg.new <- function (y, x, type = "classical", xnew = NULL, yb = NULL){
  # y: Composition (CLR transformed inside if yb is null)
  # x: Predictor (e.g., Age)
  
  if (is.null(yb)) {
    z <- as.matrix(as.data.frame(compositions::clr(y)))
  } else z <- yb
  
  if (type == "lmfit") {
    runtime <- proc.time()
    x <- model.matrix(z ~ ., data.frame(x))
    be <- solve(crossprod(x), crossprod(x, z))
    if (!is.null(xnew)) {
      xnew <- model.matrix(~., data.frame(xnew))
      est1 <- xnew %*% be
    } else est1 <- NULL
    seb <- NULL
    runtime <- proc.time() - runtime
  }
  else if (type == "classical") {
    runtime <- proc.time()
    # Compositional 패키지의 multivreg 이용
    mod <- Compositional::multivreg(z, x, plot = FALSE, xnew = xnew)
    res <- mod$suma
    di <- ncol(z)
    be <- seb <- matrix(nrow = NCOL(x) + 1, ncol = di)
    for (i in 1:di) {
      be[, i] <- res[, 1, i]
      seb[, i] <- res[, 2, i]
    }
    rownames(seb) <- rownames(be) <- rownames(res[, , 1])
    colnames(seb) <- colnames(be) <- colnames(mod$fitted)
    est1 <- mod$est
    runtime <- proc.time() - runtime
  }
  else if (type == "spatial") {
    mod <- Compositional::spatmed.reg(z, x, xnew = xnew)
    be <- mod$be
    seb <- mod$seb
    est1 <- mod$est
    runtime <- mod$runtime
  }
  
  # Back-transformation (CLR -> Simplex)
  est <- NULL
  if (!is.null(est1)) {
    est2 <- exp(est1)
    est <- est2/Rfast::rowsums(est2) # Closure operation
  }
  
  list(runtime = runtime, be = be, seb = seb, est = est)
}

# (2) CoDA Regression Wrapper Function (Pipeline 통합용)
run_coda_regression <- function(composition, covariate, target_var="target_var", scale_arrow = 1.0) {
  
  # 필요한 패키지 로드
  require(ggtern); require(MASS); require(dplyr); require(tidyr); require(phyloseq)
  
  x_comp <- composition
  y_cov  <- covariate
  taxa_names <- colnames(x_comp)
  
  # Zero Replacement & CLR Setup
  X_mat <- as.matrix(x_comp)
  min_val <- min(X_mat[X_mat > 0])
  X_mat[X_mat == 0] <- 0.5 * min_val # Simple replacement
  
  calc_clr <- function(x) { 
    if(is.vector(x)) x <- t(x)
    log_x <- log(x)
    return(log_x - rowMeans(log_x)) 
  }
  X_clr <- calc_clr(X_mat)
  
  # Regression Execution (using comp.reg.new)
  yseq <- seq(min(y_cov, na.rm=TRUE), max(y_cov, na.rm=TRUE), length.out=100)
  reg_res <- comp.reg.new(X_clr, y_cov, xnew=yseq, yb=X_clr)
  est <- reg_res$est # Predicted compositions on the grid
  coef_df <- reg_res$be # Coefficients
  
  # Data Frame Preparation for Plotting
  df <- data.frame(x_comp, y=y_cov)
  colnames(df) <- c(paste0("P", 1:ncol(x_comp)), "y")
  
  grid_df <- data.frame(est, y_pred=yseq)
  colnames(grid_df) <- c(paste0("P", 1:ncol(x_comp)), "y_pred")
  
  # Plotting Limits
  common_limits <- range(c(df$y, grid_df$y_pred), na.rm = TRUE)
  mid_point_val <- mean(df$y, na.rm = TRUE)
  
  # Ternary Plot with Trajectory
  p <- ggtern(data=df, aes(x = P1, y = P2, z = P3)) +
    
    geom_path(data = grid_df, aes(color = y_pred), color = "black", linewidth = 2.5, lineend = "round", linejoin = "round") +
    geom_path(data = grid_df, aes(color = y_pred), size = 1) + # Regression Curve
    
    geom_point(aes(fill = y), shape = 21, color = "black", size = 2.5, stroke = 0.3) +
    labs(x = taxa_names[1], y = taxa_names[2], z = taxa_names[3],
         title = paste("CoDA Regression: Effect on", target_var)) +
    scale_fill_gradient2(
      name = target_var, low = "#4575B4", mid = "#FFFFBF", high = "#D73027", 
      midpoint = mid_point_val, limits = common_limits
    ) +
    scale_color_gradient2(
      name = target_var, low = "#4575B4", mid = "#FFFFBF", high = "#D73027", 
      midpoint = mid_point_val, limits = common_limits
    ) +
    theme_bw() + theme_showarrows() +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  
  return(list(plot = p, coefficients = coef_df, model_info = "GLS on CLR coordinates"))
}


# 5-1. 데이터 준비 (상위 3개 Phylum 추출)
top3_phyla <- names(sort(taxa_sums(ps_clean), decreasing = TRUE))[1:3]
ps_top3 <- prune_taxa(top3_phyla, ps_clean) %>% 
  phyloseq::transform_sample_counts(function(x) ifelse(x==0, x+0.5, x)) %>% 
  phyloseq::transform_sample_counts(function(x) x/sum(x))

df_tern <- psmelt(ps_top3) %>%
  dplyr::select(Sample, OTU, Abundance, PER_GENDER_M2, AGE = PER_AGE_M2) %>%
  tidyr::pivot_wider(names_from = OTU, values_from = Abundance)

# 5-2. 단순 시각화 (기본 Ternary Plot)
p_tern_basic <- ggtern(data = df_tern, 
                       aes_string(x = top3_phyla[1], y = top3_phyla[2], z = top3_phyla[3])) +
  geom_point(aes(fill = AGE), shape = 21, size = 3, color = "black") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Simple Ternary Plot") +
  theme_rgbw()
print(p_tern_basic)

# 5-3. CoDA Regression 실행 (Numerical & Categorical)
comp_matrix <- df_tern[, top3_phyla] %>% as.matrix

# (A) 연속형 변수 (AGE) 회귀분석
res_num <- run_coda_regression(composition=comp_matrix, covariate=df_tern$AGE, target_var="AGE")
print(res_num$plot)

# (B) 범주형 변수 (Gender) 회귀분석 (Numeric factor로 변환 후 실행)
res_cat <- run_coda_regression(composition=comp_matrix, as.numeric(as.factor(df_tern$PER_GENDER_M2)), target_var="GENDER_Level")
print(res_cat$plot)

```


  
  ## [Section 6] Differential Abundance Analysis (ANCOM-BC2)
  
  ```r
# ------------------------------------------------------------------------------
# [Section 6] Differential Abundance Analysis (ANCOM-BC2)
# ------------------------------------------------------------------------------
message(">>> [Step 6] ANCOM-BC2 분석 수행...")

# 6-1. 그룹 변수 Factor 변환
sample_data(ps_clean)$PER_GENDER_M2 <- factor(sample_data(ps_clean)$PER_GENDER_M2, 
                                              levels = c("1", "2"), labels = c("Male", "Female"))

# 6-2. ANCOM-BC2 실행
# (Bias-corrected, Structural Zero detection included)
out_ancom <- ancombc2(
  data = ps_clean,
  fix_formula = "PER_GENDER_M2",
  p_adj_method = "holm",
  group = "PER_GENDER_M2",
  struc_zero = TRUE,
  neg_lb = TRUE,
  verbose = TRUE
)


# 6-3. 결과 정리 및 Waterfall Plot
res_df <- out_ancom$res

# 컬럼명 자동 매칭 (Reference: Male -> Female 비교)
target_col <- colnames(res_df)[grep("lfc_", colnames(res_df))][2] # 첫 번째 LFC 컬럼
se_col     <- colnames(res_df)[grep("se_", colnames(res_df))][2]
q_col      <- colnames(res_df)[grep("q_", colnames(res_df))][2]



viz_df <- res_df %>%
  mutate(nlog10p = -log10(get(q_col))) %>% 
  dplyr::select(taxon, all_of(c(target_col, se_col, q_col)), nlog10p) %>%
  rename(lfc = !!target_col, se = !!se_col, q_val = !!q_col) %>%
  mutate(
    Status = case_when(
      q_val < 0.05 & lfc > 0 ~ "Enriched in Female",
      q_val < 0.05 & lfc < 0 ~ "Enriched in Male",
      TRUE ~ "Non-significant"
    )
  )





# 상위 유의미한 Taxa 라벨링을 위해 데이터 분리
top_taxa <- viz_df %>%
  filter(Status != "Non-significant") %>%
  arrange(q_val) %>%
  head(10) # 상위 10개만 라벨 표시

ggplot(viz_df, aes(x = lfc, y = nlog10p, color = Status)) +
  geom_point(alpha = 0.7, size = 2.5) +
  
  # 가이드라인 (P-value 0.05 & LFC 1)
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  
  # 상위 Taxa 이름 표시 (ggrepel 패키지 추천)
  ggrepel::geom_text_repel(data = top_taxa, aes(label = taxon), size = 3, show.legend = FALSE) +
  
  scale_color_manual(values = c("Enriched in Female" = "#377EB8", 
                                "Enriched in Male" = "#E41A1C", 
                                "Non-significant" = "grey80")) +
  labs(
    title = "ANCOM-BC2 Volcano Plot",
    x = "Bias-corrected Log Fold Change",
    y = "-Log10 (Adjusted P-value)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "top"
  )




# Significant Taxa Plotting
sig_taxa <- viz_df %>%
  filter(Status != "Non-significant") %>%
  arrange(lfc) %>%
  mutate(taxon = factor(taxon, levels = taxon),
         ci_lower = lfc - 1.96 * se,
         ci_upper = lfc + 1.96 * se)

p_ancom <- ggplot(sig_taxa, aes(x = lfc, y = taxon, fill = Status)) +
  geom_col(width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.2) +
  geom_vline(xintercept = 0, linetype = "solid") +
  scale_fill_manual(values = c("Enriched in Male" = "#377EB8", "Enriched in Female" = "#E41A1C")) +
  labs(title = "Differential Abundance (ANCOM-BC2)",
       subtitle = "Male (Ref) vs Female",
       x = "Bias-corrected Log Fold Change (95% CI)") +
  theme_classic() +
  theme(axis.text.y = element_text(face = "italic"))

print(p_ancom)

```
