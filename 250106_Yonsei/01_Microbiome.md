# Microbiome Data Analysis Workshop: CoDA Pipeline

**Date:** 2026-01-05  
**Author:** Prof. Kipoong Kim (Dept. of Statistics, Changwon National University)  
**Topic:** Compositional Data Analysis (CoDA) & Differential Abundance Analysis

---

## 📖 Introduction
본 워크숍에서는 마이크로바이옴 데이터가 가지는 고유한 특성인 **희소성(Sparsity)**과 **조성성(Compositionality)**을 고려한 통계적 분석 파이프라인을 실습합니다. 

주요 내용은 다음과 같습니다:
1. **Phyloseq**을 이용한 데이터 핸들링
2. **CLR (Centered Log Ratio)** 변환의 이해
3. **Diversity Analysis** (Alpha/Beta w/ Aitchison Distance)
4. **CoDA Regression** 및 Ternary Plot 시각화
5. **ANCOM-BC2**를 이용한 차별 풍부도 분석

---

## [Section 0] 라이브러리 로드 및 환경 설정

분석에 필요한 CRAN 및 Bioconductor 패키지를 로드합니다. `microViz`와 같은 GitHub 패키지는 별도 설치가 필요할 수 있습니다.

```r
rm(list=ls())

# 1. 패키지 설치 확인 및 로드 함수
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

# 2. CRAN 패키지 로드
pacman::p_load(
  tidyverse, readxl, tibble, ggpubr, vegan, MASS, # 기본 통계 및 데이터 처리
  ggtern, ggrepel,                                # 시각화 확장 (Ternary plot 등)
  remotes, composition, Rfast, Compositional      # CoDA 관련 패키지
)

# 3. Bioconductor 패키지 로드
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# 필요시 설치: BiocManager::install(c("phyloseq", "microbiome", "ANCOMBC", "ComplexHeatmap"))

library(phyloseq)
library(microbiome)
library(ComplexHeatmap)
library(ANCOMBC)

# 4. GitHub 패키지 로드 (microViz)
# remotes::install_github("david-barnett/microViz")
library(microViz)

# 시각화 테마 설정
theme_set(theme_bw())

```

---

## [Section 1] 데이터 로드 및 전처리 (Data Import)

실습 데이터(`Gutmicrobiome_genus_2019.xlsx` 등)를 불러와 `phyloseq` 객체를 생성하는 과정입니다.

```r
message(">>> [Step 1] 데이터 로드 및 Phyloseq 객체 생성...")

# 1-1. Raw Data 불러오기
# (경로는 실제 데이터 위치에 맞게 수정 필요)
raw_counts   <- read_excel("../data/Gutmicrobiome_genus_2019.xlsx")
raw_metadata <- read_excel("../data/Gutmicrobiome_survey2019.xlsx") %>% as.data.frame()
id_mapping   <- read_excel("../data/Gutmicrobiome_IDmatch.xlsx")

# 1-2. 샘플 ID 매핑 및 정렬
id_dictionary <- setNames(unlist(id_mapping[,1]), unlist(id_mapping[,2]))
colnames(raw_counts)[-1] <- id_dictionary[colnames(raw_counts)[-1]]

sorted_indices <- order(as.numeric(colnames(raw_counts)[-1])) + 1
proc_counts    <- raw_counts[, c(1, sorted_indices)]

# 1-3. 메타데이터 필터링
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

tax_mat <- as.matrix(proc_counts_split[, c("kingdom", "phylum", "class", "order", "family","genus")])
otu_mat <- as.matrix(proc_counts_split[, 7:ncol(proc_counts_split)])
class(otu_mat) <- "numeric"

# 1-5. Phyloseq 객체 생성 (Main Object)
ps_gut <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(clean_metadata)
)

print(ps_gut)

```

---

## [Tutorial] Phyloseq 기본 조작 실습 (Practice)

본 섹션은 `GlobalPatterns` 예제 데이터를 사용하여 필터링 및 변환 기초를 익히는 단계입니다.

```r
# 데이터 로드
data("GlobalPatterns")
ps <- GlobalPatterns

# [실습 1] 데이터 필터링 (Pruning & Filtering)
# 1. Total Count > 0 인 Taxa만 유지
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# 2. Read depth 1,000 이상인 샘플만 유지
ps <- prune_samples(sample_sums(ps) >= 1000, ps)

# 3. (Optional) 노이즈 제거: 3개 이상 샘플에서 10번 이상 발견된 균만 유지
ps_filtered <- filter_taxa(ps, function(x) sum(x > 10) > 3, TRUE)

# [실습 2] 부분 데이터 추출 (Subsetting)
# 1. 특정 Phylum 추출
ps_firmicutes <- subset_taxa(ps, Phylum == "Firmicutes")
# 2. 오염원 제거 (Mitochondria)
ps_clean_ex <- subset_taxa(ps, Family != "Mitochondria")
# 3. 특정 샘플 타입 추출 (Feces)
ps_feces <- subset_samples(ps, SampleType == "Feces")

# [실습 3] 데이터 변환 (Transformation)
# 1. 상대 풍부도 (Relative Abundance)
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# 2. CLR (Centered Log Ratio) 변환 (Zero handling: +0.5)
# CoDA 분석의 핵심 전처리 단계입니다.
ps_clr <- transform_sample_counts(ps, function(x) compositions::clr(0.5 + x))

# [실습 4] 분류 레벨 병합 (Agglomeration)
ps_genus <- tax_glom(ps, taxrank = "Genus")

```

---

## [Section 2] Quality Control (Filtering)

실제 데이터(`ps_gut`)에 대한 QC를 수행합니다.

```r
message(">>> [Step 2] Sparsity 필터링 및 QC 수행...")

# 파라미터 설정
min_depth <- 1000        # 샘플 당 최소 Read 수
min_count <- 10          # Taxa 당 최소 Count
min_sample_prop <- 0.05  # 최소 5% 샘플에서 발견되어야 함
min_prevalence <- nsamples(ps_gut) * min_sample_prop

# 필터링 파이프라인
ps_clean <- ps_gut %>%
  subset_taxa(family != "Mitochondria" & order != "Chloroplast") %>%
  subset_taxa(phylum != "Unclassified") %>% 
  subset_taxa(!is.na(phylum)) %>%
  prune_samples(sample_sums(.) >= min_depth, .) %>%
  prune_taxa(taxa_sums(.) >= min_count, .) %>%
  prune_taxa(taxa_sums(.) >= min_prevalence, .)

print(ps_clean)

```

---

## [Section 3] Diversity Analysis

### 3-1. Alpha Diversity (Shannon)

```r
alpha_meta <- estimate_richness(ps_clean, measures = c("Shannon")) %>%
  cbind.data.frame(sample_data(ps_clean))

# Visualization
ggviolin(alpha_meta, x = "PER_GENDER_M2", y = "Shannon",
         fill = "PER_GENDER_M2", palette = "jco", add = "boxplot") +
  geom_jitter(width = 0.1, alpha = 0.5) +
  stat_compare_means(comparisons = list(c("1", "2"))) +
  labs(title = "Alpha Diversity (Shannon)", x = "Group")

```

### 3-2. Beta Diversity (PCA Biplot with CLR)

Aitchison Distance(CLR 변환 후 유클리드 거리)를 기반으로 한 PCA 분석입니다.

```r
ps_clean %>%
  tax_transform("clr", rank = "genus") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(
    color = "PER_GENDER_M2", shape = "PER_GENDER_M2",
    size = 3, alpha = 0.8,
    plot_taxa = 1:5,
    tax_lab_style = list(size = 4, fontface = "bold.italic")
  ) +
  stat_ellipse(aes(color = PER_GENDER_M2), type = "t") +
  labs(title = "PCA Biplot (Aitchison Distance)")

```

---

## [Section 4] Visualization: CoDA-based Heatmap

상위 20개 Genus에 대해 CLR 변환 후 Z-score를 계산하여 군집화를 수행합니다.

```r
# 상위 20개 추출
top20_taxa <- names(sort(taxa_sums(ps_clean), decreasing = TRUE))[1:20]
ps_heatmap <- prune_taxa(top20_taxa, ps_clean)

# Matrix 변환 및 CLR -> Z-score
mat_count <- as(otu_table(ps_heatmap), "matrix")
if(!taxa_are_rows(ps_heatmap)) mat_count <- t(mat_count)
mat_clr <- apply(mat_count + 0.5, 2, function(x) log(x) - mean(log(x)))
mat_plot <- t(scale(t(mat_clr)))

# Heatmap 그리기
ha <- HeatmapAnnotation(
  Group = factor(sample_data(ps_heatmap)$PER_GENDER_M2, levels=1:2, labels=c("F","M")),
  col = list(Group = c("F" = "#E41A1C", "M" = "#377EB8"))
)

Heatmap(
  mat_plot,
  name = "Rel. Abundance\n(Z-score of CLR)",
  top_annotation = ha,
  cluster_rows = TRUE, cluster_columns = TRUE,
  col = circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
)

```

---

## [Section 5] CoDA Regression & Ternary Plot

**Ternary Plot**은 3개의 구성요소(예: Top 3 Phyla)의 비율 변화를 시각화하는 강력한 도구입니다.

### 5-1. Custom Function 정의

```r
# 회귀분석 및 시각화를 위한 래퍼 함수
run_coda_regression <- function(composition, covariate, target_var="target_var") {
  require(ggtern); require(MASS); require(dplyr); require(phyloseq)
  
  # (Custom regression logic omitted for brevity - see full script)
  # ... [중략: 위 코드의 comp.reg.new 및 run_coda_regression 함수 포함] ...
}

```

*(참고: 실제 R 스크립트 실행 시에는 `comp.reg.new` 함수와 `run_coda_regression` 전체 코드를 포함해야 합니다.)*

### 5-2. 실행 및 시각화

```r
# 데이터 준비 (Top 3 Phyla)
top3_phyla <- names(sort(taxa_sums(ps_clean), decreasing = TRUE))[1:3]
ps_top3 <- prune_taxa(top3_phyla, ps_clean) %>% 
  transform_sample_counts(function(x) ifelse(x==0, x+0.5, x)) %>% 
  transform_sample_counts(function(x) x/sum(x))

df_tern <- psmelt(ps_top3) %>%
  dplyr::select(Sample, OTU, Abundance, PER_GENDER_M2, AGE = PER_AGE_M2) %>%
  pivot_wider(names_from = OTU, values_from = Abundance)

# Ternary Plot (AGE에 따른 변화)
comp_matrix <- df_tern[, top3_phyla] %>% as.matrix
res_num <- run_coda_regression(composition=comp_matrix, covariate=df_tern$AGE, target_var="AGE")
print(res_num$plot)

```

---

## [Section 6] Differential Abundance (ANCOM-BC2)

기존의 t-test 등이 가지는 Compositional Bias를 교정한 **ANCOM-BC2** 분석입니다.

```r
# 1. ANCOM-BC2 실행
sample_data(ps_clean)$PER_GENDER_M2 <- factor(sample_data(ps_clean)$PER_GENDER_M2, 
                                              levels = c("1", "2"), labels = c("Male", "Female"))
out_ancom <- ancombc2(
  data = ps_clean, fix_formula = "PER_GENDER_M2",
  p_adj_method = "holm", group = "PER_GENDER_M2",
  struc_zero = TRUE, neg_lb = TRUE
)

# 2. 결과 시각화 (Volcano Plot)
res_df <- out_ancom$res
# ... [시각화 코드 생략: Volcano Plot & Bar Plot 코드] ...

```

---

**End of Workshop Tutorial**
