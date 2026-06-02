RNA-seq Preprocessing and Normalization
================
Zoe Dellaert
2026-06-02

- [Preproccessing of bulk RNA-seq
  data](#preproccessing-of-bulk-rna-seq-data)
  - [0. Setup species-specific
    parameters](#0-setup-species-specific-parameters)
  - [1. Read in raw count data](#1-read-in-raw-count-data)
  - [2. Extract metadata from sample
    names](#2-extract-metadata-from-sample-names)
  - [3. Remove outliers, if
    identified](#3-remove-outliers-if-identified)
  - [4. pOverA filtering to reduce
    dataset](#4-povera-filtering-to-reduce-dataset)
    - [Note to self: maybe replace this with treatment-specific
      filtering. To get genes expressed only at one timepoint in one
      treatment](#note-to-self-maybe-replace-this-with-treatment-specific-filtering-to-get-genes-expressed-only-at-one-timepoint-in-one-treatment)
  - [5. Create DESeq object and run
    DESeq2](#5-create-deseq-object-and-run-deseq2)
  - [6. VST-Transforming count data for
    visualization](#6-vst-transforming-count-data-for-visualization)
  - [7. Visualize sample-sample
    relationships](#7-visualize-sample-sample-relationships)
    - [PCA](#pca)
    - [Hierarchical Clustering](#hierarchical-clustering)
    - [Heatmap of variable genes](#heatmap-of-variable-genes)
  - [Preprocessing Summary](#preprocessing-summary)
- [DE Analysis](#de-analysis)
  - [1. Extract results for bulk vs. OralGastro
    contrast](#1-extract-results-for-bulk-vs-oralgastro-contrast)
    - [MA Plots with Log2 Fold Change Transform
      Comparisons](#ma-plots-with-log2-fold-change-transform-comparisons)
  - [2. Extract results for adjusted p-value \< 0.05 with LFC transform
    of choice (or
    none)](#2-extract-results-for-adjusted-p-value--005-with-lfc-transform-of-choice-or-none)
    - [Join with annotation data](#join-with-annotation-data)
    - [Save csvs](#save-csvs)
  - [3. Heatmap of differentially expressed genes, with Swissprot
    annotation](#3-heatmap-of-differentially-expressed-genes-with-swissprot-annotation)
  - [Appendix](#appendix)

# Preproccessing of bulk RNA-seq data

``` r
# set up file paths so that Rmd outputs can be viewed using github markdown
knitr::opts_knit$set(base.dir = normalizePath("output_RNA/reports/"), base.url = "./")

knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE,fig.width = 10, fig.height = 8,
                      fig.path = "DE_Analysis_files/figure-gfm/")

#load packages
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(genefilter)
library(ggnewscale)
library(BiocParallel)

# set number of cores to use for parallel DESeq2 processing
register(MulticoreParam(workers = 18))

sessionInfo() #provides list of loaded packages and version of R
```

    ## R version 4.5.1 (2025-06-13)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.1 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Etc/UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] BiocParallel_1.44.0         ggnewscale_0.5.2            genefilter_1.90.0           RColorBrewer_1.1-3          pheatmap_1.0.13            
    ##  [6] DESeq2_1.50.2               SummarizedExperiment_1.40.0 Biobase_2.70.0              MatrixGenerics_1.22.0       matrixStats_1.5.0          
    ## [11] GenomicRanges_1.62.0        Seqinfo_1.0.0               IRanges_2.44.0              S4Vectors_0.48.0            BiocGenerics_0.56.0        
    ## [16] generics_0.1.4              lubridate_1.9.4             forcats_1.0.0               stringr_1.6.0               dplyr_1.2.1                
    ## [21] purrr_1.2.1                 readr_2.2.0                 tidyr_1.3.2                 tibble_3.3.1                ggplot2_4.0.3              
    ## [26] tidyverse_2.0.0            
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] DBI_1.2.3               rlang_1.2.0             magrittr_2.0.5          otel_0.2.0              compiler_4.5.1          RSQLite_3.52.0         
    ##  [7] png_0.1-9               systemfonts_1.3.2       vctrs_0.7.3             pkgconfig_2.0.3         crayon_1.5.3            fastmap_1.2.0          
    ## [13] XVector_0.50.0          labeling_0.4.3          rmarkdown_2.31          tzdb_0.5.0              UCSC.utils_1.4.0        ragg_1.5.2             
    ## [19] bit_4.6.0               xfun_0.56               cachem_1.1.0            GenomeInfoDb_1.44.3     jsonlite_2.0.0          blob_1.2.4             
    ## [25] DelayedArray_0.36.0     irlba_2.3.7             parallel_4.5.1          R6_2.6.1                stringi_1.8.7           SQUAREM_2026.1         
    ## [31] numDeriv_2016.8-1.1     Rcpp_1.1.1-1.1          knitr_1.51              Matrix_1.6-4            splines_4.5.1           timechange_0.3.0       
    ## [37] tidyselect_1.2.1        rstudioapi_0.17.1       dichromat_2.0-0.1       abind_1.4-8             yaml_2.3.12             codetools_0.2-20       
    ## [43] plyr_1.8.9              lattice_0.22-7          withr_3.0.2             KEGGREST_1.50.0         S7_0.2.2                coda_0.19-4.1          
    ## [49] evaluate_1.0.5          survival_3.8-3          Biostrings_2.78.0       pillar_1.11.1           rsconnect_1.4.2         invgamma_1.2           
    ## [55] emdbook_1.3.14          truncnorm_1.0-9         hms_1.1.4               scales_1.4.0            ashr_2.2-63             xtable_1.8-8           
    ## [61] glue_1.8.1              apeglm_1.30.0           tools_4.5.1             annotate_1.86.1         locfit_1.5-9.12         mvtnorm_1.3-3          
    ## [67] XML_3.99-0.18           grid_4.5.1              bbmle_1.0.25.1          bdsmatrix_1.3-7         AnnotationDbi_1.72.0    GenomeInfoDbData_1.2.14
    ## [73] cli_3.6.6               textshaping_1.0.5       mixsqp_0.3-54           S4Arrays_1.10.0         gtable_0.3.6            digest_0.6.39          
    ## [79] SparseArray_1.10.2      farver_2.1.2            memoise_2.0.1           htmltools_0.5.9         lifecycle_1.0.5         httr_1.4.8             
    ## [85] MASS_7.3-65             bit64_4.6.0-1

``` r
save_ggplot <- function(plot, filename, width = 10, height = 7, units = "in", dpi = 300,bg = "white") {
  png_path <- file.path(outdir_plots, paste0(filename, ".png"))
  pdf_dir <- file.path(outdir_plots, "pdf_figs")
  pdf_path <- file.path(pdf_dir, paste0(filename, ".pdf"))
  
  # Ensure the pdf_figs directory exists
  if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
  
  # Save plots
  ggsave(filename = png_path, plot = plot, width = width, height = height, units = units, dpi = dpi,bg = bg)
  ggsave(filename = pdf_path, plot = plot, width = width, height = height, units = units, dpi = dpi,bg = bg)
}
```

## 0. Setup species-specific parameters

``` r
# set up necessary output directories if they don't exist
outdir <- file.path("output_RNA/analysis")
outdir_plots <- file.path(outdir,"plots")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (!dir.exists(outdir_plots)) dir.create(outdir_plots, recursive = TRUE)

reportdir <- file.path("output_RNA/reports/DE_Analysis_files/figure-gfm/")
if (!dir.exists(reportdir)) dir.create(reportdir, recursive = TRUE)
```

## 1. Read in raw count data

``` r
# load in data
counts_raw <- read.csv(file.path("output_RNA/count_matrices/POR_Pcomp_gene_count_matrix.csv"), row.names = 1)

# make list of samples 
samples <- colnames(counts_raw)
cat("Raw counts:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples")
```

    ## Raw counts: 44130 genes x 3 samples

``` r
# read in SwissProt annotation
SwissProt <- read.delim(file.path("../HI_genome_annotations/annotation/Porites_compressa_HIv1_Swissprot_GO.tsv"))
cat("\nAnnotations:", nrow(SwissProt), "Swissprot-annotated genes")
```

    ## 
    ## Annotations: 22929 Swissprot-annotated genes

## 2. Extract metadata from sample names

``` r
# remove _S* from sample names
colnames(counts_raw) <- gsub("\\_S\\d\\d","",colnames(counts_raw))
samples <- colnames(counts_raw)

# create metadata dataframe from sample names
meta <- data.frame(
  sample = samples, 
  slide = c("POR_3_Slide1","POR_4_Slide3","POR_R72_C2_1"),
  tissue = c("bulk","bulk","Oral_Gastro")
)

# add rownames
rownames(meta) <- meta$sample

# make time and treatment factors
meta$tissue <- factor(meta$tissue)

# save metadata
write.csv(meta, paste0("output_RNA/RNA_seq_metadata.csv"))

# reorder count matrix to be in order of metadata table (should be already but just in case)
counts_raw <- counts_raw[, meta$sample]
```

## 3. Remove outliers, if identified

Nope!

## 4. pOverA filtering to reduce dataset

### Note to self: maybe replace this with treatment-specific filtering. To get genes expressed only at one timepoint in one treatment

``` r
# Keep genes expressed at 10+ counts in at least 33% of samples (aka one sample)

ffun<-filterfun(pOverA(0.33,10))
counts_filt_poa <- genefilter((counts_raw), ffun) #apply filter

filtered_counts <- counts_raw[counts_filt_poa,] #keep only rows that passed filter

paste0("Number of genes after filtering: ", sum(counts_filt_poa))
```

    ## [1] "Number of genes after filtering: 20338"

``` r
paste0("% of genes kept: ", round(100*(sum(counts_filt_poa)/nrow(counts_raw)),digits=2),"%")
```

    ## [1] "% of genes kept: 46.09%"

``` r
write.csv(filtered_counts, file = file.path(outdir, "filtered_counts.csv"))
```

## 5. Create DESeq object and run DESeq2

``` r
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = meta,
                              design= ~ tissue)

dds <- DESeq(dds, parallel = TRUE)

# Estimate size factors to determine if we can use VST
SF.dds <- estimateSizeFactors(dds) 
print(sort(sizeFactors(SF.dds))) #View size factors
```

    ##   POR_13   POR_61   POR_34 
    ## 0.760436 1.090223 1.338733

``` r
# if all are less than 4 we can use the VST transformation
all(sizeFactors(SF.dds)) < 4
```

    ## [1] TRUE

## 6. VST-Transforming count data for visualization

``` r
vst <- vst(dds, blind=FALSE)

#save the vst transformation
vst_mat <- assay(vst)
write.csv(vst_mat, file = file.path(outdir, "vst_expression_matrix.csv"))
```

## 7. Visualize sample-sample relationships

### PCA

``` r
pcaData <- plotPCA(vst, intgroup=c("tissue"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

PCA_simple <- ggplot(data = pcaData, aes(x=PC1, y=PC2, color=tissue, shape=sample)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  labs(color = "Tissue", shape = "Sample") +
  coord_fixed() + theme_bw() + ggtitle("PCA of VST-transformed counts")

print(PCA_simple)
```

![](./DE_Analysis_files/figure-gfm/pca-1.png)<!-- -->

``` r
save_ggplot(PCA_simple, "PCA", width = 8, height = 6)
```

``` r
pcaData <- plotPCA(vst, intgroup=c("tissue"), returnData=TRUE, ntop = nrow(vst))
percentVar <- round(100 * attr(pcaData, "percentVar"))

PCA_simple <- ggplot(data = pcaData, aes(x=PC1, y=PC2, color=tissue, shape=sample)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  labs(color = "Tissue", shape = "Sample") +
  coord_fixed() + theme_bw() + ggtitle("PCA of VST-transformed counts")

print(PCA_simple)
```

![](./DE_Analysis_files/figure-gfm/pca-allgene-1.png)<!-- -->

``` r
save_ggplot(PCA_simple, "PCA_allGenes", width = 8, height = 6)
```

### Hierarchical Clustering

``` r
sampleTree <- hclust(dist(t(vst_mat)), method = "average")

par(mar = c(8, 4, 2, 2))
plot(sampleTree, 
     xlab = "", sub = "", cex = 0.7)
abline(h = 100, col = "red", lty = 2)
```

![](./DE_Analysis_files/figure-gfm/cluster-1.png)<!-- -->

### Heatmap of variable genes

``` r
topVarGenes <- head(order(rowVars(vst_mat), decreasing=TRUE), 500)

png("output_RNA/analysis/plots/top500vargenes_heatmap.png", width = 2000, height = 2400, res = 300)
pheatmap(vst_mat[topVarGenes, ], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, cutree_cols = 2,
         annotation_col= meta %>% select(tissue))
```

![](./DE_Analysis_files/figure-gfm/heatmap-1.png)<!-- -->

``` r
dev.off()
```

    ## png 
    ##   3

## Preprocessing Summary

    ## Preprocessing Summary:

    ## Input

    ## ----------------------------------------

    ##   Initial genes: 44130

    ## Filtering

    ## ----------------------------------------

    ##   Low-expression genes removed: 23792

    ##   pOverA filter: >= 10 counts in >= 33 % of samples

    ## Output

    ## ----------------------------------------

    ##   Final genes: 20338

    ##   Final samples: 3

    ##   Output directory: output_RNA/analysis

    ## QC Notes

    ## ----------------------------------------

    ##   Size factors range: 0.76 - 1.34

    ##   VST appropriate: Yes

    ##   PC1 variance: 64 %

    ##   PC2 variance: 36 %

# DE Analysis

## 1. Extract results for bulk vs. OralGastro contrast

``` r
resultsNames(dds) 
```

    ## [1] "Intercept"                  "tissue_Oral_Gastro_vs_bulk"

``` r
res <- results(dds, name="tissue_Oral_Gastro_vs_bulk")
```

### MA Plots with Log2 Fold Change Transform Comparisons

``` r
resNorm <- lfcShrink(dds, coef="tissue_Oral_Gastro_vs_bulk", type="normal")
resAsh <- lfcShrink(dds, coef="tissue_Oral_Gastro_vs_bulk", type="ashr")
resLFC <- lfcShrink(dds, coef="tissue_Oral_Gastro_vs_bulk", res=res, type = "apeglm")

par(mfrow=c(1,4), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-20,20)
plotMA(res, xlim=xlim, ylim=ylim, main="no LFC transform")
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
```

![](./DE_Analysis_files/figure-gfm/res_transform-1.png)<!-- -->

## 2. Extract results for adjusted p-value \< 0.05 with LFC transform of choice (or none)

``` r
#res <- resLFC #resAsh

resOrdered <- res[order(res$pvalue),] # save differentially expressed genes

DE_05 <- as.data.frame(resOrdered) %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)
DE_05_Up <- DE_05 %>% filter(log2FoldChange > 0) #Higher in Oral Gastro
DE_05_Down <- DE_05 %>% filter(log2FoldChange < 0) #Lower in Oral Gastro

nrow(DE_05)
```

    ## [1] 1478

``` r
nrow(DE_05_Up) #Higher in Oral Gastro
```

    ## [1] 304

``` r
nrow(DE_05_Down) #Lower in Oral Gastro
```

    ## [1] 1174

### Join with annotation data

``` r
DE_05$query <- rownames(DE_05)
resOrdered$query <- rownames(resOrdered)

DESeq_SwissProt <- as.data.frame(resOrdered) %>% left_join(SwissProt) %>% select(query,everything()) 
DE_05_SwissProt <- DESeq_SwissProt %>% filter(query %in% DE_05$query)
```

### Save csvs

``` r
write.csv(DESeq_SwissProt, 
          file = file.path(outdir, "DESeq_results.csv"))

write.csv(DE_05_SwissProt, 
          file = file.path(outdir, "DEG_05.csv"))
```

## 3. Heatmap of differentially expressed genes, with Swissprot annotation

``` r
DE_05_SwissProt$short_name <- ifelse(nchar(DE_05_SwissProt$ProteinNames) > 50, 
                            paste0(substr(DE_05_SwissProt$ProteinNames, 1, 47), "..."), 
                            DE_05_SwissProt$ProteinNames)

gene_labels <- DE_05_SwissProt %>% 
  select(query,short_name) %>%
  mutate_all(~ ifelse(is.na(.), "", .)) #replace NAs with "" for labelling purposes

#view most significantly differentially expressed genes in order by p-value with labels
topDEGenes <- order(res$padj)[1:50]

png("output_RNA/analysis/plots/top50_DE_heatmap_swissprot.png", width = 2000, height = 2400, res = 300)
pheatmap(vst_mat[topDEGenes, ], 
         cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, cutree_cols = 2,
         annotation_col=(meta%>% select(tissue)),
         labels_row = gene_labels[match(rownames(res)[topDEGenes],(gene_labels$query)),2], fontsize_row = 6)
```

![](./DE_Analysis_files/figure-gfm/heatmap-swissprot-1.png)<!-- -->

``` r
dev.off()
```

    ## png 
    ##   3

``` r
png("output_RNA/analysis/plots/top50_DE_ordered_heatmap_swissprot.png", width = 2000, height = 2400, res = 300)
pheatmap(vst_mat[topDEGenes, ], 
         cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, cutree_cols = 2,
         annotation_col=(meta%>% select(tissue)),
         labels_row = gene_labels[match(rownames(res)[topDEGenes],(gene_labels$query)),2], fontsize_row = 6)
dev.off()
```

    ## png 
    ##   3

## Appendix

To knit: rmarkdown::render(“DE_Analysis.Rmd”, output_dir =
“output_RNA/reports/”)

``` r
sessionInfo()
```

    ## R version 4.5.1 (2025-06-13)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.1 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Etc/UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] BiocParallel_1.44.0         ggnewscale_0.5.2            genefilter_1.90.0           RColorBrewer_1.1-3          pheatmap_1.0.13            
    ##  [6] DESeq2_1.50.2               SummarizedExperiment_1.40.0 Biobase_2.70.0              MatrixGenerics_1.22.0       matrixStats_1.5.0          
    ## [11] GenomicRanges_1.62.0        Seqinfo_1.0.0               IRanges_2.44.0              S4Vectors_0.48.0            BiocGenerics_0.56.0        
    ## [16] generics_0.1.4              lubridate_1.9.4             forcats_1.0.0               stringr_1.6.0               dplyr_1.2.1                
    ## [21] purrr_1.2.1                 readr_2.2.0                 tidyr_1.3.2                 tibble_3.3.1                ggplot2_4.0.3              
    ## [26] tidyverse_2.0.0            
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] DBI_1.2.3               rlang_1.2.0             magrittr_2.0.5          otel_0.2.0              compiler_4.5.1          RSQLite_3.52.0         
    ##  [7] png_0.1-9               systemfonts_1.3.2       vctrs_0.7.3             pkgconfig_2.0.3         crayon_1.5.3            fastmap_1.2.0          
    ## [13] XVector_0.50.0          labeling_0.4.3          rmarkdown_2.31          tzdb_0.5.0              UCSC.utils_1.4.0        ragg_1.5.2             
    ## [19] bit_4.6.0               xfun_0.56               cachem_1.1.0            GenomeInfoDb_1.44.3     jsonlite_2.0.0          blob_1.2.4             
    ## [25] DelayedArray_0.36.0     irlba_2.3.7             parallel_4.5.1          R6_2.6.1                stringi_1.8.7           SQUAREM_2026.1         
    ## [31] numDeriv_2016.8-1.1     Rcpp_1.1.1-1.1          knitr_1.51              Matrix_1.6-4            splines_4.5.1           timechange_0.3.0       
    ## [37] tidyselect_1.2.1        rstudioapi_0.17.1       dichromat_2.0-0.1       abind_1.4-8             yaml_2.3.12             codetools_0.2-20       
    ## [43] plyr_1.8.9              lattice_0.22-7          withr_3.0.2             KEGGREST_1.50.0         S7_0.2.2                coda_0.19-4.1          
    ## [49] evaluate_1.0.5          survival_3.8-3          Biostrings_2.78.0       pillar_1.11.1           rsconnect_1.4.2         invgamma_1.2           
    ## [55] emdbook_1.3.14          truncnorm_1.0-9         hms_1.1.4               scales_1.4.0            ashr_2.2-63             xtable_1.8-8           
    ## [61] glue_1.8.1              apeglm_1.30.0           tools_4.5.1             annotate_1.86.1         locfit_1.5-9.12         mvtnorm_1.3-3          
    ## [67] XML_3.99-0.18           grid_4.5.1              bbmle_1.0.25.1          bdsmatrix_1.3-7         AnnotationDbi_1.72.0    GenomeInfoDbData_1.2.14
    ## [73] cli_3.6.6               textshaping_1.0.5       mixsqp_0.3-54           S4Arrays_1.10.0         gtable_0.3.6            digest_0.6.39          
    ## [79] SparseArray_1.10.2      farver_2.1.2            memoise_2.0.1           htmltools_0.5.9         lifecycle_1.0.5         httr_1.4.8             
    ## [85] MASS_7.3-65             bit64_4.6.0-1
