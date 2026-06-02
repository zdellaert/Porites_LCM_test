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
  - [7. Two tools to identiy potential
    outliers:](#7-two-tools-to-identiy-potential-outliers)
    - [PCA](#pca)
    - [Hierarchical Clustering](#hierarchical-clustering)
  - [Heatmap of variable genes](#heatmap-of-variable-genes)
    - [Text summary](#text-summary)

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
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Etc/UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] BiocParallel_1.44.0         ggnewscale_0.5.2            genefilter_1.90.0           RColorBrewer_1.1-3         
    ##  [5] pheatmap_1.0.13             DESeq2_1.50.2               SummarizedExperiment_1.40.0 Biobase_2.70.0             
    ##  [9] MatrixGenerics_1.22.0       matrixStats_1.5.0           GenomicRanges_1.62.0        Seqinfo_1.0.0              
    ## [13] IRanges_2.44.0              S4Vectors_0.48.0            BiocGenerics_0.56.0         generics_0.1.4             
    ## [17] SeuratWrappers_0.4.0        future_1.70.0               scCustomize_3.3.0           patchwork_1.3.2            
    ## [21] Seurat_5.5.0                SeuratObject_5.4.0          sp_2.2-1                    lubridate_1.9.4            
    ## [25] forcats_1.0.0               stringr_1.6.0               dplyr_1.2.1                 purrr_1.2.1                
    ## [29] readr_2.2.0                 tidyr_1.3.2                 tibble_3.3.1                ggplot2_4.0.3              
    ## [33] tidyverse_2.0.0             Matrix_1.6-4               
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RcppAnnoy_0.0.23        splines_4.5.1           later_1.4.8             R.oo_1.27.1             polyclip_1.10-7        
    ##   [6] janitor_2.2.1           XML_3.99-0.18           fastDummies_1.7.6       lifecycle_1.0.5         globals_0.19.1         
    ##  [11] lattice_0.22-7          hdf5r_1.3.12            MASS_7.3-65             magrittr_2.0.5          plotly_4.12.0          
    ##  [16] rmarkdown_2.31          yaml_2.3.12             remotes_2.5.0           httpuv_1.6.17           otel_0.2.0             
    ##  [21] sctransform_0.4.3       spam_2.11-4             spatstat.sparse_3.2-0   reticulate_1.46.0       DBI_1.2.3              
    ##  [26] cowplot_1.2.0           pbapply_1.7-4           pkgload_1.5.2           abind_1.4-8             Rtsne_0.17             
    ##  [31] R.utils_2.13.0          GenomeInfoDbData_1.2.14 circlize_0.4.18         ggrepel_0.9.8           irlba_2.3.7            
    ##  [36] listenv_0.10.1          spatstat.utils_3.2-3    goftest_1.2-3           RSpectra_0.16-2         annotate_1.86.1        
    ##  [41] spatstat.random_3.5-0   fitdistrplus_1.2-6      parallelly_1.47.0       codetools_0.2-20        DelayedArray_0.36.0    
    ##  [46] tidyselect_1.2.1        shape_1.4.6.1           UCSC.utils_1.4.0        farver_2.1.2            spatstat.explore_3.8-1 
    ##  [51] jsonlite_2.0.0          progressr_0.19.0        ggridges_0.5.7          survival_3.8-3          systemfonts_1.3.2      
    ##  [56] tools_4.5.1             ragg_1.5.2              ica_1.0-3               Rcpp_1.1.1-1.1          glue_1.8.1             
    ##  [61] SparseArray_1.10.2      gridExtra_2.3           xfun_0.56               GenomeInfoDb_1.44.3     withr_3.0.2            
    ##  [66] BiocManager_1.30.27     fastmap_1.2.0           mcprogress_0.1.1        digest_0.6.39           rsvd_1.0.5             
    ##  [71] timechange_0.3.0        R6_2.6.1                mime_0.13               textshaping_1.0.5       ggprism_1.0.7          
    ##  [76] colorspace_2.1-2        scattermore_1.2         tensor_1.5.1            RSQLite_3.52.0          dichromat_2.0-0.1      
    ##  [81] spatstat.data_3.1-9     R.methodsS3_1.8.2       RhpcBLASctl_0.23-42     data.table_1.18.4       httr_1.4.8             
    ##  [86] htmlwidgets_1.6.4       S4Arrays_1.10.0         uwot_0.2.4              pkgconfig_2.0.3         gtable_0.3.6           
    ##  [91] rsconnect_1.4.2         blob_1.2.4              lmtest_0.9-40           S7_0.2.2                XVector_0.50.0         
    ##  [96] htmltools_0.5.9         dotCall64_1.2           scales_1.4.0            png_0.1-9               harmony_2.0.3          
    ## [101] spatstat.univar_3.2-0   snakecase_0.11.1        knitr_1.50              rstudioapi_0.17.1       tzdb_0.5.0             
    ## [106] reshape2_1.4.5          nlme_3.1-168            cachem_1.1.0            zoo_1.8-15              GlobalOptions_0.1.4    
    ## [111] KernSmooth_2.23-26      parallel_4.5.1          miniUI_0.1.2            vipor_0.4.7             AnnotationDbi_1.72.0   
    ## [116] ggrastr_1.0.2           pillar_1.11.1           grid_4.5.1              vctrs_0.7.3             RANN_2.6.2             
    ## [121] promises_1.5.0          xtable_1.8-8            cluster_2.1.8.2         beeswarm_0.4.0          paletteer_1.7.0        
    ## [126] evaluate_1.0.5          locfit_1.5-9.12         cli_3.6.5               compiler_4.5.1          crayon_1.5.3           
    ## [131] rlang_1.2.0             future.apply_1.20.2     labeling_0.4.3          rematch2_2.1.2          plyr_1.8.9             
    ## [136] ggbeeswarm_0.7.3        stringi_1.8.7           viridisLite_0.4.3       deldir_2.0-4            Biostrings_2.78.0      
    ## [141] lazyeval_0.2.3          spatstat.geom_3.8-1     RcppHNSW_0.7.0          hms_1.1.4               bit64_4.6.0-1          
    ## [146] KEGGREST_1.50.0         shiny_1.13.0            ROCR_1.0-12             memoise_2.0.1           igraph_2.3.2           
    ## [151] bit_4.6.0

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
cat("Annotations:", nrow(SwissProt), "Swissprot-annotated genes")
```

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

## 7. Two tools to identiy potential outliers:

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
save_ggplot(PCA_simple, "PCA_simple", width = 8, height = 6)
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

## Heatmap of variable genes

``` r
topVarGenes <- head(order(rowVars(vst_mat), decreasing=TRUE), 500)

png("output_RNA/analysis/plots/top500vargenes_heatmap.png")
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

### Text summary

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

    ##   PC1 variance: 84 %

    ##   PC2 variance: 16 %

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
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Etc/UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] BiocParallel_1.44.0         ggnewscale_0.5.2            genefilter_1.90.0           RColorBrewer_1.1-3         
    ##  [5] pheatmap_1.0.13             DESeq2_1.50.2               SummarizedExperiment_1.40.0 Biobase_2.70.0             
    ##  [9] MatrixGenerics_1.22.0       matrixStats_1.5.0           GenomicRanges_1.62.0        Seqinfo_1.0.0              
    ## [13] IRanges_2.44.0              S4Vectors_0.48.0            BiocGenerics_0.56.0         generics_0.1.4             
    ## [17] SeuratWrappers_0.4.0        future_1.70.0               scCustomize_3.3.0           patchwork_1.3.2            
    ## [21] Seurat_5.5.0                SeuratObject_5.4.0          sp_2.2-1                    lubridate_1.9.4            
    ## [25] forcats_1.0.0               stringr_1.6.0               dplyr_1.2.1                 purrr_1.2.1                
    ## [29] readr_2.2.0                 tidyr_1.3.2                 tibble_3.3.1                ggplot2_4.0.3              
    ## [33] tidyverse_2.0.0             Matrix_1.6-4               
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RcppAnnoy_0.0.23        splines_4.5.1           later_1.4.8             R.oo_1.27.1             polyclip_1.10-7        
    ##   [6] janitor_2.2.1           XML_3.99-0.18           fastDummies_1.7.6       lifecycle_1.0.5         globals_0.19.1         
    ##  [11] lattice_0.22-7          hdf5r_1.3.12            MASS_7.3-65             magrittr_2.0.5          plotly_4.12.0          
    ##  [16] rmarkdown_2.31          yaml_2.3.12             remotes_2.5.0           httpuv_1.6.17           otel_0.2.0             
    ##  [21] sctransform_0.4.3       spam_2.11-4             spatstat.sparse_3.2-0   reticulate_1.46.0       DBI_1.2.3              
    ##  [26] cowplot_1.2.0           pbapply_1.7-4           pkgload_1.5.2           abind_1.4-8             Rtsne_0.17             
    ##  [31] R.utils_2.13.0          GenomeInfoDbData_1.2.14 circlize_0.4.18         ggrepel_0.9.8           irlba_2.3.7            
    ##  [36] listenv_0.10.1          spatstat.utils_3.2-3    goftest_1.2-3           RSpectra_0.16-2         annotate_1.86.1        
    ##  [41] spatstat.random_3.5-0   fitdistrplus_1.2-6      parallelly_1.47.0       codetools_0.2-20        DelayedArray_0.36.0    
    ##  [46] tidyselect_1.2.1        shape_1.4.6.1           UCSC.utils_1.4.0        farver_2.1.2            spatstat.explore_3.8-1 
    ##  [51] jsonlite_2.0.0          progressr_0.19.0        ggridges_0.5.7          survival_3.8-3          systemfonts_1.3.2      
    ##  [56] tools_4.5.1             ragg_1.5.2              ica_1.0-3               Rcpp_1.1.1-1.1          glue_1.8.1             
    ##  [61] SparseArray_1.10.2      gridExtra_2.3           xfun_0.56               GenomeInfoDb_1.44.3     withr_3.0.2            
    ##  [66] BiocManager_1.30.27     fastmap_1.2.0           mcprogress_0.1.1        digest_0.6.39           rsvd_1.0.5             
    ##  [71] timechange_0.3.0        R6_2.6.1                mime_0.13               textshaping_1.0.5       ggprism_1.0.7          
    ##  [76] colorspace_2.1-2        scattermore_1.2         tensor_1.5.1            RSQLite_3.52.0          dichromat_2.0-0.1      
    ##  [81] spatstat.data_3.1-9     R.methodsS3_1.8.2       RhpcBLASctl_0.23-42     data.table_1.18.4       httr_1.4.8             
    ##  [86] htmlwidgets_1.6.4       S4Arrays_1.10.0         uwot_0.2.4              pkgconfig_2.0.3         gtable_0.3.6           
    ##  [91] rsconnect_1.4.2         blob_1.2.4              lmtest_0.9-40           S7_0.2.2                XVector_0.50.0         
    ##  [96] htmltools_0.5.9         dotCall64_1.2           scales_1.4.0            png_0.1-9               harmony_2.0.3          
    ## [101] spatstat.univar_3.2-0   snakecase_0.11.1        knitr_1.50              rstudioapi_0.17.1       tzdb_0.5.0             
    ## [106] reshape2_1.4.5          nlme_3.1-168            cachem_1.1.0            zoo_1.8-15              GlobalOptions_0.1.4    
    ## [111] KernSmooth_2.23-26      parallel_4.5.1          miniUI_0.1.2            vipor_0.4.7             AnnotationDbi_1.72.0   
    ## [116] ggrastr_1.0.2           pillar_1.11.1           grid_4.5.1              vctrs_0.7.3             RANN_2.6.2             
    ## [121] promises_1.5.0          xtable_1.8-8            cluster_2.1.8.2         beeswarm_0.4.0          paletteer_1.7.0        
    ## [126] evaluate_1.0.5          locfit_1.5-9.12         cli_3.6.5               compiler_4.5.1          crayon_1.5.3           
    ## [131] rlang_1.2.0             future.apply_1.20.2     labeling_0.4.3          rematch2_2.1.2          plyr_1.8.9             
    ## [136] ggbeeswarm_0.7.3        stringi_1.8.7           viridisLite_0.4.3       deldir_2.0-4            Biostrings_2.78.0      
    ## [141] lazyeval_0.2.3          spatstat.geom_3.8-1     RcppHNSW_0.7.0          hms_1.1.4               bit64_4.6.0-1          
    ## [146] KEGGREST_1.50.0         shiny_1.13.0            ROCR_1.0-12             memoise_2.0.1           igraph_2.3.2           
    ## [151] bit_4.6.0
