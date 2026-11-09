# BLCA_Nerve

Code and Quarto report for the bladder cancer (BLCA) nerve blockade project.

## Repository contents

```text
code/
  01.Cell_Type.R
  02.Immune.R
  03.Stromal.R
  04.Communication.R
  05.Proximal_Label.R
  TCGA_an_vis.R
report/
  01.Cell_Typing.qmd
  02.Immune.qmd
  03.Stromal.qmd
  04.Communication.qmd
sessionInfo
```
## Analysis scripts (`code/`)

#### `01.Cell_Type.R`

Reads the main 10x scRNA-seq data, performs QC and sample-wise DoubletFinder
filtering, and generates Seurat PCA/UMAP clusters. Assigns broad and
refined cell-type labels, reclusters immune, cancer, and stromal subsets, and
performs cancer-cell GO enrichment.

#### `02.Immune.R`

Analyzes the annotated immune subset by testing treatment-versus-control
expression changes within each immune cell type. Performs GO enrichment and
GSEA, and evaluates T-cell states, markers, and signatures.

#### `03.Stromal.R`

Analyzes treatment-associated expression changes in each annotated stromal
cell type. Performs GO enrichment, compares treatment-associated genes and
pathways, and examines CGRP-receptor-related expression.

#### `04.Communication.R`

Combines the immune, cancer, and stromal annotations and infers mouse
cell-cell communication with CellChat. Summarizes global networks and
focuses on PD-L1 and IGF signaling.

#### `05.Proximal_Label.R`

Processes the separate proximal-tissue 10x dataset using the same general QC
and sample-wise doublet-removal workflow. Assigns cell labels and
summarizes cell-type composition across conditions.

#### `TCGA_an_vis.R`

Contains independent TCGA-BLCA analysis and visualization code.

## Reports (`report/`)

The report sources summarize the generated outputs.

| Source | Scope |
|---|---|
| [`01.Cell_Typing.qmd`](report/01.Cell_Typing.qmd) | QC, clustering, and broad/refined cell-type annotation |
| [`02.Immune.qmd`](report/02.Immune.qmd) | Immune-cell differential expression, enrichment, and T-cell analysis |
| [`03.Stromal.qmd`](report/03.Stromal.qmd) | Stromal expression, differential expression, and enrichment |
| [`04.Communication.qmd`](report/04.Communication.qmd) | CellChat networks, PD-L1 signaling, and IGF signaling |

After the required private data and intermediate objects are available, reports
can be rendered with:

```bash
quarto render report/01.Cell_Typing.qmd
quarto render report/02.Immune.qmd
quarto render report/03.Stromal.qmd
quarto render report/04.Communication.qmd
```

These commands are expected to be run from the repository root.

## Software environment

The recorded environment is listed in [`sessionInfo`](sessionInfo):

```text
R version 4.5.1 (2025-06-13)
Platform: aarch64-apple-darwin20
Running under: macOS Sequoia 15.8.2

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Asia/Shanghai
tzcode source: internal

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggalluvial_0.12.5      NMF_0.28               cluster_2.1.8.1        rngtools_1.5.2        
 [5] registry_0.5-1         CellChat_2.2.0         igraph_2.2.1           ProjecTILs_3.7.0      
 [9] enrichplot_1.30.4      SignatuR_0.3.0         org.Mm.eg.db_3.22.0    AnnotationDbi_1.72.0  
[13] IRanges_2.44.0         S4Vectors_0.48.0       Biobase_2.70.0         BiocGenerics_0.56.0   
[17] generics_0.1.4         clusterProfiler_4.18.4 ggraph_2.2.2           cli_3.6.5             
[21] DoubletFinder_2.0.6    RColorBrewer_1.1-3     ggplot2_4.0.2          dplyr_1.2.0           
[25] Seurat_5.4.0           SeuratObject_5.3.0     sp_2.2-0            
```

## Questions/Comments

For any comments or questions, please feel free to submit a GitHub issue or contact me via email at liuzhenyu@pumc.edu.cn.

## License

See [`LICENSE`](LICENSE) @MIT.
