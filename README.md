# CASCADE

**CASCADE** , or **C**ancer **A**ggressiveness via **S**ingle-**C**ell **A**nalysis **D**uring **E**volution, is a bioinformatics tool that quantitatively analyzes liver tumors from single-cell or bulk transcriptomic data. It distills the tumor landscape down into 3 metrics:  

- Lineage score: representation of the malignant cell community. Ranges from 0 to 1, where a 0 indicates a hepatocyte- or myeloid cell-like tumor lineage and a 1 indicates a cholangiocyte- or mesenchymal-like tumor lineage. 
- Ecological score: representation of the tumor microenvironment. Ranges from 0 to 1, where a 0 indicates a restrictive tumor microenvironment and a 1 indicates a lenient one.
- Tumor score: the product of the lineage and ecological score.

The intersection of the lineage and ecological scores separates samples into 4 quadrants based on whether they are greater than or less than 0.5 along each metric. These quadrants are termed as follows:
- A1: lineage score < 0.5 and ecological score < 0.5
- A2: lineage score < 0.5 and ecological score >= 0.5
- B1: lineage score >= 0.5 and ecological score < 0.5
- B2: lineage score >= 0.5 and ecological score >= 0.5

The tumor score is an alternative to the lineage and ecological scores that distills this quadrant system down to a single number.

## Citation

TBD

## Install

Install CASCADE directly from GitHub:

```r
library(devtools)
install_github("MaLab621/CASCADE")
library(CASCADE)
```

## Usage

CASCADE works with either single-cell or bulk transcriptomic data. For single-cell data, input either a single Seurat object containing malignancy information in the metadata or separate Seurat objects for malignant and nonmalignant cells. Both forms of input also require metadata columns for sample identity and cell type. For bulk data, input a normalized matrix with genes on the rows and samples on the columns.

Running with a single Seurat object:
```r
obj <- CreateSeuratObject(counts = gene_counts)
obj$malignancy_label <- c("malignant", "non-malignant", "non-malignant")  # your label vector would be much longer than this
obj$sample_name <- c("s1", "s1", "s2")                                    # your label vector would be much longer than this
obj$cell_type <- c("tumor cells", "B cells", "CD8-GZMH-effector T cells") # your label vector would be much longer than this

cascade_df <- sc_CASCADE(total_obj = obj, malignancy_col = "malignancy_label", 
                         ident_col = "sample_name", subtype_col = "cell_type")
```

Running with separate malignant and non-malignant objects:
```r
malignant_obj <- CreateSeuratObject(counts = malignant_gene_counts)
malignant_obj$sample_name <- c("s1", "s2", "s3")                        # your label vector would be much longer than this

nonmalignant_obj <- CreateSeuratObject(counts = nonmalignant_gene_counts)
nonmalignant_obj$sample_name <- c("s1", "s4", "s5")                     # your label vector would be much longer than this
nonmalignant_obj$cell_type <- c("MAIT", "B cells", "CD4-KLRB1-T cells") # your label vector would be much longer than this

cascade_df <- sc_CASCADE(tumor_obj = malignant_obj, nontumor_obj = nonmalignant_obj,
                         ident_col = "sample_name", subtype_col = "cell_type", min_n_cells = 10)
```

Running with bulk data:
```r
cascade_df <- bulk_CASCADE(bulk_expr = bulk_data)
```

Sample visualization plot:
```r
ggplot(cascade_df, aes(x = lineage_score, y = eco_score)) +
  geom_point() +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  xlim(c(0,1)) +
  ylim(c(0,1)) + 
  labs(x = "Lineage score", y = "Ecological score") +
  theme_bw() +
  theme(aspect.ratio = 1)
```

![Example quadrant plot](imgs/example_quadrants.tif)

## Usage notes

CASCADE is specifically designed to analyze liver tumors. As such, it contains very specific requirements on how inputted data should look. For single-cell data, only the following 25 cell type labels will be considered in the calculation of the ecological score:

B cells  
Plasma cells  
CD4-CD69-memory T cells  
CD4-FOXP3-regulatory T cells  
CD4-IL7R-central memory T cells  
CD4-KLRB1-T cells  
CD8-CD69-memory T cells  
CD8-GZMH-effector T cells  
CD8-GZMK-effector memory T cells  
MAIT  
NK-CD160-tissue resident  
NK-GNLY-circulatory  
T cells-MKI67-proliferative  
c0-LUM-inflammatory CAF  
c1-MYH11-vascular CAF  
c2-APOA1-hepatocyte like CAF  
c0-S100A8-Monocyte  
c1-CXCL10-M1 Macrophage  
c2-CCL4L2-M2 Macrophage  
c3-TPSB2-Mast cells  
c0-VWF-endothelial  
c1-ANGPT2-endothelial  
c2-CRHBP-endothelial  
c3-CCL5-endothelial  
c4-RGS5-endothelial  

Your data does not need to contain all 25 of these cell types. You can label cells as types not included in this list, but those cells will be excluded from all score calculations. Make sure that labels are an exact match in spelling and case. 

Additionally, the lineage score bases its calculations on the expression of genes related to 4 specific tumor cell lineages: Cholangiocyte-like tumor cells, Hepatocyte-like tumor cells, Myeloid cell-like tumor cells, and Mesenchymal-like tumor cells. No action is required to accomodate this functionality, but do consider if these lineages are appropriate for your dataset.

