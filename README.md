# CASCADE

**CASCADE** , or **C**ancer **A**ggressiveness via **S**ingle-**C**ell **A**nalysis **D**uring **E**volution, is a bioinformatics tool that quantitatively analyzes liver tumors from single-cell or bulk transcriptomic data. 

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

CASCADE works with either single-cell or bulk transcriptomic data. For single-cell data, input either a single Seurat object containing malignancy information in the metadata or separate Seurat objects for malignant and nonmalignant cells. For bulk data, input a normalized matrix with genes on the rows and samples on the columns.

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
malignant_obj$sample_name <- c("s1", "s2", "s2")                        # your label vector would be much longer than this

nonmalignant_obj <- CreateSeuratObject(counts = nonmalignant_gene_counts)
nonmalignant_obj$sample_name <- c("s1", "s1", "s2")                     # your label vector would be much longer than this
nonmalignant_obj$cell_type <- c("MAIT", "B cells", "CD4-KLRB1-T cells") # your label vector would be much longer than this

cascade_df <- sc_CASCADE(tumor_obj = malignant_obj, nontumor_obj = nonmalignant_obj,
                         ident_col = "sample_name", subtype_col = "cell_type")
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
