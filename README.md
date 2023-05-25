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

Single Seurat object:
```r
obj <- CreateSeuratObject(counts = gene_counts)
obj@meta.data$malignancy_label <- c("malignant", "non-malignant", malignant") # your label vector would be much longer than this
cascade_df <- sc_CASCADE(total_obj = obj, malignancy_col = "malignancy_label")
```

Separate malignant and non-malignant objects:
```r
malignant_obj <- CreateSeuratObject(counts = malignant_gene_counts)
nonmalignant_obj <- CreateSeuratObject(counts = nonmalignant_gene_counts)
cascade_df <- sc_CASCADE(tumor_obj = malignant_obj, nontumor_obj = nonmalignant_obj)
```

Bulk data:
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
