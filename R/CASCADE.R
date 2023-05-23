### Core functions for CASCADE
# @author Mahler Revsine
# @date 12 January 2023 

#' Ecological (Eco) score calculation for single cell samples.
#' Yields a hazard score, resource score, and ecological score for each sample.
#' 
#' Hazard score describes the amount of hazards facing a tumor.
#' Resource score describes the amount of resources available to a tumor.
#' Eco score describes the overall state of the TME with respect to a tumor,
#' where a low score indicates hazardous conditions and a high score indicates
#' resource-heavy conditions. This metric is uniformly distributed from 0 to 1. 
#' 
#' Hazards and resources are described from the perspective of a tumor, where
#' hazards attack it while resources aid its survival.
#' 
#' @param tme_comp Data frame containing nontumor cell subtype composition for
#' each sample. Samples are on the rows, subtypes on the columns. Rows must be
#' scaled to sum to 1. 
#' @param subtype_scores Data frame containing nontumor cell subtype weights
#' and labels. Rownames are subtypes. Must contain a `category` column that 
#' labels each subtype as a "Hazard" or "Resource" and a `score` column that
#' lists a weight for each subtype. Defaults to a built-in weight table curated
#' from a single-cell study of liver cancer.
#' 
#' The following cell types are included in the default table. If using the 
#' default, include as many of these cell types as possible in your data. 
#' Not all cell types between this table and your data must overlap, but any
#' shared types must be an exact string match.
#' 
#' B cells
#' Plasma cells
#' CD4-CD69-memory T cells
#' CD4-FOXP3-regulatory T cells
#' CD4-IL7R-central memory T cells
#' CD4-KLRB1-T cells
#' CD8-CD69-memory T cells
#' CD8-GZMH-effector T cells
#' CD8-GZMK-effector memory T cells
#' MAIT
#' NK-CD160-tissue resident
#' NK-GNLY-circulatory
#' T cells-MKI67-proliferative
#' c0-LUM-inflammatory CAF
#' c1-MYH11-vascular CAF
#' c2-APOA1-hepatocyte like CAF
#' c0-S100A8-Monocyte
#' c1-CXCL10-M1 Macrophage
#' c2-CCL4L2-M2 Macrophage
#' c3-TPSB2-Mast cells
#' c0-VWF-endothelial
#' c1-ANGPT2-endothelial
#' c2-CRHBP-endothelial
#' c3-CCL5-endothelial
#' c4-RGS5-endothelial
#' 
#' @return Data frame with three columns: "hazards" - the hazard score, 
#' "resources" - the resource score, and "eco_score" - the ecological score
#' for each sample. Samples are row names.
#' 
#' @export
#' 
#' @example 
#' samples <- c("s1", "s2", "s3")
#' tme_comp <- data.frame(cell1 = c(0.2, 0.8, 0.1),
#'                        cell2 = c(0.5, 0.2, 0.4),
#'                        cell3 = c(0.3, 0.0, 0.5))
#' rownames(tme_comp) <- samples
#' 
#' subtype_scores <- data.frame(category = c("Hazard", "Resource", "Resource"),
#'                              score = c(0.9, 0.7, 0.4))
#' rownames(subtype_scores) <- c("cell1", "cell2", "cell3")
#' 
#' sc_ECO(tme_comp, subtype_scores)
#' 
sc_ECO <- function(tme_comp, subtype_scores = eco_subtype_scores) {
  
  # assigning important data to variables 
  samples <- rownames(tme_comp)
  subtypes <- rownames(subtype_scores)
  represented_subtypes <- colnames(tme_comp)
  
  # get weights and classifications of subtypes represented in the TME
  shared_subtypes <- represented_subtypes[represented_subtypes %in% subtypes]
  shared_st_weights <- subtype_scores[shared_subtypes, "score"]
  shared_st_labels <- subtype_scores[shared_subtypes, "category"]
  
  # calculate weighted TME composition of each sample
  weighted_tme_comp <- data.frame(t(t(as.matrix(tme_comp[shared_subtypes])) * 
                                      shared_st_weights))
  colnames(weighted_tme_comp) <- shared_subtypes
  
  # calculate hazard and resource scores
  hazard_scores <- rowSums(weighted_tme_comp[shared_st_labels == "Hazard"])
  resource_scores <- rowSums(weighted_tme_comp[shared_st_labels == "Resource"])
  
  # compile metrics and calculate eco score
  eco_df <- data.frame(hazards = hazard_scores,
                       resources = resource_scores,
                       eco_score = resource_scores / 
                         (resource_scores + hazard_scores))
  rownames(eco_df) <- samples
  
  eco_df
  
}

#' Lineage score calculation for single cell samples
#' 
#' Lineage score reflects characteristics of a tumor cell population. It is
#' uniformly distributed between 0 and 1, where a low score indicates a 
#' hepatocytic or myeloid cell-like lineage and a high score indicates a
#' cholangiocytic or mesenchymal lineage.
#' 
#' @param expr_data Matrix of normalized gene counts from a single-cell dataset.
#' Genes are on the rows and samples are on the columns.
#' @param idents Vector of sample labels for the columns of `expr_data`.
#' @param lineage_genes List of vectors containing genes for each lineage. 
#' Gene set vectors must be named after the lineages that they represent. 
#' Defaults to a built-in list of genes curated from a single-cell study of 
#' liver cancer. Contains gene sets for CLTCs (cholangiocyte-like tumor cells),
#' HLTCs (hepatocyte-like tumor cells), MCTCs (myeloid cell-like tumor cells),
#' and MLTCs (mesenchymal-like tumor cells). If an independent value is given
#' for this parameter, it must contain all 4 of these lineages in order to 
#' calculate a lineage score. 
#' @return Data frame containing average gene expression for each lineage and 
#' the final lineage score in the "lineage_score" column. Samples are row names.
#' 
#' @export
#' 
#' @example 
#' expr_data: cell1 cell2 cell3 ...
#'      gene1     0   0.1   1.2
#'      gene2   0.9     0     0
#'      gene3   2.2   2.6   1.5  
#'      .
#'      .
#'      .
#' idents <- c("Sample1", "Sample1", "Sample2")
#' lineage_genes <- list(CLTC = c("gene12", "gene23", "gene45"),
#'                       HLTC = c("gene1", "gene91", "gene92", "gene93"),
#'                       MCTC = c("gene2", "gene96", "gene97", "gene98"),
#'                       MLTC = c("gene5", "gene6", "gene14"))
#'                       
#' sc_LINEAGE(expr_data, idents, lineage_genes)
#' 
sc_LINEAGE <- function(expr_data, idents, lineage_genes = top_lineage_genes) {
  
  # assigning important data to variables 
  lineages <- names(lineage_genes)
  samples <- unique(idents)
  
  # expression of each lineage in each cell
  sample_lineage_expr <- data.frame(do.call(cbind, 
    lapply(lineages, function(lin) {
      lineage_gene_set <- lineage_genes[[lin]]
      incl_lineage_gene_set <- lineage_gene_set[
        lineage_gene_set %in% rownames(expr_data)
      ]
      missing_genes <- setdiff(lineage_gene_set, incl_lineage_gene_set)
      if (length(missing_genes) == length(lineage_gene_set)) {
        stop(paste0("None of the genes needed for determining ", lin, 
                    " lineage score are included in the expression data. ",
                    "CASCADE scores cannot be calculated for this dataset."))
      }
      if (length(missing_genes) != 0) {
        warning(paste0("Not all genes needed for determining ", lin, 
                       " lineage score are included in the expression data. ",
                       length(missing_genes), "/", length(lineage_gene_set), 
                       " genes are missing: ", 
                       paste(missing_genes, collapse = ","), "."))
      }
      if (length(incl_lineage_gene_set) == 1) {
        expr_data[incl_lineage_gene_set, ]
      } else {
        colMeans(expr_data[incl_lineage_gene_set, ])
      }
    })))
  colnames(sample_lineage_expr) <- lineages
  
  # expression of each lineage in each sample, averaged across its cells
  lineage_df <- data.frame(do.call(rbind, lapply(samples, function(sample) {
    sample_idxs <- idents == sample
    if (sum(sample_idxs) == 1) { # this is possible if `min_n_cells` is 1
      sample_lineage_expr[idents == sample, ]
    } else {
      colMeans(sample_lineage_expr[idents == sample, ])
    }
  })))
  colnames(lineage_df) <- lineages
  rownames(lineage_df) <- samples
  
  # calculate lineage score
  lineage_df$lineage_score <- (lineage_df$CLTC + lineage_df$MLTC) / 
    (lineage_df$CLTC + lineage_df$HLTC + lineage_df$MCTC + lineage_df$MLTC)
  rownames(lineage_df) <- samples
  
  lineage_df
  
}

#' CASCADE calculation for single cell samples.
#'
#' @param total_obj Seurat object containing all cells. Either this parameter
#' or both `tumor_obj` and `nontumor_obj` must be supplied. If `total_obj` is
#' given, then this object must contain a metadata column with malignancy info
#' specified by the `malignancy_col` parameter.
#' @param malignancy_col Column name for metadata column containing malignancy
#' information for `total_obj`. This column must be a string vector with 
#' "malignant" and "non-malignant" labels. Not needed if `total_obj` is null.
#' @param tumor_obj Seurat object containing all malignant cells. Ignored if
#' `total_obj` is supplied. If `total_obj` is not supplied, then both 
#' `tumor_obj` and `nontumor_obj` must be supplied. 
#' @param nontumor_obj Seurat object containing all non-malignant cells. 
#' Ignored if `total_obj` is supplied. If `total_obj` is not supplied, then 
#' both `tumor_obj` and `nontumor_obj` must be supplied. 
#' @param ident_col Name of the metadata column containing a sample label for
#' each cell. Must be included in the metadata of each supplied Seurat object.
#' @param subtype_col Name of the metadata column containing a subtype label for
#' each cell. Must be included in the metadata of `total_obj` or `nontumor_obj`.
#' @param min_n_cells Integer number of the minimum required cell count for 
#' each sample. Samples with fewer than this number of cells will be excluded.
#' @param subtype_scores Data frame containing nontumor cell subtype weights
#' and labels. Rownames are subtypes. Must contain a `category` column that 
#' labels each subtype as a "Hazard" or "Resource" and a `score` column that
#' lists a weight for each subtype. Defaults to a built-in weight table curated
#' from a single-cell study of liver cancer.
#' @param lineage_genes List of vectors containing genes for each lineage. 
#' Gene set vectors must be named after the lineages that they represent. 
#' Defaults to a built-in list of genes curated from a single-cell study of 
#' liver cancer. Contains gene sets for CLTCs (cholangiocyte-like tumor cells),
#' HLTCs (hepatocyte-like tumor cells), MCTCs (myeloid cell-like tumor cells),
#' and MLTCs (mesenchymal-like tumor cells). If an independent value is given
#' for this parameter, it must contain all 4 of these lineages in order to 
#' calculate a lineage score. 
#' @return Data frame with the following columns:
#'   CLTC: expression of cholangiocyte-like tumor cell markers in tumor cells
#'   HLTC: expression of hepatocyte-like tumor cell markers in tumor cells
#'   MCTC: expression of myeloid cell-like tumor cell markers in tumor cells
#'   MLTC: expression of mesenchymal-like tumor cell markers in tumor cells
#'   additional columns for any other lineages supplied in `lineage_genes`
#'   lineage_score: lineage score in each sample
#'   hazards: amount of hazards facing the non-malignant cells
#'   resources: amount of resources available to the non-malignant cells
#'   eco_score: eco score in each sample
#'   tumor_score: tumor score in each sample
#'   quadrant: lineage-eco quadrant for each sample
#'
#' @export
#' 
#' @examples 
#' sc_CASCADE(total_obj = x, malignancy_col = "mal_anno", ident_col = "sample")
#' sc_CASCADE(tumor_obj = t, nontumor_obj = n, subtype_col = "cell_type")
#' 
sc_CASCADE <- function(total_obj = NULL, 
                       tumor_obj = NULL,
                       nontumor_obj = NULL,
                       malignancy_col = "malignancy",
                       ident_col = "ident.name", 
                       subtype_col = "subtype", 
                       min_n_cells = 15,
                       subtype_scores = eco_subtype_scores,
                       lineage_genes = top_lineage_genes) {
  
  # Make sure total_obj or both tumor_obj and nontumor_obj are supplied
  if (is.null(total_obj) & (is.null(tumor_obj) | is.null(nontumor_obj))) {
    stop(paste0("Must supply an object containing all cells to total_obj or ",
                "two objects to tumor_obj and nontumor_obj containing ",
                "all malignant and nonmalignant cells respectively."))
  }
  
  # If total_obj is supplied, we ignore tumor_obj and nontumor_obj args
  # We will divide it into tumor_obj and nontumor_obj manually
  if (!is.null(total_obj)) {
    
    # Make sure total_obj contains a proper malignancy metadata column
    if (!(malignancy_col %in% colnames(total_obj@meta.data))) {
      stop(paste0("No column named ", malignancy_col," in total_obj metadata."))
    } else if (!("malignant" %in% total_obj@meta.data[[malignancy_col]]) |
               !("non-malignant" %in% total_obj@meta.data[[malignancy_col]])) {
      stop(paste0("Malignancy column in total_obj is not formatted properly ",
                  "or does not contain both malignant and non-malignant ",
                  "labels. Make sure this column is a string vector with ",
                  "malignant cells labelled as 'malignant' and non-malignant ",
                  "cells labelled as 'non-malignant' for proper use."))
    }
    
    # Make sure total_obj contains an identity column
    if (!(ident_col %in% colnames(total_obj@meta.data))) {
      stop(paste0("No column named ", ident_col, " in total_obj metadata."))
    } 
    
    # Make sure total_obj contains a subtype column
    # Throw an error if it contains none of our subtypes
    # Give a warning if it is missing at least 1 subtype 
    if (!(subtype_col %in% colnames(total_obj@meta.data))) {
      stop(paste0("No column named ", subtype_col, " in total_obj metadata."))
    } else {
      represented_subtypes <- unique(total_obj@meta.data[[subtype_col]])
      missing_subtypes <- subtype_scores$subtype[!(subtype_scores$subtype %in% 
                                                     represented_subtypes)]
      if (length(missing_subtypes) == nrow(subtype_scores)) {
        stop(paste0("total_obj does not contain any of the needed subtypes. ",
                    "Please make sure the ", subtype_col, " metadata column ",
                    "follows the naming convention outlined in the ",
                    "documentation."))
      } else if (length(missing_subtypes) > 0) {
        warning(paste0("There are ", length(missing_subtypes), " missing ",
                       "subtypes in the ", subtype_col, " column of ",
                       "total_obj's metadata. Missing subtypes: ",
                       paste(missing_subtypes, collapse = ", ")))
      }
    }
    
    tumor_idxs <- total_obj@meta.data[[malignancy_col]] == "malignant"
    tumor_cells <- colnames(total_obj)[tumor_idxs]
    tumor_obj <- subset(total_obj, cells = tumor_cells)
    
    nontumor_idxs <- total_obj@meta.data[[malignancy_col]] == "non-malignant"
    nontumor_cells <- colnames(total_obj)[nontumor_idxs]
    nontumor_obj <- subset(total_obj, cells = nontumor_cells)
    
    # To reach this block, both tumor and nontumor obj will have been supplied
  } else {
    
    # Make sure tumor_obj contains an identity column
    if (!(ident_col %in% colnames(tumor_obj@meta.data))) {
      stop(paste0("No column named ", ident_col, " in tumor_obj metadata."))
    } 
    
    # Make sure nontumor_obj contains an identity column
    if (!(ident_col %in% colnames(nontumor_obj@meta.data))) {
      stop(paste0("No column named ", ident_col, " in nontumor_obj metadata."))
    } 
    
    # Make sure nontumor_obj contains a subtype column
    # Throw an error if it contains none of our subtypes
    # Give a warning if it is missing at least 1 subtype 
    if (!(subtype_col %in% colnames(nontumor_obj@meta.data))) {
      stop(paste0("No column named ", subtype_col," in nontumor_obj metadata."))
    } else {
      represented_subtypes <- unique(nontumor_obj@meta.data[[subtype_col]])
      missing_subtypes <- subtype_scores$subtype[!(subtype_scores$subtype %in% 
                                                     represented_subtypes)]
      if (length(missing_subtypes) == nrow(subtype_scores)) {
        stop(paste0("nontumor_obj does not contain any of the needed ", 
                    "subtypes. Please make sure the ", subtype_col, " ",
                    "metadata column follows the naming convention outlined ",
                    "in the documentation."))
      } else if (length(missing_subtypes) > 0) {
        warning(paste0("There are ", length(missing_subtypes), " missing ",
                       "subtypes in the ", subtype_col, " column of ",
                       "nontumor_obj's metadata. Missing subtypes: ",
                       paste(missing_subtypes, collapse = ", "), "."))
      }
    }
    
  }
  
  # filter tumor_obj by min_n_cells in each sample
  tumor_sample_cts <- table(tumor_obj@meta.data[[ident_col]])
  filt_tumor_samples <- names(tumor_sample_cts)[tumor_sample_cts >= min_n_cells]
  if (length(filt_tumor_samples) == 0) {
    stop(paste0("No samples contained at least ", min_n_cells, " malignant ",
                "cells. Lower the value of min_n_cells or make sure your ",
                "input data contains an adequate number of malignant cells ",
                "for each sample."))
  } else if (length(filt_tumor_samples) < length(tumor_sample_cts)) {
    message(paste0(length(tumor_sample_cts) - length(filt_tumor_samples), 
                   " sample(s) did not have at least ", min_n_cells, 
                   " malignant cells and were therefore discarded. ",
                   "Discarded sample(s): ", 
                   paste(setdiff(names(tumor_sample_cts),filt_tumor_samples),
                         collapse = ", ")), ".")
  }
  filt_tumor_idxs <- tumor_obj@meta.data[[ident_col]] %in% filt_tumor_samples
  filt_tumor_cells <- colnames(tumor_obj)[filt_tumor_idxs]
  filt_tumor_obj <- subset(tumor_obj, cells = filt_tumor_cells)
  
  # filter nontumor_obj by min_n_cells in each sample
  nontumor_sample_cts <- table(nontumor_obj@meta.data[[ident_col]])
  filt_nontumor_samples <- names(nontumor_sample_cts)[nontumor_sample_cts >= 
                                                        min_n_cells]
  if (length(filt_nontumor_samples) == 0) {
    stop(paste0("No samples contained at least ", min_n_cells, " nonmalignant ",
                "cells. Lower the value of min_n_cells or make sure your ",
                "input data contains an adequate number of nonmalignant cells ",
                "for each sample."))
  } else if (length(filt_nontumor_samples) < length(nontumor_sample_cts)) {
    message(paste0(length(nontumor_sample_cts) - length(filt_nontumor_samples), 
                   " sample(s) did not have at least ", min_n_cells, 
                   " nonmalignant cells and were therefore discarded. ",
                   "Discarded sample(s): ", 
                   paste(setdiff(names(nontumor_sample_cts),
                                 filt_nontumor_samples),
                         collapse = ", ")), ".")
  }
  filt_nontumor_idxs <- nontumor_obj@meta.data[[ident_col]] %in% 
    filt_nontumor_samples
  filt_nontumor_cells <- colnames(nontumor_obj)[filt_nontumor_idxs]
  filt_nontumor_obj <- subset(nontumor_obj, cells = filt_nontumor_cells)
  
  # Calculate lineage score
  lineage_df <- sc_LINEAGE(expr_data =as.matrix(filt_tumor_obj@assays$RNA@data), 
                           idents = filt_tumor_obj@meta.data[[ident_col]], 
                           lineage_genes = top_lineage_genes)
  
  # Calculate eco score using TME composition profiles
  nontumor_tme_comp <- table(filt_nontumor_obj@meta.data[c(ident_col,
                                                           subtype_col)])
  nontumor_tme_comp <- as.data.frame(t(apply(nontumor_tme_comp,1,function(row) {
    row / sum(row)
  })))
  eco_df <- sc_ECO(tme_comp = nontumor_tme_comp, 
                   subtype_scores = eco_subtype_scores)
  
  # Calculate tumor score
  cascade_df <- merge(lineage_df, eco_df, by = "row.names", all = TRUE)
  rownames(cascade_df) <- cascade_df$Row.names
  cascade_df <- cascade_df[2:ncol(cascade_df)]
  cascade_df <- cascade_df[gtools::mixedsort(rownames(cascade_df)), ]
  cascade_df$tumor_score <- cascade_df$lineage * cascade_df$eco
  
  # Assign quadrants
  cascade_df$quadrant <- apply(cascade_df, 1, function(row) {
    lineage <- row[["lineage_score"]]
    eco <- row[["eco_score"]]
    if (is.na(lineage) | is.na(eco)) {
      NA
    } else {
      lineage_char <- if (lineage < 0.5) "A" else "B"
      eco_char <- if (eco < 0.5) "1" else "2"
      paste0(lineage_char, eco_char)
    }
  })
  
  cascade_df
  
}

#' CASCADE calculation for bulk transcriptomic samples.
#' 
#' @param bulk_expr Matrix of normalized gene counts from bulk dataset.
#' @param degs Built-in gene set to calculate CASCADE values in bulk datasets.
#' @return Data frame with the following columns:
#'   L1: expression of low-lineage score markers
#'   L2: expression of high-lineage score markers
#'   E1: expression of low-eco score markers
#'   E2: expression of high-eco score markers
#'   lineage_score: lineage score in each sample
#'   eco_score: eco score in each sample
#'   tumor_score: tumor score in each sample
#'   quadrant: lineage-eco quadrant for each sample
#'
#' @export
#' 
#' @example 
#' bulk_expr: cell1 cell2 cell3 ...
#'      gene1   5.2   4.9   6.1
#'      gene2   6.3   5.8   7.0
#'      gene3   6.6   7.4   6.9  
#'      .
#'      .
#'      .
#' bulk_CASCADE(bulk_expr)
#' 
bulk_CASCADE <- function(bulk_expr, degs = cascade_degs) {
  
  # function to scale input vector from 0 to 1
  scale_vector <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }
  
  # filter down to genes included in the dataset
  incl_degs <- lapply(names(degs), function(geneset) {
    genes <- degs[[geneset]]
    incl_genes <- genes[genes %in% rownames(bulk_expr)]
    missing_genes <- setdiff(genes, incl_genes) 
    if (length(missing_genes) == length(genes)) {
      stop(paste0("None of the ", geneset, " genes needed for determining ",
                  "CASCADE score from bulk transcriptomics data are included ",
                  "in the supplied dataset. CASCADE scores cannot be ",
                  "calculated for this dataset."))
    }
    if (length(missing_genes) != 0) {
      warning(paste0("Not all ", geneset, " genes needed for determining ",
                     "CASCADE scores are included in the expression data. ",
                     length(missing_genes), "/", length(genes), 
                     " genes are missing: ", 
                     paste(missing_genes, collapse = ","), "."))
    }
    incl_genes
  })
  names(incl_degs) <- names(degs)
  
  # calculate expression of each gene set in each bulk sample
  cascade_df <- data.frame(do.call(cbind, lapply(incl_degs, function(genes) {
    colMeans(bulk_expr[genes, ])
  })))
  colnames(cascade_df) <- names(incl_degs)
  rownames(cascade_df) <- colnames(bulk_expr)
  
  # calculate scaled lineage, eco, and tumor scores for each sample
  cascade_df$lineage_score <- scale_vector(cascade_df$L2 - cascade_df$L1)
  cascade_df$eco_score <- scale_vector(cascade_df$E2 - cascade_df$E1)
  cascade_df$tumor_score <- cascade_df$lineage_score * cascade_df$eco_score
  
  # calculate quadrant for each sample
  cascade_df$quadrant <- apply(cascade_df, 1, function(row) {
    lineage <- row[["lineage_score"]]
    eco <- row[["eco_score"]]
    if (is.na(lineage) | is.na(eco)) {
      NA
    } else {
      lin_char <- if (lineage < median(cascade_df$lineage_score)) "A" else "B"
      eco_char <- if (eco < median(cascade_df$eco_score)) "1" else "2"
      paste0(lin_char, eco_char)
    }
  })
  
  cascade_df
  
}
