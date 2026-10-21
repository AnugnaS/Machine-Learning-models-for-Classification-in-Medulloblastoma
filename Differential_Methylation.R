#=========================================================================================================#
# ----------- DMP analysis on M-values (limma) -> significant CpGs and genes -> Random Forest
# trained on the M-values of the selected CpGs -> tested on Test NMB and Cavalli-----------------------
#=========================================================================================================#

# ------ Load Packages and Functions ----

source("~/Custom_functions.R")

cran_pgs <- c("pheatmap", "RColorBrewer", "caret", "randomForest", "MLmetrics", "pROC", "ggplot2")
bioc_pgs <- c("limma", "minfi", "missMethyl", "IlluminaHumanMethylation450kanno.ilmn12.hg19") 

sapply(cran_pgs, requireNamespace, quietly = TRUE)

sapply(bioc_pgs, requireNamespace, quietly = TRUE)

missing_cran <- cran_pgs[!sapply(cran_pgs, requireNamespace, quietly = TRUE)]

if (length(missing_cran) > 0) install.packages(missing_cran, dependencies = TRUE)

missing_bioc <- bioc_pgs[!sapply(bioc_pgs, requireNamespace, quietly = TRUE)]

if (length(missing_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(limma)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(pheatmap)
  library(caret)
  library(ggplot2)
})

# ---- Settings ----

adj_p_cut      <- 0.05    # FDR Threshold
min_abs_logFC  <- 0       # difference in DMP coefficients
n_features     <- 10000   # top 10k to be selected
n_heat         <- 100     # top number of probes to be applied for heatmaps

HISTO_LEVELS  <- c("CLA", "DN", "LCA")
HISTO_COLOURS <- c("CLA" = "purple", "LCA" = "pink", "DN" = "darkgreen")

out_dir <- "~/Thesis_models/DMP_Mval/"

dir.create(path.expand(out_dir), recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(path.expand(out_dir), "evaluation"), recursive = TRUE, showWarnings = FALSE)


# ---- Training Set ----

NMB_train_pheno <- readRDS("~/NMB/NMB_train_pheno.rds")

# DN and MBEN due to their similarities and frequent occurence in SHH-pathway activated tumours, are combined into one class as "DN"

NMB_train_pheno$CPR_Histology<-gsub("MBEN","DN", NMB_train_pheno$CPR_Histology)

NMB_train_pheno <- NMB_train_pheno[!is.na(NMB_train_pheno$Sample_ID), ]

NMB_train_pheno$CPR_Histology <- factor(NMB_train_pheno$CPR_Histology, levels = HISTO_LEVELS)

NMB_train_pheno$Mol_group     <- factor(NMB_train_pheno$Mol_group)

message("Training samples: histology (rows) by molecular group (columns)")

print(table(NMB_train_pheno$CPR_Histology, NMB_train_pheno$Mol_group))

# ---- Test Sets ----

NMB_test_pheno <- readRDS("~/NMB/NMB_test_pheno.rds")

# DN and MBEN due to their similarities and frequent occurence in SHH-pathway activated tumours, are combined into one class as "DN"

NMB_test_pheno$CPR_Histology<-gsub("MBEN","DN", NMB_test_pheno$CPR_Histology)

NMB_test_pheno <- NMB_test_pheno[!is.na(NMB_test_pheno$Sample_ID), ]

NMB_test_pheno$CPR_Histology<- factor(NMB_test_pheno$CPR_Histology, levels = HISTO_LEVELS)

Cavalli_pheno <- readRDS("~/Cavalli/Final_Cavalli_pheno.rds")

# DN and MBEN due to their similarities and frequent occurence in SHH-pathway activated tumours, are combined into one class as "DN"

Cavalli_pheno$histology<-gsub("MBEN","DN",Cavalli_pheno$histology)

Cavalli_pheno<- Cavalli_pheno[!is.na(Cavalli_pheno$Sample_ID), ]

Cavalli_pheno$histology <- factor(Cavalli_pheno$histology, levels = HISTO_LEVELS)


# ---- M-values: NMB (450k and EPIC samples combined) ----

NMB_450k_Mval <- minfi::getM(readRDS("~/Final450k/MPP_450k_mSet.rds"))
NMB_EPIC_Mval <- minfi::getM(readRDS("~/EPIC/NMBEPIC/NMB_EPIC_mSet.rds"))

shared_probes <- intersect(rownames(NMB_450k_Mval), rownames(NMB_EPIC_Mval))
NMB_Mval <- cbind(NMB_450k_Mval[shared_probes, , drop = FALSE],   
                  NMB_EPIC_Mval[shared_probes, , drop = FALSE])
rm(NMB_450k_Mval, NMB_EPIC_Mval)
stopifnot("duplicated sample names across the 450k and EPIC sets" = !anyDuplicated(colnames(NMB_Mval)))
saveRDS(NMB_Mval, "~/NMB/NMB_mval.rds")

# ---- M-values: Cavalli ----

Cavalli_Mval <- minfi::getM(readRDS("~/Cavalli/Cavalli-test/Cavalli_test_mSet.rds"))

#----- M value matrices ----

train_NMB_Mval <- NMB_Mval[, NMB_train_pheno$Sample_ID, drop = FALSE]
test_NMB_Mval  <- NMB_Mval[, NMB_test_pheno$Sample_ID,  drop = FALSE]

Cavalli_Mval   <- Cavalli_Mval[, match(as.character(Cavalli_pheno$Sample_ID),colnames(Cavalli_Mval))]

NMB_train_pheno<-NMB_train_pheno[match(colnames(train_NMB_Mval),NMB_train_pheno$Sample_ID),]

train_NMB_Mval<-train_NMB_Mval[,match(NMB_train_pheno$Sample_ID,colnames(train_NMB_Mval))]

NMB_test_pheno<-NMB_test_pheno[match(colnames(test_NMB_Mval),NMB_test_pheno$Sample_ID),]

test_NMB_Mval<-test_NMB_Mval[,match(NMB_test_pheno$Sample_ID,colnames(test_NMB_Mval))]

Cavalli_pheno<-Cavalli_pheno[match(colnames(Cavalli_Mval),Cavalli_pheno$Sample_ID),] 

stopifnot("columns of train_NMB_Mval do not follow NMB_train_pheno" =
            identical(colnames(train_NMB_Mval), as.character(NMB_train_pheno$Sample_ID)))


#-------- limma DMP analysis (train set samples only) -----------
# Model: M-value ~ histology + molecular group (+ Platform). Histology and molecular group are tangled
# so each histology difference is estimated after adjusting for
# molecular group


design_terms   <- c("CPR_Histology", "Mol_group", "Platform")                

design_formula <- as.formula(paste("~ 0 +", paste(design_terms, collapse = " + ")))

design         <- model.matrix(design_formula, data = NMB_train_pheno)

colnames(design) <- make.names(colnames(design))                           

stopifnot("samples were dropped while building the design (missing values?)" = nrow(design) == nrow(NMB_train_pheno))

not_estimable <- limma::nonEstimable(design)

if (!is.null(not_estimable)) stop("the design is not full rank; non-estimable terms: ", paste(not_estimable, collapse = ", "))

fit <- limma::lmFit(train_NMB_Mval, design)

# ---- Pairwise Histology Comparisons ----

hist_levels <- levels(NMB_train_pheno$CPR_Histology)

hist_pairs  <- t(combn(hist_levels, 2))

contrast_exprs <- paste0("CPR_Histology", hist_pairs[, 1], " - CPR_Histology", hist_pairs[, 2])

contrast.matrix <- limma::makeContrasts(contrasts = contrast_exprs, levels = design)

colnames(contrast.matrix) <- paste0(hist_pairs[, 1], "_vs_", hist_pairs[, 2])

fit2 <- limma::contrasts.fit(fit, contrast.matrix)

fit2 <- limma::eBayes(fit2)

omni_exprs <- paste0("CPR_Histology", hist_levels[-1], " - CPR_Histology", hist_levels[1])

omni_matrix <- limma::makeContrasts(contrasts = omni_exprs, levels = design)

colnames(omni_matrix) <- paste0(hist_levels[-1], "_vs_", hist_levels[1])

fit_omni <- limma::eBayes(limma::contrasts.fit(fit, omni_matrix))

omnibus  <- limma::topTable(fit_omni, coef = seq_len(ncol(omni_matrix)), number = Inf,
                            adjust.method = "BH", sort.by = "F")


#------ Significant CpGs, their genes and regions-----


# ---- annotation of the 450k probes ----

ann <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

ann_df <- as.data.frame(ann[, c("chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island")])

ann_df <- ann_df[match(rownames(train_NMB_Mval), rownames(ann_df)), ]

rownames(ann_df) <- rownames(train_NMB_Mval) 

gene_text <- ann_df$UCSC_RefGene_Name;  gene_text[is.na(gene_text)] <- ""

ann_df$gene <- vapply(strsplit(gene_text, ";", fixed = TRUE),
                      function(g) paste(unique(g[nzchar(g)]), collapse = ";"), character(1))

region_text <- ann_df$UCSC_RefGene_Group;  region_text[is.na(region_text)] <- ""

ann_df$promoter <- grepl("TSS1500|TSS200|5'UTR|1stExon", region_text)

# ---- one table per contrast ----

dmp_sig     <- list()
dmp_summary <- data.frame()

for (k in seq_len(nrow(hist_pairs))) {
  
  cn     <- colnames(contrast.matrix)[k]
  first  <- hist_pairs[k, 1]
  second <- hist_pairs[k, 2]
  
  tt  <- limma::topTable(fit2, coef = cn, number = Inf, adjust.method = "BH", sort.by = "P")
  idx <- match(rownames(tt), rownames(ann_df))
  
  tt$probe           <- rownames(tt)
  tt$higher_in       <- ifelse(tt$logFC > 0, first, second)           
  tt$chr             <- ann_df$chr[idx] 
  tt$pos             <- ann_df$pos[idx]
  tt$gene            <- ann_df$gene[idx]
  tt$region          <- ann_df$UCSC_RefGene_Group[idx]
  tt$island_relation <- ann_df$Relation_to_Island[idx]
  tt$promoter        <- ann_df$promoter[idx]
  
  sig <- tt[!is.na(tt$adj.P.Val) & tt$adj.P.Val < adj_p_cut & abs(tt$logFC) >= min_abs_logFC, ]
  
  dmp_sig[[cn]] <- sig
  
  write.csv(sig, file.path(path.expand(out_dir), paste0("DMP_", cn, "_significant.csv")), row.names = FALSE)
  
  genes_here <- unique(unlist(strsplit(sig$gene, ";", fixed = TRUE)))
  
  if (length(genes_here) > 0) {
    
    gene_counts <- sort(table(unlist(strsplit(sig$gene, ";", fixed = TRUE))), decreasing = TRUE)
    
    write.csv(data.frame(gene = names(gene_counts), n_significant_cpgs = as.integer(gene_counts)),
              file.path(path.expand(out_dir), paste0("DMP_", cn, "_genes.csv")), row.names = FALSE)
  }
  
  dmp_summary <- rbind(dmp_summary, data.frame(contrast = cn, n_tested = nrow(tt), n_significant = nrow(sig),
                                               higher_in_first_group = sum(sig$logFC > 0),
                                               higher_in_second_group = sum(sig$logFC < 0),
                                               n_in_promoters = sum(sig$promoter %in% TRUE),
                                               n_genes = length(genes_here)))
}

message("Significant DMPs per comparison (FDR < ", adj_p_cut, ", |M difference| >= ", min_abs_logFC, "):")

print(dmp_summary, row.names = FALSE)

write.csv(dmp_summary, file.path(path.expand(out_dir), "DMP_summary.csv"), row.names = FALSE)

# ---- across comparisons ----

all_sig_cpgs    <- unlist(lapply(dmp_sig, rownames))

sig_cpgs_all  <- unique(all_sig_cpgs)

sig_genes_all <- unique(unlist(lapply(dmp_sig, function(s) unlist(strsplit(s$gene, ";", fixed = TRUE)))))

message("CpGs significant in at least one comparison: ", length(sig_cpgs_all),
        " | genes: ", length(sig_genes_all))

message("CpGs significant in (number of comparisons = number of CpGs):")
print(table(table(all_sig_cpgs)))

saveRDS(list(sig_cpgs = sig_cpgs_all, sig_genes = sig_genes_all),
        file.path(path.expand(out_dir), "DMP_sig_cpgs_and_genes.rds"))

# ---- the overall test ----

omni_idx <- match(rownames(omnibus), rownames(ann_df))

omnibus$probe    <- rownames(omnibus)

omnibus$gene     <- ann_df$gene[omni_idx]

omnibus$region   <- ann_df$UCSC_RefGene_Group[omni_idx]

omnibus$promoter <- ann_df$promoter[omni_idx]

omni_effect      <- apply(abs(as.matrix(omnibus[, colnames(omni_matrix), drop = FALSE])), 1, max) 

omni_sig         <- omnibus[!is.na(omnibus$adj.P.Val) & omnibus$adj.P.Val < adj_p_cut & omni_effect >= min_abs_logFC, ]

message("CpGs differing between any histologies (overall F test): ", nrow(omni_sig))

write.csv(head(omni_sig, 50000), file.path(path.expand(out_dir), "DMP_overall_F_test_top50000.csv"), row.names = FALSE)


#-------GO enrichment of the significant CpGs--------


for (cn in names(dmp_sig)) {
   
  library(missMethyl)
    
    sig_cpg <- head(rownames(dmp_sig[[cn]]), 50000)                      
    if (length(sig_cpg) < 50) { message("  ", cn, ": fewer than 50 significant CpGs, enrichment skipped"); next }
    
    go <- tryCatch(
      missMethyl::gometh(sig.cpg = sig_cpg, all.cpg = rownames(train_NMB_Mval), collection = "GO",
                         array.type = "450K", plot.bias = FALSE),
      error = function(e) { message("  enrichment failed for ", cn, ": ", conditionMessage(e)); NULL })
    if (is.null(go)) next
    
    write.csv(missMethyl::topGSA(go, number = 50), file.path(path.expand(out_dir), paste0("GO_", cn, ".csv")))
  }


#------- Heatmaps of the top DMPs per comparison---------
  
annotation_col <- data.frame(Mol_group = NMB_train_pheno$Mol_group,
                               Histology = NMB_train_pheno$CPR_Histology,
                               MYC_status = NMB_train_pheno$`MYC Status`,
                               MYCN_status = NMB_train_pheno$`MYCN Status`,
                               row.names = as.character(NMB_train_pheno$Sample_ID))
  
mol_colors <- c("WNT-MB" = "blue", "SHH-MB" = "red", "Group3-MB" = "yellow3", "Group4-MB" = "darkgreen")
  
MYC_colours<-c("0"= "grey", "1" = "red4")
  
MYCN_colours<-c("0" = "grey", "1" = "orange")
  
annotation_colors <- list(Mol_group = mol_colors, Histology = HISTO_COLOURS, MYC_status = MYC_colours, MYCN_status = MYCN_colours)
  
heat_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdBu")))(100)
  
for (cn in names(dmp_sig)) {
  sig <- dmp_sig[[cn]] 
  for (region_name in c("promoter", "non_promoter")) {
    s <- if (region_name == "promoter") sig[sig$promoter %in% TRUE, ] else sig[!(sig$promoter %in% TRUE), ]
    if (nrow(s) < 2) {
      message(cn, " | ", region_name, ": skipped, fewer than 2 CpGs")
      next
    } 
    
    probes <- head(s$probe, n_heat)
    
    pdf(file.path(path.expand(out_dir), paste0("Heatmap_", region_name, "_DMPs_", cn, ".pdf")),
        width = 8, height = 10)
    
    tryCatch(
      pheatmap::pheatmap(train_NMB_Mval[probes, , drop = FALSE], scale = "row",
                         annotation_col = annotation_col, annotation_colors = annotation_colors,
                         show_rownames = FALSE, show_colnames = FALSE, color = heat_colors,
                         main = paste0(gsub("_", " ", region_name), " DMPs: ", gsub("_", " ", cn),
                                       " (top ", length(probes), ")")),
      error = function(e) message("   heatmap failed: ", conditionMessage(e)),
      finally = dev.off())
  }
}

# With M-values adjusted for Molecular-group effects

for (cn in names(dmp_sig)) {
    sig <- dmp_sig[[cn]]
      for (region_name in c("promoter", "non_promoter")) {
      s <- if (region_name == "promoter") sig[sig$promoter %in% TRUE, ] else sig[!(sig$promoter %in% TRUE), ]
      if (nrow(s) < 2) next
      
      probes <- head(s$probe, n_heat) 
      
      mat_adj <- limma::removeBatchEffect(
        train_NMB_Mval[probes, , drop = FALSE],
        batch  = NMB_train_pheno$Mol_group,
        design = model.matrix(~ CPR_Histology, data = NMB_train_pheno))
      
      
      pdf(file.path(path.expand(out_dir), paste0("Heatmap_", region_name, "_DMPs_", cn, "adj.pdf")), width = 8, height = 10)
      
      pheatmap::pheatmap(mat_adj, scale = "row",
                         annotation_col = annotation_col, annotation_colors = annotation_colors,
                         show_rownames = FALSE, show_colnames = FALSE, color = heat_colors,
                         main = paste0(gsub("_", " ", region_name), " DMPs: ", gsub("_", " ", cn), " (top ", length(probes) , ")"))
      dev.off()
    }
  } 


# ------- The CpG feature set for the ML model-------
# The top n_features CpGs of the overall test (ranked by the F statistic) among those that pass the thresholds.

feature_cpgs <- head(omni_sig$probe, n_features)

if (length(feature_cpgs) < n_features) message("only ", length(feature_cpgs), " CpGs pass the thresholds (asked for ", n_features, ")")
stopifnot("too few CpGs for a model" = length(feature_cpgs) >= 50)

saveRDS(feature_cpgs, file.path(path.expand(out_dir), "DMP_feature_cpgs.rds"))
message("Feature set: ", length(feature_cpgs), " CpGs, ", sum(ann_df[feature_cpgs, "promoter"]), " in promoters")



#------- Train a Random Forest on the M-values of the feature set --------------------

histo   <- NMB_train_pheno$CPR_Histology                                           
train_x <- as.data.frame(t(train_NMB_Mval[feature_cpgs, , drop = FALSE]))          
stopifnot(identical(rownames(train_x), as.character(NMB_train_pheno$Sample_ID)))

ctrl <- caret::trainControl(method = "cv", 
                            number = 5, 
                            classProbs = TRUE,
                            summaryFunction = caret::multiClassSummary,
                            savePredictions = "final")

mtry_grid <- data.frame(mtry = unique(pmax(1, round(c(0.5, 1, 2, 4) * sqrt(ncol(train_x))))))

set.seed(42)

rf_Mval <- caret::train(x = train_x, 
                        y = histo, 
                        method = "rf",
                        metric = "logLoss",                        
                        trControl = ctrl, 
                        tuneGrid = mtry_grid, 
                        ntree = 500, 
                        importance = TRUE)


saveRDS(rf_Mval, file.path(path.expand(out_dir), "RF_DMP_Mval.rds"))



#-------- Test on Test NMB and Cavalli -----------

test_sets <- list(
  "Test NMB" = list(mval = test_NMB_Mval, ids = NMB_test_pheno$Sample_ID, labels = NMB_test_pheno$CPR_Histology),
  "Cavalli"  = list(mval = Cavalli_Mval,  ids = Cavalli_pheno$Sample_ID,  labels = Cavalli_pheno$histology))

test_results <- list()

for (nm in names(test_sets)) {
  ts <- test_sets[[nm]]
  
  test_x      <- as.data.frame(t(ts$mval[feature_cpgs, , drop = FALSE]))
  true_labels <- factor(ts$labels, levels = levels(histo))                        
  stopifnot(identical(rownames(test_x), as.character(ts$ids)))
  
  message("Testing the M-value Random Forest on ", nm, " (", nrow(test_x), " samples)")
  
  res <- tryCatch(
    test_and_evaluate(
      model = rf_Mval,
      test_data = test_x,
      true_test_labels = true_labels,
      test_data_name = nm,
      model_name = "Random Forest",
      method = "caret",
      file_id = paste0("RF_DMP_Mval_", gsub(" ", "_", nm)),
      feature_selection_method = paste0("top ", length(feature_cpgs), " DMPs (M-values)"),
      balanced_method = NULL,
      save_in = file.path(path.expand(out_dir), "evaluation")),
    error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL })
  
  if (is.null(res)) next
  res$summary$test_set  <- nm
  res$summary$input     <- "M-values"
  res$classwise$test_set <- nm
  test_results[[nm]] <- res
}

final_summary   <- do.call(rbind, lapply(test_results, function(r) r$summary))
final_classwise <- do.call(rbind, lapply(test_results, function(r) r$classwise))
print(final_summary)

write.csv(final_summary,   file.path(path.expand(out_dir), "RF_DMP_Mval_summary_test_sets.csv"),   row.names = FALSE)

write.csv(final_classwise, file.path(path.expand(out_dir), "RF_DMP_Mval_classwise_test_sets.csv"), row.names = FALSE)

