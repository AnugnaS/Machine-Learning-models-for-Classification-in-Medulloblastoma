
# ---- Custom Functions to test the ML models and produce corresponding evaluation metrics and plots ------------------------------------------

# for title-wrapping to avoid the plots getting cut off over long descriptive titles.

if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")

wrap_title <- function(x, width = 55) { str_wrap(x, width = width) } 


# Wrapper function to test and evaulate models trained with both caret and methylClass (Liu,2024)
# Evaluation metrics (overall and/or classwise)- Confusion Matrix, accuracy, balanced accuracy,brier scores,log loss, sensitivity,specificity,precision, recall, MCC, 
# ROC-AUC, PR-AUC,F1 score and kappa scores
# Plots for Confusion Matrix, ROC and PR Curves. 

# IMPORTANT: Change Train set and Test set names,class levels and labels/titles accordingly

test_and_evaluate <- function(model, 
                              method,
                              test_data, 
                              true_test_labels, 
                              test_data_name, 
                              model_name,
                              file_id, 
                              train_set_name,
                              feature_selection_method, 
                              balanced_method, 
                              positive_class,
                              calibrated = FALSE, 
                              calibration_method, 
                              save_in) {
  
  
  # model : model object
  
  # method: 'caret' or 'methylClass' ; "character"
  
  # test_data : test data beta values ; "matrix" or "data.frame"
  
  # true_test_labels : Given labels before prediction ; "vector"
  
  # test_data_name (To set labels of plots); "character"
  
  # model_name (To set labels of plots): "Random Forest", "Support Vector Machine", "XGBoost" ; "character"
  
  # file_id ( To set the names of the saved files) ; "character" 
  
  # train_set_name (To set labels of plots and names of the saved files) : "Train NMB" ; "character"
  
  # feature_selection_method : "10k", "LASSO", "all_probes"; "character
  
  # balanced_method : "Under-sampled", "Over-sampled", "SMOTE", "ADASYN", "Tomek", "SMOTE-Tomek"; "character"
  
  # positive class (for binary models only); "character", case.sensitive
  
  # calibrated (To calibrated prediction probabilities) : TRUE
  
  # calibration_method: "LR","FLR" and/or "MR" (applicable for method = "methylClass") ; "character"
  
  # save_in (directory to save the outputs in) ; "character"
  
  
  dir.create(save_in, recursive = TRUE, showWarnings = FALSE)
  
  # ---- Names used in file names and plot titles ------------------------------------------
  
  details    <- c(feature_selection_method, balanced_method)
  details    <- details[!is.na(details) & nzchar(details)]
  detail_str <- if (length(details) > 0) paste0(" (", paste(details, collapse = ", "), ")") else ""
  title_core <- paste0(model_name, " trained on ", train_set_name, detail_str,
                       ", tested on ", test_data_name)
  
  
  if (is.null(file_id)) {
    file_id <- paste(c(model_name, train_set_name, feature_selection_method, balanced_method, test_data_name),
                     collapse = "_")
  }
  
  file_id  <- gsub("[^A-Za-z0-9_]+", "_", file_id)
  out_file <- function(suffix) file.path(save_in, paste0(file_id, "_", suffix))
  
  # ---- Packages -----------------------------------------------------------------------------
  
  required_pgs <- c("caret", "pROC", "ggplot2", "tibble", "stringr", "MLmetrics", "dplyr", "scales")
  missing_pgs  <- required_pgs[!sapply(required_pgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pgs) > 0) {
    message("Installing missing packages for evaluation methods: ", paste(missing_pgs, collapse = ", "))
    install.packages(missing_pgs, dependencies = TRUE)
  }
  invisible(lapply(required_pgs, function(pkg) library(pkg, character.only = TRUE)))
  
  # ---- Input checks ---------------------------------------------------------------------------
  
  if (!is.data.frame(test_data) && !is.matrix(test_data)) test_data <- as.data.frame(test_data)
  if (!is.factor(true_test_labels)) true_test_labels <- as.factor(true_test_labels)
  if (nrow(test_data) != length(true_test_labels)) {
    stop("Check the match between the rownames of the test data and of the ground truth labels")
  }
  
  classes <- levels(true_test_labels)
  
  detect_engine <- function(model) if (inherits(model, "train")) "caret" else "methylClass"
  
  if (is.null(method)) method <- detect_engine(model)
  
  model_labels <- if (method == "caret") model$levels else levels(as.factor(model$mod[[1]]$classes))
  
  if (!setequal(model_labels, classes)) {
    stop("The histology levels of the model (", paste(model_labels, collapse = ", "),
         ") and the test data (", paste(classes, collapse = ", "), ") DO NOT match")
  }
  
  # ---- Predictions ------------------------------------------------------------------------------
  
  if (method == "caret") {
    pred_labels <- tryCatch(predict(model, newdata = test_data),
                            error = function(e) stop("predict() failed: ", conditionMessage(e)))
    pred_probs  <- tryCatch(predict(model, newdata = test_data, type = "prob"),
                            error = function(e) stop("predict(prob) failed: ", conditionMessage(e)))
    pred_probs  <- as.data.frame(pred_probs)
    
  } else if (calibrated) {
    calfit <- switch(as.character(calibration_method),
                     MR  = model$glmnet.calfit,
                     LR  = model$platt.calfits,
                     FLR = model$platt.brglm.calfits,
                     stop("calibration_method must be one of 'MR', 'LR' or 'FLR' when calibrated = TRUE"))
    
    if (is.null(calfit)) stop("This model has no '", calibration_method,
                              "' calibration fit (was calibrationmethod passed to maintrain()?)")
    pred <- tryCatch(methylClass::mainpredict(newdat = test_data, mod = model$mod, calibratemod = calfit),
                     error = function(e) stop("methylClass::predict error: ", conditionMessage(e)))
    
    pred_probs  <- as.data.frame(pred$probs)
    
    pred_labels <- as.factor(pred$pres)
    
  } else {
    pred <- tryCatch(methylClass::mainpredict(newdat = test_data, mod = model$mod),
                     error = function(e) stop("methylClass::predict error: ", conditionMessage(e)))
    pred_labels <- as.factor(pred$rawpres)
    pred_probs  <- as.data.frame(pred$scores)
  }
  
  
  pred_labels <- factor(pred_labels, levels = classes)
  
  # ---- Confusion matrix + plot ---------------------------------------------------------------------
  
  cm <- tryCatch(
    caret::confusionMatrix(pred_labels, true_test_labels, mode = "everything", positive = positive_class),
    error = function(e) stop("confusion matrix encountered error: ", conditionMessage(e)))
  
  cm_df <- as.data.frame(cm$table)
  
  levls <- levels(as.factor(true_test_labels))
  
  
  molgrp_pattern    <- "MB$"                                         # SHH-MB, WNT-MB, Group3-MB ...
  subgroup_pattern  <- "^(SHH[_-][0-9]+[A-C]?|G34[_-][IVX]+)$"       # SHH_1, SHH_3A, G34_I ... G34_VIII
  histology_pattern <- "^(DN|CLA|LCA|MBEN|non[_-]DN)$"               # histology labels
  
  # IMPORTANT: CHANGE LABELS ACCORDINGLY
  
  COLOURS <- NULL
  
  if (any(grepl(molgrp_pattern, levls, ignore.case = TRUE))) {
    COLOURS <- c("WNT-MB" = "#4DBBD5", "SHH-MB" = "#E64B35", "Group3-MB" = "yellow", "Group4-MB" = "green")
    
  } else if (any(grepl(subgroup_pattern, levls, ignore.case = TRUE))) {
    
    COLOURS <- c("SHH_1" = "#FFB6C1", "SHH_2" = "#FF69B4", "SHH_3" = "#C71585", "SHH_4" = "#FF1493",
                 "SHH_3A" = "#9B1B6E", "SHH_3B" = "#D4449A", "SHH_3C" = "#F0A8D0",
                 "G34_I" = "#800080", "G34_II" = "#C71585", "G34_III" = "#FF8C00", "G34_IV" = "#FFD700",
                 "G34_V" = "#ADFF2F", "G34_VI" = "#90EE90", "G34_VII" = "#ADD8E6", "G34_VIII" = "#006400",
                 "WNT" = "darkblue", "MBNOS" = "black", "NA" = "grey")
    
  } else if (any(grepl(histology_pattern, levls, ignore.case = TRUE))) {
    
    COLOURS <- c("DN" = "darkgreen", "non_DN" = "purple", "non-DN" = "purple",
                 "CLA" = "purple", "LCA" = "pink", "MBEN" = "green")
  }
  
  
  
  cm_plot <- ggplot(cm_df, aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = Reference, alpha = Freq), color = "white") +
    geom_text(aes(label = Freq), size = 4,
              color = ifelse(cm_df$Freq > max(cm_df$Freq) * 0.5, "white", "black"), fontface = "bold") +
    scale_fill_manual(values = COLOURS, guide = "none") +
    scale_alpha_continuous(range = c(0.15, 1), guide = "none") +
    labs(title = wrap_title(paste0("Confusion Matrix for ", title_core)),
         x = "True Class", y = "Predicted Class") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          axis.title = element_text(hjust = 0.5, face = "bold"),
          plot.margin = ggplot2::margin(t = 15, r = 15, b = 15, l = 15)) +
    coord_cartesian(clip = "off")
  
  ggsave(out_file("confusion_matrix.png"), cm_plot, width = 9, height = 6, dpi = 300, bg = "white")
  
  # ---- Accuracy, kappa ---------------------------------------------------------------------------------
  
  accuracy <- unname(cm$overall["Accuracy"])
  kappa    <- unname(cm$overall["Kappa"])
  
  # ---- Brier + log-loss  -----------------------------------------
  
  P  <- as.matrix(pred_probs[, classes, drop = FALSE])
  storage.mode(P) <- "double"
  rs <- rowSums(P)
  P  <- P / ifelse(is.finite(rs) & rs > 0, rs, NA_real_)
  lab     <- as.character(true_test_labels)
  Y       <- outer(lab, classes, "==") * 1
  brier   <- mean(rowSums((P - Y)^2), na.rm = TRUE)
  p_true  <- P[cbind(seq_len(nrow(P)), match(lab, classes))]
  logloss <- -mean(log(pmin(pmax(p_true, 1e-15), 1)), na.rm = TRUE)
  
  # ---- Sensitivity / specificity ----------------------------------------------------------------------------
  
  byclass <- cm$byClass
  if (is.null(dim(byclass))) {     
    
    positive_cls <- cm$positive
    other_class  <- setdiff(classes, positive_cls)
    sensitivity <- setNames(c(unname(byclass["Sensitivity"]), unname(byclass["Specificity"])),
                            c(positive_cls, other_class))
    specificity <- setNames(c(unname(byclass["Specificity"]), unname(byclass["Sensitivity"])),
                            c(positive_cls, other_class))
  } else {                                           
    sensitivity <- byclass[, "Sensitivity"]
    specificity <- byclass[, "Specificity"]
    names(sensitivity) <- gsub("^Class: ", "", rownames(byclass))
    names(specificity) <- gsub("^Class: ", "", rownames(byclass))
  }
  
  # balanced accuracy 
  
  balanced_accuracy <- mean(sensitivity, na.rm = TRUE)
  
  # ---- Matthews correlation coefficient  -------------------
  
  tab     <- as.matrix(cm$table)[classes, classes, drop = FALSE]
  s       <- sum(tab)
  correct <- sum(diag(tab))
  row_tot <- rowSums(tab)
  col_tot <- colSums(tab)
  mcc_den <- sqrt((s^2 - sum(col_tot^2)) * (s^2 - sum(row_tot^2)))
  mcc     <- if (is.finite(mcc_den) && mcc_den > 0) (correct * s - sum(row_tot * col_tot)) / mcc_den else 0
  
  # ---- Precision, recall, F1 ------------------------------------------------------------------------------------
  
  stopifnot(identical(names(dimnames(cm$table)), c("Prediction", "Reference")))
  cm_table <- t(cm$table)                           
  
  classwise_precision <- sapply(classes, function(cls) {
    tp <- cm_table[cls, cls]; fp <- sum(cm_table[, cls]) - tp; tp / (tp + fp)
  })
  classwise_recall <- sapply(classes, function(cls) {
    tp <- cm_table[cls, cls]; fn <- sum(cm_table[cls, ]) - tp; tp / (tp + fn)
  })
  names(classwise_precision) <- classes
  names(classwise_recall)    <- classes
  classwise_f1 <- 2 * (classwise_precision * classwise_recall) / (classwise_precision + classwise_recall)
  names(classwise_f1) <- classes
  
  precision <- mean(classwise_precision, na.rm = TRUE)   
  recall    <- mean(classwise_recall,    na.rm = TRUE)
  f1        <- mean(classwise_f1,        na.rm = TRUE)
  
  # ---- ROC / AUC ------------------------------------------------------------------------------------------------
  
  roc_list <- tryCatch({
    setNames(lapply(classes, function(cls) {
      pROC::roc(response = as.numeric(true_test_labels == cls), predictor = pred_probs[[cls]], quiet = TRUE)
    }), classes)
  }, error = function(e) stop("ROC Curve error: ", conditionMessage(e)))
  
  class_auc <- sapply(roc_list, function(r) as.numeric(pROC::auc(r)))
  
  handtill_obj <- tryCatch(
    pROC::multiclass.roc(response = true_test_labels, predictor = pred_probs[, classes, drop = FALSE], quiet = TRUE),
    error = function(e) stop("Hand-Till AUC Error: ", conditionMessage(e)))
  macro_auc <- as.numeric(handtill_obj$auc)
  
  # DeLong 95% CI for the AUC (2 classes only)
  
  macro_auc_ci <- c(NA_real_, NA_real_)
  if (length(classes) == 2 && !is.null(roc_list)) {
    ci <- tryCatch(as.numeric(pROC::ci.auc(roc_list[[classes[1]]], method = "delong")), error = function(e) NULL)
    if (!is.null(ci) && length(ci) == 3) macro_auc_ci <- ci[c(1, 3)]
  }
  
  roc_df <- do.call(rbind, lapply(names(roc_list), function(cls) {
    r <- roc_list[[cls]]
    data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities, class = cls,
               class_label = paste0(cls, " (AUC = ", round(pROC::auc(r), 3), ")"))
  }))
  
  roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = class)) +
    geom_line(linewidth = 1) +
    geom_abline(linetype = "dashed", color = "grey50") +
    scale_color_manual(values = COLOURS,
                       labels = setNames(roc_df$class_label, roc_df$class)[unique(roc_df$class)]) +
    labs(title = wrap_title(paste0("ROC Curves for ", title_core)),
         subtitle = paste0("Macro AUC (Hand-Till) = ", round(macro_auc, 3)),
         x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.margin = ggplot2::margin(t = 15, r = 15, b = 15, l = 15),
          axis.title = element_text(hjust = 0.5, face = "bold")) +
    coord_cartesian(clip = "off")
  
  ggsave(out_file("roc_curves.png"), roc_plot, width = 9, height = 6, dpi = 300, bg = "white")
  
  # ---- Precision-recall curves / PR-AUC ----------------------------------------------------------------------------
  
  pr_list <- tryCatch({
    setNames(lapply(classes, function(cls) {
      prs <- pROC::coords(roc_list[[cls]], x = "all", ret = c("recall", "precision"), transpose = FALSE)
      prs <- prs[!is.na(prs$precision) & !is.na(prs$recall), ]
      prs[order(prs$recall), ]
    }), classes)
  }, error = function(e) stop("PR Curve Error: ", conditionMessage(e)))
  
  pr_auc <- tryCatch({
    sapply(classes, function(cls) {
      MLmetrics::PRAUC(y_pred = pred_probs[[cls]], y_true = as.numeric(true_test_labels == cls))
    })
  }, error = function(e) stop("PRAUC Error: ", conditionMessage(e)))
  names(pr_auc) <- classes
  
  pr_df <- do.call(rbind, lapply(classes, function(cls) {
    d <- pr_list[[cls]]
    data.frame(recall = d$recall, precision = d$precision, class = cls,
               class_label = paste0(cls, " (PR-AUC =", round(pr_auc[[cls]], 3), ")"))
  }))
  
  macro_pr_auc <- mean(pr_auc, na.rm = TRUE)
  
  pr_plot <- ggplot(pr_df, aes(x = recall, y = precision, color = class)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = COLOURS,
                       labels = setNames(pr_df$class_label, pr_df$class)[unique(pr_df$class)]) +
    labs(title = wrap_title(paste0("Precision-Recall Curves for ", title_core)),
         subtitle = paste0("Macro PR-AUC = ", round(macro_pr_auc, 3)),
         x = "Recall", y = "Precision", color = "Class") +
    theme_minimal() + ylim(0, 1) + xlim(0, 1) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(hjust = 0.5, face = "bold"),
          plot.margin = ggplot2::margin(t = 15, r = 15, b = 15, l = 15)) +
    coord_cartesian(clip = "off")
  
  ggsave(out_file("pr_curves.png"), pr_plot, width = 9, height = 6, dpi = 300, bg = "white")
  
  # ---- Per-sample probabilities -----------------------------------------------------------------------------------------
  
  pred_probs$pred <- as.factor(pred_labels)
  pred_probs %>%
    tibble::rownames_to_column(var = "Sample_ID") %>%
    write.csv(out_file("probabilities.csv"), row.names = FALSE)
  
  
  # ---- Summary tables -------------------------------------------------------------------------------------------------------
  
  meta <- data.frame(model = model_name, train_set = train_set_name, test_data = test_data_name,
                     stringsAsFactors = FALSE)
  if (!is.null(feature_selection_method)) meta$feature_selection_method <- feature_selection_method
  if (!is.null(balanced_method))         meta$balanced_method         <- balanced_method
  
  summary_df <- cbind(meta, data.frame(
    accuracy           = accuracy,
    balanced_accuracy  = balanced_accuracy,
    macro_precision    = precision,
    macro_recall       = recall,
    macro_f1           = f1,
    mcc                = mcc,
    kappa              = kappa,
    brier              = brier,
    log_loss           = logloss,
    macro_auc          = macro_auc,
    macro_auc_ci_lower = macro_auc_ci[1],
    macro_auc_ci_upper = macro_auc_ci[2],
    macro_pr_auc       = macro_pr_auc))
  
  # Binary only
  
  if (!is.null(positive_class)) {
    summary_df$positive_class <- positive_class
    summary_df$sensitivity    <- unname(sensitivity[positive_class])
    summary_df$specificity    <- unname(specificity[positive_class])
    summary_df$precision      <- unname(classwise_precision[positive_class])
    summary_df$f1             <- unname(classwise_f1[positive_class])
    summary_df$auc            <- unname(class_auc[positive_class])
    summary_df$pr_auc         <- unname(pr_auc[positive_class])
  }
  
  summary_df <- summary_df[, !duplicated(names(summary_df)), drop = FALSE]
  
  classwise_meta <- meta[rep(1, length(classes)), , drop = FALSE]
  
  rownames(classwise_meta) <- NULL
  
  classwise_df <- cbind(classwise_meta, data.frame(
    class       = classes,
    sensitivity = as.numeric(sensitivity[classes]),
    specificity = as.numeric(specificity[classes]),
    precision   = as.numeric(classwise_precision[classes]),
    recall      = as.numeric(classwise_recall[classes]),
    f1          = as.numeric(classwise_f1[classes]),
    auc         = as.numeric(class_auc[classes]),
    pr_auc      = as.numeric(pr_auc[classes]),
    youden_j    = as.numeric(sensitivity[classes]) + as.numeric(specificity[classes]) - 1))
  
  
  classwise_df <- classwise_df[, !duplicated(names(classwise_df)), drop = FALSE]
  
  write.csv(summary_df,   out_file("summary_metrics.csv"),   row.names = FALSE)
  write.csv(classwise_df, out_file("classwise_metrics.csv"), row.names = FALSE)
  
  list(confusion_matrix = cm, cm_plot = cm_plot,
       roc_list = roc_list, roc_plot = roc_plot,
       pr_list = pr_list, pr_plot = pr_plot,
       class_auc = class_auc, pr_auc = pr_auc, macro_pr_auc = macro_pr_auc,
       sensitivity = sensitivity, specificity = specificity,
       classwise_precision = classwise_precision, classwise_recall = classwise_recall,
       classwise_f1 = classwise_f1, kappa = kappa,
       macro_auc = macro_auc, handtill_obj = handtill_obj,
       summary = summary_df, classwise = classwise_df, pred_probs = pred_probs)
} 



# Wrapper function to produce probability score heatmaps with the output from test_and_evaluate() 
# and to also produce complementing barplots/scatter plots/density plots for
# probability score distribution/correlations in clinical correlate groups

# IMPORTANT: Change Train set and Test set names, and labels/titles accordingly

prob_heatmaps <- function(prob_df,
                          pred, 
                          test_pheno, 
                          sample_id_column, 
                          heatmap_label,
                          true_test_labels, 
                          mol_group, 
                          sub_group, 
                          MYC_status, 
                          MYCN_status, 
                          OS, 
                          follow_up, 
                          Age,
                          continuum = "NULL", 
                          heatmap_panel = TRUE, 
                          correlate_plot = NULL,
                          column_for_panel = "DN", 
                          file_id = NULL, 
                          save_in ) { 
  
  # prob_df: probability score matrix from test_and_evaluate()$pred_probs ; "vector-data.frame"
  
  # pred: prediction labels from test_and_evaluate()$pred_probs$pred , "vector"
  
  # test_pheno : Phenotype table for test data ; "data.frame"
  
  # sample_id_column : Sample IDs of test data from test_pheno; "character"
  
  # heatmap_label : Title of the output heatmap; "character" 
  
  # true_test_labels : Given labels of the test data before prediction from test_pheno; "vector"
  
  # mol_group : Principal molecular group column from test_pheno;"character"
  
  # sub_group : molecular sub group column from test_pheno; "character"
  
  # MYC_status : MYC Amplification column from test_pheno ; "character"
  
  # MYCN_status :  MYCN Amplification column from test_pheno ; "character"
  
  # OS : Overall Survival (Event : 1, No Event : 0) column from test_pheno ; "character"
  
  # follow_up : test_pheno column for Follow up time for Overall Survival; "character"
  
  # Age : Chronological Age column from test_pheno ; "character"
  
  # continuum : To choose which level of true_test_label would be plotted for its probability score enrichment across all the test samples ; "character"
  
  # heatmap_panel : TRUE if other clinical correlate should be plotted against probability scores and displayed next to probability heatmap. 
  
  # correlate_plot : To choose with clincal correlate plot to be saved, along with saving standalone probability score heatmap: "Heatmap", "Age", "MYC/MYCN", "Mol/Subgroup", "Survival" ; "character"
  
  # column_for_panel : level of true_test_label whose probability score is to be plotted in the correlate plots "DN", "non-DN", "CLA","LCA", "MBEN"; "character"
  
  # file_id : file name prefix for the outputs to be saved with; "character"
  
  dir.create(save_in, recursive = TRUE, showWarnings = FALSE)
  
  # ---- Packages -----------------------------------------------------------------------------
  
  required_pgs <- c("ggplot2", "ComplexHeatmap", "circlize", "cowplot", "ggpubr", "survival", "dplyr", "survminer", "tibble")
  
  bioc_pgs <- c("ComplexHeatmap")
  
  missing_pgs <- required_pgs[!sapply(required_pgs, requireNamespace, quietly = TRUE)]
  
  if (length(missing_pgs) > 0) {
    missing_cran <- setdiff(missing_pgs, bioc_pgs)
    missing_bioc <- intersect(missing_pgs, bioc_pgs)
    
    if (length(missing_cran) > 0) {
      message("Installing missing CRAN packages for Heatmaps: ", paste(missing_cran, collapse = ", "))
      install.packages(missing_cran, dependencies = TRUE)
    }
    
    if (length(missing_bioc) > 0) {
      message("Installing missing Bioconductor packages for Heatmaps: ", paste(missing_bioc, collapse = ", "))
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
    }
  }
  
  invisible(lapply(required_pgs, function(pkg) { 
    library(pkg, character.only = TRUE)
  }))
  
  
  
  # IMPORTANT: CHANGE LABELS ACCORDINGLY
  
  heatmap_label_display <- wrap_title(heatmap_label, width = 70)
  
  if (all(rownames(prob_df) == test_pheno[[sample_id_column]])) { 
    if (continuum == "DN/multi") {
      prob_df$DN_continuum <- prob_df$DN / (prob_df$CLA + prob_df$DN)
      continuum_score  <- prob_df$DN_continuum
      range_label <- "DN_cont"
    } else if (continuum == "LCA") {
      prob_df$LCA_continuum <- prob_df$LCA / (prob_df$CLA + prob_df$LCA)
      continuum_score  <- prob_df$LCA_continuum
      range_label <- "LCA_cont"
    } else if (continuum == "CLA") {
      prob_df$CLA_continuum <- prob_df$CLA / (prob_df$CLA + prob_df$LCA)
      continuum_score  <- prob_df$CLA_continuum
      range_label <- "CLA_cont"
    } else if (continuum == "DN/bi") {
      prob_df$DN_continuum <- prob_df$DN / (prob_df$non_DN + prob_df$DN)
      continuum_score  <- prob_df$DN_continuum
      range_label <- "DN_cont"
    } else {
      stop("Unknown continuum '", continuum, "'. Use one of: DN/bi, DN/multi, LCA, CLA")
    }
  } else {
    stop("Probability score rownames don't match Sample IDs in the Phenotype table")
  }
  
  # ---- Prediction Class Labels -----------------------------------------------------------------------------
  
  if (all(rownames(prob_df) == test_pheno[[sample_id_column]])) {
    
    prob_df$True_label <- as.factor(test_pheno[[true_test_labels]])
    
  } else {
    
    test_pheno<-test_pheno[match(rownames(prob_df),test_pheno[[sample_id_column]]),]
    
    prob_df$True_label <- as.factor(test_pheno[[true_test_labels]])
  }
  
  if(binary){
    
    prob_df <- prob_df %>% mutate(bi_histo = case_when( 
      True_label %in% c("DN", "MBEN") ~ "DN",
      True_label %in% c("CLA", "LCA") ~ "non_DN"
    ))
    
    prob_df$bi_histo <- as.factor(prob_df$bi_histo)
    
  }
  
  classify_outcome <- function(true_label, pred) {
    truth <- as.character(true_label)
    pred  <- as.character(pred)
    ifelse(truth == pred, paste0("Correct_", truth), paste0(truth, "->", pred))
  }
  
  if(binary){
    
    prob_df$Class <- classify_outcome(prob_df$bi_histo, prob_df$pred, binary = binary)
    
  } else {
    
    prob_df$Class <- classify_outcome(prob_df$True_label, prob_df$pred, binary = binary) } 
  
  
  prob_cols   <- if (binary) c("DN", "non_DN") else intersect(c("CLA", "DN", "LCA", "MBEN"), colnames(prob_df))
  
  prob_matrix <- as.matrix(prob_df[, prob_cols, drop = FALSE])
  
  stopifnot(
    sample_id_column %in% colnames(test_pheno),
    nrow(prob_matrix) == nrow(test_pheno),
    all(rownames(prob_matrix) == test_pheno[[sample_id_column]])
  )
  
  
  # ---- Annotation Labels-----------------------------------------------------------------------------
  
  test_pheno[[MYC_status]]  <- gsub("0", "Neutral", test_pheno[[MYC_status]])
  test_pheno[[MYC_status]]  <- gsub("1", "AMP", test_pheno[[MYC_status]])
  test_pheno[[MYC_status]]  <- ifelse(test_pheno[[MYC_status]] %in% c("Neutral", "AMP"), test_pheno[[MYC_status]], "NoData")
  test_pheno[[MYC_status]]  <- as.factor(test_pheno[[MYC_status]])
  test_pheno[[MYC_status]]  <- droplevels(test_pheno[[MYC_status]])
  test_pheno[[MYCN_status]] <- gsub("0", "Neutral", test_pheno[[MYCN_status]])
  test_pheno[[MYCN_status]] <- gsub("1", "AMP", test_pheno[[MYCN_status]])
  test_pheno[[MYCN_status]] <- ifelse(test_pheno[[MYCN_status]] %in% c("Neutral", "AMP"), test_pheno[[MYCN_status]], "NoData")
  test_pheno[[MYCN_status]] <- as.factor(test_pheno[[MYCN_status]])
  test_pheno[[MYCN_status]] <- droplevels(test_pheno[[MYCN_status]])
  
  test_pheno[[OS]] <- ifelse(as.character(test_pheno[[OS]]) %in% c("0", "1"), as.character(test_pheno[[OS]]), "NoData")
  test_pheno[[OS]] <- factor(test_pheno[[OS]], levels = c("0", "1", "NoData"))
  
  HISTO_COLOURS<-c("CLA" = "purple2", "LCA" = "pink",
                   "DN"  = "darkgreen", "MBEN" = "lightgreen")
  
  PRED_COLOURS<-c("non_DN" = "purple", "DN" = "darkgreen", "CLA" = "purple2", "LCA" = "pink", "MBEN" = "lightgreen")
  
  CLASS_COLOURS<-c("Correct_DN"    = "darkgreen", "Correct_non_DN" = "purple4",
                   "non_DN->DN"    = "pink3",     "DN->non_DN"     = "orange", "Correct_CLA" = "purple3", "Correct_LCA" = "pink",
                   "CLA->LCA" = "#D4449A", "LCA->CLA" = "darkblue", "CLA->DN" = "green4","DN->CLA" = "orange",
                   "LCA->DN" = "pink3","DN->LCA" = "#9B1B6E", "Correct_MBEN" = "lightgreen", "DN->MBEN" = "#90EE90", "LCA->MBEN" = "#FF8C00",
                   "CLA->MBEN" = "#FFD700" , "MBEN->CLA" = "#ADFF2F", "MBEN->LCA" = "orange", "MBEN->DN" = "green3" )
  
  MOL_GRP_COLOURS<-c("WNT-MB" = "#4DBBD5", "SHH-Inf-MB" = "#E64B35",
                     "SHH-Old-MB" = "orange", "SHH-MB" = "#E64B35",
                     "Group3-MB" = "yellow", "Group4-MB" = "green")
  
  SUB_GRP_COLOURS<- c(
    "SHH_1" = "#FFB6C1", "SHH_2" = "#FF69B4", "SHH_3" = "#C71585", "SHH_4" = "#FF1493",
    "SHH_3A" = "#9B1B6E", "SHH_3B" = "#D4449A", "SHH_3C" = "#F0A8D0",
    "G34_I"   = "#800080", "G34_II"  = "#C71585", "G34_III" = "#FF8C00", "G34_IV" = "#FFD700",
    "G34_V"   = "#ADFF2F", "G34_VI"  = "#90EE90", "G34_VII" = "#ADD8E6", "G34_VIII" = "#006400",
    "WNT" = "darkblue", "MBNOS" = "black", "NA" = "grey"
  )
  
  
  OS_COLOURS <- c("0" = "grey70", "1" = "red", "NoData" = "grey40")
  
  prob_df$True_label  <- factor(as.character(prob_df$True_label), 
                                levels = intersect(names(HISTO_COLOURS), 
                                                   unique(as.character(prob_df$True_label))))
  
  prob_df$pred        <- factor(as.character(prob_df$pred), 
                                levels = intersect(names(PRED_COLOURS), 
                                                   unique(as.character(prob_df$pred))))
  prob_df$Class      <- factor(prob_df$Class, 
                               levels = intersect(names(CLASS_COLOURS), 
                                                  unique(prob_df$Class)))
  
  if (binary) {class_display <- gsub("non_DN", "nonDN", as.character(prob_df$Class))}
  
  prob_df$Age <- test_pheno[[Age]] 
  
  prob_df <- prob_df %>%
    mutate(age_bin = case_when(
      Age >= 0  & Age < 5   ~ "0-4",
      Age >= 5  & Age <= 10 ~ "5-10",
      Age > 10  & Age <= 15 ~ "11-15",
      Age > 15              ~ "15+"
    ))
  
  col_fun <- colorRamp2(
    c(0, 0.25, 0.5, 0.75, 1),
    c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B")
  )
  
  col_ha <- HeatmapAnnotation(
    Pred       = prob_df$pred,
    true_label = test_pheno[[true_test_labels]],
    mol_grp    = test_pheno[[mol_group]],
    sub_grp    = test_pheno[[sub_group]],
    MYC        = test_pheno[[MYC_status]],
    MYCN       = test_pheno[[MYCN_status]],
    OS         = test_pheno[[OS]],
    Range      = continuum_score,
    Age_bin    = factor(prob_df$age_bin, levels = c("0-4", "5-10", "11-15", "15+")),
    
    annotation_label = c(
      Pred       = "Pred",
      true_label = "True label",
      mol_grp    = "Mol.group",
      sub_grp    = "Sub group",
      MYC        = "MYC",
      MYCN       = "MYCN",
      OS         = "OS",
      Range      = range_label,
      Age_bin    = "Age Range"
    ),
    
    col = list(
      true_label = HISTO_COLOURS,
      Pred       = PRED_COLOURS,
      mol_grp    = MOL_GRP_COLOURS,
      sub_grp    = SUB_GRP_COLOURS,
      
      MYC        = c("Neutral" = "white", "AMP" = "darkred",  "NoData" = "grey"),
      MYCN       = c("Neutral" = "white", "AMP" = "orange",   "NoData" = "grey"),
      OS         = c("0" = "white", "1" = "red", "NoData" = "grey"),
      
      Range = colorRamp2(
        breaks = c(0.21, 0.52, 0.83),
        colors = c("#FEF9E7", "#F7DC6F", "#1A5276")),
      
      Age_bin = c(
        "0-4"   = "#2166AC",
        "5-10"  = "#92C5DE",
        "11-15" = "#F4A582",
        "15+"   = "#B2182B"
      )
    ),
    
    simple_anno_size     = unit(4, "mm"),
    annotation_name_gp   = gpar(fontsize = 9),
    annotation_name_side = "left",
    border               = FALSE,
    
    annotation_legend_param = list(
      Age_bin = list(title = "Age (years)",
                     at =  c("0-4", "5-10", "11-15", "15+"), 
                     labels =  c("0-4", "5-10", "11-15", "15+")),
      
      OS  = list(title = "Death", at = c("0", "1"), labels = c("No", "Yes")),
      
      MYC  = list(title  = "MYC",
                  at     = c("Neutral", "AMP", "NoData"),
                  labels = c("WT", "Amp", "No data")),
      
      MYCN = list(title  = "MYCN",
                  at     = c("Neutral", "AMP", "NoData"),
                  labels = c("WT", "Amp", "No data")),
      
      Range = list(title  = range_label,
                   at     = c(0.21, 0.52, 0.83),
                   labels = c("0.21", "0.52", "0.83")),
      
      border = FALSE,
      gp = gpar(col = "black", lwd = 0.5))
  )
  
  stopifnot(
    nrow(prob_matrix) == nrow(prob_df),
    all(rownames(prob_matrix) == test_pheno[[sample_id_column]])
  )
  
  # ---- Heatmap -----------------------------------------------------------------------------
  
  ht <- Heatmap(
    t(prob_matrix),                     
    
    name      = "Probability",
    col       = col_fun,
    
    column_order = order(continuum_score),
    column_split = class_display,
    column_title_rot = if (binary) 0 else 90,        
    column_title_gp  = gpar(fontsize = 6, fontface = "bold"),
    column_gap       = unit(3, "mm"),
    
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    show_column_dend  = FALSE,
    
    top_annotation = col_ha,
    
    height = unit(1.5 * length(prob_cols), "cm"), 
    
    show_column_names = FALSE,
    show_row_names    = TRUE,
    row_labels        = paste0("P(", rownames(t(prob_matrix)), ")"),
    row_names_side    = "left",
    
    row_title    = NULL,
    rect_gp      = gpar(col = "grey85", lwd = 0.3),
    border       = TRUE,
    
    heatmap_legend_param = list(
      title      = "Probability",
      direction  = "vertical",
      title_gp   = gpar(fontsize = 9, fontface = "bold"),
      labels_gp  = gpar(fontsize = 8),
      grid_width = unit(4, "mm")
    )
  )
  
  
  ht_grob <- grid.grabExpr(
    draw(ht,
         column_title            = heatmap_label_display,
         column_title_gp         = gpar(fontsize = 12, fontface = "bold"),
         heatmap_legend_side     = "right",
         annotation_legend_side  = "right",
         merge_legend            = TRUE,
         padding                 = unit(c(2.5, 2.5, 2.5, 2.5), "mm")),
    width  = 15,
    height = 5 
  )
  
  ht_titled <- ht_grob
  
  right_theme <- theme_bw() +
    theme(
      plot.margin     = unit(c(0.5, 1, 0.5, 0.5), "mm"),
      legend.key.size = unit(3, "mm"),
      legend.text     = element_text(size = 6),
      legend.title    = element_text(size = 7, face = "bold"),
      axis.text       = element_text(size = 7),
      axis.title      = element_text(size = 6, hjust = 0.5, face = "bold"),
      plot.title           = element_text(hjust = 0.5, face = "bold") 
    )  
  
  prob_df$follow_up_time <- test_pheno[[follow_up]] 
  prob_df$event          <- as.numeric(as.character(test_pheno[[OS]]))
  
  class_palette      <- unname(CLASS_COLOURS[levels(prob_df$Class)])
  pred_palette       <- unname(PRED_COLOURS[levels(prob_df$pred)])
  true_label_palette <- unname(HISTO_COLOURS[levels(prob_df$True_label)])
  
  km_fit_class      <- survfit(Surv(follow_up_time, event) ~ Class,      data = prob_df)
  km_fit_pred       <- survfit(Surv(follow_up_time, event) ~ pred,       data = prob_df)
  km_fit_true_label <- survfit(Surv(follow_up_time, event) ~ True_label, data = prob_df)
  
  # ---- Inbuilt function for Kaplan Meier Curve Plot  -----------------------------------------------------------------------------
  
  build_km_plot <- function(fit, plot_title,palette = NULL, theme_to_apply = TRUE) {
    km <- ggsurvplot(
      fit,
      data       = prob_df,
      pval       = TRUE,
      conf.int   = FALSE,
      risk.table = FALSE,
      xlab       = "Follow-up time",
      title      = wrap_title(plot_title, width = 45),
      palette    = palette
    ) 
    
    km_plot <- km$plot + coord_cartesian(clip = "off") 
    
    if (theme_to_apply == TRUE) km_plot <- km_plot + theme (axis.title = element_text ( hjust = 0.5,face = "bold"),
                                                            plot.title = element_text(hjust = 0.5, face = "bold"))  
    
    km_plot     
  } 
  
  # ---- Correlate Plots within Heatmap Panel  -----------------------------------------------------------------------------
  
  if (heatmap_panel) {
    if (all(rownames(prob_df) == test_pheno[[sample_id_column]])) {
      
      prob_df$Age <- test_pheno[[Age]]
      
      p_age <- ggplot(prob_df, aes(x = Age, y = .data[[column_for_panel]], colour = Class)) +
        geom_point(size = 1, alpha = 0.7) +
        scale_x_continuous(limits = c(0, 45), breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45)) +
        scale_y_continuous(limits = c(0, 1),  breaks = c(0, 0.25, 0.5, 0.75, 1)) +
        scale_color_manual(values = CLASS_COLOURS) +
        stat_cor(method      = "pearson",
                 label.x     = 12,
                 label.y     = 1,
                 size        = 2.5, 
                 colour      = "black",
                 inherit.aes = FALSE,
                 data        = prob_df,
                 mapping     = aes(x = Age, y = .data[[column_for_panel]])) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
        labs(x = "Age (years)", y = paste0("P(", column_for_panel, ")"), colour = "Class") +
        right_theme 
      
      prob_df <- prob_df %>%
        mutate(age_bin = case_when(
          Age >= 0  & Age < 5   ~ "0-4",
          Age >= 5  & Age <= 10 ~ "5-10",
          Age > 10  & Age <= 15 ~ "11-15",
          Age > 15              ~ "15+"
        ))
      
      p_age_density <- ggplot(prob_df, aes(x = .data[[column_for_panel]], fill = age_bin)) +
        geom_density(alpha = 0.6) +
        scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
        scale_fill_manual(values = c(
          "0-4"   = "#2166AC",
          "5-10"  = "#92C5DE",
          "11-15" = "#F4A582",
          "15+"   = "#B2182B"
        ),
        breaks = c("0-4", "5-10", "11-15", "15+")) +
        geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey30") +
        labs(x = paste0("P(", column_for_panel, ")"), y = "Density", fill = "Age",
             title = wrap_title(paste0("P(", column_for_panel, ") ", "in age-ranges"))) +
        right_theme 
      
      prob_df$MYC_status <- as.factor(test_pheno[[MYC_status]])
      
      p_myc_density <- ggplot(prob_df,
                              aes(x = .data[[column_for_panel]], fill = MYC_status)) +
        geom_density(alpha = 0.6) +
        scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
        scale_fill_manual(values = c("AMP" = "orange", "Neutral" = "grey70")) +
        geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey30") +
        coord_flip(clip = "off") +
        labs(x = paste0("P(", column_for_panel, ")"), y = "Density", fill = "MYC",
             title = wrap_title(paste0("P(", column_for_panel, ") ", "by MYC Amplification Status"))) +
        right_theme 
      
      prob_df$MYCN_status <- as.factor(test_pheno[[MYCN_status]])
      
      p_mycn_density <- ggplot(prob_df,
                               aes(x = .data[[column_for_panel]], fill = MYCN_status)) +
        geom_density(alpha = 0.6) +
        scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
        scale_fill_manual(values = c("AMP" = "orange", "Neutral" = "grey70")) +
        geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey30") +
        coord_flip(clip = "off") +
        labs(x = paste0("P(", column_for_panel, ")"), y = "Density", fill = "MYCN",
             title = wrap_title(paste0("P(", column_for_panel, ") ", "by MYCN Amplification Status"))) +
        right_theme
      
      prob_df$OS <- as.factor(test_pheno[[OS]])
      
      p_os_density <- ggplot(prob_df,
                             aes(x = .data[[column_for_panel]], fill = OS)) +
        geom_density(alpha = 0.6) +
        scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
        scale_fill_manual(values = OS_COLOURS) +                      
        geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey30") +
        coord_flip(clip = "off") +
        labs(x = paste0("P(", column_for_panel, ")"), y = "Density", fill = "OS",
             title = wrap_title(paste0("P(", column_for_panel, ") ", "by Survival Group"))) +
        right_theme 
      
      prob_df$mol_group <- as.factor(test_pheno[[mol_group]])
      
      p_mol_group <- ggplot(prob_df,
                            aes(x = mol_group, y = .data[[column_for_panel]], fill = mol_group)) +
        geom_boxplot(outlier.size = 0.5, linewidth = 0.4, alpha = 0.8) +
        geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
        scale_fill_manual(values = MOL_GRP_COLOURS) +
        stat_compare_means(method = "kruskal.test",
                           label.y = 1,
                           size    = 2.5) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
        scale_y_continuous(limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
        labs(x = "Molecular group", y = paste0("P(", column_for_panel, ")"),
             title = wrap_title(paste0("P(", column_for_panel, ") ", "by Principal Molecular Group"))) + 
        theme(legend.position = "none",
              axis.text.x     = element_text(angle = 45, hjust = 0.5, size = 6),
              axis.title = element_text ( hjust = 0.5,face = "bold"),
              plot.title = element_text(hjust = 0.5, face = "bold")) + coord_cartesian(clip = "off") +
        right_theme 
      
      prob_df$sub_group <- as.factor(test_pheno[[sub_group]])
      
      p_subgroup <- ggplot(prob_df,
                           aes(x = sub_group, y = .data[[column_for_panel]], fill = sub_group)) +
        geom_boxplot(outlier.size = 0.5, linewidth = 0.4, alpha = 0.8) +
        geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
        scale_fill_manual(values = SUB_GRP_COLOURS) +
        stat_compare_means(method = "kruskal.test",
                           label.y = 1,
                           size    = 2.5) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
        scale_y_continuous(limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
        labs(x = "Sub group", y = paste0("P(", column_for_panel, ")"),
             title = wrap_title(paste0("P(", column_for_panel, ") ", "by Molecular Sub group"))) +
        theme(legend.position = "none",
              axis.text.x     = element_text(angle = 45, hjust = 0.5, size = 6),
              axis.title = element_text ( hjust = 0.5,face = "bold"),
              plot.title = element_text(hjust = 0.5, face = "bold")) +
        coord_cartesian(clip = "off") + 
        right_theme
      
      strip_margin <- theme(plot.margin = unit(c(0, 1, 0, 0), "mm"))
      
      p_age          <- p_age          + strip_margin
      p_age_density  <- p_age_density  + strip_margin
      p_myc_density  <- p_myc_density  + strip_margin
      p_mycn_density <- p_mycn_density + strip_margin
      p_os_density   <- p_os_density   + strip_margin
      p_mol_group    <- p_mol_group    + strip_margin
      p_subgroup     <- p_subgroup     + strip_margin
      
      right_panel <- plot_grid(
        p_age,
        p_age_density,
        p_myc_density,
        p_mycn_density,
        p_os_density,
        p_mol_group,
        p_subgroup,
        ncol  = 4,
        align = "hv",
        axis  = "tblr"
      )
      
      full_figure <- plot_grid(
        ht_titled,
        right_panel,
        ncol        = 1,
        rel_heights = c(0.45, 0.55)
      )
      
      pdf(file.path(save_in, paste0(file_id, "_.pdf")),
          width  = 15,
          height = 10,
          bg     = "white")
      
      print(full_figure)
      
      dev.off()
      
    } else {
      
      stop("Probability score rownames don't match Sample IDs in the Phenotype table to plot clinical correlates")
    }
    
  } else {
    
    ht_titled
    
    pdf(file.path(save_in, paste0(file_id, "_.pdf")),
        width  = 15,
        height = 10,
        bg     = "white")
    
    
  if (all(rownames(prob_df) == test_pheno[[sample_id_column]])) {
      
# ---- Correlate Plots to be returned individually  -----------------------------------------------------------------------------
      
    available_plots <- c("Heatmap", "Age", "MYC/MYCN", "Mol/Subgroup", "Survival")
      
    library(dplyr)
      
    cohort_name <- case_when(
      grepl("Test NMB", heatmap_label)  ~ "Test NMB",
      grepl("'NMB", heatmap_label) ~ "Train set",
      grepl("Cavalli", heatmap_label)   ~ "Cavalli",
      TRUE ~ NA_character_
      )
      
      if (any(is.na(cohort_name))) {
        print(unique(heatmap_label[is.na(cohort_name)]))
        warning("Unmatched heatmap_label values above - update the case_when mapping")
      }
      
      plots_to_run <- if (is.null(correlate_plot) || identical(correlate_plot, "All")) {
        
        available_plots
        
      } else {
        
        plots_to_run <- correlate_plot
      }
      
      unknown_plots <- setdiff(plots_to_run, available_plots)
      
      if (length(unknown_plots) > 0) {
        stop("Unknown correlate_plot value(s): ", paste(unknown_plots, collapse = ", "),
             ". Valid options are: ", paste(available_plots, collapse = ", "), ", or \"All\".")
      }
      
      ggsave(file.path(save_in, paste0(label_shorthand, "_heatmap.png")),
             ht_titled, width = 15, height = 5, dpi = 300, bg = "white")
      
      if ("Age" %in% plots_to_run) {
        
        prob_df$Age <- test_pheno[[Age]]
        
        p_age <- ggplot(prob_df, aes(x = Age, y = .data[[column_for_panel]], colour = Class)) +
          geom_point(size = 1, alpha = 0.7) +
          scale_x_continuous(limits = c(0, 45), breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45)) +
          scale_y_continuous(limits = c(0, 1),  breaks = c(0, 0.25, 0.5, 0.75, 1)) +
          scale_color_manual(values = CLASS_COLOURS) +
          stat_cor(method      = "pearson",
                   label.x     = quantile(prob_df$Age, 0.55, na.rm = TRUE),
                   label.y     = 1,
                   size        = 2.5,
                   colour      = "black",
                   inherit.aes = FALSE,
                   data        = prob_df,
                   mapping     = aes(x = Age, y = .data[[column_for_panel]])) +
          geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
          labs(x = "Age (years)", y = paste0("P(", column_for_panel, ")"), colour = "Class", 
               title = wrap_title(paste0("P(", column_for_panel, ") by Chronological Age in ", cohort_name))) +
          theme(axis.title = element_text ( hjust = 0.5,face = "bold"),
                plot.title = element_text(hjust = 0.5, face = "bold")) +
          coord_cartesian(clip = "off")
        
        ggsave(file.path(save_in, paste0(label_shorthand, "_age_scatter.png")),
               p_age, width = 9, height = 6, dpi = 300, bg = "white")
        
        prob_df <- prob_df %>%
          mutate(age_bin = case_when(
            Age >= 0  & Age < 5   ~ "0-4",
            Age >= 5  & Age <= 10 ~ "5-10",
            Age > 10  & Age <= 15 ~ "11-15",
            Age > 15              ~ "15+"
          ))
        
        p_age_density <- ggplot(prob_df, aes(x = .data[[column_for_panel]], fill = age_bin)) +
          geom_density(alpha = 0.6) +
          scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
          scale_fill_manual(values = c(
            "0-4"   = "#2166AC",
            "5-10"  = "#92C5DE",
            "11-15" = "#F4A582",
            "15+"   = "#B2182B"
          ),
          breaks = c("0-4", "5-10", "11-15", "15+")) +
          geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey30") +
          labs(x = paste0("P(", column_for_panel, ")"), y = "Density", fill = "Age",
               title = wrap_title(paste0("P(", column_for_panel, ") ", "in age-ranges in ", cohort_name))) +
          theme(axis.title = element_text ( hjust = 0.5,face = "bold"),
                plot.title = element_text(hjust = 0.5, face = "bold")) +
          coord_cartesian(clip = "off")
        
        ggsave(file.path(save_in, paste0(label_shorthand, "_age_density.png")),
               p_age_density, width = 9, height = 6, dpi = 300, bg = "white")
      }
      
      if ("MYC/MYCN" %in% plots_to_run) {
        
        prob_df$MYC_status <- as.factor(test_pheno[[MYC_status]])
        
        p_myc_density <- ggplot(prob_df,
                                aes(x = .data[[column_for_panel]], fill = MYC_status)) +
          geom_density(alpha = 0.6) +
          scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
          scale_fill_manual(values = c("AMP" = "orange", "Neutral" = "grey70")) +
          geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey30") +
          labs(x = paste0("P(", column_for_panel, ")"), y = "Density", fill = "MYC",
               title = wrap_title(paste0("P(", column_for_panel, ") ", "by MYC Amplification Status in ", cohort_name))) +
          theme(axis.title = element_text(hjust = 0.5, face = "bold"),
                plot.title      = element_text(hjust = 0.5, face = "bold")) +
          coord_cartesian(clip = "off")
        
        ggsave(file.path(save_in, paste0(label_shorthand, "_MYC_density.png")),
               p_myc_density, width = 9, height = 6, dpi = 300, bg = "white")
        
        prob_df$MYCN_status <- as.factor(test_pheno[[MYCN_status]])
        
        p_mycn_density <- ggplot(prob_df,
                                 aes(x = .data[[column_for_panel]], fill = MYCN_status)) +
          geom_density(alpha = 0.6) +
          scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
          scale_fill_manual(values = c("AMP" = "orange", "Neutral" = "grey70")) +
          geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey30") +
          labs(x = paste0("P(", column_for_panel, ")"), y = "Density", fill = "MYCN",
               title = wrap_title(paste0("P(", column_for_panel, ") ", "by MYCN Amplification Status in ", cohort_name))) +
          theme(axis.title = element_text (hjust = 0.5,face = "bold"),
                plot.title = element_text(hjust = 0.5, face = "bold")) +
          coord_cartesian(clip = "off")
        
        
        ggsave(file.path(save_in, paste0(label_shorthand, "_MYCN_density.png")),
               p_mycn_density, width = 9, height = 6, dpi = 300, bg = "white")
      }
      
      if ("Mol/Subgroup" %in% plots_to_run) {
        
        prob_df$mol_group <- as.factor(test_pheno[[mol_group]])
        
        p_mol_group <- ggplot(prob_df,
                              aes(x = mol_group, y = .data[[column_for_panel]], fill = mol_group)) +
          geom_boxplot(outlier.size = 0.5, linewidth = 0.4, alpha = 0.8) +
          geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
          scale_fill_manual(values = MOL_GRP_COLOURS) +
          stat_compare_means(method = "kruskal.test",
                             label.y = 1,
                             size    = 2.5) +
          geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
          scale_y_continuous(limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
          labs(x = "Molecular group", y = paste0("P(", column_for_panel, ")"),
               title = wrap_title(paste0("P(", column_for_panel, ") ", "by Principal Molecular Group in ", cohort_name))) +
          theme(legend.position = "none",
                axis.text.x     = element_text(angle = 45, hjust = 0.5, size = 6),
                axis.title = element_text ( hjust = 0.5,face = "bold"),
                plot.title = element_text(hjust = 0.5, face = "bold")) +
          coord_cartesian(clip = "off")
        
        ggsave(file.path(save_in, paste0(label_shorthand, "_mol_grp_boxplot.png")),
               p_mol_group, width = 9, height = 6, dpi = 300, bg = "white")
        
        prob_df$sub_group <- as.factor(test_pheno[[sub_group]])
        
        p_subgroup <- ggplot(prob_df,
                             aes(x = sub_group, y = .data[[column_for_panel]], fill = sub_group)) +
          geom_boxplot(outlier.size = 0.5, linewidth = 0.4, alpha = 0.8) +
          geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
          scale_fill_manual(values = SUB_GRP_COLOURS) +
          stat_compare_means(method = "kruskal.test",
                             label.y = 1,
                             size    = 2.5) +
          geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
          scale_y_continuous(limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
          labs(x = "Sub group", y = paste0("P(", column_for_panel, ")"),
               title = wrap_title(paste0("P(", column_for_panel, ") ", "by Molecular Sub group in ", cohort_name))) +
          theme(legend.position = "none",
                axis.text.x     = element_text(angle = 45, hjust = 0.5, size = 6),
                axis.title = element_text ( hjust = 0.5,face = "bold"),
                plot.title = element_text(hjust = 0.5, face = "bold")) +
          coord_cartesian(clip = "off")
        
        ggsave(file.path(save_in, paste0(label_shorthand, "_sub_grp_boxplot.png")),
               p_subgroup, width = 9, height = 6, dpi = 300, bg = "white")
      }
      
      if ("Survival" %in% plots_to_run) {
        
        p_km_class <- build_km_plot(km_fit_class, plot_title = paste0("Survival by Class in ", cohort_name), palette = class_palette, theme_to_apply = TRUE)
        ggsave(file.path(save_in, paste0(label_shorthand, "_KM_class.png")),
               p_km_class, width = 9, height = 6, dpi = 300, bg = "white")
        
        p_km_pred <- build_km_plot(km_fit_pred, plot_title = paste0("Survival by Predicted Label in ", cohort_name),palette = pred_palette, theme_to_apply = TRUE)
        ggsave(file.path(save_in, paste0(label_shorthand, "_KM_pred.png")),
               p_km_pred, width = 9, height = 6, dpi = 300, bg = "white")
        
        p_km_true <- build_km_plot(km_fit_true_label, plot_title = paste0("Survival by True Histology Label in ", cohort_name), palette = true_label_palette, theme_to_apply = TRUE)
        ggsave(file.path(save_in, paste0(label_shorthand, "_KM_true_label.png")),
               p_km_true, width = 9, height = 6, dpi = 300, bg = "white")
      }
      
    } else {
      stop("Probability score rownames don't match Sample IDs in the Phenotype table to plot clinical correlates")
    }
  }
}