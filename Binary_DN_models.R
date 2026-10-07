#=========================================================================================================================================#
# Training Binary Random Forest Models for Classification of Desmoplastic Nodular Histology in Medulloblastoma.

# This script includes loading of samplesheets and pre-processed methylation data (in-house),
# Initial training of Binary Random Forest models for classification of DN and non-DN MB groups across all principal Molecular groups
# Training Binary Random Forest models for DN histology prediction within SHH-pathway activated MB patient group
# Horvath (Horvath, 2013) Epigenetic Age prediction and Epigenetic Age acceleration
# Statistical testing for significant difference and correlation 
# Differential Gene Expression Analysis 
# Heatmaps, Bar plots, Scatter plots and density plots to visualise output probability score distribution and other correlates.
# Sample-size projection on SHH-MB pool : a log-linear learning curve of AUC against total n
# and a DN / non-DN count grid whose coefficients show which class's samples matter more,
# each evaluated by repeated stratified CV and by training on subsets of the pool and testing on Cavalli-SHH.
#=========================================================================================================================================#

# IMPORTANT: Change default labels and titles in the defined functions according to your application. 

# ------ Start from a clean session ------ #

rm(list = ls())
 
# ------ Load packages ------ #

require(dplyr) 
require(caret)
require(MLmetrics) 
require(pROC)
require(circlize) 
require(ComplexHeatmap)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(stringr)
require(ggrepel)
require(tibble)
require(tidyr)  

# ------------------------------------------------------------------------------------------------------------------- #
# Load forked methylClass - do NOT add require(methylClass)/library(methylClass)
# anywhere in this script
# ------------------------------------------------------------------------------------------------------------------- #

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::load_all("~/methylClass-fork")


# ------ Custom Functions ------ #

# for title-wrapping to avoid the plots getting cut off over long descriptive titles.

wrap_title <- function(x, width = 55) {
  stringr::str_wrap(x, width = width)
}


# Wrapper function to test and evaulate models trained with both caret and methylClass (Liu,2024)
# Evaluation metrics (overall and/or classwise)- Confusion Matrix, accuracy, balanced accuracy, sensitivity,specificity,precision, recall, MCC, 
# ROC-AUC, PR-AUC,F1 score and kappa scores
# Plots for Confusion Matrix, ROC and PR Curves. 

# IMPORTANT: Change Train set and Test set names, and labels/titles accordingly

test_and_evaluate_model <- function(model, method, test_data, true_test_labels, test_data_name, model_name,
                                    train_set_name,calibrated = TRUE, calibration_method = NULL,
                                    save_in = "~/Thesis_models/DN_models/") { 
# model : model object

# method: 'caret' or 'methylClass' ; "character"
  
# test_data : test beta value matrix
  
# true_test_labels : Given labels before prediction, "vector"
  
# model_name (To set labels of plots and names of the saved files): "Random Forest", "Support Vector Machine", "XGBoost"; "character"
  
# train_set_name (To set labels of plots and names of the saved files) : "Train NMB" ; "character"
  
# calibrated (To calibrated prediction probabilities) : TRUE
  
# calibration_method: "LR","FLR" and/or "MR" (applicable for method = "methylClass") ; "character"

 
  dir.create(save_in, recursive = TRUE, showWarnings = FALSE)
  
  # Shortform to save files
  
  model_shorthand<-ifelse(grepl("Random Forest", model_name),"RF",
                          ifelse(grepl("Support Vector Machine", model_name), "SVM", "XGB")) 
  
  #IMPORTANT: PLEASE CHANGE ACCORDINGLY
  
  train_set_shorthand <- dplyr::case_when(
    grepl("Under-sampled", train_set_name, ignore.case = TRUE) ~ "US_NMB",
    grepl("Train NMB", train_set_name, ignore.case = TRUE) & !grepl("SHH", train_set_name, ignore.case = TRUE) ~ "NMB",
    grepl("SHH-Inf-MB", train_set_name, ignore.case = TRUE) & grepl("extended", train_set_name, ignore.case = TRUE) ~ "SHH_Inf_sansage_ext",
    grepl("SHH-Inf-MB", train_set_name, ignore.case = TRUE) & grepl("Age-related", train_set_name, ignore.case = TRUE) ~ "SHH_Inf_sansage",
    grepl("SHH-Inf-MB", train_set_name, ignore.case = TRUE) ~ "SHH_Inf",
    grepl("SHH", train_set_name, ignore.case = TRUE) & grepl("extended", train_set_name, ignore.case = TRUE) & grepl("Age-related", train_set_name, ignore.case = TRUE) ~ "SHH_sansage_ext",
    grepl("SHH", train_set_name, ignore.case = TRUE) & grepl("Age-related", train_set_name, ignore.case = TRUE) ~ "SHH_sansage",
    grepl("SHH", train_set_name, ignore.case = TRUE) ~ "SHH",
    TRUE ~ make.names(train_set_name)
  )
  
  # To check required packages
  
  required_pgs <- c("caret", "pROC", "ggplot2", "tibble", "stringr", "MLmetrics", "dplyr")
  
  missing_pgs <- required_pgs[!sapply(required_pgs, requireNamespace, quietly = TRUE)]
  
  if (length(missing_pgs) > 0) {
    message("Installing missing packages for evaluation methods: ", paste(missing_pgs, collapse = ", "))
    install.packages(missing_pgs, dependencies = TRUE)
  }
  
  invisible(lapply(required_pgs, function(pkg) {
    library(pkg, character.only = TRUE)
  }))
  
  # Check data types
  
  if (!is.data.frame(test_data) && !is.matrix(test_data)) {
    test_data <- as.data.frame(test_data)
  }
  
  if (!is.factor(true_test_labels)) {
    true_test_labels <- as.factor(true_test_labels)
  }
  
  if (nrow(test_data) != length(true_test_labels)) {
    stop("Check the match between the rownames of the test data and of the ground truth labels")
  }
  
  classes <- levels(true_test_labels)
  
  if (method == "caret") {
    model_labels <- model$levels
    test_labels <- levels(true_test_labels)
  } else {
    model_labels <- levels(as.factor(model$mod[[1]]$classes))
    test_labels <- levels( as.factor(true_test_labels)) 
  }
  
  if (!setequal(model_labels, test_labels)) {
    stop("The histology levels of the model and the test data DO NOT match")
  }
  
  # ---Predictions---
  
  if (method == "caret") {
    
    pred_labels <- tryCatch(
      predict(model, newdata = test_data),
      error = function(e) {
        stop("predict() failed: ", conditionMessage(e))
      })
    
    pred_probs <- tryCatch(
      predict(model, newdata = test_data, type = "prob"),
      error = function(e) {
        stop("predict(prob) failed: ", conditionMessage(e))
      })
    
    pred_probs <- as.data.frame(pred_probs) 
    
  } else {
    if (calibrated){
      if(calibration_method == "MR"){
        pred <- tryCatch(
          methylClass::mainpredict(newdat = test_data, mod = model$mod, calibratemod = model$glmnet.calfit),
          error = function(e) stop("methylClass::predict error: ", conditionMessage(e))
        )
        pred_probs <- as.data.frame(pred$probs)
        pred_labels <- as.factor(pred$pres)
      } else if (calibration_method == "LR"){
        pred <- tryCatch(
          methylClass::mainpredict(newdat = test_data, mod = model$mod, calibratemod = model$platt.calfits),
          error = function(e) stop("methylClass::predict error: ", conditionMessage(e))
        )
        pred_probs <- as.data.frame(pred$probs)
        pred_labels <- as.factor(pred$pres)
      } else if ( calibration_method == "FLR"){
        pred <- tryCatch(
          methylClass::mainpredict(newdat = test_data, mod = model$mod, calibratemod = model$platt.brglm.calfits),
          error = function(e) stop("methylClass::predict error: ", conditionMessage(e))
        )
        pred_probs <- as.data.frame(pred$probs)
        pred_labels <- as.factor(pred$pres)
        
      } else {
        stop("calibration_method must be one of 'MR', 'LR' or 'FLR' when calibrated = TRUE")
      }
    } else {
      pred <- tryCatch(
        methylClass::mainpredict(newdat = test_data, mod = model$mod),
        error = function(e) stop("methylClass::predict error: ", conditionMessage(e))
      )
      pred_labels <- as.factor(pred$rawpres)
      pred_probs <- as.data.frame(pred$scores)
    }
  }
  
  # ---Confusion Matrix
  
  cm <- tryCatch(
    caret::confusionMatrix(pred_labels, true_test_labels, mode = "everything"),
    error = function(e) {
      stop("confusion matrix encountered error: ", conditionMessage(e))
    })
  
  cm_df <- as.data.frame(cm$table)
  
  levls <- levels(as.factor(true_test_labels))
  
  molgrp_pattern    <- "MB$"                                   
  subgroup_pattern  <- "^SHH[_-][0-9]+$"
  histology_pattern <- "^(DN|non_DN)$"
  
  # IMPORTANT: CHANGE LABELS ACCORDINGLY 
  
  if (any(grepl(molgrp_pattern, levls, ignore.case = TRUE))) {
    
    COLOURS <- c("WNT-MB" = "#4DBBD5", "SHH-MB" = "#E64B35",
                 "Group3-MB" = "yellow", "Group4-MB" = "green")
    
  } else if (any(grepl(subgroup_pattern, levls, ignore.case = TRUE))) {
    
    COLOURS <- c(
      "SHH_1" = "#FFB6C1", "SHH_2" = "#FF69B4", "SHH_3" = "#C71585", "SHH_4" = "#FF1493",
      "SHH_3A" = "#9B1B6E", "SHH_3B" = "#D4449A", "SHH_3C" = "#F0A8D0",
      "G34_I"   = "#800080", "G34_II"  = "#C71585", "G34_III" = "#FF8C00", "G34_IV" = "#FFD700",
      "G34_V"   = "#ADFF2F", "G34_VI"  = "#90EE90", "G34_VII" = "#ADD8E6", "G34_VIII" = "#006400",
      "WNT" = "darkblue", "MBNOS" = "black", "NA" = "grey"
    )
    
  } else if (any(grepl(histology_pattern, levls, ignore.case = TRUE))) {
    
    COLOURS <- c("DN" = "darkgreen", "non_DN" = "purple", "CLA" = "purple","LCA" = "pink","MBEN" = "green")
    
  } else {
    warning("No known colour mapping for levels:", paste(levls, collapse = ", "))
    return(NULL) 
  }
  
  
  cm_plot <- ggplot(cm_df, aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = Reference, alpha = Freq), color = "white") +
    geom_text(aes(label = Freq), size = 4,
              color = ifelse(cm_df$Freq > max(cm_df$Freq) * 0.5, "white", "black"), fontface = "bold") +
    scale_fill_manual(values = COLOURS, guide = "none") +
    scale_alpha_continuous(range = c(0.15, 1), guide = "none") +
    labs( title = wrap_title(
      paste0("Confusion Matrix for ", model_name, " trained on ", train_set_name, ", tested on ", test_data_name)
    ),
    x = "True Class", y = "Predicted Class"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          axis.title = element_text ( hjust = 0.5,face = "bold"),
          plot.margin = ggplot2::margin(t = 15, r = 15, b = 15, l = 15)) +
    coord_cartesian(clip = "off")
  
  ggsave(file.path(save_in, paste0(model_shorthand, "_",train_set_shorthand,"_", test_data_name, "_confusion_matrix.png")),
         cm_plot, width = 9, height = 6, dpi = 300, bg = "white")
  
  #--- Overall Accuracy
  
  accuracy <- unname(cm$overall["Accuracy"])
  
  
  #---Kappa
  
  kappa <- unname(cm$overall["Kappa"])
  
  # Brier and Log-loss
  
  P  <- as.matrix(pred_probs[, classes, drop = FALSE])
  storage.mode(P) <- "double"
  rs <- rowSums(P)
  P  <- P / ifelse(is.finite(rs) & rs > 0, rs, NA_real_)
  lab    <- as.character(true_labels)
  Y      <- outer(lab, classes, "==") * 1
  brier  <- mean(rowSums((P - Y)^2), na.rm = TRUE)
  p_true <- P[cbind(seq_len(nrow(P)), match(lab, classes))]
  logloss <- -mean(log(pmin(pmax(p_true, 1e-15), 1)), na.rm = TRUE)
  
  
  #--- Matthew's Correlation Coefficient 
  
  tab     <- as.matrix(cm$table)[classes, classes, drop = FALSE]
  s       <- sum(tab)
  correct <- sum(diag(tab))
  row_tot <- rowSums(tab)
  col_tot <- colSums(tab)
  sens    <- as.numeric(sensitivity[classes])
  spec    <- as.numeric(specificity[classes])
  
  mcc_den <- sqrt((s^2 - sum(col_tot^2)) * (s^2 - sum(row_tot^2)))
  mcc     <- if (is.finite(mcc_den) && mcc_den > 0) (correct * s - sum(row_tot * col_tot)) / mcc_den else 0
  
  #---Sensitivity and Specificity ---
  
  byclass <- cm$byClass
  
  if (is.null(dim(byclass))) {
    positive_class <- cm$positive
    other_class <- setdiff(classes, positive_class)
    
    sensitivity <- setNames(
      c(unname(byclass["Sensitivity"]), unname(byclass["Specificity"])),
      c(positive_class, other_class)
    )
    specificity <- setNames(
      c(unname(byclass["Specificity"]), unname(byclass["Sensitivity"])),
      c(positive_class, other_class)
    )
    
  } else {
    sensitivity <- byclass[, "Sensitivity"]
    specificity <- byclass[, "Specificity"]
    names(sensitivity) <- gsub("^Class: ", "", rownames(byclass))
    names(specificity) <- gsub("^Class: ", "", rownames(byclass))
  }
  
  #--- Balanced Accuracy 
  
  balanced_accuracy = mean(sensitivity, na.rm = TRUE)
  
  #---Precision, Recall and F1--- 
  
  stopifnot(identical(names(dimnames(cm$table)), c("Prediction", "Reference")))
  
  cm_table <- t(cm$table)   # rows = Reference (actual/truth), columns = Prediction
  
  classwise_precision <- sapply(classes, function(cls) {
    tp <- cm_table[cls, cls]
    fp <- sum(cm_table[, cls]) - tp  
    tp / (tp + fp)
  })
  
  classwise_recall <- sapply(classes, function(cls) {
    tp <- cm_table[cls, cls]
    fn <- sum(cm_table[cls, ]) - tp  
    tp / (tp + fn)
  })
  
  names(classwise_precision) <- classes
  names(classwise_recall)    <- classes
  
  classwise_f1 <- 2 * (classwise_precision * classwise_recall) / (classwise_precision + classwise_recall)
  names(classwise_f1) <- classes
  
  precision <- mean(classwise_precision, na.rm = TRUE)  # macro precision
  recall    <- mean(classwise_recall,    na.rm = TRUE)  # macro recall
  f1        <- mean(classwise_f1,        na.rm = TRUE)  # macro F1
  
  #---ROC and ROC-AUC
  
  # Class wise ROC and AUC
  
  roc_list <- tryCatch({
    setNames(lapply(classes, function(cls) {
      pROC::roc(response = as.numeric(true_test_labels == cls),
                predictor = pred_probs[[cls]], quiet = TRUE)
    }), classes)
  }, error = function(e) {
    stop("ROC Curve error: ", conditionMessage(e))
  })
  
  class_auc <- sapply(roc_list, function(r) as.numeric(pROC::auc(r)))
  
  # Hand and Till (2001) Multi-class AUC
  
  handtill_obj <- tryCatch(
    pROC::multiclass.roc(response = true_test_labels, predictor = pred_probs, quiet = TRUE),
    error = function(e) stop("Hand-Till AUC Error: ", conditionMessage(e))
  )
  
  macro_auc <- as.numeric(handtill_obj$auc)
  
  # --- DeLong 95% CI for the AUC (binary only)
  
  macro_auc_ci <- c(NA_real_, NA_real_)
  if (length(classes) == 2 && !is.null(roc_list)) {
    ci <- tryCatch(as.numeric(pROC::ci.auc(roc_list[[classes[1]]], method = "delong")), error = function(e) NULL)
    if (!is.null(ci) && length(ci) == 3) macro_auc_ci <- ci[c(1, 3)]
  }
  
  # --- ROC curve plot
  
  roc_df <- do.call(rbind, lapply(names(roc_list), function(cls) {
    r <- roc_list[[cls]]
    data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities,
               class = cls,
               class_label = paste0(cls, " (AUC = ", round(pROC::auc(r), 3), ")"))
  }))
  
  roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = class)) +
    geom_line(linewidth = 1) +
    geom_abline(linetype = "dashed", color = "grey50") +
    scale_color_manual(values = COLOURS,
                       labels = setNames(roc_df$class_label, roc_df$class)[unique(roc_df$class)]) +
    labs(title = wrap_title(paste("ROC Curves for", model_name, "trained on", train_set_name,
                                  "tested on", test_data_name)),
         subtitle = paste0("Macro AUC (Hand-Till) = ", round(macro_auc, 3)),
         x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate (Sensitivity)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5,face = "bold"),
          plot.margin = ggplot2::margin(t = 15, r = 15, b = 15, l = 15),
          axis.title = element_text ( hjust = 0.5,face = "bold")) +
    coord_cartesian(clip = "off")
  
  ggsave(file.path(save_in, paste0(model_shorthand, "_", train_set_shorthand, "_", test_data_name, "_roc_curves.png")),
         roc_plot, width = 9, height = 6, dpi = 300, bg = "white")
  
  #--- Precision Recall Curves and PR-AUC
  
  pr_list <- tryCatch({
    setNames(lapply(classes, function(cls) {
      prs <- pROC::coords(roc_list[[cls]], x = "all",
                          ret = c("recall", "precision"), transpose = FALSE)
      prs <- prs[!is.na(prs$precision) & !is.na(prs$recall), ]
      prs <- prs[order(prs$recall), ]
      prs
    }), classes)
  }, error = function(e) {
    stop("PR Curve Error: ", conditionMessage(e))
  })
  
  pr_auc <- tryCatch({
    sapply(classes, function(cls) {
      MLmetrics::PRAUC(y_pred = pred_probs[[cls]],
                       y_true = as.numeric(true_test_labels == cls))
    })
  }, error = function(e) {
    stop("PRAUC Error: ", conditionMessage(e))
  })
  
  names(pr_auc) <- classes
  
  pr_df <- do.call(rbind, lapply(classes, function(cls) {
    d <- pr_list[[cls]]
    data.frame(recall = d$recall, precision = d$precision,
               class = cls,
               class_label = paste0(cls, " (PR-AUC =", round(pr_auc[[cls]], 3), ")"))
  }))
  
  #--Macro PR-AUC
  
  macro_pr_auc <- mean(pr_auc, na.rm = TRUE)  
  
  pr_plot <- ggplot(pr_df, aes(x = recall, y = precision, color = class)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = COLOURS,
                       labels = setNames(pr_df$class_label,pr_df$class)[unique(pr_df$class)]) +
    labs(title = wrap_title(paste("Precision-Recall Curves for", model_name, "trained on",
                                  train_set_name, "tested on", test_data_name)),
         subtitle = paste0("Macro PR-AUC = ", round(macro_pr_auc, 3)),
         x = "Recall", y = "Precision", color = "Class") +
    theme_minimal() + ylim(0, 1) + xlim(0, 1) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text ( hjust = 0.5,face = "bold"),
          plot.margin = ggplot2::margin(t = 15, r = 15, b = 15, l = 15)) +
    coord_cartesian(clip = "off")
  
  ggsave(file.path(save_in, paste0(model_shorthand, "_",train_set_shorthand,"_", test_data_name, "_pr_curves.png")),
         pr_plot, width = 9, height = 6, dpi = 300, bg = "white")
  
  #---Summary
  
  pred_probs$pred <- as.factor(pred_labels)
  
  pred_probs %>% 
    tibble::rownames_to_column(var = "Sample_ID") %>%
    write.csv(file.path(save_in, paste0(model_shorthand, "_", train_set_shorthand, "_", 
                                        test_data_name, "_probabilities_.csv")), 
              row.names = FALSE)
  
  extra <- extra_eval_metrics(cm_table = cm_table, pred_probs = pred_probs, true_labels = true_test_labels,
                              classes = classes, roc_list = roc_list,
                              sensitivity = sensitivity, specificity = specificity)
  
  summary_df <- data.frame(
    model = model_shorthand,
    test_data = test_data_name,
    macro_precision = precision,
    mcc             = mcc,
    brier           = brier,
    log_loss        = logloss,
    macro_recall = recall,
    macro_f1 = f1,
    macro_auc = macro_auc,
    macro_auc_ci = macro_auc_ci,
    macro_pr_auc = macro_pr_auc,
    kappa = kappa
  )
  summary_df <- cbind(summary_df, extra$summary)
  
  #---Class-wise metrics
  
  classwise_df <- data.frame(
    model = model_shorthand,
    test_data = test_data_name,
    class = classes, 
    sensitivity = as.numeric(sensitivity[classes]),
    specificity = as.numeric(specificity[classes]),
    precision   = as.numeric(classwise_precision[classes]),
    recall      = as.numeric(classwise_recall[classes]),
    f1          = as.numeric(classwise_f1[classes]),
    auc         = as.numeric(class_auc[classes]),
    pr_auc      = as.numeric(pr_auc[classes]),
    youden_j    = sensitivity + specificity- 1,
    
  )
  classwise_df <- cbind(classwise_df, extra$classwise)
  
  #---Save as .csv files
  
  write.csv(summary_df, file.path(save_in, paste0(model_shorthand, "_", train_set_shorthand, "_", test_data_name, "_summary_metrics.csv")), row.names = FALSE)
  write.csv(classwise_df, file.path(save_in, paste0(model_shorthand, "_", train_set_shorthand, "_", test_data_name, "_classwise_metrics.csv")), row.names = FALSE)
  
  
  #---Return to object 
  
  list(
    confusion_matrix    = cm,
    cm_plot             = cm_plot,
    roc_list            = roc_list,
    roc_plot            = roc_plot,
    pr_list             = pr_list,
    pr_plot             = pr_plot,
    class_auc           = class_auc,
    pr_auc              = pr_auc,
    macro_pr_auc        = macro_pr_auc,
    sensitivity         = sensitivity,
    specificity         = specificity,
    classwise_precision = classwise_precision,
    classwise_recall    = classwise_recall,
    classwise_f1        = classwise_f1,
    kappa               = kappa,
    macro_auc           = macro_auc,   
    handtill_obj        = handtill_obj,
    summary             = summary_df,
    classwise           = classwise_df,
    pred_probs          = pred_probs
  )
}


# Wrapper function to produce probability score heatmaps with the output from test_and_evaluate() 
# and to also produce complementing barplots/scatter plots/density plots for
# probability score distribution/correlations in clinical correlate groups

# IMPORTANT: Change Train set and Test set names, and labels/titles accordingly

prob_heatmaps <- function(prob_df, pred, test_pheno, sample_id_column, heatmap_label,
                          true_test_labels, mol_group, sub_group, MYC_status, MYCN_status, OS, follow_up, Age,
                          continuum = "DN/bi", heatmap_panel = TRUE, correlate_plot = NULL,
                          column_for_panel = "DN", file_prefix = NULL, save_in = "~/Thesis_models/DN_models/") { 

  # prob_df: probability score matrix from test_and_evaluate()$pred_probs
  
  # pred: prediction labels from test_and_evaluate()$pred_probs$pred
  
  # test_pheno : Phenotype table for test data 
  
  # true_test_labels : Given labels of the test data before prediction from test_pheno;"character"
  
  # mol_group : Principal molecular group column from test_pheno;"character"
  
  # sub_group : molecular sub group column from test_pheno; "character"
  
  # MYC_status : MYC Amplification column from test_pheno ; "character"
  
  # MYCN_status :  MYCN Amplification column from test_pheno ; "character"
  
  # OS : Overall Survival (Event : 1, No Event : 0) column from test_pheno ; "character"
  
  # follow_up : test_pheno column for Follow up time for Overall Survival; "character"
  
  # Age : Chronological Age column from test_pheno ; "character"
  
  # continuum : To choose which level of true_test_label would be plotted for its probability score enrichment across all the test samples ; "character"
  
  # heatmap_panel : TRUE if other clinical correlate should be plotted against probability scores and displayed next to probability heatmap. 
  
  # correlate_plot : To choose with clincal correlate plot to be saved, along with saving standalone probability score heatmap: "character"
  
  # column_for_panel : level of true_test_label whose probability score is to be plotted in the correlate plots; "character"
  
  # file_prefix : file name prefix for the outputs to be saved with; "character"
  
  
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
  
  if (!is.null(file_prefix)) {
    label_shorthand <- file_prefix
  } else if (grepl("Random", heatmap_label) && grepl("Test NMB", heatmap_label)) {
    label_shorthand <- "RF_test_NMB"
  } else if (grepl("Random", heatmap_label) && grepl("Cavalli-SHH-Inf-MB", heatmap_label) && grepl("Age", heatmap_label, ignore.case = TRUE)) {
    label_shorthand <- "RF_Cavalli_SHH_Inf_sansage"
  } else if (grepl("Random", heatmap_label) && grepl("Cavalli-SHH-Inf-MB", heatmap_label)) {
    label_shorthand <- "RF_Cavalli_SHH_Inf"
  } else if (grepl("Random", heatmap_label) && grepl("Cavalli-SHH-MB", heatmap_label) && grepl("Age", heatmap_label, ignore.case = TRUE)) {
    label_shorthand <- "RF_Cavalli_SHH_sansage"
  } else if (grepl("Random", heatmap_label) && grepl("Cavalli-SHH", heatmap_label)) {
    label_shorthand <- "RF_Cavalli_SHH"
  } else if (grepl("Random", heatmap_label) && grepl("Cavalli", heatmap_label)) {
    label_shorthand <- "RF_Cavalli"
  } else if (grepl("Support", heatmap_label) && grepl("Test NMB", heatmap_label)) {
    label_shorthand <- "SVM_test_NMB"
  } else if (grepl("Support", heatmap_label) && grepl("Cavalli", heatmap_label)) {
    label_shorthand <- "SVM_Cavalli"
  } else if (grepl("Random", heatmap_label) && grepl("'NMB+Infant'", heatmap_label, fixed = TRUE)) {
    label_shorthand <- "RF_train" 
  } else if (grepl("Age", heatmap_label, ignore.case = TRUE) && grepl("'NMB+Infant'", heatmap_label, fixed = TRUE)) {
    label_shorthand <- "RF_train_sansage" 
  } else {
    stop("Check labels for filenames") 
  }  
  
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
    }
  } else {
    stop("Probability score rownames don't match Sample IDs in the Phenotype table")
  }
  
  if (all(rownames(prob_df) == test_pheno[[sample_id_column]])) {
    prob_df$True_label <- as.factor(test_pheno[[true_test_labels]])
    
    prob_df <- prob_df %>% mutate(bi_histo = case_when( 
      True_label %in% c("DN", "MBEN") ~ "DN",
      True_label %in% c("CLA", "LCA") ~ "non_DN"
    ))
    prob_df$bi_histo <- as.factor(prob_df$bi_histo)
    
    prob_df <- prob_df %>% mutate(Class = case_when(
      True_label == "CLA"  & pred == "non_DN" ~ "Correct_non_DN",
      True_label == "LCA"  & pred == "non_DN" ~ "Correct_non_DN",
      True_label == "DN"   & pred == "DN"     ~ "Correct_DN",
      True_label == "MBEN" & pred == "DN"     ~ "Correct_DN",
      True_label == "DN"   & pred == "non_DN" ~ "DN->non_DN",
      True_label == "MBEN" & pred == "non_DN" ~ "DN->non_DN",
      True_label == "CLA"  & pred == "DN"     ~ "non_DN->DN",
      True_label == "LCA"  & pred == "DN"     ~ "non_DN->DN",
      TRUE ~ paste0(True_label, "->", pred) 
    ))
  } else {
    stop("Probability score rownames don't match Sample IDs in the Phenotype table")
  } 
  
  prob_matrix <- as.matrix(prob_df[, c("DN", "non_DN")])
  
  stopifnot(
    sample_id_column %in% colnames(test_pheno),
    nrow(prob_matrix) == nrow(test_pheno),
    all(rownames(prob_matrix) == test_pheno[[sample_id_column]])
  )
  
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
  
  prob_df$True_label  <- factor(as.character(prob_df$True_label), 
                                levels = intersect(names(HISTO_COLOURS), 
                                                   unique(as.character(prob_df$True_label))))
  
  prob_df$pred        <- factor(as.character(prob_df$pred), 
                                levels = intersect(names(PRED_COLOURS), 
                                                   unique(as.character(prob_df$pred))))
  prob_df$Class      <- factor(prob_df$Class, 
                               levels = intersect(names(CLASS_COLOURS), 
                                                  unique(prob_df$Class)))
  
  class_display <- gsub("non_DN", "nonDN", as.character(prob_df$Class))
  
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
  
  ht <- Heatmap(
    t(prob_matrix),                     
    
    name      = "Probability",
    col       = col_fun,
    
    column_order = order(continuum_score),
    column_split = class_display,
    column_title_rot = 0,
    column_title_gp  = gpar(fontsize = 6, fontface = "bold"),
    column_gap       = unit(3, "mm"),
    
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    show_column_dend  = FALSE,
    
    top_annotation = col_ha,
    
    height = unit(3, "cm"),            
    
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
    
  } 
  
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
        scale_fill_manual(values = c("AMP" = "orange", "Neutral" = "grey70")) +
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
      
      pdf(file.path(save_in, paste0(label_shorthand, "_.pdf")),
          width  = 15,
          height = 10,
          bg     = "white")
      
      print(full_figure)
      
      dev.off()
      
    } else {
      stop("Probability score rownames don't match Sample IDs in the Phenotype table to plot clinical correlates")
    }
  } else {
    if (all(rownames(prob_df) == test_pheno[[sample_id_column]])) {
      
      available_plots <- c("Heatmap", "Age", "MYC/MYCN", "Mol/Subgroup", "Survival")
      
      library(dplyr)
      
      cohort_name <- case_when(
        grepl("Test NMB", heatmap_label)  ~ "Test NMB",
        grepl("'NMB' and 'Infant'", heatmap_label) ~ "Train SHH-MB",
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
        correlate_plot
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

# =========================================================================================================================================#
# ---- Balancing ablation helper functions ----
#
# Used throughout this script to compare ways RF training within maintrain can handle class imbalance:
#   - "balanced_both"                : bal_feature_sel = TRUE,  bal_training = TRUE   
#   - "balanced_featsel_imbal_train"  : bal_feature_sel = TRUE,  bal_training = FALSE  
#   - "neither"                      : bal_feature_sel = FALSE, bal_training = FALSE
# imp_replace = TRUE throughout (to sample with replacement when under sampling the majority class, to leave OOB samples in the minority class).
#
# =========================================================================================================================================#

ablation_configs <- list(
  balanced_both = list(
    feature_sel = TRUE, bal_feature_sel = TRUE, imp_replace = TRUE, bal_training = TRUE
  ),
  balanced_featsel_imbal_train = list(
    feature_sel = TRUE, bal_feature_sel = TRUE, imp_replace = TRUE, bal_training = FALSE
  ),
  neither = list(
    feature_sel = TRUE, bal_feature_sel = FALSE, imp_replace = TRUE, bal_training = FALSE
  ),
  
  imbal_featsel_bal_train = list(
    feature_sel = TRUE, bal_feature_sel = FALSE, imp_replace = TRUE, bal_training = TRUE
  ),
  no_feature_selection = list(
    feature_sel = FALSE, bal_feature_sel = FALSE, imp_replace = TRUE, bal_training = TRUE
  ),
  no_feature_selection_imbal_train = list(
    feature_sel = FALSE, bal_feature_sel = FALSE, imp_replace = TRUE, bal_training = FALSE
  )
)


# Trains one model per configuration, then evaluates each trained model against every test set
# supplied in `test_sets` 

# test_sets: a named list, each element itself a list with elements test_data, test_data_name, test_y
#            (and optionally heatmap_args for a full prob_heatmaps() output per config x test set)
#

run_rf_ablation <- function(train_betas, train_y, test_sets,
                            train_set_name, save_in,
                            seed = 42, ntrees = 100, p = 200,
                            configs = ablation_configs) {
  
  dir.create(save_in, recursive = TRUE, showWarnings = FALSE)
  
  models       <- list()
  evals        <- list()
  summary_rows <- list()
  
  for (cfg_name in names(configs)) {
    
    cfg <- configs[[cfg_name]]
    cfg_save_in <- file.path(save_in, cfg_name)
    dir.create(cfg_save_in, recursive = TRUE, showWarnings = FALSE)
    
    model_i <- tryCatch(
      maintrain(
        betas..              = train_betas,
        y..                  = train_y,
        method               = "RF",
        ntrees               = ntrees,
        p                    = p,
        topfeaturenumber     = NULL,
        subset.CpGs          = NULL,
        seed                 = seed,
        feature_sel          = cfg$feature_sel,
        bal_feature_sel      = cfg$bal_feature_sel,
        imp_replace          = cfg$imp_replace,
        bal_training         = cfg$bal_training,
        calibrationmethod    = c("MR", "FLR", "LR")
      ),
      error = function(e) {
        message(sprintf("FAILED training [%s / %s]: %s", train_set_name, cfg_name, conditionMessage(e)))
        NULL
      }
    )
    
    if (is.null(model_i)) next
    
    saveRDS(model_i, file.path(cfg_save_in, paste0("model_", cfg_name, ".rds")))
    models[[cfg_name]] <- model_i
    
    for (ts_name in names(test_sets)) {
      
      ts <- test_sets[[ts_name]]
      
      eval_i <- tryCatch(
        test_and_evaluate_model(
          model               = model_i,
          method              = "methylClass",
          calibrated          = TRUE,
          calibration_method  = "LR",
          train_set_name      = train_set_name,
          test_data           = ts$test_data,
          test_data_name      = ts$test_data_name,
          true_test_labels    = ts$test_y,
          model_name          = "Binary Random Forest",
          save_in             = cfg_save_in
        ),
        error = function(e) {
          message(sprintf("FAILED evaluation [%s / %s / %s]: %s", train_set_name, cfg_name, ts_name, conditionMessage(e)))
          NULL
        }
      )
      
      if (is.null(eval_i)) next
      
      evals[[paste0(cfg_name, "__", ts_name)]] <- eval_i
      
      if (!is.null(ts$heatmap_args)) {
        tryCatch(
          do.call(prob_heatmaps, c(
            list(
              prob_df       = eval_i$pred_probs,
              pred          = eval_i$pred_probs$pred,
              heatmap_label = paste0(train_set_name, "tested on ", ts$test_data_name),
              heatmap_panel = FALSE,
              file_prefix   = paste0(cfg_name, "_", ts_name),
              save_in       = cfg_save_in
            ),
            ts$heatmap_args
          )),
          error = function(e) {
            message(sprintf("FAILED prob_heatmaps [%s / %s / %s]: %s",
                            train_set_name, cfg_name, ts_name, conditionMessage(e)))
          }
        )
      }
      
    }
  }
  
  summary_df <- dplyr::bind_rows(summary_rows)
  
  write.csv(summary_df, file.path(save_in, "ablation_summary.csv"), row.names = FALSE)
  
  list(models = models, evals = evals, summary = summary_df)
}

# Comparison bar chart: one panel per test set, one bar per configuration. Default metric is macro AUC, with the
# DeLong 95% CI over the test samples as error bars and a dashed line at chance (0.5).

plot_ablation_comparison <- function(summary_df, title, save_path, metric = "macro_auc", drop_sets = "TrainSHH") {
  
  if (nrow(summary_df) == 0) {
    message("plot_ablation_comparison(): summary_df is empty, skipping plot for '", title, "'")
    return(invisible(NULL))
  }
  if (!(metric %in% colnames(summary_df))) {
    message("plot_ablation_comparison(): column '", metric, "' not found - use rebuild_ablation_summary() first. Skipping '", title, "'")
    return(invisible(NULL))
  }
  
  dropped <- intersect(drop_sets, unique(summary_df$test_set))
  plot_df <- summary_df[!(summary_df$test_set %in% drop_sets), , drop = FALSE]
  if (nrow(plot_df) == 0) { plot_df <- summary_df; dropped <- character(0) }
  
  plot_df$value  <- plot_df[[metric]]
  plot_df$config <- factor(plot_df$config, levels = unique(summary_df$config))
  has_ci  <- metric == "macro_auc" && all(c("auc_ci_lower", "auc_ci_upper") %in% colnames(plot_df))
  bounded <- metric %in% c("macro_auc", "macro_pr_auc", "accuracy", "balanced_accuracy")
  
  p <- ggplot(plot_df, aes(x = config, y = value)) +
    geom_col(fill = "#4C78A8", width = 0.7)
  
  if (metric == "macro_auc") p <- p + geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40")
  if (has_ci) p <- p + geom_errorbar(aes(ymin = auc_ci_lower, ymax = auc_ci_upper), width = 0.25, na.rm = TRUE)
  
  p <- p +
    geom_text(aes(label = ifelse(is.na(value), "NA", sprintf("%.3f", value))), y = 0.04, colour = "white", size = 3.2) +
    facet_wrap(~test_set) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(gsub("_", " ", x), 14)) +
    labs(x = "Feature-selection / training configuration", y = gsub("_", " ", metric),
         title = wrap_title(title),
         subtitle = paste(c(if (has_ci) "Error bars: DeLong 95% CI over test samples.",
                            if (metric == "macro_auc") "Dashed line: chance (0.5).",
                            if (length(dropped) > 0) paste0("Omitted (resubstitution): ", paste(dropped, collapse = ", "), ".")),
                          collapse = " ")) +
    theme_minimal() +
    theme(axis.text.x   = element_text(size = 8),
          plot.title    = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          strip.text    = element_text(face = "bold"))
  
  if (bounded) p <- p + scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.25, 0.5, 0.75, 1))
  
  ggsave(save_path, p, width = 10, height = 6, dpi = 300, bg = "white")
  p
}

# Pick the config with the best value of `metric` on `test_set`, from an ablation result.

# NOTE: choosing on Cavalli and then reporting Cavalli performance for the chosen config is optimistic for that config.

choose_best_config <- function(ablation_result, metric = "macro_auc", test_set = "Cavalli",
                               candidates = names(ablation_configs_extended), alpha = 0.05) {
  
  s <- ablation_result$summary
  if (is.null(s) || !(metric %in% colnames(s)) || !("auc_ci_lower" %in% colnames(s))) {
    message("choose_best_config(): summary lacks '", metric, "' / CI columns - rebuilding it from the stored evaluations.")
    s <- rebuild_ablation_summary(ablation_result)
  }
  
  s <- s[s$test_set == test_set & s$config %in% candidates, , drop = FALSE]
  if (nrow(s) == 0) stop("choose_best_config(): no rows for test set '", test_set, "' among the candidate configs.")
  
  s    <- s[order(-s[[metric]]), , drop = FALSE]
  best <- s$config[1]
  
  get_roc <- function(cfg) {
    ev <- ablation_result$evals[[paste0(cfg, "__", test_set)]]
    if (is.null(ev)) NULL else ev$roc_list[[1]]
  }
  roc_best <- get_roc(best)
  
  s$p_vs_best <- vapply(s$config, function(cfg) {
    if (metric != "macro_auc" || cfg == best || is.null(roc_best)) return(NA_real_)
    r <- get_roc(cfg)
    if (is.null(r)) return(NA_real_)
    tryCatch(as.numeric(pROC::roc.test(roc_best, r, method = "delong", paired = TRUE)$p.value),
             error = function(e) NA_real_)
  }, numeric(1))
  
  show_cols <- intersect(c("config", metric, "auc_ci_lower", "auc_ci_upper", "p_vs_best",
                           "DN_sensitivity", "nonDN_sensitivity", "DN_pred_rate", "macro_pr_auc"), colnames(s))
  message(sprintf("\nConfig ranking by %s on '%s' (%d candidate configs):", metric, test_set, nrow(s)))
  print(s[, show_cols, drop = FALSE], row.names = FALSE, digits = 3)
  
  message(sprintf("\nSelected: '%s'  (%s = %.3f)", best, metric, s[[metric]][1]))
  
  tied <- s$config[!is.na(s$p_vs_best) & s$p_vs_best >= alpha]
  if (length(tied) > 0) {
    message(sprintf("NOT clearly separated from: %s (paired DeLong p >= %.2f). The ranking among these is within test-set noise.",
                    paste(tied, collapse = ", "), alpha))
  } else if (metric == "macro_auc" && nrow(s) > 1) {
    message("The best config is separated from every other candidate (paired DeLong p < ", alpha, ").")
  }
  message(sprintf("NOTE: chosen on '%s'. Reporting performance for this config on the same set is optimistic.", test_set))
  
  list(best = best, table = s)
}

# Training binary Random Forest (RF) model on Train NMB for classification 
# of DN and non-DN histology

# ------ Load data ------ #

#Train and Test split of the 'NMB' Cohort is done with caret::createDataPartition 
#in the previous script "Distribution_charts.R". "NOS" histology labels have been removed before the split

train_NMB_pheno<-readRDS("~/NMB/NMB_train_pheno.rds")

if (any(is.na(train_NMB_pheno$Sample_ID))){
  train_NMB_pheno<-train_NMB_pheno[!is.na(train_NMB_pheno$Sample_ID),]
}

test_NMB_pheno<-readRDS("~/NMB/NMB_test_pheno.rds")

if (any(is.na(test_NMB_pheno$Sample_ID))){
  test_NMB_pheno<-test_NMB_pheno[!is.na(test_NMB_pheno$Sample_ID),]
} 

# Beta value matrix post pre processing with minfi- 450k + EPIC combined

NMB_betas<-readRDS("~/NMB/Final_NMB_betas.rds")

train_NMB<-NMB_betas[match(train_NMB_pheno$Sample_ID, rownames(NMB_betas)),]

if(is.null(train_NMB$bi_histo)){
  stopifnot(all(rownames(train_NMB) == train_NMB_pheno$Sample_ID))
  train_NMB$bi_histo<-as.factor(train_NMB_pheno$CPR_Histology)
}

train_NMB<-train_NMB %>% mutate(bi_histo = case_when( 
  bi_histo %in% c("CLA","LCA") ~ "non_DN", 
  bi_histo %in% c("DN", "MBEN") ~ "DN"))


train_NMB_pheno<-train_NMB_pheno %>% mutate(bi_histo = case_when( 
  CPR_Histology %in% c("CLA","LCA") ~ "non_DN", 
  CPR_Histology %in% c("DN", "MBEN") ~ "DN"))

# Test sets 

# Test NMB

test_NMB<-NMB_betas[match(test_NMB_pheno$Sample_ID, rownames(NMB_betas)),]

if(all(test_NMB_pheno$Sample_ID == rownames(test_NMB))){
  test_NMB_pheno<-test_NMB_pheno %>% mutate(bi_histo = case_when( 
    CPR_Histology  %in% c("CLA","LCA") ~ "non_DN", 
    CPR_Histology %in% c("DN", "MBEN") ~ "DN"))
}

# Independent Test dataset- Cavalli (Cavalli et al., 2017)

Cavalli_pheno<- readRDS("~/Cavalli/Final_Cavalli_pheno.rds")

# Matching Cavalli's Age distribution to train_NMB's

Cavalli_pheno<-Cavalli_pheno %>% filter(Age <= max(train_NMB_pheno$Age))

if (any(is.na(Cavalli_pheno$Sample_ID))){
  Cavalli_pheno<-Cavalli_pheno[!is.na(Cavalli_pheno$Sample_ID),]
}

if (any(is.na(Cavalli_pheno$histology))){
  Cavalli_pheno<-Cavalli_pheno[!is.na(Cavalli_pheno$histology),]
}


Cavalli_betas<-readRDS("~/Cavalli/Cavalli-test/Cavalli_test_betas.rds")

Cavalli_betas<-as.data.frame(t(Cavalli_betas)) 

Cavalli_betas<-Cavalli_betas[match(Cavalli_pheno$Sample_ID, rownames(Cavalli_betas)),] 

if(all(Cavalli_pheno$Sample_ID == rownames(Cavalli_betas))){
  Cavalli_pheno<-Cavalli_pheno %>% mutate(bi_histo = case_when( 
    histology  %in% c("CLA","LCA") ~ "non_DN", 
    histology  %in% c("DN", "MBEN") ~ "DN"))
} else {
  message("Mismatch in Sample IDs in Pheno and betas")
  
}

# ------ Balancing Train NMB with informed under-sampling (To match histology type distribution among MB molecular groups in the parent Train set) ------- #

counts    <- table(train_NMB_pheno$CPR_Histology, train_NMB_pheno$Mol_group)
DN_counts <- sum(train_NMB_pheno$bi_histo == "DN")

non_DN_total <- sum(counts[c("CLA", "LCA"), ])

train_NMB_pheno_US <- bind_rows(
  train_NMB_pheno %>% filter(CPR_Histology == "CLA", Mol_group == "Group4-MB") %>%
    slice_sample(n = round(DN_counts * counts["CLA", "Group4-MB"] / non_DN_total)),
  
  train_NMB_pheno %>% filter(CPR_Histology == "LCA", Mol_group == "Group4-MB") %>%
    slice_sample(n = round(DN_counts * counts["LCA", "Group4-MB"] / non_DN_total)),
  
  train_NMB_pheno %>% filter(CPR_Histology == "CLA", Mol_group == "Group3-MB") %>%
    slice_sample(n = round(DN_counts * counts["CLA", "Group3-MB"] / non_DN_total)),
  
  train_NMB_pheno %>% filter(CPR_Histology == "LCA", Mol_group == "Group3-MB") %>%
    slice_sample(n = round(DN_counts * counts["LCA", "Group3-MB"] / non_DN_total)),
  
  train_NMB_pheno %>% filter(bi_histo == "non_DN", Mol_group == "SHH-MB") %>%
    slice_sample(n = round(DN_counts * sum(counts[c("CLA", "LCA"), "SHH-MB"]) / non_DN_total)),
  
  train_NMB_pheno %>% filter(bi_histo == "non_DN", Mol_group == "WNT-MB") %>%
    slice_sample(n = round(DN_counts * sum(counts[c("CLA", "LCA"), "WNT-MB"]) / non_DN_total)),
  
  train_NMB_pheno %>% filter(bi_histo == "DN")
)

US_train_NMB<-train_NMB[match(train_NMB_pheno_US$Sample_ID,rownames(train_NMB)),]

if(!all(rownames(US_train_NMB) == train_NMB_pheno_US$Sample_ID)){
  message("CHECK the match between Undersampled Phenotype table and Beta value data frame")
}

# ------ Training methylClass (Liu, 2024) models with all 441870 input probes- 
# methylClass for high-dimensional methylation datasets is computationally inexpensive in comparison to caret ------ #


RF_bi_NMB <- maintrain(betas.. = train_NMB[,!colnames(train_NMB) %in% "bi_histo"],
                       y.. = as.factor(train_NMB$bi_histo),
                       method = "RF",
                       ntrees = 100, #default
                       p = 200, #default
                       topfeaturenumber = NULL,
                       subset.CpGs = NULL,
                       calibrationmethod = c("MR","FLR","LR")
)

# Train an RF model on under-sampled dataset

# NOTE: bal_feature_sel/bal_training explicitly FALSE here - US_train_NMB is already balanced by
# construction (informed under-sampling, stratified by molecular group, above).

RF_bi_US_NMB <- maintrain(betas.. = as.matrix(US_train_NMB[,!colnames(US_train_NMB) %in% "bi_histo"]),
                          y.. = as.factor(US_train_NMB$bi_histo),
                          method = "RF",
                          ntrees = 100, #default
                          p = 200, #default
                          topfeaturenumber = NULL,
                          subset.CpGs = NULL,
                          feature_sel     = TRUE,
                          bal_feature_sel = FALSE,
                          imp_replace     = TRUE,
                          bal_training    = FALSE,
                          calibrationmethod = c("MR","FLR","LR")
)


# -------- Testing of the models and Evaluation metrics -------------- #

# Binary model trained on Imbalanced data

message("Testing and evaluating on Test NMB")

if (all(rownames(test_NMB) == test_NMB_pheno$Sample_ID)) {
  test_NMB_results_imb <- tryCatch(
    test_and_evaluate_model(
      model = RF_bi_NMB,
      method = "methylClass",
      calibrated = TRUE,
      calibration_method = "LR",
      train_set_name = "Train NMB",
      test_data = test_NMB,
      test_data_name = "Test NMB",
      true_test_labels = as.factor(test_NMB_pheno$bi_histo),
      model_name = "Binary Random Forest",
      save_in = "~/Thesis_models/DN_models/Imb_Train_NMB/"
    ),
    error = function(e) {
      message(sprintf("FAILED [%s / %s]: %s",
                      "Binary Random Forest trained on Train NMB", "Test NMB", conditionMessage(e)))
      NULL
    }
  )
} else {
  stop("Mismatch in test set and test pheno Sample IDs")
}


test_NMB_imb_heatmap <- prob_heatmaps(
  prob_df           = test_NMB_results_imb$pred_probs,
  pred              = test_NMB_results_imb$pred_probs$pred,
  test_pheno        = test_NMB_pheno,
  sample_id_column  = "Sample_ID",
  true_test_labels  = "CPR_Histology",
  mol_group         = "Mol_grp",
  sub_group         = "sub_group",
  MYC_status        = "MYC Status",
  MYCN_status       = "MYCN Status",
  OS                = "OS",
  follow_up         = "Follow up years",  
  Age               = "Age",
  continuum         = "DN/bi",
  heatmap_panel     = FALSE,
  column_for_panel  = "DN",
  heatmap_label     = "Binary Random Forest Model for classification of DN histopathology,\n trained on Imbalanced Train NMB and tested on Test NMB",
  file_prefix       = "RF_imbNMB_vs_testNMB",
  save_in           = "~/Thesis_models/DN_models/Imb_Train_NMB/"
)

message("Testing and evaluating on Cavalli")

if(all(rownames(Cavalli_betas) == Cavalli_pheno$Sample_ID)){
  tryCatch(Cavalli_results_imb<-test_and_evaluate_model(
    model = RF_bi_NMB,
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "Train NMB",
    test_data = Cavalli_betas,
    test_data_name = "Cavalli",
    true_test_labels = as.factor(Cavalli_pheno$bi_histo),
    model_name = "Binary Random Forest",
    save_in = "~/Thesis_models/DN_models/Imb_Train_NMB/"
  ), error = function(e){
    message(sprintf("FAILED [%s / %s]: %s",
                    "Binary Random Forest trained on Train NMB", "Cavalli", conditionMessage(e)))
    NULL
  } )
  
} else {
  stop("Mismatch in Cavalli test set and test pheno Sample IDs")
}

Cavalli_imb_heatmap <- prob_heatmaps(
  prob_df           = Cavalli_results_imb$pred_probs,
  pred              = Cavalli_results_imb$pred_probs$pred,
  test_pheno        = Cavalli_pheno,
  sample_id_column  = "Sample_ID",
  true_test_labels  = "histology",
  mol_group         = "mol_grp",
  sub_group         = "sub_grp",
  MYC_status        = "MYC_Status",
  MYCN_status       = "MYCN_Status",
  OS                = "Dead",
  follow_up         = "OS (years)",
  Age               = "Age",
  continuum         = "DN/bi",
  heatmap_label     = "Binary Random Forest Model for classification of DN histopathology,\n trained on Imbalanced Train NMB and tested on Cavalli",
  heatmap_panel     = FALSE,
  file_prefix       = "RF_imbNMB_vs_Cavalli",
  save_in           = "~/Thesis_models/DN_models/Imb_Train_NMB/"
)

# ---- Balancing ablation: Train NMB (imbalanced), tested on both Test NMB and Cavalli ---- #

heatmap_args_testnmb <- list(
  test_pheno = test_NMB_pheno, sample_id_column = "Sample_ID", true_test_labels = "CPR_Histology",
  mol_group = "Mol_grp", sub_group = "sub_group", MYC_status = "MYC Status", MYCN_status = "MYCN Status",
  OS = "OS", follow_up = "Follow up years", Age = "Age", continuum = "DN/bi"
)

heatmap_args_cavalli_nmb <- list(
  test_pheno = Cavalli_pheno, sample_id_column = "Sample_ID", true_test_labels = "histology",
  mol_group = "mol_grp", sub_group = "sub_grp", MYC_status = "MYC_Status", MYCN_status = "MYCN_Status",
  OS = "Dead", follow_up = "OS (years)", Age = "Age", continuum = "DN/bi"
)

NMB_imb_ablation <- run_rf_ablation(
  train_betas    = train_NMB[, !colnames(train_NMB) %in% "bi_histo"],
  train_y        = as.factor(train_NMB$bi_histo),
  test_sets      = list(
    TestNMB = list(test_data = test_NMB, test_data_name = "Test NMB",
                   test_y = as.factor(test_NMB_pheno$bi_histo), heatmap_args = heatmap_args_testnmb),
    Cavalli = list(test_data = Cavalli_betas, test_data_name = "Cavalli",
                   test_y = as.factor(Cavalli_pheno$bi_histo), heatmap_args = heatmap_args_cavalli_nmb)
  ),
  train_set_name = "Train NMB (imbalanced)",
  save_in        = "~/Thesis_models/DN_models/Imb_Train_NMB/Ablation/",
  seed = 42, ntrees = 100, p = 200,
  configs        = ablation_configs_extended
)

print(NMB_imb_ablation$summary)

plot_ablation_comparison(
  NMB_imb_ablation$summary,
  title     = "Macro AUC Across RF Balancing Configurations - Train NMB (imbalanced)",
  save_path = "~/Thesis_models/DN_models/Imb_Train_NMB/Ablation/ablation_comparison.png"
)

# Binary model trained on balanced (under sampled) data

message("Testing and evaluating on Test NMB")

if (all(rownames(test_NMB) == test_NMB_pheno$Sample_ID)) {
  test_NMB_results_bal <- tryCatch(
    test_and_evaluate_model(
      model = RF_bi_US_NMB,
      method = "methylClass",
      calibrated = TRUE,
      calibration_method = "LR",
      train_set_name = "Under-sampled Train NMB",
      test_data = test_NMB,
      test_data_name = "Test NMB",
      true_test_labels = as.factor(test_NMB_pheno$bi_histo),
      model_name = "Binary Random Forest",
      save_in = "~/Thesis_models/DN_models/Balanced_train_NMB/"
    ),
    error = function(e) {
      message(sprintf("FAILED [%s / %s]: %s",
                      "Binary Random Forest trained on Train NMB", "Test NMB", conditionMessage(e)))
      NULL
    }
  )
} else {
  stop("Mismatch in test set and test pheno Sample IDs")
}


test_NMB_bal_heatmap<-prob_heatmaps(
  prob_df           = test_NMB_results_bal$pred_probs,
  pred              = test_NMB_results_bal$pred_probs$pred,
  test_pheno        = test_NMB_pheno,
  sample_id_column  = "Sample_ID",
  true_test_labels  = "CPR_Histology",
  mol_group         = "Mol_grp",
  sub_group         = "sub_group",
  MYC_status        = "MYC Status",
  MYCN_status       = "MYCN Status",
  OS                = "OS",
  follow_up         = "Follow up years",  
  Age               = "Age",
  continuum         = "DN/bi",
  heatmap_panel     = FALSE,
  column_for_panel  = "DN",
  heatmap_label     = "Binary Random Forest Model for classification of DN histopathology,\n trained on Under-sampled Train NMB and tested on Test NMB",
  file_prefix       = "RF_balNMB_vs_testNMB",
  save_in           = "~/Thesis_models/DN_models/Balanced_train_NMB/"
)


message("Testing and evaluating on Cavalli")

if(all(rownames(Cavalli_betas) == Cavalli_pheno$Sample_ID)){
  tryCatch(Cavalli_results_bal<-test_and_evaluate_model(
    model = RF_bi_US_NMB,
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "Under-sampled Train NMB",
    test_data = Cavalli_betas,
    test_data_name = "Cavalli",
    true_test_labels = as.factor(Cavalli_pheno$bi_histo),
    model_name = "Binary Random Forest",
    save_in = "~/Thesis_models/DN_models/Balanced_train_NMB/"
  ), error = function(e){
    message(sprintf("FAILED [%s / %s]: %s",
                    "Binary Random Forest trained on Train NMB", "Cavalli", conditionMessage(e)))
    NULL
  } )
  
} else {
  
  stop("Mismatch in Cavalli test set and test pheno Sample IDs")
}



Cavalli_bal_heatmap <- prob_heatmaps(
  prob_df           = Cavalli_results_bal$pred_probs,
  pred              = Cavalli_results_bal$pred_probs$pred,
  test_pheno        = Cavalli_pheno,
  sample_id_column  = "Sample_ID",
  true_test_labels  = "histology",
  mol_group         = "mol_grp",
  sub_group         = "sub_grp",
  MYC_status        = "MYC_Status",
  MYCN_status       = "MYCN_Status",
  OS                = "Dead",
  follow_up         = "OS (years)",
  Age               = "Age",
  continuum         = "DN/bi",
  heatmap_label     = "Binary Random Forest Model for classification of DN histopathology,\ntrained on Under-sampled Train NMB and tested on Cavalli",
  heatmap_panel     = FALSE,
  file_prefix       = "RF_balNMB_vs_Cavalli",
  save_in           = "~/Thesis_models/DN_models/Balanced_train_NMB/"
)

# ---- Balancing ablation: Train NMB (manually under-sampled), tested on both Test NMB and Cavalli ---- #
# This dataset is already balanced by construction, so EVERY config with bal_training=TRUE is
# excluded here. Only the three bal_training=FALSE configs are run: balanced_featsel_imbal_train,
# neither, and no_feature_selection_imbal_train.

US_NMB_ablation <- run_rf_ablation(
  train_betas    = US_train_NMB[, !colnames(US_train_NMB) %in% "bi_histo"],
  train_y        = as.factor(US_train_NMB$bi_histo),
  test_sets      = list(
    TestNMB = list(test_data = test_NMB, test_data_name = "Test NMB",
                   test_y = as.factor(test_NMB_pheno$bi_histo), heatmap_args = heatmap_args_testnmb),
    Cavalli = list(test_data = Cavalli_betas, test_data_name = "Cavalli",
                   test_y = as.factor(Cavalli_pheno$bi_histo), heatmap_args = heatmap_args_cavalli_nmb)
  ),
  train_set_name = "Train NMB (under-sampled)",
  save_in        = "~/Thesis_models/DN_models/Balanced_train_NMB/Ablation/",
  seed = 42, ntrees = 100, p = 200,
  configs        = Filter(function(cfg) !cfg$bal_training, ablation_configs_extended)
)

print(US_NMB_ablation$summary)

plot_ablation_comparison(
  US_NMB_ablation$summary,
  title     = "Macro AUC Across RF Balancing Configurations - Train NMB (under-sampled)",
  save_path = "~/Thesis_models/DN_models/Balanced_train_NMB/Ablation/ablation_comparison.png"
)

#-----Training Binary Random Forest Models on SHH-pathway activated MB patient group------ 

# Wrapper function to combine SHH samples from 'Infant' data with NMB MB-SHH samples, 
#to form the train set 'NMB+Infant-SHH'and to load Cavalli SHH test set. 


shh_pool <- function(age_probes = NULL,
                     nmb_pheno_path = "~/NMB/NMB_Final_Cohort_26.rds",
                     nmb_betas_path = "~/NMB/Final_NMB_betas.rds",
                     Inf_pheno_path = "~/NMB/updatedSampleSheet.csv",
                     Inf_betas_path = "~/NMB/Final_Inf_betas.rds",
                     Cavalli_betas_path = "~/Cavalli/Cavalli-test/Cavalli_test_betas.rds",
                     Cavalli_pheno_path = "~/Cavalli/Final_Cavalli_pheno.rds",
                     infant = FALSE) {
  
  # ---- NMB ----
  
  NMB_pheno <- readRDS(nmb_pheno_path)
  NMB_pheno <- NMB_pheno[NMB_pheno$CPR_Histology != "NOS", ]
  NMB_pheno <- NMB_pheno[!is.na(NMB_pheno$Sample_ID), ]
  NMB_pheno <- NMB_pheno %>%
    mutate(bi_histo = case_when(
      CPR_Histology %in% c("CLA", "LCA") ~ "non_DN",
      CPR_Histology %in% c("DN", "MBEN") ~ "DN"
    ))
  
  NMB_betas <- readRDS(nmb_betas_path)
  NMB_betas <- NMB_betas[match(NMB_pheno$Sample_ID, rownames(NMB_betas)), ]
  
  SHH_NMB_pheno <- NMB_pheno[NMB_pheno$Mol_group == "SHH-MB",]
  
  SHH_NMB_table<-data.frame("Sample_ID" = SHH_NMB_pheno$Sample_ID,
                            "CPR_Histology" = as.factor(SHH_NMB_pheno$CPR_Histology),
                            "Mol_group" = as.factor(SHH_NMB_pheno$Mol_grp),
                            "Sub_group" = as.factor(SHH_NMB_pheno$sub_group),
                            "Age" = SHH_NMB_pheno$Age,
                            "MYC_Status" = as.factor(SHH_NMB_pheno$`MYC Status`),
                            "MYCN_Status" = as.factor(SHH_NMB_pheno$`MYCN Status`),
                            "Bi_histo" = as.factor(SHH_NMB_pheno$bi_histo),
                            "OS" = as.factor(SHH_NMB_pheno$OS),
                            "Follow_up" = SHH_NMB_pheno$`Follow up years`,
                            "Source"   = "NMB")
  
  # ---- Cavalli ----
  
  Cavalli_pheno<-readRDS(Cavalli_pheno_path)
  
  Cavalli_betas<-readRDS(Cavalli_betas_path)
  
  Cavalli_betas<-as.data.frame(t(Cavalli_betas)) 
  
  if (any(is.na(Cavalli_pheno$Sample_ID))){
    Cavalli_pheno<-Cavalli_pheno[!is.na(Cavalli_pheno$Sample_ID),]
  }
  
  if (any(is.na(Cavalli_pheno$histology))){
    Cavalli_pheno<-Cavalli_pheno[!is.na(Cavalli_pheno$histology),]
  }
  
  # Age-matching with NMB Cohort
  
  Cavalli_pheno<-Cavalli_pheno %>% filter(Age <= max(NMB_pheno$Age))
  
  Cavalli_betas<-Cavalli_betas[match(Cavalli_pheno$Sample_ID, rownames(Cavalli_betas)),] 
  
  if(all(Cavalli_pheno$Sample_ID == rownames(Cavalli_betas))){
    Cavalli_pheno<-Cavalli_pheno %>% mutate(bi_histo = case_when( 
      histology  %in% c("CLA","LCA") ~ "non_DN", 
      histology  %in% c("DN", "MBEN") ~ "DN"))
  } else {
    message("Mismatch in Sample IDs in Pheno and betas")
    
  }
  
  Cavalli_SHH_pheno<-Cavalli_pheno[Cavalli_pheno$Mol_grp == "SHH-MB",]
  
  Cavalli_SHH_pheno<-Cavalli_SHH_pheno[!is.na(Cavalli_SHH_pheno$Mol_grp),]
  
  Cavalli_SHH_pheno<-droplevels(Cavalli_SHH_pheno)  # drop unused factor levels left over from subsetting the full Cavalli cohort down to SHH-MB only
  
  Cavalli_SHH_betas<-Cavalli_betas[match(Cavalli_SHH_pheno$Sample_ID,rownames(Cavalli_betas)),]
  
  
  # ---- Infant -----
  
  Inf_betas <- readRDS(Inf_betas_path)
  Inf_betas <- Inf_betas[, !colnames(Inf_betas) %in% c("mol.grp", "histo")]
  
  Inf_pheno <- read.csv(Inf_pheno_path)
  Inf_pheno <- Inf_pheno[Inf_pheno$Path_master != "NOS", ]
  
  Inf_pheno<-Inf_pheno[Inf_pheno$Nature_of_pathology_review == "Central",] 
  
  Inf_betas <- Inf_betas[match(Inf_pheno$Sentrix_ID, rownames(Inf_betas)), ]
  
  NMB_betas<-NMB_betas[,colnames(Inf_betas) %in% colnames(NMB_betas)]
  
  Inf_betas<-Inf_betas[,colnames(Inf_betas) %in% colnames(NMB_betas)] 
  
  
  if (all(rownames(Inf_betas) == Inf_pheno$Sentrix_ID)) {
    Inf_pheno$bi_histo <- as.factor(Inf_pheno$Path_master)
    Inf_pheno <- Inf_pheno %>%
      mutate(bi_histo = case_when(
        Path_master %in% c("CLA", "LCA") ~ "non_DN",
        Path_master %in% c("DN", "MBEN") ~ "DN"
      ))
  } else {
    Inf_pheno <- Inf_pheno[match(rownames(Inf_betas), Inf_pheno$Sentrix_ID), ]   # FIX: was lowercase inf_pheno (undefined)
    Inf_pheno <- Inf_pheno[!is.na(Inf_pheno$Sentrix_ID), ]
    Inf_betas <- Inf_betas[match(Inf_pheno$Sentrix_ID, rownames(Inf_betas)), ]
    
    if (all(rownames(Inf_betas) == Inf_pheno$Sentrix_ID)) {
      Inf_pheno$bi_histo <- as.factor(Inf_pheno$Path_master)
      Inf_pheno <- Inf_pheno %>%
        mutate(bi_histo = case_when(
          Path_master %in% c("CLA", "LCA") ~ "non_DN",
          Path_master %in% c("DN", "MBEN") ~ "DN"
        ))
    } else {
      stop("Mismatch in Inf_betas and Inf_pheno Sentrix/Sample_IDs")
    }
  }
  
  SHH_Inf_pheno <- Inf_pheno[Inf_pheno$Subgroup == "SHH", ]
  
  SHH_Inf_table<-data.frame("Sample_ID" = SHH_Inf_pheno$Sentrix_ID,
                            "CPR_Histology" = as.factor(SHH_Inf_pheno$Path_master),
                            "Sub_group" = as.factor(SHH_Inf_pheno$MNP_12.5_call),
                            "Age" = SHH_Inf_pheno$Age_at_diagnosis,
                            "MYC_Status" = as.factor(SHH_Inf_pheno$MYC_Amplified),
                            "MYCN_Status" = as.factor(SHH_Inf_pheno$MYCN_Amplified),
                            "Bi_histo" = as.factor(SHH_Inf_pheno$bi_histo),
                            "OS" = as.factor(SHH_Inf_pheno$Overall.Survival),
                            "Follow_up" = SHH_Inf_pheno$Overall.survival.time,
                            "Source" = "Infant")
  
  
  SHH_Inf_table$Mol_group <- "SHH-Inf-MB"
  
  
  
  # Remove any NMB-cohort and Cavalli-cohort samples already represented among the Infant samples
  if (any(SHH_Inf_pheno$Sentrix_ID %in% SHH_NMB_pheno$Sentrix_ID)) {
    SHH_NMB_pheno <- SHH_NMB_pheno[!SHH_NMB_pheno$Sentrix_ID %in% SHH_Inf_pheno$Sentrix_ID, ]
    SHH_NMB_pheno <- SHH_NMB_pheno[!is.na(SHH_NMB_pheno$Sentrix_ID), ]
    SHH_NMB_table <- SHH_NMB_table[SHH_NMB_table$Sample_ID %in% SHH_NMB_pheno$Sample_ID, ]
  }
  
  if (any(Cavalli_pheno$Sentrix_ID %in% SHH_Inf_pheno$Sentrix_ID)) {
    SHH_Inf_pheno <- SHH_Inf_pheno[!SHH_Inf_pheno$Sentrix_ID %in% Cavalli_pheno$Sentrix_ID, ]
    SHH_Inf_pheno <- SHH_Inf_pheno[!is.na(SHH_Inf_pheno$Sentrix_ID), ]
    SHH_Inf_table <- SHH_Inf_table[match(SHH_Inf_pheno$Sentrix_ID, SHH_Inf_table$Sample_ID),]
  }
  
  
  SHH_NMB_betas <- NMB_betas[match(SHH_NMB_table$Sample_ID, rownames(NMB_betas)), ]
  
  SHH_Inf_betas <- Inf_betas[match(SHH_Inf_table$Sample_ID, rownames(Inf_betas)),] 
  
  if (all(colnames(SHH_NMB_betas) == colnames(SHH_Inf_betas))) {
    Train_SHH <- rbind(SHH_NMB_betas, SHH_Inf_betas)
  } else {
    SHH_Inf_betas <- SHH_Inf_betas[, match(colnames(SHH_NMB_betas), colnames(SHH_Inf_betas))]
    Train_SHH <- rbind(SHH_NMB_betas, SHH_Inf_betas)
  }
  
  Train_SHH_pheno <- rbind(SHH_NMB_table, SHH_Inf_table)
  Train_SHH_pheno <- Train_SHH_pheno[match(rownames(Train_SHH), Train_SHH_pheno$Sample_ID), ]
  Train_SHH_pheno <- Train_SHH_pheno[!is.na(Train_SHH_pheno$Sample_ID), ]
  
  if (all(rownames(Train_SHH) == Train_SHH_pheno$Sample_ID)) {
    Train_SHH$bi_histo <- as.factor(Train_SHH_pheno$Bi_histo)
    Train_SHH <- Train_SHH[!is.na(Train_SHH$bi_histo), ]
  }
  
  Train_SHH_pheno <- Train_SHH_pheno[match(rownames(Train_SHH), Train_SHH_pheno$Sample_ID), ]
  Train_SHH_pheno <- droplevels(Train_SHH_pheno)  # drop unused factor levels left over from rbind() / row filtering
  
  # Removing age-correlated probes (using the previously-screened set)
  
  if (!is.null(age_probes)){
    Train_SHH_mat <- Train_SHH[, !colnames(Train_SHH) %in% "bi_histo"]
    Train_SHH_mat <- Train_SHH_mat[, !colnames(Train_SHH_mat) %in% age_probes] 
    
  } else {
    
    Train_SHH_mat <- Train_SHH[, !colnames(Train_SHH) %in% "bi_histo"]
    
  }
  
  if (infant == TRUE){
    
    Train_SHH_Inf_pheno <- Train_SHH_pheno[Train_SHH_pheno$Mol_group == "SHH-Inf-MB", ]
    match_idx <- match(Train_SHH_Inf_pheno$Sample_ID, rownames(Train_SHH_mat))
    Train_SHH_Inf_pheno <- Train_SHH_Inf_pheno[!is.na(match_idx), ]
    Train_SHH_Inf_pheno <- droplevels(Train_SHH_Inf_pheno)  # drop levels (e.g. "SHH-Old-MB") that no longer appear after subsetting to SHH-Inf-MB only
    Train_SHH_inf <- Train_SHH_mat[match_idx[!is.na(match_idx)], ]
    
    
    Cavalli_inf_pheno<-Cavalli_SHH_pheno[Cavalli_SHH_pheno$mol_grp == "SHH-Inf-MB",]
    
    Cavalli_inf_pheno<-Cavalli_inf_pheno[!is.na(Cavalli_inf_pheno$mol_grp),]
    Cavalli_inf_pheno<-droplevels(Cavalli_inf_pheno)  # same reasoning, for the Cavalli-side subset
    
    Cavalli_inf_betas<-Cavalli_SHH_betas[match(Cavalli_inf_pheno$Sample_ID, rownames(Cavalli_SHH_betas)),]
    
    list(Train_inf_pheno = Train_SHH_Inf_pheno, Train_inf_betas = Train_SHH_inf, Cavalli_inf_pheno = Cavalli_inf_pheno, Cavalli_inf_betas = Cavalli_inf_betas)
    
  } else {
    
    list (Train_SHH_pheno = Train_SHH_pheno, Train_SHH_betas = Train_SHH_mat, Cavalli_SHH_betas = Cavalli_SHH_betas, Cavalli_SHH_pheno = Cavalli_SHH_pheno)
  }
}

# ---Prepare SHH Data---- 

SHH_sets <- shh_pool(age_probes = NULL, infant = FALSE)

Train_SHH_betas<-SHH_sets$Train_SHH_betas

Train_SHH_pheno<-SHH_sets$Train_SHH_pheno

Train_SHH_pheno<-Train_SHH_pheno[match(rownames(Train_SHH_betas), Train_SHH_pheno$Sample_ID),]

Cavalli_SHH_betas<-SHH_sets$Cavalli_SHH_betas

Cavalli_SHH_pheno<-SHH_sets$Cavalli_SHH_pheno


# ------ Training methylClass (Liu, 2024) models with all 441870 input probes for DN Classification within SHH-MB ----- 

RF_comb_SHH <- maintrain(betas.. = as.data.frame(Train_SHH_betas[,!colnames(Train_SHH_betas) %in% "bi_histo"]),
                         y.. = as.factor(Train_SHH_pheno$Bi_histo),
                         method = "RF",
                         ntrees = 100, #default
                         p = 200, #default
                         topfeaturenumber = NULL,
                         subset.CpGs = NULL,
                         seed = 42,
                         calibrationmethod = c("MR","FLR","LR")
)


# Test the model on Train_SHH (itself)

message("Testing and evaluating on Train-SHH")

# (RF_comb_SHH is already in memory - trained above in this run.)

if(all(rownames(Train_SHH_betas) == Train_SHH_pheno$Sample_ID)){
  tryCatch(train_results_SHH<-test_and_evaluate_model(
    model = RF_comb_SHH,
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "'NMB'+'Infant' SHH-MB",
    test_data = Train_SHH_betas,
    test_data_name = "Train-SHH-MB",
    true_test_labels = as.factor(Train_SHH_pheno$Bi_histo),
    model_name = "Binary Random Forest",
    save_in = "~/Thesis_models/DN_models/SHH_NMB_Infant/"
  ), error = function(e){
    message(sprintf("FAILED [%s / %s]: %s",
                    "Binary Random Forest trained on Train NMB", "Train-SHH", conditionMessage(e)))
    NULL
  } )
  
} else {
  
  stop("Mismatch in Train SHH betas and pheno Sample IDs")
}


train_SHH_heatmap <- prob_heatmaps(
  prob_df           = train_results_SHH$pred_probs,
  pred              = train_results_SHH$pred_probs$pred,
  test_pheno        = Train_SHH_pheno,
  sample_id_column  = "Sample_ID",
  true_test_labels  = "CPR_Histology",
  mol_group         = "Mol_group",
  sub_group         = "Sub_group",
  MYC_status        = "MYC_Status",
  MYCN_status       = "MYCN_Status",
  OS                = "OS",
  follow_up         = "Follow_up",
  Age               = "Age",
  continuum         = "DN/bi",
  heatmap_label     = "Binary Random Forest Model for classification of DN histopathology,\ntrained on SHH-MB in 'NMB' and 'Infant',and tested on 'NMB+Infant' SHH-MB",
  heatmap_panel     = FALSE,
  file_prefix       = "RF_SHH_all_vs_TrainSHHall",
  save_in           = "~/Thesis_models/DN_models/SHH_NMB_Infant/"
)

# Test the model on Cavalli

message("Testing and evaluating on Cavalli")

if(all(rownames(Cavalli_SHH_betas) == Cavalli_SHH_pheno$Sample_ID)){
  tryCatch(Cavalli_results_SHH<-test_and_evaluate_model(
    model = RF_comb_SHH,
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "'NMB'+'Infant' SHH-MB",
    test_data = Cavalli_SHH_betas,
    test_data_name = "Cavalli-SHH-MB",
    true_test_labels = as.factor(Cavalli_SHH_pheno$bi_histo),
    model_name = "Binary Random Forest",
    save_in = "~/Thesis_models/DN_models/SHH_NMB_Infant/"
  ), error = function(e){
    message(sprintf("FAILED [%s / %s]: %s",
                    "Binary Random Forest trained on Train NMB", "Cavalli", conditionMessage(e)))
    NULL
  } )
  
} else {
  
  stop("Mismatch in Cavalli SHH test set and test pheno Sample IDs")
}



Cavalli_SHH_heatmap <- prob_heatmaps(
  prob_df           = Cavalli_results_SHH$pred_probs,
  pred              = Cavalli_results_SHH$pred_probs$pred,
  test_pheno        = Cavalli_SHH_pheno,
  sample_id_column  = "Sample_ID",
  true_test_labels  = "histology",
  mol_group         = "mol_grp",
  sub_group         = "sub_grp",
  MYC_status        = "MYC_Status",
  MYCN_status       = "MYCN_Status",
  OS                = "Dead",
  follow_up         = "OS (years)",
  Age               = "Age",
  continuum         = "DN/bi",
  heatmap_label     = "Binary Random Forest Model for classification of DN histopathology,\ntrained on SHH-MB in 'NMB' and 'Infant',and tested on Cavalli-SHH-MB",
  heatmap_panel     = FALSE,
  file_prefix       = "RF_SHHall_vs_CavalliSHHall",
  save_in           = "~/Thesis_models/DN_models/SHH_NMB_Infant/"
)

# ---- Balancing ablation: combined SHH-MB ('NMB'+'Infant', all ages, with age probes) ---- #

heatmap_args_trainshh <- list(
  test_pheno = Train_SHH_pheno, sample_id_column = "Sample_ID", true_test_labels = "CPR_Histology",
  mol_group = "Mol_group", sub_group = "Sub_group", MYC_status = "MYC_Status", MYCN_status = "MYCN_Status",
  OS = "OS", follow_up = "Follow_up", Age = "Age", continuum = "DN/bi"
)

heatmap_args_cavalli_shh <- list(
  test_pheno = Cavalli_SHH_pheno, sample_id_column = "Sample_ID", true_test_labels = "histology",
  mol_group = "mol_grp", sub_group = "sub_grp", MYC_status = "MYC_Status", MYCN_status = "MYCN_Status",
  OS = "Dead", follow_up = "OS (years)", Age = "Age", continuum = "DN/bi"
)

SHH_comb_ablation <- run_rf_ablation(
  train_betas    = as.data.frame(Train_SHH_betas[, !colnames(Train_SHH_betas) %in% "bi_histo"]),
  train_y        = as.factor(Train_SHH_pheno$Bi_histo),
  test_sets      = list(
    TrainSHH = list(test_data = Train_SHH_betas, test_data_name = "Train-SHH-MB",
                    test_y = as.factor(Train_SHH_pheno$Bi_histo), heatmap_args = heatmap_args_trainshh),
    Cavalli  = list(test_data = Cavalli_SHH_betas, test_data_name = "Cavalli-SHH-MB",
                    test_y = as.factor(Cavalli_SHH_pheno$bi_histo), heatmap_args = heatmap_args_cavalli_shh)
  ),
  train_set_name = "'NMB'+'Infant' SHH-MB",
  save_in        = "~/Thesis_models/DN_models/SHH_NMB_Infant/Ablation/",
  seed = 42, ntrees = 100, p = 200,
  configs        = ablation_configs_extended
)

print(SHH_comb_ablation$summary)

plot_ablation_comparison(
  SHH_comb_ablation$summary,
  title     = "Macro AUC Across RF Balancing Configurations - Combined SHH-MB",
  save_path = "~/Thesis_models/DN_models/SHH_NMB_Infant/Ablation/ablation_comparison.png"
)

# ---- Age distribution in DN and non-DN, and the age-related probe screen used to build the sans-age model ----

Train_SHH_pheno <- Train_SHH_pheno %>%
  group_by(Bi_histo) %>%
  mutate(q1 = quantile(Age, 0.25), q3 = quantile(Age, 0.75), iqr = q3 - q1,
         outliers =  Age > (q3 + 1.5 * iqr)) %>%
  ungroup()

p_age_train <- ggplot(Train_SHH_pheno,
                      aes(x = Bi_histo, y = Age, fill = Bi_histo)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.4, alpha = 0.8) +
  geom_jitter(aes(alpha = outliers), width = 0.15, size = 0.5) +  
  scale_fill_manual(values = c("DN" = "darkgreen", "non_DN" = "purple")) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.25), guide = "none") +
  stat_compare_means(method = "kruskal.test",
                     label.y = max(Train_SHH_pheno$Age) - 5,
                     size    = 2.5) +
  labs(x = "Histology group", y = "Age (in years)",
       title = wrap_title("Age Distribution among Histology Groups in Train 'NMB+Infant' SHH")) +
  theme(legend.position = "none",
        axis.text.x     = element_text(hjust = 0.5, size = 6),
        axis.title      = element_text(hjust = 0.5, face = "bold"),
        plot.title      = element_text(hjust = 0.5, face = "bold")) +
  coord_cartesian(clip = "off")


kw_excl <- kruskal.test(Age ~ Bi_histo, data = subset(Train_SHH_pheno, !outliers))

outlier_ids <- Train_SHH_pheno$Sample_ID[Train_SHH_pheno$outliers]

for (id in outlier_ids) {
  kw <- kruskal.test(Age ~ Bi_histo, data = subset(Train_SHH_pheno, Sample_ID != id))
  cat(sprintf("removing %s -> p = %.3g\n", id, kw$p.value))
} #Still significant

# ------ Training the RF model on feature set without chronological age related probes 

if (!all(rownames(Train_SHH_betas) == Train_SHH_pheno$Sample_ID)) {
  Train_SHH_pheno <- Train_SHH_pheno[match(rownames(Train_SHH_betas), Train_SHH_pheno$Sample_ID), ]
}

age_cor <- apply(Train_SHH_betas, 2, function(x) {
  ct <- cor.test(x, Train_SHH_pheno$Age, method = "spearman")
  c(rho = ct$estimate, p = ct$p.value)
})

age_cor_df <- as.data.frame(t(age_cor)) |>
  tibble::rownames_to_column("probe") |>
  dplyr::rename(rho = `rho.rho`) |>
  dplyr::mutate(p_adj = p.adjust(p, method = "BH")) |>
  dplyr::arrange(p_adj)

age_probes <- age_cor_df$probe[age_cor_df$p_adj < 0.05 & abs(age_cor_df$rho) > 0.3]

histo_probes <- rownames(RF_comb_SHH$mod[[1]]$importance)

message(sprintf("%d of the top-importance histology probes are also significantly age-correlated probes",
                length(which(histo_probes %in% age_probes))))


n_probes_before_age_removal <- ncol(Train_SHH_betas)

# Removing chronological age related probes 

Train_SHH_betas<-Train_SHH_betas[,!colnames(Train_SHH_betas) %in% age_probes]

#check 

length(age_probes)

n_probes_before_age_removal 

ncol(Train_SHH_betas) 


# Model

RF_comb_SHH_sansage <- maintrain(betas.. = Train_SHH_betas[,!colnames(Train_SHH_betas) %in% "bi_histo"],
                                 y.. = as.factor(Train_SHH_pheno$Bi_histo),
                                 method = "RF",
                                 ntrees = 100, #default
                                 p = 200, #default
                                 topfeaturenumber = NULL,
                                 subset.CpGs = NULL,
                                 seed = 42,
                                 calibrationmethod = c("MR","FLR","LR")
)

# Test the model on Train_SHH (itself)

message("Testing and evaluating on Train-SHH")

if(all(rownames(Train_SHH_betas) == Train_SHH_pheno$Sample_ID)){
  tryCatch(train_results_SHH_sa<-test_and_evaluate_model(
    model = RF_comb_SHH_sansage,
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "'NMB'+'Infant' SHH-MB with removed Age-related probes",
    test_data = Train_SHH_betas,
    test_data_name = "Train-SHH-MB",
    true_test_labels = as.factor(Train_SHH_pheno$Bi_histo),
    model_name = "Binary Random Forest",
    save_in = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/"
  ), error = function(e){
    message(sprintf("FAILED [%s / %s]: %s",
                    "Binary Random Forest trained on Train SHH", "Train-SHH", conditionMessage(e)))
    NULL
  } )
  
} else {
  
  stop("Mismatch in Train SHH betas and pheno Sample IDs")
} 


train_SHH_heatmap_sa <- prob_heatmaps(
  prob_df           = train_results_SHH_sa$pred_probs,
  pred              = train_results_SHH_sa$pred_probs$pred,
  test_pheno        = Train_SHH_pheno,
  sample_id_column  = "Sample_ID",
  true_test_labels  = "CPR_Histology",
  mol_group         = "Mol_group",
  sub_group         = "Sub_group",
  MYC_status        = "MYC_Status",
  MYCN_status       = "MYCN_Status",
  OS                = "OS",
  follow_up         = "Follow_up",
  Age               = "Age",
  continuum         = "DN/bi",
  heatmap_label     = "Binary Random Forest Model for classification of DN histopathology,\ntrained on SHH-MB in 'NMB' and 'Infant' with removed age-related probes,and tested on 'NMB+Infant' SHH-MB",
  heatmap_panel     = FALSE,
  file_prefix       = "RF_SHHall_sansage_vs_TrainSHHall",
  save_in           = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/"
)

# Test the model on Cavalli

message("Testing and evaluating on Cavalli")

if(all(rownames(Cavalli_SHH_betas) == Cavalli_SHH_pheno$Sample_ID)){
  tryCatch(Cavalli_results_SHH_sa<-test_and_evaluate_model(
    model = RF_comb_SHH_sansage,
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "'NMB'+'Infant' SHH-MB with removed Age-related probes",
    test_data = Cavalli_SHH_betas,
    test_data_name = "Cavalli-SHH-MB",
    true_test_labels = as.factor(Cavalli_SHH_pheno$bi_histo),
    model_name = "Binary Random Forest",
    save_in = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/"
  ), error = function(e){
    message(sprintf("FAILED [%s / %s]: %s",
                    "Binary Random Forest trained on SHH-MB", "Cavalli", conditionMessage(e)))
    NULL
  } )
  
} else {
  
  stop("Mismatch in Cavalli SHH test set and test pheno Sample IDs")
}


Cavalli_SHH_heatmap_sa <- prob_heatmaps(
  prob_df           = Cavalli_results_SHH_sa$pred_probs,
  pred              = Cavalli_results_SHH_sa$pred_probs$pred,
  test_pheno        = Cavalli_SHH_pheno,
  sample_id_column  = "Sample_ID",
  true_test_labels  = "histology",
  mol_group         = "mol_grp",
  sub_group         = "sub_grp",
  MYC_status        = "MYC_Status",
  MYCN_status       = "MYCN_Status",
  OS                = "Dead",
  follow_up         = "OS (years)",
  Age               = "Age",
  continuum         = "DN/bi",
  heatmap_label     = "Binary Random Forest Model for classification of DN histopathology,\ntrained on SHH-MB in 'NMB' and 'Infant' with removed age-related probes,and tested on Cavalli-SHH-MB",
  heatmap_panel     = FALSE,
  file_prefix       = "RF_SHHall_sansage_vs_CavalliSHHall",
  save_in           = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/"
) 

# ---- Balancing ablation: combined SHH-MB ('NMB'+'Infant', all ages, sans age-related probes) ---- #

SHH_comb_sansage_ablation <- run_rf_ablation(
  train_betas    = Train_SHH_betas[, !colnames(Train_SHH_betas) %in% "bi_histo"],
  train_y        = as.factor(Train_SHH_pheno$Bi_histo),
  test_sets      = list(
    TrainSHH = list(test_data = Train_SHH_betas, test_data_name = "Train-SHH-MB",
                    test_y = as.factor(Train_SHH_pheno$Bi_histo), heatmap_args = heatmap_args_trainshh),
    Cavalli  = list(test_data = Cavalli_SHH_betas, test_data_name = "Cavalli-SHH-MB",
                    test_y = as.factor(Cavalli_SHH_pheno$bi_histo), heatmap_args = heatmap_args_cavalli_shh)
  ),
  train_set_name = "'NMB'+'Infant' SHH-MB with removed Age-related probes",
  save_in        = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Ablation/",
  seed = 42, ntrees = 100, p = 200,
  configs        = ablation_configs_extended
)

print(SHH_comb_sansage_ablation$summary)

plot_ablation_comparison(
  SHH_comb_sansage_ablation$summary,
  title     = "Macro AUC Across RF Balancing Configurations - Combined SHH-MB (sans-age)",
  save_path = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Ablation/ablation_comparison.png"
)

# ----- Horvath (Horvath, 2013) Methylation Age Estimation and correlation -----

#BiocManager::install("methylclock")

require(methylclock) 


compute_horvath_accel <- function(betas_mat, pheno, sample_id_col, age_col) {
  
  probes_t   <- t(betas_mat)
  clocks     <- methylclock::DNAmAge(probes_t)
  clock_df   <- as.data.frame(clocks)
  clock_df$join_id <- as.character(clock_df$id)
  
  pheno$join_id <- as.character(pheno[[sample_id_col]])
  clock_merged  <- merge(pheno, clock_df, by = "join_id")
  
  clock_merged$horvath_accel     <- clock_merged$Horvath - clock_merged[[age_col]]
  lm_fit                          <- lm(Horvath ~ get(age_col), data = clock_merged)
  clock_merged$horvath_accel_res <- resid(lm_fit)
  
  clock_merged <- clock_merged[match(pheno$join_id, clock_merged$join_id), ]
  
  if (!all(pheno$join_id == clock_merged$join_id)) {
    warning("Check mismatch in Sample IDs when computing Horvath epigenetic clock")
  }
  
  pheno$Horvath_Age        <- clock_merged$Horvath
  pheno$Horvath_accel      <- clock_merged$horvath_accel
  pheno$Horvath_accel_res  <- clock_merged$horvath_accel_res
  pheno$join_id            <- NULL
  
  pheno
}

attach_clock_to_probs <- function(prob_df, pheno, sample_id_col, true_label_col,
                                  age_col, pred_col = "pred") {
  
  stopifnot(all(prob_df$Sample_ID == pheno[[sample_id_col]]))
  
  prob_df$True_label       <- as.factor(pheno[[true_label_col]])
  prob_df$Age               <- pheno[[age_col]]
  prob_df$Horvath_Age       <- pheno$Horvath_Age
  prob_df$Horvath_accel     <- pheno$Horvath_accel
  prob_df$Horvath_accel_res <- pheno$Horvath_accel_res
  
  prob_df %>%
    dplyr::mutate(Class = dplyr::case_when(
      True_label == "CLA"  & .data[[pred_col]] == "non_DN" ~ "Correct_non_DN",
      True_label == "LCA"  & .data[[pred_col]] == "non_DN" ~ "Correct_non_DN",
      True_label == "DN"   & .data[[pred_col]] == "DN"     ~ "Correct_DN",
      True_label == "MBEN" & .data[[pred_col]] == "DN"     ~ "Correct_DN",
      True_label == "DN"   & .data[[pred_col]] == "non_DN" ~ "DN->non_DN",
      True_label == "MBEN" & .data[[pred_col]] == "non_DN" ~ "DN->non_DN",
      True_label == "CLA"  & .data[[pred_col]] == "DN"     ~ "non_DN->DN",
      True_label == "LCA"  & .data[[pred_col]] == "DN"     ~ "non_DN->DN"
    ))
}


plot_prob_by_accel <- function(df, accel_col, prob_col = "DN", group_col = "Class",
                               group_colors, cohort_name, title_text, save_path) {
  
  label_x <- stats::quantile(df[[accel_col]], 0.8, na.rm = TRUE)
  
  p <- ggplot(df, aes(x = .data[[accel_col]], y = .data[[prob_col]], colour = .data[[group_col]])) +
    geom_point(size = 1, alpha = 0.7) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_color_manual(values = group_colors) +
    stat_cor(method      = "spearman",
             size        = 2.5,
             colour      = "black",
             label.x     = label_x,
             label.y     = 1,
             inherit.aes = FALSE,
             data        = df,
             mapping     = aes(x = .data[[accel_col]], y = .data[[prob_col]])) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = 0,   linetype = "dashed", colour = "grey50") +
    labs(x = "Horvath Epigenetic Age Acceleration (years)",
         y = paste0("P(", prob_col, ")"), colour = group_col,
         title = wrap_title(title_text)) +
    theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(hjust = 0.5, face = "bold")) +
    coord_cartesian(clip = "off")
  
  ggsave(save_path, p, width = 12, height = 7, dpi = 300, bg = "white")
  p
}


plot_accel_by_group <- function(df, accel_col, group_col, group_colors,
                                y_label, title_text, save_path, ref_line = 0) {
  
  accel_range <- range(df[[accel_col]], na.rm = TRUE)
  label_y     <- accel_range[2] + diff(accel_range) * 0.08
  
  p <- ggplot(df, aes(x = .data[[group_col]], y = .data[[accel_col]], fill = .data[[group_col]])) +
    geom_boxplot(outlier.size = 0.5, linewidth = 0.4, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
    scale_fill_manual(values = group_colors) +
    stat_compare_means(method = "kruskal.test", label.y = label_y, size = 2.5) +
    geom_hline(yintercept = ref_line, linetype = "dashed", colour = "grey50") +
    labs(x = "Histology group", y = y_label, title = wrap_title(title_text)) +
    theme(legend.position = "none",
          axis.text  = element_text(hjust = 0.5, face = "bold"),
          plot.title      = element_text(hjust = 0.5, face = "bold")) +
    coord_cartesian(clip = "off")
  
  ggsave(save_path, p, width = 12, height = 7, dpi = 300, bg = "white")
  p
}


plot_accel_vs_age <- function(df, accel_col, age_col, group_col, group_colors,
                              y_label, title_text, save_path, ref_line = 0) {
  
  label_x <- stats::quantile(df[[age_col]], 0.75, na.rm = TRUE)
  label_y <- max(df[[accel_col]], na.rm = TRUE) * 0.9
  
  p <- ggplot(df, aes(x = .data[[age_col]], y = .data[[accel_col]], colour = .data[[group_col]])) +
    geom_point(size = 1, alpha = 0.7) +
    geom_smooth(aes(group = 1), method = "lm", color = "#D85A30", se = TRUE) +
    scale_color_manual(values = group_colors) +
    stat_cor(method      = "spearman",
             label.x     = label_x,
             label.y     = label_y,
             size        = 2.5,
             colour      = "black",
             inherit.aes = FALSE,
             data        = df,
             mapping     = aes(x = .data[[age_col]], y = .data[[accel_col]])) +
    geom_hline(yintercept = ref_line, linetype = "dashed", colour = "grey50") +
    labs(x = "Chronological Age (years)", y = y_label, colour = "Histology Group",
         title = wrap_title(title_text)) +
    theme(axis.title = element_text(hjust = 0.5, face = "bold"),
          plot.title      = element_text(hjust = 0.5, face = "bold")) +
    coord_cartesian(clip = "off")
  
  ggsave(save_path, p, width = 12, height = 7, dpi = 300, bg = "white")
  p
}

CLASS_COLOURS_HORVATH <- c("Correct_DN" = "darkgreen", "Correct_non_DN" = "purple4",
                           "non_DN->DN" = "pink3",     "DN->non_DN"     = "orange")
HISTO_BIHISTO_COLOURS    <- c("DN" = "darkgreen", "non_DN" = "purple")

save_dir <- "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/"

# ---- Train-SHH ('NMB'+'Infant') ----

Train_SHH_pheno <- compute_horvath_accel(
  betas_mat     = Train_SHH_betas[, !colnames(Train_SHH_betas) %in% "bi_histo"],
  pheno         = Train_SHH_pheno,
  sample_id_col = "Sample_ID",
  age_col       = "Age"
)


Train_SHH_prob_df <- read_probs_csv(file.path(save_dir, "RF_SHH_sansage_Train-SHH-MB_probabilities_.csv"))

Train_SHH_prob_df <- attach_clock_to_probs(
  prob_df        = Train_SHH_prob_df,
  pheno          = Train_SHH_pheno,
  sample_id_col  = "Sample_ID",
  true_label_col = "CPR_Histology",
  age_col        = "Age"
)


plot_prob_by_accel(
  df           = Train_SHH_prob_df,
  accel_col    = "Horvath_accel",
  group_colors = CLASS_COLOURS_HORVATH,
  title_text   = "P(DN) by Horvath Epigenetic Age Acceleration in Train-SHH-MB",
  save_path    = file.path(save_dir, "Train_SHH_Horvath_accel.png")
)

plot_prob_by_accel(
  df           = Train_SHH_prob_df,
  accel_col    = "Horvath_accel_res",
  group_colors = CLASS_COLOURS_HORVATH,
  title_text   = "P(DN) by Residual Horvath Epigenetic Age Acceleration in Train-SHH-MB",
  save_path    = file.path(save_dir, "Train_SHH_Horvath_accel_res.png")
)

plot_accel_by_group(
  df           = Train_SHH_prob_df,
  accel_col    = "Horvath_accel",
  group_col    = "pred",
  group_colors = HISTO_BIHISTO_COLOURS,
  y_label      = "Horvath Methylation Age Acceleration (in years)",
  title_text   = "Horvath Epigenetic Age Acceleration in Histology Groups of Train 'NMB+Infant' SHH-MB",
  save_path    = file.path(save_dir, "Train_SHH_Horvath_accel_in_histo.png")
)

plot_accel_by_group(
  df           = Train_SHH_prob_df,
  accel_col    = "Horvath_accel_res",
  group_col    = "pred",
  group_colors = HISTO_BIHISTO_COLOURS,
  y_label      = "Horvath Methylation Age Acceleration (in years)",
  title_text   = "Residual Horvath Epigenetic Age Acceleration in Histology Groups of Train 'NMB+Infant' SHH-MB",
  save_path    = file.path(save_dir, "Train_SHH_Horvath_accel_res_in_histo.png")
)

plot_accel_vs_age(
  df           = Train_SHH_prob_df,
  accel_col    = "Horvath_accel",
  age_col      = "Age",
  group_col    = "pred",
  group_colors = HISTO_BIHISTO_COLOURS,
  y_label      = "Horvath Methylation Age Acceleration (in years)",
  title_text   = "Horvath Epigenetic Age Acceleration vs Chronological Age in Train 'NMB+Infant' SHH-MB",
  save_path    = file.path(save_dir, "Train_SHH_Horvath_accel_vs_age.png")
)


plot_accel_vs_age(
  df           = Train_SHH_prob_df,
  accel_col    = "Horvath_accel_res",
  age_col      = "Age",
  group_col    = "pred",
  group_colors = HISTO_BIHISTO_COLOURS,
  y_label      = "Horvath Methylation Age Acceleration (in years)",
  title_text   = "Residual Epigenetic Age Acceleration vs Chronological Age in Train 'NMB+Infant' SHH-MB",
  save_path    = file.path(save_dir, "Train_SHH_Horvath_accel_res_vs_age.png")
)

# ---- Cavalli-SHH ----

Cavalli_SHH_pheno$Sample_ID <- as.factor(Cavalli_SHH_pheno$Sample_ID)

Cavalli_SHH_pheno <- compute_horvath_accel(
  betas_mat     = Cavalli_SHH_betas,
  pheno         = Cavalli_SHH_pheno,
  sample_id_col = "Sample_ID",
  age_col       = "Age"
)

Cavalli_SHH_prob_df <- read_probs_csv(file.path(save_dir, "RF_SHH_sansage_Cavalli-SHH-MB_probabilities_.csv"))

Cavalli_SHH_prob_df <- attach_clock_to_probs(
  prob_df        = Cavalli_SHH_prob_df,
  pheno          = Cavalli_SHH_pheno,
  sample_id_col  = "Sample_ID",
  true_label_col = "histology",
  age_col        = "Age"
)

Cavalli_SHH_prob_df$bi_histo <- as.factor(Cavalli_SHH_pheno$bi_histo)


plot_prob_by_accel(
  df           = Cavalli_SHH_prob_df,
  accel_col    = "Horvath_accel",
  group_colors = CLASS_COLOURS_HORVATH,
  title_text   = "P(DN) by Horvath Epigenetic Age Acceleration in Cavalli SHH-MB",
  save_path    = file.path(save_dir, "Cavalli_Horvath_accel.png")
)

plot_prob_by_accel(
  df           = Cavalli_SHH_prob_df,
  accel_col    = "Horvath_accel_res",
  group_colors = CLASS_COLOURS_HORVATH,
  title_text   = "P(DN) by Residual Horvath Epigenetic Age Acceleration in Cavalli SHH-MB",
  save_path    = file.path(save_dir, "Cavalli_Horvath_accel_res.png")
)

plot_accel_by_group(
  df           = Cavalli_SHH_prob_df,
  accel_col    = "Horvath_accel",
  group_col    = "bi_histo",
  group_colors = HISTO_BIHISTO_COLOURS,
  y_label      = "Horvath Methylation Age Acceleration (in years)",
  title_text   = "Horvath Epigenetic Age Acceleration in Histology Groups of Cavalli SHH-MB",
  save_path    = file.path(save_dir, "Cavalli_Horvath_accel_in_histo.png")
)


plot_accel_by_group(
  df           = Cavalli_SHH_prob_df,
  accel_col    = "Horvath_accel_res",
  group_col    = "bi_histo",
  group_colors = HISTO_BIHISTO_COLOURS,
  y_label      = "Horvath Methylation Age Acceleration (in years)",
  title_text   = "Residual Horvath Epigenetic Age Acceleration in Histology Groups of Cavalli SHH-MB",
  save_path    = file.path(save_dir, "Cavalli_Horvath_accel_res_in_histo.png")
)

plot_accel_vs_age(
  df           = Cavalli_SHH_prob_df,
  accel_col    = "Horvath_accel",
  age_col      = "Age",
  group_col    = "bi_histo",
  group_colors = HISTO_BIHISTO_COLOURS,
  y_label      = "Horvath Methylation Age Acceleration (in years)",
  title_text   = "Horvath Epigenetic Age Acceleration vs Chronological Age in Cavalli SHH-MB",
  save_path    = file.path(save_dir, "Cavalli_Horvath_accel_vs_age.png")
)

plot_accel_vs_age(
  df           = Cavalli_SHH_prob_df,
  accel_col    = "Horvath_accel_res",
  age_col      = "Age",
  group_col    = "bi_histo",
  group_colors = HISTO_BIHISTO_COLOURS,
  y_label      = "Residual Horvath Methylation Age Acceleration (in years)",
  title_text   = "Residual Epigenetic Age Acceleration vs Chronological Age in Cavalli SHH-MB",
  save_path    = file.path(save_dir, "Cavalli_Horvath_accel_res_vs_age.png")
)

# ----- Training RF model on SHH-Infant (Ages 0- 4.99) ---- 

SHH_inf_sets <- shh_pool(age_probes = age_probes, infant = TRUE)

# shh_pool(infant=TRUE) returns fields named Train_inf_betas/Train_inf_pheno

Train_SHH_inf_betas<-SHH_inf_sets$Train_inf_betas

Train_SHH_inf_pheno<-SHH_inf_sets$Train_inf_pheno

Cavalli_inf_pheno<-SHH_inf_sets$Cavalli_inf_pheno

Cavalli_inf_betas<-SHH_inf_sets$Cavalli_inf_betas


Train_SHH_inf_betas <- Train_SHH_inf_betas[match(Train_SHH_inf_pheno$Sample_ID, rownames(Train_SHH_inf_betas)),]

RF_comb_SHH_inf <- maintrain(betas.. = as.data.frame(Train_SHH_inf_betas[,!colnames(Train_SHH_inf_betas) %in% "bi_histo"]),
                             y.. = as.factor(Train_SHH_inf_pheno$Bi_histo),
                             method = "RF", ntrees = 100, p = 200, topfeaturenumber = NULL,
                             seed = 42, subset.CpGs = NULL,
                             calibrationmethod = c("MR","FLR","LR"))

RF_Inf_NMB <- RF_comb_SHH_inf     

message("Testing and evaluating RF_Inf_MB on Train SHH (including all Mol. grps)")

if (all(rownames(Train_SHH_betas) == Train_SHH_pheno$Sample_ID)) {
  Train_SHH_inf_results <- tryCatch(
    test_and_evaluate_model(
      model = RF_Inf_NMB,
      method = "methylClass",
      calibrated = TRUE,
      calibration_method = "LR",
      train_set_name = "'NMB'+'Infant' SHH-Inf-MB with removed Age-related probes",
      test_data = Train_SHH_betas,
      test_data_name = "Train-SHH-MB",
      true_test_labels = as.factor(Train_SHH_pheno$Bi_histo),
      model_name = "Binary Random Forest",
      save_in = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/" 
    ),
    error = function(e) {
      message(sprintf("FAILED [%s / %s]: %s",
                      "Binary Random Forest trained on SHH-Inf", "Train SHH", conditionMessage(e)))
      NULL
    }
  )
} else {
  stop("Mismatch in test set and test pheno Sample IDs")
}


if (!is.null(Train_SHH_inf_results)) {
  Train_SHH_inf_heatmap <- prob_heatmaps(
    prob_df           = Train_SHH_inf_results$pred_probs,
    pred              = Train_SHH_inf_results$pred_probs$pred,
    test_pheno        = Train_SHH_pheno,
    sample_id_column  = "Sample_ID",
    true_test_labels  = "CPR_Histology",
    mol_group         = "Mol_group",
    sub_group         = "Sub_group",
    MYC_status        = "MYC_Status",
    MYCN_status       = "MYCN_Status",
    OS                = "OS",
    follow_up         = "Follow_up",  
    Age               = "Age",
    continuum         = "DN/bi",
    heatmap_panel     = FALSE,
    column_for_panel  = "DN",
    heatmap_label     = "Binary Random Forest Model for classification of DN histopathology,\ntrained on SHH-Inf-MB in 'NMB' and 'Infant' with removed age-related probes,and tested on 'NMB+Infant' SHH-MB",
    file_prefix       = "RF_SHHInf_sansage_vs_TrainSHHall",
    save_in           = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/"
  )
} else {
  message("Skipping Train_SHH_inf_heatmap - Train_SHH_inf_results evaluation failed upstream (see FAILED message above)")
}

message("Testing and evaluating RF_Inf_MB on Cavalli-SHH")

if(all(rownames(Cavalli_SHH_betas) == Cavalli_SHH_pheno$Sample_ID)){
  tryCatch(inf_Cavalli_results <-test_and_evaluate_model(
    model = RF_Inf_NMB,
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "'NMB'+'Infant' SHH-Inf-MB with removed Age-related probes",
    test_data = Cavalli_SHH_betas,
    test_data_name = "Cavalli-SHH-MB",
    true_test_labels = as.factor(Cavalli_SHH_pheno$bi_histo),
    model_name = "Binary Random Forest",
    save_in = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/"
  ), error = function(e){
    message(sprintf("FAILED [%s / %s]: %s",
                    "Binary Random Forest trained on SHH-Inf-NMB", "Cavalli_SHH", conditionMessage(e)))
    NULL
  } )
  
} else {
  stop("Mismatch in Cavalli test set and test pheno Sample IDs")
}


inf_Cavalli_heatmap <- prob_heatmaps(
  prob_df           = inf_Cavalli_results$pred_probs,
  pred              = inf_Cavalli_results$pred_probs$pred,
  test_pheno        = Cavalli_SHH_pheno,
  sample_id_column  = "Sample_ID",
  true_test_labels  = "histology",
  mol_group         = "mol_grp",
  sub_group         = "sub_grp",
  MYC_status        = "MYC_Status",
  MYCN_status       = "MYCN_Status",
  OS                = "Dead",
  follow_up         = "OS (years)",
  Age               = "Age",
  continuum         = "DN/bi",
  heatmap_label     = "Binary Random Forest Model for classification of DN histopathology,\ntrained on SHH-Inf-MB in 'NMB' and 'Infant' with removed age-related probes,and tested on Cavalli-SHH-MB",
  heatmap_panel     = FALSE,
  file_prefix       = "RF_SHHInf_sansage_vs_CavalliSHHall",
  save_in           = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/"
)



message("Testing and evaluating RF_Inf_MB on Cavalli-SHH-Inf")

if(all(rownames(Cavalli_inf_betas) == Cavalli_inf_pheno$Sample_ID)){
  tryCatch(Cavalli_inf_results <-test_and_evaluate_model(
    model = RF_Inf_NMB,
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "'NMB'+'Infant' SHH-Inf-MB with removed Age-related probes",
    test_data = Cavalli_inf_betas,
    test_data_name = "Cavalli-SHH-Inf-MB",
    true_test_labels = as.factor(Cavalli_inf_pheno$bi_histo),
    model_name = "Binary Random Forest",
    save_in = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/"
  ), error = function(e){
    message(sprintf("FAILED [%s / %s]: %s",
                    "Binary Random Forest trained on SHH-Inf-NMB", "Cavalli_SHH_Inf", conditionMessage(e)))
    NULL
  } )
  
} else {
  stop("Mismatch in Cavalli test set and test pheno Sample IDs")
}


Cavalli_inf_heatmap <- prob_heatmaps(
  prob_df           = Cavalli_inf_results$pred_probs,
  pred              = Cavalli_inf_results$pred_probs$pred,
  test_pheno        = Cavalli_inf_pheno,
  sample_id_column  = "Sample_ID",
  true_test_labels  = "histology",
  mol_group         = "mol_grp",
  sub_group         = "sub_grp",
  MYC_status        = "MYC_Status",
  MYCN_status       = "MYCN_Status",
  OS                = "Dead",
  follow_up         = "OS (years)",
  Age               = "Age",
  continuum         = "DN/bi",
  heatmap_label     = "Binary Random Forest Model for classification of DN histopathology,\ntrained on SHH-Inf-MB in 'NMB' and 'Infant' with removed age-related probes,and tested on Cavalli-SHH-Inf-MB",
  heatmap_panel     = FALSE,
  file_prefix       = "RF_SHHInf_sansage_vs_CavalliSHH_inf",
  save_in           = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/"
)


# Model Trained on SHH from all age groups 


if(all(rownames(Cavalli_inf_betas) == Cavalli_inf_pheno$Sample_ID)){
  tryCatch(inf_Cavalli_SHH_sa_results <-test_and_evaluate_model(
    model = RF_comb_SHH_sansage,
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "'NMB'+'Infant' SHH-MB with removed Age-related probes",
    test_data = Cavalli_inf_betas,
    test_data_name = "Cavalli-SHH-Inf-MB",
    true_test_labels = as.factor(Cavalli_inf_pheno$bi_histo),
    model_name = "Binary Random Forest",
    save_in = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/"
  ), error = function(e){
    message(sprintf("FAILED [%s / %s]: %s",
                    "Binary Random Forest trained on SHH-Inf-NMB", "Cavalli_SHH", conditionMessage(e)))
    NULL
  } )
  
} else {
  stop("Mismatch in Cavalli test set and test pheno Sample IDs")
}


inf_Cavalli_SHH_sa_heatmap <- prob_heatmaps(
  prob_df           = inf_Cavalli_SHH_sa_results$pred_probs,
  pred              = inf_Cavalli_SHH_sa_results$pred_probs$pred,
  test_pheno        = Cavalli_inf_pheno,
  sample_id_column  = "Sample_ID",
  true_test_labels  = "histology",
  mol_group         = "mol_grp",
  sub_group         = "sub_grp",
  MYC_status        = "MYC_Status",
  MYCN_status       = "MYCN_Status",
  OS                = "Dead",
  follow_up         = "OS (years)",
  Age               = "Age",
  continuum         = "DN/bi",
  heatmap_label     = "Binary Random Forest Model for classification of DN histopathology,\ntrained on SHH-MB in 'NMB' and 'Infant' with removed age-related probes,and tested on Cavalli-SHH-Inf-MB",
  heatmap_panel     = FALSE,
  file_prefix       = "RF_SHHall_sansage_vs_CavalliSHHInf",
  save_in           = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/"
)

# ---- Balancing ablation on the SHH-Inf-MB (sans-age) model ----

heatmap_args_cavalli_shhinf <- list(
  test_pheno = Cavalli_inf_pheno, sample_id_column = "Sample_ID", true_test_labels = "histology",
  mol_group = "mol_grp", sub_group = "sub_grp", MYC_status = "MYC_Status", MYCN_status = "MYCN_Status",
  OS = "Dead", follow_up = "OS (years)", Age = "Age", continuum = "DN/bi"
)

titration_save_dir <- "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/"

shh_inf_ablation <- run_rf_ablation(
  train_betas    = Train_SHH_inf_betas[, !colnames(Train_SHH_inf_betas) %in% "bi_histo"],
  train_y        = as.factor(Train_SHH_inf_pheno$Bi_histo),
  test_sets      = list(
    Cavalli_Inf = list(test_data = Cavalli_inf_betas, test_data_name = "Cavalli-SHH-Inf-MB",
                       test_y = as.factor(Cavalli_inf_pheno$bi_histo), heatmap_args = heatmap_args_cavalli_shhinf)
  ),
  train_set_name = "SHH-Inf-MB (sans-age)",
  save_in        = file.path(titration_save_dir, "Ablation_main_model"),
  seed = 42, ntrees = 100, p = 200,
  configs        = ablation_configs_extended
)

print(shh_inf_ablation$summary)

plot_ablation_comparison(
  shh_inf_ablation$summary,
  title     = "Macro AUC Across RF Balancing Configurations - SHH-Inf-MB",
  save_path = file.path(titration_save_dir, "Ablation_main_model", "ablation_comparison.png")
)


# =========================================================================================================================================#
# SAMPLE-SIZE PROJECTION WITH MACRO-AUC

#   1. Train the chosen RF configuration on many sub-sampled training sets (different sizes / class mixes).
#   2. Measure the AUC of each on held-out data. Two kinds of held-out data are used everywhere:
#        - repeated stratified k-fold CV inside the SHH pool
#        - the external Cavalli-SHH cohort
#   3. Fit a straight line of AUC against the LOG of the training count(s):
#        Curve 1  : AUC = b0 + b1 * ln(n_total)                               (DN:non-DN kept at the pool's natural ratio)
#        Curve 2  : AUC = b0 + b_DN * ln(n_DN) + b_nonDN * ln(n_nonDN)         (DN and non-DN counts varied independently, 4 x 4 grid)
#      b * ln(2) is the AUC gained each time that count doubles. Curve 2 only estimates these two coefficients (with bootstrap
#      intervals): which class's samples matter more for AUC. It does not project sample sizes.
#   4. Uncertainty = bootstrap over the TEST samples (resampled within class, same rows for every design), refitting each time.
## =========================================================================================================================================#


proj_dir <- "~/Thesis_models/DN_models/SHH_NMB_Infant/Sample_size_projection/"
dir.create(proj_dir, recursive = TRUE, showWarnings = FALSE)

# Which RF configuration(s) to project. NULL = pick the best by macro AUC on Cavalli from the sans-age SHH ablation

PROJECTION_CONFIGS_TO_RUN <- NULL

TARGET_AUCS      <- c(0.80, 0.85, 0.90, 0.93)   # AUCs to project a required sample size for (Curve 1)
PLAUSIBLE_EXTRAP <- 10                          # a projected n is only reported if it is <= this x the largest n tested
N_BOOTSTRAP      <- 1000
PROJ_SEED        <- 2026

# Random forest settings (same as the ablation runs)

RF_NTREES   <- 100
RF_P        <- 200
CALIBRATION <- "LR"

# Curve 1: fractions of the available training pool used (DN and non-DN scaled together)

SIZE_FRACS <- c(0.20, 0.30, 0.40, 0.50, 0.65, 0.80, 1.00)
CV_FOLDS   <- 5
CV_REPEATS <- 5      # repeats of the k-fold CV
CAV_REPS   <- 10     # random sub-samples per design when testing on Cavalli
STRATIFY_FOLDS_BY_SOURCE <- TRUE   # folds balanced on histology x cohort (NMB / Infant)

# Curve 2: every combination of DN fraction x non-DN fraction (4 x 4 = 16 designs)

RUN_CURVE2 <- TRUE
GRID_FRACS      <- c(0.25, 0.50, 0.75, 1.00)
GRID_CV_REPEATS <- 3
GRID_CAV_REPS   <- 5 

SRC_CV  <- paste0("CV (", CV_FOLDS, "-fold)")
SRC_CAV <- "Cavalli (external)"


# ---------------------------------------------------------------------------------------------------------------------------------------- #
# 2. DATA USED BY THE PROJECTION
# ---------------------------------------------------------------------------------------------------------------------------------------- #

stopifnot(all(rownames(Train_SHH_betas)  == Train_SHH_pheno$Sample_ID),
          all(rownames(Cavalli_SHH_betas) == as.character(Cavalli_SHH_pheno$Sample_ID)),
          !("bi_histo" %in% colnames(Train_SHH_betas)))

proj_y   <- factor(Train_SHH_pheno$Bi_histo, levels = c("DN", "non_DN"))
proj_ids <- as.character(Train_SHH_pheno$Sample_ID)
cav_y    <- factor(Cavalli_SHH_pheno$bi_histo, levels = c("DN", "non_DN"))
cav_ids  <- as.character(Cavalli_SHH_pheno$Sample_ID)

pool_n        <- length(proj_y)
pool_n_dn     <- sum(proj_y == "DN")
pool_n_non    <- sum(proj_y == "non_DN")
pool_share_dn <- pool_n_dn / pool_n
message(sprintf("Projection pool: %d samples = %d DN + %d non-DN (%.1f%% DN). Cavalli test: %d samples.",
                pool_n, pool_n_dn, pool_n_non, 100 * pool_share_dn, length(cav_y)))

# --- Strata used to sub-sample: histology. The DN and non-DN counts are controlled separately. ---

hist_strata <- proj_y

# --- Strata for the CV folds (so every fold has both histologies, and both cohorts where possible) ---

fold_strata <- proj_y
if (STRATIFY_FOLDS_BY_SOURCE && "Source" %in% colnames(Train_SHH_pheno)) {
  s <- interaction(proj_y, as.character(Train_SHH_pheno$Source), drop = TRUE)
  if (all(table(s) >= CV_FOLDS)) fold_strata <- s else
    message("Some histology x cohort cells have < ", CV_FOLDS, " samples - CV folds stratified by histology only.")
}

# --- Which config(s) to run ---
if (is.null(PROJECTION_CONFIGS_TO_RUN)) {
  PROJECTION_CONFIGS_TO_RUN <- choose_best_config(SHH_comb_sansage_ablation, metric = "macro_auc", test_set = "Cavalli")$best
}
projection_configs <- ablation_configs_extended[PROJECTION_CONFIGS_TO_RUN]
stopifnot(!anyNA(names(projection_configs)))


# ---------------------------------------------------------------------------------------------------------------------------------------- #
# 3. TRAIN ONE MODEL ON A SUB-SAMPLE AND SCORE THE TEST SAMPLES
# ---------------------------------------------------------------------------------------------------------------------------------------- #

draw_training_set <- function(strata, pool_idx, fractions, seed) {
  set.seed(seed)
  chosen <- lapply(levels(strata), function(s) {
    idx    <- pool_idx[strata[pool_idx] == s]
    f      <- if (s %in% names(fractions)) fractions[[s]] else 1
    n_take <- min(length(idx), max(1, round(f * length(idx))))
    idx[sample.int(length(idx), n_take)]
  })
  sort(unlist(chosen))
}

# Train one RF with the given config and return P(DN) for every row of x_test.

fit_and_score <- function(cfg, x_train, y_train, x_test, seed) {
  
  if (any(table(y_train) < 3)) stop("fewer than 3 training samples in one class")
  
  model <- maintrain(
    betas..           = as.data.frame(x_train),
    y..               = y_train,
    method            = "RF",
    ntrees            = RF_NTREES,
    p                 = RF_P,
    topfeaturenumber  = NULL,
    subset.CpGs       = NULL,
    seed              = seed,
    feature_sel       = cfg$feature_sel,
    bal_feature_sel   = cfg$bal_feature_sel,
    imp_replace       = cfg$imp_replace,
    bal_training      = cfg$bal_training,
    calibrationmethod = CALIBRATION
  )
  
  pred <- methylClass::mainpredict(newdat = x_test, mod = model$mod, calibratemod = model$platt.calfits)
  as.numeric(as.data.frame(pred$probs)[["DN"]])
}


# ---------------------------------------------------------------------------------------------------------------------------------------- #
# 4. RUN AN EXPERIMENT = every design x every split x every repeat
# ---------------------------------------------------------------------------------------------------------------------------------------- #

# designs   : data.frame, one row per training design, one column per stratum, values = fraction of that stratum to use
# strata    : factor over the pool saying which stratum each pool sample belongs to (column names of `designs` must be its levels)
# test_on   : "cv"      - k-fold CV inside the pool (the training sub-sample is drawn from each fold's training part)
#             "cavalli" - the sub-sample is drawn from the whole pool and tested on Cavalli-SHH
#
# Returns
#   n_info  : for each design, the mean number of training samples (total and per stratum) actually used
#   p_array : P(DN) for [test sample x design x repeat]  (NA where a fit failed)
#   is_dn   : TRUE for the DN test samples

run_experiment <- function(cfg, designs, strata, test_on = c("cv", "cavalli"), n_reps, seed0, label) {
  
  test_on      <- match.arg(test_on)
  stratum_cols <- colnames(designs)
  n_designs    <- nrow(designs)
  N            <- length(proj_y)
  test_truth   <- if (test_on == "cv") proj_y else cav_y
  
  p_array  <- array(NA_real_, dim = c(length(test_truth), n_designs, n_reps))
  size_log <- list()
  
  for (rep_i in seq_len(n_reps)) {
    
    message(sprintf("  [%s] repeat %d / %d", label, rep_i, n_reps))
    
    if (test_on == "cv") {
      set.seed(seed0 + rep_i)
      folds  <- caret::createFolds(fold_strata, k = CV_FOLDS, returnTrain = FALSE)
      splits <- lapply(folds, function(test_idx) list(test = sort(test_idx), train = setdiff(seq_len(N), test_idx)))
    } else {
      splits <- list(list(test = seq_along(cav_y), train = seq_len(N)))
    }
    
    for (s in seq_along(splits)) {
      for (d in seq_len(n_designs)) {
        
        seed_i    <- seed0 + rep_i * 10000 + s * 100 + d
        train_idx <- draw_training_set(strata, splits[[s]]$train, unlist(designs[d, ]), seed_i)
        test_rows <- splits[[s]]$test
        x_test    <- if (test_on == "cv") Train_SHH_betas[test_rows, , drop = FALSE] else Cavalli_SHH_betas
        
        p <- tryCatch(
          fit_and_score(cfg, Train_SHH_betas[train_idx, , drop = FALSE], proj_y[train_idx], x_test, seed_i),
          error = function(e) {
            message(sprintf("    FAILED [%s repeat=%d split=%d design=%d]: %s", label, rep_i, s, d, conditionMessage(e)))
            NULL
          })
        if (is.null(p)) next
        
        p_array[test_rows, d, rep_i] <- p
        
        n_used <- as.numeric(table(factor(strata[train_idx], levels = stratum_cols)))
        names(n_used) <- stratum_cols
        size_log[[length(size_log) + 1]] <- c(design = d, n_train = length(train_idx), n_used)
      }
    }
  }
  
  if (length(size_log) == 0) stop("Every fit failed in '", label, "' - see the FAILED messages above.")
  
  n_info <- as.data.frame(do.call(rbind, size_log)) %>%
    dplyr::group_by(design) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), mean), .groups = "drop") %>%
    dplyr::arrange(design)
  
  list(n_info  = n_info,
       p_array = p_array[, n_info$design, , drop = FALSE],
       is_dn   = test_truth == "DN")
}


run_both_tests <- function(cfg, designs, strata, cv_reps, cav_reps, seed_base, label) {
  out <- list(
    run_experiment(cfg, designs, strata, "cv",      cv_reps,  seed_base,          paste(label, "- CV")),
    run_experiment(cfg, designs, strata, "cavalli", cav_reps, seed_base + 500000, paste(label, "- Cavalli"))
  )
  names(out) <- c(SRC_CV, SRC_CAV)
  out
}


# ---------------------------------------------------------------------------------------------------------------------------------------- #
# 5. AUC, LINE FITTING AND BOOTSTRAP
# ---------------------------------------------------------------------------------------------------------------------------------------- #

auc_rank <- function(score, is_pos) {
  ok <- !is.na(score); score <- score[ok]; is_pos <- is_pos[ok]
  n1 <- sum(is_pos); n0 <- sum(!is_pos)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[is_pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# One AUC per design: computed for each repeat, then averaged over repeats.
auc_by_design <- function(p_array, is_dn) {
  a <- apply(p_array, c(2, 3), function(score) auc_rank(score, is_dn))
  a <- matrix(a, nrow = dim(p_array)[2])
  m <- rowMeans(a, na.rm = TRUE)
  m[is.nan(m)] <- NA_real_
  m
}

# Least-squares fit of AUC on the log counts. log_counts: one row per design, one column per term. Returns (intercept, slopes).
fit_log_linear <- function(auc, log_counts) {
  X   <- cbind("(Intercept)" = 1, log_counts)
  out <- setNames(rep(NA_real_, ncol(X)), colnames(X))
  ok  <- is.finite(auc)
  if (sum(ok) > ncol(X)) out[] <- lm.fit(X[ok, , drop = FALSE], auc[ok])$coefficients
  out
}

# Bootstrap over TEST samples: resample within class (same rows for every design), recompute every design's AUC, refit the line.

bootstrap_fits <- function(p_array, is_dn, log_counts, n_boot = N_BOOTSTRAP) {
  set.seed(PROJ_SEED)
  by_class <- list(which(is_dn), which(!is_dn))
  coefs <- matrix(NA_real_, n_boot, ncol(log_counts) + 1, dimnames = list(NULL, c("(Intercept)", colnames(log_counts))))
  aucs  <- matrix(NA_real_, n_boot, nrow(log_counts))
  for (b in seq_len(n_boot)) {
    rows       <- unlist(lapply(by_class, function(ix) ix[sample.int(length(ix), length(ix), replace = TRUE)]))
    aucs[b, ]  <- auc_by_design(p_array[rows, , , drop = FALSE], is_dn[rows])
    coefs[b, ] <- fit_log_linear(aucs[b, ], log_counts)
  }
  list(coef = coefs, auc = aucs)
}

# Everything about one fitted learning curve: observed AUCs, the fitted line, and its bootstrap.

fit_learning_curve <- function(exp, log_counts) {
  auc <- auc_by_design(exp$p_array, exp$is_dn)
  list(auc = auc, log_counts = log_counts, point = fit_log_linear(auc, log_counts),
       boot = bootstrap_fits(exp$p_array, exp$is_dn, log_counts))
}

# log of the mean training count(s) in the given n_info columns (rows = designs)
log_matrix <- function(n_info, cols) log(pmax(as.matrix(n_info[, cols, drop = FALSE]), 1))

ci_or_na <- function(v, min_valid = 100) {
  v <- v[!is.na(v)]
  if (length(v) >= min_valid) unname(quantile(v, c(0.025, 0.975))) else c(NA_real_, NA_real_)
}

fmt_num <- function(x, digits = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))


# ---------------------------------------------------------------------------------------------------------------------------------------- #
# 6. CURVE 1: tables
# ---------------------------------------------------------------------------------------------------------------------------------------- #

# Total n at which the fitted line reaches `target`. NA unless the slope is positive.

n_for_target <- function(b0, b1, target) ifelse(is.finite(b1) & b1 > 0, exp((target - b0) / b1), NA_real_)

curve1_fit_table <- function(fit, source) {
  n      <- exp(fit$log_counts[, 1])
  slopes <- fit$boot$coef[, 2]
  ci     <- ci_or_na(slopes)
  data.frame(source = source,
             intercept = fit$point[[1]], slope = fit$point[[2]],
             slope_ci_lower = ci[1], slope_ci_upper = ci[2],
             pct_boot_slope_positive = 100 * mean(slopes > 0, na.rm = TRUE),
             auc_per_doubling = fit$point[[2]] * log(2),
             auc_at_smallest_n = fit$auc[which.min(n)], auc_at_largest_n = fit$auc[which.max(n)],
             largest_n_tested = max(n), row.names = NULL)
}

# Projected total training size for each target AUC. Reported only if within PLAUSIBLE_EXTRAP x the largest n tested.

curve1_projection_table <- function(fit, source) {
  limit <- PLAUSIBLE_EXTRAP * max(exp(fit$log_counts[, 1]))
  b0 <- fit$point[[1]]; b1 <- fit$point[[2]]
  dplyr::bind_rows(lapply(TARGET_AUCS, function(target) {
    raw   <- n_for_target(b0, b1, target)
    acc   <- if (!is.na(raw) && raw <= limit) raw else NA_real_
    raw_b <- n_for_target(fit$boot$coef[, 1], fit$boot$coef[, 2], target)
    acc_b <- ifelse(raw_b <= limit, raw_b, NA_real_)
    ci    <- ci_or_na(acc_b)
    data.frame(source = source, target_auc = target,
               n_needed_uncapped = raw, n_needed = acc, ci_lower = ci[1], ci_upper = ci[2],
               pct_boot_valid = 100 * mean(!is.na(acc_b)),
               times_pool = acc / pool_n,
               dn_at_pool_mix = acc * pool_share_dn, non_dn_at_pool_mix = acc * (1 - pool_share_dn))
  }))
}

curve1_points_table <- function(fit, source) {
  ci <- t(apply(fit$boot$auc, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
  data.frame(source = source, n = exp(fit$log_counts[, 1]), auc = fit$auc, ci_lower = ci[, 1], ci_upper = ci[, 2])
}

report_curve1 <- function(fit_tab, proj_tab) {
  for (i in seq_len(nrow(fit_tab))) {
    f <- fit_tab[i, ]
    message(sprintf("\n  [%s] AUC %s -> %s over n = up to %d. Slope %s per ln(n) = %s AUC per doubling (95%% CI of slope %s to %s)",
                    f$source, fmt_num(f$auc_at_smallest_n), fmt_num(f$auc_at_largest_n), round(f$largest_n_tested),
                    fmt_num(f$slope), fmt_num(f$auc_per_doubling), fmt_num(f$slope_ci_lower), fmt_num(f$slope_ci_upper)))
    message("    ", if (!is.na(f$slope_ci_lower) && f$slope_ci_lower > 0)
      "AUC rises clearly with training size in the range tested; anything beyond the pool is extrapolation."
      else "slope is NOT clearly positive (CI includes 0): more data has not been shown to help in this range.")
    print(proj_tab[proj_tab$source == f$source,
                   c("target_auc", "n_needed", "ci_lower", "ci_upper", "times_pool", "dn_at_pool_mix", "non_dn_at_pool_mix")],
          row.names = FALSE, digits = 3)
  }
  message("  (NA = not estimable: slope not positive, or n beyond ", PLAUSIBLE_EXTRAP, "x the largest size tested)")
}


# ---------------------------------------------------------------------------------------------------------------------------------------- #
# 7. CURVE 2: coefficients and partial-effect curves
# ---------------------------------------------------------------------------------------------------------------------------------------- #

# One row per term (DN count and non-DN count): coefficient, bootstrap CI, and AUC gained per doubling.

coefficient_table <- function(fit, analysis, source) {
  terms <- colnames(fit$log_counts)
  est   <- unname(fit$point[terms])
  boot  <- fit$boot$coef[, terms, drop = FALSE]
  ci    <- t(sapply(terms, function(tm) ci_or_na(boot[, tm])))
  data.frame(analysis = analysis, source = source, term = terms,
             estimate = est, ci_lower = ci[, 1], ci_upper = ci[, 2],
             pct_boot_positive = 100 * colMeans(boot > 0, na.rm = TRUE),
             auc_per_doubling = est * log(2),
             ci_lower_per_doubling = ci[, 1] * log(2), ci_upper_per_doubling = ci[, 2] * log(2),
             row.names = NULL)
}

# "Partial effect" of each term = the fitted line for that term with all other terms adjusted away.
#   points : observed AUC of each design, adjusted for the other terms (line + residual)
#   lines  : fitted line (slope = the coefficient) with a band from the bootstrap slopes

partial_effects <- function(fit, analysis, source) {
  
  if (anyNA(fit$point)) return(list(points = NULL, lines = NULL))
  
  X     <- fit$log_counts
  b     <- fit$point[colnames(X)]
  fitted <- as.vector(cbind(1, X) %*% fit$point)
  resid  <- fit$auc - fitted
  centre <- mean(fit$auc, na.rm = TRUE)
  pts <- list(); lns <- list()
  
  for (j in seq_len(ncol(X))) {
    x     <- X[, j]
    xs    <- seq(min(x), max(x), length.out = 60)
    slope <- fit$boot$coef[, colnames(X)[j]]
    band  <- sapply(xs, function(v) quantile(slope * (v - mean(x)), c(0.025, 0.975), na.rm = TRUE))
    pts[[j]] <- data.frame(analysis = analysis, source = source, term = colnames(X)[j], count = exp(x),
                           partial_auc = centre + b[[j]] * (x - mean(x)) + resid)
    lns[[j]] <- data.frame(analysis = analysis, source = source, term = colnames(X)[j], count = exp(xs),
                           line = centre + b[[j]] * (xs - mean(x)),
                           lo = centre + band[1, ], hi = centre + band[2, ])
  }
  list(points = dplyr::bind_rows(pts), lines = dplyr::bind_rows(lns))
}


analyse_counts <- function(exps, cols, analysis) {
  fits  <- lapply(exps, function(e) fit_learning_curve(e, log_matrix(e$n_info, cols)))
  coefs <- dplyr::bind_rows(lapply(names(fits), function(s) coefficient_table(fits[[s]], analysis, s)))
  pe    <- lapply(names(fits), function(s) partial_effects(fits[[s]], analysis, s))
  list(fits = fits, coefs = coefs,
       partial_points = dplyr::bind_rows(lapply(pe, `[[`, "points")),
       partial_lines  = dplyr::bind_rows(lapply(pe, `[[`, "lines")))
}

report_coefficients <- function(coefs) {
  message("\n  AUC gained per doubling of each count (positive = more of that class raises AUC):")
  print(coefs[, c("analysis", "source", "term", "auc_per_doubling", "ci_lower_per_doubling",
                  "ci_upper_per_doubling", "pct_boot_positive")], row.names = FALSE, digits = 3)
  message("  (a CI that includes 0 means the data cannot show that count helps; pct_boot_positive = % of bootstrap fits with a positive coefficient)")
}


# ---------------------------------------------------------------------------------------------------------------------------------------- #
# 8. PLOTS
# ---------------------------------------------------------------------------------------------------------------------------------------- #

# Curve 1 learning curves. Points = observed AUC (bar = bootstrap 95% CI), blue line + band = the fitted log-linear curve,
# dashed vertical = current pool size, grey area = beyond the largest n tested (extrapolation), dotted horizontals = AUC targets.

plot_curve1 <- function(fits, title, save_path) {
  
  pts <- list(); lns <- list(); rects <- list()
  
  for (src in names(fits)) {
    f      <- fits[[src]]
    n      <- exp(f$log_counts[, 1])
    n_seq  <- exp(seq(log(min(n)), log(PLAUSIBLE_EXTRAP * max(n)), length.out = 200))
    boot_l <- outer(f$boot$coef[, 1], rep(1, 200)) + outer(f$boot$coef[, 2], log(n_seq))
    band   <- apply(boot_l, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    
    pts[[src]]   <- curve1_points_table(f, src)
    lns[[src]]   <- data.frame(source = src, n = n_seq, auc = pmin(1, f$point[[1]] + f$point[[2]] * log(n_seq)),
                               lo = pmin(1, band[1, ]), hi = pmin(1, band[2, ]))
    rects[[src]] <- data.frame(source = src, xmin = max(n))
  }
  pts <- dplyr::bind_rows(pts); lns <- dplyr::bind_rows(lns); rects <- dplyr::bind_rows(rects)
  
  p <- ggplot() +
    geom_rect(data = rects, aes(xmin = xmin, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "grey85", alpha = 0.4) +
    geom_hline(yintercept = TARGET_AUCS, colour = "grey60", linetype = "dotted") +
    geom_vline(xintercept = pool_n, colour = "grey30", linetype = "dashed") +
    geom_ribbon(data = lns, aes(x = n, ymin = lo, ymax = hi), fill = "#378ADD", alpha = 0.2) +
    geom_line(data = lns, aes(x = n, y = auc), colour = "#378ADD", linewidth = 0.8) +
    geom_linerange(data = pts, aes(x = n, ymin = ci_lower, ymax = ci_upper), alpha = 0.6) +
    geom_point(data = pts, aes(x = n, y = auc), size = 2.8) +
    facet_wrap(~source) +
    scale_x_log10(labels = scales::comma) +
    coord_cartesian(ylim = c(0.4, 1)) +
    labs(x = "Training samples (natural DN:non-DN ratio, log scale)", y = "AUC", title = wrap_title(title),
         subtitle = paste0("Points = observed (bar: bootstrap 95% CI); line = fitted AUC = b0 + b1 ln(n); dashed = current pool (n = ",
                           pool_n, "); grey = beyond largest n tested")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 8),
          strip.text = element_text(face = "bold"))
  
  ggsave(save_path, p, width = 11, height = 5.5, dpi = 300, bg = "white")
  p
}

# Curve 2 partial-effect plots: one panel per term. x = that count, y = AUC adjusted for the other terms.
# The slope of the line IS the coefficient reported in the table.


plot_partial_effects <- function(res, title, save_path) {
  
  if (nrow(res$partial_lines) == 0) { message("  (partial-effect plot skipped: fit not estimable)"); return(invisible(NULL)) }
  
  pts <- res$partial_points; lns <- res$partial_lines
  pts$term_label <- gsub("_", "-", pts$term); lns$term_label <- gsub("_", "-", lns$term)
  
  p <- ggplot() +
    geom_ribbon(data = lns, aes(x = count, ymin = lo, ymax = hi), fill = "#378ADD", alpha = 0.2) +
    geom_line(data = lns, aes(x = count, y = line), colour = "#378ADD", linewidth = 0.8) +
    geom_point(data = pts, aes(x = count, y = partial_auc), alpha = 0.7, size = 1.8) +
    facet_grid(source ~ term_label, scales = "free_x") +
    scale_x_log10() +
    labs(x = "Training samples of this class (log scale)", y = "AUC adjusted for the other terms",
         title = wrap_title(title),
         subtitle = "Points = designs; line = fitted effect of this count alone (slope = coefficient); band = bootstrap 95% CI of the slope") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 8),
          strip.text = element_text(face = "bold"))
  
  ggsave(save_path, p, width = max(8, 2.6 * length(unique(pts$term)) + 2), height = 6, dpi = 300, bg = "white")
  p
}

# Coefficients of the DN and non-DN counts, as AUC gained per doubling.

plot_coefficients <- function(coefs, title, save_path) {
  
  d <- coefs[!is.na(coefs$estimate), ]
  if (nrow(d) == 0) { message("  (coefficient plot skipped: nothing estimable)"); return(invisible(NULL)) }
  d$term_label <- gsub("_", "-", d$term)
  
  p <- ggplot(d, aes(x = auc_per_doubling, y = term_label, colour = source)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_pointrange(aes(xmin = ci_lower_per_doubling, xmax = ci_upper_per_doubling), position = position_dodge(width = 0.5)) +
    facet_wrap(~analysis, scales = "free_y", ncol = 1) +
    labs(x = "AUC gained per doubling of the count (bootstrap 95% CI)", y = NULL, colour = "Test data", title = wrap_title(title)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), strip.text = element_text(face = "bold"))
  
  ggsave(save_path, p, width = 8, height = 3 + 1.4 * length(unique(d$term_label)), dpi = 300, bg = "white")
  p
}


# The grid itself: observed AUC for every DN count x non-DN count (what the regression is fitted to).

plot_grid_heatmap <- function(res, title, save_path) {
  
  cells <- dplyr::bind_rows(lapply(names(res$fits), function(src) {
    f <- res$fits[[src]]
    data.frame(source = src, n_dn = exp(f$log_counts[, "DN"]), n_non = exp(f$log_counts[, "non_DN"]), auc = f$auc)
  }))
  
  p <- ggplot(cells, aes(x = factor(round(n_dn)), y = factor(round(n_non)), fill = auc)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = ifelse(is.na(auc), "NA", sprintf("%.2f", auc))), size = 3) +
    facet_wrap(~source) +
    scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0.5, 1), name = "AUC") +
    labs(x = "Mean DN training samples", y = "Mean non-DN training samples", title = wrap_title(title)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), strip.text = element_text(face = "bold"))
  
  ggsave(save_path, p, width = 11, height = 5.5, dpi = 300, bg = "white")
  p
}

# ---------------------------------------------------------------------------------------------------------------------------------------- #
# 9. RUN
# ---------------------------------------------------------------------------------------------------------------------------------------- #

run_projection <- function(cfg_name, cfg) {
  
  out_dir <- file.path(proj_dir, cfg_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ---------------- Curve 1: total n, natural DN:non-DN ratio ---------------- #
  
  message("\nCurve 1: AUC vs total training size")
  designs1 <- data.frame(DN = SIZE_FRACS, non_DN = SIZE_FRACS)
  exps1    <- run_both_tests(cfg, designs1, hist_strata, CV_REPEATS, CAV_REPS, seed_base = PROJ_SEED, label = "Curve 1")
  fits1    <- lapply(exps1, function(e) fit_learning_curve(e, log_matrix(e$n_info, "n_train")))
  
  points1 <- dplyr::bind_rows(lapply(names(fits1), function(s) curve1_points_table(fits1[[s]], s)))
  fit1    <- dplyr::bind_rows(lapply(names(fits1), function(s) curve1_fit_table(fits1[[s]], s)))
  proj1   <- dplyr::bind_rows(lapply(names(fits1), function(s) curve1_projection_table(fits1[[s]], s)))
  
  write.csv(points1, file.path(out_dir, "curve1_points.csv"),     row.names = FALSE)
  write.csv(fit1,    file.path(out_dir, "curve1_fit.csv"),        row.names = FALSE)
  write.csv(proj1,   file.path(out_dir, "curve1_projection.csv"), row.names = FALSE)
  plot_curve1(fits1, paste0("AUC vs training size - overall SHH-MB [", cfg_name, "]"), file.path(out_dir, "curve1_learning_curves.png"))
  report_curve1(fit1, proj1)
  
  result <- list(config = cfg_name, curve1 = list(fits = fits1, points = points1, fit = fit1, projection = proj1))
  
  # ---------------- Curve 2: which counts drive the AUC ---------------- #
  
  if (RUN_CURVE2) {
    
    message("\nCurve 2: DN and non-DN counts varied independently")
    designs2a <- expand.grid(DN = GRID_FRACS, non_DN = GRID_FRACS)
    exps2a    <- run_both_tests(cfg, designs2a, hist_strata, GRID_CV_REPEATS, GRID_CAV_REPS,
                                seed_base = PROJ_SEED + 1000000, label = "Curve 2")
    res2a     <- analyse_counts(exps2a, cols = c("DN", "non_DN"), analysis = "DN vs non-DN counts")
    plot_partial_effects(res2a, paste0("Effect of DN and non-DN training counts on AUC [", cfg_name, "]"),
                         file.path(out_dir, "curve2_partial_effects.png"))
    plot_grid_heatmap(res2a, paste0("AUC across DN / non-DN training counts [", cfg_name, "]"),
                      file.path(out_dir, "curve2_grid_heatmap.png"))
    
    all_coefs <- res2a$coefs
    result$curve2 <- res2a
    
    write.csv(all_coefs, file.path(out_dir, "curve2_coefficients.csv"), row.names = FALSE)
    plot_coefficients(all_coefs, paste0("AUC gained per doubling of each training count [", cfg_name, "]"),
                      file.path(out_dir, "curve2_coefficients.png"))
    report_coefficients(all_coefs)
    result$curve2_coefficients <- all_coefs
  }
  
  saveRDS(result, file.path(out_dir, "projection_result.rds"))   # output for re-plotting; never read back by this script
  result
}

projection_results <- lapply(names(projection_configs), function(nm) run_projection(nm, projection_configs[[nm]]))
names(projection_results) <- names(projection_configs)

write.csv(dplyr::bind_rows(lapply(names(projection_results), function(nm) cbind(config = nm, projection_results[[nm]]$curve1$projection))),
          file.path(proj_dir, "curve1_projection_all_configs.csv"), row.names = FALSE)
