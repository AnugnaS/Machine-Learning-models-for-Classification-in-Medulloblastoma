#=========================================================================================================================================#
# Training Binary Random Forest Models for Classification of Desmoplastic Nodular Histology in Medulloblastoma.

# This script includes loading of samplesheets and pre-processed methylation data (in-house),
# Initial training of Binary Random Forest models for classification of DN and non-DN MB groups across all principal Molecular groups
# Training Binary Random Forest models for DN histology prediction within SHH-pathway activated MB patient group
# Horvath (Horvath, 2013) Epigenetic Age prediction and Epigenetic Age acceleration
# Statistical testing for significant difference and correlation 
# Gene and Pathway Enrichment Analyses 
# Differential Gene Expression Analysis 
# Heatmaps, Bar plots, Scatter plots and density plots to visualise output probability score distribution
#==========================================================================================================================================#

# ------ Load packages ------ #

require(dplyr)

require(methylClass)

require(caret)

require(MLmetrics)

require(pROC)

require(circlize)

require(ComplexHeatmap)

require(ggplot2)

require(ggpubr)

require(cowplot)

require(ComplexHeatmap) 


# ------ Custom Functions ------ #

# Wrapper function for caret models - confusion matrix, Receiver Operating Characteristic (ROC) curve plots
# Area under the Curve (AUC) scores, F1 scores and Kappa scores

test_and_evaluate_model <- function(model, method, test_data, true_test_labels, test_data_name, model_name,
                                    train_set_name,calibrated = TRUE, calibration_method = NULL,
                                    save_in = "~/Thesis_models/DN_models/") {
  
  dir.create(save_in, recursive = TRUE, showWarnings = FALSE)
  
  # Shortform to save files
  
  model_shorthand<-ifelse(grepl("Random Forest", model_name),"RF",
                          ifelse(grepl("Support Vector Machine", model_name), "SVM", "XGB")) 
  
  train_set_shorthand <- ifelse(
    grepl("SHH", train_set_name) & grepl("Age-related", train_set_name, ignore.case = TRUE),
    "SHH_sansage",
    ifelse(grepl("SHH", train_set_name), "SHH", train_set_name)
  ) 
  
  # To check required packages
  
  required_pgs <- c("caret", "pROC", "ggplot2","tibble","stringr") 
  
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
        
      }
    } else {
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
  
  HISTO_COLOURS <- c("DN" = "darkgreen", "non_DN" = "purple") # "CLA" = "purple", "LCA" = "pink", "MBEN" = "lightgreen"
  
  cm_plot <- ggplot(cm_df, aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = Reference, alpha = Freq), color = "white") +
    geom_text(aes(label = Freq), size = 4,
              color = ifelse(cm_df$Freq > max(cm_df$Freq) * 0.5, "white", "black"), fontface = "bold") +
    scale_fill_manual(values = HISTO_COLOURS, guide = "none") +
    scale_alpha_continuous(range = c(0.15, 1), guide = "none") +
    labs( title = str_wrap(
      paste0("Confusion Matrix for ", model_name, " trained on ", train_set_name, ", tested on ", test_data_name),
      width = 50 ),
      x = "True Class", y = "Predicted Class"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          title = element_text(face = "bold"),
          axis.title = element_text ( hjust = 0.5,face = "bold"),
          plot.margin = ggplot2::margin(t = 15, r = 15, b = 15, l = 15)) +
    coord_cartesian(clip = "off")
  
  ggsave(file.path(save_in, paste0(model_shorthand, "_",train_set_shorthand,"_", test_data_name, "_confusion_matrix.png")),
         cm_plot, width = 9, height = 6, dpi = 300, bg = "white")
  
  #--- Overall Accuracy
  
  accuracy <- unname(cm$overall["Accuracy"])
  
  #---Kappa
  
  kappa <- unname(cm$overall["Kappa"])
  
  #---Precision and Recall- Class wise and Macro
  
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
    precision <- unname(byclass["Precision"])  
    recall <- unname(byclass["Recall"])
    
  } else {
    precision <- mean(byclass[, "Precision"], na.rm = TRUE)
    recall <- mean(byclass[, "Recall"], na.rm = TRUE)
    sensitivity <- byclass[, "Sensitivity"]
    specificity <- byclass[, "Specificity"]
    names(sensitivity) <- gsub("^Class: ", "", rownames(byclass))
    names(specificity) <- gsub("^Class: ", "", rownames(byclass))
  }
  
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
    scale_color_manual(values = HISTO_COLOURS,
                       labels = setNames(roc_df$class_label, roc_df$class)[unique(roc_df$class)]) +
    labs(title = str_wrap(paste("ROC Curves for", model_name, " trained on ",train_set_name, ",", "\n", 
                                "tested on", test_data_name), width = 50),
         subtitle = paste0("Macro AUC (Hand-Till) = ", round(macro_auc, 3)),
         x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate (Sensitivity)") +
    theme_minimal() +
    theme(title = element_text(face = "bold"),
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
    scale_color_manual(values = HISTO_COLOURS,
                       labels = setNames(pr_df$class_label,pr_df$class)[unique(pr_df$class)]) +
    labs(title = str_wrap(paste("Precision-Recall Curves for", model_name, " trained on ",
                                train_set_name, ",", "\n", "tested on", test_data_name), width = 50),
         subtitle = paste0("Macro PR-AUC = ", round(macro_pr_auc, 3)),
         x = "Recall", y = "Precision", color = "Class") +
    theme_minimal() + ylim(0, 1) + xlim(0, 1) + 
    theme(title = element_text(face = "bold"),
          axis.title = element_text ( hjust = 0.5,face = "bold"),
          plot.margin = ggplot2::margin(t = 15, r = 15, b = 15, l = 15)) +
    coord_cartesian(clip = "off")
  
  ggsave(file.path(save_in, paste0(model_shorthand, "_",train_set_shorthand,"_", test_data_name, "_pr_curves.png")),
         pr_plot, width = 9, height = 6, dpi = 300, bg = "white")
  
  #---Summary
  
  pred_probs$pred <- as.factor(pred_labels)
  
  pred_probs %>% 
    rownames_to_column(var = "Sample_ID") %>%
    write.csv(file.path(save_in, paste0(model_shorthand, "_", train_set_shorthand, "_", 
                                        test_data_name, "_probabilities_.csv")), 
              row.names = FALSE)
  
  summary_df <- data.frame(
    model = model_shorthand,
    test_data = test_data_name,
    macro_precision = precision,
    macro_recall = recall,
    macro_auc = macro_auc,
    macro_pr_auc = macro_pr_auc
  )
  
  #---Class-wise metrics
  
  classwise_df <- data.frame(
    model = model_shorthand,
    test_data = test_data_name,
    class = classes, 
    sensitivity = as.numeric(sensitivity[classes]),
    specificity = as.numeric(specificity[classes]),
    auc = as.numeric(class_auc[classes]),
    pr_auc = as.numeric(pr_auc[classes])
  )
  
  #---Save as .csv files
  
  write.csv(summary_df, file.path(save_in, paste0(model_shorthand, "_", train_set_shorthand, "_", test_data_name, "_summary_metrics.csv")), row.names = FALSE)
  write.csv(classwise_df, file.path(save_in, paste0(model_shorthand, "_", train_set_shorthand, "_", test_data_name, "_classwise_metrics.csv")), row.names = FALSE)
  
  
  #---Return to object 
  
  list(
    confusion_matrix = cm,
    cm_plot          = cm_plot,
    roc_list         = roc_list,
    roc_plot         = roc_plot,
    pr_list          = pr_list,
    pr_plot          = pr_plot,
    class_auc        = class_auc,
    pr_auc           = pr_auc,
    macro_pr_auc     = macro_pr_auc,
    sensitivity      = sensitivity,
    specificity      = specificity,
    macro_auc        = macro_auc,   
    handtill_obj     = handtill_obj,
    summary          = summary_df,
    classwise        = classwise_df,
    pred_probs       = pred_probs
  )
} 


# Heat maps for prediction probabilities 

prob_heatmaps <- function(prob_df, pred, test_pheno, sample_id_column, heatmap_label,
                          true_test_labels, mol_group, sub_group, MYC_status, MYCN_status, OS, follow_up, Age,
                          continuum = "DN/bi", heatmap_panel = TRUE, correlate_plot = NULL,
                          column_for_panel = "DN", file_prefix = NULL, save_in = "~/Thesis_models/DN_models/") { 
  
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
  
  if (!is.null(file_prefix)) {
    label_shorthand <- file_prefix
  } else if (grepl("Random", heatmap_label) && grepl("Test NMB", heatmap_label)) {
    label_shorthand <- "RF_test_NMB"
  } else if (grepl("Random", heatmap_label) && grepl("Cavalli-SHH-Inf-MB", heatmap_label) && grepl("Age", heatmap_label)) {
    label_shorthand <- "RF_Cavalli_SHH_Inf_sansage"
  } else if (grepl("Random", heatmap_label) && grepl("Cavalli-SHH-Inf-MB", heatmap_label)) {
    label_shorthand <- "RF_Cavalli_SHH_Inf"
  } else if (grepl("Random", heatmap_label) && grepl("Cavalli-SHH-MB", heatmap_label) && grepl("Age", heatmap_label)) {
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
  } else if (grepl("Age", heatmap_label) && grepl("'NMB+Infant'", heatmap_label, fixed = TRUE)) {
    label_shorthand <- "RF_train_sansage" 
  } else {
    stop("Check labels for filenames") 
  }  
  
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
  test_pheno[[MYC_status]]  <- as.factor(test_pheno[[MYC_status]])
  test_pheno[[MYC_status]]  <- droplevels(test_pheno[[MYC_status]])
  test_pheno[[MYCN_status]] <- gsub("0", "Neutral", test_pheno[[MYCN_status]])
  test_pheno[[MYCN_status]] <- gsub("1", "AMP", test_pheno[[MYCN_status]])
  test_pheno[[MYCN_status]] <- as.factor(test_pheno[[MYCN_status]])
  test_pheno[[MYCN_status]] <- droplevels(test_pheno[[MYCN_status]])
  
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
      
      OS  = list(title = "Death", at = c(0, 1), labels = c("No", "Yes")),
      
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
         column_title            = heatmap_label,
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
      title           = element_text(hjust = 0.5, face = "bold") 
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
      title      = plot_title,
      palette    = palette
    ) 
    
    km_plot <- km$plot + coord_cartesian(clip = "off") 
    
    if (theme_to_apply == TRUE) km_plot <- km_plot + theme (axis.title = element_text ( hjust = 0.5,face = "bold"),
                                                             title = element_text(hjust = 0.5, face = "bold"))  
                                                             
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
             title = paste0("P(", column_for_panel, ") ", "in age-ranges")) +
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
             title = paste0("P(", column_for_panel, ") ", "by MYC Amplification Status")) +
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
             title = paste0("P(", column_for_panel, ") ", "by MYCN Amplification Status")) +
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
             title = paste0("P(", column_for_panel, ") ", "by Survival Group")) +
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
             title = paste0("P(", column_for_panel, ") ", "by Principal Molecular Group")) + 
        theme(legend.position = "none",
              axis.text.x     = element_text(angle = 45, hjust = 0.5, size = 6),
              axis.title = element_text ( hjust = 0.5,face = "bold"),
              title = element_text(hjust = 0.5, face = "bold")) + coord_cartesian(clip = "off") +
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
             title = paste0("P(", column_for_panel, ") ", "by Molecular Sub group")) +
        theme(legend.position = "none",
              axis.text.x     = element_text(angle = 45, hjust = 0.5, size = 6),
              axis.title = element_text ( hjust = 0.5,face = "bold"),
              title = element_text(hjust = 0.5, face = "bold")) +
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
      
      pdf(paste0(save_in, label_shorthand, "_.pdf"),
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
                   label.x     = 12,
                   label.y     = 1,
                   size        = 2.5,
                   colour      = "black",
                   inherit.aes = FALSE,
                   data        = prob_df,
                   mapping     = aes(x = Age, y = .data[[column_for_panel]])) +
          geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
          labs(x = "Age (years)", y = paste0("P(", column_for_panel, ")"), colour = "Class", 
               title = paste0("P(", column_for_panel, ") ", " by Chronological Age in ", cohort_name)) +
          theme(axis.title = element_text ( hjust = 0.5,face = "bold"),
                title = element_text(hjust = 0.5, face = "bold")) +
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
               title = paste0("P(", column_for_panel, ") ", "in age-ranges in ", cohort_name)) +
          theme(axis.title = element_text ( hjust = 0.5,face = "bold"),
                title = element_text(hjust = 0.5, face = "bold")) +
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
               title = paste0("P(", column_for_panel, ") ", "by MYC Amplification Status in ", cohort_name)) +
          theme(axis.title = element_text(hjust = 0.5, face = "bold"),
                title      = element_text(hjust = 0.5, face = "bold")) + 
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
               title = paste0("P(", column_for_panel, ") ", "by MYCN Amplification Status in ", cohort_name)) +
          theme(axis.title = element_text (hjust = 0.5,face = "bold"),
                title = element_text(hjust = 0.5, face = "bold")) +
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
               title = paste0("P(", column_for_panel, ") ", "by Principal Molecular Group in ", cohort_name)) +
          theme(legend.position = "none",
                axis.text.x     = element_text(angle = 45, hjust = 0.5, size = 6),
                axis.title = element_text ( hjust = 0.5,face = "bold"),
                title = element_text(hjust = 0.5, face = "bold")) +
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
               title = paste0("P(", column_for_panel, ") ", "by Molecular Sub group in ",cohort_name)) +
          theme(legend.position = "none",
                axis.text.x     = element_text(angle = 45, hjust = 0.5, size = 6),
                axis.title = element_text ( hjust = 0.5,face = "bold"),
                title = element_text(hjust = 0.5, face = "bold")) +
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

# Training binary Random Forest (RF) model on Train NMB for classification 
# of DN and non-DN histology

# ------ Load data ------ #

#Train and Test split of the 'NMB' Cohort is done with caret::createDataPartition 
#in the previous script "Distribution_charts.R". "NOS" histology labels have been removed before the split

train_NMB_pheno<-readRDS("~/NMB/NMB_train_pheno.rds")

if (any(is.na(train_NMB_pheno$Sample_ID))){
  train_NMB_pheno<-train_NMB_pheno[!is.na(train_NMB_pheno$Sample_ID),]
}

# train_NMB_pheno<-train_NMB_pheno %>% mutate(mol_group = case_when( 
#   Age >= 0 & Age < 4.99 ~ "SHH_Inf", 
#   TRUE ~ "SHH_Old"))

test_NMB_pheno<-readRDS("~/NMB/NMB_test_pheno.rds")

# test_NMB_pheno<-test_NMB_pheno %>% mutate(Mol_group = case_when( 
#   Age >= 0 & Age < 4.99 ~ "SHH_Inf", 
#   TRUE ~ "SHH_Old"))


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

if (any(is.na(Cavalli_pheno$Study_ID...12))){
  Cavalli_pheno<-Cavalli_pheno[!is.na(Cavalli_pheno$Study_ID...12),]
}

if (any(is.na(Cavalli_pheno$histology))){
  Cavalli_pheno<-Cavalli_pheno[!is.na(Cavalli_pheno$histology),]
}


Cavalli_betas<-readRDS("~/Cavalli/Cavalli-test/Cavalli_test_betas.rds")

Cavalli_betas<-as.data.frame(t(Cavalli_betas)) 

Cavalli_betas<-Cavalli_betas[match(Cavalli_pheno$Study_ID...12, rownames(Cavalli_betas)),] 

if(all(Cavalli_pheno$Study_ID...12 == rownames(Cavalli_betas))){
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
  train_NMB_pheno %>% filter(CPR_Histology == "CLA", Mol_group == "Group3-MB") %>%
    slice_sample(n = round(DN_counts * counts["CLA", "Group3-MB"] / non_DN_total)),
  
  train_NMB_pheno %>% filter(CPR_Histology == "LCA", Mol_group == "Group4-MB") %>%
    slice_sample(n = round(DN_counts * counts["LCA", "Group4-MB"] / non_DN_total)),
  
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


RF_bi_NMB<-maintrain(betas.. = train_NMB[,!colnames(train_NMB) %in% "bi_histo"],
                     y.. = as.factor(train_NMB$bi_histo),
                     method = "RF",
                     ntrees = 500, #default
                     p = 200, #default
                     topfeaturenumber = NULL,
                     subset.CpGs = NULL,
                     calibrationmethod = c("MR","FLR","LR")
)

saveRDS(RF_bi_NMB, "~/Thesis_models/DN_models/Imb_Train_NMB/RF_imb_NMB.rds") 

# Train an RF model on under-sampled dataset

RF_bi_US_NMB<-methylClass::maintrain(betas.. = as.matrix(US_train_NMB[,!colnames(US_train_NMB) %in% "bi_histo"]),
                                     y.. = as.factor(US_train_NMB$bi_histo),
                                     method = "RF",
                                     ntrees = 500, #default
                                     p = 200, #default
                                     topfeaturenumber = NULL,
                                     subset.CpGs = NULL,
                                     calibrationmethod = c("MR","FLR","LR")
)

saveRDS(RF_bi_US_NMB, "~/Thesis_models/DN_models/Balanced_train_NMB/RF_bi_US_NMB.rds") 


# -------- Testing of the models and Evaluation metrics -------------- #

# Binary model trained on Imbalanced data

# Load the model if not in the env.

RF_bi_NMB<-readRDS("~/Thesis_models/DN_models/Imb_Train_NMB/RF_imb_NMB.rds") 

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

if(all(rownames(Cavalli_betas) == Cavalli_pheno$Study_ID...12)){
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
  sample_id_column  = "Study_ID...12",
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

# Binary model trained on balanced (under sampled) data

# Load the model if not in the env. 

RF_bi_US_NMB<-readRDS("~/Thesis_models/DN_models/Balanced_train_NMB/RF_bi_US_NMB.rds") 

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

if(all(rownames(Cavalli_betas) == Cavalli_pheno$Study_ID...12)){
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
  sample_id_column  = "Study_ID...12",
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


# ----- Training Binary Random Forest Models on SHH-pathway activated MB patient group ------ #

# Combining SHH samples from 'Infant' data with NMB MB-SHH samples to form the train set 'NMB+Infant-SHH'

# --- Load the Data --- #

# Load  NMB if not in the env.

NMB_pheno<-readRDS("~/NMB/NMB_Final_Cohort_26.rds")

NMB_pheno<-NMB_pheno[NMB_pheno$CPR_Histology != "NOS",] 

NMB_pheno<-NMB_pheno[!is.na(NMB_pheno$Sample_ID),] 

NMB_pheno<- NMB_pheno %>% mutate(bi_histo = case_when( 
  CPR_Histology %in% c("CLA","LCA") ~ "non_DN", 
  CPR_Histology %in% c("DN", "MBEN") ~ "DN")) 

# Beta value matrix post pre processing with minfi- 450k + EPIC combined

NMB_betas<-readRDS("~/NMB/Final_NMB_betas.rds")

NMB_betas<-NMB_betas[match(NMB_pheno$Sample_ID, rownames(NMB_betas)),]


# Test sets 

# Cavalli

Cavalli_pheno<- readRDS("~/Cavalli/Final_Cavalli_pheno.rds")

# Matching Cavalli's Age distribution to train_NMB's

Cavalli_pheno<-Cavalli_pheno %>% filter(Age <= max(NMB_pheno$Age))

if (any(is.na(Cavalli_pheno$Study_ID...12))){
  Cavalli_pheno<-Cavalli_pheno[!is.na(Cavalli_pheno$Study_ID...12),]
}

if (any(is.na(Cavalli_pheno$histology))){
  Cavalli_pheno<-Cavalli_pheno[!is.na(Cavalli_pheno$histology),]
}

Cavalli_betas<-readRDS("~/Cavalli/Cavalli-test/Cavalli_test_betas.rds")

Cavalli_betas<-as.data.frame(t(Cavalli_betas)) 

Cavalli_betas<-Cavalli_betas[match(Cavalli_pheno$Study_ID...12, rownames(Cavalli_betas)),] 

if(all(Cavalli_pheno$Study_ID...12 == rownames(Cavalli_betas))){
  Cavalli_pheno<-Cavalli_pheno %>% mutate(bi_histo = case_when( 
    histology  %in% c("CLA","LCA") ~ "non_DN", 
    histology  %in% c("DN", "MBEN") ~ "DN"))
} else {
  message("Mismatch in Sample IDs in Pheno and betas")
  
}

Cavalli_SHH_pheno<-Cavalli_pheno[Cavalli_pheno$mol_grp %in% c("SHH-Inf-MB", "SHH-Old-MB", "SHH-MB"),]

Cavalli_SHH_betas<-Cavalli_betas[match(Cavalli_SHH_pheno$Study_ID...12, rownames(Cavalli_betas)),] 

# SHH-MB

SHH_NMB_pheno<-NMB_pheno[NMB_pheno$Mol_grp %in% c("SHH-Inf-MB", "SHH-Old-MB", "SHH-MB"),] 

SHH_NMB_table<-data.frame("Sample_ID" = SHH_NMB_pheno$Sample_ID,
                          "CPR_Histology" = as.factor(SHH_NMB_pheno$CPR_Histology),
                          "Mol_group" = as.factor(SHH_NMB_pheno$Mol_grp),
                          "Sub_group" = as.factor(SHH_NMB_pheno$sub_group),
                          "Age" = SHH_NMB_pheno$Age,
                          "MYC_Status" = as.factor(SHH_NMB_pheno$`MYC Status`),
                          "MYCN_Status" = as.factor(SHH_NMB_pheno$`MYCN Status`),
                          "Bi_histo" = as.factor(SHH_NMB_pheno$bi_histo),
                          "OS" = as.factor(SHH_NMB_pheno$OS),
                          "Follow_up" = SHH_NMB_pheno$`Follow up years`)

# In-house Infant Cohort (Richardson et.al, 2024)

inf_betas<-readRDS("~/NMB/Final_Inf_betas.rds")

inf_betas<-inf_betas[,!colnames(inf_betas) %in% c("mol.grp","histo")]

inf_pheno<-read.csv("~/NMB/updatedSampleSheet.csv")

# Extracting only CPR-labelled samples 

inf_pheno<-inf_pheno[inf_pheno$Nature_of_pathology_review == "Central",]

inf_pheno<-inf_pheno[inf_pheno$Path_master != "NOS",]

inf_betas<-inf_betas[match(inf_pheno$Sentrix_ID, rownames(inf_betas)),]

if(all(rownames(inf_betas) == inf_pheno$Sentrix_ID)){
  inf_pheno$bi_histo<-as.factor(inf_pheno$Path_master)
  inf_pheno<-inf_pheno %>% mutate ( bi_histo = case_when( bi_histo %in% c("CLA","LCA") ~ "non_DN",
                                                          bi_histo %in% c("DN","MBEN") ~ "DN"))
} else {
  inf_pheno<-inf_pheno[match(rownames(inf_betas), inf_pheno$Sentrix_ID),]
  
  inf_pheno<-inf_pheno[!is.na(inf_pheno$Sentrix_ID),]
  
  inf_betas<-inf_betas[match(inf_pheno$Sentrix_ID, rownames(inf_betas)),] 
  
  if(all(rownames(inf_betas) == inf_pheno$Sentrix_ID)) { 
    
    inf_pheno$bi_histo<-as.factor(inf_pheno$Path_master)
    
    inf_pheno<-inf_pheno %>% mutate ( bi_histo = case_when( bi_histo %in% c("CLA","LCA") ~ "non_DN",
                                                            bi_histo %in% c("DN","MBEN") ~ "DN"))
    
  } else { 
    
    message("Mismatch in inf_betas and inf_pheno Sentrix/Sample_IDs") 
  }
} 

SHH_inf_pheno<-inf_pheno[inf_pheno$Subgroup == "SHH",]

SHH_Inf_table<-data.frame("Sample_ID" = SHH_inf_pheno$Sentrix_ID,
                          "CPR_Histology" = as.factor(SHH_inf_pheno$Path_master),
                          "Sub_group" = as.factor(SHH_inf_pheno$MNP_12.5_call),
                          "Age" = SHH_inf_pheno$Age_at_diagnosis,
                          "MYC_Status" = as.factor(SHH_inf_pheno$MYC_Amplified),
                          "MYCN_Status" = as.factor(SHH_inf_pheno$MYCN_Amplified),
                          "Bi_histo" = as.factor(SHH_inf_pheno$bi_histo),
                          "OS" = as.factor(SHH_inf_pheno$Overall.Survival),
                          "Follow_up" = SHH_inf_pheno$Overall.survival.time)

SHH_Inf_table$Mol_group<-"SHH-Inf-MB"

SHH_inf_betas<-inf_betas[match(SHH_inf_pheno$Sentrix_ID, rownames(inf_betas)),]



# Check for Overlap in SHH-NMB and Infant 
# (Sentrix_ID confirmed as the correct join key for this de-duplication step)

if(any(SHH_inf_pheno$Sentrix_ID %in% SHH_NMB_pheno$Sentrix_ID)){
  SHH_NMB_pheno<-SHH_NMB_pheno[!SHH_NMB_pheno$Sentrix_ID %in% SHH_inf_pheno$Sentrix_ID,]
  
  SHH_NMB_pheno<-SHH_NMB_pheno[!is.na(SHH_NMB_pheno$Sentrix_ID),] 
} 

SHH_NMB_betas<-NMB_betas[match(SHH_NMB_pheno$Sample_ID, rownames(NMB_betas)),] 

if(all(colnames(SHH_NMB_betas) == colnames(SHH_inf_betas))){
  
  Train_SHH<-rbind(SHH_NMB_betas,SHH_inf_betas)
  
} else {
  tryCatch({
    SHH_inf_betas <- SHH_inf_betas[, match(colnames(SHH_NMB_betas),colnames(SHH_inf_betas))]
  }, error = function(e) {
    message("Column matching failed: ", conditionMessage(e))
  })
  Train_SHH<-rbind(SHH_NMB_betas,SHH_inf_betas) 
}


Train_SHH_pheno<-rbind(SHH_NMB_table,SHH_Inf_table)

Train_SHH_pheno<-Train_SHH_pheno[match(rownames(Train_SHH), Train_SHH_pheno$Sample_ID),]

Train_SHH_pheno<-Train_SHH_pheno[!is.na(Train_SHH_pheno$Sample_ID),]

if(all(rownames(Train_SHH) == Train_SHH_pheno$Sample_ID)){
  
  Train_SHH$bi_histo<-as.factor(Train_SHH_pheno$Bi_histo)
  
  Train_SHH<-Train_SHH[!is.na(Train_SHH$bi_histo),]
  
} else {
  
  message("No histology labels in Train-SHH") 
} 

Train_SHH_pheno<-Train_SHH_pheno[match(rownames(Train_SHH),Train_SHH_pheno$Sample_ID),] 

# ------ Training methylClass (Liu, 2024) models with all 441870 input probes for DN Classification within SHH-MB ----- 

RF_comb_SHH<-maintrain(betas.. = as.data.frame(Train_SHH[,!colnames(Train_SHH) %in% "bi_histo"]),
                       y.. = as.factor(Train_SHH$bi_histo),
                       method = "RF",
                       ntrees = 500, #default
                       p = 200, #default
                       topfeaturenumber = NULL,
                       subset.CpGs = NULL,
                       seed = 42,
                       calibrationmethod = c("MR","FLR","LR")
)

saveRDS(RF_comb_SHH, "~/Thesis_models/DN_models/SHH_NMB_Infant/RF_comb_SHH.rds")


# Test the model on Train_SHH (itself)

message("Testing and evaluating on Train-SHH")

RF_comb_SHH<-readRDS("~/Thesis_models/DN_models/SHH_NMB_Infant/RF_comb_SHH.rds")

if(all(rownames(Train_SHH) == Train_SHH_pheno$Sample_ID)){
  tryCatch(train_results_SHH<-test_and_evaluate_model(
    model = RF_comb_SHH,
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "'NMB'+'Infant' SHH-MB",
    test_data = Train_SHH,
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

if(all(rownames(Cavalli_SHH_betas) == Cavalli_SHH_pheno$Study_ID...12)){
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
  sample_id_column  = "Study_ID...12",
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

# ------ Training the RF model on feature set without chronological age related probes 

Train_SHH_mat<-Train_SHH[,!colnames(Train_SHH) %in% "bi_histo"]

if(all(rownames(Train_SHH_mat) == Train_SHH_pheno$Sample_ID)){
  
  library("tibble")
  library("dplyr")
  
  age_cor <- apply(Train_SHH_mat, 2, function(x) {
    ct <- cor.test(x, Train_SHH_pheno$Age, method = "spearman")
    c(rho = ct$estimate, p = ct$p.value)
  })
  
  age_cor_df <- as.data.frame(t(age_cor)) |>
    tibble::rownames_to_column("probe") |>
    dplyr::rename(rho = `rho.rho`) |>
    dplyr::mutate(p_adj = p.adjust(p, method = "BH")) |>
    dplyr::arrange(p_adj)
  
  age_probes <- age_cor_df$probe[age_cor_df$p_adj < 0.05 & abs(age_cor_df$rho) > 0.3] } else {
    
    Train_SHH_pheno<-Train_SHH_pheno[match(rownames(Train_SHH), Train_SHH_pheno$Sample_ID),] 
    
    library("tibble")
    library("dplyr")
    
  
    age_cor <- apply(Train_SHH_mat, 2, function(x) {
      ct <- cor.test(x, Train_SHH_pheno$Age, method = "spearman")
      c(rho = ct$estimate, p = ct$p.value)
    })
    
    age_cor_df <- as.data.frame(t(age_cor)) |>
      tibble::rownames_to_column("probe") |>
      dplyr::rename(rho = `rho.rho`) |>
      dplyr::mutate(p_adj = p.adjust(p, method = "BH")) |>
      dplyr::arrange(p_adj)
    
    age_probes <- age_cor_df$probe[age_cor_df$p_adj < 0.05 & abs(age_cor_df$rho) > 0.3]
    
  }

saveRDS(age_probes, "~/Thesis_models/DN_models/SHH_NMB_Infant/age_probes.rds") 

# Removing chronological age related probes 

Train_SHH_mat<-Train_SHH_mat[,!colnames(Train_SHH_mat) %in% age_probes]

#check 

length(age_probes)

ncol(Train_SHH) #should be 441870 CpG probes + 1 histo col.

ncol(Train_SHH_mat) 

if(all(rownames(Train_SHH_mat) == Train_SHH_pheno$Sample_ID)) {
  Train_SHH_mat$bi_histo<-as.factor(Train_SHH_pheno$Bi_histo)
}

# Model

RF_comb_SHH_sansage<-maintrain(betas.. = Train_SHH_mat[,!colnames(Train_SHH_mat) %in% "bi_histo"],
                               y.. = as.factor(Train_SHH_mat$bi_histo),
                               method = "RF",
                               ntrees = 500, #default
                               p = 200, #default
                               topfeaturenumber = NULL,
                               subset.CpGs = NULL,
                               seed = 42,
                               calibrationmethod = c("MR","FLR","LR")
)

saveRDS(RF_comb_SHH_sansage, "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/RF_comb_SHH_sansage.rds")

# Test the model on Train_SHH (itself)

RF_comb_SHH_sansage<-readRDS("~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/RF_comb_SHH_sansage.rds")

message("Testing and evaluating on Train-SHH")

if(all(rownames(Train_SHH_mat) == Train_SHH_pheno$Sample_ID)){
  tryCatch(train_results_SHH_sa<-test_and_evaluate_model(
    model = RF_comb_SHH_sansage,
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "'NMB+'Infant' SHH-MB with removed Age-related probes",
    test_data = Train_SHH_mat,
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

if(all(rownames(Cavalli_SHH_betas) == Cavalli_SHH_pheno$Study_ID...12)){
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
  sample_id_column  = "Study_ID...12",
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

# ----- Horvath (Horvath, 2013) Methylation Age Estimation and correlation -----

#BioCmanager::install("methylclock")

require(methylclock) 

# In Train-SHH

Train_SHH_probes<-t(Train_SHH_mat[,!colnames(Train_SHH_mat) %in% "bi_histo"])

Train_SHH_clocks<-DNAmAge(Train_SHH_probes) 

Train_SHH_clock_df<-as.data.frame(Train_SHH_clocks)

Train_SHH_clock_df$Sample_ID<-as.factor(Train_SHH_clock_df$id)

Train_SHH_clock<-merge(Train_SHH_pheno, Train_SHH_clock_df, by = "Sample_ID")

# Simple Epigenetic Age Acceleration ( as in previous Medulloblastoma papers)

Train_SHH_clock$horvath_accel  <- Train_SHH_clock$Horvath - Train_SHH_clock$Age

# Residuals from Linear Regression - to control for chronological age confounding 

train_lm_fit <- lm(Horvath ~ Age, data = Train_SHH_clock)

Train_SHH_clock$horvath_accel_res <- resid(train_lm_fit) 

Train_SHH_clock<-Train_SHH_clock[match(Train_SHH_pheno$Sample_ID, Train_SHH_clock$Sample_ID),] 

if(all(Train_SHH_pheno$Sample_ID == Train_SHH_clock$Sample_ID)) {
  Train_SHH_pheno$Horvath_Age   <- Train_SHH_clock$Horvath
  Train_SHH_pheno$Horvath_accel <- Train_SHH_clock$horvath_accel
  Train_SHH_pheno$Horvath_accel_res<-Train_SHH_clock$horvath_accel_res
} else {
  warning("Check mismatch in Sample IDs for Train SHH Epigenetic Clock")
}

# Load saved probability scores

Train_SHH_prob_df<-read.csv("~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/RF_SHH_sansage_Train-SHH-MB_probabilities_.csv")

if(all(Train_SHH_prob_df$Sample_ID == Train_SHH_pheno$Sample_ID)){
  
  Train_SHH_prob_df$True_label<-as.factor(Train_SHH_pheno$CPR_Histology)
  
  Train_SHH_prob_df$Age<-Train_SHH_pheno$Age
  
  Train_SHH_prob_df$Horvath_Age<-Train_SHH_pheno$Horvath_Age
  
  Train_SHH_prob_df$Horvath_accel<-Train_SHH_pheno$Horvath_accel
  
  Train_SHH_prob_df$Horvath_accel_res<-Train_SHH_pheno$Horvath_accel_res
} else {
  
  message ("Check for match in Sample IDs of the loaded Prob_df and Phenotype table")
}  

Train_SHH_prob_df <- Train_SHH_prob_df %>% mutate (Class = case_when (True_label == 'CLA' & pred == "non_DN" ~ "Correct_non_DN",
                                                                      True_label == 'LCA' & pred == "non_DN" ~ "Correct_non_DN",
                                                                      True_label == "DN" & pred == "DN" ~ "Correct_DN",
                                                                      True_label == "MBEN"& pred == "DN" ~ "Correct_DN",
                                                                      True_label == 'DN' & pred == "non_DN" ~ "DN->non_DN",
                                                                      #true_label == 'DN' & pred == "MBEN" ~ "DN->MBEN",
                                                                      True_label == 'MBEN' & pred == "non_DN" ~ "DN->non_DN",
                                                                      #true_label == 'MBEN' & pred == "DN" ~ "MBEN->DN",
                                                                      True_label == 'CLA' & pred == "DN" ~ "non_DN->DN",
                                                                      #true_label == 'CLA' & pred == "MBEN" ~ "non_DN->MBEN",
                                                                      True_label == 'LCA' & pred == "DN" ~ "non_DN->DN"))

p_horvath_train_SHH<- ggplot(Train_SHH_prob_df, aes(x = Horvath_accel, y = DN, colour = Class)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_color_manual(values = c("Correct_DN"    = "darkgreen", "Correct_non_DN" = "purple4",
                                "non_DN->DN"    = "pink3",     "DN->non_DN"     = "orange")) +
  stat_cor(method      = "pearson",
           size        = 2.5,
           colour      = "black",
           label.x = 20,
           label.y = 1,
           inherit.aes = FALSE,
           data        = Train_SHH_prob_df,
           mapping     = aes(x = Horvath_accel, y = DN)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0,   linetype = "dashed", colour = "grey50") +
  labs(x = "Horvath Epigenetic Age acceleration (years)", y = "P(DN)", colour = "Class", 
       title = "P (DN) by Horvath Epigenetic Age Acceleration in 'Infant+NMB' SHH-MB") + 
  theme(legend.position = "none",
         title = element_text(hjust = 0.5, face = "bold"), 
         axis.title = element_text(hjust = 0.5, face = "bold")) +
  coord_cartesian(clip = "off")

ggsave(file.path("~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Train_SHH_Horvath_accel.png"),
       p_horvath_train_SHH, width = 9, height = 6, dpi = 300, bg = "white")


# Simple Epigenetic Age Acceleration in Histology groups 

p_horvath_accel_histo <- ggplot(Train_SHH_prob_df,
                     aes(x = pred, y = Horvath_accel, fill = pred)) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.4, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
  scale_fill_manual(values = c("DN" = "darkgreen", "non_DN" = "purple")) +
  stat_compare_means(method = "kruskal.test",
                     label.y = 1,
                     size    = 2.5) +
  labs(x = "Histology group", y = "Horvath Methylation Age Acceleration (in years)",
       title = "Horvath Epigenetic Age Acceleration in Histology Groups of Train 'NMB+Infant' SHH-MB") +
  theme(legend.position = "none",
        axis.text     = element_text(hjust = 0.5, face = "bold"),
        title = element_text( hjust = 0.5, face = "bold")) +
  coord_cartesian(clip = "off")

ggsave(file.path("~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Train_SHH_Horvath_accel_in_histo.png"),
       p_horvath_accel_histo, width = 9, height = 6, dpi = 300, bg = "white")


# Correlation between Simple Epigenetic Age Acceleration and Chronological Age- with Spearman correlation coefficient 

p_accel_vs_age <- ggplot(Train_SHH_prob_df, aes(x = Age, y = Horvath_accel, colour = pred)) +
         geom_point(size = 1, alpha = 0.7) +
         scale_x_continuous(limits = c(0, 45), breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45)) +
         geom_smooth(method = "lm", color = "#D85A30", se = TRUE) +
         scale_color_manual(values = c("DN" = "darkgreen", "non_DN" = "purple")) +
         stat_cor(method      = "spearman",
                                 label.x     = 35,
                                label.y     = 30,
                                 size        = 2.5,
                                 colour      = "black",
                                inherit.aes = FALSE,
                                 data        = Train_SHH_prob_df,
                                 mapping     = aes(x = Age, y = Horvath_accel)) + 
         geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
         labs(x = "Chronological Age (years)", y = "Horvath Methylation Age Acceleration (in years)", colour = "Histology Group", 
              title = "Horvath Epigenetic Age Acceleration vs Chronological Age in Train 'NMB+Infant' SHH-MB") + 
  theme(axis.title = element_text(hjust = 0.5, face = "bold"),title = element_text(hjust = 0.5, face = "bold")) +
  coord_cartesian(clip = "off")

ggsave(file.path("~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Train_SHH_Horvath_accel_vs_age.png"),
       p_accel_vs_age, width = 9, height = 6, dpi = 300, bg = "white")



# Residual Epigenetic Age Acceleration in Histology groups 

p_horvath_accel_res_histo <- ggplot(Train_SHH_prob_df,
                                aes(x = pred, y = Horvath_accel_res, fill = pred)) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.4, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
  scale_fill_manual(values = c("DN" = "darkgreen", "non_DN" = "purple")) +
  stat_compare_means(method = "kruskal.test",
                     label.y = 1,
                     size    = 2.5) +
  labs(x = "Histology group", y = "Residual Horvath Methylation Age Acceleration (in years)",
       title = " Residual Horvath Epigenetic Age Acceleration in Histology Groups of Train 'NMB+Infant' SHH-MB") +
  theme(legend.position = "none",
        axis.text     = element_text(hjust = 0.5, face = "bold"),
        title = element_text( hjust = 0.5, face = "bold")) +
  coord_cartesian(clip = "off")

ggsave(file.path("~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Train_SHH_Horvath_accel_res_in_histo.png"),
       p_horvath_accel_res_histo, width = 9, height = 6, dpi = 300, bg = "white")

# Correlation between residuals Epigenetic Age Acceleration and Chronological Age- with Spearman correlation coefficient 

p_accel_res_vs_age <- ggplot(Train_SHH_prob_df, aes(x = Age, y = Horvath_accel_res, colour = pred)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_x_continuous(limits = c(0, 45), breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45)) +
  geom_smooth(method = "lm", color = "#D85A30", se = TRUE) +
  scale_color_manual(values = c("DN" = "darkgreen", "non_DN" = "purple")) +
  stat_cor(method      = "spearman",
           label.x     = 35,
           label.y     = 30,
           size        = 2.5,
           colour      = "black",
           inherit.aes = FALSE,
           data        = Train_SHH_prob_df,
           mapping     = aes(x = Age, y = Horvath_accel)) + 
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(x = "Chronological Age (years)", y = " Residual Horvath Methylation Age Acceleration (in years)", colour = "Histology Group", 
       title = "Residual Epigenetic Age Acceleration vs Chronological Age in Train 'NMB+Infant' SHH-MB") + 
  theme(axis.title = element_text(hjust = 0.5, face = "bold"),title = element_text(hjust = 0.5, face = "bold")) +
  coord_cartesian(clip = "off")

ggsave(file.path("~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Train_SHH_Horvath_accel_res_vs_age.png"),
       p_accel_res_vs_age, width = 9, height = 6, dpi = 300, bg = "white")

# In Cavalli

Cavalli_SHH_pheno$Sample_ID<-as.factor(Cavalli_SHH_pheno$Study_ID...12)

Cavalli_probes<-t(Cavalli_SHH_betas)

Cavalli_clocks<-DNAmAge(Cavalli_probes) 

Cavalli_clock_df<-as.data.frame(Cavalli_clocks)

Cavalli_clock_df$Sample_ID<-as.factor(Cavalli_clock_df$id)

Cavalli_clock<-merge(Cavalli_SHH_pheno, Cavalli_clock_df, by = "Sample_ID")


# Simple Epigenetic Age Acceleration 

Cavalli_clock$horvath_accel  <- Cavalli_clock$Horvath - Cavalli_clock$Age

# Residuals from Linear Regression - to control for chronological age confounding 

cavalli_lm_fit <- lm(Horvath ~ Age, data = Cavalli_clock)

Cavalli_clock$horvath_accel_res <- resid(cavalli_lm_fit)

Cavalli_clock<-Cavalli_clock[match(Cavalli_SHH_pheno$Study_ID...12, Cavalli_clock$Sample_ID),] 

if(all(Cavalli_SHH_pheno$Sample_ID == Cavalli_clock$Sample_ID)) {
  Cavalli_SHH_pheno$Horvath_Age   <- Cavalli_clock$Horvath
  Cavalli_SHH_pheno$Horvath_accel <- Cavalli_clock$horvath_accel
  Cavalli_SHH_pheno$Horvath_accel_res <- Cavalli_clock$horvath_accel_res
} else {
  warning("Check mismatch in Sample IDs for Cavalli Epigenetic Clock")
} 

Cavalli_SHH_prob_df<-read.csv("~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/RF_SHH_sansage_Cavalli-SHH-MB_probabilities_.csv")

if(all(Cavalli_SHH_prob_df$Sample_ID == Cavalli_SHH_pheno$Study_ID...12)){
  
  Cavalli_SHH_prob_df$True_label<-as.factor(Cavalli_SHH_pheno$histology)
  
  Cavalli_SHH_prob_df$Age<-Cavalli_SHH_pheno$Age
  
  Cavalli_SHH_prob_df$Horvath_Age<-Cavalli_SHH_pheno$Horvath_Age
  
  Cavalli_SHH_prob_df$Horvath_accel<-Cavalli_SHH_pheno$Horvath_accel
  
  Cavalli_SHH_prob_df$Horvath_accel_res <- Cavalli_SHH_pheno$Horvath_accel_res
  
  Cavalli_SHH_prob_df$bi_histo <- as.factor(Cavalli_SHH_pheno$bi_histo)
  
} else {
  
  message ("Check for match in Sample IDs of the loaded Prob_df and Phenotype table")
}  

Cavalli_SHH_prob_df <- Cavalli_SHH_prob_df %>% mutate (Class = case_when (True_label == 'CLA' & pred == "non_DN" ~ "Correct_non_DN",
                                                                          True_label == 'LCA' & pred == "non_DN" ~ "Correct_non_DN",
                                                                          True_label == "DN" & pred == "DN" ~ "Correct_DN",
                                                                          True_label == "MBEN"& pred == "DN" ~ "Correct_DN",
                                                                          True_label == 'DN' & pred == "non_DN" ~ "DN->non_DN",
                                                                          #true_label == 'DN' & pred == "MBEN" ~ "DN->MBEN",
                                                                          True_label == 'MBEN' & pred == "non_DN" ~ "DN->non_DN",
                                                                          #true_label == 'MBEN' & pred == "DN" ~ "MBEN->DN",
                                                                          True_label == 'CLA' & pred == "DN" ~ "non_DN->DN",
                                                                          #true_label == 'CLA' & pred == "MBEN" ~ "non_DN->MBEN",
                                                                          True_label == 'LCA' & pred == "DN" ~ "non_DN->DN"))

p_horvath_Cavalli<- ggplot(Cavalli_SHH_prob_df, aes(x = Horvath_accel, y = DN, colour = Class)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_color_manual(values = c("Correct_DN"    = "darkgreen", "Correct_non_DN" = "purple4",
                                "non_DN->DN"    = "pink3",     "DN->non_DN"     = "orange")) +
  stat_cor(method      = "pearson",
           size        = 2.5,
           colour      = "black",
           label.x = 20,
           label.y = 1,
           inherit.aes = FALSE,
           data        = Cavalli_SHH_prob_df,
           mapping     = aes(x = Horvath_accel, y = DN)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0,   linetype = "dashed", colour = "grey50") +
  labs(x = "Horvath Epigenetic Age acceleration (years)", y = "P(DN)", colour = "Class", 
       title = "P (DN) by Horvath Epigenetic Age Acceleration in 'Infant+NMB' SHH-MB")+ 
  theme( legend.position = "none",
   title = element_text(hjust = 0.5, face = "bold"), 
    axis.title = element_text(hjust = 0.5, face = "bold")) +
  coord_cartesian(clip = "off")

ggsave(file.path("~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Cavalli_Horvath_accel.png"),
       p_horvath_Cavalli, width = 9, height = 6, dpi = 300, bg = "white")


# Simple Epigenetic Age Acceleration in Histology groups 

p_horvath_accel_histo_cavalli <- ggplot(Cavalli_SHH_prob_df,
                                aes(x = pred, y = Horvath_accel, fill = bi_histo)) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.4, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
  scale_fill_manual(values = c("DN" = "darkgreen", "non_DN" = "purple")) +
  stat_compare_means(method = "kruskal.test",
                     label.y = 1,
                     size    = 2.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  scale_y_continuous(limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "Histology group", y = "Horvath Methylation Age Acceleration (in years)",
       title = "Horvath Epigenetic Age Acceleration in Histology Groups of Cavalli SHH-MB") +
  theme(legend.position = "none",
        axis.text     = element_text(hjust = 0.5, face = "bold"),
        title = element_text( hjust = 0.5, face = "bold")) +
  coord_cartesian(clip = "off")

ggsave(file.path("~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Cavalli_Horvath_accel_in_histo.png"),
       p_horvath_accel_histo_cavalli, width = 9, height = 6, dpi = 300, bg = "white")


# Correlation between Simple Epigenetic Age Acceleration and Chronological Age- with Spearman correlation coefficient 

p_accel_vs_age_cavalli <- ggplot(Cavalli_SHH_prob_df, aes(x = Age, y = Horvath_accel, colour = bi_histo)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_x_continuous(limits = c(0, 45), breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45)) +
  geom_smooth(method = "lm", color = "#D85A30", se = TRUE) +
  scale_color_manual(values = c("DN" = "darkgreen", "non_DN" = "purple")) +
  stat_cor(method      = "spearman",
           label.x     = 35,
           label.y     = 30,
           size        = 2.5,
           colour      = "black",
           inherit.aes = FALSE,
           data        = Cavalli_SHH_prob_df,
           mapping     = aes(x = Age, y = Horvath_accel)) + 
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(x = "Chronological Age (years)", y = "Horvath Methylation Age Acceleration (in years)", colour = "Histology Group", 
       title = "Horvath Epigenetic Age Acceleration vs Chronological Age in Cavalli SHH-MB") + 
  theme(axis.title = element_text(hjust = 0.5, face = "bold"),
        title = element_text(hjust = 0.5, face = "bold")) +
  coord_cartesian(clip = "off")

ggsave(file.path("~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Cavalli_Horvath_accel_vs_age.png"),
       p_accel_vs_age_cavalli, width = 9, height = 6, dpi = 300, bg = "white")



# Residual Epigenetic Age Acceleration in Histology groups 

p_horvath_accel_res_histo_cavalli <- ggplot(Cavalli_SHH_prob_df,
                                    aes(x = pred, y = Horvath_accel_res, fill = bi_histo)) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.4, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
  scale_fill_manual(values = c("DN" = "darkgreen", "non_DN" = "purple")) +
  stat_compare_means(method = "kruskal.test",
                     label.y = 1,
                     size    = 2.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  scale_y_continuous(limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "Histology group", y = "Residual Horvath Methylation Age Acceleration (in years)",
       title = " Residual Horvath Epigenetic Age Acceleration in Histology Groups of Cavalli SHH-MB") +
  theme(legend.position = "none",
        axis.text     = element_text(hjust = 0.5, face = "bold"),
        title = element_text( hjust = 0.5, face = "bold")) +
  coord_cartesian(clip = "off")

ggsave(file.path("~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Cavalli_Horvath_accel_res_in_histo.png"),
       p_horvath_accel_res_histo_cavalli, width = 9, height = 6, dpi = 300, bg = "white")

# Correlation between residuals Epigenetic Age Acceleration and Chronological Age- with Spearman correlation coefficient 

p_accel_res_vs_age_cavalli <- ggplot(Cavalli_SHH_prob_df, aes(x = Age, y = Horvath_accel_res, colour = bi_histo)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_x_continuous(limits = c(0, 45), breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45)) +
  geom_smooth(method = "lm", color = "#D85A30", se = TRUE) +
  scale_color_manual(values = c("DN" = "darkgreen", "non_DN" = "purple")) +
  stat_cor(method      = "spearman",
           label.x     = 35,
           label.y     = 30,
           size        = 2.5,
           colour      = "black",
           inherit.aes = FALSE,
           data        = Cavalli_SHH_prob_df,
           mapping     = aes(x = Age, y = Horvath_accel)) + 
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(x = "Chronological Age (years)", y = " Residual Horvath Methylation Age Acceleration (in years)", colour = "Histology Group", 
       title = "Residual Epigenetic Age Acceleration vs Chronological Age in Cavalli SHH-MB") + 
  theme(axis.title = element_text(hjust = 0.5, face = "bold"),title = element_text(hjust = 0.5, face = "bold")) +
  coord_cartesian(clip = "off")

ggsave(file.path("~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Cavalli_Horvath_accel_res_vs_age.png"),
       p_accel_res_vs_age_cavalli, width = 9, height = 6, dpi = 300, bg = "white")



# ----- Training RF model on SHH-Infant (Ages 0- 4.99) ---- 

Train_SHH_inf_pheno<-Train_SHH_pheno[Train_SHH_pheno$Mol_group == "SHH-Inf-MB",] 

summary(Train_SHH_inf_pheno$Age) # Check the ages

Train_SHH_inf<-Train_SHH_mat[match(Train_SHH_inf_pheno$Sample_ID, rownames(Train_SHH_mat)),] #without age-probes 

Train_SHH_inf<-Train_SHH_inf[!is.na(rownames(Train_SHH_inf)),] 

# Cavalli SHH-Infant test

Cavalli_SHH_inf_pheno<-Cavalli_SHH_pheno[Cavalli_SHH_pheno$mol_grp == "SHH-Inf-MB",]

summary(Cavalli_SHH_inf_pheno$Age) # Check the ages

Cavalli_SHH_inf_betas<-Cavalli_SHH_betas[match(Cavalli_SHH_inf_pheno$Sample_ID, rownames(Cavalli_SHH_betas)),]

Cavalli_SHH_inf_betas<-Cavalli_SHH_inf_betas[!is.na(rownames(Cavalli_SHH_inf_betas)),]

# Model

if(all(rownames(Train_SHH_inf) == Train_SHH_inf_pheno$Sample_ID)){
  
  Train_SHH_inf$bi_histo<-as.factor(Train_SHH_inf_pheno$Bi_histo)
} 


RF_comb_SHH_inf<-maintrain(betas.. = as.data.frame(Train_SHH_inf[,!colnames(Train_SHH_inf) %in% "bi_histo"]),
                           y.. = as.factor(Train_SHH_inf$bi_histo),
                           method = "RF",
                           ntrees = 500, #default
                           p = 200, #default
                           topfeaturenumber = NULL,
                           seed = 42,
                           subset.CpGs = NULL,
                           calibrationmethod = c("MR","FLR","LR")
)

saveRDS(RF_comb_SHH_inf, "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/RF_comb_SHH_inf.rds") 


# Load the model if not in the env.

RF_Inf_NMB<-readRDS("~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/RF_comb_SHH_inf.rds") 

message("Testing and evaluating on Train SHH (including all Mol. grps)")

if (all(rownames(Train_SHH_mat) == Train_SHH_pheno$Sample_ID)) {
  Train_SHH_inf_results <- tryCatch(
    test_and_evaluate_model(
      model = RF_Inf_NMB,
      method = "methylClass",
      calibrated = TRUE,
      calibration_method = "LR",
      train_set_name = "'NMB'+'Infant' SHH-Inf-MB with removed Age-related probes",
      test_data = Train_SHH_mat,
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

message("Testing and evaluating on Cavalli-SHH")

if(all(rownames(Cavalli_betas) == Cavalli_pheno$Study_ID...12)){
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
  sample_id_column  = "Study_ID...12",
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



message("Testing and evaluating RF_SHH_MB on Cavalli-SHH-Inf")

# Model Trained on SHH from all age groups 

RF_comb_SHH_sansage<-readRDS("~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/RF_comb_SHH_sansage.rds")


if(all(rownames(Cavalli_betas) == Cavalli_pheno$Study_ID...12)){
  tryCatch(inf_Cavalli_SHH_sa_results <-test_and_evaluate_model(
    model = RF_comb_SHH_sansage,
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "'NMB'+'Infant' SHH-MB with removed Age-related probes",
    test_data = Cavalli_SHH_inf_betas,
    test_data_name = "Cavalli-SHH-Inf-MB",
    true_test_labels = as.factor(Cavalli_SHH_inf_pheno$bi_histo),
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
  test_pheno        = Cavalli_SHH_inf_pheno,
  sample_id_column  = "Study_ID...12",
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



message("Testing and evaluating RF_inf_MB on Cavalli-SHH-Inf")

if(all(rownames(Cavalli_betas) == Cavalli_pheno$Study_ID...12)){
  tryCatch(Cavalli_inf_results <-test_and_evaluate_model(
    model = RF_Inf_NMB,
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "'NMB'+'Infant' SHH-Inf-MB with removed Age-related probes",
    test_data = Cavalli_SHH_inf_betas,
    test_data_name = "Cavalli-SHH-Inf-MB",
    true_test_labels = as.factor(Cavalli_SHH_inf_pheno$bi_histo),
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

Cavalli_inf_heatmap <- prob_heatmaps(
  prob_df           = Cavalli_inf_results$pred_probs,
  pred              = Cavalli_inf_results$pred_probs$pred,
  test_pheno        = Cavalli_SHH_inf_pheno,
  sample_id_column  = "Study_ID...12",
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
  file_prefix       = "RF_SHHInf_sansage_vs_CavalliSHHInf",
  save_in           = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/"
) 


# ------- Testing with the models trained on all of the infant samples which include Local Path calls (with more samples) ------ 


# Train datasets 

NMB_pheno<-readRDS("~/NMB/NMB_Final_Cohort_26.rds")

NMB_pheno<-NMB_pheno[NMB_pheno$CPR_Histology != "NOS",] 

NMB_pheno<-NMB_pheno[!is.na(NMB_pheno$Sample_ID),] 

NMB_pheno<- NMB_pheno %>% mutate(bi_histo = case_when( 
  CPR_Histology %in% c("CLA","LCA") ~ "non_DN", 
  CPR_Histology %in% c("DN", "MBEN") ~ "DN")) 

# Beta value matrix post pre processing with minfi- 450k + EPIC combined

NMB_betas<-readRDS("~/NMB/Final_NMB_betas.rds")

NMB_betas<-NMB_betas[match(NMB_pheno$Sample_ID, rownames(NMB_betas)),]

# SHH-MB

SHH_NMB_pheno<-NMB_pheno[NMB_pheno$Mol_grp %in% c("SHH-Inf-MB", "SHH-Old-MB", "SHH-MB"),] 

SHH_NMB_table<-data.frame("Sample_ID" = SHH_NMB_pheno$Sample_ID,
                          "CPR_Histology" = as.factor(SHH_NMB_pheno$CPR_Histology),
                          "Mol_group" = as.factor(SHH_NMB_pheno$Mol_grp),
                          "Sub_group" = as.factor(SHH_NMB_pheno$sub_group),
                          "Age" = SHH_NMB_pheno$Age,
                          "MYC_Status" = as.factor(SHH_NMB_pheno$`MYC Status`),
                          "MYCN_Status" = as.factor(SHH_NMB_pheno$`MYCN Status`),
                          "Bi_histo" = as.factor(SHH_NMB_pheno$bi_histo),
                          "OS" = as.factor(SHH_NMB_pheno$OS),
                          "Follow_up" = SHH_NMB_pheno$`Follow up years`)

# In-house Infant Cohort (Richardson et.al, 2024)

inf_betas<-readRDS("~/NMB/Final_Inf_betas.rds")

inf_betas<-inf_betas[,!colnames(inf_betas) %in% c("mol.grp","histo")]

inf_pheno<-read.csv("~/NMB/updatedSampleSheet.csv")

inf_pheno<-inf_pheno[inf_pheno$Path_master != "NOS",]

inf_betas<-inf_betas[match(inf_pheno$Sentrix_ID, rownames(inf_betas)),]

if(all(rownames(inf_betas) == inf_pheno$Sentrix_ID)){
  inf_pheno$bi_histo<-as.factor(inf_pheno$Path_master)
  inf_pheno<-inf_pheno %>% mutate ( bi_histo = case_when( bi_histo %in% c("CLA","LCA") ~ "non_DN",
                                                          bi_histo %in% c("DN","MBEN") ~ "DN"))
} else {
  inf_pheno<-inf_pheno[match(rownames(inf_betas), inf_pheno$Sentrix_ID),]
  
  inf_pheno<-inf_pheno[!is.na(inf_pheno$Sentrix_ID),]
  
  inf_betas<-inf_betas[match(inf_pheno$Sentrix_ID, rownames(inf_betas)),] 
  
  if(all(rownames(inf_betas) == inf_pheno$Sentrix_ID)) { 
    
    inf_pheno$bi_histo<-as.factor(inf_pheno$Path_master)
    
    inf_pheno<-inf_pheno %>% mutate ( bi_histo = case_when( bi_histo %in% c("CLA","LCA") ~ "non_DN",
                                                            bi_histo %in% c("DN","MBEN") ~ "DN"))
    
  } else { 
    
    message("Mismatch in inf_betas and inf_pheno Sentrix/Sample_IDs") 
  }
} 

SHH_inf_pheno<-inf_pheno[inf_pheno$Subgroup == "SHH",]

SHH_Inf_table<-data.frame("Sample_ID" = SHH_inf_pheno$Sentrix_ID,
                          "CPR_Histology" = as.factor(SHH_inf_pheno$Path_master),
                          "Sub_group" = as.factor(SHH_inf_pheno$MNP_12.5_call),
                          "Age" = SHH_inf_pheno$Age_at_diagnosis,
                          "MYC_Status" = as.factor(SHH_inf_pheno$MYC_Amplified),
                          "MYCN_Status" = as.factor(SHH_inf_pheno$MYCN_Amplified),
                          "Bi_histo" = as.factor(SHH_inf_pheno$bi_histo),
                          "OS" = as.factor(SHH_inf_pheno$Overall.Survival),
                          "Follow_up" = SHH_inf_pheno$Overall.survival.time)

SHH_Inf_table$Mol_group<-"SHH-Inf-MB"

SHH_inf_betas<-inf_betas[match(SHH_inf_pheno$Sentrix_ID, rownames(inf_betas)),]


# Check for Overlap in SHH-NMB and Infant 

if(any(SHH_inf_pheno$Sentrix_ID %in% SHH_NMB_pheno$Sentrix_ID)){
  SHH_NMB_pheno<-SHH_NMB_pheno[!SHH_NMB_pheno$Sentrix_ID %in% SHH_inf_pheno$Sentrix_ID,]
  
  SHH_NMB_pheno<-SHH_NMB_pheno[!is.na(SHH_NMB_pheno$Sentrix_ID),] 
} 

SHH_NMB_betas<-NMB_betas[match(SHH_NMB_pheno$Sample_ID, rownames(NMB_betas)),] 

if(all(colnames(SHH_NMB_betas) == colnames(SHH_inf_betas))){
  
  Train_SHH<-rbind(SHH_NMB_betas,SHH_inf_betas)
  
} else {
  tryCatch({
    SHH_inf_betas <- SHH_inf_betas[, match(colnames(SHH_NMB_betas),colnames(SHH_inf_betas))]
  }, error = function(e) {
    message("Column matching failed: ", conditionMessage(e))
  })
  Train_SHH<-rbind(SHH_NMB_betas,SHH_inf_betas) 
}


Train_SHH_pheno<-rbind(SHH_NMB_table,SHH_Inf_table)

Train_SHH_pheno<-Train_SHH_pheno[match(rownames(Train_SHH), Train_SHH_pheno$Sample_ID),]

Train_SHH_pheno<-Train_SHH_pheno[!is.na(Train_SHH_pheno$Sample_ID),]

if(all(rownames(Train_SHH) == Train_SHH_pheno$Sample_ID)){
  
  Train_SHH$bi_histo<-as.factor(Train_SHH_pheno$Bi_histo)
  
  Train_SHH<-Train_SHH[!is.na(Train_SHH$bi_histo),]
  
} else {
  
  message("No histology labels in Train-SHH") 
} 

Train_SHH_pheno<-Train_SHH_pheno[match(rownames(Train_SHH),Train_SHH_pheno$Sample_ID),] 


# Removing age-correlated probes

Train_SHH_mat<-Train_SHH[,!colnames(Train_SHH) %in% "bi_histo"]

age_probes<-readRDS("~/Thesis_models/DN_models/SHH_NMB_Infant/age_probes.rds")

Train_SHH_mat<-Train_SHH[,!colnames(Train_SHH_mat) %in% age_probes] 

if(all(rownames(Train_SHH_mat) == Train_SHH_pheno$Sample_ID)){
  
  RF_comb_SHH_sansage_ext<-maintrain(betas.. = Train_SHH_mat[,!colnames(Train_SHH_mat) %in% "bi_histo"],
                                     y.. = as.factor(Train_SHH_pheno$Bi_histo),
                                     method = "RF",
                                     ntrees = 500, #default
                                     p = 200, #default
                                     topfeaturenumber = NULL,
                                     subset.CpGs = NULL,
                                     seed = 42,
                                     calibrationmethod = c("MR","FLR","LR")
  ) } else {
    
    
    Train_SHH_mat<-Train_SHH_mat[match(Train_SHH_pheno$Sample_ID, rownames(Train_SHH_mat)),]
    
    RF_comb_SHH_sansage_ext<-maintrain(betas.. = Train_SHH_mat[,!colnames(Train_SHH_mat) %in% "bi_histo"],
                                       y.. = as.factor(Train_SHH_pheno$Bi_histo),
                                       method = "RF",
                                       ntrees = 500, #default
                                       p = 200, #default
                                       topfeaturenumber = NULL,
                                       subset.CpGs = NULL,
                                       seed = 42,
                                       calibrationmethod = c("MR","FLR","LR")
    )
    
  }


saveRDS(RF_comb_SHH_sansage_ext, "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Extended/RF_comb_SHH_sansage_ext.rds")


# Testing the model

# On Train-SHH (itself)

RF_comb_SHH_sansage_ext<-readRDS("~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Extended/RF_comb_SHH_sansage_ext.rds")

message("Testing and evaluating on Train-SHH")

if(all(rownames(Train_SHH_mat) == Train_SHH_pheno$Sample_ID)){
  tryCatch(train_results_SHH_ext<-test_and_evaluate_model(
    model = RF_comb_SHH_sansage_ext,
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "'NMB'+'Infant'SHH-MB with removed Age-related probes,extended",
    test_data = Train_SHH_mat,
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


train_SHH_heatmap_ext <- prob_heatmaps(
  prob_df           = train_results_SHH_ext$pred_probs,
  pred              = train_results_SHH_ext$pred_probs$pred,
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
  heatmap_label     = "Binary Random Forest Model for classification of DN histopathology,\ntrained on Extended SHH-MB in 'NMB' and 'Infant' with removed age-related probes,and tested on 'NMB+Infant' SHH-MB",
  heatmap_panel     = FALSE,
  save_in           = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Extended/"
)

# Testing on Cavalli

message("Testing and evaluating on Cavalli")

if(all(rownames(Cavalli_SHH_betas) == Cavalli_SHH_pheno$Study_ID...12)){
  tryCatch(Cavalli_results_SHH_ext<-test_and_evaluate_model(
    model = RF_comb_SHH_sansage_ext,
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "'NMB'+'Infant' SHH-MB with removed Age-related probes,extended",
    test_data = Cavalli_SHH_betas,
    test_data_name = "Cavalli-SHH-MB",
    true_test_labels = as.factor(Cavalli_SHH_pheno$bi_histo),
    model_name = "Binary Random Forest",
    save_in = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Extended/"
  ), error = function(e){
    message(sprintf("FAILED [%s / %s]: %s",
                    "Binary Random Forest trained on SHH-MB", "Cavalli", conditionMessage(e)))
    NULL
  } )
  
} else {
  
  stop("Mismatch in Cavalli SHH test set and test pheno Sample IDs")
}


Cavalli_SHH_heatmap_ext <- prob_heatmaps(
  prob_df           = Cavalli_results_SHH_ext$pred_probs,
  pred              = Cavalli_results_SHH_ext$pred_probs$pred,
  test_pheno        = Cavalli_SHH_pheno,
  sample_id_column  = "Study_ID...12",
  true_test_labels  = "histology",
  mol_group         = "mol_grp",
  sub_group         = "sub_grp",
  MYC_status        = "MYC_Status",
  MYCN_status       = "MYCN_Status",
  OS                = "Dead",
  follow_up         = "OS (years)",
  Age               = "Age",
  continuum         = "DN/bi",
  heatmap_label     = "Binary Random Forest Model for classification of DN histopathology,\ntrained on Extended SHH-MB in 'NMB' and 'Infant' with removed age-related probes,and tested on Cavalli-SHH-MB",
  heatmap_panel     = FALSE,
  save_in           = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Extended/"
) 


# On Cavalli SHH-Inf-MB

message("Testing and evaluating on Cavalli")

if(all(rownames(Cavalli_SHH_betas) == Cavalli_SHH_pheno$Study_ID...12)){
  tryCatch(SHH_Cavalli_inf_results_ext<-test_and_evaluate_model(
    model = RF_comb_SHH_sansage_ext, 
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "'NMB'+'Infant' SHH-MB with removed Age-related probes,extended",
    test_data = Cavalli_SHH_inf_betas,
    test_data_name = "Cavalli-SHH-Inf-MB",
    true_test_labels = as.factor(Cavalli_SHH_inf_pheno$bi_histo),
    model_name = "Binary Random Forest",
    save_in = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Extended/"
  ), error = function(e){
    message(sprintf("FAILED [%s / %s]: %s",
                    "Binary Random Forest trained on SHH-MB-ext", "Cavalli_inf", conditionMessage(e)))
    NULL
  } )
  
} else {
  
  stop("Mismatch in Cavalli SHH test set and test pheno Sample IDs")
}


SHH_Cavalli_inf_heatmap_ext <- prob_heatmaps(
  prob_df           = SHH_Cavalli_inf_results_ext$pred_probs,
  pred              = SHH_Cavalli_inf_results_ext$pred_probs$pred,
  test_pheno        = Cavalli_SHH_inf_pheno,
  sample_id_column  = "Study_ID...12",
  true_test_labels  = "histology",
  mol_group         = "mol_grp",
  sub_group         = "sub_grp",
  MYC_status        = "MYC_Status",
  MYCN_status       = "MYCN_Status",
  OS                = "Dead",
  follow_up         = "OS (years)",
  Age               = "Age",
  continuum         = "DN/bi",
  heatmap_label     = "Binary Random Forest Model for classification of DN histopathology,\ntrained on Extended SHH-MB in 'NMB' and 'Infant' with removed age-related probes,and tested on Cavalli-SHH-Inf-MB",
  heatmap_panel     = FALSE,
  save_in           = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/Extended/"
) 


# Infant models

Train_SHH_inf_pheno<-Train_SHH_pheno[Train_SHH_pheno$Mol_group == "SHH-Inf-MB",] 

summary(Train_SHH_inf_pheno$Age) # Check the ages

Train_SHH_inf<-Train_SHH_mat[match(Train_SHH_inf_pheno$Sample_ID, rownames(Train_SHH_mat)),] #without age-probes 

Train_SHH_inf<-Train_SHH_inf[!is.na(rownames(Train_SHH_inf)),] 

if(all(rownames(Train_SHH_inf) == Train_SHH_inf_pheno$Sample_ID)){
  
  RF_comb_Inf_sansage_ext<-maintrain(betas.. = Train_SHH_inf[,!colnames(Train_SHH_inf) %in% "bi_histo"],
                                     y.. = as.factor(Train_SHH_inf_pheno$Bi_histo),
                                     method = "RF",
                                     ntrees = 500, #default
                                     p = 200, #default
                                     topfeaturenumber = NULL,
                                     subset.CpGs = NULL,
                                     seed = 42,
                                     calibrationmethod = c("MR","FLR","LR")
  ) } else {
    
    
    Train_SHH_inf<-Train_SHH_inf[match(Train_SHH_inf_pheno$Sample_ID, rownames(Train_SHH_inf)),]
    
    RF_comb_Inf_sansage_ext<-maintrain(betas.. = Train_SHH_inf[,!colnames(Train_SHH_inf) %in% "bi_histo"],
                                       y.. = as.factor(Train_SHH_inf_pheno$Bi_histo),
                                       method = "RF",
                                       ntrees = 500, #default
                                       p = 200, #default
                                       topfeaturenumber = NULL,
                                       subset.CpGs = NULL,
                                       seed = 42,
                                       calibrationmethod = c("MR","FLR","LR")
    )
    
  }


saveRDS(RF_comb_Inf_sansage_ext, "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/Extended/RF_comb_Inf_sansage_ext.rds")

# Testing the model

# On Train-SHH (itself)

RF_comb_Inf_sansage_ext<-readRDS("~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/Extended/RF_comb_Inf_sansage_ext.rds")

message("Testing and evaluating on Train-SHH")

if(all(rownames(Train_SHH_mat) == Train_SHH_pheno$Sample_ID)){
  tryCatch(train_results_inf_ext<-test_and_evaluate_model(
    model = RF_comb_Inf_sansage_ext,
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "'NMB'+'Infant'SHH-Inf-MB with removed Age-related probes,extended",
    test_data = Train_SHH_mat,
    test_data_name = "Train-SHH-MB",
    true_test_labels = as.factor(Train_SHH_pheno$Bi_histo),
    model_name = "Binary Random Forest",
    save_in = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/Extended/"
  ), error = function(e){
    message(sprintf("FAILED [%s / %s]: %s",
                    "Binary Random Forest trained on SHH-Inf-MB-ext", "Train-SHH", conditionMessage(e)))
    NULL
  } )
  
} else {
  
  stop("Mismatch in Train SHH betas and pheno Sample IDs")
} 


train_inf_heatmap_ext <- prob_heatmaps(
  prob_df           = train_results_SHH_ext$pred_probs,
  pred              = train_results_SHH_ext$pred_probs$pred,
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
  heatmap_label     = "Binary Random Forest Model for classification of DN histopathology,\ntrained on Extended SHH-Inf-MB in 'NMB' and 'Infant' with removed age-related probes,and tested on 'NMB+Infant' SHH-MB",
  heatmap_panel     = FALSE,
  save_in           = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/Extended/"
)

# Testing on Cavalli

message("Testing and evaluating on Cavalli")

if(all(rownames(Cavalli_SHH_betas) == Cavalli_SHH_pheno$Study_ID...12)){
  tryCatch(inf_Cavalli_results_SHH_ext<-test_and_evaluate_model(
    model = RF_comb_Inf_sansage_ext,
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "'NMB'+'Infant' SHH-MB with removed Age-related probes,extended",
    test_data = Cavalli_SHH_betas,
    test_data_name = "Cavalli-SHH-MB",
    true_test_labels = as.factor(Cavalli_SHH_pheno$bi_histo),
    model_name = "Binary Random Forest",
    save_in = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/Extended/"
  ), error = function(e){
    message(sprintf("FAILED [%s / %s]: %s",
                    "Binary Random Forest trained on SHH-Inf-ext", "Cavalli", conditionMessage(e)))
    NULL
  } )
  
} else {
  
  stop("Mismatch in Cavalli SHH test set and test pheno Sample IDs")
}


inf_Cavalli_SHH_heatmap_ext <- prob_heatmaps(
  prob_df           = inf_Cavalli_results_SHH_ext$pred_probs,
  pred              = inf_Cavalli_results_SHH_ext$pred_probs$pred,
  test_pheno        = Cavalli_SHH_pheno,
  sample_id_column  = "Study_ID...12",
  true_test_labels  = "histology",
  mol_group         = "mol_grp",
  sub_group         = "sub_grp",
  MYC_status        = "MYC_Status",
  MYCN_status       = "MYCN_Status",
  OS                = "Dead",
  follow_up         = "OS (years)",
  Age               = "Age",
  continuum         = "DN/bi",
  heatmap_label     = "Binary Random Forest Model for classification of DN histopathology,\ntrained on Extended SHH-Inf-MB in 'NMB' and 'Infant' with removed age-related probes,and tested on Cavalli-SHH-Inf-MB",
  heatmap_panel     = FALSE,
  save_in           = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/Extended/"
) 

# On Cavalli SHH-Inf-MB

message("Testing and evaluating on Cavalli")

if(all(rownames(Cavalli_SHH_betas) == Cavalli_SHH_pheno$Study_ID...12)){
  tryCatch(Cavalli_inf_results_ext<-test_and_evaluate_model(
    model = RF_comb_Inf_sansage_ext, 
    method = "methylClass",
    calibrated = TRUE,
    calibration_method = "LR",
    train_set_name = "'NMB'+'Infant' SHH-MB with removed Age-related probes,extended",
    test_data = Cavalli_SHH_inf_betas,
    test_data_name = "Cavalli-SHH-Inf-MB",
    true_test_labels = as.factor(Cavalli_SHH_inf_pheno$bi_histo),
    model_name = "Binary Random Forest",
    save_in = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/Extended/"
  ), error = function(e){
    message(sprintf("FAILED [%s / %s]: %s",
                    "Binary Random Forest trained on SHH-Inf-ext", "Cavalli_Inf", conditionMessage(e)))
    NULL
  } )
  
} else {
  
  stop("Mismatch in Cavalli SHH test set and test pheno Sample IDs")
}


Cavalli_inf_heatmap_ext <- prob_heatmaps(
  prob_df           = Cavalli_inf_results_ext$pred_probs,
  pred              = Cavalli_inf_results_ext$pred_probs$pred,
  test_pheno        = Cavalli_SHH_inf_pheno,
  sample_id_column  = "Study_ID...12",
  true_test_labels  = "histology",
  mol_group         = "mol_grp",
  sub_group         = "sub_grp",
  MYC_status        = "MYC_Status",
  MYCN_status       = "MYCN_Status",
  OS                = "Dead",
  follow_up         = "OS (years)",
  Age               = "Age",
  continuum         = "DN/bi",
  heatmap_label     = "Binary Random Forest Model for classification of DN histopathology,\ntrained on Extended SHH-Inf-MB in 'NMB' and 'Infant' with removed age-related probes,and tested on Cavalli-SHH-Inf-MB",
  heatmap_panel     = FALSE,
  save_in           = "~/Thesis_models/DN_models/SHH_NMB_Infant/Sans_age/SHH_Inf/Extended/"
) 

# Warning message:Computation failed in `stat_compare_means()`.
# Caused by error in `kruskal.test.default()`: 
# ! all observations are in the same group as samples only from SHH-Inf-MB are applied for mol. group comparison of P(DN)




