#================================================================================================================================#
# ----------------------- RNAseq ML models------------------------------------------
#===============================================================================================================================#

# This script includes loading of samplesheets and RNA count data, including pre split batch-corrected, pre-processed train and test sets and selected feature sets
# Training RF, SVM , XGB models and testing them with custom functions (Custom_functions.R)
# Splitting CPR NMB Cohort to train and test models on CPR cohorts
# Integration of RNAseq count data with methylation beta values with SNF and training ML models on network features
# 
#==================================================================================================================#

# ---- Load packages -----

require(matrixStats) 

require(dplyr) 

require(doParallel)

#IMPORTANT - Load older version of xgboost before loading caret as the version used by caret has issues

xgb_lib <- path.expand("~/R_xgb_old")

.libPaths(c(xgb_lib, .libPaths()))

library(xgboost)

stopifnot(packageVersion("xgboost") == "1.7.8.1") 

require(caret)

require(MLmetrics)

require(randomForest) 

require(kernlab)

require(pROC) 

require(DESeq2)

# Please make sure the updated version of methylClass is loaded on every session

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::document("~/methylClass-fork") 
devtools::load_all("~/methylClass-fork")


# ---- Custom Wrapper Functions -----

source("~/Custom_functions.R",echo = TRUE)


#------------Load pre-processed train and test sets (in previous RNAseq_processing.R script)-------------------

RNA_dir<-"~/Thesis_models/RNA_dir/"

# Batch corrected with Combat-ref

train_NMB_sd_10k<-readRDS(file.path(path.expand(RNA_dir), "Train_RNAseq_sd10k.rds"))

train_NMB_mad_10k<-readRDS(file.path(path.expand(RNA_dir), "Train_RNAseq_mad10k.rds"))

train_NMB<-readRDS(file.path(path.expand(RNA_dir), "Train_RNAseq.rds"))

test_mat<-readRDS(file.path(path.expand(RNA_dir), "Test_RNAseq.rds")) 

train_pheno<-colData(readRDS(file.path(path.expand(RNA_dir),"Train_RNAseq_vst.rds")))

test_pheno<-colData(readRDS(file.path(path.expand(RNA_dir),"Test_RNAseq_vst.rds")))

# Batch corrected with limma

train_NMB_sd_10k_limma<-readRDS(file.path(path.expand(RNA_dir), "Train_RNAseq_sd10k_limma.rds"))

train_NMB_mad_10k_limma<-readRDS(file.path(path.expand(RNA_dir), "Train_RNAseq_mad10k_limma.rds"))

train_mat_limma<-readRDS(file.path(path.expand(RNA_dir), "Train_RNAseq_limma.rds"))

test_RNAseq_raw<-readRDS(file.path(path.expand(RNA_dir), "vst_test_sans_batch_corr.rds"))

test_mat_raw<-as.data.frame(t(assay(test_RNAseq_raw)))

test_pheno<-colData(test_RNAseq_raw)

#-------Caret Models with cross validation, on top 10k genes - on Train NMB corrected with combat-ref------------------

set.seed(42)

train_NMB_sd_10k<-train_NMB_sd_10k[match(rownames(train_NMB),rownames(train_NMB_sd_10k)),]

train_NMB_mad_10k<-train_NMB_mad_10k[match(rownames(train_NMB),rownames(train_NMB_mad_10k)),]

if (identical(rownames(train_NMB), rownames(train_NMB_sd_10k)) && identical(rownames(train_NMB), rownames(train_NMB_mad_10k))){
  
  folds<-createMultiFolds(train_NMB$histo, k = 10, times = 3)
  
} else { 
  message("IMPORTANT- mismatch in row order, to be consistent when training models on same n samples with different p features- MATCHING NOW")
  
  train_NMB_sd_10k<-train_NMB_sd_10k[match(rownames(train_NMB),rownames(train_NMB_sd_10k)),]
  
  train_NMB_mad_10k<-train_NMB_mad_10k[match(rownames(train_NMB),rownames(train_NMB_mad_10k)),]
  
  stopifnot(identical(rownames(train_NMB), rownames(train_NMB_sd_10k)))
  stopifnot(identical(rownames(train_NMB), rownames(train_NMB_mad_10k)))
  
  folds<-createMultiFolds(train_NMB$histo, k = 10, times = 3)
}

cv <- trainControl(
  method = "repeatedcv",
  index = folds,
  classProbs = TRUE,
  summaryFunction = multiClassSummary,
  savePredictions = "final", 
  allowParallel = TRUE,
) 


# Hyper-parameter tuning

# For Random Forests

rf_tunegrid <- expand.grid(mtry = c(50,75,100,150,200))

#For Support Vector Machines 

svm_grid <- expand.grid(
  C = 2^(-1:3),
  sigma = c(0.001,0.01,0.05,0.1)
) 

#For eXtreme Gradient Boosting

# xgb_grid <- expand.grid(
#   nrounds = 100,
#   max_depth = c(3, 5, 9),
#   eta = c(0.01,0.1,0.3),
#   gamma = 0,
#   colsample_bytree = c(0.2,0.5,0.8,1.0),
#   min_child_weight = 1,
#   subsample = 1
# )

# Parallel Processing 

n_cores<-20

clusters <-makeForkCluster(n_cores)

registerDoParallel(clusters)

# Training models on top 10k probes by Standard deviation 

#Random Forest

tryCatch({ 
  RF_sd_10k<- caret::train(histo~.,
                    data = train_NMB_sd_10k,
                    method = "rf",
                    metric = "logLoss",
                    trControl = cv, 
                    tuneGrid = rf_tunegrid,
                    maximize = FALSE
  ) 
  
  
  saveRDS(RF_sd_10k,"~/Thesis_models/RNA_dir/RNA_models/Combat-ref/Imb_RF_sd_10k_ref.rds")
  
  # Support Vector Machine
  
  SVM_sd_10k<- caret::train(histo~.,
                     data = train_NMB_sd_10k,
                     method = "svmRadial",
                     metric = "logLoss",
                     tuneGrid = svm_grid,
                     trControl = cv,
                     maximise = FALSE
  )
  
  
  saveRDS(SVM_sd_10k,"~/Thesis_models/RNA_dir/RNA_models/Combat-ref/Imb_SVM_sd_10k_ref.rds")
  
  
  # Training models on Top 10k probes by Median Absolute Deviation
  
  # Random Forest
  
  RF_mad_10k<- caret::train(histo~.,
                     data = train_NMB_mad_10k,
                     method = "rf",
                     metric = "logLoss",
                     tuneGrid = rf_tunegrid,
                     trControl = cv,
                     maximize = FALSE
  )
  
  saveRDS(RF_mad_10k,"~/Thesis_models/RNA_dir/RNA_models/Combat-ref/Imb_RF_mad_10k_ref.rds")
  
  #Support Vector Machine
  
  SVM_mad_10k<- caret::train(histo~.,
                      data = train_NMB_mad_10k,
                      method = "svmRadial",
                      metric = "logLoss",
                      tuneGrid = svm_grid,
                      trControl = cv,
                      maximise = FALSE
  )
  
  saveRDS(SVM_mad_10k,"~/Thesis_models/RNA_dir/RNA_models/Combat-ref/Imb_SVM_mad_10k_ref.rds")
  
  
  # Terminate parallel processing 
  
}, error = function(e) {
  
  message("Error pushed:", conditionMessage(e))
  
}, finally = {
  
  stopCluster(clusters) 
  
  registerDoSEQ()
  
  message("Parallel clusters terminated")
}) 


# ------------------- Training methylClass (Liu, 2024) models with all 441870 and/or top 10k input probes (Train NMB batch corrected with Combat-ref)- methylClass for high-dimensional methylation datasets
# is computationally inexpensive in comparison to caret as it applies parallel processing for each nested fit -----------------------------------------------------------------------------------------#

# Random Forest

RF_train_NMB<-maintrain(betas.. = train_NMB[,!colnames(train_NMB) %in% "histo"],
                        y.. = train_NMB$histo,
                        seed = 42,
                        method = "RF",
                        ntrees = 100, 
                        feature_sel = FALSE,
                        bal_training = FALSE,
                        topfeaturenumber = NULL, 
                        subset.CpGs = NULL,# can be set to 10000 for top variable probes
                        calibrationmethod = c("MR","FLR","LR")
)

saveRDS(RF_train_NMB,"~/Thesis_models/RNA_dir/RNA_models/Combat-ref/Imb_RF_NMB_allprobes_ref.rds") 

# Support Vector Machine - 5 fold cross validation

SVM_train_NMB<-maintrain( betas.. = train_NMB[,!colnames(train_NMB) %in% "histo"],
                          y.. = train_NMB$histo,
                          seed = 42,
                          method = "SVM",
                          modelcv = 5, #default
                          gridsearch = TRUE,
                          C.base = 10, #default
                          C.min = -3, #default
                          C.max = -2, #default
                          topfeaturenumber = NULL,
                          subset.CpGs = NULL,
                          calibrationmethod = c("MR","FLR","LR"))


saveRDS(SVM_train_NMB,"~/Thesis_models/RNA_dir/RNA_models/Combat-ref/Imb_SVM_NMB_allprobes_ref.rds")



#Training XGB models with methylClass- wrapper for caret XGB models

## There's an issue with caret's xgbTree's implementation that pushes the error "Something is wrong; all the Accuracy metric values are missing"
## It is an issue with the mismatch of the type of the object being created for each fit (the latest xgboost version; ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost)
## methylClass also applies caret::xgbTree and hence,a downgraded version of xgboost ("1.7.8.1") is to be applied


#IMPORTANT: Refresh the R session before installation 

xgb_lib <- path.expand("~/R_xgb_old")
dir.create(xgb_lib, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

remotes::install_version(
  "xgboost",
  version = "1.7.8.1",
  lib     = xgb_lib,
  repos   = "https://cloud.r-project.org"
)


xgb_lib <- path.expand("~/R_xgb_old")
.libPaths(c(xgb_lib, .libPaths()))

library(xgboost)
stopifnot(packageVersion("xgboost") == "1.7.8.1") 

message("xgboost ", as.character(packageVersion("xgboost")), " loaded from: ",
        find.package("xgboost"))


#IMPORTANT: Only load the forked methylClass repo to be able to apply recent commits

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::load_all("~/methylClass-fork") 


library(doParallel)

clusters <-makeForkCluster(n_cores)

registerDoParallel(clusters)

n_cores  <- 15
clusters <- makeForkCluster(n_cores)
registerDoParallel(clusters)

## Helper

run_xgb <- function(betas_df, y, out_file) {
  tryCatch({
    fit <- maintrain(
      betas.. = as.matrix(betas_df[, !colnames(betas_df) %in% "histo"]),
      y.. = y,
      seed = 1234,
      method = "XGB",
      modelcv = 5,                          # default
      gridsearch = TRUE,
      max_depth = 6,                        # default
      eta = c(0.1, 0.3),                    # default
      gamma = c(0, 0.01),                   # default
      min.chwght = 1,                       # default
      nrounds = 100,                        # default
      colsample_bytree = c(0.2, 0.5, 0.8, 1.0),   # same as caret XGB models
      topfeaturenumber = NULL,
      subset.CpGs = NULL,
      calibrationmethod = c("LR","FLR","MR") 
    )
    saveRDS(fit, out_file)
    message("Saved: ", out_file)
    fit
  }, error = function(e) {
    message("Error in ", basename(out_file), ": ", conditionMessage(e))
    NULL
  })
}


tryCatch({
  
  XGB_sd_10k <- run_xgb(train_NMB_sd_10k,
                        train_NMB_sd_10k$histo,
                        "~/Thesis_models/RNA_dir/RNA_models/Combat-ref/Imb_XGB_sd_10k_ref.rds")
  
  XGB_mad_10k <- run_xgb(train_NMB_mad_10k,
                         train_NMB_mad_10k$histo,     
                         "~/Thesis_models/RNA_dir/RNA_models/Combat-ref/Imb_XGB_mad_10k_ref.rds")
  
  # all probes
  
  XGB_train_NMB <- run_xgb(train_NMB,
                           train_NMB$histo,
                           "~/Thesis_models/RNA_dir/RNA_models/Combat-ref/Imb_XGB_NMB_allprobes_ref.rds")
  
}, finally = {
  stopCluster(clusters)
  registerDoSEQ()
  message("Parallel clusters terminated")
}) 


# ---- Testing and comparing of all the models  -----

# IMPORTANT: Load the old version of xgboost before running to avoid clash with caret loading the newer version 

xgb_lib <- path.expand("~/R_xgb_old")
.libPaths(c(xgb_lib, .libPaths()))

library(xgboost)
stopifnot(packageVersion("xgboost") == "1.7.8.1") 

# Load models if not in the env.

algos<-c("RF", "SVM", "XGB")

algo_names<-c(RF = "Random Forest", SVM = "Support Vector Machine", XGB = "eXtreme Gradient Boosting")


load_models<-function(algo,model_dir){
  files<-list.files(path.expand(model_dir),
                    pattern = paste0("^Imb_",algo,"_.*\\.rds$"),
                    full.names = TRUE)
  if(length(files) == 0) message (paste0("No files found for", algo, "in",path.expand(model_dir)))
  feature_selection <-gsub(paste0("^Imb_",algo,"_(NMB_)?|\\.rds$"),"",basename(files))
  setNames(lapply(files,readRDS), feature_selection)
}

models<-setNames(lapply(algos,load_models,model_dir = "~/Thesis_models/RNA_dir/RNA_models/Combat-ref/"),algos) 

# Testing on test NMB 

save_dir <- "~/Thesis_models/RNA_dir/Imb_Train_NMB/Combat-ref/"

dir.create(path.expand(save_dir), recursive = TRUE, showWarnings = FALSE)

test_mat_results<-list()

model_dict<-data.frame(model = character(), model_label = character())

for (algo in names(models)){
  for (feature_sel in names(models[[algo]])){
    
    model_key<-paste(algo, feature_sel, sep = "_")
    
    feature_label <- if (grepl("^sd", feature_sel)) {
      "top 10k by SD"
    } else if (grepl("^mad", feature_sel)) {
      "top 10k by MAD"
    } else {
      "all probes"
    }
    
    label  <- paste(algo_names[[algo]], feature_label, sep = " | ")
    
    message("Testing and evaluating ", label, " on Test NMB")
    
    model<-models[[algo]][[feature_sel]]
    
    test_data<-test_mat
    
    if (algo == "XGB") {
      fn <- model$mod[[1]]$feature_names              
      if (!is.null(fn)) {
        stopifnot("some probes the XGB model was trained on are missing from the test data" =
                    all(fn %in% colnames(test_mat)))
        test_data <- test_mat[, fn, drop = FALSE]
      }
    }
    
    res<-tryCatch(
      test_and_evaluate(
        model = model,
        test_data = test_data, 
        test_data_name = "Test RNAseq",
        calibrated = FALSE,
        train_set_name = "Train NMB RNAseq",
        true_test_labels = test_pheno$Histology,
        model_name = algo_names[[algo]],
        feature_selection_method = feature_label,
        file_id = paste(algo, feature_sel, sep = "_"),
        save_in = save_dir 
      ),
      error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL })
    
    if (is.null(res)) next
    
    res$summary$model_label<-label
    
    res$classwise$model_label<-label
    
    res$summary$algorithm<-algo                          
    
    res$summary$features<-feature_label      
    
    test_mat_results[[model_key]]<-res
    
    model_dict<-rbind(model_dict, data.frame(model=model_key, model_label = label))
    
    test_pheno<-as.data.frame(test_pheno)
    
    pred_probs<-res$pred_probs
    
    pred_probs<-pred_probs[match(test_pheno$Sample_ID, rownames(pred_probs)),] 
    
    test_RNAseq_heatmap <- prob_heatmaps(
      prob_df           = pred_probs,
      pred              = pred_probs$pred,
      test_pheno        = test_pheno,
      sample_id_column  = "Sample_ID",
      true_test_labels  = "Histology",
      mol_group         = "Group",
      sub_group         = "Subgroup",
      MYC_status        = "MYC",
      MYCN_status       = "MYCN",
      OS                = "OS",
      follow_up         = "Follow_up",
      Age               = "Age",
      cohort            = "Cohort",   
      column_for_panel  = "LCA",
      continuum         = "LCA",
      heatmap_label     = paste0( algo_names[[algo]], "for classification of histopathology,\ntrained on NMB RNAseq_",feature_label, "_and tested on non-CPR RNAseq cohorts"),
      heatmap_panel     = FALSE,
      file_id           = paste0(algo, "_com_ref_cal"),
      save_in           = file.path(path.expand(save_dir))
    )
  }
}

stopifnot("No model could be evaluated on Test NMB" = length(test_mat_results) > 0)

final_test_mat_summary<-do.call(rbind, lapply(test_mat_results, function(x) x$summary))

final_test_mat_classwise_summary<-do.call(rbind, lapply(test_mat_results, function(x) x$classwise))

print(final_test_mat_summary) 

write.csv(final_test_mat_summary,
          file.path(path.expand(save_dir), "imbalanced_summary_all_models.csv"), row.names = FALSE)
write.csv(final_test_mat_classwise_summary,
          file.path(path.expand(save_dir), "imbalanced_classwise_all_models.csv"), row.names = FALSE)

# Comparison Bar plot 

p_grouped_recall <- ggplot(final_test_mat_summary, aes(x = algorithm, y = macro_recall, fill = features)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(macro_recall, 3)),
            position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
  labs(title = wrap_title("Macro Recall Comparison by Machine Learning Model and Feature Set in Test RNAseq", width = 60),
       x = "Algorithm", y = "Macro Recall score", fill = "Feature Set") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
  theme_minimal()

ggsave("~/Thesis_models/RNA_dir/Imb_Train_NMB/Combat-ref/recall_comparision.png", plot = p_grouped_recall, width = 7, height = 5, dpi = 300)



p_grouped_auc <- ggplot(final_test_mat_summary, aes(x = algorithm, y = macro_auc, fill = features)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(macro_auc, 3)),
            position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
  labs(title = wrap_title("Area Under the Curve (AUC) score Comparison by Machine Learning Model and Feature Set in Test RNAseq", width = 60),
       x = "Algorithm", y = "AUC score", fill = "Feature Set") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
  theme_minimal() 

ggsave("~/Thesis_models/RNA_dir/Imb_Train_NMB/Combat-ref/auc_comparision.png", plot = p_grouped_auc, width = 7, height = 5, dpi = 300)




#-------Caret Models with cross validation, on top 10k genes - on Train NMB corrected with limma------------------

set.seed(42)

train_NMB_sd_10k_limma<-train_NMB_sd_10k_limma[match(rownames(train_mat_limma),rownames(train_NMB_sd_10k_limma)),]

train_NMB_mad_10k_limma<-train_NMB_mad_10k_limma[match(rownames(train_mat_limma),rownames(train_NMB_mad_10k_limma)),]

if (identical(rownames(train_mat_limma), rownames(train_NMB_sd_10k_limma)) && identical(rownames(train_mat_limma), rownames(train_NMB_mad_10k_limma))){
  
  folds<-createMultiFolds(train_mat_limma$histo, k = 10, times = 3)
  
} else { 
  message("IMPORTANT- mismatch in row order, to be consistent when training models on same n samples with different p features- MATCHING NOW")
  
  train_NMB_sd_10k_limma<-train_NMB_sd_10k_limma[match(rownames(train_mat_limma),rownames(train_NMB_sd_10k_limma)),]
  
  train_NMB_mad_10k_limma<-train_NMB_mad_10k_limma[match(rownames(train_mat_limma),rownames(train_NMB_mad_10k_limma)),]
  
  stopifnot(identical(rownames(train_mat_limma), rownames(train_NMB_sd_10k_limma)))
  stopifnot(identical(rownames(train_mat_limma), rownames(train_NMB_mad_10k_limma)))
  
  folds<-createMultiFolds(train_mat_limma$histo, k = 10, times = 3)
}

cv <- trainControl(
  method = "repeatedcv",
  index = folds,
  classProbs = TRUE,
  summaryFunction = multiClassSummary,
  savePredictions = "final", 
  allowParallel = TRUE,
) 


# Hyper-parameter tuning

# For Random Forests

rf_tunegrid <- expand.grid(mtry = c(50,75,100,150,200))

#For Support Vector Machines 

svm_grid <- expand.grid(
  C = 2^(-1:3),
  sigma = c(0.001,0.01,0.05,0.1)
) 

#For eXtreme Gradient Boosting

# xgb_grid <- expand.grid(
#   nrounds = 100,
#   max_depth = c(3, 5, 9),
#   eta = c(0.01,0.1,0.3),
#   gamma = 0,
#   colsample_bytree = c(0.2,0.5,0.8,1.0),
#   min_child_weight = 1,
#   subsample = 1
# )

# Parallel Processing 

n_cores<-20

clusters <-makeForkCluster(n_cores)

registerDoParallel(clusters)

# Training models on top 10k probes by Standard deviation 

#Random Forest

tryCatch({ 
  RF_sd_10k<- caret::train(histo~.,
                    data = train_NMB_sd_10k_limma,
                    method = "rf",
                    metric = "logLoss",
                    trControl = cv, 
                    tuneGrid = rf_tunegrid,
                    maximize = FALSE
  ) 
  
  
  saveRDS(RF_sd_10k,"~/Thesis_models/RNA_dir/RNA_models/limma/Imb_RF_sd_10k_limma.rds")
  
  # Support Vector Machine
  
  SVM_sd_10k<- caret::train(histo~.,
                     data = train_NMB_sd_10k_limma,
                     method = "svmRadial",
                     metric = "logLoss",
                     tuneGrid = svm_grid,
                     trControl = cv,
                     maximise = FALSE
  )
  
  
  saveRDS(SVM_sd_10k,"~/Thesis_models/RNA_dir/RNA_models/limma/Imb_SVM_sd_10k_limma.rds")
  
  
  # Training models on Top 10k probes by Median Absolute Deviation
  
  # Random Forest
  
  RF_mad_10k<- caret::train(histo~.,
                     data = train_NMB_mad_10k_limma,
                     method = "rf",
                     metric = "logLoss",
                     tuneGrid = rf_tunegrid,
                     trControl = cv,
                     maximize = FALSE
  )
  
  saveRDS(RF_mad_10k,"~/Thesis_models/RNA_dir/RNA_models/limma/Imb_RF_mad_10k_limma.rds")
  
  #Support Vector Machine
  
  SVM_mad_10k<- caret::train(histo~.,
                      data = train_NMB_mad_10k_limma,
                      method = "svmRadial",
                      metric = "logLoss",
                      tuneGrid = svm_grid,
                      trControl = cv,
                      maximise = FALSE
  )
  
  saveRDS(SVM_mad_10k,"~/Thesis_models/RNA_dir/RNA_models/limma/Imb_SVM_mad_10k_limma.rds")
  
  
  # Terminate parallel processing 
  
}, error = function(e) {
  
  message("Error pushed:", conditionMessage(e))
  
}, finally = {
  
  stopCluster(clusters) 
  
  registerDoSEQ()
  
  message("Parallel clusters terminated")
}) 


# ------------------- Training methylClass (Liu, 2024) models with all 441870 and/or top 10k input probes (Train NMB batch corrected with limma)- methylClass for high-dimensional methylation datasets
# is computationally inexpensive in comparison to caret as it applies parallel processing for each nested fit -----------------------------------------------------------------------------------------#

# Random Forest

RF_train_NMB<-maintrain(betas.. = train_mat_limma[,!colnames(train_mat_limma) %in% "histo"],
                        y.. = train_mat_limma$histo,
                        seed = 42,
                        method = "RF",
                        ntrees = 100, 
                        feature_sel = FALSE,
                        bal_training = FALSE,
                        topfeaturenumber = NULL, 
                        subset.CpGs = NULL,# can be set to 10000 for top variable probes
                        calibrationmethod = c("MR","FLR","LR")
)

saveRDS(RF_train_NMB,"~/Thesis_models/RNA_dir/RNA_models/limma/Imb_RF_NMB_allprobes_limma.rds") 

# Support Vector Machine - 5 fold cross validation

SVM_train_NMB<-maintrain( betas.. = train_mat_limma[,!colnames(train_mat_limma) %in% "histo"],
                          y.. = train_mat_limma$histo,
                          seed = 42,
                          method = "SVM",
                          modelcv = 5, #default
                          gridsearch = TRUE,
                          C.base = 10, #default
                          C.min = -3, #default
                          C.max = -2, #default
                          topfeaturenumber = NULL,
                          subset.CpGs = NULL,
                          calibrationmethod = c("MR","FLR","LR"))


saveRDS(SVM_train_NMB,"~/Thesis_models/RNA_dir/RNA_models/limma/Imb_SVM_NMB_allprobes_limma.rds")



#Training XGB models with methylClass- wrapper for caret XGB models

## There's an issue with caret's xgbTree's implementation that pushes the error "Something is wrong; all the Accuracy metric values are missing"
## It is an issue with the mismatch of the type of the object being created for each fit (the latest xgboost version; ALTLIST classes must provide a Set_elt method [class: XGBAltrepPointerClass, pkg: xgboost)
## methylClass also applies caret::xgbTree and hence,a downgraded version of xgboost ("1.7.8.1") is to be applied


#IMPORTANT: Refresh the R session before installation 

xgb_lib <- path.expand("~/R_xgb_old")
dir.create(xgb_lib, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

remotes::install_version(
  "xgboost",
  version = "1.7.8.1",
  lib     = xgb_lib,
  repos   = "https://cloud.r-project.org"
)


xgb_lib <- path.expand("~/R_xgb_old")
.libPaths(c(xgb_lib, .libPaths()))

library(xgboost)
stopifnot(packageVersion("xgboost") == "1.7.8.1") 

message("xgboost ", as.character(packageVersion("xgboost")), " loaded from: ",
        find.package("xgboost"))


#IMPORTANT: Only load the forked methylClass repo to be able to apply recent commits

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::load_all("~/methylClass-fork") 


library(doParallel)

clusters <-makeForkCluster(n_cores)

registerDoParallel(clusters)

n_cores  <- 15
clusters <- makeForkCluster(n_cores)
registerDoParallel(clusters)

## Helper

run_xgb <- function(betas_df, y, out_file) {
  tryCatch({
    fit <- maintrain(
      betas.. = as.matrix(betas_df[, !colnames(betas_df) %in% "histo"]),
      y.. = y,
      seed = 1234,
      method = "XGB",
      modelcv = 5,                          # default
      gridsearch = TRUE,
      max_depth = 6,                        # default
      eta = c(0.1, 0.3),                    # default
      gamma = c(0, 0.01),                   # default
      min.chwght = 1,                       # default
      nrounds = 100,                        # default
      colsample_bytree = c(0.2, 0.5, 0.8, 1.0),   # same as caret XGB models
      topfeaturenumber = NULL,
      subset.CpGs = NULL,
      calibrationmethod = c("LR","FLR","MR") 
    )
    saveRDS(fit, out_file)
    message("Saved: ", out_file)
    fit
  }, error = function(e) {
    message("Error in ", basename(out_file), ": ", conditionMessage(e))
    NULL
  })
}


tryCatch({
  
  # XGB_sd_10k <- run_xgb(train_NMB_sd_10k_limma,
  #                       train_NMB_sd_10k_limma$histo,
  #                       "~/Thesis_models/RNA_dir/RNA_models/limma/Imb_XGB_sd_10k_limma.rds")
  # 
  # XGB_mad_10k <- run_xgb(train_NMB_mad_10k_limma,
  #                        train_NMB_mad_10k_limma$histo,     
  #                        "~/Thesis_models/RNA_dir/RNA_models/limma/Imb_XGB_mad_10k_limma.rds")
  # 
  # all probes
  
  XGB_train_NMB <- run_xgb(train_NMB_mat_limma,
                           train_NMB_mat_limma$histo,
                           "~/Thesis_models/RNA_dir/RNA_models/limma/Imb_XGB_NMB_allprobes_limma.rds")
  
}, finally = {
  stopCluster(clusters)
  registerDoSEQ()
  message("Parallel clusters terminated")
}) 


# ---- Testing and comparing of all the models  -----

# IMPORTANT: Load the old version of xgboost before running to avoid clash with caret loading the newer version 

xgb_lib <- path.expand("~/R_xgb_old")
.libPaths(c(xgb_lib, .libPaths()))

library(xgboost)
stopifnot(packageVersion("xgboost") == "1.7.8.1") 

# Load models if not in the env.

algos<-c("RF", "SVM", "XGB")

algo_names<-c(RF = "Random Forest", SVM = "Support Vector Machine", XGB = "eXtreme Gradient Boosting")


load_models<-function(algo,model_dir){
  files<-list.files(path.expand(model_dir),
                    pattern = paste0("^Imb_",algo,"_.*\\.rds$"),
                    full.names = TRUE)
  if(length(files) == 0) message (paste0("No files found for", algo, "in",path.expand(model_dir)))
  feature_selection <-gsub(paste0("^Imb_",algo,"_(NMB_)?|\\.rds$"),"",basename(files))
  setNames(lapply(files,readRDS), feature_selection)
}

models<-setNames(lapply(algos,load_models,model_dir = "~/Thesis_models/RNA_dir/RNA_models/limma/"),algos) 

# Testing on test NMB 

save_dir <- "~/Thesis_models/RNA_dir/Imb_Train_NMB/limma/"

dir.create(path.expand(save_dir), recursive = TRUE, showWarnings = FALSE)

test_mat_raw<-test_mat_raw[match(test_pheno$Sample_ID, rownames(test_mat_raw)),]

test_mat_results<-list()

model_dict<-data.frame(model = character(), model_label = character())

for (algo in names(models)){
  for (feature_sel in names(models[[algo]])){
    
    model_key<-paste(algo, feature_sel, sep = "_")
    
    feature_label <- if (grepl("^sd", feature_sel)) {
      "top 10k by SD"
    } else if (grepl("^mad", feature_sel)) {
      "top 10k by MAD"
    } else {
      "all probes"
    }
    
    label  <- paste(algo_names[[algo]], feature_label, sep = " | ")
    
    message("Testing and evaluating ", label, " on Test NMB")
    
    model<-models[[algo]][[feature_sel]]
    
    test_data<-test_mat_raw
    
    if (algo == "XGB") {
      fn <- model$mod[[1]]$feature_names              
      if (!is.null(fn)) {
        stopifnot("some probes the XGB model was trained on are missing from the test data" =
                    all(fn %in% colnames(test_mat_raw)))
        test_data <- test_mat_raw[, fn, drop = FALSE]
      }
    }
    
    res<-tryCatch(
      test_and_evaluate(
        model = model,
        test_data = test_data, 
        test_data_name = "Test RNAseq",
        calibrated = FALSE,
        train_set_name = "Train NMB RNAseq-limma",
        true_test_labels = test_pheno$Histology,
        model_name = algo_names[[algo]],
        feature_selection_method = feature_label,
        file_id = paste(algo, feature_sel, sep = "_"),
        save_in = save_dir 
      ),
      error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL })
    
    if (is.null(res)) next
    
    res$summary$model_label<-label
    
    res$classwise$model_label<-label
    
    res$summary$algorithm<-algo                          
    
    res$summary$features<-feature_label      
    
    test_mat_results[[model_key]]<-res
    
    model_dict<-rbind(model_dict, data.frame(model=model_key, model_label = label))
    
    pred_probs<-res$pred_probs
    
    test_pheno<-as.data.frame(test_pheno)
    
    pred_probs<-pred_probs[match(test_pheno$Sample_ID, rownames(pred_probs)),] 
    
    test_RNAseq_heatmap <- prob_heatmaps(
      prob_df           = pred_probs,
      pred              = pred_probs$pred,
      test_pheno        = test_pheno,
      sample_id_column  = "Sample_ID",
      true_test_labels  = "Histology",
      mol_group         = "Group",
      sub_group         = "Subgroup",
      MYC_status        = "MYC",
      MYCN_status       = "MYCN",
      OS                = "OS",
      follow_up         = "Follow_up",
      Age               = "Age",
      cohort            = "Cohort",   
      column_for_panel  = "LCA",
      continuum         = "LCA",
      heatmap_label     = paste0( algo_names[[algo]], "for classification of histopathology,\ntrained on NMB RNAseq batch corrected with limma_",feature_label, "_and tested on non-CPR RNAseq cohorts"),
      heatmap_panel     = FALSE,
      file_id           = paste0(algo, "_limma_cal"),
      save_in           = file.path(path.expand(save_dir))
    )
  }
}

stopifnot("No model could be evaluated on Test NMB" = length(test_mat_results) > 0)

final_test_mat_summary<-do.call(rbind, lapply(test_mat_results, function(x) x$summary))

final_test_mat_classwise_summary<-do.call(rbind, lapply(test_mat_results, function(x) x$classwise))

print(final_test_mat_summary) 

write.csv(final_test_mat_summary,
          file.path(path.expand(save_dir), "imbalanced_summary_all_models.csv"), row.names = FALSE)
write.csv(final_test_mat_classwise_summary,
          file.path(path.expand(save_dir), "imbalanced_classwise_all_models.csv"), row.names = FALSE)

# Comparison Bar plot 

p_grouped_recall <- ggplot(final_test_mat_summary, aes(x = algorithm, y = macro_recall, fill = features)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(macro_recall, 3)),
            position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
  labs(title = wrap_title("Macro Recall Comparison by Machine Learning Model and Feature Set  (BC with limma) in Test RNAseq", width = 60),
       x = "Algorithm", y = "Macro Recall score", fill = "Feature Set") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
  theme_minimal()

ggsave("~/Thesis_models/RNA_dir/Imb_Train_NMB/limma/recall_comparision.png", plot = p_grouped_recall, width = 7, height = 5, dpi = 300)



p_grouped_auc <- ggplot(final_test_mat_summary, aes(x = algorithm, y = macro_auc, fill = features)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(macro_auc, 3)),
            position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
  labs(title = wrap_title("Area Under the Curve (AUC) score Comparison by Machine Learning Model and Feature Set(BC with limma) in Test RNAseq", width = 60),
       x = "Algorithm", y = "AUC score", fill = "Feature Set") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
  theme_minimal() 

ggsave("~/Thesis_models/RNA_dir/Imb_Train_NMB/limma/auc_comparision.png", plot = p_grouped_auc, width = 7, height = 5, dpi = 300)


#-----------------Splitting CPR NMB Cohort to train and test models on CPR cohorts---------------------------

NMB_Final_Cohort<-readRDS("~/NMB/NMB_Final_Cohort_26.rds")

NMB_Final_Cohort<-NMB_Final_Cohort[match(train_pheno$Sample_ID, NMB_Final_Cohort$Sample_ID2),]

stopifnot(all(train_pheno$Sample_ID == NMB_Final_Cohort$Sample_ID2))

train_pheno$FS<-as.factor(NMB_Final_Cohort$`Feasibility Study`)

train_pheno_non_FS<-train_pheno[train_pheno$FS == "0",]

test_pheno_FS<-train_pheno[train_pheno$FS == "1",] 

train_non_FS<-train_NMB[match(train_pheno_non_FS$Sample_ID, rownames(train_NMB)),]

test_FS<-train_NMB[match(test_pheno_FS$Sample_ID, rownames(train_NMB)),]

# Training RF and SVM models on non-FS

RF_train_non_FS<-maintrain(betas.. = train_non_FS[,!colnames(train_non_FS) %in% "histo"],
                        y.. = train_non_FS$histo,
                        seed = 42,
                        method = "RF",
                        ntrees = 100, 
                        feature_sel = FALSE,
                        bal_training = FALSE,
                        topfeaturenumber = NULL, 
                        subset.CpGs = NULL,# can be set to 10000 for top variable probes
                        calibrationmethod = c("MR","FLR","LR")
)

saveRDS(RF_train_non_FS,"~/Thesis_models/RNA_dir/RNA_models/Imb_RF_non_FS_all_probes.rds") 



RF_batch_corr<-readRDS(file.path(path.expand(batch_corr_dir), "RF_mol_batch_corr.rds")) 

message("Testing and evaluating on Test NMB")

if (all(rownames(test_NMB) == NMB_test_pheno$Sample_ID)) {
  test_NMB_results_adj<- tryCatch(
    test_and_evaluate (
      model = RF_batch_corr,
      method = "methylClass",
      calibrated = TRUE,
      calibration_method = "MR",
      train_set_name = "Train NMB adjusted for Molecular Group",
      test_data = test_NMB_adj,
      test_data_name = "Test NMB",
      true_test_labels = as.factor(NMB_test_pheno$CPR_Histology),
      model_name = "Random Forest",
      save_in = file.path(path.expand(batch_corr_dir_2))
    ),
    error = function(e) {
      message(sprintf("FAILED [%s / %s]: %s",
                      " Batch corr Random Forest trained on Train NMB", "Test NMB", conditionMessage(e)))
      NULL
    }
  )
} else {
  stop("Mismatch in test set and test pheno Sample IDs")
}

Test_NMB_heatmap_batch_corr <- prob_heatmaps(
  prob_df           = test_NMB_results_adj$pred_probs,
  pred              = test_NMB_results_adj$pred_probs$pred,
  test_pheno        = NMB_test_pheno,
  sample_id_column  = "Sample_ID",
  true_test_labels  = "CPR_Histology",
  mol_group         = "Mol_group",
  sub_group         = "sub_group",
  MYC_status        = "MYC Status",
  MYCN_status       = "MYCN Status",
  OS                = "OS",
  follow_up         = "Follow up years",
  Age               = "Age",
  continuum         = "LCA",
  heatmap_label     = "Random Forest Model for classification of histopathology,\ntrained on 'Train NMB' after adjusting for Molecular group, tested on 'Test_NMB",
  heatmap_panel     = FALSE,
  column_for_panel  = "LCA",
  file_id           = "RF_corr_test_NMB",
  save_in           = file.path(path.expand(batch_corr_dir))
)

# On Cavalli


RF_batch_corr<-readRDS(file.path(path.expand(batch_corr_dir), "RF_mol_batch_corr.rds")) 

if (all(rownames(Cavalli_betas) == Cavalli_pheno$Sample_ID)) {
  test_Cavalli_results_adj<- tryCatch(
    test_and_evaluate (
      model = RF_batch_corr,
      method = "methylClass",
      calibrated = TRUE,
      calibration_method = "MR",
      train_set_name = "Train NMB Adjusted for Molecular Group",
      test_data = Cavalli_adj,
      test_data_name = "Cavalli",
      true_test_labels = as.factor(Cavalli_pheno$histology),
      model_name = "Random Forest",
      save_in = file.path(path.expand(batch_corr_dir_2))
    ),
    error = function(e) {
      message(sprintf("FAILED [%s / %s]: %s",
                      " Batch corr Random Forest trained on Train NMB", "Cavalli", conditionMessage(e)))
      NULL
    }
  )
} else {
  stop("Mismatch in test set and test pheno Sample IDs")
}

Test_Cavalli_heatmap_batch_corr <- prob_heatmaps(
  prob_df           = test_Cavalli_results_adj$pred_probs,
  pred              = test_Cavalli_results_adj$pred_probs$pred,
  test_pheno        = Cavalli_pheno,
  sample_id_column  = "Sample_ID",
  true_test_labels  = "histology",
  mol_group         = "Mol_grp",
  sub_group         = "sub_grp",
  MYC_status        = "MYC_Status",
  MYCN_status       = "MYCN_Status",
  OS                = "Dead",
  follow_up         = "OS (years)",
  Age               = "Age",
  continuum         = "LCA",
  heatmap_label     = "Random Forest Model for classification of histopathology,\ntrained on 'Train NMB' after adjusting for Molecular group, tested on 'Cavalli",
  heatmap_panel     = FALSE,
  column_for_panel  = "LCA",
  file_id           = "RF_corr_Cavalli",
  save_in           = file.path(path.expand(batch_corr_dir))
)




# Support Vector Machine - 5 fold cross validation

SVM_train_non_FS<-maintrain( betas.. = train_non_FS[,!colnames(train_non_FS) %in% "histo"],
                          y.. = train_non_FS$histo,
                          seed = 42,
                          method = "SVM",
                          modelcv = 5, #default
                          gridsearch = TRUE,
                          C.base = 10, #default
                          C.min = -3, #default
                          C.max = -2, #default
                          topfeaturenumber = NULL,
                          subset.CpGs = NULL,
                          calibrationmethod = c("MR","FLR","LR"))


saveRDS(SVM_train_NMB,"~/Thesis_models/RNA_dir/RNA_models/Imb_SVM_non_FS_allprobes.rds") 


