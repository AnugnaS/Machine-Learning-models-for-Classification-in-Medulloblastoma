#===================================================================================================================#
#Train Machine Learning Models for Classification of Histopathology in Medulloblastoma

# This script includes loading of samplesheets and pre processed methylation data (in-house),
# Feature selection of Top most variable probes,
# Initial training of models on imbalanced data, 
# Training and optimising of caret and methylClass models- Random Forest(RF),Support Vector Machine (SVM) and eXtreme Gradient Boosting (XGB,
# Testing of models on independent test sets and subsequent evaluation metrics
# Balancing the train set with Sampling Methods
# RF trained on balanced data and evaluated on test sets 
# Optimising the Train Feature Set with Principal Component Selection and Feature selection with LASSO
# Stratifying NMB to form the non-Feasibility Train set and Feasibility test set as in (Crosier et al., 2021). 
# Training RF models on train set balanced with molecular group=proportion-aware under-sampling, with ablation of custom methylClass hyperparameters
# SHAP for overlap features for histology and Principal Molecular group of MB.
# Cox models for survival association of RF model's probability estimates
# Hazard Ratio Forest Plots
# Youden's J statistic Probability score thresholds and sequential application
# Definition of Uncertainty zone to depict CLA-LCA continuum reflected in probability scores of assignments.

# PLEASE NOTE: Every section loads the datasets and functions, and updates hyperparameter settings as for a fresh session. 
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

# Please make sure the updated version of methylClass is loaded on every session

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::document("~/methylClass-fork") 
devtools::load_all("~/methylClass-fork")


# ---- Custom Wrapper Functions -----

source("~/Custom_functions.R",echo = TRUE)

# ---- Load Data -----

#Train and Test split of the 'NMB' Cohort is done with caret::createDataPartition 
#in the previous script "Distribution_charts.R". "NOS" histology labels have been removed before the split

NMB_train_pheno<-readRDS("~/NMB/NMB_train_pheno.rds")

if (any(is.na(NMB_train_pheno$Sample_ID))){
  NMB_train_pheno<-NMB_train_pheno[!is.na(NMB_train_pheno$Sample_ID),]
}

NMB_test_pheno<-readRDS("~/NMB/NMB_test_pheno.rds")

##NMB_test_pheno<-NMB_test_pheno[NMB_test_pheno != "NMB_411_450k",] #missing

if (any(is.na(NMB_test_pheno$Sample_ID))){
  NMB_test_pheno<-NMB_test_pheno[!is.na(NMB_test_pheno$Sample_ID),]
}

#TODO: Save the above as .csv sample sheets

# Beta value matrix post pre processing with minfi- 450k + EPIC combined

NMB_betas<-readRDS("~/NMB/Final_NMB_betas.rds")

train_NMB<-NMB_betas[match(NMB_train_pheno$Sample_ID, rownames(NMB_betas)),]

test_NMB<-NMB_betas[match(NMB_test_pheno$Sample_ID, rownames(NMB_betas)),]

# Independent Test dataset- Cavalli (Cavalli et al., 2017)

Cavalli_pheno<- readRDS("~/Cavalli/Final_Cavalli_pheno.rds")

#TODO: Save the above as .csv sample sheet

if (any(is.na(Cavalli_pheno$Sample_ID))){
  Cavalli_pheno<-Cavalli_pheno[!is.na(Cavalli_pheno$Sample_ID),] 
}

Cavalli_betas<-readRDS("~/Cavalli/Cavalli-test/Cavalli_test_betas.rds") 

Cavalli_betas<-as.data.frame(t(Cavalli_betas)) 

Cavalli_betas<-Cavalli_betas[match(Cavalli_pheno$Sample_ID, rownames(Cavalli_betas)),] 

# ------ Prepare data ------ #

# Top 10k probes by Standard Deviation

train_NMB_mat<- as.matrix(train_NMB)

if (any(is.na(train_NMB_mat))) { 
  message("NAs in train_NMB:using na.rm = TRUE")
  train_NMB_sds<-colSds(train_NMB_mat, na.rm = TRUE)
} else {
  train_NMB_sds<-colSds(train_NMB_mat, na.rm = FALSE)
} 

train_NMB_sd_10k<-train_NMB[,order(train_NMB_sds,decreasing = TRUE)[1:10000]]

if (all(rownames(train_NMB_sd_10k) == NMB_train_pheno$Sample_ID)){
  train_NMB_sd_10k$histo<-as.factor(NMB_train_pheno$CPR_Histology)
} else { 
  NMB_train_pheno<-NMB_train_pheno[match(rownames(train_NMB_sd_10k), NMB_train_pheno$Sample_ID),]
  train_NMB_sd_10k$histo<-as.factor(NMB_train_pheno$CPR_Histology)
}

#train_NMB_10k<-train_NMB[,order(apply(train_NMB,2,sd), decreasing = TRUE)[1:10000]]

#Top 10k probes by Median Absolute Deviation

if (any(is.na(train_NMB_mat))) { 
  message("NAs in train_NMB:using na.rm = TRUE")
  train_NMB_mads<-colMads(train_NMB_mat, na.rm = TRUE)
} else {
  train_NMB_mads<-colMads(train_NMB_mat, na.rm = FALSE)
} 

train_NMB_mad_10k<-train_NMB[,order(train_NMB_mads,decreasing = TRUE)[1:10000]]

#train_NMB_10k<-train_NMB[,order(apply(train_NMB,2,mad), decreasing = TRUE)[1:10000]]


if (all(rownames(train_NMB_mad_10k) == NMB_train_pheno$Sample_ID)){
  train_NMB_mad_10k$histo<-as.factor(NMB_train_pheno$CPR_Histology)
} else { 
  NMB_train_pheno<-NMB_train_pheno[match(rownames(train_NMB_mad_10k), NMB_train_pheno$Sample_ID),]
  train_NMB_mad_10k$histo<-as.factor(NMB_train_pheno$CPR_Histology)
}


# Add Histology label annotation to NMB train set. 

if (all(rownames(train_NMB) == NMB_train_pheno$Sample_ID)){
  train_NMB$histo<-as.factor(NMB_train_pheno$CPR_Histology)
} else { 
  NMB_train_pheno<-NMB_train_pheno[match(rownames(train_NMB), NMB_train_pheno$Sample_ID),]
  train_NMB$histo<-as.factor(NMB_train_pheno$CPR_Histology)
}

#Test sets

if (any(is.na(NMB_test_pheno$CPR_Histology) | NMB_test_pheno$CPR_Histology == "NOS")){
  message("Check Test NMB Histology Labels")
}

NMB_test_pheno<-NMB_test_pheno[match(rownames(test_NMB),NMB_test_pheno$Sample_ID),]

Cavalli_pheno<-Cavalli_pheno[!is.na(Cavalli_pheno$histology) & Cavalli_pheno$histology != "NOS",]

Cavalli_betas<-Cavalli_betas[match(Cavalli_pheno$Sample_ID, rownames(Cavalli_betas)),]

if (all(Cavalli_pheno$Sample_ID == rownames(Cavalli_betas))){
  Cavalli_pheno<-Cavalli_pheno[match(rownames(Cavalli_betas), Cavalli_pheno$Sample_ID),]
} 

#----------- Training caret (Kuhn, 2007) models with top 10k probes with 5 fold cross validation
# with parallel processing ----------------------------------------------------------------------#

# Folds for Cross Validation - 5 folds and 3 repeats

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
  stopifnot(identical(ownames(train_NMB), rownames(train_NMB_mad_10k)))
  
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
RF_sd_10k<- train(histo~.,
                             data = train_NMB_sd_10k,
                             method = "rf",
                             metric = "logLoss",
                             trControl = cv, 
                             tuneGrid = rf_tunegrid,
                             maximize = FALSE
) 


saveRDS(RF_sd_10k,"~/Thesis_models/Imb_RF_sd_10k.rds")

# Support Vector Machine

SVM_sd_10k<- train(histo~.,
                   data = train_NMB_sd_10k,
                   method = "svmRadial",
                   metric = "logLoss",
                   tuneGrid = svm_grid,
                   trControl = cv,
                   maximise = FALSE
)


saveRDS(SVM_sd_10k,"~/Thesis_models/Imb_SVM_sd_10k.rds")

  
# Training models on Top 10k probes by Median Absolute Deviation

# Random Forest

RF_mad_10k<- train(histo~.,
                  data = train_NMB_mad_10k,
                  method = "rf",
                  metric = "logLoss",
                  tuneGrid = rf_tunegrid,
                  trControl = cv,
                  maximize = FALSE
)

saveRDS(RF_mad_10k,"~/Thesis_models/Imb_RF_mad_10k.rds")

#Support Vector Machine

SVM_mad_10k<- train(histo~.,
                   data = train_NMB_mad_10k,
                   method = "svmRadial",
                   metric = "logLoss",
                   tuneGrid = svm_grid,
                   trControl = cv,
                   maximise = FALSE
)

saveRDS(SVM_mad_10k,"~/Thesis_models/Imb_SVM_mad_10k.rds")


# Terminate parallel processing 

}, error = function(e) {
  
  message("Error pushed:", conditionMessage(e))
  
}, finally = {
  
  stopCluster(clusters) 
  
  registerDoSEQ()
  
  message("Parallel clusters terminated")
}) 


# ------------------- Training methylClass (Liu, 2024) models with all 441870 and/or top 10k input probes- methylClass for high-dimensional methylation datasets
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

saveRDS(RF_train_NMB,"~/Thesis_models/Imb_RF_NMB_allprobes.rds") 

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


saveRDS(SVM_train_NMB,"~/Thesis_models/Imb_SVM_NMB_allprobes.rds")



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
      subset.CpGs = NULL
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
                        "~/Thesis_models/Imb_XGB_sd_10k.rds")
  
  XGB_mad_10k <- run_xgb(train_NMB_mad_10k,
                         train_NMB_mad_10k$histo,     
                         "~/Thesis_models/Imb_XGB_mad_10k.rds")
  
  # all probes
  
  XGB_train_NMB <- run_xgb(train_NMB,
                           train_NMB$histo,
                           "~/Thesis_models/Imb_XGB_NMB_allprobes.rds")
  
}, finally = {
  stopCluster(clusters)
  registerDoSEQ()
  message("Parallel clusters terminated")
}) 


# ---- Testing and comparing of all the models -----

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

models<-setNames(lapply(algos,load_models,model_dir = "~/Thesis_models/Imb_Train_NMB/"),algos) 

# Testing on test NMB 

save_dir <- "~/Thesis_models/test_NMB_imbal/"

dir.create(path.expand(save_dir), recursive = TRUE, showWarnings = FALSE)

test_NMB_results<-list()

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
    
    test_data<-test_NMB
    
    if (algo == "XGB") {
      fn <- model$mod[[1]]$feature_names              
      if (!is.null(fn)) {
        stopifnot("some probes the XGB model was trained on are missing from the test data" =
                    all(fn %in% colnames(test_NMB)))
        test_data <- test_NMB[, fn, drop = FALSE]
      }
    }
    
    res<-tryCatch(
      test_and_evaluate(
        model = model,
        test_data = test_data, 
        test_data_name = "Test NMB",
        train_set_name = "Train NMB",
        true_test_labels = NMB_test_pheno$CPR_Histology,
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
    
    test_NMB_results[[model_key]]<-res
    
    model_dict<-rbind(model_dict, data.frame(model=model_key, model_label = label))
  }
}
 
stopifnot("No model could be evaluated on Test NMB" = length(test_NMB_results) > 0)

final_test_NMB_summary<-do.call(rbind, lapply(test_NMB_results, function(x) x$summary))

final_test_NMB_classwise_summary<-do.call(rbind, lapply(test_NMB_results, function(x) x$classwise))

print(final_test_NMB_summary) 

write.csv(final_test_NMB_summary,
          file.path(path.expand(save_dir), "imbalanced_summary_all_models.csv"), row.names = FALSE)
write.csv(final_test_NMB_classwise_summary,
          file.path(path.expand(save_dir), "imbalanced_classwise_all_models.csv"), row.names = FALSE)

# Comparison Bar plot 

p_grouped_accuracy <- ggplot(final_test_NMB_summary, aes(x = algorithm, y = accuracy, fill = features)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(accuracy, 3)),
            position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
  labs(title = wrap_title("Accuracy Comparison by Machine Learning Model and Feature Set in Test NMB", width = 60),
       x = "Algorithm", y = "Accuracy score", fill = "Feature Set") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
  theme_minimal()

ggsave("~/Thesis_models/test_NMB_imbal/accuracy_comparision.png", plot = p_grouped_accuracy, width = 7, height = 5, dpi = 300)



p_grouped_auc <- ggplot(final_test_NMB_summary, aes(x = algorithm, y = macro_auc, fill = features)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(macro_auc, 3)),
            position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
  labs(title = wrap_title("Area Under the Curve (AUC) score Comparison by Machine Learning Model and Feature Set in Test NMB", width = 60),
       x = "Algorithm", y = "AUC score", fill = "Feature Set") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
  theme_minimal() 

ggsave("~/Thesis_models/test_NMB_imbal/auc_comparision.png", plot = p_grouped_auc, width = 7, height = 5, dpi = 300)


# Testing on Cavalli 

dir.create(path.expand("~/Thesis_models/Cavalli_imbal/"), recursive = TRUE, showWarnings = FALSE)

Cavalli_results<-list()

model_dict<-data.frame(model = character(), model_label = character())

stopifnot("Mismatch in test set and test pheno Sample IDs" = all(rownames(Cavalli_betas) == Cavalli_pheno$Study_ID...12))

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
    
    model<-models[[algo]][[feature_sel]]
    
    test_data<-Cavalli_betas
    
    message("Testing and evaluating ", label, " on Cavalli")
    
    if (algo == "XGB") {
      fn <- model$mod[[1]]$feature_names              
      if (!is.null(fn)) {
        stopifnot("some probes the XGB model was trained on are missing from the test data" =
                    all(fn %in% colnames(Cavalli_betas)))
        test_data <- Cavalli_betas[, fn, drop = FALSE]
      }
    }
    
    res<-tryCatch(
      test_and_evaluate (
        model = model,
        test_data = test_data,
        test_data_name = "Cavalli",
        train_set_name = "Train NMB",
        true_test_labels = Cavalli_pheno$histology,
        model_name = algo_names[[algo]],
        feature_selection_method = feature_label,
        calibrated = FALSE,
        calibration_method = NULL,
        
        file_id = paste(algo, feature_sel, sep = "_"),
        save_in = "~/Thesis_models/Cavalli_imbal/"
      ),
      error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL })
    
    if (is.null(res)) next
    
    res$summary$model_label<-label
    
    res$classwise$model_label<-label
    
    res$summary$algorithm<-algo                          
    
    res$summary$features<-feature_label      
    
    Cavalli_results[[model_key]]<-res
    
    model_dict<-rbind(model_dict, data.frame(model=model_key, model_label = label))
  }
}

stopifnot("No model could be evaluated on Cavalli" = length(Cavalli_results) > 0)

final_Cavalli_summary<-do.call(rbind, lapply(Cavalli_results, function(x) x$summary))

final_Cavalli_classwise_summary<-do.call(rbind, lapply(Cavalli_results, function(x) x$classwise))

print(final_Cavalli_summary)

write.csv(final_Cavalli_summary, "~/Thesis_models/Cavalli_imbal/all_probes_summary_all_models.csv", row.names = FALSE)
write.csv(final_Cavalli_classwise_summary, "~/Thesis_models/Cavalli_imbal/all_probes_classwise_all_models.csv", row.names = FALSE)

# Comparison Bar plot 

p_grouped_accuracy <- ggplot(final_Cavalli_summary, aes(x = algorithm, y = accuracy, fill = features)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(accuracy, 3)),
            position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
  labs(title = wrap_title("Accuracy Comparison by Machine Learning Model and Feature set in Cavalli", width = 60),
       x = "Algorithm", y = "Accuracy score", fill = "Feature Set") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
  theme_minimal()

ggsave("~/Thesis_models/Cavalli_imbal/accuracy_comparision.png", plot = p_grouped_accuracy, width = 7, height = 5, dpi = 300)



p_grouped_auc <- ggplot(final_Cavalli_summary, aes(x = algorithm, y = macro_auc, fill = features)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(macro_auc, 3)),
            position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
  labs(title = wrap_title("Area Under the Curve (AUC) score Comparison by Machine Learning Model and Feature Set in Cavalli", width = 60),
       x = "Algorithm", y = "AUC score", fill = "Feature Set") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
  theme_minimal() 

ggsave("~/Thesis_models/Cavalli_imbal/auc_comparision.png", plot = p_grouped_auc, width = 7, height = 5, dpi = 300)


#=====================================================================================#
#----Balancing the Train set ------
#=====================================================================================#

#Custom 'balance' function for sampling and synthetic sampling methods

balance <- function(train_df, 
                    outcome_labels, 
                    method = c("none","up","down","methylSMOTE",
                               "smote","adasyn","tomek","smote_tomek"),
                    over_ratio = 1, 
                    minority_threshold = NULL, 
                    neighbours = 5, 
                    downsampling = FALSE, 
                    seed = 123,
                    save_in = NULL) {    
  
  # train_df: Train beta value matrix as a data frame; "data.frame"
  
  # outcome_labels: name of the column containing Train class' labels; "character"
  
  # method : "none", "up" - oversampling of the minority classes, "down" - under-sampling of the majority class,
  # "methylSMOTE"- SMOTE with methylClass::balancesampling, "smote" - themis recipe , "adasyn"- Adaptive Synthetic sample generation,
  # "tomek" - Removal of majority class instances in 'Tomek' links, "smote_tomek" - SMOTE with removal of majority class instances in 'Tomek' links ; "character"
  
  # over_ratio : themis based sampling methods ; "numeric"
  
  # minority_threshold (cutoff in methylClass::balancesampling) : the number of samples in the majority class, 
  # all the classes with sample numbers below this will be over sampled with SMOTE ; "numeric"
  # NULL (default) = size of the largest class.
  
  # neighbours: Number of nearest neighbours for SMOTE to interpolate with ; "numeric"
  
  # downsampling: for methylClass::balancesampling, set to TRUE in case the minority_threshold != n in majority class
  # to reduce the sample numbers down to the number set with minority_threshold
  
  # save_in : folder to save the balanced set in, as <method>_train_NMB.rds; NULL = don't save
  
  # ---- Packages -----------------------------------------------------------------------------
  
  required_pgs <- c("caret","themis","recipes")
  
  missing_pgs <- required_pgs[!sapply(required_pgs, requireNamespace, quietly = TRUE)]
  
  if (length(missing_pgs) > 0) {
    message("Installing missing packages for sampling methods:", paste(missing_pgs, collapse = ", "))
    install.packages(missing_pgs, dependencies = TRUE)
  } 
  
  invisible(lapply(required_pgs, library, character.only = TRUE))
  
  method <- match.arg(method)
  
  stopifnot(is.data.frame(train_df), is.factor(train_df[[outcome_labels]]))
  
  set.seed(seed)
  
  x <- train_df[, setdiff(names(train_df), outcome_labels), drop = FALSE]
  y <- train_df[[outcome_labels]]
  
  if (is.null(minority_threshold)) minority_threshold <- max(table(y))
  
  if (method == "none") {
    
    out <- train_df
    
    # ---- Random Oversampling/Under-sampling -----------------------------------------------------------------------------
    # caret::upSample()/downSample() split the table into one wide table per class and glue the
    # pieces back with plyr::rbind.fill, which takes a lot of time
    #  Here the row indexes are chosen first and the rows are then picked from a matrix in one step.
    
  } else if (method %in% c("up", "down")) {
    
    counts <- table(y)
    counts <- counts[counts > 0]                                  
    target <- if (method == "up") max(counts) else min(counts)    
    
    idx <- unlist(lapply(names(counts), function(l) {
      i <- which(y == l)                                          
      if (method == "up") c(i, i[sample.int(length(i), target - length(i), replace = TRUE)])  
      else                i[sample.int(length(i), target)]                                     
    }), use.names = FALSE)
    
    out <- as.data.frame(as.matrix(x)[idx, , drop = FALSE])       
    rownames(out) <- make.unique(rownames(train_df)[idx])       
    out[[outcome_labels]] <- y[idx]
    
    # ---- Synthetic Minority Oversampling (SMOTE) -----------------------------------------------------------------------------
    
  } else if (method == "methylSMOTE") {
    
    out_methylsmote <- methylClass::balancesampling(x, labels = y, 
                                                    topfeaturenumber = NULL,
                                                    cutoff = minority_threshold,
                                                    downsampling = FALSE,
                                                    sampleseed = seed,
                                                    k = neighbours,
                                                    adjustbetas = TRUE) 
    
    if (!is.null(save_in)) {                                      
      dir.create(path.expand(save_in), recursive = TRUE, showWarnings = FALSE)
      saveRDS(out_methylsmote, file.path(path.expand(save_in), "methylSMOTE_raw_output.rds"))
    }
    
    
    if (is.null(out_methylsmote[["dat"]]) || is.null(out_methylsmote[["labels"]])) {
      utils::str(out_methylsmote, max.level = 1)
      stop("balancesampling() did not return elements named 'dat' and 'labels'; see str() above")
    }
    
    out <- as.data.frame(out_methylsmote[["dat"]])
    
    if (nrow(out) != length(out_methylsmote[["labels"]])) {
      stop("rows in 'dat' (", nrow(out), ") do not match 'labels' (", length(out_methylsmote[["labels"]]),
           ")")
    }
    
    out[[outcome_labels]] <- out_methylsmote[["labels"]]
    
    # ---- Themis application -----------------------------------------------------------------------------
    
  } else {
    
    rec <- recipes::recipe(stats::as.formula(paste(outcome_labels, "~ .")), data = train_df)
    
    rec <- switch(method, 
                  smote = rec %>%
                    step_smote(all_outcomes(), over_ratio = over_ratio,
                               neighbors = neighbours, seed = seed),
                  
                  adasyn = rec %>% 
                    step_adasyn(all_outcomes(), over_ratio = over_ratio,
                                neighbors = neighbours, seed = seed),
                  
                  tomek = rec %>%
                    step_tomek(all_outcomes(), seed = seed),
                  
                  smote_tomek = rec %>%
                    step_smote(all_outcomes(), over_ratio = over_ratio,
                               neighbors = neighbours, seed = seed) %>%
                    step_tomek(all_outcomes(), seed = seed)
    )
    
    rec_prep <- prep(rec, training = train_df, retain = TRUE)
    
    out <- as.data.frame(bake(rec_prep, new_data = NULL))
  }
  
  
  stopifnot("a column of the input is missing from the balanced set" = all(names(train_df) %in% names(out)))
  
  out <- out[, names(train_df), drop = FALSE]
  
  out[[outcome_labels]] <- factor(as.character(out[[outcome_labels]]), levels = levels(y))
  
  stopifnot("the balanced set has labels that are not among the input's classes" = !anyNA(out[[outcome_labels]]))
  
  if (!is.null(save_in) && method != "none") {
    
    dir.create(path.expand(save_in), recursive = TRUE, showWarnings = FALSE)
    
    saveRDS(out, file.path(path.expand(save_in), paste0(method, "_train_NMB.rds"))) 
  }
  
  out     
} 


#==========================================================================================#
# ---- Applying balancing methods, training ML models on balanced_sets and testing them----
#==========================================================================================#

library(ggplot2)

# Save_in

file_dir <- "~/Thesis_models/balanced_sets/"                              

dir.create(path.expand(file_dir), recursive = TRUE, showWarnings = FALSE)

balancing_methods<-c("up","smote","adasyn","tomek","smote_tomek") #down,"methylSMOTE"

balanced_sets<-list()

for (method in balancing_methods) {
  
  message("Balancing with ", method, " ...")
  
  balanced_train_NMB<-tryCatch(
    balance (train_df = train_NMB,
             outcome_labels = "histo",
             method = method,
             minority_threshold = max(table(train_NMB$histo)),
             neighbours = 5,
             over_ratio = 1,
             seed = 42 ),
    error = function(e) { message("  FAILED (", method, "): ", conditionMessage(e)); NULL })
  
  if (is.null(balanced_train_NMB)) next
  
  saveRDS(balanced_train_NMB, file.path(path.expand(file_dir), paste0(method,"_train_NMB.rds")))
  
  balanced_sets[[method]]<-balanced_train_NMB
} 

#------ Training ML models on balanced sets -----------

# Load balanced data sets if not in the env.

get_balanced_set <- function(method, dir = "~/Thesis_models/balanced_sets/"){
  readRDS(file.path(path.expand(dir), paste0(method,"_train_NMB.rds"))) }     
  
  
balancing_methods<-c("up","down",
                       "smote","adasyn","tomek","smote_tomek") # "methylSMOTE"

balanced_sets<-list() 
  
for (method in balancing_methods) {
    
    if (file.exists(file.path(path.expand("~/Thesis_models/balanced_sets/"), paste0(method,"_train_NMB.rds")))) { 
      
      balanced_sets[[method]]<-get_balanced_set(method, dir = "~/Thesis_models/balanced_sets/")
      
    } else {
      
      message("No saved set for ", method)
    }
  }
  
  # Training and testing methylClass models on each of the balanced Train NMB set with all probes
  
  balanced_models<-list()
  
  for (set in names(balanced_sets)){
    
    train_set<-balanced_sets[[set]]
    
    x <- train_set[,setdiff(names(train_set), "histo"), drop = FALSE]
    
    y<-as.factor(train_set[["histo"]])
    
    RF_model<-methylClass::maintrain(betas.. = x, 
                                     y.. = y ,
                                     method = "RF",
                                     feature_sel = FALSE,
                                     bal_training = FALSE,
                                     subset.CpGs = NULL,
                                     topfeaturenumber = NULL,
                                     seed = 42,
                                     ntrees = 500, #default
                                     p = 200 )
    
    saveRDS(RF_model, paste0("~/Thesis_models/balanced_sets/","RF_",set,"_all_probes.rds"))
    
    
    # SVM_model <-methylClass::maintrain( betas.. = x,
    #                                     y.. = y,
    #                                     method = "SVM",
    #                                     modelcv = 5, #default
    #                                     gridsearch = TRUE,
    #                                     C.base = 10, #default
    #                                     C.min = -3, #default
    #                                     C.max = -2, #default
    #                                     topfeaturenumber = NULL,
    #                                     subset.CpGs = NULL )
    # 
    # 
    # saveRDS(SVM_model, paste0("~/Thesis_models/balanced_sets/","SVM_",set,"_all_probes.rds"))
    # 
    # #eXtreme Gradient Boosting (XGB)
    # 
    # XGB_model<-methylClass::maintrain( betas.. = as.matrix(x),
    #                                    y.. = y,
    #                                    method = "XGB",
    #                                    modelcv = 5, #default
    #                                    gridsearch = TRUE,
    #                                    max_depth = 6, #default
    #                                    eta = c(0.1,0.3), #default
    #                                    gamma = c(0,0.01), #default
    #                                    min.chwght = 1, #default
    #                                    nrounds = 100, 
    #                                    colsample_bytree = c(0.2,0.5,0.8,1.0), #same as caret XGB models
    #                                    topfeaturenumber = NULL,
    #                                    subset.CpGs = NULL )
    # 
    # 
    # saveRDS(XGB_model, paste0("~/Thesis_models/balanced_sets/","XGB_",set,"_all_probes.rds"))
    
    balanced_models[[set]]<-list(rf = RF_model) #,svm = SVM_model, xgb = XGB_model)
  }
  
  #---Load models if not in the env.
  
  
  algos<-c("RF","SVM","XGB")
  
  load_balanced_models<-function(algo, model_dir){
    files<-list.files(path.expand(model_dir),
                      pattern = paste0("^", algo,"_.*_all_probes\\.rds$"),
                      full.names = TRUE)
    method_names<-gsub(paste0("^", algo, "_|_all_probes.rds$"), "", basename(files))
    setNames(lapply(files,readRDS), method_names)                      
  }
  
  
  balanced_models<-setNames(lapply(algos, load_balanced_models, model_dir = "~/Thesis_models/balanced_sets/"), algos)
  
  
  # ---- Names and probe lists used by every testing section below -----------------------------------------
  
  algo_names<-c(RF = "Random Forest") #SVM = "Support Vector Machine", XGB = "eXtreme Gradient Boosting")
  
  balancing_labels<-c(up = "Up-sampled", down = "Down-sampled", smote = "SMOTE", adasyn = "ADASYN",
                      tomek = "Tomek", smote_tomek = "SMOTE-Tomek") #methylSMOTE = "methylSMOTE
  
  
  all_probes<-setdiff(colnames(train_NMB), "histo")
  #sd_probes <-setdiff(colnames(train_NMB_sd_10k), "histo")
  #mad_probes<-setdiff(colnames(train_NMB_mad_10k), "histo")
  
  # Testing on Test NMB 
  
  dir.create(path.expand("~/Thesis_models/test_NMB_bal/"), recursive = TRUE, showWarnings = FALSE)
  
  test_NMB_results<-list()
  
  model_dict<-data.frame(model = character(), model_label = character())
  
  stopifnot("Mismatch in test set and test pheno Sample IDs" = all(rownames(test_NMB) == NMB_test_pheno$Sample_ID))
  
  for (algo in intersect(names(balanced_models), names(algo_names))){
    for (bal in names(balanced_models[[algo]])){
      
      model_key<-paste(algo, bal, sep = "_")
      label<-paste(algo_names[[algo]], "all probes", balancing_labels[[bal]], sep = " | ")  
      
      
      message("Testing and evaluating ", label, " on Test NMB")
      
      res<-tryCatch(
        test_and_evaluate (
          model = balanced_models[[algo]][[bal]],
          test_data = test_NMB[, all_probes, drop = FALSE],
          test_data_name = "Test NMB",
          true_test_labels = NMB_test_pheno$CPR_Histology,
          model_name = algo_names[[algo]],
          feature_selection_method = "all probes",
          balanced_method = balancing_labels[[bal]],
          file_id = paste(algo, "all", bal, sep = "_"),
          save_in = "~/Thesis_models/test_NMB_bal/"
        ),
        error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL })
      
      if (is.null(res)) next
      
      res$summary$model_label<-label
      
      res$classwise$model_label<-label
      
      res$summary$algorithm<-algo                         
      res$summary$balanced<-balancing_labels[[bal]]        
      
      test_NMB_results[[model_key]]<-res
      
      model_dict<-rbind(model_dict, data.frame(model=model_key, model_label = label))
    }
  }
  
  stopifnot("No model could be evaluated on Test NMB" = length(test_NMB_results) > 0)
  
  final_test_NMB_summary<-do.call(rbind, lapply(test_NMB_results, function(x) x$summary))
  
  final_test_NMB_classwise_summary<-do.call(rbind, lapply(test_NMB_results, function(x) x$classwise))
  
  print(final_test_NMB_summary)
  
  write.csv(final_test_NMB_summary, "~/Thesis_models/test_NMB_bal/all_probes_summary_all_models.csv", row.names = FALSE)
  
  write.csv(final_test_NMB_classwise_summary, "~/Thesis_models/test_NMB_bal/all_probes_classwise_all_models.csv", row.names = FALSE)
  
  # Comparison Bar plot 
  
  p_grouped_accuracy <- ggplot(final_test_NMB_summary, aes(x = balanced, y = accuracy, fill = balanced)) +
    geom_col(position = "dodge", width = 0.6) +
    geom_text(aes(label = round(accuracy, 3)),
              position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
    labs(title = wrap_title("Accuracy Comparison for Random Forest Model by Balancing Method used on Train NMB- Test NMB", width = 60),
         x = "Algorithm", y = "Accuracy score", fill = "Balancing Method") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
    theme_minimal()
  
  ggsave("~/Thesis_models/test_NMB_bal/accuracy_comparision.png", plot = p_grouped_accuracy, width = 7, height = 5, dpi = 300)
  
  
  
  p_grouped_auc <- ggplot(final_test_NMB_summary, aes(x = balanced, y = macro_auc, fill = balanced)) +
    geom_col(position = "dodge", width = 0.6) +
    geom_text(aes(label = round(macro_auc, 3)),
              position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
    labs(title = wrap_title("Area Under the Curve (AUC) score Comparison for Random Forest Model by Balancing Method used on Train NMB- Test NMB", width = 60),
         x = "Algorithm", y = "AUC score", fill = "Balancing Method") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
    theme_minimal() 
  
  ggsave("~/Thesis_models/test_NMB_bal/auc_comparision.png", plot = p_grouped_auc, width = 7, height = 5, dpi = 300)
  
  
  # Testing on Cavalli 
  
  dir.create(path.expand("~/Thesis_models/Cavalli_bal/"), recursive = TRUE, showWarnings = FALSE)
  
  Cavalli_results<-list()
  
  model_dict<-data.frame(model = character(), model_label = character())
  
  stopifnot("Mismatch in test set and test pheno Sample IDs" = all(rownames(Cavalli_betas) == Cavalli_pheno$Sample_ID))
  
  for (algo in intersect(names(balanced_models), names(algo_names))){
    for (bal in names(balanced_models[[algo]])){
      
      model_key<-paste(algo, bal, sep = "_")
      label<-paste(algo_names[[algo]], "all probes", balancing_labels[[bal]], sep = " | ")  
      
      message("Testing and evaluating ", label, " on Cavalli")
      
      
      res<-tryCatch(
        test_and_evaluate (
          model = balanced_models[[algo]][[bal]],
          test_data = Cavalli_betas[, all_probes, drop = FALSE],
          test_data_name = "Cavalli",
          true_test_labels = Cavalli_pheno$histology,
          model_name = algo_names[[algo]],
          feature_selection_method = "all probes",
          balanced_method = balancing_labels[[bal]],
          file_id = paste(algo, "all", bal, sep = "_"),
          save_in = "~/Thesis_models/Cavalli_bal/"
        ),
        error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL })
      
      if (is.null(res)) next
      
      res$summary$model_label<-label
      
      res$classwise$model_label<-label
      
      res$summary$algorithm<-algo                          
      res$summary$balanced<-balancing_labels[[bal]]        
      
      Cavalli_results[[model_key]]<-res
      
      model_dict<-rbind(model_dict, data.frame(model=model_key, model_label = label))
    }
  }
  
  stopifnot("No model could be evaluated on Cavalli" = length(Cavalli_results) > 0)
  
  final_Cavalli_summary<-do.call(rbind, lapply(Cavalli_results, function(x) x$summary))
  
  final_Cavalli_classwise_summary<-do.call(rbind, lapply(Cavalli_results, function(x) x$classwise))
  
  print(final_Cavalli_summary)
  
  write.csv(final_Cavalli_summary, "~/Thesis_models/Cavalli_bal/all_probes_summary_all_models.csv", row.names = FALSE)
  write.csv(final_Cavalli_classwise_summary, "~/Thesis_models/Cavalli_bal/all_probes_classwise_all_models.csv", row.names = FALSE)
  
  # Comparison Bar plot 
  
  p_grouped_accuracy <- ggplot(final_Cavalli_summary, aes(x = balanced, y = accuracy, fill = balanced)) +
    geom_col(position = "dodge", width = 0.6) +
    geom_text(aes(label = round(accuracy, 3)),
              position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
    labs(title = wrap_title("Accuracy Comparison for Random Forest Model by Balancing Method used on Train NMB- Cavalli", width = 60),
         x = "Algorithm", y = "Accuracy score", fill = "Balancing Method") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
    theme_minimal()
  
  ggsave("~/Thesis_models/Cavalli_bal/accuracy_comparision.png", plot = p_grouped_accuracy, width = 7, height = 5, dpi = 300)
  
  
  
  p_grouped_auc <- ggplot(final_Cavalli_summary, aes(x = balanced, y = macro_auc, fill = balanced)) +
    geom_col(position = "dodge", width = 0.6) +
    geom_text(aes(label = round(macro_auc, 3)),
              position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
    labs(title = wrap_title("Area Under the Curve (AUC) score Comparison for Random Forest Model by Balancing Method used on Train NMB- Cavalli", width = 60),
         x = "Algorithm", y = "AUC score", fill = "Balancing Method") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
    theme_minimal() 
  
  ggsave("~/Thesis_models/Cavalli_bal/auc_comparision.png", plot = p_grouped_auc, width = 7, height = 5, dpi = 300)
  
  
  # Training and Testing methylClass models on each of the balanced Train NMB set with top 10k probes (sub setted before applying balancing methods)
  # Using methylClass to allow quick run in comparision to computationally expensive caret models.
  
  # Standard Deviation
  
  sd_cols<-unique(c(colnames(train_NMB_sd_10k), "histo"))     
  
  sd_10k_balanced_sets<-list()
  
  for (set in names(balanced_sets)) {
    
    train_set<-as.data.frame(balanced_sets[[set]])
    
    x<-train_set[, sd_cols, drop = FALSE]
    
    if (all(colnames(x) == sd_cols)) {
      
      saveRDS(x, paste0("~/Thesis_models/balanced_sets/",set,"_10k_sd.rds"))
      
      sd_10k_balanced_sets[[set]]<-x                              
      
    } else {
      
      stop("Error in creating SD 10k balanced sets")
    }
    
  }
  
  # sd_balanced_models<-list()
  # 
  # for (set in names(sd_10k_balanced_sets)){
  #   
  #   train_set<-sd_10k_balanced_sets[[set]]                        
  #   
  #   x <- train_set[,setdiff(names(train_set), "histo"), drop = FALSE]
  #   
  #   y<-as.factor(train_set[["histo"]])
  #   
  #   RF_model<-methylClass::maintrain(betas.. = x, 
  #                                    y.. = y ,
  #                                    method = "RF",
  #                                    subset.CpGs = NULL,
  #                                    topfeaturenumber = NULL,
  #                                    feature_sel = FALSE,
  #                                    bal_training = FALSE,
  #                                    seed = 42,
  #                                    ntrees = 500, #default
  #                                    p = 200 )
  #   
  #   saveRDS(RF_model, paste0("~/Thesis_models/balanced_sets/","RF_",set,"_sd_10k.rds"))
  #   
  #   
  #   SVM_model <-methylClass::maintrain( betas.. = x, 
  #                                       y.. = y,
  #                                       method = "SVM",
  #                                       modelcv = 5, #default
  #                                       gridsearch = TRUE,
  #                                       C.base = 10, #default
  #                                       C.min = -3, #default
  #                                       C.max = -2, #default
  #                                       topfeaturenumber = NULL,
  #                                       subset.CpGs = NULL )
  #   
  #   
  #   saveRDS(SVM_model, paste0("~/Thesis_models/balanced_sets/","SVM_",set,"_sd_10k.rds"))
  #   
  #   #eXtreme Gradient Boosting (XGB)
  #   
  #   XGB_model<-methylClass::maintrain( betas.. = x,
  #                                      y.. = y,
  #                                      method = "XGB",
  #                                      modelcv = 5, #default
  #                                      gridsearch = TRUE,
  #                                      max_depth = 6, #default
  #                                      eta = c(0.1,0.3), #default
  #                                      gamma = c(0,0.01), #default
  #                                      min.chwght = 1, #default
  #                                      nrounds = 100, #default
  #                                      colsample_bytree = c(0.2,0.5,0.8,1.0), #same as caret XGB models
  #                                      topfeaturenumber = NULL,
  #                                      subset.CpGs = NULL )
  #   
  #   
  #   saveRDS(XGB_model, paste0("~/Thesis_models/balanced_sets/","XGB_",set,"_sd_10k.rds"))
  #   
  #   sd_balanced_models[[set]]<-list(rf = RF_model, svm = SVM_model, xgb = XGB_model)
  # }
  # 
  # #--Load the models if not in the env.
  # 
  # load_balanced_models<-function(algo, model_dir){
  #   files<-list.files( path.expand(model_dir),
  #                      pattern = paste0("^", algo,"_.*_sd_10k\\.rds$"),
  #                      full.names = TRUE)
  #   method_names<-gsub(paste0("^", algo, "_|_sd_10k\\.rds$"), "", basename(files))
  #   setNames(lapply(files,readRDS), method_names)                     
  # }
  # 
  # sd_balanced_models<-setNames(lapply(algos, load_balanced_models, model_dir = "~/Thesis_models/balanced_sets/"), algos)  
  # 
  # # Testing on Test NMB 
  # 
  # dir.create(path.expand("~/Thesis_models/test_NMB_bal/"), recursive = TRUE, showWarnings = FALSE)
  # 
  # test_NMB_results<-list()
  # 
  # model_dict<-data.frame(model = character(), model_label = character())
  # 
  # stopifnot("Mismatch in test set and test pheno Sample IDs" = all(rownames(test_NMB) == NMB_test_pheno$Sample_ID))
  # 
  # 
  # for (algo in names(sd_balanced_models)){
  #   for (bal in names(sd_balanced_models[[algo]])){
  #     
  #     model_key<-paste(algo, bal, sep = "_")
  #     label<-paste(algo_names[[algo]], "top 10k by SD", balancing_labels[[bal]], sep = " | ")   
  #     
  #     message("Testing and evaluating ", label, " on Test NMB")
  #     
  #     res<-tryCatch(
  #       test_and_evaluate(
  #         model = sd_balanced_models[[algo]][[bal]],
  #         test_data = test_NMB[, sd_probes, drop = FALSE],
  #         test_data_name = "Test NMB",
  #         true_test_labels = NMB_test_pheno$CPR_Histology,
  #         model_name = algo_names[[algo]],
  #         feature_selection_method = "top 10k by SD",
  #         balancing_method = balancing_labels[[bal]],
  #         file_id = paste(algo, "sd", bal, sep = "_"),
  #         save_in = "~/Thesis_models/test_NMB_bal/"
  #       ),
  #       error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL })
  #     
  #     if (is.null(res)) next
  #     
  #     res$summary$model_label<-label
  #     
  #     res$classwise$model_label<-label
  #     
  #     res$summary$algorithm<-algo                          
  #     res$summary$balanced<-balancing_labels[[bal]]        
  #     
  #     test_NMB_results[[model_key]]<-res
  #     
  #     model_dict<-rbind(model_dict, data.frame(model=model_key, model_label = label))
  #   }
  # }
  # 
  # stopifnot("No model could be evaluated on Test NMB" = length(test_NMB_results) > 0)
  # 
  # final_test_NMB_summary<-do.call(rbind, lapply(test_NMB_results, function(x) x$summary))
  # 
  # final_test_NMB_classwise_summary<-do.call(rbind, lapply(test_NMB_results, function(x) x$classwise))
  # 
  # print(final_test_NMB_summary)
  # 
  # write.csv(final_test_NMB_summary, "~/Thesis_models/test_NMB_bal/sd_10k_summary_all_models.csv", row.names = FALSE)
  # write.csv(final_test_NMB_classwise_summary, "~/Thesis_models/test_NMB_bal/sd_10k_classwise_all_models.csv", row.names = FALSE)
  # 
  # # Comparison Bar plot 
  # 
  # p_grouped_accuracy <- ggplot(final_test_NMB_summary, aes(x = algorithm, y = accuracy, fill = balanced)) +
  #   geom_col(position = "dodge", width = 0.6) +
  #   geom_text(aes(label = round(accuracy, 3)),
  #             position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
  #   labs(title = wrap_title("Accuracy Comparison by Machine Learning Model trained on Top 10k probes (Standard Deviation) and Balancing Method in Test NMB",width = 60),
  #        x = "Algorithm", y = "Accuracy score", fill = "Balancing Method") +
  #   scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
  #   theme_minimal()
  # 
  # ggsave("~/Thesis_models/test_NMB_bal/sd_10k_accuracy_comparision.png", plot = p_grouped_accuracy, width = 7, height = 5, dpi = 300)
  # 
  # 
  # p_grouped_auc <- ggplot(final_test_NMB_summary, aes(x = algorithm, y = macro_auc, fill = balanced)) +
  #   geom_col(position = "dodge", width = 0.6) +
  #   geom_text(aes(label = round(macro_auc, 3)),
  #             position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
  #   labs(title = wrap_title("Area Under the Curve (AUC) score Comparison by Machine Learning Model trained on Top 10k probes (Standard Deviation) and Balancing Method in Test NMB", width = 60),
  #        x = "Algorithm", y = "AUC score", fill = "Balancing Method") +
  #   scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
  #   theme_minimal() 
  # 
  # ggsave("~/Thesis_models/test_NMB_bal/sd_10k_auc_comparision.png", plot = p_grouped_auc, width = 7, height = 5, dpi = 300)
  # 
  # 
  # # Testing on Cavalli 
  # 
  # dir.create(path.expand("~/Thesis_models/Cavalli_bal/"), recursive = TRUE, showWarnings = FALSE)
  # 
  # Cavalli_results<-list()
  # 
  # model_dict<-data.frame(model = character(), model_label = character())
  # 
  # stopifnot("Mismatch in test set and test pheno Sample IDs" = all(rownames(Cavalli_betas) == Cavalli_pheno$Study_ID...12))
  # 
  # for (algo in names(sd_balanced_models)){
  #   for (bal in names(sd_balanced_models[[algo]])){
  #     
  #     model_key<-paste(algo, bal, sep = "_")
  #     label<-paste(algo_names[[algo]], "top 10k by SD", balancing_labels[[bal]], sep = " | ") 
  #     
  #     message("Testing and evaluating ", label, " on Cavalli")
  #     
  #     res<-tryCatch(
  #       test_and_evaluate(
  #         model = sd_balanced_models[[algo]][[bal]],
  #         test_data = Cavalli_betas[, sd_probes, drop = FALSE],
  #         test_data_name = "Cavalli",
  #         true_test_labels = Cavalli_pheno$histology,
  #         model_name = algo_names[[algo]],
  #         feature_selection_method = "top 10k by SD",
  #         balancing_method = balancing_labels[[bal]],
  #         file_id = paste(algo, "sd", bal, sep = "_"),
  #         save_in = "~/Thesis_models/Cavalli_bal/"
  #       ),
  #       error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL })
  #     
  #     if (is.null(res)) next
  #     
  #     res$summary$model_label<-label
  #     
  #     res$classwise$model_label<-label
  #     
  #     res$summary$algorithm<-algo                          
  #     res$summary$balanced<-balancing_labels[[bal]]        
  #     
  #     Cavalli_results[[model_key]]<-res
  #     
  #     model_dict<-rbind(model_dict, data.frame(model=model_key, model_label = label))
  #   }
  # }
  # 
  # stopifnot("No model could be evaluated on Cavalli" = length(Cavalli_results) > 0)
  # 
  # final_Cavalli_summary<-do.call(rbind, lapply(Cavalli_results, function(x) x$summary))
  # 
  # final_Cavalli_classwise_summary<-do.call(rbind, lapply(Cavalli_results, function(x) x$classwise))
  # 
  # print(final_Cavalli_summary)
  # 
  # write.csv(final_Cavalli_summary, "~/Thesis_models/Cavalli_bal/sd_10k_summary_all_models.csv", row.names = FALSE)
  # write.csv(final_Cavalli_classwise_summary, "~/Thesis_models/Cavalli_bal/sd_10k_classwise_all_models.csv", row.names = FALSE)
  # 
  # # Comparison Bar plot 
  # 
  # p_grouped_accuracy <- ggplot(final_Cavalli_summary, aes(x = algorithm, y = accuracy, fill = balanced)) +
  #   geom_col(position = "dodge", width = 0.6) +
  #   geom_text(aes(label = round(accuracy, 3)),
  #             position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
  #   labs(title = wrap_title("Accuracy Comparison by Machine Learning Model trained on Top 10k probes (Standard Deviation) and Balancing Method in Cavalli", width = 60),
  #        x = "Algorithm", y = "Accuracy score", fill = "Balancing Method") +
  #   scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
  #   theme_minimal()
  # 
  # ggsave("~/Thesis_models/Cavalli_bal/sd_10k_accuracy_comparision.png", plot = p_grouped_accuracy, width = 7, height = 5, dpi = 300)
  # 
  # 
  # p_grouped_auc <- ggplot(final_Cavalli_summary, aes(x = algorithm, y = macro_auc, fill = balanced)) +
  #   geom_col(position = "dodge", width = 0.6) +
  #   geom_text(aes(label = round(macro_auc, 3)),
  #             position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
  #   labs(title = wrap_title("Area Under the Curve (AUC) score Comparison by Machine Learning Model trained on Top 10k probes (Standard Deviation) and Balancing Method in Cavalli", width = 60),
  #        x = "Algorithm", y = "AUC score", fill = "Balancing Method") +
  #   scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
  #   theme_minimal() 
  # 
  # ggsave("~/Thesis_models/Cavalli_bal/sd_10k_auc_comparision.png", plot = p_grouped_auc, width = 7, height = 5, dpi = 300)
  # 
  # 
  # Median Absolute Deviation
  
  mad_cols<-unique(c(colnames(train_NMB_mad_10k), "histo"))    
  
  mad_10k_balanced_sets<-list()
  
  for (set in names(balanced_sets)) {
    
    train_set<-as.data.frame(balanced_sets[[set]])
    
    x<-train_set[, mad_cols, drop = FALSE]
    
    if (all(colnames(x) == mad_cols)) {
      saveRDS(x, paste0("~/Thesis_models/balanced_sets/",set,"_10k_mad.rds"))
      mad_10k_balanced_sets[[set]]<-x                              
    } else {
      stop("Error in creating MAD 10k balanced sets")
    }
    
  }
  
  mad_balanced_models<-list()
  
  for (set in names(mad_10k_balanced_sets)){
    
    train_set<-mad_10k_balanced_sets[[set]]                       
    
    x <- train_set[,setdiff(names(train_set), "histo"), drop = FALSE]
    
    y<-as.factor(train_set[["histo"]])
    
    RF_model<-methylClass::maintrain(betas.. = x, y.. = y ,
                                     method = "RF",
                                     subset.CpGs = NULL,
                                     topfeaturenumber = NULL,
                                     feature_sel = FALSE,
                                     bal_training = FALSE,
                                     seed = 42,
                                     ntrees = 500, #default
                                     p = 200 )
    
    saveRDS(RF_model, paste0("~/Thesis_models/balanced_sets/","RF_",set,"_mad_10k.rds"))
    
    
    # SVM_model <-methylClass::maintrain( betas.. = x,
    #                                     y.. = y,
    #                                     method = "SVM",
    #                                     modelcv = 5, #default
    #                                     gridsearch = TRUE,
    #                                     C.base = 10, #default
    #                                     C.min = -3, #default
    #                                     C.max = -2, #default
    #                                     topfeaturenumber = NULL,
    #                                     subset.CpGs = NULL )
    # 
    # 
    # saveRDS(SVM_model, paste0("~/Thesis_models/balanced_sets/","SVM_",set,"_mad_10k.rds"))
    # 
    # #eXtreme Gradient Boosting (XGB)
    # 
    # XGB_model<-methylClass::maintrain( betas.. = x,
    #                                    y.. = y,
    #                                    method = "XGB",
    #                                    modelcv = 5, #default
    #                                    gridsearch = TRUE,
    #                                    max_depth = 6, #default
    #                                    eta = c(0.1,0.3), #default
    #                                    gamma = c(0,0.01), #default
    #                                    min.chwght = 1, #default
    #                                    nrounds = 100, #default
    #                                    colsample_bytree = c(0.2,0.5,0.8,1.0), #same as caret XGB models
    #                                    topfeaturenumber = NULL,
    #                                    subset.CpGs = NULL )
    # 
    # 
    # saveRDS(XGB_model, paste0("~/Thesis_models/balanced_sets/","XGB_",set,"_mad_10k.rds"))
    
    mad_balanced_models[[set]]<-list(rf = RF_model) #,svm = SVM_model, xgb = XGB_model)
  }
  
  #--Load the models if not in the env.
  
  load_balanced_models<-function(algo, model_dir){
    files<-list.files( path.expand(model_dir),
                       pattern = paste0("^", algo,"_.*_mad_10k\\.rds$"),
                       full.names = TRUE)
    method_names<-gsub(paste0("^", algo, "_|_mad_10k\\.rds$"), "", basename(files))
    setNames(lapply(files,readRDS), method_names)                     
  }
  
  mad_balanced_models<-setNames(lapply(algos, load_balanced_models, model_dir = "~/Thesis_models/balanced_sets/"), algos)  
  
  # Testing on Test NMB 
  
  dir.create(path.expand("~/Thesis_models/test_NMB_bal/"), recursive = TRUE, showWarnings = FALSE)
  
  test_NMB_results<-list()
  
  model_dict<-data.frame(model = character(), model_label = character())
  
  stopifnot("Mismatch in test set and test pheno Sample IDs" = all(rownames(test_NMB) == NMB_test_pheno$Sample_ID))
  
  for (algo in names(mad_balanced_models)){
    for (bal in names(mad_balanced_models[[algo]])){
      
      model_key<-paste(algo, bal, sep = "_")
      label<-paste(algo_names[[algo]], "top 10k by MAD", balancing_labels[[bal]], sep = " | ")   
      
      message("Testing and evaluating ", label, " on Test NMB")
      
      res<-tryCatch(
        test_and_evaluate (
          model = mad_balanced_models[[algo]][[bal]],
          test_data = test_NMB[, mad_probes, drop = FALSE],
          test_data_name = "Test NMB",
          true_test_labels = NMB_test_pheno$CPR_Histology,
          model_name = algo_names[[algo]],
          feature_selection_method = "top 10k by MAD",
          balancing_method = balancing_labels[[bal]],
          file_id = paste(algo, "mad", bal, sep = "_"),
          save_in = "~/Thesis_models/test_NMB_bal/"
        ),
        error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL })
      
      if (is.null(res)) next
      
      res$summary$model_label<-label
      
      res$classwise$model_label<-label
      
      res$summary$algorithm<-algo                          
      res$summary$balanced<-balancing_labels[[bal]]        
      
      test_NMB_results[[model_key]]<-res
      
      model_dict<-rbind(model_dict, data.frame(model=model_key, model_label = label))
    }
  }
  
  stopifnot("No model could be evaluated on Test NMB" = length(test_NMB_results) > 0)
  
  final_test_NMB_summary<-do.call(rbind, lapply(test_NMB_results, function(x) x$summary))
  
  final_test_NMB_classwise_summary<-do.call(rbind, lapply(test_NMB_results, function(x) x$classwise))
  
  print(final_test_NMB_summary)
  
  write.csv(final_test_NMB_summary, "~/Thesis_models/test_NMB_bal/mad_10k_summary_all_models.csv", row.names = FALSE)
  write.csv(final_test_NMB_classwise_summary, "~/Thesis_models/test_NMB_bal/mad_10k_classwise_all_models.csv", row.names = FALSE)
  
  # Comparison Bar plot 
  
  p_grouped_accuracy <- ggplot(final_test_NMB_summary, aes(x = balanced, y = accuracy, fill = balanced)) +
    geom_col(position = "dodge", width = 0.6) +
    geom_text(aes(label = round(accuracy, 3)),
              position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
    labs(title = wrap_title("Accuracy Comparison by Random Forest Model trained on Top 10k probes (Median Absolute Deviation) and Balancing Method in Test NMB", width = 60),
         x = "Algorithm", y = "Accuracy score", fill = "Balancing Method") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
    theme_minimal()
  
  ggsave("~/Thesis_models/test_NMB_bal/mad_10k_accuracy_comparision.png", plot = p_grouped_accuracy, width = 7, height = 5, dpi = 300)
  
  
  
  p_grouped_auc <- ggplot(final_test_NMB_summary, aes(x = balanced, y = macro_auc, fill = balanced)) +
    geom_col(position = "dodge", width = 0.6) +
    geom_text(aes(label = round(macro_auc, 3)),
              position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
    labs(title = wrap_title("Area Under the Curve (AUC) score Comparison by Random Forest Model trained on Top 10k probes (Median Absolute Deviation) and Balancing Method in Test NMB", width = 60),
         x = "Algorithm", y = "AUC score", fill = "Balancing Method") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
    theme_minimal() 
  
  ggsave("~/Thesis_models/test_NMB_bal/mad_10k_auc_comparision.png", plot = p_grouped_auc, width = 7, height = 5, dpi = 300)
  
  
  # Testing on Cavalli 
  
  dir.create(path.expand("~/Thesis_models/Cavalli_bal/"), recursive = TRUE, showWarnings = FALSE)
  
  Cavalli_results<-list()
  
  model_dict<-data.frame(model = character(), model_label = character())
  
  stopifnot("Mismatch in test set and test pheno Sample IDs" = all(rownames(Cavalli_betas) == Cavalli_pheno$Study_ID...12))
  
  for (algo in names(mad_balanced_models)){
    for (bal in names(mad_balanced_models[[algo]])){
      
      model_key<-paste(algo, bal, sep = "_")
      label<-paste(algo_names[[algo]], "top 10k by MAD", balancing_labels[[bal]], sep = " | ")
      
      message("Testing and evaluating ", label, " on Cavalli")
      
      res<-tryCatch(
        test_and_evaluate (
          model = mad_balanced_models[[algo]][[bal]],
          test_data = Cavalli_betas[, mad_probes, drop = FALSE],
          test_data_name = "Cavalli",
          true_test_labels = Cavalli_pheno$histology,
          model_name = algo_names[[algo]],
          feature_selection_method = "top 10k by MAD",
          balancing_method = balancing_labels[[bal]],
          file_id = paste(algo, "mad", bal, sep = "_"),
          save_in = "~/Thesis_models/Cavalli_bal/"
        ),
        error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL })
      
      if (is.null(res)) next
      
      res$summary$model_label<-label
      
      res$classwise$model_label<-label
      
      res$summary$algorithm<-algo                          
      res$summary$balanced<-balancing_labels[[bal]]        
      
      Cavalli_results[[model_key]]<-res
      
      model_dict<-rbind(model_dict, data.frame(model=model_key, model_label = label))
    }
  }
  
  stopifnot("No model could be evaluated on Cavalli" = length(Cavalli_results) > 0)
  
  final_Cavalli_summary<-do.call(rbind, lapply(Cavalli_results, function(x) x$summary))
  
  final_Cavalli_classwise_summary<-do.call(rbind, lapply(Cavalli_results, function(x) x$classwise))
  
  print(final_Cavalli_summary)
  
  write.csv(final_Cavalli_summary, "~/Thesis_models/Cavalli_bal/mad_10k_summary_all_models.csv", row.names = FALSE)
  write.csv(final_Cavalli_classwise_summary, "~/Thesis_models/Cavalli_bal/mad_10k_classwise_all_models.csv", row.names = FALSE)
  
  # Comparison Bar plot 
  
  p_grouped_accuracy <- ggplot(final_Cavalli_summary, aes(x = balanced, y = accuracy, fill = balanced)) +
    geom_col(position = "dodge", width = 0.6) +
    geom_text(aes(label = round(accuracy, 3)),
              position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
    labs(title = wrap_title("Accuracy Comparison by Random Forest Model trained on Top 10k probes (Median Absolute Deviation) and Balancing Method in Cavalli", width = 60),
         x = "Algorithm", y = "Accuracy score", fill = "Balancing Method") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
    theme_minimal()
  
  ggsave("~/Thesis_models/Cavalli_bal/mad_10k_accuracy_comparision.png", plot = p_grouped_accuracy, width = 7, height = 5, dpi = 300)
  
  
  
  p_grouped_auc <- ggplot(final_Cavalli_summary, aes(x = balanced, y = macro_auc, fill = balanced)) +
    geom_col(position = "dodge", width = 0.6) +
    geom_text(aes(label = round(macro_auc, 3)),
              position = position_dodge(width = 0.6), angle = 90, hjust = -0.15, size = 2.6) +
    labs(title = wrap_title("Area Under the Curve (AUC) score Comparison by Random Forest Model trained on Top 10k probes (Median Absolute Deviation) and Balancing Method in Cavalli", width = 60),
         x = "Algorithm", y = "AUC score", fill = "Balancing Method") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.2))) +
    theme_minimal() 
  
  ggsave("~/Thesis_models/Cavalli_bal/mad_10k_auc_comparision.png", plot = p_grouped_auc, width = 7, height = 5, dpi = 300)
  
  
  #----PCA plots to show default clsutering of MB's methylation groups and where synthetic samples generated by SMOTE sit -----
  
  require(ggplot2)
  
  bal_train_NMB<-readRDS("~/Thesis_models/balanced_sets/smote_train_NMB.rds") 
  
  
  #PCA for train NMB
  
  stopifnot("NAs in train_NMB: impute or remove them before PCA" = !anyNA(train_NMB))
  
  pca_dir<-"~/Thesis_models/PCA_plots/"
  
  dir.create(path.expand(pca_dir), recursive = TRUE, showWarnings = FALSE)
  
  #train_NMB<-train_NMB[, matrixStats::colSds(as.matrix(train_NMB)) > 0, drop = FALSE]
  
  NMB_pca<-prcomp(train_NMB,center = FALSE, scale. = FALSE)
  
  saveRDS(NMB_pca, file.path(path.expand(pca_dir), "NMB_pca.rds")) 
  
  if (all(rownames(NMB_pca$x) == NMB_train_pheno$Sample_ID)){
    NMB_pca_df<-data.frame("PC1" = NMB_pca$x[,1],
                           "PC2" = NMB_pca$x[,2],
                           "Histology" = NMB_train_pheno$CPR_Histology,
                           "Molecular_group" = NMB_train_pheno$Mol_group)
    
  } else {
    NMB_train_pheno<-NMB_train_pheno[match(rownames(NMB_pca$x), NMB_train_pheno$Sample_ID),]
    NMB_pca_df<-data.frame("PC1" = NMB_pca$x[,1],
                           "PC2" = NMB_pca$x[,2],
                           "Histology" = NMB_train_pheno$CPR_Histology,
                           "Molecular_group" = NMB_train_pheno$Mol_group)   
    
  }
  
  saveRDS(NMB_pca_df, file.path(path.expand(pca_dir), "NMB_pca_df.rds"))
  
  # PCA plot
  
  NMB_PCA_plot<-ggplot(NMB_pca_df,aes(x = PC1, y = PC2, fill = Molecular_group, shape = Histology)) +
    geom_point(size = 3, alpha = 0.7, stroke = 1 , colour = "black") + 
    labs(title = paste(wrap_title(paste0("Patient Samples in Train_NMB in Methylation Feature Space - PCA plot"), width = 60), collapse = "\n"), x = "PC1", y = "PC2",
         subtitle = "PC1 & PC2:Principal Components along which the samples vary the most", caption = wrap_title(" Clusters are mostly homogenous with MB molecular groups and histology types are not exclusive to a particular molecular group,resulting in heterogenous clusters for Histoloogy types of MB"),fill = "Principle Molecular Group", shape = "Histology_Type") + 
    scale_fill_manual(values = c("SHH-MB" = "red", "WNT-MB" =  "blue","Group4-MB" = "darkgreen", "Group3-MB" = "yellow")) +
    
    scale_shape_manual(values = c("CLA" = 21, "DN" = 24 , "LCA" = 22,"MBEN" = 23)) +
    
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 4)),
           shape = guide_legend(override.aes = list(size = 4))) +
    theme_minimal (base_size = 12) +
    
    theme (plot.title = element_text (face = "bold", size = 14,hjust = 0.5),
           legend.position = "bottom",
           panel.grid.minor = element_blank(),
           plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30", 
                                        margin = ggplot2::margin(b = 10)),
           plot.caption = element_text(size = 9, hjust = 0.5, color = "gray40",
                                       face = "italic", margin = ggplot2::margin(t = 10))) 
  
  ggsave(file.path(path.expand(pca_dir), "PCA_Train_NMB_mol_grp.png"), NMB_PCA_plot,
         width = 15, height = 6.5, dpi = 300, bg = "white")
  
  # Projecting synthetic samples generated by SMOTE onto methylation feature space with PCA fit generated above
  
  if (!"histo" %in% colnames(bal_train_NMB)){
    stop("No histo column in the balanced set")
  } else {
    bal_train_NMB$histo<-as.factor(bal_train_NMB$histo)
  }
  
  bal_train_histo<-bal_train_NMB$histo
  
  bal_train_NMB<-bal_train_NMB[,!colnames(bal_train_NMB) %in% "histo"] 
  
  bal_train_pca<-predict(NMB_pca,newdata = bal_train_NMB)
  
  # To label synthetic samples apart from actual samples, we are using digest 
  # to create hash representations to compare/match
  
  HISTO_COLOURS<-c("CLA" = "purple","LCA" = "pink","DN" = "darkgreen", "MBEN" = "lightgreen")
  
  require(digest)
  
  pca_cols<-rownames(NMB_pca$rotation)
  
  var_exp <- round(100 * NMB_pca$sdev^2 / sum(NMB_pca$sdev^2), 1)
  
  og_hash<-apply(as.matrix(train_NMB[, pca_cols]), 1, function(r) digest::digest(unname(round(r, 8))))
  
  bal_hash<-apply(as.matrix(bal_train_NMB[, pca_cols]), 1, function(r) digest::digest(unname(round(r, 8))))
  
  syn_type<-ifelse(bal_hash %in% og_hash, "Given", "Synthetic")
  
  bal_pca_df <- data.frame(PC1 = bal_train_pca[, 1],
                           PC2 = bal_train_pca[, 2],
                           Histology = bal_train_histo,
                           Sample_type = factor(syn_type, levels = c("Given", "Synthetic")))
  
  bal_pca_df <- bal_pca_df[order(bal_pca_df$Sample_type), ]
  
  syn_PCA_plot <- ggplot(bal_pca_df, aes(x = PC1, y = PC2, colour = Histology, shape = Sample_type)) +
    geom_point(size = 2.5, alpha = 0.7) +
    scale_colour_manual(values = HISTO_COLOURS) +
    scale_shape_manual(values = c("Given" = 16, "Synthetic" = 17)) +
    labs(title = paste(wrap_title(paste0("Given patient samples and synthetic samples (generated by SMOTE) of Train NMB in methylation feature space - PCA"), width = 60), collapse = "\n"),
         subtitle = "PCA fitted on the given samples only; the synthetic samples are projected onto it",
         x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"), 
         colour = "Histology", shape = "Sample") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "gray30"),
          legend.position = "bottom", panel.grid.minor = element_blank())
  
  ggsave(file.path(path.expand(pca_dir), "PCA_given_vs_synthetic.png"), syn_PCA_plot,
         width = 8, height = 6.5, dpi = 300, bg = "white")
  
#----- Optimising the Train Feature Set ------

# Checking which pair of PCA components allow for maximal visual separation of the histology groups 
  
library(ggplot2)

if(!requireNamespace(GGally, quitely = TRUE)) install.packages("GGally")
  
pca_dir<-"~/Thesis_models/PCA_plots/"

NMB_pca<-readRDS(file.path(path.expand(pca_dir), "NMB_pca.rds"))

stopifnot(identical(rownames(NMB_pca$x), as.character(NMB_train_pheno$Sample_ID)))

pca20_df <- data.frame(NMB_pca$x[, 1:20],
                       Histology = factor(NMB_train_pheno$CPR_Histology, levels = c("CLA", "DN", "LCA", "MBEN")))
pca20_df <- pca20_df[order(pca20_df$Histology), ]

HISTO_COLOURS<-c("CLA" = "purple","LCA" = "pink","DN" = "darkgreen", "MBEN" = "lightgreen")
  

g_20 <- ggpairs(pca20_df, columns = 1:20,
             aes(colour = Histology, fill = Histology),
             upper = list(continuous = "blank"),
             lower = list(continuous = wrap("points", size = 0.4, alpha = 0.6)),
             diag  = list(continuous = wrap("densityDiag", alpha = 0.4)),
             legend = c(2, 1)) +
  scale_colour_manual(values = HISTO_COLOURS) +
  scale_fill_manual(values = HISTO_COLOURS) +
  theme_bw(base_size = 6) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), legend.position = "bottom")

png(file.path(path.expand(pca_dir), "PCA_PC1_PC20_ggpairs.png"), width = 3600, height = 3600, res = 200)
print(g_20)
dev.off() 


# -------- LASSO Feature Selection with glmnet ---------

# ------ Load Packages and Functions ----

source("~/Custom_functions.R")

lasso_dir<-"~/Thesis_models/LASSO/"

dir.create(path.expand(lasso_dir), recursive = TRUE, showWarnings = FALSE)

required_pgs <- c("caret", "glmnet", "MLmetrics", "pROC", "ggplot2")

sapply(required_pgs, library, character.only = TRUE, quietly = TRUE)

missing_pgs <- required_pgs[!sapply(required_pgs, requireNamespace, quietly = TRUE)]

if (length(missing_pgs) > 0) install.packages(missing_pgs, dependencies = TRUE)

HISTO_LEVELS<-c("CLA","LCA","DN")

NMB_train_pheno <- readRDS("~/NMB/NMB_train_pheno.rds")

# DN and MBEN due to their similarities and frequent occurence in SHH-pathway activated tumours, are combined into one class as "DN"

NMB_train_pheno$CPR_Histology<-gsub("MBEN","DN", NMB_train_pheno$CPR_Histology)

NMB_train_pheno <- NMB_train_pheno[!is.na(NMB_train_pheno$Sample_ID), ]

NMB_train_pheno$CPR_Histology <- factor(NMB_train_pheno$CPR_Histology, levels = HISTO_LEVELS)

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

#----- Beta Matrices ----

NMB_betas<-readRDS("~/NMB/Final_NMB_betas.rds")

train_NMB<-NMB_betas[match(NMB_train_pheno$Sample_ID, rownames(NMB_betas)),]

test_NMB<-NMB_betas[match(NMB_test_pheno$Sample_ID, rownames(NMB_betas)),]

Cavalli_betas<-readRDS("~/Cavalli/Cavalli-test/Cavalli_test_betas.rds") 

Cavalli_betas<-as.data.frame(t(Cavalli_betas)) 

Cavalli_betas<-Cavalli_betas[match(Cavalli_pheno$Sample_ID, rownames(Cavalli_betas)),] 

# ----- GLMNET to optimise regularisation factor   -------

set.seed(42)

x_train <- as.matrix(train_NMB)

if(all(rownames(x_train) == NMB_train_pheno$Sample_ID)) { y_train<- as.factor(NMB_train_pheno$CPR_Histology) 

} else {
  
  NMB_train_pheno<-NMB_train_pheno[match(rownames(x_train), NMB_train_pheno$Sample_ID),]
  
  y_train<- as.factor(NMB_train_pheno$CPR_Histology) }

cv.lasso <- cv.glmnet(
  x = x_train,
  y = y_train,
  nfolds = 10, #default
  family = "multinomial",
  alpha = 1,
  type.measure = "class"
)

selected_features <- unique(unlist(
  lapply(coef(cv.lasso, s = "lambda.min"), function(m) {
    rownames(m)[which(m[, 1] != 0)]
  })
))

selected_features <- setdiff(selected_features, "(Intercept)") 

saveRDS(selected_features,file.path(path.expand(lasso_dir), "Train_NMB_Lasso_features.rds"))

cat("Number of features selected by LASSO:", length(selected_features), "\n")

train_df <- data.frame(x_train[, selected_features, drop = FALSE])


# ----- Training and testing an RF model on LASSO selected features --------

RF_LASSO<-maintrain(betas.. = train_df,
                    y.. = if(all(rownames(train_df) == NMB_train_pheno$Sample_ID)) NMB_train_pheno$CPR_Histology,
                    seed = 42,
                    method = "RF",
                    ntrees = 100, 
                    p = "89",
                    feature_sel = FALSE,
                    bal_feature_sel = FALSE,
                    bal_training = FALSE,
                    stratified_bal_group = NULL,
                    topfeaturenumber = NULL, 
                    imp_replace = TRUE,
                    subset.CpGs = NULL,
                    calibrationmethod = c("MR","FLR","LR") 
)

saveRDS(RF_LASSO,file.path(path.expand(lasso_dir), "RF_Train_NMB_LASSO.rds")) 

RF_LASSO<-readRDS("~/Thesis_models/LASSO/RF_Train_NMB_LASSO.rds")

message("Testing and evaluating on Test NMB")

if (all(rownames(test_NMB) == NMB_test_pheno$Sample_ID)) {
  test_NMB_results_lasso <- tryCatch(
    test_and_evaluate (
      model = RF_LASSO,
      method = "methylClass",
      calibrated = FALSE,
      calibration_method = NULL,
      train_set_name = "Train NMB",
      test_data = test_NMB,
      test_data_name = "Test NMB",
      true_test_labels = as.factor(NMB_test_pheno$CPR_Histology),
      model_name = "LASSO Random Forest",
      save_in = file.path(path.expand(lasso_dir))
    ),
    error = function(e) {
      message(sprintf("FAILED [%s / %s]: %s",
                      "LASSO Random Forest trained on Train NMB", "Test NMB", conditionMessage(e)))
      NULL
    }
  )
} else {
  stop("Mismatch in test set and test pheno Sample IDs")
}

Test_NMB_heatmap_lasso <- prob_heatmaps(
  prob_df           = test_NMB_results_lasso$pred_probs,
  pred              = test_NMB_results_lasso$pred_probs$pred,
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
  column_for_panel  = "LCA",
  continuum         = "LCA",
  heatmap_label     = "Random Forest Model for classification of histopathology,\ntrained on LASSO selected features in 'Train NMB',and tested on 'Test NMB'",
  heatmap_panel     = FALSE,
  file_id           = "RF_lasso_test_NMB",
  save_in           = file.path(path.expand(lasso_dir))
)


message("Testing and evaluating on Cavalli")

if (all(rownames(Cavalli_betas) == Cavalli_pheno$Sample_ID)) {
  test_Cavalli_results_lasso <- tryCatch(
    test_and_evaluate(
      model = RF_LASSO,
      method = "methylClass",
      calibrated = FALSE,
      calibration_method = NULL,
      train_set_name = "Train NMB",
      test_data = Cavalli_betas,
      test_data_name = "Cavalli",
      true_test_labels = as.factor(Cavalli_pheno$histology),
      model_name = "LASSO Random Forest",
      save_in = file.path(path.expand(lasso_dir)) 
    ),
    error = function(e) {
      message(sprintf("FAILED [%s / %s]: %s",
                      "LASSO Random Forest trained on Train NMB", "Cavalli", conditionMessage(e)))
      NULL
    }
  )
} else {
  stop("Mismatch in test set and test pheno Sample IDs")
} 


Cavalli_heatmap_lasso <- prob_heatmaps(
  prob_df           = test_Cavalli_results_lasso$pred_probs,
  pred              = test_Cavalli_results_lasso$pred_probs$pred,
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
  column_for_panel  = "LCA",
  continuum         = "LCA",
  heatmap_label     = "Random Forest Model for classification of histopathology,\ntrained on LASSO selected features in 'Train NMB',and tested on 'Cavalli'",
  heatmap_panel     = FALSE,
  file_id           = "RF_lasso_Cavalli",
  save_in           = file.path(path.expand(lasso_dir))
)


# --------- Treating the composition bias of molecular group as a batch effect -------

# ----- Load Data if not in the env.-----

devtools::document("~/methylClass-fork") 
devtools::load_all("~/methylClass-fork")

HISTO_LEVELS<-c("CLA","LCA","DN")

NMB_train_pheno <- readRDS("~/NMB/NMB_train_pheno.rds")

NMB_train_pheno$CPR_Histology<-gsub("MBEN","DN", NMB_train_pheno$CPR_Histology)

NMB_train_pheno <- NMB_train_pheno[!is.na(NMB_train_pheno$Sample_ID), ]

NMB_train_pheno$CPR_Histology <- factor(NMB_train_pheno$CPR_Histology, levels = HISTO_LEVELS)

# ---- Test Sets ----

NMB_test_pheno <- readRDS("~/NMB/NMB_test_pheno.rds")

NMB_test_pheno$CPR_Histology<-gsub("MBEN","DN", NMB_test_pheno$CPR_Histology)

NMB_test_pheno <- NMB_test_pheno[!is.na(NMB_test_pheno$Sample_ID), ]

NMB_test_pheno$CPR_Histology<- factor(NMB_test_pheno$CPR_Histology, levels = HISTO_LEVELS)

Cavalli_pheno <- readRDS("~/Cavalli/Final_Cavalli_pheno.rds")

Cavalli_pheno$histology<-gsub("MBEN","DN",Cavalli_pheno$histology)

Cavalli_pheno<- Cavalli_pheno[!is.na(Cavalli_pheno$Sample_ID), ]

Cavalli_pheno$histology <- factor(Cavalli_pheno$histology, levels = HISTO_LEVELS)

#----- Beta Matrices ----

NMB_betas<-readRDS("~/NMB/Final_NMB_betas.rds")

train_NMB<-NMB_betas[match(NMB_train_pheno$Sample_ID, rownames(NMB_betas)),]

test_NMB<-NMB_betas[match(NMB_test_pheno$Sample_ID, rownames(NMB_betas)),]

Cavalli_betas<-readRDS("~/Cavalli/Cavalli-test/Cavalli_test_betas.rds") 

Cavalli_betas<-as.data.frame(t(Cavalli_betas)) 

Cavalli_betas<-Cavalli_betas[match(Cavalli_pheno$Sample_ID, rownames(Cavalli_betas)),] 

# Using M-values for batch correction 

batch_corr_dir<-"~/Thesis_models/mol_batch_corr"

dir.create(path.expand(batch_corr_dir), recursive = TRUE, showWarnings = FALSE)

library(limma)

library(minfi)

NMB_Mval<-readRDS("~/NMB/NMB_mval.rds")

NMB_Mval<-NMB_Mval[match(colnames(test_NMB), rownames(NMB_Mval)),] 

train_NMB_Mval <- NMB_Mval[, NMB_train_pheno$Sample_ID, drop = FALSE]

stopifnot(all(rownames(as.data.frame(t(train_NMB_Mval))) == NMB_train_pheno$Sample_ID))

adj_train_NMB <- removeBatchEffect(train_NMB_Mval,
batch  = NMB_train_pheno$Mol_group,
design = model.matrix(~ CPR_Histology, data = NMB_train_pheno))


train_NMB_adj <- as.data.frame(t(minfi::ilogit2(adj_train_NMB)))

if(all(rownames(train_NMB_adj) == NMB_train_pheno$Sample_ID)) {
  train_NMB_adj$histo<-as.factor(NMB_train_pheno$CPR_Histology) 
} else {
  train_NMB_adj<-train_NMB_adj[match(NMB_train_pheno$Sample_ID, rownames(train_NMB_adj)),]
  
  train_NMB_adj$histo<-as.factor(NMB_train_pheno$CPR_Histology) 
  
  }

saveRDS(train_NMB_adj,file.path(path.expand(batch_corr_dir), "Train_NMB_mol_adj.rds"))


RF_batch_corr<-maintrain(betas.. = train_NMB_adj[,!colnames(train_NMB_adj) %in% "histo"],
                    y.. = train_NMB_adj$histo,
                    seed = 42,
                    method = "RF",
                    ntrees = 100, 
                    p = "200", #default
                    feature_sel = FALSE,
                    bal_feature_sel = FALSE,
                    bal_training = FALSE,
                    stratified_bal_group = NULL,
                    topfeaturenumber = NULL, 
                    imp_replace = TRUE,
                    subset.CpGs = NULL,
                    calibrationmethod = c("MR","FLR","LR") 
)

saveRDS(RF_batch_corr,file.path(path.expand(batch_corr_dir), "RF_mol_batch_corr.rds")) 

# Testing on NMB

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
      test_data = test_NMB,
      test_data_name = "Test NMB",
      true_test_labels = as.factor(NMB_test_pheno$CPR_Histology),
      model_name = "Random Forest",
      save_in = file.path(path.expand(batch_corr_dir))
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
      test_data = Cavalli_betas,
      test_data_name = "Cavalli",
      true_test_labels = as.factor(Cavalli_pheno$histology),
      model_name = "Random Forest",
      save_in = file.path(path.expand(batch_corr_dir))
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



# ----------- Stratifying NMB Cohort into FS and non-FS to allow an other CPR cohort for testing (Feasibility Study samples from Crosier et.al, 2021) ---------

# ----- Load datasets ----

FS_dir<-"~/Thesis_models/FS/" 

dir.create(path.expand(FS_dir), recursive = TRUE, showWarnings = FALSE)

# The methylation array sets representing both of the sets have been pre-processed separately to avoid data leakage 

# Train set - non_FS

NMB_betas<-readRDS("~/NMB/Final_NMB_betas.rds")

non_FS_betas<-as.data.frame(t(readRDS("~/FS/non_FS/non_FSnon_FS_meth_betas.rds")))

non_FS_betas<-non_FS_betas[,match(colnames(NMB_betas),colnames(non_FS_betas))]

NMB_pheno<-readRDS("~/NMB/NMB_Final_Cohort_26.rds")

NMB_pheno$CPR_Histology<-gsub("MBEN","DN",NMB_pheno$CPR_Histology)

NMB_pheno<-NMB_pheno[NMB_pheno$CPR_Histology != "NOS",]

NMB_pheno$CPR_Histology<-as.factor(NMB_pheno$CPR_Histology) 

NMB_pheno$CPR_Histology<-droplevels(NMB_pheno$CPR_Histology)

non_FS_pheno<-NMB_pheno[match(rownames(non_FS_betas), NMB_pheno$Sample_ID),]

non_FS_pheno<-non_FS_pheno[!is.na(non_FS_pheno$Sample_ID),]

if(all(!rownames(non_FS_betas) == non_FS_pheno$Sample_ID)){
  
  non_FS_betas<-non_FS_betas[match(non_FS_pheno$Sample_ID, rownames(non_FS_betas)),]
}
# Test sets 

FS_betas<-as.data.frame(t(readRDS("~/FS/FS_meth_betas.rds")))

FS_pheno<-NMB_pheno[match(rownames(FS_betas), NMB_pheno$Sample_ID),]

FS_pheno<-FS_pheno[!is.na(FS_pheno$Sample_ID),]

FS_betas<-FS_betas[match(FS_pheno$Sample_ID, rownames(FS_betas)),]

Cavalli_pheno <- readRDS("~/Cavalli/Final_Cavalli_pheno.rds")

Cavalli_pheno$histology<-gsub("MBEN","DN",Cavalli_pheno$histology)

Cavalli_pheno<- Cavalli_pheno[!is.na(Cavalli_pheno$Sample_ID), ]

Cavalli_pheno$histology <- factor(Cavalli_pheno$histology, levels = HISTO_LEVELS)

Cavalli_betas<-readRDS("~/Cavalli/Cavalli-test/Cavalli_test_betas.rds") 

Cavalli_betas<-as.data.frame(t(Cavalli_betas)) 

Cavalli_betas<-Cavalli_betas[match(Cavalli_pheno$Sample_ID, rownames(Cavalli_betas)),] 

# Training a model on non-FS


RF_non_FS<-maintrain(betas.. = non_FS_betas,
                         y.. = non_FS_pheno$CPR_Histology,
                         seed = 42,
                         method = "RF",
                         ntrees = 100, 
                         p = "200", #default
                         feature_sel = FALSE,
                         bal_feature_sel = FALSE,
                         bal_training = FALSE,
                         stratified_bal_group = NULL,
                         topfeaturenumber = NULL, 
                         imp_replace = TRUE,
                         subset.CpGs = NULL,
                         calibrationmethod = c("MR","FLR","LR") 
)

saveRDS(RF_non_FS,file.path(path.expand(FS_dir), "RF_non_FS.rds")) 

# Testing on NMB

RF_non_FS<-readRDS(file.path(path.expand(FS_dir), "RF_non_FS.rds")) 

message("Testing and evaluating on Test NMB")

if (all(rownames(FS_betas) == FS_pheno$Sample_ID)) {
  FS_results <- tryCatch(
    test_and_evaluate (
      model = RF_non_FS,
      method = "methylClass",
      calibrated = TRUE,
      calibration_method = "MR",
      train_set_name = "Non FS",
      test_data = FS_betas,
      test_data_name = "FS",
      true_test_labels = as.factor(FS_pheno$CPR_Histology),
      model_name = "Random Forest",
      save_in = file.path(path.expand(FS_dir))
    ),
    error = function(e) {
      message(sprintf("FAILED [%s / %s]: %s",
                      "Random Forest trained on non-FS", "FS", conditionMessage(e)))
      NULL
    }
  )
} else {
  stop("Mismatch in test set and test pheno Sample IDs")
}

FS_heatmap  <- prob_heatmaps(
  prob_df           = FS_results$pred_probs,
  pred              = FS_results$pred_probs$pred,
  test_pheno        = FS_pheno,
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
  heatmap_label     = "Random Forest Model for classification of histopathology,\ntrained on 'non-FS', tested on 'FS",
  heatmap_panel     = FALSE,
  column_for_panel  = "LCA",
  file_id           = "RF_non_FS_FS",
  save_in           = file.path(path.expand(FS_dir))
)

# On Cavalli

RF_non_FS<-readRDS(file.path(path.expand(FS_dir), "RF_non_FS.rds")) 

if (all(rownames(Cavalli_betas) == Cavalli_pheno$Sample_ID)) {
  test_Cavalli_results_non_FS<- tryCatch(
    test_and_evaluate (
      model = RF_non_FS,
      method = "methylClass",
      calibrated = TRUE,
      calibration_method = "MR",
      train_set_name = "Non FS",
      test_data = Cavalli_betas,
      test_data_name = "Cavalli",
      true_test_labels = as.factor(Cavalli_pheno$histology),
      model_name = "Random Forest",
      save_in = file.path(path.expand(FS_dir))
    ),
    error = function(e) {
      message(sprintf("FAILED [%s / %s]: %s",
                      " Random Forest trained on Non-FS", "Cavalli", conditionMessage(e)))
      NULL
    }
  )
} else {
  stop("Mismatch in test set and test pheno Sample IDs")
}

Test_Cavalli_heatmap_non_FS <- prob_heatmaps(
  prob_df           = test_Cavalli_results_non_FS$pred_probs,
  pred              = test_Cavalli_results_non_FS$pred_probs$pred,
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
  heatmap_label     = "Random Forest Model for classification of histopathology,\ntrained on 'non-FS', tested on 'Cavalli",
  heatmap_panel     = FALSE,
  column_for_panel  = "LCA",
  file_id           = "RF_non_FS_Cavalli",
  save_in           = file.path(path.expand(FS_dir))
) 

# --------------- Training Random Forest Models on Datasets Balanced with 'Informed' Undersampling of the Majority Class ---------------

# 'Informed' here refers to MB molecular group 'composition-aware' sampling of the train set for bootstrap of the tree training within Random Forest model, 
#  to ensure all of the molecular groups, among the histology classes are represented, also reproducing biological incidence.

# To apply the above, custom functionalities have been added to trainRF function within methylClass-fork

# ---- Custom Wrapper Functions -----

source("~/Custom_functions.R",echo = TRUE) 

# ----- Load Packages and Data ------

# Load the package every new session 

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::document("~/methylClass-fork") 
devtools::load_all("~/methylClass-fork")

US_dir<-"~/Thesis_models/Undersampled"

dir.create(US_dir, showWarnings = FALSE, recursive = TRUE)

US_dir_uncal<-"~/Thesis_models/Undersampled/Uncal/"

dir.create(US_dir_uncal, showWarnings = FALSE , recursive = TRUE)

HISTO_LEVELS<-c("CLA","LCA","DN")

NMB_train_pheno <- readRDS("~/NMB/NMB_train_pheno.rds")

NMB_train_pheno$CPR_Histology<-gsub("MBEN","DN", NMB_train_pheno$CPR_Histology)

NMB_train_pheno <- NMB_train_pheno[!is.na(NMB_train_pheno$Sample_ID), ]

NMB_train_pheno$CPR_Histology <- factor(NMB_train_pheno$CPR_Histology, levels = HISTO_LEVELS)

# ---- Test Sets ----

NMB_test_pheno <- readRDS("~/NMB/NMB_test_pheno.rds")

NMB_test_pheno$CPR_Histology<-gsub("MBEN","DN", NMB_test_pheno$CPR_Histology)

NMB_test_pheno <- NMB_test_pheno[!is.na(NMB_test_pheno$Sample_ID), ]

NMB_test_pheno$CPR_Histology<- factor(NMB_test_pheno$CPR_Histology, levels = HISTO_LEVELS)

Cavalli_pheno <- readRDS("~/Cavalli/Final_Cavalli_pheno.rds")

Cavalli_pheno$histology<-gsub("MBEN","DN",Cavalli_pheno$histology)

Cavalli_pheno<- Cavalli_pheno[!is.na(Cavalli_pheno$Sample_ID), ]

Cavalli_pheno$histology <- factor(Cavalli_pheno$histology, levels = HISTO_LEVELS)

# Also using train NMB adjusted for molecular group

train_NMB_adj<-readRDS("~/Thesis_models/mol_batch_corr/Train_NMB_mol_adj.rds")

train_NMB_adj<-train_NMB_adj[,!colnames(train_NMB_adj) %in% "histo"]

#----- Beta Matrices ----

NMB_betas<-readRDS("~/NMB/Final_NMB_betas.rds")

train_NMB<-NMB_betas[match(NMB_train_pheno$Sample_ID, rownames(NMB_betas)),]

NMB_betas<-readRDS("~/NMB/Final_NMB_betas.rds")

train_NMB<-NMB_betas[match(NMB_train_pheno$Sample_ID, rownames(NMB_betas)),]

test_NMB<-NMB_betas[match(NMB_test_pheno$Sample_ID, rownames(NMB_betas)),]

Cavalli_betas<-readRDS("~/Cavalli/Cavalli-test/Cavalli_test_betas.rds") 

Cavalli_betas<-as.data.frame(t(Cavalli_betas)) 

Cavalli_betas<-Cavalli_betas[match(Cavalli_pheno$Sample_ID, rownames(Cavalli_betas)),]

# Training the RF model 

Train_sets<-list("Train NMB " = train_NMB, "Train NMB Adjusted for Molecular Group" = train_NMB_adj )

for (train_set in names(Train_sets)){ 
  
  stopifnot(all(rownames(Train_sets[[train_set]]) == NMB_train_pheno$Sample_ID))
  
  # Training with feature selection (p = 500)
  
  RF_500<- tryCatch(maintrain(betas.. = Train_sets[[train_set]],
                              y.. = NMB_train_pheno$CPR_Histology,
                              method = "RF",
                              subset.CpGs = NULL,
                              topfeaturenumber = NULL,
                              seed = 42,
                              feature_sel = TRUE,
                              bal_feature_sel = TRUE,
                              p = 500,
                              ntrees = 500,
                              bal_training = TRUE,
                              stratified_bal_group = as.factor(NMB_train_pheno$Mol_group),
                              calibrationmethod = "MR")
                    
                    ,error = function(e) {
                      message(sprintf("FAILED [%s / %s]: %s",
                                      "Random Forest on US_Train_NMB_f500", conditionMessage(e))) })
  
  
  saveRDS(RF_500, file.path(path.expand(US_dir), paste0("RF_",train_set,"_f500",".rds")))
  
  message("Testing and evaluating on Test NMB")
  
  if (all(rownames(test_NMB) == NMB_test_pheno$Sample_ID)) {
    test_NMB_results_f500<- tryCatch(
      test_and_evaluate (
        model = RF_500,
        method = "methylClass",
        calibrated = FALSE,
        calibration_method = NULL,
        train_set_name = train_set,
        test_data = test_NMB,
        test_data_name = "Test NMB",
        true_test_labels = as.factor(NMB_test_pheno$CPR_Histology),
        model_name = "Random Forest_top_500_features",
        save_in = file.path(path.expand(US_dir_uncal))
      ),
      error = function(e) {
        message(sprintf("FAILED [%s / %s]: %s",
                        "US_f500 Random Forest trained on Train NMB", "Test NMB", conditionMessage(e)))
        NULL
      }
    )
  } else {
    stop("Mismatch in test set and test pheno Sample IDs")
  }
  
  Test_NMB_heatmap_f500 <- prob_heatmaps(
    prob_df           = test_NMB_results_f500$pred_probs,
    pred              = test_NMB_results_f500$pred_probs$pred,
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
    column_for_panel  = "LCA",
    continuum         = "LCA",
    heatmap_label     = paste0("Random Forest Model for classification of histopathology,\ntrained on", train_set,"with feature selection (500) and Informed Undersampling",",and tested on 'Test NMB'"),
    heatmap_panel     = FALSE,
    file_id           = paste0("f500_", train_set,"_test_NMB"),
    save_in           = file.path(path.expand(US_dir_uncal))
  )
  
  
  message("Testing and evaluating on Cavalli")
  
  if (all(rownames(Cavalli_betas) == Cavalli_pheno$Sample_ID)) {
    test_Cavalli_results_f500<- tryCatch(
      test_and_evaluate(
        model = RF_500,
        method = "methylClass",
        calibrated = FALSE,
        calibration_method = NULL,
        train_set_name = train_set,
        test_data = Cavalli_betas,
        test_data_name = "Cavalli",
        true_test_labels = as.factor(Cavalli_pheno$histology),
        model_name = "Random Forest_top_500_features",
        save_in = file.path(path.expand(US_dir_uncal)) 
      ),
      error = function(e) {
        message(sprintf("FAILED [%s / %s]: %s",
                        "US_f500 Random Forest trained on Train NMB", "Cavalli", conditionMessage(e)))
        NULL
      }
    )
  } else {
    stop("Mismatch in test set and test pheno Sample IDs")
  } 
  
  
  Cavalli_heatmap_f500 <- prob_heatmaps(
    prob_df           = test_Cavalli_results_f500$pred_probs,
    pred              = test_Cavalli_results_f500$pred_probs$pred,
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
    column_for_panel  = "LCA",
    continuum         = "LCA",
    heatmap_label     = paste0("Random Forest Model for classification of histopathology,\ntrained on", train_set,"with feature selection and Informed Undersampling",",and tested on 'Cavalli'"),
    heatmap_panel     = FALSE,
    file_id           = paste0("f500_", train_set,"_Cavalli"),
    save_in           = file.path(path.expand(US_dir_uncal))
  )
  
  # Without feature selection but balanced training 
  
  RF <- tryCatch(maintrain(betas.. = Train_sets[[train_set]],
                           y.. = NMB_train_pheno$CPR_Histology,
                           method = "RF",
                           subset.CpGs = NULL,
                           topfeaturenumber = NULL,
                           seed = 42,
                           feature_sel = FALSE,
                           bal_feature_sel = FALSE,
                           p = 500,
                           ntrees = 500,
                           bal_training = TRUE,
                           stratified_bal_group = as.factor(NMB_train_pheno$Mol_group),
                           calibrationmethod = "MR")
                 ,
                 error = function(e) {
                   message(sprintf("FAILED [%s / %s]: %s",
                                   "Random Forest on US_Train_NMB", conditionMessage(e))) })
  
  saveRDS(RF, file.path(path.expand(US_dir), paste0("RF_",train_set,".rds")))
  
  message("Testing and evaluating on Test NMB")
  
  if (all(rownames(test_NMB) == NMB_test_pheno$Sample_ID)) {
    test_NMB_results<- tryCatch(
      test_and_evaluate (
        model = RF,
        method = "methylClass",
        calibrated = FALSE,
        calibration_method = NULL,
        train_set_name = train_set,
        test_data = test_NMB,
        test_data_name = "Test NMB",
        true_test_labels = as.factor(NMB_test_pheno$CPR_Histology),
        model_name = "Random Forest",
        save_in = file.path(path.expand(US_dir_uncal))
      ),
      error = function(e) {
        message(sprintf("FAILED [%s / %s]: %s",
                        "US Random Forest trained on Train NMB", "Test NMB", conditionMessage(e)))
        NULL
      }
    )
  } else {
    stop("Mismatch in test set and test pheno Sample IDs")
  }
  
  Test_NMB_heatmap_US <- prob_heatmaps(
    prob_df           = test_NMB_results$pred_probs,
    pred              = test_NMB_results$pred_probs$pred,
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
    column_for_panel  = "LCA",
    continuum         = "LCA",
    heatmap_label     = paste0("Random Forest Model for classification of histopathology,\ntrained on", train_set,"with Informed Undersampling",",and tested on 'Test NMB'"),
    heatmap_panel     = FALSE,
    file_id           = paste0("US_", train_set,"_test_NMB"),
    save_in           = file.path(path.expand(US_dir_uncal))
  )
  
  
  message("Testing and evaluating on Cavalli")
  
  if (all(rownames(Cavalli_betas) == Cavalli_pheno$Sample_ID)) {
    test_Cavalli_results<- tryCatch(
      test_and_evaluate(
        model = RF,
        method = "methylClass",
        calibrated = FALSE,
        calibration_method = NULL,
        train_set_name = train_set,
        test_data = Cavalli_betas,
        test_data_name = "Cavalli",
        true_test_labels = as.factor(Cavalli_pheno$histology),
        model_name = "Random Forest",
        save_in = file.path(path.expand(US_dir_uncal)) 
      ),
      error = function(e) {
        message(sprintf("FAILED [%s / %s]: %s",
                        "US Random Forest trained on Train NMB", "Cavalli", conditionMessage(e)))
        NULL
      }
    )
  } else {
    stop("Mismatch in test set and test pheno Sample IDs")
  } 
  
  
  Cavalli_heatmap_US <- prob_heatmaps(
    prob_df           = test_Cavalli_results$pred_probs,
    pred              = test_Cavalli_results$pred_probs$pred,
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
    heatmap_label     = paste0("Random Forest Model for classification of histopathology,\ntrained on", train_set,"with Informed Undersampled",",and tested on 'Cavalli'"),
    heatmap_panel     = FALSE,
    column_for_panel  = "LCA",
    file_id           = paste0("US_", train_set,"_Cavalli"),
    save_in           = file.path(path.expand(US_dir_uncal))
  )
  
}

# ------------- Shapley Analysis ---------------------------------------------------------------------------

# ----- Load data and models if not in the env -------- 

SHAP_dir<-"~/Thesis_models/Undersampled/Shapley/"

dir.create(SHAP_dir, showWarnings = FALSE, recursive = TRUE)

US_dir<-"~/Thesis_models/Undersampled/"

HISTO_LEVELS<-c("CLA","LCA","DN")

NMB_train_pheno <- readRDS("~/NMB/NMB_train_pheno.rds")

NMB_train_pheno$CPR_Histology<-gsub("MBEN","DN", NMB_train_pheno$CPR_Histology)

NMB_train_pheno <- NMB_train_pheno[!is.na(NMB_train_pheno$Sample_ID), ]

NMB_train_pheno$CPR_Histology <- factor(NMB_train_pheno$CPR_Histology, levels = HISTO_LEVELS)

# ---- Test Sets ----

NMB_test_pheno <- readRDS("~/NMB/NMB_test_pheno.rds")

NMB_test_pheno$CPR_Histology<-gsub("MBEN","DN", NMB_test_pheno$CPR_Histology)

NMB_test_pheno <- NMB_test_pheno[!is.na(NMB_test_pheno$Sample_ID), ]

NMB_test_pheno$CPR_Histology<- factor(NMB_test_pheno$CPR_Histology, levels = HISTO_LEVELS)

Cavalli_pheno <- readRDS("~/Cavalli/Final_Cavalli_pheno.rds")

Cavalli_pheno$histology<-gsub("MBEN","DN",Cavalli_pheno$histology)

Cavalli_pheno<- Cavalli_pheno[!is.na(Cavalli_pheno$Sample_ID), ]

Cavalli_pheno$histology <- factor(Cavalli_pheno$histology, levels = HISTO_LEVELS)

NMB_betas<-readRDS("~/NMB/Final_NMB_betas.rds")

train_NMB<-NMB_betas[match(NMB_train_pheno$Sample_ID, rownames(NMB_betas)),]

NMB_betas<-readRDS("~/NMB/Final_NMB_betas.rds")

train_NMB<-NMB_betas[match(NMB_train_pheno$Sample_ID, rownames(NMB_betas)),]

test_NMB<-NMB_betas[match(NMB_test_pheno$Sample_ID, rownames(NMB_betas)),]

Cavalli_betas<-readRDS("~/Cavalli/Cavalli-test/Cavalli_test_betas.rds") 

Cavalli_betas<-as.data.frame(t(Cavalli_betas)) 

Cavalli_betas<-Cavalli_betas[match(Cavalli_pheno$Sample_ID, rownames(Cavalli_betas)),]

# Training a model on Train set that has been undersampled 'outside' of the training loop and with prior feature selection

library(dplyr)

set.seed(42)

target_n <- min(table(NMB_train_pheno$CPR_Histology))   

NMB_train_pheno_US <- NMB_train_pheno %>%
  group_by(CPR_Histology) %>%
  mutate(class_n = n()) %>%                      
  group_by(CPR_Histology, Mol_group) %>%
  mutate(cell_n = n(),                           
         keep_n = round(target_n * cell_n / class_n),   
         rnd    = sample(n())) %>%               
  filter(rnd <= keep_n) %>%                     
  ungroup() %>%
  select(-class_n, -cell_n, -keep_n, -rnd)

saveRDS(NMB_train_pheno_US, "~/Thesis_models/Undersampled/Prior_US_NMB.rds")

table(NMB_train_pheno_US$CPR_Histology, NMB_train_pheno_US$Mol_group)

train_NMB_US<-train_NMB[match(NMB_train_pheno_US$Sample_ID,rownames(train_NMB)),]

RF_US_NMB<-tryCatch(maintrain(betas.. = train_NMB_US,
                              y.. = NMB_train_pheno_US$CPR_Histology,
                              method = "RF",
                              subset.CpGs = NULL,
                              topfeaturenumber = NULL,
                              seed = 42,
                              feature_sel = TRUE,
                              bal_feature_sel = TRUE,
                              p = 500, 
                              ntrees = 500,
                              bal_training = TRUE,
                              stratified_bal_group = NULL,
                              calibrationmethod = "MR",
                              imp_replace = TRUE)
                    
                    ,error = function(e) {
                      message(sprintf("FAILED [%s]: %s",
                                      "Random Forest for Undersampled Train NMB", conditionMessage(e))) })

saveRDS(RF_US_NMB, "~/Thesis_models/Undersampled/RF_prior_US_NMB_500.rds")


message("Testing and evaluating on Test NMB")

RF_US_NMB<-readRDS("~/Thesis_models/Undersampled/RF_prior_US_NMB_500.rds")

if (all(rownames(test_NMB) == NMB_test_pheno$Sample_ID)) {
  test_NMB_results_US<- tryCatch(
    test_and_evaluate (
      model = RF_US_NMB,
      method = "methylClass",
      calibrated = TRUE,
      calibration_method = "MR",
      train_set_name = "Undersampled Train NMB",
      test_data = test_NMB,
      test_data_name = "Test NMB",
      true_test_labels = as.factor(NMB_test_pheno$CPR_Histology),
      model_name = "Random Forest",
      save_in = file.path(path.expand(US_dir))
    ),
    error = function(e) {
      message(sprintf("FAILED [%s / %s]: %s",
                      "Random Forest trained on US Train NMB", "Test NMB", conditionMessage(e)))
      NULL
    }
  )
} else {
  stop("Mismatch in test set and test pheno Sample IDs")
}

Test_NMB_heatmap_US <- prob_heatmaps(
  prob_df           = test_NMB_results_US$pred_probs,
  pred              = test_NMB_results_US$pred_probs$pred,
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
  column_for_panel  = "LCA",
  continuum         = "LCA",
  heatmap_label     = "Random Forest Model for classification of histopathology,\ntrained on Informed-Undersampled Train NMB, and tested on 'Test NMB'",
  heatmap_panel     = FALSE,
  file_id           = "prior_US_test_NMB",
  save_in           = file.path(path.expand(US_dir))
)


message("Testing and evaluating on Cavalli")

if (all(rownames(Cavalli_betas) == Cavalli_pheno$Sample_ID)) {
  test_Cavalli_results_US<- tryCatch(
    test_and_evaluate(
      model = RF_US_NMB,
      method = "methylClass",
      calibrated = TRUE,
      calibration_method = "MR",
      train_set_name ="Undersampled Train NMB",
      test_data = Cavalli_betas,
      test_data_name = "Cavalli",
      true_test_labels = as.factor(Cavalli_pheno$histology),
      model_name = "Random Forest",
      save_in = file.path(path.expand(US_dir)) 
    ),
    error = function(e) {
      message(sprintf("FAILED [%s / %s]: %s",
                      "Random Forest trained on US Train NMB", "Cavalli", conditionMessage(e)))
      NULL
    }
  )
} else {
  stop("Mismatch in test set and test pheno Sample IDs")
} 


Cavalli_heatmap_US <- prob_heatmaps(
  prob_df           = test_Cavalli_results_US$pred_probs,
  pred              = test_Cavalli_results_US$pred_probs$pred,
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
  column_for_panel  = "LCA",
  continuum         = "LCA",
  heatmap_panel     =  FALSE,
  heatmap_label     = "Random Forest Model for classification of histopathology,\ntrained on Informed-Undersampled Train NMB, and tested on 'Cavalli'",
  file_id           = "US_prior_Cavalli",
  save_in           = file.path(path.expand(US_dir))
) 


# Training an RF model to classify MB molecular groups to check for the overlap in molecular-group-specific classification features
# and histology-group-specific-features

RF_mol_grp<- tryCatch(maintrain(betas.. = train_NMB,
                                y.. = NMB_train_pheno$Mol_group,
                                method = "RF",
                                subset.CpGs = NULL,
                                topfeaturenumber = NULL,
                                seed = 42,
                                feature_sel = FALSE,
                                bal_feature_sel = FALSE,
                                p = 500, 
                                ntrees = 500,
                                bal_training = FALSE,
                                calibrationmethod = "MR")
                      
                      ,error = function(e) {
                        message(sprintf("FAILED [%s]: %s",
                                        "Random Forest for Molecular group", conditionMessage(e))) })

saveRDS(RF_mol_grp, "~/Thesis_models/Undersampled/RF_mol_grp_NMB_unbal.rds")

message("Testing and evaluating on Test NMB")

RF_mol_grp<-readRDS("~/Thesis_models/Undersampled/RF_mol_grp_NMB_unbal.rds")

if (all(rownames(test_NMB) == NMB_test_pheno$Sample_ID)) {
  test_NMB_results_mol<- tryCatch(
    test_and_evaluate (
      model = RF_mol_grp,
      method = "methylClass",
      calibrated = TRUE,
      calibration_method = "MR",
      train_set_name = "Train NMB",
      test_data = test_NMB,
      test_data_name = "Test NMB",
      true_test_labels = as.factor(NMB_test_pheno$Mol_group),
      model_name = "Random Forest",
      save_in = file.path(path.expand(US_dir))
    ),
    error = function(e) {
      message(sprintf("FAILED [%s:%s]: %s",
                      "Random Forest trained on Train NMB Mol. grp.", "Test NMB", conditionMessage(e)))
      NULL
    }
  )
} else {
  stop("Mismatch in test set and test pheno Sample IDs")
}


message("Testing and evaluating on Cavalli")

if (all(rownames(Cavalli_betas) == Cavalli_pheno$Sample_ID)) {
  test_Cavalli_results_mol<- tryCatch(
    test_and_evaluate(
      model = RF_mol_grp,
      method = "methylClass",
      calibrated = TRUE,
      calibration_method = "MR",
      train_set_name = "Train NMB",
      test_data = Cavalli_betas,
      test_data_name = "Cavalli",
      true_test_labels = as.factor(Cavalli_pheno$Mol_grp),
      model_name = "Random Forest",
      save_in = file.path(path.expand(US_dir)) 
    ),
    error = function(e) {
      message(sprintf("FAILED [%s : %s ]: %s",
                      "Random Forest trained on Train NMB", "Cavalli", conditionMessage(e)))
      NULL
    }
  )
} else {
  stop("Mismatch in test set and test pheno Sample IDs")
} 


#------- Overlap in molecular-group-specific classification features and histology-group-specific-features --------

Histo_features<-as.factor(rownames(RF_US_NMB$mod[[1]]$importance)) #500

Mol_grp_features<-as.factor(rownames(RF_mol_grp$mod[[1]]$importance)) #500

overlap<-intersect(Mol_grp_features,Histo_features)

length(overlap)

# Applying fastshap with Monte-Carlo algorithm to extract class-wise features important for classification

#------- Custom function for Shapley Analysis within multi-class classification ---------

shap_multi <- function(model,
                       train_set,
                       test_set,
                       classes,
                       calibration_method = "MR",
                       nsim = 20,
                       save_in,
                       file_id,
                       seed = 42,
                       cores = 1,                  
                       subset_features = FALSE) { 
  
  # model - A trained ML model: list methylClass object
  
  # train_set - Train set used to train the model (samples x probes): data.frame
  
  # test_set - The test set to explain (samples x probes): data.frame
  
  # classes - levels of the classification classes: character vector
  
  # calibration_method - "MR", "LR", "FLR", or NULL for no calibration
  
  # nsim - Number of Monte-Carlo repetitions for Shapley values: numeric
  
  # save_in - directory to save the class-wise output in: character
  
  # file_id - name to save the class-wise outputs with: character
  
  # seed - for reproducibility: numeric
  
  # cores - number of cores to share the work across: numeric
  
  # subset_features - if TRUE, only the model's used probes are tested for
  
  required_pgs <- c("dplyr", "ggplot2", "tidyr", "doParallel", "doRNG","tidytext")
  missing_pgs  <- required_pgs[!sapply(required_pgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pgs) > 0) install.packages(missing_pgs, dependencies = TRUE)
  invisible(sapply(required_pgs, library, character.only = TRUE, quietly = TRUE))
  
  if (!require(fastshap)) {
    install.packages("fastshap",
                     repos = c("https://bgreenwell.r-universe.dev",
                               "https://cloud.r-project.org"))
    library(fastshap)
  }
  
  calfit <- if (!is.null(calibration_method)) {
    switch(as.character(calibration_method),
           MR  = model$glmnet.calfit,
           LR  = model$platt.calfits,
           FLR = model$platt.brglm.calfits)
  } else {
    NULL
  }
  
  X_train <- as.data.frame(train_set)
  X_test  <- as.data.frame(test_set)
  
  if (subset_features) {
    feats   <- rownames(model$mod[[1]]$importance)
    feats   <- intersect(feats, colnames(X_train))
    X_train <- X_train[, feats, drop = FALSE]
    X_test  <- X_test[,  feats, drop = FALSE]
  }
  
  
  
  if (cores > 1) {
    cl <- parallel::makeCluster(cores, type = "FORK")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    doParallel::registerDoParallel(cl)
  }
  
  pred_fun <- function(class_name) {
    force(class_name)
    function(object, newdata) {
      p <- methylClass::mainpredict(mod          = object$mod,
                                    newdat       = newdata,   
                                    calibratemod = calfit)
      if (!is.null(calfit)) {
        as.numeric(p$probs[, class_name])
      } else {
        as.numeric(p$scores[, class_name])
      }
    }
  }
  
  shap_results <- list()
  
  for (class_name in classes) {
    
    pred_wrapper <- pred_fun(class_name)
    
    if (cores > 1) {
      doRNG::registerDoRNG(seed)
    } else {
      set.seed(seed)
    }
    
    shap_values <- explain(
      object       = model,
      X            = X_train,
      pred_wrapper = pred_wrapper,
      nsim         = nsim, 
      newdata      = X_test,
      adjust       = TRUE,
      parallel     = cores > 1
    )
    
    shap_results[[class_name]] <- shap_values
    
    saveRDS(shap_values,
            file.path(path.expand(save_in),
                      paste0(class_name, "_features_", file_id, ".rds")))
  }
  
  return(shap_results)
}

# ------- Histology SHAP -------

RF_US_NMB <- readRDS("~/Thesis_models/Undersampled/RF_prior_US_NMB_500.rds")

shap_histo <- shap_multi(model = RF_US_NMB,
                         train_set = train_NMB,
                         test_set = test_NMB,
                         classes = HISTO_LEVELS,
                         calibration_method = "MR",
                         nsim = 20,
                         save_in = SHAP_dir,
                         file_id = "US_train_NMB_test_NMB",
                         seed = 42,
                         cores = 5,
                         subset_features = TRUE)

saveRDS(shap_histo, file.path(path.expand(SHAP_dir), "RF_US_train_NMB_histo_features.rds"))

# ------- Molecular Group SHAP -------

RF_mol_grp <- readRDS("~/Thesis_models/Undersampled/RF_mol_grp_NMB_500.rds")

mol_classes <- c("SHH-MB", "Group3-MB", "Group4-MB", "WNT-MB")

shap_mol <- shap_multi(model = RF_mol_grp,
                       train_set = train_NMB,
                       test_set = test_NMB,
                       classes = mol_classes,
                       calibration_method = "MR",
                       nsim = 20,
                       save_in = SHAP_dir,
                       file_id = "Train_NMB_test_NMB",
                       seed = 42,
                       cores = 5,
                       subset_features = TRUE) 

saveRDS(shap_mol, file.path(path.expand(SHAP_dir), "RF_train_NMB_mol_grp_features.rds"))

# ------- Top features by SHAP values -------

top_shap <- function(shap_matrix, top_n = 100) {
  
  mean_abs_shap <- colMeans(abs(shap_matrix), na.rm = TRUE)
  
  top_features <- sort(mean_abs_shap, decreasing = TRUE)[1:min(top_n, length(mean_abs_shap))]
  
  data.frame(
    feature    = names(top_features),
    importance = as.numeric(top_features),
    rank       = seq_along(top_features)
  )
}

# Extracting for histology

histo_importance <- list()

for(class_name in names(shap_histo)) {
  histo_importance[[class_name]] <- top_shap(shap_histo[[class_name]], top_n = 500)
  histo_importance[[class_name]]$class <- class_name
  histo_importance[[class_name]]$type <- "Histology"
}

# Extracting for molecular groups

mol_importance <- list()
for(class_name in names(shap_mol)) {
  mol_importance[[class_name]] <- top_shap(shap_mol[[class_name]],top_n = 500)
  mol_importance[[class_name]]$class <- class_name
  mol_importance[[class_name]]$type <- "Molecular"
}

# Combining all

all_importance <- rbind(
  do.call(rbind, histo_importance),
  do.call(rbind, mol_importance)
)

saveRDS(all_importance, file.path(path.expand(SHAP_dir),"US_train_NMB_test_NMB_all_shap_importance.rds"))
saveRDS(histo_importance,file.path(path.expand(SHAP_dir), "US_train_NMB_test_NMB_histo_importance.rds"))
saveRDS(mol_importance,file.path(path.expand(SHAP_dir), "US_train_NMB_test_NMB_mol_grp_importance.rds"))


#-------- Overlap of top Histology SHAP features and Molecular Group SHAP features ------------------

classes_all <- c(names(shap_histo), names(shap_mol))

overlap_matrix <- matrix(0, nrow = length(classes_all), ncol = length(classes_all))

rownames(overlap_matrix) <- classes_all

colnames(overlap_matrix) <- classes_all

n_compare <- 200

for(i in 1:length(classes_all)) {
  for(j in 1:length(classes_all)) {
    
    class_i <- classes_all[i]
    class_j <- classes_all[j]
    
    if(class_i %in% names(shap_histo)) {
      features_i <- histo_importance[[class_i]]$feature[1:n_compare]
    } else {
      features_i <- mol_importance[[class_i]]$feature[1:n_compare]
    }
    
    if(class_j %in% names(shap_histo)) {
      
      features_j <- histo_importance[[class_j]]$feature[1:n_compare]
    } else {
      features_j <- mol_importance[[class_j]]$feature[1:n_compare]
    }
    
    overlap_matrix[i, j] <- length(intersect(features_i, features_j))
  }
}

saveRDS(overlap_matrix,file.path(path.expand(SHAP_dir), "US_train_NMB_test_NMB_mol_grp_histo_overlap.rds"))


# --------- Visualizing the Overlap --------------

overlap_long <- as.data.frame(as.table(overlap_matrix))

names(overlap_long) <- c("Class1", "Class2", "Overlap")

p1 <- ggplot(overlap_long, aes(x = Class1, y = Class2, fill = Overlap)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = Overlap), color = "black", size = 3) +
  scale_fill_gradient2(low = "white", mid = "lightblue", high = "darkblue",
                       midpoint = n_compare/2,
                       name = "Shared\nFeatures") +
  labs(title = "Top 200 SHAP Feature Overlap Between Histology Classes and Molecular Group Classes",
       x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        panel.grid = element_blank())

ggsave(file.path(path.expand(SHAP_dir),"shap_overlap_heatmap_200.pdf"), p1, width = 10, height = 8)


# ------------- Survival Analysis for P(LCA) --------------------

#---------- Load Data and packages if not in the env. --------

required_pgs <- c("dplyr", "ggplot2", "survival", "broom", "survminer")

missing_pgs  <- required_pgs[!sapply(required_pgs, requireNamespace, quietly = TRUE)]

if (length(missing_pgs) > 0) install.packages(missing_pgs, dependencies = TRUE)

invisible(sapply(required_pgs, library, character.only = TRUE, quietly = TRUE)) 

Cox_dir<-"~/Thesis_models/Undersampled/Survival/"

dir.create(Cox_dir, showWarnings = FALSE, recursive = TRUE)

HISTO_LEVELS<-c("CLA","LCA","DN")

NMB_test_pheno <- readRDS("~/NMB/NMB_test_pheno.rds")

NMB_test_pheno$CPR_Histology<-gsub("MBEN","DN", NMB_test_pheno$CPR_Histology)

NMB_test_pheno <- NMB_test_pheno[!is.na(NMB_test_pheno$Sample_ID), ]

NMB_test_pheno$CPR_Histology<- factor(NMB_test_pheno$CPR_Histology, levels = HISTO_LEVELS)

Cavalli_pheno <- readRDS("~/Cavalli/Final_Cavalli_pheno.rds")

Cavalli_pheno$histology<-gsub("MBEN","DN",Cavalli_pheno$histology)

Cavalli_pheno<- Cavalli_pheno[!is.na(Cavalli_pheno$Sample_ID), ]

Cavalli_pheno$histology <- factor(Cavalli_pheno$histology, levels = HISTO_LEVELS)

# ----------- Cox Regression for Survival Analysis of 'LCAness' as assigned by the model ---------------

# In the previous steps, None of the probe features important for classification of Histology were seen to 
# be important for classification of MB molecular group, suggesting weak confounding effect of Mol.grp.

#In LCA however, the P(LCA) scores in True and False Positive LCA were seen to be relatively higher in the poor survival group

# and hence, with Cox regression models, we are testing if P(LCA) could be a predictor variable for Cox survival models and can 

# and could alter hazard ratios.

# THE MODEL BEING TESTED FOR IS RF TRAINED ON UNDERSAMPLED TRAIN NMB (Undersampled outside the RF training loop)

test_NMB_probs<-read.csv("~/Thesis_models/Undersampled/Random_Forest_Undersampled_Train_NMB_Test_NMB_probabilities.csv")

#Scaling it to measure HR per 10% increase in P(LCA)

test_NMB_probs$LCAness<-(test_NMB_probs$LCA)*10

stopifnot("Rownames of Test NMB probabilities doesn't match Pheno Samples' ID" = test_NMB_probs$Sample_ID == NMB_test_pheno$Sample_ID)

test_NMB_probs$OS<-as.numeric(NMB_test_pheno$OS)

test_NMB_probs$Follow_up<-as.numeric(NMB_test_pheno$`Follow up years`)

test_NMB_probs$Mol_grp<-as.factor(NMB_test_pheno$Mol_group)

Cavalli_probs<-read.csv("~/Thesis_models/Undersampled/Random_Forest_Undersampled_Train_NMB_Cavalli_probabilities.csv")

Cavalli_probs$LCAness<-(Cavalli_probs$LCA)*10

stopifnot("Rownames of Cavalli probabilities doesn't match Pheno Samples' ID" = Cavalli_probs$Sample_ID == Cavalli_pheno$Sample_ID)

Cavalli_probs$OS<-as.numeric(Cavalli_pheno$Dead)

Cavalli_probs$Follow_up<-as.numeric(Cavalli_pheno$`OS (years)`)

Cavalli_probs$Mol_grp<-as.factor(Cavalli_pheno$Mol_grp)

#---------------- Cox models for Test_NMB --------------------------------

cox_NMB_unadj <- coxph(Surv(Follow_up, OS) ~ LCAness , 
                       data = test_NMB_probs, id = rownames(test_NMB_probs))

summary(cox_NMB_unadj)

saveRDS(cox_NMB_unadj, file.path(path.expand(Cox_dir),"Test_NMB_Cox_model_unadj.rds"))

# Adjusting for Molecular group to check if any significance present is irrespective of the confounding effect of Molecular group

cox_NMB_adj <- coxph(Surv(Follow_up, OS) ~ LCAness + Mol_grp, 
                     data = test_NMB_probs, id = rownames(test_NMB_probs))

summary(cox_NMB_adj)

saveRDS(cox_NMB_adj, file.path(path.expand(Cox_dir),"Test_NMB_Cox_model_adj.rds"))

#---------------- Cox models for Cavalli --------------------------------

cox_cavalli_unadj <- coxph(Surv(Follow_up, OS) ~ LCAness, 
                           data = Cavalli_probs)

summary(cox_cavalli_unadj)

saveRDS(cox_cavalli_unadj, file.path(path.expand(Cox_dir),"Test_cavalli_Cox_model_unadj.rds"))

# Adjusting for Molecular group to check if any significance present is irrespective of the confounding effect of Molecular group

cox_cavalli_adj <- coxph(Surv(Follow_up, OS) ~ LCAness + Mol_grp, 
                         data = Cavalli_probs)

saveRDS(cox_cavalli_adj, file.path(path.expand(Cox_dir),"Test_cavalli_Cox_model_adj.rds"))


summary(cox_cavalli_adj)


lca_results <- function(cox_model, cohort, adjustment) {
  
  result <- tidy(cox_model, exponentiate = TRUE, conf.int = TRUE)
  
  lca_result <- result %>% 
    filter(grepl("LCAness", term)) %>%
    mutate(
      cohort = cohort,
      adjustment = adjustment,
      n = cox_model$n,
      events = cox_model$nevent
    )
  
  return(lca_result)
}


forest_data <- bind_rows(
  lca_results(cox_NMB_unadj,  "Test-NMB",     "Unadjusted"),
  lca_results(cox_NMB_adj,    "Test-NMB",     "Adjusted"),
  lca_results(cox_cavalli_unadj,   "Test-Cavalli", "Unadjusted"),
  lca_results(cox_cavalli_adj,     "Test-Cavalli", "Adjusted")
)

forest_data <- forest_data %>%
  mutate(
    cohort     = factor(cohort,     levels = c("Test-NMB", "Test-Cavalli")),
    adjustment = factor(adjustment, levels = c("Unadjusted", "Adjusted")),
    label = sprintf("%s - %s\n(N=%d, Events=%d)",
                    as.character(cohort), as.character(adjustment),
                    n, events)
  ) %>%
  arrange(cohort, adjustment) %>%                 
  mutate(label = factor(label, levels = rev(label)))


forest_data %>%
  select(cohort, adjustment, estimate, conf.low, conf.high, p.value)


forest_plot <- ggplot(forest_data, aes(x = estimate, y = label)) +
  
  # Null effect line
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "#D32F2F", linewidth = 1) +
  
  # 95% CI
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high, colour = adjustment),
                orientation = "y", width = 0.3, linewidth = 1.5) +
  
  # Hazard ratio
  geom_point(aes(colour = adjustment, shape = adjustment), size = 5) +
  
  # HR and p-value
  geom_text(aes(x = 3.5,
                label = sprintf("HR: %.2f (%.2f-%.2f)\n%s",
                                estimate, conf.low, conf.high,
                                ifelse(p.value < 0.001, "p < 0.001",
                                       sprintf("p = %.3f", p.value)))),
            hjust = 1, size = 3, fontface = "bold", lineheight = 0.9) +
  
  scale_colour_manual(values = c(Unadjusted = "#0277BD", Adjusted = "#F57C00"),
                      name = "Model") +
  scale_shape_manual(values = c(Unadjusted = 15, Adjusted = 17),
                     name = "Model") +
  scale_x_log10(breaks = c(0.5, 1, 1.5, 2, 3,4), limits = c(0.5, 4))+
  
  labs(
    title = "P(LCA) and Overall Survival with and without Molecular Group Adjustment",  
    x = "Hazard ratio (95% CI) per 10% increase",   
    y = NULL
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "top",
    legend.title       = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 11, hjust = 1, face = "bold"),
    axis.text.x        = element_text(size = 10),
    plot.title         = element_text(hjust = 0.5, face = "bold", size = 14)
  )


ggsave(file.path(path.expand(Cox_dir),"forest_plot_adjusted_all_cohorts.png"), forest_plot,
       width = 14, height = 8, dpi = 300, bg = "white")


# ------------------ Applying Class-wise Probability Score thresholds ----------------------------

# ------------- Custom function to calculate thresholds and apply them for predictions -----------


youden_thresholds <- function(prob_matrix, true_labels, classes) {
  
  # prob_matrix: out-of-fold/out-of-bag train probabilities
  
  thr <- sapply(classes, function(cls) {
    binary_labels <- ifelse(true_labels == cls, 1, 0)      
    roc_class <- roc(binary_labels, prob_matrix[, cls],   
                     levels = c(0, 1), direction = "<", quiet = TRUE)
    coords(roc_class, "best", best.method = "youden",
           ret = "threshold")$threshold[1]                
  })
  names(thr) <- classes
  thr                                                      
}


predict_with_thresholds <- function(prob_matrix, thresholds, classes,
                                    allow_uncertain = FALSE) {
  
  # To address mis-classifications as a result of 'lapse in judgement' in the terms of probability scores
  
  # of the classes involved varying by <0.2, when sequential = TRUE, the class with most False Positives is 
  
  # predicted first with its specific threshold (LCA here), so once the True Positives are predicted, the rest
  
  # are to be followed up by other classes. 
  
  prob_matrix<-prob_matrix %>% mutate(seq_pred = case_when(
    
    # Priority to getting LCA right - can change accordingly
    
    LCA >= thresholds[["LCA"]] ~ "LCA",
    
    CLA >= thresholds[["CLA"]] ~ "CLA",
    
    DN >= thresholds[["DN"]] ~ "DN",
  )
  )
  
  prob_matrix$seq_pred<-as.factor(prob_matrix$seq_pred) 
  
  # Applying class-wise thresholds simultaneously
  
  scaled <- sweep(as.matrix(prob_matrix[, classes]), 2, thresholds[classes], "/")
  
  youden_pred <- classes[max.col(scaled, ties.method = "first")]
  
  if (allow_uncertain) {
    youden_pred[apply(scaled, 1, max) < 1] <- "Uncertain"        
  }
  
  prob_matrix$youden_pred<-as.factor(youden_pred)
  
  return(prob_matrix)
}

# --------- Load packages and data if not in the env.------------------

source("~/Custom_functions.R",echo = TRUE) 

required_pgs <- c("dplyr", "ggplot2", "caret", "pROC", "patchwork", "stringr")

missing_pgs  <- required_pgs[!sapply(required_pgs, requireNamespace, quietly = TRUE)]

if (length(missing_pgs) > 0) install.packages(missing_pgs, dependencies = TRUE)

invisible(sapply(required_pgs, library, character.only = TRUE, quietly = TRUE)) 

youden_dir<-"~/Thesis_models/Undersampled/Youden/"

dir.create(youden_dir, showWarnings = FALSE, recursive = TRUE)

# OOB Train NMB Probabilities

HISTO_LEVELS<-c("CLA","LCA","DN")

RF_US_NMB<-readRDS("~/Thesis_models/Undersampled/RF_prior_US_NMB_500.rds")

US_NMB_pheno<-readRDS("~/Thesis_models/Undersampled/Prior_US_NMB.rds")

train_probs<-as.data.frame(RF_US_NMB$probs.mr)

train_probs$pred<-as.factor(colnames(train_probs)[max.col(train_probs, ties.method = "first")])

stopifnot(rownames(train_probs) == US_NMB_pheno$Sample_ID)

# Youden J's class-wise thresholds from train probabilities (calibrated)

thresholds<-youden_thresholds(prob_matrix = train_probs,
                              true_labels = US_NMB_pheno$CPR_Histology,
                              classes = HISTO_LEVELS)

# Predicting with class-wise thresholds 

# Test NMB

NMB_test_pheno <- readRDS("~/NMB/NMB_test_pheno.rds")

NMB_test_pheno$CPR_Histology<-gsub("MBEN","DN", NMB_test_pheno$CPR_Histology)

NMB_test_pheno <- NMB_test_pheno[!is.na(NMB_test_pheno$Sample_ID), ]

NMB_test_pheno$CPR_Histology<- factor(NMB_test_pheno$CPR_Histology, levels = HISTO_LEVELS)

#Probabilities

test_NMB_probs<-read.csv("~/Thesis_models/Undersampled/Random_Forest_Undersampled_Train_NMB_Test_NMB_probabilities.csv")

rownames(test_NMB_probs)<-as.factor(test_NMB_probs$Sample_ID)

test_NMB_probs<-predict_with_thresholds(prob_matrix = test_NMB_probs,thresholds = thresholds, classes = HISTO_LEVELS,allow_uncertain = FALSE)


# Confusion Matrix for Class-wise thresholds predictions

stopifnot(test_NMB_probs$Sample_ID == NMB_test_pheno$Sample_ID)

cm <- confusionMatrix(test_NMB_probs$youden_pred, NMB_test_pheno$CPR_Histology, mode = "everything")

cm_df <- as.data.frame(cm$table) 

COLOURS <- c("DN" = "darkgreen", "non_DN" = "purple", "non-DN" = "purple",
             "CLA" = "purple", "LCA" = "pink", "MBEN" = "green")

cm_plot <- ggplot(cm_df, aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = Reference, alpha = Freq), color = "white") +
  geom_text(aes(label = Freq), size = 4,
            color = ifelse(cm_df$Freq > max(cm_df$Freq) * 0.5, "white", "black"), fontface = "bold") +
  scale_fill_manual(values = COLOURS, guide = "none") +
  scale_alpha_continuous(range = c(0.15, 1), guide = "none") +
  labs(title = str_wrap(paste0("Confusion Matrix for RF trained on Undersampled Train NMB and tested on Test NMB,with class-wise Youden's thresholds"),width = 60),
       x = "True Class", y = "Predicted Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title = element_text(hjust = 0.5, face = "bold"),
        plot.margin = ggplot2::margin(t = 15, r = 15, b = 15, l = 15)) +
  coord_cartesian(clip = "off")

ggsave(file.path(path.expand(youden_dir), "Test_NMB_youden_pred.png"),cm_plot, width = 9, height = 6, dpi = 300, bg = "white")


Test_NMB_heatmap_yd <- prob_heatmaps(
  prob_df           = test_NMB_probs,
  pred              = test_NMB_probs$youden_pred,
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
  column_for_panel  = "LCA",
  continuum         = "LCA",
  heatmap_label     = "Random Forest Model for classification of histopathology,\ntrained on with Informed Undersampling,and tested on 'Test NMB', with Youden's Class-wise thresholds",
  heatmap_panel     = FALSE,
  file_id           = "youden_pred_test_NMB",
  save_in           = file.path(path.expand(youden_dir))
)

# Confusion Matrix for Sequential application of thresholds 

stopifnot(test_NMB_probs$Sample_ID == NMB_test_pheno$Sample_ID)

cm <- confusionMatrix(test_NMB_probs$seq_pred, NMB_test_pheno$CPR_Histology, mode = "everything")

cm_df <- as.data.frame(cm$table) 

COLOURS <- c("DN" = "darkgreen", "non_DN" = "purple", "non-DN" = "purple",
             "CLA" = "purple", "LCA" = "pink", "MBEN" = "green")

cm_plot <- ggplot(cm_df, aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = Reference, alpha = Freq), color = "white") +
  geom_text(aes(label = Freq), size = 4,
            color = ifelse(cm_df$Freq > max(cm_df$Freq) * 0.5, "white", "black"), fontface = "bold") +
  scale_fill_manual(values = COLOURS, guide = "none") +
  scale_alpha_continuous(range = c(0.15, 1), guide = "none") +
  labs(title = str_wrap(paste0("Confusion Matrix for RF trained on Undersampled Train NMB and tested on Test NMB,with sequential class-wise Youden's thresholds"),width = 60),
       x = "True Class", y = "Predicted Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title = element_text(hjust = 0.5, face = "bold"),
        plot.margin = ggplot2::margin(t = 15, r = 15, b = 15, l = 15)) +
  coord_cartesian(clip = "off")

ggsave(file.path(path.expand(youden_dir), "Test_NMB_seq_pred.png"),cm_plot, width = 9, height = 6, dpi = 300, bg = "white")


Test_NMB_heatmap_seq <- prob_heatmaps(
  prob_df           = test_NMB_probs,
  pred              = test_NMB_probs$seq_pred,
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
  column_for_panel  = "LCA",
  continuum         = "LCA",
  heatmap_label     = "Random Forest Model for classification of histopathology,\ntrained on with Informed Undersampling,and tested on 'Test NMB', with sequential Youden's Class-wise thresholds",
  heatmap_panel     = FALSE,
  file_id           = "seq_pred_test_NMB",
  save_in           = file.path(path.expand(youden_dir))
)



# Cavalli 

Cavalli_pheno <- readRDS("~/Cavalli/Final_Cavalli_pheno.rds")

Cavalli_pheno$histology<-gsub("MBEN","DN",Cavalli_pheno$histology)

Cavalli_pheno<- Cavalli_pheno[!is.na(Cavalli_pheno$Sample_ID), ]

Cavalli_pheno$histology <- factor(Cavalli_pheno$histology, levels = HISTO_LEVELS)

#Probabilities

Cavalli_probs<-read.csv("~/Thesis_models/Undersampled/Random_Forest_Undersampled_Train_NMB_Cavalli_probabilities.csv")

rownames(Cavalli_probs)<-as.factor(Cavalli_probs$Sample_ID)

Cavalli_probs<-predict_with_thresholds(prob_matrix = Cavalli_probs,thresholds = thresholds,classes = HISTO_LEVELS,allow_uncertain = FALSE)

# Confusion Matrix for Class-wise thresholds predictions

stopifnot(Cavalli_probs$Sample_ID == Cavalli_pheno$Sample_ID)

cm <- confusionMatrix(Cavalli_probs$youden_pred, Cavalli_pheno$histology, mode = "everything")

cm_df <- as.data.frame(cm$table) 

COLOURS <- c("DN" = "darkgreen", "non_DN" = "purple", "non-DN" = "purple",
             "CLA" = "purple", "LCA" = "pink", "MBEN" = "green")

cm_plot <- ggplot(cm_df, aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = Reference, alpha = Freq), color = "white") +
  geom_text(aes(label = Freq), size = 4,
            color = ifelse(cm_df$Freq > max(cm_df$Freq) * 0.5, "white", "black"), fontface = "bold") +
  scale_fill_manual(values = COLOURS, guide = "none") +
  scale_alpha_continuous(range = c(0.15, 1), guide = "none") +
  labs(title = str_wrap(paste0("Confusion Matrix for RF trained on Undersampled Train NMB and tested on Cavalli,with class-wise Youden's thresholds"),width = 60),
       x = "True Class", y = "Predicted Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title = element_text(hjust = 0.5, face = "bold"),
        plot.margin = ggplot2::margin(t = 15, r = 15, b = 15, l = 15)) +
  coord_cartesian(clip = "off")

ggsave(file.path(path.expand(youden_dir), "Cavalli_youden_pred.png"),cm_plot, width = 9, height = 6, dpi = 300, bg = "white")


Cavalli_heatmap_yd <- prob_heatmaps(
  prob_df           = Cavalli_probs,
  pred              = Cavalli_probs$youden_pred,
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
  column_for_panel  = "LCA",
  continuum         = "LCA",
  heatmap_label     = "Random Forest Model for classification of histopathology,\ntrained on with Informed Undersampling,and tested on 'Cavalli', with Youden's Class-wise thresholds",
  heatmap_panel     = FALSE,
  file_id           = "youden_pred_Cavalli",
  save_in           = file.path(path.expand(youden_dir))
)

# Confusion Matrix for Sequential application of thresholds 

stopifnot(Cavalli_probs$Sample_ID == Cavalli_pheno$Sample_ID)

cm <- confusionMatrix(Cavalli_probs$seq_pred, Cavalli_pheno$histology, mode = "everything")

cm_df <- as.data.frame(cm$table) 

COLOURS <- c("DN" = "darkgreen", "non_DN" = "purple", "non-DN" = "purple",
             "CLA" = "purple", "LCA" = "pink", "MBEN" = "green")

cm_plot <- ggplot(cm_df, aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = Reference, alpha = Freq), color = "white") +
  geom_text(aes(label = Freq), size = 4,
            color = ifelse(cm_df$Freq > max(cm_df$Freq) * 0.5, "white", "black"), fontface = "bold") +
  scale_fill_manual(values = COLOURS, guide = "none") +
  scale_alpha_continuous(range = c(0.15, 1), guide = "none") +
  labs(title = str_wrap(paste0("Confusion Matrix for RF trained on Undersampled Train NMB and tested on Cavalli,with sequential class-wise Youden's thresholds"),width = 60),
       x = "True Class", y = "Predicted Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title = element_text(hjust = 0.5, face = "bold"),
        plot.margin = ggplot2::margin(t = 15, r = 15, b = 15, l = 15)) +
  coord_cartesian(clip = "off")

ggsave(file.path(path.expand(youden_dir), "Cavalli_seq_pred.png"),cm_plot, width = 9, height = 6, dpi = 300, bg = "white")


Cavalli_heatmap_seq <- prob_heatmaps(
  prob_df           = Cavalli_probs,
  pred              = Cavalli_probs$seq_pred,
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
  column_for_panel  = "LCA",
  continuum         = "LCA",
  heatmap_label     = "Random Forest Model for classification of histopathology,\ntrained on with Informed Undersampling,and tested on 'Cavalli', with sequential Youden's Class-wise thresholds",
  heatmap_panel     = FALSE,
  file_id           = "seq_pred_Cavalli",
  save_in           = file.path(path.expand(youden_dir))
)

# Checking for the change in classification with class-wise threshold-based predictions (especially "CLA->LCA")

stopifnot(test_NMB_probs$Sample_ID == NMB_test_pheno$Sample_ID)

test_NMB_probs$true_label<-as.factor(NMB_test_pheno$CPR_Histology)

classify_outcome <- function(true_label, pred) {
  truth <- as.character(true_label)
  pred  <- as.character(pred)
  ifelse(truth == pred, paste0("Correct_", truth), paste0(truth, "->", pred))
}

test_NMB_probs$Class<-classify_outcome(test_NMB_probs$true_label,test_NMB_probs$pred)

test_NMB_probs$Class_youden<-classify_outcome(test_NMB_probs$true_label,test_NMB_probs$youden_pred)

test_NMB_probs$Class_seq<-classify_outcome(test_NMB_probs$true_label,test_NMB_probs$seq_pred)

tradeoff_data <- test_NMB_probs %>%
  filter(true_label %in% c("CLA", "LCA")) %>%
  mutate(
    Status_youden = case_when(
      
      true_label == "CLA" & Class == "CLA->LCA" & Class_youden == "Correct_CLA" ~ "CLA Recovered",
      
      true_label == "LCA" & Class == "Correct_LCA" & Class_youden== "LCA->CLA" ~ "LCA Now Wrong",
      
      true_label == "CLA" & Class == "CLA->LCA" & Class_youden == "CLA->LCA" ~ "CLA Still Wrong",
      
      TRUE ~ "No Change"
    )
  ) %>% mutate (
    Status_seq = case_when(
      
      true_label == "CLA" & Class == "CLA->LCA" & Class_seq == "Correct_CLA" ~ "CLA Recovered",
      
      true_label == "LCA" & Class == "Correct_LCA" & Class_seq == "LCA->CLA" ~ "LCA Now Wrong",
      
      true_label == "CLA" & Class == "CLA->LCA" & Class_seq  == "CLA->LCA" ~ "CLA Still Wrong",
      
      TRUE ~ "No Change"
    )
  )

#-------------------- Plot to visualize the effect of applying class-wise thresholds-----------------------

pred_change_plot <- function(data,
                             thresholds,
                             pred_col    = "youden_pred",
                             cohort_name,
                             sequential  = FALSE,      
                             lims        = c(0.05, 0.65),
                             save_in     = youden_dir,
                             file_id     = "pred_change_NMB.png") {
  
  dn_thr  <- thresholds[["DN"]]
  data$pred_plot <- data[[pred_col]]
  
  type_pred_title <- if (sequential) {
    sprintf("P(LCA) vs P(CLA) for predictions with Youden's thresholds, applied sequentially, in %s", cohort_name)
  } else {
    sprintf("P(LCA) vs P(CLA) for predictions with Youden's thresholds, in %s", cohort_name)
  }
  
  p <- ggplot(data, aes(x = CLA, y = LCA)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "blue", alpha = 0.6) +
    geom_vline(xintercept = thresholds[["CLA"]], linetype = "dashed",
               colour = "purple", alpha = 0.5) +
    geom_hline(yintercept = thresholds[["LCA"]], linetype = "dashed",
               colour = "hotpink", alpha = 0.6) +
    geom_point(aes(fill = pred_plot, shape = true_label),
               size = 3.5, stroke = 1, alpha = 0.9) +
    annotate("label", x = thresholds[["CLA"]], y = lims[2] - 0.01,
             label = sprintf("CLA threshold = %.2f", thresholds[["CLA"]]),
             colour = "purple", fill = "white", alpha = 0.8,
             linewidth = 0, size = 3, hjust = 0.5, vjust = 1) +
    annotate("label", x = lims[2] - 0.01, y = thresholds[["LCA"]],
             label = sprintf("LCA threshold = %.2f", thresholds[["LCA"]]),
             colour = "hotpink3", fill = "white", alpha = 0.8,
             linewidth = 0, size = 3, hjust = 1, vjust = -0.3) +
    annotate("label", x = lims[1] + 0.02, y = lims[1] + 0.05,      
             label = "P(CLA) = P(LCA)",
             colour = "blue4", fill = "white", alpha = 0.8, 
             linewidth = 0, size = 3, hjust = 0) +
    scale_fill_manual(
      values = c(CLA = "purple", LCA = "pink", DN = "darkgreen"),
      name   = "Predicted label",
      guide  = guide_legend(override.aes = list(shape = 21, colour = "grey30"), order = 1)) +
    scale_shape_manual(
      values = c(CLA = 21, LCA = 24, DN = 23),
      name   = "True label",
      guide  = guide_legend(override.aes = list(fill = "white", colour = "grey30"), order = 2)) +
    labs(title = type_pred_title, x = "P(CLA)", y = "P(LCA)") +
    theme_minimal() +
    theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5)) +
    coord_cartesian(xlim = lims, ylim = lims, expand = FALSE)
  
  ggsave(file.path(path.expand(save_in), file_id), p,
         width = 12, height = 8, dpi = 300, bg = "white")
  return(p)
}

p_seq <- pred_change_plot(
  test_NMB_probs_CLA_LCA,
  thresholds  = thresholds,
  pred_col    = "youden_pred",
  cohort_name = "Test NMB",
  sequential  = FALSE,
  save_in     = youden_dir,
  file_id     = "youden_pred_NMB.png"
)

p_seq <- pred_change_plot(
  test_NMB_probs_CLA_LCA,
  thresholds  = thresholds,
  pred_col    = "seq_pred",
  cohort_name = "Test NMB",
  sequential  = TRUE,
  save_in     = youden_dir,
  file_id     = "seq_pred_NMB.png"
)


# -------------- Defining Intermediate/Uncertainty Zone ---------------------------------------------

# if |P(CLA) -P(LCA)| < 0.15, it could be suggested that the model is genuinely confused about either of the classes

# which could explain the hypothesis of CLA and LCA exisiting on a continuum

uncertainty <- function(data, dn_thr, diff_cut = 0.15) {
  sum_min <- 1 - dn_thr
  data %>%
    dplyr::mutate(
      prob_sum  = CLA + LCA,
      prob_diff = abs(CLA - LCA),
      uncertainty_status = ifelse(
        prob_diff < diff_cut & prob_sum > sum_min,
        "Uncertain", "Confident")
    )
}

uncertainty_plot <- function(data,
                             thresholds,
                             pred_col    = "youden_pred",
                             cohort_name,
                             diff_cut    = 0.2,
                             lims        = c(0.05, 0.75),
                             save_in     = youden_dir,
                             file_id     = "uncertainty_NMB.png") {
  
  dn_thr  <- thresholds[["DN"]]
  sum_min <- 1 - dn_thr
  data    <- uncertainty(data, dn_thr, diff_cut)
  data$pred_plot <- data[[pred_col]]
  
  n_unc   <- sum(data$uncertainty_status == "Uncertain", na.rm = TRUE)
  pct_unc <- 100 * n_unc / nrow(data)
  
  
  zone <- data.frame(x = seq(lims[1] - 0.1, lims[2] + 0.1, length.out = 1000)) %>%
    dplyr::mutate(lo = pmax(x - diff_cut, sum_min - x),
                  hi = x + diff_cut) %>%
    dplyr::filter(lo < hi)
  
  p <- ggplot(data, aes(x = CLA, y = LCA)) +
    geom_ribbon(data = zone, aes(x = x, ymin = lo, ymax = hi),
                fill = "red", alpha = 0.12, inherit.aes = FALSE) +
    geom_abline(slope = 1, intercept = c(-diff_cut, diff_cut),
                linetype = "dotted", colour = "blue", alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "blue", alpha = 0.6) +
    geom_abline(slope = -1, intercept = sum_min,
                linetype = "dotted", colour = "orange", alpha = 0.9) +
    geom_vline(xintercept = thresholds[["CLA"]], linetype = "dashed",
               colour = "purple", alpha = 0.5) +
    geom_hline(yintercept = thresholds[["LCA"]], linetype = "dashed",
               colour = "hotpink", alpha = 0.6) +
    geom_point(aes(fill = pred_plot, shape = true_label, colour = uncertainty_status),
               size = 3.5, stroke = 1, alpha = 0.9) +
    annotate("label", x = thresholds[["CLA"]], y = lims[2] - 0.01,
             label = sprintf("CLA threshold = %.2f", thresholds[["CLA"]]),
             colour = "purple", fill = "white", alpha = 0.8,
             linewidth = 0, size = 3, hjust = 0.5, vjust = 1) +
    
    annotate("label", x = lims[2] - 0.01, y = thresholds[["LCA"]],
             label = sprintf("LCA threshold = %.2f", thresholds[["LCA"]]),
             colour = "hotpink3", fill = "white", alpha = 0.8,
             linewidth = 0, size = 3, hjust = 1, vjust = -0.3) +
    
    annotate("label", x = sum_min - (lims[2] - 0.03), y = lims[2] - 0.03,
             label = sprintf("%.2f (1 - DN thr)", sum_min),
             colour = "darkorange3", fill = "white", alpha = 0.8,
             linewidth = 0, size = 3, hjust = 0) +
    
    annotate("label", x = lims[1] + 0.02, y = lims[1] + 0.02 + 0.03,
             label = "P(CLA) = P(LCA)",
             colour = "blue4", fill = "white", alpha = 0.8,
             linewidth = 0, size = 3, hjust = 0) +
    
    scale_fill_manual(
      values = c(CLA = "purple", LCA = "pink", DN = "darkgreen"),
      name   = "Predicted label",
      guide  = guide_legend(override.aes = list(shape = 21, colour = "grey30"), order = 1)) +
    scale_shape_manual(
      values = c(CLA = 21, LCA = 24, DN = 23),
      name   = "True label",
      guide  = guide_legend(override.aes = list(fill = "white", colour = "grey30"), order = 2)) +
    scale_colour_manual(
      values = c(Confident = "grey30", Uncertain = "red"),
      name   = "Outline",
      guide  = guide_legend(override.aes = list(shape = 21, fill = "white"), order = 3)) +
    labs(
      title = sprintf("P(LCA) vs P(CLA) with uncertainty zone in %s", cohort_name),
      subtitle = sprintf(
        "Uncertain: |P(LCA) - P(CLA)| < %.2f and P(CLA) + P(LCA) > %.2f [1-P(DN) Threshold] | %d samples (%.0f%%)",
        diff_cut, sum_min, n_unc, pct_unc),
      x = "P(CLA)", y = "P(LCA)") +
    theme_minimal() +
    theme(plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5)) +
    coord_cartesian(xlim = lims, ylim = lims, expand = FALSE)
  
  ggsave(file.path(path.expand(save_in), file_id), p,
         width = 12, height = 8, dpi = 300, bg = "white")
  return(p)
}


uncertainty_summary <- function(data, dn_thr, diff_cut = 0.2) {
  data <- uncertainty(data, dn_thr, diff_cut)
  list(
    by_group = data %>%
      dplyr::count(uncertainty_status, pred, true_label) %>%
      tidyr::pivot_wider(names_from = uncertainty_status,
                         values_from = n, values_fill = 0),
    overall = data %>%
      dplyr::summarise(total          = dplyr::n(),
                       in_uncertainty = sum(uncertainty_status == "Uncertain"),
                       pct_uncertain  = 100 * in_uncertainty / total))
}

test_NMB_probs_CLA_LCA <- test_NMB_probs %>%
  dplyr::filter(true_label %in% c("CLA", "LCA"))

uncertainty_summary(test_NMB_probs_CLA_LCA,
                    dn_thr   = thresholds[["DN"]],
                    diff_cut = 0.2)


p_youden_uncertain <- uncertainty_plot(
  test_NMB_probs_CLA_LCA,
  thresholds  = thresholds,
  pred_col    = "youden_pred",
  cohort_name = "Test NMB",
  diff_cut    = 0.15,
  save_in     = youden_dir,
  file_id     = "youden_pred_uncertainty_NMB.png"
)

p_seq_uncertain <- uncertainty_plot(
  test_NMB_probs_CLA_LCA,
  thresholds  = thresholds,
  pred_col    = "seq_pred",
  cohort_name = "Test NMB",
  diff_cut    = 0.15,
  save_in     = youden_dir,
  file_id     = "seq_pred_uncertainty_NMB.png"
)







