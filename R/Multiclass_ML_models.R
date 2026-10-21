#===================================================================================================================#
#Train Machine Learning Models for Classification of Histopathology in Medulloblastoma

# This script includes loading of samplesheets and pre processed methylation data (in-house),
# Feature selection of Top most variable probes,
# Initial training of models on imbalanced data, 
# Training and optimising of caret and methylClass models- Random Forest(RF),Support Vector Machine (SVM) and eXtreme Gradient Boosting (XGB,
# Testing of models on independent test sets and subsequent evaluation metrics
# Balancing the train set with Synthetic Minority Oversampling Technique (SMOTE)
# RF trained on balanced data and evaluated on test sets 
# The 'NMB' train cohort stratified to train and test models on Central Pathology Review (CPR) cohorts

#==================================================================================================================#

# ------ Load packages ------ #

require(matrixStats)

require(dplyr) 

require(doParallel)

require(caret)

require(MLmetrics)

require(randomForest)

require(kernlab)

require(xgboost)

require(methylClass) 

require(pROC)

# ------ Load data ------ #

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

if (any(is.na(Cavalli_pheno$Study_ID...12))){
  Cavalli_pheno<-Cavalli_pheno[!is.na(Cavalli_pheno$Study_ID...12),] 
}

Cavalli_betas<-readRDS("~/Cavalli/Cavalli-test/Cavalli_test_betas.rds") 

Cavalli_betas<-as.data.frame(t(Cavalli_betas)) 

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

Cavalli_betas<-Cavalli_betas[match(Cavalli_pheno$Study_ID...12, rownames(Cavalli_betas)),]

if (all(Cavalli_pheno$Study_ID...12 == rownames(Cavalli_betas))){
  Cavalli_pheno<-Cavalli_pheno[match(rownames(Cavalli_betas), Cavalli_pheno$Study_ID...12),]
} 

#----------- Training caret (Kuhn, 2007) models with top 10k probes with 5 fold cross validation
# with parallel processing ----------------------------------------------------------------------#

# Folds for Cross Validation - 5 folds and 3 repeats

set.seed(42)

train_NMB_sd_10k<-train_NMB_sd_10k[match(rownames(train_NMB),rownames(train_NMB_sd_10k)),]

train_NMB_mad_10k<-train_NMB_sd_10k[match(rownames(train_NMB),rownames(train_NMB_mad_10k)),]

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
# is computationally inexpensive in comparison to caret -----------------------------------------------------------------------------------------#

# Random Forest


RF_train_NMB<-maintrain(betas.. = train_NMB[,!colnames(train_NMB) %in% "histo"],
                        y.. = train_NMB$histo,
                        seed = 42,
                        method = "RF",
                        ntrees = 500, #default
                        p = 200, #default
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
## It is apparently an issue with how objects are called in the nested predict function (https://github.com/topepo/caret/issues/1412)
## methylClass also applies caret::xgbTree and hence, a forked repo version of caret with a fix described by (https://github.com/CeresBarros/caret/commit/badfb66ada4e95bd00bc5b6511b7107f95144387)
## is being used instead of the native caret package 

require(remotes)

remotes::install_github("CeresBarros/caret", subdir = "pkg/caret") 

#IMPORTANT: Restart R, load the required packages again

library(caret)

library(methylClass)

library(doParallel)

# Check the version

packageVersion("caret") #7.0.2.9001

#eXtreme Gradient Boosting (XGB)-10k

n_cores<-20

clusters <-makeForkCluster(n_cores)

registerDoParallel(clusters)

tryCatch({ XGB_sd_10k<-maintrain( betas.. = as.matrix(train_NMB_sd_10k[,!colnames(train_NMB_sd_10k) %in% "histo"]),
                          y.. = train_NMB_sd_10k$histo,
                          seed = 1234,
                          method = "XGB",
                          modelcv = 5, #default
                          gridsearch = TRUE,
                          max_depth = 6, #default
                          eta = c(0.1,0.3), #default
                          gamma = c(0,0.01), #default
                          min.chwght = 1, #default
                          nrounds = 100, #default
                          colsample_bytree = c(0.2,0.5,0.8,1.0), #same as caret XGB models
                          topfeaturenumber = NULL,
                          subset.CpGs = NULL )

saveRDS(XGB_sd_10k,"~/Thesis_models/Imb_XGB_sd_10k.rds")


XGB_mad_10k<-maintrain( betas.. = as.matrix(train_NMB_mad_10k[,!colnames(train_NMB_mad_10k) %in% "histo"]),
                       y.. = train_NMB$histo,
                       seed = 1234,
                       method = "XGB",
                       modelcv = 5, #default
                       gridsearch = TRUE,
                       max_depth = 6, #default
                       eta = c(0.1,0.3), #default
                       gamma = c(0,0.01), #default
                       min.chwght = 1, #default
                       nrounds = 100, #default
                       colsample_bytree = c(0.2,0.5,0.8,1.0), #same as caret XGB models
                       topfeaturenumber = NULL,
                       subset.CpGs = NULL )

saveRDS(XGB_mad_10k,"~/Thesis_models/Imb_XGB_mad_10k.rds")


#eXtreme Gradient Boosting (XGB) - all probes 

XGB_train_NMB<-maintrain( betas.. = train_NMB[,!colnames(train_NMB) %in% "histo"],
                         y.. = train_NMB$histo,
                         seed = 1234,
                         method = "XGB",
                         modelcv = 5, #default
                         gridsearch = TRUE,
                         max_depth = 6, #default
                         eta = c(0.1,0.3), #default
                         gamma = c(0,0.01), #default
                         min.chwght = 1, #default
                         nrounds = 100, #default
                         colsample_bytree = c(0.2,0.5,0.8,1.0), #same as caret XGB models
                         topfeaturenumber = NULL,
                         subset.CpGs = NULL )


saveRDS(XGB_train_NMB,"~/Thesis_models/Imb_XGB_NMB_allprobes.rds") 

}, error = function(e) {
  
  message("Error pushed:", conditionMessage(e))
  
}, finally = {
  
  stopCluster(clusters)
  
  registerDoSEQ() 
  
  message("Parallel clusters terminated")
}) 


# -------- Testing of the models and Evaluation metrics -------------- #

# Wrapper function for caret models - confusion matrix, Receiver Operating Characteristic (ROC) curve plots
# Area under the Curve (AUC) scores, F1 scores and Kappa scores

test_and_evaluate_model<- function ( model, method, test_data, true_test_labels, test_data_name, model_name,feature_selection_method,
                        save_in = "~/Thesis_models/Imbal_caret_results/") {
  
  dir.create(save_in, recursive = TRUE , showWarnings = FALSE)
  
  # To check required packages
  
  required_pgs<- c("caret","pROC","ggplot2")
  
  missing_pgs<-required_pgs[!sapply(required_pgs, requireNamespace, quietly = TRUE)]
  
  if(length(missing_pgs) > 0){
    message("Installing missing packages for evaluation methods", paste(missing_pgs, collapse = ", "))
    install.packages(missing_pgs, dependencies = TRUE)
  } 
  
  # Check data types
  
  if (!is.data.frame(test_data) && !is.matrix (test_data)) {
    
    test_data<-as.data.frame(test_data)
  }
  
  if (!is.factor(true_test_labels)){
    
    true_test_labels <- as.factor(true_test_labels)
  }
  
  if (nrow(test_data) != length(true_test_labels)){
    
    stop("Check the match between the rownames of the test data and of the ground truth labels")
  }
  
  if (method == "caret") {
  
  model_labels<-model$levels
  test_labels<-levels(true_test_labels)  
  
  } else {
    
    model_labels<- levels(model$mod[[1]]$classes)
    test_labels<-levels(true_test_labels)
  }
  
  if (!setequal(model_labels, test_labels)){
    stop ("The histology levels of the model and the test data DO NOT match")
  }
  
  #---Predictions---
  
  classes<-levels(true_test_labels)
  
  if ( method == "caret") {
  
  pred_labels <- tryCatch(
    predict(model, newdata = test_data)
  ,error = function (e) {
    stop("predict() failed:",conditionMessage(e))
  })
  
  pred_probs <- tryCatch(
    predict(model, newdata = test_data, type = "prob")
  , error = function (e) {
     stop( "predict(prob) failed:", conditionMessage(e))
  })
  
  } else {
    
    pred<- tryCatch(methylClass::mainpredict(newdat = test_data,
                                              mod = model$mod), error = function(e)
                                                stop("methylClass::predict error"))
    
    pred_labels<-as.factor(pred$rawpres)
    
    pred_probs<-as.data.frame(pred$scores)
  }
  
  #---Confusion Matrix---
  
  cm<-tryCatch(
  confusionMatrix(pred_labels, true_test_labels, mode = "everything")
  , error = function(e) {
    stop("confusion matrix encountered error:", conditionMessage(e))
  })
  
  cm_df<-as.data.frame(cm$table)
  
  HISTO_COLOURS<-  c("CLA" = "purple", "LCA" = "pink", "MBEN" = "lightgreen", "DN" = "darkgreen")
  
  cm_plot<-ggplot(cm_df, aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = Reference, alpha = Freq), color = "white") +
    geom_text(aes(label = Freq), size = 4, 
              color = ifelse(cm_df$Freq > max(cm_df$Freq) * 0.5, "white", "black")) +
    scale_fill_manual(values = HISTO_COLOURS, guide = "none") +
    scale_alpha_continuous(range = c(0.15, 1), guide = "none") +
    labs(
      title = paste0("Confusion Matrix for ", model_name,
                     "\n(", ") tested on ", test_data_name),
      x = "True Class", y = "Predicted Class"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          plot.title = element_text(hjust = 0.5)) 
  
  ggsave(file.path(save_in, paste0(model_name,"_", test_data_name, "_confusion_matrix.png")),
         cm_plot, width = 6, height = 5, dpi = 300)
  

  #--- Overall Accuracy
  
  accuracy<-unname(cm$overall["Accuracy"])

  
  #---Kappa
  
  kappa<-unname(cm$overall["Kappa"])
  
  #---Precision and Recall- Class wise and Macro
  
   byclass<-cm$byClass
   
   if(is.null(dim(byclass))){
     precision <- unname(byclass["Precision"])
     recall<- unname(byclass["Recall"])
     sensitivity <-setNames(unname(byclass["Sensitivity"]), classes[2])
     specificity <- setNames(unname(byclass["Specificity"]),classes[2])
   } else {
     precision<-mean(byclass[,"Precision"], na.rm = TRUE)
     recall <- mean( byclass[, "Recall"], na.rm = TRUE)
     sensitivity <- byclass[, "Sensitivity"]
     specificity <- byclass[, "Specificity"]
     names(sensitivity) <- gsub("^Class: ", "", rownames(byclass))
     names(specificity) <- gsub("^Class: ", "", rownames(byclass))
   }
  
  #---ROC and ROC-AUC
   
   # One-vs-Rest class wise ROC and AUC
   
   classes<-levels(true_test_labels)
   
   roc_list<-tryCatch({
     setNames(lapply(classes, function(cls){
       pROC::roc(response = as.numeric(true_test_labels == cls),
                 predictor = pred_prob[[cls]], quiet = TRUE)
     }),classes)
   }, error = function (e) { stop ("ROC Curve error:", conditionMessage(e))
     
     })
   
   class_auc <- sapply(roc_list, function(r) as.numeric(pROC::auc(r)))
   
   # Hand and Till (2001) Multi-class AUC
   
   handtill_obj<-tryCatch(pROC::multiclass.roc(response = true_test_labels,
                                               predictors = pred_probs, quiet = TRUE),
                          error = function(e) stop ("Hand-Till AUC Error:", conditionMessage(e)))
   
   macro_auc<- as.numeric(handtill_obj$auc)
   
  
   #Multi Class ROC curve plot
   
   roc_df<-do.call(rbind,lapply(names(roc_list), function(cls){ 
     
     r<-roc_list[[cls]]
     
     data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities,
                class = paste0(cls, " (AUC", round (pROC::auc(r), 3), ")"))
     }))
   
   roc_plot <- ggplot(roc_df , aes( x = fpr, y = tpr , colour = class)) +
     geom_line(linewidth = 1) +
     geom_abline(linetype = "dashed", color = "grey50") +
     scale_color_manual(values = HISTO_COLOURS) +
     labs(title = paste("One-vs-Rest ROC Curves for ", 
                        model_name,"\n(", feature_selection_method, ") tested on ",test_data_name),
          subtitle = paste0("Macro AUC (Hand-Till) = ", round(macro_auc, 3))) +
     theme_minimal() +
     theme(plot.title = element_text(hjust = 0.5))
   
     ggsave (file.path(save_in, paste0(model_name,"_", test_data_name, "_roc_curves.png")),
             roc_plot, width = 7, height = 5, dpi = 300)
     
     #--- Precision Recall Curves and PR-AUC

     pr_list<-tryCatch({
       setNames(lapply(classes, function(cls) {
         prs<-pROC::coords(roc_list[[cls]], x = "all",
                           ret = c("recall","precision"), transpose = FALSE)
         prs<-prs[!is.na(prs$precision) & !is.na(prs$recall),]
         prs<-prs[order(prs$recall),]
         prs}), classes)
     }, error = function (e) {
       
       stop ("PR Curve Error:", conditionMessage(e))
     })

     pr_auc <- tryCatch({
       sapply(classes, function(cls){
         MLmetrics::PRAUC(y_pred = pred_probs[[cls]],
                          y_true = as.numeric(true_test_labels == cls))
       })
     }, error = function(e){ stop ("PRAUC Error:", conditionMessage(e))
       })
     
     names(pr_auc)<-classes
     
     pr_df <- do.call(rbind, lapply(classes, function(cls) {
       d <- pr_list[[cls]]
       data.frame(recall = d$recall, precision = d$precision,
                  class = paste0(cls, " (PR-AUC=", round(pr_auc[[cls]], 3), ")"))
     }))
     
     #--Macro PR-AUC
     
     macro_pr_auc<-mean(pr_auc, na.rm = TRUE)
     
     pr_plot <- ggplot(pr_df, aes(x = recall, y = precision, color = class)) +
       geom_line(linewidth = 1) +
       scale_color_manual(values = HISTO_COLOURS) +
       labs(title = paste("One-vs-Rest Precision-Recall Curves for", model_name, "\n(", feature_selection_method, ") tested on " , test_data_name),
            subtitle = paste0("Macro PR-AUC = ", round(macro_pr_auc, 3)),
            x = "Recall", y = "Precision", color = "Class") +
       theme_minimal() + ylim(0, 1) + xlim(0, 1) + theme(plot.title = element_text(hjust = 0.5))
     
     ggsave(file.path(out_dir, paste0(model_name, "_pr_curves.png")),
            pr_plot, width = 7, height = 5, dpi = 300)
    


    #---Summary
     
     summary_df<-data.frame(
       model = model_name ,
       test_data = test_data_name,
       macro_precision = precision ,
       macro_recall = recall,
       macro_auc = macro_auc,
       macro_pr_auc = macro_pr_auc
     )


     #---Class-wise metrics
     
     classwise_df <- data.frame(
       model = model_name,
       test_data = test_data_name,
       sensitivity = as.numeric(sensitivity[classes]),
       specificity = as.numeric(specificity[classes]),
       auc = as.numeric(class_auc[classes]),
       pr_auc = as.numeric(pr_auc[classes])
     )
     
     #---Save as .csv files
     
     write.csv(summary_df, file.path(save_in, paste0(model_name,"_", test_data_name, "_summary_metrics.csv")), row.names = FALSE)
     write.csv(classwise_df, file.path(save_in, paste0(model_name, "_", test_data_name,"_classwise_metrics.csv")), row.names = FALSE)
     
     #---To Console
  
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
       classwise        = classwise_df
     )
}

# Load models if not in the env.

RF_sd<-readRDS("~/Thesis_models/Imb_RF_sd_10k.rds")

SVM_sd<-readRDS("~/Thesis_models/Imb_SVM_sd_10k.rds")

XGB_sd<-readRDS("~/Thesis_models/Imb_XGB_sd_10k.rds") 

RF_mad<-readRDS("~/Thesis_models/Imb_RF_mad_10k.rds")

SVM_mad<-readRDS("~/Thesis_models/Imb_SVM_mad_10k.rds")

XGB_mad<-readRDS("~/Thesis_models/Imb_XGB_mad_10k.rds") 

RF_all<-readRDS("~/Thesis_models/Imb_RF_NMB_allprobes.rds")

SVM_all<-readRDS("~/Thesis_models/Imb_SVM_NMB_allprobesrds")

XGB_all<-readRDS("~/Thesis_models/Imb_XGB_NMB_allprobes.rds") 


models<-list(RF_sd = RF_sd,
             SVM_sd = SVM_sd,
             #XGB_sd = XGB_sd,
             RF_mad = RF_mad,
             SVM_mad = SVM_mad,
             #XGB_mad = XGB_mad,
             RF_all = RF_all,
             SVM_all = SVM_all,
             XGB_all = XGB_all,
             )

# Setting titles for plots and results

get_model_desc<-function(model_key){
  
  algo<-if ( grepl("^RF", model_key, ignore.case = TRUE)){
    model_name<-"Random Forest"
  } else if (grepl("^SVM", model_key, ignore.case = TRUE)) {
      model_name<-"Support Vector Machine"
    } else if ( grepl("^XGB", model_key, ignore.case = TRUE)){
      model_name<-"eXtreme Gradient Boosting"
    } else {
      model_key
    }

feature_selection<-if(grepl("sd",model_key, ignore.case = TRUE)){
  "Top 10k variable CpG probes by Standard Deviation"
} else if (grepl("mad", model_key, ignore.case = TRUE)){
  "Top 10k variable CpG probes by Median Absolute Deviation"
} else {
  "all_probes"
}

balancing_method<- if (grepl("up", model_key,ignore.case = TRUE)){
  "Up-sampled minority classes"
} else if (grepl("down",model_key, ignore.case = TRUE)){
  "Down-sampled CLA"
} else if (grepl("SMOTE", model_key, ignore.case = FALSE)){
  "SMOTE"
} else {
  NA
}

list(algo = algo, feature_selection = feature_selection, balancing_method = balancing_method)

} 


# Testing on test NMB 

test_NMB_results<-list()

model_dict<-data.frame(model = character(), model_label = character())

for (model_key in names(models)){
  
  label<-get_model_desc(model_key)
  
  message("Testing and evaluating", label, "on Test NMB")
  
 if(all(rownames(test_NMB) == NMB_test_pheno$Sample_ID)){
   tryCatch( test_NMB_results[[model_key]]<-test_and_evaluate_model(
    model = models[[model_key]],
    test_data = test_NMB,
    test_data_name = "Test NMB",
    true_test_labels = NMB_test_pheno$CPR_Histology,
    model_name = label[["algo"]],
    feature_selection_method = label[["feature_selection"]],
    save_in = "~/Thesis_models/test_NMB_imbal/"
    ),
   error = function(e) {
     message(sprintf("FAILED [%s / %s]: %s",
                     "Binary Random Forest trained on Train NMB", "Test NMB", conditionMessage(e)))
     NULL
   } )
   
 } else {
   stop("Mismatch in test set and test pheno Sample IDs")
 }
 
  test_NMB_results[[model_key]]$summary$model_label <- label
  
  test_NMB_results[[model_key]]$classwise$model_label<-label
  
  model_dict<-rbind(model_dict, data.frame(model=model_key, model_label = label))
}

final_test_NMB_summary<-do.call(rbind, lapply(test_NMB_results, function(x) x$summary))

final_test_NMB_classwise_summary<-do.call(rbind, lapply(test_NMB_results, function(x) x$classwise))

print(final_test_NMB_summary)

# Comparison Bar plot 

final_test_NMB_summary$algorithm <- ifelse(grepl("^RF", final_test_NMBsummary$model), "RF",
                                           ifelse(grepl("^SVM", final_test_NMB_summary$model), "SVM", "XGB"))
final_test_NMB_summary$feature_selection <- ifelse(grepl("sd", final_test_NMB_ummary$model, ignore.case = TRUE),
                                                   "SD", ifelse(grepl("^mad", final_test_NMB_summary), ignore.case = TRUE), "MAD", "All_probes")


p_grouped_accuracy <- ggplot(final_test_NMB_summary, aes(x = algorithm, y = accuracy, fill = feature_selection)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(accuracy, 3)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 3.5) +
  labs(title = "Accuracy Comparison by Machine Learning Model and Feature Selection Method in Test NMB",
       x = "Algorithm", y = "Accuracy score", fill = "Feature Selection Method") +
  ylim(0, 1) +
  theme_minimal()

ggsave("~/Thesis_models/test_NMB_imbal/accuracy_comparision.png", width = 7, height = 5, dpi = 300)


p_grouped_auc <- ggplot(final_test_NMB_summary, aes(x = algorithm, y = auc, fill = feature_selection)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(accuracy, 3)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 3.5) +
  labs(title = "Area Under the Curve (AUC) scores Comparison by Machine Learning Model and Feature Selection Method in Test NMB",
       x = "Algorithm", y = "AUC score", fill = "Feature Selection Method") +
  ylim(0, 1) +
  theme_minimal() 

ggsave("~/Thesis_models/test_NMB_imbal/auc_comparision.png", width = 7, height = 5, dpi = 300)

# Testing on Cavalli

Cavalli_results<-list()

model_dict<-data.frame(model = character(), model_label = character())

for (model_key in names(models)){
  
  label<-get_model_desc(model_key)
  
  message("Testing and evaluating", label, "on Cavalli")
  
  if(all(rownames(Cavalli_betas) == Cavalli_pheno$Study_ID...12)){
    Cavalli_results[[model_key]]<-test_and_evaluate_model(
      model = models[[model_key]],
      test_data = Cavalli_betas,
      test_data_name = "Cavalli",
      true_test_labels = Cavalli_pheno$Study_ID...12,
      model_name = label[["algo"]],
      feature_selection_method = label[["feature_selection"]],
      save_in = "~/Thesis_models/Cavalli_imbal/"
    )
  } else {
    stop("Mismatch in test set and test pheno Sample IDs")
  }
  
 Cavalli_results[[model_key]]$summary$model_label <- label
  
  Cavalli_results[[model_key]]$classwise$model_label<-label
  
  model_dict<-rbind(model_dict, data.frame(model=model_key, model_label = label))
}

final_Cavalli_summary<-do.call(rbind, lapply(Cavalli_results, function(x) x$summary))

final_Cavalli_classwise_summary<-do.call(rbind, lapply(Cavalli_results, function(x) x$classwise))

print(final_Cavalli_summary)

# Comparison Bar plot 

final_Cavalli_summary$algorithm <- ifelse(grepl("^RF", final_Cavalli_summary$model, ignore.case = TRUE), "RF",
                                           ifelse(grepl("^SVM", final_Cavalli_summary$model, ignore.case = TRUE), "SVM", "XGB"))
final_Cavalli_summary$feature_selection <- ifelse(grepl("sd", final_Cavalli_summary$model, ignore.case = TRUE),
                                                   "SD", ifelse(grepl("^mad", final_Cavalli_summary), ignore.case = TRUE), "MAD", "All_probes")


p_grouped_accuracy <- ggplot(final_Cavalli_summary, aes(x = algorithm, y = accuracy, fill = feature_selection)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(accuracy, 3)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 3.5) +
  labs(title = "Accuracy Comparison by Machine Learning Model and Feature Selection Method in Test-Cavalli",
       x = "Algorithm", y = "Accuracy score", fill = "Feature Selection Method") +
  ylim(0, 1) +
  theme_minimal()

ggsave("~/Thesis_models/Cavalli_imbal/accuracy_comparision.png", width = 7, height = 5, dpi = 300)


p_grouped_auc <- ggplot(final_Cavalli_summary, aes(x = algorithm, y = auc, fill = feature_selection)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(accuracy, 3)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 3.5) +
  labs(title = "Area Under the Curve (AUC) scores Comparison by Machine Learning Model and Feature Selection Method in Test-Cavalli",
       x = "Algorithm", y = "AUC score", fill = "Feature Selection Method") +
  ylim(0, 1) +
  theme_minimal() 

ggsave("~/Thesis_models/Cavalli_imbal/auc_comparision.png", width = 7, height = 5, dpi = 300)


# ------- Balancing the train set --------------- #

required_pgs<-c("caret","themis","methylClass","recipes")

missing_pgs<-required_pgs[!sapply(required_pgs, requireNamespace, quietly = TRUE)]

if(length(missing_pgs) > 0){
  message("Installing missing packages for sampling methods:", paste(missing_pgs, collapse = ", "))
  install.packages(missing_pgs, dependencies = TRUE)
} 


# ---- Wrapper function to apply sampling methods ---- #

balance<-function(train_df, outcome_labels, method = c("none","up","down","methylSMOTE",
                                                       "smote","adasyn","tomek","smote_tomek"),
                  over_ratio = 1, 
                  minority_threshold = 1, 
                  neighbours = 5, 
                  downsampling = FALSE, 
                  seed = 123) {

#MethylSMOTE- SMOTE with methylClass::balancesampling
  
#smote - themis recipe

#over_ratio for themis based sampling methods 
  
#minority_threshold (cutoff in methylClass::balancesampling) is the number of samples in the majority class, 
#all the classes with sample numbers below this will be over sampled with SMOTE 
  
#Number of nearest neighbours for SMOTE to interpolate with. 
  
#downsampling for methylClass::balancesampling is set to TRUE in case the minority_threshold != n in majority class
#to reduce the sample numbers down to the number set with minority_threshold

  method<-match.arg(method)
  
  stopifnot(is.factor(train_df[[outcome_labels]]))
  
  set.seed(seed)
  
  if (method == "none") {
    
    return(train_df)
    
    #---Random Oversampling/Undersampling
    
    if(method %in% c("up","down")){
      
      x<-train_df[,setdiff(names(train_df), outcome_labels), drop = FALSE]
      y<-train_df[[outcome_labels]]
      
      out<-if (method == "up"){
        
        caret::upSample(x = x, y = y, yname = outcome_labels)
      } else {
        caret::downSample(x = x, y = y, yname = outcome_labels)
      }
      
      return(out)
    }
    
  }
  
  if (method == "methylSMOTE") {
    
    x<-train_df[,setdiff(names(train_df), outcome_labels), drop = FALSE]
    y<-train_df[[outcome_labels]]
    
    out<-methylClass::balancesampling(x, labels = y, 
                                      topfeaturenumber = NULL, # change if selecting top features
                                      cutoff = minority_threshold,
                                      downsampling = FALSE,
                                      sampleseed = seed,
                                      k = neighbours,
                                      adjustbetas = TRUE ) 
    
    
  }
  
  if(method %in% c("smote","adasyn","tomek","smote_tomek")){
  
  rec<-switch( method, 
               smote = rec %>%
                 step_smote(all_of(train_df[[outcome_labels]]), over_ratio = over_ratio,
                            neighbors = neighbours, seed = seed),
               
               adasyn = rec %>% 
                 step_adasyn(all_of(train_df[[outcome_labels]]), over_ratio = over_ratio,
                             neighbors = neighbours, seed = seed),
               
               
               tomek = rec %>%
                 step_tomek (all_of(train_df[[outcome_labels]]),
               neighbors = neighbours, seed = seed),
               
               smote_tomek = rec %>%
                 step_smote(all_of(train_df[[outcome_labels]]),
                            neighbors = neighbors, seed = seed) %>%
                 step_tomek(all_of(train_df[[outcome_labels]]))
  )
  
  rec_prep<-prep(rec, training = train_df, retain = TRUE)
  
  bake(rec_prep, new_data = NULL) }
      
} 

#Save_in

dir.create("~/Thesis_models/balanced_sets/", showWarnings = TRUE, recursive = TRUE)

balancing_methods<-c("methylSMOTE",
           "smote","adasyn","tomek","smote_tomek") #"up","down"

balanced_sets<-list()

for (method in balancing_methods) {
  
  balanced_train_NMB<-balance ( train_df = train_NMB,
                              outcome_labels = "histo",
                              method = method,
                              minority_threshold = max(table(train_NMB$histo)),
                              neighbours = 5, #default in most smote applications
                              over_ratio = 1,
                              seed = 42 ) 
  
  saveRDS(balanced_train_NMB,paste0("~/Thesis_models/balanced_sets/", method,"_train_NMB.rds"))
  
  balanced_sets[[method]]<-balanced_train_NMB
 } 

#---------Training ML models on balanced sets -----------#

# Load balanced data sets if not in the env.

get_balanced_set <- function(method, dir = "~/Thesis_models/balanced_sets/"){
  readRDS(file.path(dir,paste0("_",method,"_train_NMB.rds")))
}

balancing_methods<-list("up","down","methylSMOTE",
              "smote","adasyn","tomek","smote_tomek")

balanced_sets<-list()

for (method in names(balancing_methods)){
  balanced_sets[[method]]<-get_balanced_set(method,
                          dir = "~/Thesis_models/balanced_sets/")
}

# Training and testing methylClass models on each of the balanced Train NMB set with all probes

balanced_models<-list()

for (set in names(balanced_sets)){
    
    train_set<-balanced_sets[[set]]
    
    x <- train_set[,setdiff(names(train_set), "histo"), drop = FALSE]
    
    y<-as.factor(train_set[["histo"]])
    
    RF_model<-methylClass::maintrain(betas.. = x, y.. = y ,
                                     method = "RF",
                                     subset.CpGs = NULL,
                                     topfeaturenumber = NULL,
                                     seed = 42,
                                     ntrees = 500, #default
                                     p = 200 )
    
    saveRDS(RF_model, paste0("~/Thesis_models/balanced_sets/","RF_",set,"_all_probes.rds"))
  
    
    SVM_model <-maintrain( betas.. = x,
                              y.. = y,
                              method = "SVM",
                              modelcv = 5, #default
                              gridsearch = TRUE,
                              C.base = 10, #default
                              C.min = -3, #default
                              C.max = -2, #default
                              topfeaturenumber = NULL,
                              subset.CpGs = NULL )
    
    
    saveRDS(SVM_model, paste0("~/Thesis_models/balanced_sets/","SVM_",set,"_all_probes.rds"))
    
    #eXtreme Gradient Boosting (XGB)
    
    XGB_model<-maintrain( betas.. = x,
                              y.. = y,
                              method = "XGB",
                              modelcv = 5, #default
                              gridsearch = TRUE,
                              max_depth = 6, #default
                              eta = c(0.1,0.3), #default
                              gamma = c(0,0.01), #default
                              min.chwght = 1, #default
                              nrounds = 100, #default
                              colsample_bytree = c(0.2,0.5,0.8,1.0), #same as caret XGB models
                              topfeaturenumber = NULL,
                              subset.CpGs = NULL )
    
    
    saveRDS(XGB_model, paste0("~/Thesis_models/balanced_sets/","XGB_",set,"_all_probes.rds"))
    
  balanced_models[[set]]<-list(rf = RF_model, svm = SVM_model, xgb = XGB_model)
}

#---Load models if not in the env.


algos<-c("RF","SVM","XGB")

load_balanced_models<-function(algo, model_dir){
  files<-list.files( path.expand(model_dir),
                     pattern = paste0("^", algo,"_.*_all_probes\\.rds$"),
                     full.names = TRUE)
  method_names<-gsub(paste0("^", algo, "_|_all_probes.rds$"), "", basename(files))
  setnames(lapply(files,readRDS), method_names)
}

balanced_models<-setNames(lapply(algos, load_balanced_models, model_dir = "~Thesis_models/balanced_sets/"), algos)


# Testing on test NMB 

test_NMB_results<-list()

model_dict<-data.frame(model = character(), model_label = character())

for (model_key in names(balanced_models)){
  
  label<-get_model_desc(model_key)
  
  message("Testing and evaluating", label, "on Test NMB")
  
  if(all(rownames(test_NMB) == NMB_test_pheno$Sample_ID)){
    test_NMB_results[[model_key]]<-test_and_evaluate_model(
      model = models[[model_key]],
      test_data = test_NMB,
      test_data_name = "Test NMB",
      true_test_labels = NMB_test_pheno$CPR_Histology,
      model_name = label[["algo"]],
      feature_selection_method = label[["feature_selection"]],
      save_in = "~/Thesis_models/test_NMB_bal/"
    )
  } else {
    stop("Mismatch in test set and test pheno Sample IDs")
  }
  
  test_NMB_results[[model_key]]$summary$model_label <- label
  
  test_NMB_results[[model_key]]$classwise$model_label<-label
  
  model_dict<-rbind(model_dict, data.frame(model=model_key, model_label = label))
}

final_test_NMB_summary<-do.call(rbind, lapply(test_NMB_results, function(x) x$summary))

final_test_NMB_classwise_summary<-do.call(rbind, lapply(test_NMB_results, function(x) x$classwise))

print(final_test_NMB_summary)

# Comparison Bar plot 

final_test_NMB_summary$algorithm <- ifelse(grepl("^RF", final_test_NMBsummary$model, ignore.case = TRUE), "RF",
                                           ifelse(grepl("^SVM", final_test_NMB_summary$model, ignore.case = TRUE), "SVM", "XGB"))

# final_test_NMB_summary$feature_selection <- ifelse(grepl("^sd", final_test_NMB_summary$model, ignore.case = TRUE),
#                                                    "SD", ifelse(grepl("^mad", final_test_NMB_summary), ignore.case = TRUE), "MAD", "All_probes")

final_test_NMB_summary$balanced <- ifelse(
  grepl("smote_tomek", final_test_NMB_summary$model_label, ignore.case = TRUE), "SMOTE-Tomek",
  ifelse(grepl("^up",    final_test_NMB_summary$model_label, ignore.case = TRUE), "Up-sampled",
         ifelse(grepl("^down",  final_test_NMB_summary$model_label, ignore.case = TRUE), "Down-sampled",
                ifelse(grepl("^smote", final_test_NMB_summary$model_label, ignore.case = TRUE), "SMOTE",
                       ifelse(grepl("tomek",  final_test_NMB_summary$model_label, ignore.case = TRUE), "Tomek",
                              "Unbalanced"))))) 

p_grouped_accuracy <- ggplot(final_test_NMB_summary, aes(x = algorithm, y = accuracy, fill = balanced)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(accuracy, 3)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 3.5) +
  labs(title = "Accuracy Comparison by Machine Learning Model and Balancing Method in Test NMB",
       x = "Algorithm", y = "Accuracy score", fill = "Balancing Method") +
  ylim(0, 1) +
  theme_minimal()

ggsave("~/Thesis_models/test_NMB_bal/accuracy_comparision.png", width = 7, height = 5, dpi = 300)


p_grouped_auc <- ggplot(final_test_NMB_summary, aes(x = algorithm, y = auc, fill = balanced)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(accuracy, 3)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 3.5) +
  labs(title = "Area Under the Curve (AUC) scores Comparison by Machine Learning Model and Balancing Method in Test NMB",
       x = "Algorithm", y = "AUC score", fill = "Balancing Method") +
  ylim(0, 1) +
  theme_minimal() 

ggsave("~/Thesis_models/test_NMB_bal/auc_comparision.png", width = 7, height = 5, dpi = 300)

# Testing on Cavalli

Cavalli_results<-list()

model_dict<-data.frame(model = character(), model_label = character())

for (model_key in names(models)){
  
  label<-get_model_desc(model_key)
  
  message("Testing and evaluating", label, "on Cavalli")
  
  if(all(rownames(Cavalli_betas) == Cavalli_pheno$Study_ID...12)){
    Cavalli_results[[model_key]]<-test_and_evaluate_model(
      model = models[[model_key]],
      test_data = Cavalli_betas,
      test_data_name = "Cavalli",
      true_test_labels = Cavalli_pheno$Study_ID...12,
      model_name = label[["algo"]],
      feature_selection_method = label[["feature_selection"]],
      save_in = "~/Thesis_models/Cavalli_bal/"
    )
  } else {
    stop("Mismatch in test set and test pheno Sample IDs")
  }
  
  Cavalli_results[[model_key]]$summary$model_label <- label
  
  Cavalli_results[[model_key]]$classwise$model_label<-label
  
  model_dict<-rbind(model_dict, data.frame(model=model_key, model_label = label))
}

final_Cavalli_summary<-do.call(rbind, lapply(Cavalli_results, function(x) x$summary))

final_Cavalli_classwise_summary<-do.call(rbind, lapply(Cavalli_results, function(x) x$classwise))

print(final_Cavalli_summary)

# Comparison Bar plot 

final_Cavalli_summary$algorithm <- ifelse(grepl("^RF", final_Cavalli_summary$model, ignore.case = TRUE), "RF",
                                          ifelse(grepl("^SVM", final_Cavalli_summary$model, ignore.case = TRUE), "SVM", "XGB"))
# final_Cavalli_summary$feature_selection <- ifelse(grepl("sd", final_Cavalli_summary$model, ignore.case = TRUE),
#                                                   "SD", ifelse(grepl("^mad", final_Cavalli_summary), ignore.case = TRUE), "MAD", "All_probes")

final_Cavalli_summary$balanced <- ifelse(
  grepl("smote_tomek", final_Cavalli_summary$model_label, ignore.case = TRUE), "SMOTE-Tomek",
  ifelse(grepl("^up",    final_Cavalli_summary$model_label, ignore.case = TRUE), "Up-sampled",
         ifelse(grepl("^down",  final_Cavalli_summary$model_label, ignore.case = TRUE), "Down-sampled",
                ifelse(grepl("^smote", final_Cavalli_summary$model_label, ignore.case = TRUE), "SMOTE",
                       ifelse(grepl("tomek",  final_Cavalli_summary$model_label, ignore.case = TRUE), "Tomek",
                              "Unbalanced"))))) 

p_grouped_accuracy <- ggplot(final_Cavalli_summary, aes(x = algorithm, y = accuracy, fill = balanced)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(accuracy, 3)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 3.5) +
  labs(title = "Accuracy Comparison by Machine Learning Model and Balancing Method in Test-Cavalli",
       x = "Algorithm", y = "Accuracy score", fill = "Balancing Method") +
  ylim(0, 1) +
  theme_minimal()

ggsave("~/Thesis_models/Cavalli_bal/accuracy_comparision.png", width = 7, height = 5, dpi = 300)


p_grouped_auc <- ggplot(final_Cavalli_summary, aes(x = algorithm, y = auc, fill = balanced)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(accuracy, 3)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 3.5) +
  labs(title = "Area Under the Curve (AUC) scores Comparison by Machine Learning Model and Balancing Method in Test-Cavalli",
       x = "Algorithm", y = "AUC score", fill = "Balancing Method") +
  ylim(0, 1) +
  theme_minimal() 

ggsave("~/Thesis_models/Cavalli_bal/auc_comparision.png", width = 7, height = 5, dpi = 300)

# Training and Testing methylClass models on each of the balanced Train NMB set with top 10k probes (sub setted before applying balancing methods)
# Using methylClass to allow quick run in comparision to computationally expensive caret models.

# Standard Deviation

sd_10k_balanced_sets<-list()

for (set in names(balanced_sets)) {
  
  train_set<-as.data.frame(balanced_sets[[set]])
  
  x<-train_set[,colnames(train_set) %in% colnames(train_NMB_sd_10k)]  #train_NMB_sd_10k generated earlier in the script
  
  if (all(colnames(x) == colnames(train_NMB_sd_10k))) {
    saveRDS(x, paste0("~/Thesis_models/balanced_sets/",set,"_10k_sd.rds"))
  } else {
    stop("Error in creating SD 10k balanced sets")
  }
  
}

sd_balanced_models<-list()

for (set in names(sd_10k_balanced_sets)){
  
  train_set<-balanced_sets[[set]]
  
  x <- train_set[,setdiff(names(train_set), "histo"), drop = FALSE]
  
  y<-as.factor(train_set[["histo"]])
  
  RF_model<-methylClass::maintrain(betas.. = x, y.. = y ,
                                   method = "RF",
                                   subset.CpGs = NULL,
                                   topfeaturenumber = NULL,
                                   seed = 42,
                                   ntrees = 500, #default
                                   p = 200 )
  
  saveRDS(RF_model, paste0("~/Thesis_models/balanced_sets/","RF_",set,"_sd_10k.rds"))
  
  
  SVM_model <-maintrain( betas.. = x,
                         y.. = y,
                         method = "SVM",
                         modelcv = 5, #default
                         gridsearch = TRUE,
                         C.base = 10, #default
                         C.min = -3, #default
                         C.max = -2, #default
                         topfeaturenumber = NULL,
                         subset.CpGs = NULL )
  
  
  saveRDS(SVM_model, paste0("~/Thesis_models/balanced_sets/","SVM_",set,"_sd_10k.rds"))
  
  #eXtreme Gradient Boosting (XGB)
  
  XGB_model<-maintrain( betas.. = x,
                        y.. = y,
                        method = "XGB",
                        modelcv = 5, #default
                        gridsearch = TRUE,
                        max_depth = 6, #default
                        eta = c(0.1,0.3), #default
                        gamma = c(0,0.01), #default
                        min.chwght = 1, #default
                        nrounds = 100, #default
                        colsample_bytree = c(0.2,0.5,0.8,1.0), #same as caret XGB models
                        topfeaturenumber = NULL,
                        subset.CpGs = NULL )
  
  
  saveRDS(XGB_model, paste0("~/Thesis_models/balanced_sets/","XGB_",set,"_sd_10k.rds"))
  
 sd_balanced_models[[set]]<-list(rf = RF_model, svm = SVM_model, xgb = XGB_model)
}

#--Load the models if not in the env.

algos<-c("RF","SVM","XGB")

load_balanced_models<-function(algo, model_dir){
  files<-list.files( path.expand(model_dir),
                     pattern = paste0("^", algo,"_.*_sd_10k\\.rds$"),
                     full.names = TRUE)
  method_names<-gsub(paste0("^", algo, "_|_sd_10k\\.rds$"), "", basename(files))
  setnames(lapply(files,readRDS), method_names)
}

sd_balanced_models<-setNames(lapply(algos, load_balanced_models, model_dir = "~Thesis_models/balanced_sets/"), algos)


# Testing on Test NMB

test_NMB_results<-list()

model_dict<-data.frame(model = character(), model_label = character())

for (model_key in names(sd_balanced_models)){
  
  label<-get_model_desc(model_key)
  
  message("Testing and evaluating", label, "on Test NMB")
  
  if(all(rownames(test_NMB) == NMB_test_pheno$Sample_ID)){
    test_NMB_results[[model_key]]<-test_and_evaluate_model(
      model = models[[model_key]],
      test_data = test_NMB,
      test_data_name = "Test NMB",
      true_test_labels = NMB_test_pheno$CPR_Histology,
      model_name = label[["algo"]],
      feature_selection_method = label[["feature_selection"]],
      save_in = "~/Thesis_models/test_NMB_bal/"
    )
  } else {
    stop("Mismatch in test set and test pheno Sample IDs")
  }
  
  test_NMB_results[[model_key]]$summary$model_label <- label
  
  test_NMB_results[[model_key]]$classwise$model_label<-label
  
  model_dict<-rbind(model_dict, data.frame(model=model_key, model_label = label))
}

final_test_NMB_summary<-do.call(rbind, lapply(test_NMB_results, function(x) x$summary))

final_test_NMB_classwise_summary<-do.call(rbind, lapply(test_NMB_results, function(x) x$classwise))

print(final_test_NMB_summary)

# Comparison Bar plot 

final_test_NMB_summary$algorithm <- ifelse(grepl("^RF", final_test_NMBsummary$model, ignore.case = TRUE), "RF",
                                           ifelse(grepl("^SVM", final_test_NMB_summary$model, ignore.case = TRUE), "SVM", "XGB"))

# final_test_NMB_summary$feature_selection <- ifelse(grepl("^sd", final_test_NMB_summary$model, ignore.case = TRUE),
#                                                    "SD", ifelse(grepl("^mad", final_test_NMB_summary), ignore.case = TRUE), "MAD", "All_probes")

final_test_NMB_summary$balanced <- ifelse(
  grepl("smote_tomek", final_test_NMB_summary$model_label, ignore.case = TRUE), "SMOTE-Tomek",
  ifelse(grepl("^up",    final_test_NMB_summary$model_label, ignore.case = TRUE), "Up-sampled",
         ifelse(grepl("^down",  final_test_NMB_summary$model_label, ignore.case = TRUE), "Down-sampled",
                ifelse(grepl("^smote", final_test_NMB_summary$model_label, ignore.case = TRUE), "SMOTE",
                       ifelse(grepl("tomek",  final_test_NMB_summary$model_label, ignore.case = TRUE), "Tomek",
                              "Unbalanced"))))) 

p_grouped_accuracy <- ggplot(final_test_NMB_summary, aes(x = algorithm, y = accuracy, fill = balanced)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(accuracy, 3)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 3.5) +
  labs(title = "Accuracy Comparison by Machine Learning Model trained on Top 10k probes (Standard Deviation) and Balancing Method in Test NMB",
       x = "Algorithm", y = "Accuracy score", fill = "Balancing Method") +
  ylim(0, 1) +
  theme_minimal()

ggsave("~/Thesis_models/test_NMB_bal/sd_10k_accuracy_comparision.png", width = 7, height = 5, dpi = 300)


p_grouped_auc <- ggplot(final_test_NMB_summary, aes(x = algorithm, y = auc, fill = balanced)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(accuracy, 3)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 3.5) +
  labs(title = "Area Under the Curve (AUC) score Comparison by Machine Learning Model trained on Top 10k probes (Standard Deviation) and Balancing Method in Test NMB",
       x = "Algorithm", y = "Accuracy score", fill = "Balancing Method") +
  ylim(0, 1) +
  theme_minimal()

ggsave("~/Thesis_models/test_NMB_bal/sd_10k_auc_comparision.png", width = 7, height = 5, dpi = 300)


# Testing on Cavalli

Cavalli_results<-list()

model_dict<-data.frame(model = character(), model_label = character())

for (model_key in names(sd_balanced_models)){
  
  label<-get_model_desc(model_key)
  
  message("Testing and evaluating", label, "on Test NMB")
  
  if(all(rownames(Cavalli_betas) == Cavalli_pheno$Study_ID...12)){
    test_NMB_results[[model_key]]<-test_and_evaluate_model(
      model = models[[model_key]],
      test_data = Cavalli_betas,
      test_data_name = "Cavalli",
      true_test_labels = Cavalli_pheno$histology,
      model_name = label[["algo"]],
      feature_selection_method = label[["feature_selection"]],
      save_in = "~/Thesis_models/test_NMB_bal/"
    )
  } else {
    stop("Mismatch in test set and test pheno Sample IDs")
  }
  
  Cavalli_results[[model_key]]$summary$model_label <- label
  
  Cavalli_results[[model_key]]$classwise$model_label<-label
  
  model_dict<-rbind(model_dict, data.frame(model = model_key, model_label = label))
}

final_Cavalli_summary<-do.call(rbind, lapply(Cavalli_results, function(x) x$summary))

final_Cavalli_classwise_summary<-do.call(rbind, lapply(Cavalli_results, function(x) x$classwise))

print(final_Cavalli_summary)

# Comparison Bar plot 

final_Cavalli_summary$algorithm <- ifelse(grepl("^RF", final_Cavallisummary$model, ignore.case = TRUE), "RF",
                                           ifelse(grepl("^SVM", final_Cavalli_summary$model, ignore.case = TRUE), "SVM", "XGB"))

# final_Cavalli_summary$feature_selection <- ifelse(grepl("^sd", final_Cavalli_summary$model, ignore.case = TRUE),
#                                                    "SD", ifelse(grepl("^mad", final_Cavalli_summary), ignore.case = TRUE), "MAD", "All_probes")

final_Cavalli_summary$balanced <- ifelse(
  grepl("smote_tomek", final_Cavalli_summary$model_label, ignore.case = TRUE), "SMOTE-Tomek",
  ifelse(grepl("^up",    final_Cavalli_summary$model_label, ignore.case = TRUE), "Up-sampled",
         ifelse(grepl("^down",  final_Cavalli_summary$model_label, ignore.case = TRUE), "Down-sampled",
                ifelse(grepl("^smote", final_Cavalli_summary$model_label, ignore.case = TRUE), "SMOTE",
                       ifelse(grepl("tomek",  final_Cavalli_summary$model_label, ignore.case = TRUE), "Tomek",
                              "Unbalanced"))))) 

p_grouped_accuracy <- ggplot(final_Cavalli_summary, aes(x = algorithm, y = accuracy, fill = balanced)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(accuracy, 3)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 3.5) +
  labs(title = "Accuracy Comparison by Machine Learning Model trained on Top 10k probes (Standard Deviation) and Balancing Method in Cavalli",
       x = "Algorithm", y = "Accuracy score", fill = "Balancing Method") +
  ylim(0, 1) +
  theme_minimal()

ggsave("~/Thesis_models/Cavalli_bal/sd_10k_accuracy_comparision.png", width = 7, height = 5, dpi = 300)


p_grouped_auc <- ggplot(final_Cavalli_summary, aes(x = algorithm, y = auc, fill = balanced)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(accuracy, 3)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 3.5) +
  labs(title = "Area Under the Curve (AUC) score Comparison by Machine Learning Model trained on Top 10k probes (Standard Deviation) and Balancing Method in Cavalli",
       x = "Algorithm", y = "Accuracy score", fill = "Balancing Method") +
  ylim(0, 1) +
  theme_minimal()

ggsave("~/Thesis_models/Cavalli_bal/sd_10k_auc_comparision.png", width = 7, height = 5, dpi = 300)




# Median Absolute Deviation

mad_10k_balanced_sets<-list()

for (set in names(balanced_sets)) {
  
  train_set<-balanced_sets[[set]]
  
  x<-train_set[,colnames(train_set) %in% colnames(train_NMB_mad_10k)]  #train_NMB_mad_10k generated earlier in the script
  
  if (all(colnames(x) == colnames(train_NMB_mad_10k))) {
    saveRDS(x, paste0("~/Thesis_models/balanced_sets/",set,"_10k_mad.rds"))
  } else {
    stop("Error in creating MAD 10k balanced sets")
  }
  
}


mad_balanced_models<-list()

for (set in names(mad_10k_balanced_sets)){
  
  train_set<-balanced_sets[[set]]
  
  x <- train_set[,setdiff(names(train_set), "histo"), drop = FALSE]
  
  y<-as.factor(train_set[["histo"]])
  
  RF_model<-methylClass::maintrain(betas.. = x, y.. = y ,
                                   method = "RF",
                                   subset.CpGs = NULL,
                                   topfeaturenumber = NULL,
                                   seed = 42,
                                   ntrees = 500, #default
                                   p = 200 )
  
  saveRDS(RF_model, paste0("~/Thesis_models/balanced_sets/","RF_",set,"_mad_10k.rds"))
  
  
  SVM_model <-maintrain( betas.. = x,
                         y.. = y,
                         method = "SVM",
                         modelcv = 5, #default
                         gridsearch = TRUE,
                         C.base = 10, #default
                         C.min = -3, #default
                         C.max = -2, #default
                         topfeaturenumber = NULL,
                         subset.CpGs = NULL )
  
  
  saveRDS(SVM_model, paste0("~/Thesis_models/balanced_sets/","SVM_",set,"_mad_10k.rds"))
  
  #eXtreme Gradient Boosting (XGB)
  
  XGB_model<-maintrain( betas.. = x,
                        y.. = y,
                        method = "XGB",
                        modelcv = 5, #default
                        gridsearch = TRUE,
                        max_depth = 6, #default
                        eta = c(0.1,0.3), #default
                        gamma = c(0,0.01), #default
                        min.chwght = 1, #default
                        nrounds = 100, #default
                        colsample_bytree = c(0.2,0.5,0.8,1.0), #same as caret XGB models
                        topfeaturenumber = NULL,
                        subset.CpGs = NULL )
  
  
  saveRDS(XGB_model, paste0("~/Thesis_models/balanced_sets/","XGB_",set,"_mad_10k.rds"))
  
  mad_balanced_models[[set]]<-list(rf = RF_model, svm = SVM_model, xgb = XGB_model)
}

#--Load the models if not in the env.

algos<-c("RF","SVM","XGB")

load_balanced_models<-function(algo, model_dir){
  files<-list.files( path.expand(model_dir),
                     pattern = paste0("^", algo,"_.*_mad_10k\\.rds$"),
                     full.names = TRUE)
  method_names<-gsub(paste0("^", algo, "_|_mad_10k\\.rds$"), "", basename(files))
  setnames(lapply(files,readRDS), method_names)
}

mad_balanced_models<-setNames(lapply(algos, load_balanced_models, model_dir = "~Thesis_models/balanced_sets/"), algos)

# Testing on Test NMB

test_NMB_results<-list()

model_dict<-data.frame(model = character(), model_label = character())

for (model_key in names(mad_balanced_models)){
  
  label<-get_model_desc(model_key)
  
  message("Testing and evaluating", label, "on Test NMB")
  
  if(all(rownames(test_NMB) == NMB_test_pheno$Sample_ID)){
    test_NMB_results[[model_key]]<-test_and_evaluate_model(
      model = models[[model_key]],
      test_data = test_NMB,
      test_data_name = "Test NMB",
      true_test_labels = NMB_test_pheno$CPR_Histology,
      model_name = label[["algo"]],
      feature_selection_method = label[["feature_selection"]],
      save_in = "~/Thesis_models/test_NMB_bal/"
    )
  } else {
    stop("Mismatch in test set and test pheno Sample IDs")
  }
  
  test_NMB_results[[model_key]]$summary$model_label <- label
  
  test_NMB_results[[model_key]]$classwise$model_label<-label
  
  model_dict<-rbind(model_dict, data.frame(model=model_key, model_label = label))
}

final_test_NMB_summary<-do.call(rbind, lapply(test_NMB_results, function(x) x$summary))

final_test_NMB_classwise_summary<-do.call(rbind, lapply(test_NMB_results, function(x) x$classwise))

print(final_test_NMB_summary)

# Comparison Bar plot 

final_test_NMB_summary$algorithm <- ifelse(grepl("^RF", final_test_NMBsummary$model, ignore.case = TRUE), "RF",
                                           ifelse(grepl("^SVM", final_test_NMB_summary$model, ignore.case = TRUE), "SVM", "XGB"))

# final_test_NMB_summary$feature_selection <- ifelse(grepl("^sd", final_test_NMB_summary$model, ignore.case = TRUE),
#                                                    "SD", ifelse(grepl("^mad", final_test_NMB_summary), ignore.case = TRUE), "MAD", "All_probes")

final_test_NMB_summary$balanced <- ifelse(
  grepl("smote_tomek", final_test_NMB_summary$model_label, ignore.case = TRUE), "SMOTE-Tomek",
  ifelse(grepl("^up",    final_test_NMB_summary$model_label, ignore.case = TRUE), "Up-sampled",
         ifelse(grepl("^down",  final_test_NMB_summary$model_label, ignore.case = TRUE), "Down-sampled",
                ifelse(grepl("^smote", final_test_NMB_summary$model_label, ignore.case = TRUE), "SMOTE",
                       ifelse(grepl("tomek",  final_test_NMB_summary$model_label, ignore.case = TRUE), "Tomek",
                              "Unbalanced"))))) 

p_grouped_accuracy <- ggplot(final_test_NMB_summary, aes(x = algorithm, y = accuracy, fill = balanced)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(accuracy, 3)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 3.5) +
  labs(title = "Accuracy Comparison by Machine Learning Model trained on Top 10k probes (Median Absolute Deviation) and Balancing Method in Test NMB",
       x = "Algorithm", y = "Accuracy score", fill = "Balancing Method") +
  ylim(0, 1) +
  theme_minimal()

ggsave("~/Thesis_models/test_NMB_bal/mad_10k_accuracy_comparision.png", width = 7, height = 5, dpi = 300)


p_grouped_auc <- ggplot(final_test_NMB_summary, aes(x = algorithm, y = auc, fill = balanced)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(accuracy, 3)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 3.5) +
  labs(title = "Area Under the Curve (AUC) score Comparison by Machine Learning Model trained on Top 10k probes (Median Absolute Deviation) and Balancing Method in Test NMB",
       x = "Algorithm", y = "Accuracy score", fill = "Balancing Method") +
  ylim(0, 1) +
  theme_minimal()

ggsave("~/Thesis_models/test_NMB_bal/mad_10k_auc_comparision.png", width = 7, height = 5, dpi = 300)


# Testing on Cavalli

Cavalli_results<-list()

model_dict<-data.frame(model = character(), model_label = character())

for (model_key in names(mad_balanced_models)){
  
  label<-get_model_desc(model_key)
  
  message("Testing and evaluating", label, "on Test NMB")
  
  if(all(rownames(Cavalli_betas) == Cavalli_pheno$Study_ID...12)){
    test_NMB_results[[model_key]]<-test_and_evaluate_model(
      model = models[[model_key]],
      test_data = Cavalli_betas,
      test_data_name = "Cavalli",
      true_test_labels = Cavalli_pheno$histology,
      model_name = label[["algo"]],
      feature_selection_method = label[["feature_selection"]],
      save_in = "~/Thesis_models/test_NMB_bal/"
    )
  } else {
    stop("Mismatch in test set and test pheno Sample IDs")
  }
  
  Cavalli_results[[model_key]]$summary$model_label <- label
  
  Cavalli_results[[model_key]]$classwise$model_label<-label
  
  model_dict<-rbind(model_dict, data.frame(model=model_key, model_label = label))
}

final_Cavalli_summary<-do.call(rbind, lapply(Cavalli_results, function(x) x$summary))

final_Cavalli_classwise_summary<-do.call(rbind, lapply(Cavalli_results, function(x) x$classwise))

print(final_Cavalli_summary)

# Comparison Bar plot 

final_Cavalli_summary$algorithm <- ifelse(grepl("^RF", final_Cavallisummary$model, ignore.case = TRUE), "RF",
                                          ifelse(grepl("^SVM", final_Cavalli_summary$model, ignore.case = TRUE), "SVM", "XGB"))

# final_Cavalli_summary$feature_selection <- ifelse(grepl("^sd", final_Cavalli_summary$model, ignore.case = TRUE),
#                                                    "SD", ifelse(grepl("^mad", final_Cavalli_summary), ignore.case = TRUE), "MAD", "All_probes")

final_Cavalli_summary$balanced <- ifelse(
  grepl("smote_tomek", final_Cavalli_summary$model_label, ignore.case = TRUE), "SMOTE-Tomek",
  ifelse(grepl("^up",    final_Cavalli_summary$model_label, ignore.case = TRUE), "Up-sampled",
         ifelse(grepl("^down",  final_Cavalli_summary$model_label, ignore.case = TRUE), "Down-sampled",
                ifelse(grepl("^smote", final_Cavalli_summary$model_label, ignore.case = TRUE), "SMOTE",
                       ifelse(grepl("tomek",  final_Cavalli_summary$model_label, ignore.case = TRUE), "Tomek",
                              "Unbalanced"))))) 

p_grouped_accuracy <- ggplot(final_Cavalli_summary, aes(x = algorithm, y = accuracy, fill = balanced)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(accuracy, 3)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 3.5) +
  labs(title = "Accuracy Comparison by Machine Learning Model trained on Top 10k probes (Median Absolute Deviation) and Balancing Method in Cavalli",
       x = "Algorithm", y = "Accuracy score", fill = "Balancing Method") +
  ylim(0, 1) +
  theme_minimal()

ggsave("~/Thesis_models/Cavalli_bal/mad_10k_accuracy_comparision.png", width = 7, height = 5, dpi = 300)


p_grouped_auc <- ggplot(final_Cavalli_summary, aes(x = algorithm, y = auc, fill = balanced)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(accuracy, 3)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 3.5) +
  labs(title = "Area Under the Curve (AUC) score Comparison by Machine Learning Model trained on Top 10k probes (Median Absolute Deviation) and Balancing Method in Cavalli",
       x = "Algorithm", y = "Accuracy score", fill = "Balancing Method") +
  ylim(0, 1) +
  theme_minimal()

ggsave("~/Thesis_models/Cavalli_bal/mad_10k_auc_comparision.png", width = 7, height = 5, dpi = 300)


#----PCA plots to show default clsutering of MB's methylation groups and where synthetic samples generated by SMOTE sit ----#

require(ggplot2)

# Load datasets if not in the env.

train_NMB<-readRDS("~/DN Classifiers/train_NMB.rds")

#TODO: Save the above in the NMB folder

NMB_train_pheno<-readRDS("~/NMB/NMB_train_pheno.rds")

if (any(is.na(NMB_train_pheno$Sample_ID))){
  NMB_train_pheno<-NMB_train_pheno[!is.na(NMB_train_pheno$Sample_ID),]
}

train_NMB<-train_NMB[,!colnames(train_NMB) %in% "histo"]

NMB_train_pheno<-NMB_train_pheno[match(rownames(train_NMB), NMB_train_pheno$Sample_ID),]

bal_train_NMB<-readRDS("~/Thesis_models/balanced_sets/methylSMOTE_train_NMB.rds") 


#PCA for train NMB

NMB_pca<-prcomp(train_NMB,center = TRUE, scale. = TRUE)

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
                         "Molecular_group" = NMB_train_pheno$mol_group)
  
}

# PCA plot

NMB_PCA_plot<-ggplot(NMB_pca_df,aes(x = PC1, y = PC2,color = Molecular_group, shape = Histology)) +
  geom_point(size = 3, alpha = 0.7, stroke = 1 , color = "black") + 
  abs(title = "Patient Samples in Train_NMB in Methylation Feature Space - PCA plot", x = "PC1", y = "PC2",
      subtitle = "PC1 & PC2:Principal Components along which the samples vary the most", caption = " Clusters are mostly homogenous with MB molecular groups and histology types are not exclusive to a particular molecular group,resulting in heterogenous clusters for Histoloogy types of MB",fill = "Principle Molecular Group", shape = "Histology_Type") + 
  scale_fill_manual(values = #c("Correct_CLA" = "purple", "Correct_LCA" = "pink", 
                      #"Correct_DN" = "darkgreen", "CLA->LCA" = "#F40B7F", 
                      #"LCA->CLA" = "#F39B7F" , "CLA->DN" ="#8491B4","LCA->DN" = "#E08B35",
                      #"DN->LCA" = "#E98B88",
                      #"DN->CLA" = "#E87B90")) + 
                      c("SHH-MB" = "red", "WNT-MB" =  "blue","Group4-MB" = "darkgreen", "Group3-MB" = "yellow")) +
  
  scale_shape_manual(values = c("CLA" = 21, "DN" = 24 , "LCA" = 22,"MBEN" = 23)) +
  
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 4)),
         shape = guide_legend(override.aes = list(size = 4))) +
  theme_minimal (base_size = 12) +
  
  theme (plot.title = element_text ( face = "bold", size = 14,hjust = 0.5),
         legend.position = "bottom",
         panel.grid.minor = element_blank(),
         plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30", 
                                      margin = margin(b = 10)),
         plot.caption = element_text(size = 9, hjust = 0.5, color = "gray40",
                                     face = "italic", margin = margin(t = 10))) 


# Projecting synthetic samples generated by SMOTE onto methylation feature space with PCA fit generated above

if ( is.na(bal_train_NMB$histo)){
  message("No Histo column in the balanced set")
} else {
  bal_train_NMB$histo<-as.factor(bal_train_NMB$histo)
}

bal_train_histo<-bal_train_NMB$histo

bal_train_NMB<-bal_train_NMB[,!colnames(bal_train_NMB) %in% "histo"] 

bal_train_pca<-predict(NMB_pca,newdata = bal_train_NMB)

# To label synthetic samples apart from actual samples, we are using digest 
# to create hash representations to compare/match

require(digest)

og_hash<-hash_row(train_NMB)

bal_hash<-hash_row(bal_train_NMB)

syn_type<-ifelse(bal_hash %in% og_hash, "Given", "Synthetic")





  
