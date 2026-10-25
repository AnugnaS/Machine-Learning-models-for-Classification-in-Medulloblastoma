
# ------ Probability Calibration ---------- 

# ----- Load packages, functions and data ------ 

source("~/Custom_functions.R") 

required_pgs <- c("caret","doParallel", "probably", "tidymodels", "glmnet", "ggplot2", "matrixStats")

sapply(required_pgs, library, character.only = TRUE, quietly = TRUE)

missing_pgs <- required_pgs[!sapply(required_pgs, requireNamespace, quietly = TRUE)]

if (length(missing_pgs) > 0) install.packages(missing_pgs, dependencies = TRUE)

cal_dir<-"~/Thesis_models/Prob_cal"

dir.create(cal_dir, showWarnings = FALSE, recursive = TRUE)

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

# Using top10k by mad 

train_NMB_mat<- as.matrix(train_NMB)

if (any(is.na(train_NMB_mat))) { 
  message("NAs in train_NMB:using na.rm = TRUE")
  train_NMB_mads<-colMads(train_NMB_mat, na.rm = TRUE)
} else {
  train_NMB_mads<-colMads(train_NMB_mat, na.rm = FALSE)
} 

train_NMB_mad_10k<-train_NMB[,order(train_NMB_mads,decreasing = TRUE)[1:10000]]

train_NMB_mad_10k<-train_NMB_mad_10k[match(rownames(train_NMB),rownames(train_NMB_mad_10k)),]


if(all(rownames(train_NMB_mad_10k) == NMB_train_pheno$Sample_ID)) {
  
  train_NMB_mad_10k$histo<-as.factor(NMB_train_pheno$CPR_Histology)
  
} else {
  
  NMB_train_pheno<-NMB_train_pheno[match(rownames(train_NMB), NMB_train_pheno$Sample_ID),]
  
  train_NMB_mad_10k$histo<-as.factor(NMB_train_pheno$CPR_Histology)
  
}

# Test sets 

test_NMB<-NMB_betas[match(NMB_test_pheno$Sample_ID, rownames(NMB_betas)),]

Cavalli_betas<-readRDS("~/Cavalli/Cavalli-test/Cavalli_test_betas.rds") 

Cavalli_betas<-as.data.frame(t(Cavalli_betas)) 

Cavalli_betas<-Cavalli_betas[match(Cavalli_pheno$Sample_ID, rownames(Cavalli_betas)),] 


# ------- Probability calibration within methylClass for RF uses OOB probabilities,
# Here, we are training a caret RF model on top 10k probes (less computational cost; methylClass doesn't apply cv in RF) 
# by median absolute deviation
# with cross validation (cv) folds (k = 10) to calibrate probabilities with OOF probabilities ------------- 

#(Previous mad RF model was trained on levels CLA,LCA,MBEN whereas this iteration would only have 3 levels CLA,LCA,DN)

n_cores<-20

clusters <-makeForkCluster(n_cores)

registerDoParallel(clusters)

tryCatch({ 
  
  set.seed(42)

  cv <- trainControl(method = "cv",
                     number = 10,
                     repeats = 3,
                     classProbs = TRUE,
                     summaryFunction = multiClassSummary,
                     savePredictions = "final",
                     verboseIter = TRUE)
  
  
  RF_cal_mad10k <- train(histo ~ .,
                    data = train_NMB_mad_10k,
                    method = "rf",
                    trControl = cv,
                    metric = "logLoss")
  
  saveRDS(RF_cal_mad10k, file.path(path.expand(cal_dir),"RF_cal_mad10k.rds"))
  
  
  }, error = function(e) {
    
    message("Error pushed:", conditionMessage(e))
    
  }, finally = {
    
    stopCluster(clusters) 
    
    registerDoSEQ()
    
    message("Parallel clusters terminated")
  }) 


#----- Extract OOF probabilities and labels-----

RF_cal_mad10k<-readRDS(file.path(path.expand(cal_dir), "RF_cal_mad10k.rds"))

cv_probs <-RF_cal_mad10k$pred

cv_probs <- cv_probs[order(cv_probs$rowIndex), ]

y_true <- cv_probs$obs

X_probs <- as.matrix(cv_probs[, levels(y_true)])

# ------ Penalised multinomial logistic regression for calibration -----

set.seed(42)

calib_model <- cv.glmnet(X_probs, y_true,
                         family = "multinomial",
                         alpha = 0,  # ridge penalty
                         type.measure = "deviance",
                         nfolds = 10) 

saveRDS(calib_model, file.path(path.expand(cal_dir), "CV_ML_calibration_model.rds"))


#-------- Predict calibrated probabilities of test set -------

if (all(rownames(test_NMB) == NMB_test_pheno$Sample_ID)) {
  test_NMB_results_cal <- tryCatch(
    test_and_evaluate (
      model =RF_cal_mad10k,
      method = "caret",
      calibrated = FALSE,
      calibration_method = NULL,
      train_set_name = "Train NMB",
      test_data = test_NMB,
      test_data_name = "Test NMB",
      true_test_labels = as.factor(NMB_test_pheno$CPR_Histology),
      model_name = "Caret Random Forest",
      save_in = file.path(path.expand(cal_dir))
    ),
    error = function(e) {
      message(sprintf("FAILED [%s / %s]: %s",
                      "Random Forest trained on Train NMB", "Test NMB", conditionMessage(e)))
      NULL
    }
  )
} else {
  stop("Mismatch in test set and test pheno Sample IDs")
}

test_NMB_pred<-test_NMB_results_cal$pred_probs$pred

test_NMB_probs_uncal_mat<-as.matrix(test_NMB_results_cal$pred_probs[,!colnames(test_NMB_results_cal$pred_probs) %in% "pred"])

calib_model<-readRDS(file.path(path.expand(cal_dir), "CV_ML_calibration_model.rds"))

test_NMB_probs_cal <- predict(calib_model, 
                              newx = test_NMB_probs_uncal_mat, 
                              s = "lambda.min", 
                              type = "response") 

test_NMB_probs_cal <- test_NMB_probs_cal[, , 1] 

rownames(test_NMB_probs_cal) <- rownames(test_NMB_probs_uncal_mat)

colnames(test_NMB_probs_cal) <- colnames(test_NMB_probs_uncal_mat)

test_NMB_probs_cal_df <- as.data.frame(test_NMB_probs_cal)


NMB_predicted_histology <- colnames(test_NMB_probs_cal_df)[
  apply(test_NMB_probs_cal_df, 1, which.max)
]

test_NMB_probs_cal_df$pred<-as.factor(NMB_predicted_histology)

saveRDS(test_NMB_probs_cal_df,file.path(path.expand(cal_dir),"test_NMB_cal_probs.rds"))

test_NMB_probs_uncal_df<-as.data.frame(read.csv(file.path(path.expand(cal_dir), "Caret_Random_Forest_Train_NMB_Test_NMB_probabilities.csv")))

# ----- Reliability diagram and Brier scores to check for Calibration --------- 

stopifnot("Sample IDs differ between probs and pheno" =
            all(rownames(test_NMB_probs_cal) == NMB_test_pheno$Sample_ID))

test_NMB_probs_cal_df$true<-as.factor(NMB_test_pheno$CPR_Histology)

Test_NMB_reliab<-cal_plot_breaks(
  test_NMB_probs_cal_df,
  truth = true,
  estimate = c("CLA", "LCA", "DN"),
  num_breaks = 5
) + ggtitle("Calibrated Class Probabilities in 'Test NMB'")  + 
  theme(plot.title = element_text(size = 14,
                                  face = "bold",
                                  hjust = 0.5)) 

ggsave(filename = file.path(path.expand(cal_dir), "Test_NMB_reliab.png"),
       plot     = Test_NMB_reliab,
       width = 7, height = 5, dpi = 300)

# ------ Brier scores -----

stopifnot("Sample IDs differ between probs and pheno" =
            all(test_NMB_probs_uncal_df$Sample_ID == NMB_test_pheno$Sample_ID))

true_labels <- factor(NMB_test_pheno$CPR_Histology)

true_onehot <- model.matrix(~ true_labels - 1)

brier_uncal <- mean(rowSums((as.matrix(test_NMB_probs_uncal_df[, levels(true_labels)]) - true_onehot)^2))

brier_cal <- mean(rowSums((as.matrix(test_NMB_probs_cal_df[,!colnames(test_NMB_probs_cal_df) %in% c("pred","true")]) - true_onehot)^2))

cat("\nBrier Score (Uncalibrated) for Test NMB:", round(brier_uncal, 4), "\n")

cat("\nBrier Score (Calibrated) for Test NMB:", round(brier_cal, 4), "\n") 


# ----- Heatmap with calibrated probabilities ------

stopifnot("Mismatch in test_cal_prob_df and test pheno Sample IDs" = all(rownames(test_NMB_probs_cal_df) == NMB_test_pheno$Sample_ID))

Test_NMB_heatmap_cal <- prob_heatmaps(
  prob_df           = test_NMB_probs_cal_df,
  pred              = test_NMB_probs_cal_df$pred,
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
  heatmap_label     = "Random Forest Model for classification of histopathology,\ntrained on 'Train NMB' ,and tested on 'Test NMB' with calibrated probabilities",
  heatmap_panel     = FALSE,
  column_for_panel = "LCA",
  file_id           = "RF_cal_test_NMB",
  save_in           = file.path(path.expand(cal_dir))
)

# ------ Confusion Matrix for predictions with calibrated probabilities -------

test_NMB_cm <- tryCatch(
  caret::confusionMatrix(test_NMB_probs_cal_df$pred, NMB_test_pheno$histology, mode = "everything"),
  error = function(e) stop("confusion matrix encountered error: ", conditionMessage(e)))

test_NMB_cm_df <- as.data.frame(test_NMB_cm$table)

HISTO_COLOURS <- c("DN" = "darkgreen", "non_DN" = "purple", "non-DN" = "purple",
                   "CLA" = "purple", "LCA" = "pink", "MBEN" = "green")


test_NMB_cm_plot <- ggplot(test_NMB_cm_df, aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = Reference, alpha = Freq), color = "white") +
  geom_text(aes(label = Freq), size = 4,
            color = ifelse(test_NMB_cm_df$Freq > max(test_NMB_cm_df$Freq) * 0.5, "white", "black"), fontface = "bold") +
  scale_fill_manual(values = HISTO_COLOURS, guide = "none") +
  scale_alpha_continuous(range = c(0.15, 1), guide = "none") +
  labs(title = wrap_title(paste0("Confusion Matrix for 'Test NMB' with Calibrated Probabilities")),
       x = "True Class", y = "Predicted Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title = element_text(hjust = 0.5, face = "bold"),
        plot.margin = ggplot2::margin(t = 15, r = 15, b = 15, l = 15)) +
  coord_cartesian(clip = "off")

ggsave(file.path(path.expand(cal_dir),"test_NMB_cal_confusion_matrix.png"), 
       test_NMB_cm_plot, width = 9, height = 6, dpi = 300, bg = "white")


# ------ On Cavalli ------ 

if (all(rownames(Cavalli_betas) == Cavalli_pheno$Sample_ID)) {
  test_Cavalli_results_cal <- tryCatch(
    test_and_evaluate(
      model = RF_cal_mad10k,
      method = "caret",
      calibrated = FALSE,
      calibration_method = NULL,
      train_set_name = "Train NMB",
      test_data = Cavalli_betas,
      test_data_name = "Cavalli",
      true_test_labels = as.factor(Cavalli_pheno$histology),
      model_name = "Random Forest",
      save_in = file.path(path.expand(cal_dir)) 
    ),
    error = function(e) {
      message(sprintf("FAILED [%s / %s]: %s",
                      "Random Forest trained on Train NMB", "Cavalli", conditionMessage(e)))
      NULL
    }
  )
} else {
  stop("Mismatch in test set and test pheno Sample IDs")
} 


test_Cavalli_pred<-test_Cavalli_results_cal$pred_probs$pred

test_Cavalli_probs_uncal_mat<-as.matrix(test_Cavalli_results_cal$pred_probs[,!colnames(test_Cavalli_results_cal$pred_probs) %in% "pred"])

calib_model<-readRDS(file.path(path.expand(cal_dir), "CV_ML_calibration_model.rds"))

test_Cavalli_probs_cal <- predict(calib_model, 
                              newx = test_Cavalli_probs_uncal_mat, 
                              s = "lambda.min", 
                              type = "response") 

test_Cavalli_probs_cal <- test_Cavalli_probs_cal[, , 1] 

rownames(test_Cavalli_probs_cal) <- rownames(test_Cavalli_probs_uncal_mat)

colnames(test_Cavalli_probs_cal) <- colnames(test_Cavalli_probs_uncal_mat)

test_Cavalli_probs_cal_df <- as.data.frame(test_Cavalli_probs_cal)

Cavalli_predicted_histology <- colnames(test_Cavalli_probs_cal_df)[
  apply(test_Cavalli_probs_cal_df, 1, which.max)
]

test_Cavalli_probs_cal_df$pred<-as.factor(Cavalli_predicted_histology)

saveRDS(test_Cavalli_probs_cal_df, file.path(path.expand(cal_dir),"Cavalli_cal_probs.rds"))

test_Cavalli_probs_uncal_df<-as.data.frame(read.csv(file.path(path.expand(cal_dir), "Random_Forest_Train_NMB_Cavalli_probabilities.csv")))

# ----- Reliability diagram and Brier Scores to check for Calibration --------- 

stopifnot("Sample IDs differ between probs and pheno" =
            all(rownames(test_Cavalli_probs_cal) == Cavalli_pheno$Sample_ID))

test_Cavalli_probs_cal_df$true<-as.factor(Cavalli_pheno$histology)

Cavalli_reliab<-cal_plot_breaks(
  test_Cavalli_probs_cal_df,
  truth = true,
  estimate = c("CLA", "LCA", "DN"),
  num_breaks = 5
) + ggtitle("Calibrated Class Probabilities in 'Cavalli'") + 
  theme(plot.title = element_text(size = 14,
                                 face = "bold",
                                 hjust = 0.5))

ggsave(file.path(path.expand(cal_dir), "Cavalli_reliab.png"), plot = Cavalli_reliab, width = 7, height = 5, dpi = 300)


# ------ Brier scores -----

stopifnot("Sample IDs differ between probs and pheno" =
            all(test_Cavalli_probs_uncal_df$Sample_ID == Cavalli_pheno$Sample_ID))
  
true_labels <- factor(Cavalli_pheno$histology)

true_onehot <- model.matrix(~ true_labels - 1)

brier_uncal <- mean(rowSums((as.matrix(test_Cavalli_probs_uncal_df[, levels(true_labels)]) - true_onehot)^2))

brier_cal <- mean(rowSums((as.matrix(test_Cavalli_probs_cal_df[,!colnames(test_Cavalli_probs_cal_df) %in% c("pred","true")]) - true_onehot)^2))

cat("\nBrier Score (Uncalibrated) for Cavalli:", round(brier_uncal, 4), "\n")

cat("\nBrier Score (Calibrated) for Cavalli:", round(brier_cal, 4), "\n") 

# ----- Heatmap with calibrated probabilities ------

Cavalli_heatmap_cal <- prob_heatmaps(
  prob_df           = test_Cavalli_results_cal$pred_probs,
  pred              = test_Cavalli_results_cal$pred_probs$pred,
  test_pheno        = Cavalli_pheno,
  sample_id_column  = "Sample_ID",
  true_test_labels  = "histology",
  mol_group         = "Mol_grp",
  sub_group         = "sub_grp",
  MYC_status        = "MYC_Status",
  MYCN_status       = "MYCN_Status",
  OS                = "OS (years)",
  follow_up         = "Dead",
  Age               = "Age",
  continuum         = "LCA",
  heatmap_label     = "Random Forest Model for classification of histopathology,\ntrained on 'Train NMB',and tested on 'Cavalli' with calibrated probabilities",
  heatmap_panel     = FALSE,
  column_for_panel  = "LCA",
  file_id           = "RF_cal_Cavalli",
  save_in           = file.path(path.expand(cal_dir))
)

# ------ Confusion Matrix for predictions with calibrated probabilities -------

cavalli_cm <- tryCatch(
  caret::confusionMatrix(test_Cavalli_probs_cal_df$pred, Cavalli_pheno$histology, mode = "everything"),
  error = function(e) stop("confusion matrix encountered error: ", conditionMessage(e)))

cavalli_cm_df <- as.data.frame(cavalli_cm$table)

HISTO_COLOURS <- c("DN" = "darkgreen", "non_DN" = "purple", "non-DN" = "purple",
               "CLA" = "purple", "LCA" = "pink", "MBEN" = "green")


cavalli_cm_plot <- ggplot(cavalli_cm_df, aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = Reference, alpha = Freq), color = "white") +
  geom_text(aes(label = Freq), size = 4,
            color = ifelse(cavalli_cm_df$Freq > max(cavalli_cm_df$Freq) * 0.5, "white", "black"), fontface = "bold") +
  scale_fill_manual(values = HISTO_COLOURS, guide = "none") +
  scale_alpha_continuous(range = c(0.15, 1), guide = "none") +
  labs(title = wrap_title(paste0("Confusion Matrix for 'Cavalli' with Calibrated Probabilities")),
       x = "True Class", y = "Predicted Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title = element_text(hjust = 0.5, face = "bold"),
        plot.margin = ggplot2::margin(t = 15, r = 15, b = 15, l = 15)) + 
  coord_cartesian(clip = "off")

ggsave(file.path(path.expand(cal_dir),"Cavalli_cal_confusion_matrix.png"), cavalli_cm_plot, width = 9, height = 6, dpi = 300, bg = "white")


