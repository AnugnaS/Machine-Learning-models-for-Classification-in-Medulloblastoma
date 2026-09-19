
********THE CODE IS IN PROGRESS, is being complied into ready-to-run workflow scripts (alongside thesis write-up)*********

Machine Learning (ML) and Deep Learning (DL) models for Histopathology and Molecular Biomarker prediction in Medulloblastoma (MB).

Main Objectives:

1. Multi-class ML models with Methylation array data for Histology Classification
2. Multi-class ML models with RNAseq data for Histology Classification.
3. Multi-class ML models with Copy-Number Variation (CNV) Data for Histology Classification.
4. Iterative integration of the above described omic data and training ML models on fused datasets for Histology Classification.
5. Binary ML models with Methylation data for Histology Classification in SHH-pathway activated MB.
6. Application of DL architectures and toolkits for preprocessing, Tissue segmentation and patching (TRIDENT toolkit, Mahmood Lab), Feature extraction (UNI v-2, Mahmood Lab) and Nuclei [Instanseg (Goldsborough et al., 2024)]  and whole cell segmentation with watershed cell expansion and with DL pipeline described by Gu, Wang, Rong, Zhao, Wu, Zhou, Wen, et al., (2025).
7. Extracted feature vector embeddings for unsupervised analysis and training ML and DL models to predict the presence of MYC amplification. 
8. Segmented cell morphometrics with HistomicsTK (Pourakpour et al., 2025), for unsupervised analysis and to classify MYC amplified and non-MYC amplified patient samples. 
9. Attention-based Multiple Instance Learning (ABMIL) with LOOCV (Leave One out Cross Validation), with the TRIDENT toolkit to extract high-attention patches. 



The R scripts include:
  
- Training an ML model with  either caret (Kuhn,2008) and/or methylClass (Liu, 2024) R package functionalities.
- Wrapper functions for the following:
     - Test and evaluate models with custom ML models, test set, labels and filenames. The evaluation methods include - Confusion matrices, accuracy scores ROC-AUC and PR-AUC curve plots and scores, F1 score and Cohen's kappa scores
     - Balancing datasets with sampling methods, Oversampling and Undersampling, and with SMOTE,ADASYN SMOTE, Tomek-link removal combined with SMOTE.
     - Annotated Heatmaps for Probability score distribution, which could also be combined with scatter plots and box plots with statistical testing for correlation/difference.
- Application of the above described wrapper functions to test the trained ML models and visualise comparison across multiple iterations with the functionalities part of the wrapper functions. 


The Jupyter notebooks include:

- Tissue segmentation and Patching
- Feature Extraction followed up with UNIv2 , PCA, UMAP, PACMAP and HDBSCAN projected onto PACMAP components.
- Nuclei segmentation, overlay of segmented nuclei masks, expansion of nuclei masks with watershed.
- Cell membrane segmentation and Whole cell segmentation with previous nuclei masks.
- Patch sampling, Stain deconvolution and morphometric measurements of masked nuclei and cell,followed up with PCA, UMAP, PACMAP and HDBSCAN projected onto PACMAP components.
- ABMIL with LOOCV and subsequent attention heatmaps.
- Training of ML and DL models respectively. 
