#================================================================================================================================#
# ----------------------- RNAseq data Preprocessing and Differential Gene Expression------------------------------------------
#===============================================================================================================================#

# This script includes loading of samplesheets and RNA count data
# Setting colData and DESeqFromMatrix object
# Train and test split for ML models
# Normalisation with size factors and variance stablising transformation- applying train set's factors onto test set 
# Batch effect correction with limma, Combat-seq and Combat-ref
# Top variant genes -Train sets to train ML models
# Differential Gene Expression
# Feature Selection with LASSO
#==================================================================================================================#

#------------- RNAseq DESeqDataset object from RNAseq count table ----------------------------------

RNA_dir<-"~/Thesis_models/RNA_dir/"

dir.create(RNA_dir, showWarnings = FALSE, recursive = TRUE)

required_pgs <- c("dplyr", "limma","ggplot2","multcompView","ggpubr","scales")

required_bioc_pgs<-c("sva","DESeq2","edgeR","BatchQC")

missing_pgs  <- required_pgs[!sapply(required_pgs, requireNamespace, quietly = TRUE)]

missing_bioc_pgs  <- required_bioc_pgs[!sapply(required_bioc_pgs, requireNamespace, quietly = TRUE)]

if (length(missing_pgs) > 0) install.packages(missing_pgs, dependencies = TRUE)

if (length(missing_bioc_pgs) > 0) BiocManager::install(missing_bioc_pgs)

invisible(sapply(required_pgs, library, character.only = TRUE, quietly = TRUE))

invisible(sapply(required_bioc_pgs, library, character.only = TRUE, quietly = TRUE))


#------RNA-seq raw count files and samplesheet-----------

count_directory<-"/home/deant93/counts/"

RNA_counts.files<-list.files(count_directory)

RNA_seq<-read.csv("~/RNAseq/SampleSheet_NMB_PNET5Loading.csv") 

#--------- Subsetting by CPR Histology/Histology labels availability-----------

#DNA Methylation Cohort comprises all of in-house CPR samples, hence applying that to form the 'NMB RNAseq Train set'. 

NMB_Final_Cohort<-readRDS("~/NMB/NMB_Final_Cohort_26.rds")

NMB_RNA_seq<-RNA_seq[ RNA_seq$Sample_ID %in% NMB_Final_Cohort$Sample_ID2,]

NMB_Final_Cohort<-NMB_Final_Cohort[match(NMB_RNA_seq$Sample_ID, NMB_Final_Cohort$Sample_ID2),]

stopifnot(all(NMB_RNA_seq$Sample_ID == NMB_Final_Cohort$Sample_ID2))

NMB_RNA_seq$Histology<-as.factor(NMB_Final_Cohort$CPR_Histology)

NMB_RNA_seq$MYC_Amplification_Status<-as.factor(NMB_Final_Cohort$`MYC Status`)

NMB_RNA_seq$MYCN_Amplification_Status<-as.factor(NMB_Final_Cohort$`MYCN Status`)

NMB_RNA_seq$Overall.Survival.Status<-as.factor(NMB_Final_Cohort$OS)

NMB_RNA_seq$Overall.Survival.Time<-as.numeric(NMB_Final_Cohort$`Follow up years`)

NMB_RNA_seq$Age<-as.numeric(NMB_Final_Cohort$Age)

RNA_seq<-RNA_seq[!RNA_seq$Sample_ID %in% NMB_RNA_seq$Sample_ID,]

RNA_seq<-RNA_seq[!is.na(RNA_seq$Histology),]

RNA_seq<-RNA_seq[!RNA_seq$Histology == "",]

RNA_seq$Histology<-as.factor(RNA_seq$Histology)

RNA_seq_samplesheet<-rbind(NMB_RNA_seq, RNA_seq)

# Removing NOS 

RNA_seq_samplesheet <- RNA_seq_samplesheet %>% mutate(Histology = case_when( Histology %in% c("classic","Classic","CLA") ~ "CLA",
                                                                             Histology %in% c( "DN", "dnmb","MBEN") ~ "DN",
                                                                             Histology %in% c("lca","LCA") ~ "LCA"))

RNA_seq_samplesheet$Histology<-as.factor(RNA_seq_samplesheet$Histology)

RNA_seq_samplesheet<-RNA_seq_samplesheet[!is.na(RNA_seq_samplesheet$Histology),]

RNA_seq_samplesheet <- RNA_seq_samplesheet %>% mutate(Group = case_when( Group %in% c("group_3","Grp3") ~ "Group3-MB",
                                                                         Group %in% c( "group_4", "Grp4") ~ "Group4-MB",
                                                                         Group %in% c("shh","SHH") ~ "SHH-MB",
                                                                         Group == "WNT" ~ "WNT-MB",
                                                                         TRUE ~ Group
))


RNA_seq_samplesheet$Group<-as.factor(RNA_seq_samplesheet$Group)

RNA_seq_samplesheet <- RNA_seq_samplesheet %>% mutate(Subgroup = case_when( Subgroup %in% c("i","I") ~ "G34-I",
                                                                            Subgroup %in% c("ii","II") ~ "G34-II",
                                                                            Subgroup %in% c("iii","III") ~ "G34-III",
                                                                            Subgroup %in% c("iv","IV") ~ "G34-IV",
                                                                            Subgroup %in% c("v","V") ~ "G34-V",
                                                                            Subgroup %in% c("vi","VI") ~ "G34-VI",
                                                                            Subgroup %in% c("vii","VII") ~ "G34-VII",
                                                                            Subgroup %in% c("viii","VIII") ~ "G34-VIII",
                                                                            TRUE ~ Subgroup
))


RNA_seq_samplesheet$Subgroup<-as.factor(RNA_seq_samplesheet$Subgroup)


RNA_seq_samplesheet <- RNA_seq_samplesheet %>% mutate(MYC_Amplification_Status = case_when( MYC_Amplification_Status %in% c("1","AMP") ~ "1",
                                                                                            MYC_Amplification_Status %in% c("0","wt") ~ "0",
                                                                                            TRUE ~ NA
                                                                                            
))

RNA_seq_samplesheet$MYC_Amplification_Status<-as.factor(RNA_seq_samplesheet$MYC_Amplification_Status)


RNA_seq_samplesheet <- RNA_seq_samplesheet %>% mutate(MYCN_Amplification_Status = case_when( MYCN_Amplification_Status %in% c("1","AMP") ~ "1",
                                                                                             MYCN_Amplification_Status %in% c("0","wt") ~ "0",
                                                                                             TRUE ~ NA
                                                                                             
))

RNA_seq_samplesheet$MYCN_Amplification_Status<-as.factor(RNA_seq_samplesheet$MYCN_Amplification_Status)

saveRDS(RNA_seq_samplesheet,file = file.path(path.expand(RNA_dir), "RNA_seq_combined_samplesheet.rds"))

#------------Count Table and ColData for NMB--------------

sample.list<-RNA_seq_samplesheet$FileName

sample.list<-sample.list[sample.list %in% RNA_counts.files]

RNA_counts.files<-RNA_counts.files[RNA_counts.files %in% sample.list] 

sampleTable<- as.data.frame(sample.list) 

RNA_seq_samplesheet<-RNA_seq_samplesheet[match(sampleTable$sample.list, RNA_seq_samplesheet$FileName),]

sampleTable$fileName<-sampleTable$sample.list

sampleTable$Sample_ID<-as.factor(RNA_seq_samplesheet$Sample_ID)

sampleTable$SampleName <- as.factor(RNA_seq_samplesheet$Sample_ID)

sampleTable$sample.list <- NULL 

sampleTable$Cohort<-RNA_seq_samplesheet$Cohort

sampleTable$Batch <- as.factor(RNA_seq_samplesheet$Batch)

sampleTable$DetailedBatch <- as.factor(RNA_seq_samplesheet$Detailed_Batch)

sampleTable$Sex <- as.factor(RNA_seq_samplesheet$Sex)

sampleTable$Age <- as.numeric(RNA_seq_samplesheet$Age)

sampleTable$Histology <- as.factor(RNA_seq_samplesheet$Histology)

sampleTable$M_Status <- as.factor(RNA_seq_samplesheet$M_Status)

sampleTable$Group <- as.factor(RNA_seq_samplesheet$Group)

sampleTable$Subgroup <- as.factor(RNA_seq_samplesheet$Subgroup)

sampleTable$MYC <- as.factor(RNA_seq_samplesheet$MYC_Amplification_Status)

sampleTable$MYCN <- as.factor(RNA_seq_samplesheet$MYCN_Amplification_Status)

sampleTable$TP53 <-as.factor(RNA_seq_samplesheet$TP53_Status)

sampleTable$Diag_Relap<-as.factor(RNA_seq_samplesheet$DiagnosticOrRelapse)

sampleTable$OS<-as.factor(RNA_seq_samplesheet$Overall.Survival.Status)

sampleTable$Follow_up<-RNA_seq_samplesheet$Overall.Survival.Time

#sampleTable$CTNNB1 <- RNA_seq_samplesheet$CTNNB1_Status

#sampleTable$ImmuneCluster <- RNA_seq_samplesheet$Consensus_Immune_Cluster

sampleTable<-sampleTable[,!colnames(sampleTable) %in% "sample.list"]

sampleTable<-sampleTable %>% dplyr::select(SampleName, SampleName, everything()) 

#-----------DESeq object from HTseq files-------------------

RNA_seq_NMB<- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = count_directory, design = ~ DetailedBatch + Group + Histology) 

# Ignore Warning Message Note: levels of factors in the design contain characters other than
#letters, numbers, '_' and '.'. It is recommended (but not required) to use
#only letters, numbers, and delimiters '_' or '.', as these are safe characters
#for column names in R. [This is a message, not a warning or an error]- this is because Group contains levels like "Group3-MB" 

saveRDS(RNA_seq_NMB, file = file.path(path.expand(RNA_dir),"RNAseq_NMB.rds"))

#------------Count Table and ColData for other Cohorts--------------

RNA_seq_samplesheet<-readRDS(file = file.path(path.expand(RNA_dir), "RNA_seq_combined_samplesheet.rds"))

RNA_counts_sans_NMB<-read.table("~/RNAseq/all.combined.txt")

RNA_seq_samplesheet<-RNA_seq_samplesheet[RNA_seq_samplesheet$Cohort != "NMB",]

sample.list<-RNA_seq_samplesheet$Sample_ID

sample.list<-sample.list[sample.list %in% colnames(RNA_counts_sans_NMB)]

RNA_counts_sans_NMB<-RNA_counts_sans_NMB[colnames(RNA_counts_sans_NMB) %in% sample.list] 

sampleTable<- as.data.frame(sample.list) 

RNA_seq_samplesheet<-RNA_seq_samplesheet[match(sample.list, RNA_seq_samplesheet$Sample_ID),]

sampleTable$Sample_ID<-as.factor(RNA_seq_samplesheet$Sample_ID)

sampleTable$SampleName <- as.factor(RNA_seq_samplesheet$Sample_ID)

sampleTable$sample.list <- NULL 

sampleTable$Cohort<-RNA_seq_samplesheet$Cohort

sampleTable$Batch <- as.factor(RNA_seq_samplesheet$Batch)

sampleTable$DetailedBatch <- as.factor(RNA_seq_samplesheet$Detailed_Batch)

sampleTable$Sex <- as.factor(RNA_seq_samplesheet$Sex)

sampleTable$Age <- as.numeric(RNA_seq_samplesheet$Age)

sampleTable$Histology <- as.factor(RNA_seq_samplesheet$Histology)

sampleTable$M_Status <- as.factor(RNA_seq_samplesheet$M_Status)

sampleTable$Group <- as.factor(RNA_seq_samplesheet$Group)

sampleTable$Subgroup <- as.factor(RNA_seq_samplesheet$Subgroup)

sampleTable$MYC <- as.factor(RNA_seq_samplesheet$MYC_Amplification_Status)

sampleTable$MYCN <- as.factor(RNA_seq_samplesheet$MYCN_Amplification_Status)

sampleTable$TP53 <-as.factor(RNA_seq_samplesheet$TP53_Status)

sampleTable$Diag_Relap<-as.factor(RNA_seq_samplesheet$DiagnosticOrRelapse)

sampleTable$OS<-as.factor(RNA_seq_samplesheet$Overall.Survival.Status)

sampleTable$Follow_up<-RNA_seq_samplesheet$Overall.Survival.Time

#sampleTable$CTNNB1 <- RNA_seq_samplesheet$CTNNB1_Status

#sampleTable$ImmuneCluster <- RNA_seq_samplesheet$Consensus_Immune_Cluster

sampleTable<-sampleTable %>% dplyr::select(SampleName, SampleName, everything()) 


#-----------DESeq object from raw counts -------------------

RNA_seq_combined<- DESeqDataSetFromMatrix(colData = sampleTable,
                                          countData = RNA_counts_sans_NMB, design = ~ DetailedBatch + Group + Histology) 

# Ignore Warning Message Note: levels of factors in the design contain characters other than
#letters, numbers, '_' and '.'. It is recommended (but not required) to use
#only letters, numbers, and delimiters '_' or '.', as these are safe characters
#for column names in R. [This is a message, not a warning or an error]- this is because Group contains levels like "Group3-MB" 


saveRDS(RNA_seq_combined, file = file.path(path.expand(RNA_dir),"RNAseq_combined.rds"))

#------------Removing low count genes--------------------

RNA_seq_NMB<-readRDS(file.path(path.expand(RNA_dir),"RNAseq_NMB.rds"))

RNA_seq_combined<-readRDS(file.path(path.expand(RNA_dir),"RNAseq_combined.rds"))

train_RNAseq<-RNA_seq_NMB

test_RNAseq<-RNA_seq_combined

shared <- intersect(rownames(train_RNAseq), rownames(test_RNAseq))

train_RNAseq <- train_RNAseq[shared, ]

test_RNAseq  <- test_RNAseq[shared, ]

keep <- rowSums(counts(train_RNAseq)) > 0

count_data_removed<-train_RNAseq[!keep,]

saveRDS(count_data_removed, file.path(path.expand(RNA_dir), "NMB_count_data_removed.rds"))

train_RNAseq <- train_RNAseq[keep,]

test_RNAseq  <- test_RNAseq[keep,]

saveRDS(test_RNAseq,file.path(path.expand(RNA_dir),"Test_RNAseq_keep_genes.rds"))

stopifnot(all(rownames(train_RNAseq) == rownames(test_RNAseq))) 


# --------- Normalisation and Variance Stabilizing Transformation --------------------

dds_train <- estimateSizeFactors(train_RNAseq)

dds_train <- estimateDispersions(dds_train)

vst_train <- varianceStabilizingTransformation(dds_train, blind = FALSE) 

geoMeans <- exp(rowMeans(log(counts(dds_train)))) 

dds_test <- estimateSizeFactors(test_RNAseq, geoMeans = geoMeans)

# Ignore Warning Message Note: levels of factors in the design contain characters other than
#letters, numbers, '_' and '.'. It is recommended (but not required) to use
#only letters, numbers, and delimiters '_' or '.', as these are safe characters
#for column names in R. [This is a message, not a warning or an error]- this is because Group contains levels like "Group3-MB" 

dispersionFunction(dds_test) <- dispersionFunction(dds_train)  

vst_test <- varianceStabilizingTransformation(dds_test, blind = FALSE)

saveRDS(vst_test,file.path(path.expand(RNA_dir), "vst_test_sans_batch_corr.rds"))

#------------- Raw counts and processed count and metadata tables -------------------------------

train_pheno<-colData(vst_train)

test_pheno<-colData(vst_test)

train_counts<-counts(train_RNAseq)

test_counts<-counts(test_RNAseq)

train_vst_counts<-assay(vst_train)

test_vst_counts<-assay(vst_test)


# ------- PCA plots---------------------------------

GROUP_COLS <- c("WNT-MB"    = "blue",
                "SHH-MB"    = "red3",
                "Group3-MB" = "yellow3",
                "Group4-MB" = "darkgreen")

HISTO_COLS <- c("CLA"    = "purple",
                "LCA"    = "pink3",
                "DN" = "darkgreen")


pca_data <- plotPCA(vst_train,
                    intgroup  = c("Group", "Histology", "DetailedBatch"),
                    returnData = TRUE)

percentVar <- round(100 * attr(pca_data, "percentVar"))


p_histo<-ggplot(pca_data, aes(PC1, PC2, colour = Histology)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_colour_manual(values = HISTO_COLS) +
  labs(title = strwrap("PCA for NMB RNAseq count table after Normalisation and VST - Histology", width = 60),
       x = xlab(paste0("PC1: ", percentVar[1], "% variance")), 
       y = ylab(paste0("PC2: ", percentVar[2], "% variance")) , fill = "Histology") +
  theme_minimal() +
  coord_fixed() 

ggsave(filename = file.path(path.expand(RNA_dir),"PCA_Histo_sans_batch_corr.png"),plot = p_histo,width = 7, height = 5, dpi = 300)


p_mol_grp<-ggplot(pca_data, aes(PC1, PC2, colour = Group)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_colour_manual(values = GROUP_COLS) +
  labs(title = strwrap("PCA for NMB RNAseq count table after Normalisation and VST - Principal Molecular Group", width = 60),
       x = xlab(paste0("PC1: ", percentVar[1], "% variance")), 
       y = ylab(paste0("PC2: ", percentVar[2], "% variance")) , fill = "Molecular Group") +
  theme_minimal() +
  coord_fixed() 

ggsave(filename = file.path(path.expand(RNA_dir),"PCA_group_sans_batch_corr.png"),plot = p_mol_grp,width = 7, height = 5, dpi = 300)


p_batch<-ggplot(pca_data, aes(PC1, PC2, colour = Batch)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = strwrap("PCA for NMB RNAseq count table after Normalisation and VST - Batch", width = 60),
       x = xlab(paste0("PC1: ", percentVar[1], "% variance")), 
       y = ylab(paste0("PC2: ", percentVar[2], "% variance")) , fill = "Batch") +
  theme_minimal() +
  coord_fixed() 

ggsave(filename = file.path(path.expand(RNA_dir),"PCA_batch_2_sans_batch_corr.png"),plot = p_batch,width = 7, height = 5, dpi = 300)


# Check in Test set 

test_pca_data <- plotPCA(vst_test,
                    intgroup  = c("Group", "Histology", "DetailedBatch"),
                    returnData = TRUE)

percentVar_test <- round(100 * attr(test_pca_data, "percentVar"))


p_batch_test<-ggplot(test_pca_data, aes(PC1, PC2, colour = DetailedBatch)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = strwrap("PCA for Test RNAseq count table after Normalisation and VST - Batch", width = 60),
       x = xlab(paste0("PC1: ", percentVar_test[1], "% variance")), 
       y = ylab(paste0("PC2: ", percentVar_test[2], "% variance")) , fill = "Batch") +
  theme_minimal() +
  coord_fixed() 

ggsave(filename = file.path(path.expand(RNA_dir),"PCA_batch_sans_batch_corr_test.png"),plot = p_batch_test,width = 7, height = 5, dpi = 300)


# ------------ Batch effect correction ----------------

#--------With limma-------------

stopifnot(colnames(train_vst_counts) == train_pheno$Sample_ID)

design_train<-model.matrix(~Group + Histology, data = train_pheno)

batch_corr_limma<-removeBatchEffect(x = train_vst_counts,batch = train_pheno$Batch,design = design_train)

vst_train_batch_corr<-vst_train 

assay(vst_train_batch_corr)<-batch_corr_limma

saveRDS(vst_train_batch_corr,file.path(path.expand(RNA_dir), "Train_NMB_limma.rds"))

# --------- With ComBat-seq-----------------

train_cs_histo <- ComBat_seq(
  counts    = as.matrix(train_counts),
  batch     = train_pheno$Batch,
  group     = train_pheno$Histology,
  covar_mod = model.matrix(~ Group, data = train_pheno)[, -1, drop = FALSE]
)

stopifnot(all(colnames(train_cs_histo) == train_pheno$Sample_ID))


train_cs_histo_RNAseq<-DESeqDataSetFromMatrix(countData = round(train_cs_histo), 
                                                  colData = train_pheno,
                                                  design = ~ Group + Histology)


#--------With ComBat-ref----------------------------

source("~/Combat-ref/Combat_ref.R")

# functions not defined in the source file 

vec2mat <- function(vec, n_times) {
  matrix(rep(vec, n_times), ncol = n_times, byrow = FALSE)
}

vec2mat <- function(vec, n_times){
  return(matrix(rep(vec, n_times), ncol=n_times, byrow=FALSE))
}

####  Monte Carlo integration functions
monte_carlo_int_NB <- function(dat, mu, gamma, phi, gene.subset.n){
  weights <- pos_res <- list()
  for(i in 1:nrow(dat)){
    m <- mu[-i,!is.na(dat[i,])]
    x <- dat[i,!is.na(dat[i,])]
    gamma_sub <- gamma[-i]
    phi_sub <- phi[-i]
    
    # take a subset of genes to do integration - save time
    if(!is.null(gene.subset.n) & is.numeric(gene.subset.n) & length(gene.subset.n)==1){
      if(i==1){cat(sprintf("Using %s random genes for Monte Carlo integration\n", gene.subset.n))}
      mcint_ind <- sample(1:(nrow(dat)-1), gene.subset.n, replace=FALSE)
      m <- m[mcint_ind, ]; gamma_sub <- gamma_sub[mcint_ind]; phi_sub <- phi_sub[mcint_ind]
      G_sub <- gene.subset.n
    }else{
      if(i==1){cat("Using all genes for Monte Carlo integration; the function runs very slow for large number of genes\n")}
      G_sub <- nrow(dat)-1
    }
    
    #LH <- sapply(1:G_sub, function(j){sum(log2(dnbinom(x, mu=m[j,], size=1/phi_sub[j])+1))})  
    LH <- sapply(1:G_sub, function(j){prod(dnbinom(x, mu=m[j,], size=1/phi_sub[j]))})
    LH[is.nan(LH)]=0; 
    if(sum(LH)==0 | is.na(sum(LH))){
      pos_res[[i]] <- c(gamma.star=as.numeric(gamma[i]), phi.star=as.numeric(phi[i]))
    }else{
      pos_res[[i]] <- c(gamma.star=sum(gamma_sub*LH)/sum(LH), phi.star=sum(phi_sub*LH)/sum(LH))
    }
    
    weights[[i]] <- as.matrix(LH/sum(LH))
  }
  pos_res <- do.call(rbind, pos_res)
  weights <- do.call(cbind, weights)
  res <- list(gamma_star=pos_res[, "gamma.star"], phi_star=pos_res[, "phi.star"], weights=weights)	
  return(res)
} 


####  Match quantiles
# keep_zero: zero values in the original counts don't change
match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi, keep_zero=TRUE){
  new_counts_sub <- matrix(NA, nrow=nrow(counts_sub), ncol=ncol(counts_sub))
  for(a in 1:nrow(counts_sub)){
    for(b in 1:ncol(counts_sub)){
      if(keep_zero) {
        if(counts_sub[a, b] <= 1){
          new_counts_sub[a,b] <- counts_sub[a, b]
        }else{
          tmp_p <- pnbinom(counts_sub[a, b]-1, mu=old_mu[a, b], size=1/old_phi[a])
          if(abs(tmp_p-1)<1e-4){
            new_counts_sub[a,b] <- counts_sub[a, b]  
            # for outlier count, if p==1, will return Inf values -> use original count instead
          }else{
            new_counts_sub[a,b] <- 1+qnbinom(tmp_p, mu=new_mu[a, b], size=1/new_phi[a])
          }
        }
      } else {
        tmp_p = pnbinom(counts_sub[a, b], mu=old_mu[a, b], size=1/old_phi[a])
        
        if(abs(tmp_p-1) < 1e-6) {
          new_counts_sub[a,b] <- counts_sub[a, b]  
          # for outlier count, if p==1, will return Inf values -> use original count instead
        } else {
          new_counts_sub[a,b] <- qnbinom(tmp_p, mu=new_mu[a, b], size=1/new_phi[a])
        }
      }
    }
  }
  return(new_counts_sub)
}



mapDisp <- function(old_mu, new_mu, old_phi, divider){
  new_phi <- matrix(NA, nrow=nrow(old_mu), ncol=ncol(old_mu))
  for(a in 1:nrow(old_mu)){
    for(b in 1:ncol(old_mu)){
      old_var <- old_mu[a, b] + old_mu[a, b]^2 * old_phi[a, b]
      new_var <- old_var / (divider[a, b]^2)
      new_phi[a, b] <- (new_var - new_mu[a, b]) / (new_mu[a, b]^2)
    }
  }
  return(new_phi)
}     

train_cs_ref_histo <- ComBat_ref(counts = as.matrix(train_counts), 
                            batch = train_pheno$Batch, 
                            group = train_pheno$Histology)


stopifnot(all(colnames(train_cs_ref_histo) == train_pheno$Sample_ID))

train_cs_ref_histo_RNAseq<-DESeqDataSetFromMatrix(countData = round(train_cs_ref_histo), 
                                      colData = train_pheno,
                                      design = ~ Group + Histology)

saveRDS(train_cs_ref_histo_RNAseq,file.path(path.expand(RNA_dir), "Train_NMB_combat_ref_histo.rds"))



train_cs_ref_group <- ComBat_ref(counts = as.matrix(train_counts), 
                                 batch = train_pheno$Batch, 
                                 group = train_pheno$Group,
                                 genewise.disp = FALSE) 

saveRDS(train_cs_ref_group,file.path(path.expand(RNA_dir), "Train_NMB_combat_ref_group.rds"))


stopifnot(all(colnames(train_cs_ref_group) == train_pheno$Sample_ID))

train_cs_ref_group_RNAseq<-DESeqDataSetFromMatrix(countData = round(train_cs_ref_group), 
                                                  colData = train_pheno,
                                                  design = ~ Group + Histology)


corrected_sets<- list(
  Uncorrected = vst_train,
  limma       = vst_train_batch_corr,
  ComBat_ref_group  = vst(train_cs_ref_group_RNAseq, blind = FALSE),
  ComBat_seq_histo  = vst(train_cs_histo_RNAseq, blind = FALSE),
  ComBat_ref_histo  = vst(train_cs_ref_histo_RNAseq, blind = FALSE)
)

for (set in names(corrected_sets)){
  
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  
  
  pca_data <- plotPCA(corrected_sets[[set]],
                      intgroup  = c("Group", "Histology", "DetailedBatch"),
                      returnData = TRUE)
  
  
  p_histo<-ggplot(pca_data, aes(PC1, PC2, colour = Histology)) +
    geom_point(size = 3, alpha = 0.85) +
    scale_colour_manual(values = HISTO_COLS) +
    labs(title = strwrap(paste0("PCA for NMB","(", set,")", "- Histology"), width = 60),
         x = xlab(paste0("PC1: ", percentVar[1], "% variance")), 
         y = ylab(paste0("PC2: ", percentVar[2], "% variance")) , fill = "Histology") +
    theme_minimal() +
    coord_fixed() 
  
  ggsave(filename = file.path(path.expand(RNA_dir),paste0("PCA_Histo_",set,"_batch_corr.png")),plot = p_histo,width = 7, height = 5, dpi = 300)
  
  
  p_mol_grp<-ggplot(pca_data, aes(PC1, PC2, colour = Group)) +
    geom_point(size = 3, alpha = 0.85) +
    scale_colour_manual(values = GROUP_COLS) +
    labs(title = strwrap(paste0("PCA for NMB","(", set,")","-Principal Molecular Group"), width = 60),
         x = xlab(paste0("PC1: ", percentVar[1], "% variance")), 
         y = ylab(paste0("PC2: ", percentVar[2], "% variance")) , fill = "Molecular Group") +
    theme_minimal() +
    coord_fixed() 
  
  ggsave(filename = file.path(path.expand(RNA_dir),paste0("PCA_Mol_grp_",set,"_batch_corr.png")),plot = p_mol_grp,width = 7, height = 5, dpi = 300)

  p_dbatch<-ggplot(pca_data, aes(PC1, PC2, colour = DetailedBatch)) +
    geom_point(size = 3, alpha = 0.85) +
    labs(title = strwrap(paste0("PCA for NMB","(", set,")","-Detailed Batch"), width = 60),
         x = xlab(paste0("PC1: ", percentVar[1], "% variance")), 
         y = ylab(paste0("PC2: ", percentVar[2], "% variance")) , fill = "Batch") +
    theme_minimal() +
    coord_fixed() 
  
  ggsave(filename = file.path(path.expand(RNA_dir),paste0("PCA_DBatch_",set,"_batch_corr.png")),plot = p_dbatch,width = 7, height = 5, dpi = 300)
  
  
  p_batch<-ggplot(pca_data, aes(PC1, PC2, colour = Batch)) +
    geom_point(size = 3, alpha = 0.85) +
    labs(title = strwrap(paste0("PCA for NMB","(", set,")","-Batch"), width = 60),
         x = xlab(paste0("PC1: ", percentVar[1], "% variance")), 
         y = ylab(paste0("PC2: ", percentVar[2], "% variance")) , fill = "Batch") +
    theme_minimal() +
    coord_fixed() 
  
  ggsave(filename = file.path(path.expand(RNA_dir),paste0("PCA_Batch_",set,"_batch_corr.png")),plot = p_batch,width = 7, height = 5, dpi = 300)
  
}

#----Train VST (after batch correction with ComBat-ref) onto Test VST--------

RNA_dir<-"~/Thesis_models/RNA_dir/"

train_NMB_RNAseq<-readRDS(file.path(path.expand(RNA_dir), "Train_NMB_combat_ref_histo.rds"))

test_RNAseq<-readRDS(file.path(path.expand(RNA_dir),"Test_RNAseq_keep_genes.rds"))

dds_train <- estimateSizeFactors(train_NMB_RNAseq)

dds_train <- estimateDispersions(dds_train)

vst_train <- varianceStabilizingTransformation(dds_train, blind = FALSE) 

saveRDS(vst_train,file.path(path.expand(RNA_dir), "Train_RNAseq_vst.rds")) 

geoMeans <- exp(rowMeans(log(counts(dds_train)))) 

dds_test <- estimateSizeFactors(test_RNAseq, geoMeans = geoMeans)

dispersionFunction(dds_test) <- dispersionFunction(dds_train)  

vst_test <- varianceStabilizingTransformation(dds_test, blind = FALSE)

saveRDS(vst_test,file.path(path.expand(RNA_dir), "Test_RNAseq_vst.rds")) 

test_mat<-as.matrix(t(assay(vst_test)))

saveRDS(test_mat,file.path(path.expand(RNA_dir), "Test_RNAseq.rds")) 

train_NMB_mat<-as.matrix(t(assay(vst_train))) 

train_pheno<-as.data.frame(colData(vst_train))

#---------- Top variant genes-------------------------

#SDs

if (any(is.na(train_NMB_mat))) { 
  message("NAs in train_NMB:using na.rm = TRUE")
  train_NMB_sds<-colSds(train_NMB_mat, na.rm = TRUE)
} else {
  train_NMB_sds<-colSds(train_NMB_mat, na.rm = FALSE)
} 

train_NMB_sd_10k<-as.data.frame(train_NMB_mat[,order(train_NMB_sds,decreasing = TRUE)[1:10000]])

if (all(rownames(train_NMB_sd_10k) == train_pheno$Sample_ID)){
  train_NMB_sd_10k$histo<-as.factor(train_pheno$Histology)
} else { 
  train_pheno<-train_pheno[match(rownames(train_NMB_sd_10k), train_pheno$Sample_ID),]
  train_NMB_sd_10k$histo<-as.factor(train_pheno$Histology)
}

saveRDS(train_NMB_sd_10k, file.path(path.expand(RNA_dir), "Train_RNAseq_sd10k.rds"))


#MAD

if (any(is.na(train_NMB_mat))) { 
  message("NAs in train_NMB:using na.rm = TRUE") 
  train_NMB_mads<-colMads(train_NMB_mat, na.rm = TRUE)
} else {
  train_NMB_mads<-colMads(train_NMB_mat, na.rm = FALSE)
} 

train_NMB_mad_10k<-as.data.frame(train_NMB_mat[,order(train_NMB_mads,decreasing = TRUE)[1:10000]])

if (all(rownames(train_NMB_mad_10k) == train_pheno$Sample_ID)){
  train_NMB_mad_10k$histo<-as.factor(train_pheno$Histology)
} else { 
  train_pheno<-train_pheno[match(rownames(train_NMB_mat_mad_10k), train_pheno$Sample_ID),]
  train_NMB_mad_10k$histo<-as.factor(train_pheno$Histology)
} 

saveRDS(train_NMB_mad_10k, file.path(path.expand(RNA_dir), "Train_RNAseq_mad10k.rds"))

train_NMB_mat<-as.data.frame(train_NMB_mat)

if (all(rownames(train_NMB_mat) == train_pheno$Sample_ID)){
  train_NMB_mat$histo<-as.factor(train_pheno$Histology)
} else { 
  train_pheno<-train_pheno[match(rownames(train_NMB_mat), train_pheno$Sample_ID),]
  train_NMB_mat$histo<-as.factor(train_pheno$Histology)
} 

saveRDS(train_NMB_mat, file.path(path.expand(RNA_dir), "Train_RNAseq.rds"))


#----------Train VST (after batch correction with limma)-------------------------

vst_train_limma<-readRDS(file.path(path.expand(RNA_dir), "Train_NMB_limma.rds"))

train_pheno<-as.data.frame(colData(vst_train_limma))

train_NMB_mat_limma<-as.matrix(t(assay(vst_train_limma))) 

#---------- Top variant genes-------------------------

#SDs

if (any(is.na(train_NMB_mat))) { 
  message("NAs in train_NMB:using na.rm = TRUE")
  train_NMB_sds_limma<-colSds(train_NMB_mat_limma, na.rm = TRUE)
} else {
  train_NMB_sds_limma<-colSds(train_NMB_mat_limma, na.rm = FALSE)
} 

train_NMB_sd_10k_limma<-as.data.frame(train_NMB_mat_limma[,order(train_NMB_sds_limma,decreasing = TRUE)[1:10000]])

if (all(rownames(train_NMB_sd_10k_limma) == train_pheno$Sample_ID)){
  train_NMB_sd_10k_limma$histo<-as.factor(train_pheno$Histology)
} else { 
  train_pheno<-train_pheno[match(rownames(train_NMB_sd_10k_limma), train_pheno$Sample_ID),]
  train_NMB_sd_10k_limma$histo<-as.factor(train_pheno$Histology)
}

saveRDS(train_NMB_sd_10k_limma, file.path(path.expand(RNA_dir), "Train_RNAseq_sd10k_limma.rds"))


#MAD

if (any(is.na(train_NMB_mat_limma))) { 
  message("NAs in train_NMB:using na.rm = TRUE") 
  train_NMB_mads_limma<-colMads(train_NMB_mat_limma, na.rm = TRUE)
} else {
  train_NMB_mads_limma<-colMads(train_NMB_mat_limma, na.rm = FALSE)
} 

train_NMB_mad_10k_limma<-as.data.frame(train_NMB_mat_limma[,order(train_NMB_mads_limma,decreasing = TRUE)[1:10000]])

if (all(rownames(train_NMB_mad_10k_limma) == train_pheno$Sample_ID)){
  train_NMB_mad_10k_limma$histo<-as.factor(train_pheno$Histology)
} else { 
  train_pheno<-train_pheno[match(rownames(train_NMB_mad_10k_limma), train_pheno$Sample_ID),]
  train_NMB_mad_10k_limma$histo<-as.factor(train_pheno$Histology)
} 

saveRDS(train_NMB_mad_10k_limma, file.path(path.expand(RNA_dir), "Train_RNAseq_mad10k_limma.rds"))

train_NMB_mat_limma<-as.data.frame(train_NMB_mat_limma)

if (all(rownames(train_NMB_mat_limma) == train_pheno$Sample_ID)){
  train_NMB_mat_limma$histo<-as.factor(train_pheno$Histology)
} else { 
  train_pheno<-train_pheno[match(rownames(train_NMB_mat_limma), train_pheno$Sample_ID),]
  train_NMB_mat_limma$histo<-as.factor(train_pheno$Histology)
} 

saveRDS(train_NMB_mat_limma, file.path(path.expand(RNA_dir), "Train_RNAseq_limma.rds"))


# ------------------Differential Gene Expression ------------------------------------------------

#-------------Load packages and data if not in the env.-----------------------------------

require(DESeq2)

RNA_seq_NMB<-readRDS(file.path(path.expand(RNA_dir),"RNAseq_NMB.rds"))

#RNA_seq_combined<-readRDS(file.path(path.expand(RNA_dir),"RNAseq_combined.rds"))

table(RNA_seq_NMB$Group, RNA_seq_NMB$Histology)

NMB_pheno<-as.data.frame(colData(RNA_seq_NMB))

NMB_pheno<-NMB_pheno[NMB_pheno$Group != "MBNOS",]

NMB_pheno$Group<-as.factor(make.names(as.character(NMB_pheno$Group))) 

# Applying CLA histology as the baseline/reference group

NMB_pheno$Histology<-relevel(factor(NMB_pheno$Histology), ref = "CLA")

NMB_pheno$Batch<-as.factor(NMB_pheno$Batch)

NMB_pheno <- droplevels(NMB_pheno)

ids <- as.character(NMB_pheno$Sample_ID)

rownames(NMB_pheno) <- ids

NMB_counts <- as.matrix(counts(RNA_seq_NMB))[, ids]       

stopifnot(identical(colnames(NMB_counts), ids),
          identical(colnames(NMB_counts), rownames(NMB_pheno)))


design_DGE<- ~ Batch + Group + Histology  

mm<-model.matrix(design_DGE, NMB_pheno)

qr(mm)$rank == ncol(mm)

dds<-DESeqDataSetFromMatrix(countData = NMB_counts, colData = NMB_pheno , design = design_DGE)

keep <-rowSums(counts(dds) >= 10 ) >= min(table(dds$Histology)) 

dds<-dds[keep,]

nrow(dds)

dds<-DESeq(dds)

resultsNames(dds)

plotDispEsts(dds)

#--------LFC Shrinkage--------

contrasts<-list(DN_vs_CLA = c("Histology", "DN", "CLA"),
                LCA_vs_CLA = c("Histology", "LCA", "CLA"),
                LCA_vs_DN = c("Histology", "LCA", "DN"))


raw_res<-lapply(contrasts, function(c) results(dds, contrast = c , alpha = 0.05))

shrunk_res<-Map(function(c,r) 
  lfcshrink(dds,contrast = ct , res = r,type = "ashr", quite = TRUE), contrasts, raw_res)

lapply(raw_res, summary)

# To check p values

hist(raw_res$CLA_vs_LCA$p_value, breaks = 50)

sig_hits<-function(res,padj_thres = 0.05, lfc_thres = 1){
  
 sig_genes<-as.data.frame(res)
 
 sig_genes$gene<-rownames(sig_genes)
 
 sig_genes<-sig_genes[!is.na(sig_genes$p_adj) & sig_genes$padj < padj_thres & abs(sig_genes$log2FoldChange) > lfc_thres,]
 
 sig_genes[order(sig_genes$padj),]

 }

hits<-lapply(shrunk_res, sig_hits)

sapply(hits, nrow)

# ------------With interaction term in the design matrix - to see if histology effects could be modelled specific to each molecular group - Likelihood Ratio Test -------------

# Applying interactions isnt full rank as a few Molecular groups do not carry all histology types- so dropping the zero count/less count cells

int<-paste(NMB_pheno$Group, NMB_pheno$Histology , sep = "_")

int_tab<-table(int)

int <- int %>% names(int_table)[int_tab > 4]

pheno_i<-droplevels(NMB_pheno[int,])

pheno_i$int<-droplevels(as.factor(paste(pheno_i$Group, pheno_i$Histology, sep = "_")))

table(pheno_i$Group,pheno_i$Histology)

full_model <- ~Batch + int

reduced_model <- ~Batch + Group + Histology

m_full<-model.matrix(full, pheno_i)

m_red<-model.matrix(reduced, pheno_i)

stopifnot(qr(m_full)$rank == ncol(m_full), qr(m_red)$rank == ncol(m_red))

#-----Likelihood Ratio Test--------------

genes_i<-rownames(dds)

dds_i<-DESeqDataSetFromMatrix(counData = NMB_counts[genes_i, pheno_i$Sample_ID],
                              colData = pheno_i,
                              design = full_model)

dds_lrt<-DESeq(dds_i, test = "LRT", reduced = reduced_model)

res_int<-results(dds_lrt)

# Checking p-values

hist(res_int$pvalue, breaks = 50)

sum(res_int$padj < 0.1, na.rm = TRUE)

res_int<-res_int[order(res_int$padj),]

head(res_int$stat)


#-------------DGE and LRT on non-CPR cohorts---------------------

RNA_seq_combined<-readRDS(file.path(path.expand(RNA_dir),"RNAseq_combined.rds"))


table(RNA_seq_combined$Group, RNA_seq_combined$Histology)

test_pheno<-as.data.frame(colData(RNA_seq_combined))

test_pheno<-test_pheno[test_pheno$Group != "MBNOS",]

test_pheno$Group<-as.factor(make.names(as.character(test_pheno$Group))) 

# Applying CLA histology as the baseline/reference group

test_pheno$Histology<-relevel(factor(test_pheno$Histology), ref = "CLA")

test_pheno$Batch<-as.factor(test_pheno$Batch)

test_pheno <- droplevels(test_pheno)

ids <- as.character(test_pheno$Sample_ID)

rownames(test_pheno) <- ids

test_counts <- as.matrix(counts(RNA_seq_combined))[, ids]       

stopifnot(identical(colnames(test_counts), ids),
          identical(colnames(test_counts), rownames(test_pheno)))


design_DGE<- ~ Batch + Group + Histology  

mm<-model.matrix(design_DGE, test_pheno)

qr(mm)$rank == ncol(mm)

dds_non_CPR<-DESeqDataSetFromMatrix(countData = test_counts, colData = test_pheno , design = design_DGE)

keep <-rowSums(counts(dds_non_CPR) >= 10 ) >= min(table(dds_non_CPR$Histology)) 

dds_non_CPR<-dds_non_CPR[keep,]

nrow(dds_non_CPR)

dds_non_CPR<-DESeq(dds_non_CPR)

resultsNames(dds_non_CPR)

plotDispEsts(dds_non_CPR)

#--------LFC Shrinkage--------

contrasts<-list(DN_vs_CLA = c("Histology", "DN", "CLA"),
                LCA_vs_CLA = c("Histology", "LCA", "CLA"),
                LCA_vs_DN = c("Histology", "LCA", "DN"))


non_CPR_raw_res<-lapply(contrasts, function(c) results(dds_non_CPR, contrast = c , alpha = 0.05))

non_CPR_shrunk_res<-Map(function(c,r) 
  lfcshrink(dds_non_CPR,contrast = ct , res = r,type = "ashr", quite = TRUE), contrasts, raw_res)

lapply(non_CPR_shrunk_res, summary)

# To check p values

hist(non_COR_shrunk_res$CLA_vs_LCA$p_value, breaks = 50)

sig_hits<-function(res,padj_thres = 0.05, lfc_thres = 1){
  
  sig_genes<-as.data.frame(res)
  
  sig_genes$gene<-rownames(sig_genes)
  
  sig_genes<-sig_genes[!is.na(sig_genes$p_adj) & sig_genes$padj < padj_thres & abs(sig_genes$log2FoldChange) > lfc_thres,]
  
  sig_genes[order(sig_genes$padj),]
  
}

non_CPR_hits<-lapply(non_CPR_shrunk_res, sig_hits)

saveRDS(non_CPR_hits,file.path(path.expand(DGE_dir), "non_CPR_sig_genes.rds"))

sapply(non_CPR_hits, nrow)

# ------------With interaction term in the design matrix - to see if histology effects could be modelled specific to each molecular group - Likelihood Ratio Test -------------

# Applying interactions isnt full rank as a few Molecular groups do not carry all histology types- so dropping the zero count/less count cells

int<-paste(test_pheno$Group, test_pheno$Histology , sep = "_")

int_tab<-table(int)

int <- int %>% names(int_table)[int_tab > 4]

pheno_i<-droplevels(test_pheno[int,])

pheno_i$int<-droplevels(as.factor(paste(pheno_i$Group, pheno_i$Histology, sep = "_")))

table(pheno_i$Group,pheno_i$Histology)

full_model <- ~Batch + int

reduced_model <- ~Batch + Group + Histology

m_full<-model.matrix(full, pheno_i)

m_red<-model.matrix(reduced, pheno_i)

stopifnot(qr(m_full)$rank == ncol(m_full), qr(m_red)$rank == ncol(m_red))

#-----Likelihood Ratio Test--------------

genes_i<-rownames(dds)

dds_i_non_CPR<-DESeqDataSetFromMatrix(counData = test_counts[genes_i, pheno_i$Sample_ID],
                              colData = pheno_i,
                              design = full_model)

dds_lrt<-DESeq(dds_i, test = "LRT", reduced = reduced_model)

res_int<-results(dds_lrt)

# Checking p-values

hist(res_int$pvalue, breaks = 50)

sum(res_int$padj < 0.1, na.rm = TRUE)

res_int<-res_int[order(res_int$padj),]

head(res_int$stat)


dds_w <- DESeq(dds_i)

tabc  <- table(pheno_i$Group, pheno_i$Histology)
pairs <- list(c("DN", "CLA"), c("LCA", "CLA"), c("LCA", "DN"))

cts <- list()
for (g in rownames(tabc))
  for (p in pairs)
    if (all(tabc[g, p] >= min_cell))
      cts[[paste(g, p[1], "vs", p[2], sep = "_")]] <-
  c("GH", paste(g, p[1], sep = "_"), paste(g, p[2], sep = "_"))
names(cts)

res_gh <- lapply(cts, function(ct)
  lfcShrink(dds_w, contrast = ct,
            res = results(dds_w, contrast = ct, alpha = 0.05),
            type = "ashr", quiet = TRUE))
sapply(res_gh, function(r) sum(r$padj < 0.05 & abs(r$log2FoldChange) > 1, na.rm = TRUE))















