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

required_pgs <- c("dplyr", "limma","ggplot2","multcompView","ggpubr","scales", "ggplot2","tidyr","ashr")

required_bioc_pgs<-c("sva","DESeq2","edgeR","BatchQC", "biomaRt","org.Hs.eg.db")

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

NMB_RNA_seq$Sex<-as.factor(NMB_Final_Cohort$Sex)

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

RNA_seq_samplesheet <- RNA_seq_samplesheet %>% mutate(Sex = case_when( Sex %in% c("m","M","m_") ~ "M",
                                                                                            Sex %in% c("f","F") ~ "F",
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

sampleTable$Sex<-as.factor(RNA_seq_samplesheet$Sex)

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

#------------Removing low count and sex chr genes--------------------

RNA_seq_NMB<-readRDS(file.path(path.expand(RNA_dir),"RNAseq_NMB.rds"))

RNA_seq_combined<-readRDS(file.path(path.expand(RNA_dir),"RNAseq_combined.rds"))

train_RNAseq<-RNA_seq_NMB

test_RNAseq<-RNA_seq_combined

keep <- rowSums(counts(train_RNAseq)) > 0

count_data_removed <- train_RNAseq[!keep, ]

train_RNAseq       <- train_RNAseq[keep, ]

shared <- intersect(rownames(train_RNAseq), rownames(test_RNAseq))

train_RNAseq <- train_RNAseq[shared, ]

test_RNAseq  <- test_RNAseq[shared, ]

stopifnot(identical(rownames(train_RNAseq), rownames(test_RNAseq)))

length(shared)   

# Removing genes on sex chromosomes 

require(org.Hs.eg.db)

gene_symbols <- rownames(train_RNAseq) 

chr <- mapIds(org.Hs.eg.db,
              keys      = gene_symbols,
              column    = "CHR",
              keytype   = "SYMBOL",
              multiVals = "first")

table(chr, useNA = "ifany")   

autosomal <- !is.na(chr) & chr %in% as.character(1:22)

sex_genes <- !is.na(chr) & chr %in% c("X", "Y")

sex_count_data<-train_RNAseq[sex_genes,]

train_RNAseq <- train_RNAseq[autosomal, ]

saveRDS(train_RNAseq,file.path(path.expand(RNA_dir),"Train_RNAseq_keep_genes.rds"))

count_data_removed<-rbind(count_data_removed,sex_count_data)

saveRDS(count_data_removed, file.path(path.expand(RNA_dir), "NMB_count_data_removed.rds"))

test_RNAseq <- test_RNAseq[autosomal, ]

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

# Batch-correcting non-CPR cohorts test set with ComBat-ref

test_RNAseq<-readRDS(file.path(path.expand(RNA_dir),"Test_RNAseq_keep_genes.rds"))

test_pheno<-as.data.frame(colData(test_RNAseq))

test_counts<-counts(test_RNAseq)

test_cs_ref_histo <- ComBat_ref(counts = as.matrix(test_counts), 
                                 batch = test_pheno$Batch, 
                                 group = test_pheno$Histology)

saveRDS(test_cs_ref_histo, file.path(path.expand(RNA_dir), "Test_RNAseq_keep_genes_bc.rds"))

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

#-------------Load data if not in the env.-----------------------------------

RNA_dir<-"~/Thesis_models/RNA_dir/"

DGE_dir<-"~/Thesis_models/RNA_dir/DGE/"

dir.create(DGE_dir, showWarnings = FALSE, recursive = TRUE)

train_RNAseq<-readRDS(file.path(path.expand(RNA_dir),"Train_RNAseq_keep_genes.rds"))

test_RNAseq<-readRDS(file.path(path.expand(RNA_dir),"Test_RNAseq_keep_genes.rds"))

test_RNAseq_bc<-readRDS(file.path(path.expand(RNA_dir),"Test_RNAseq_keep_genes_bc.rds"))


contrast_list <- list(
  DN_vs_CLA  = c("Histology", "DN",  "CLA"),
  LCA_vs_CLA = c("Histology", "LCA", "CLA"),
  LCA_vs_DN  = c("Histology", "LCA", "DN")
)

status_cols <- c(Up = "firebrick", Down = "steelblue", NS = "grey75")

GROUP_COLS <- c("WNT-MB"    = "blue",
                "SHH-MB"    = "red3",
                "Group3-MB" = "yellow3",
                "Group4-MB" = "darkgreen")

HISTO_COLS <- c("CLA"    = "purple",
                "LCA"    = "pink3",
                "DN" = "darkgreen")

#---------Helper Functions for DGE--------------------------------

res_table <- function(raw, shrunk) {
  df <- as.data.frame(raw)[, c("baseMean", "log2FoldChange", "pvalue", "padj")]
  
  names(df)[names(df) == "log2FoldChange"] <- "log2FC_raw"
  
  df$log2FC <- shrunk$log2FoldChange
  
  df$gene   <- rownames(df) 
  
  df
}

DGE_change <- function(df, lfc_threshold, padj_cutoff) {
  
  df$status <- "NS"
  
  ok <- !is.na(df$padj) & !is.na(df$log2FC) & df$padj < padj_cutoff
  
  df$status[ok & df$log2FC >  lfc_threshold] <- "Up"
  
  df$status[ok & df$log2FC < -lfc_threshold] <- "Down"
  
  df$status <- factor(df$status, levels = c("Up", "Down", "NS"))
  
  df
}

sig_hits <- function(df,padj_cutoff, lfc_threshold) {
  
  df <- DGE_change(df, padj_cutoff = padj_cutoff, lfc_threshold = lfc_threshold)
  
  df <- df[df$status != "NS", ]
  
  df[order(df$padj), ]
}


volcano_plot <- function(df, title, lfc_threshold, padj_cutoff, n_label = 10) {
  
  df <- DGE_change(df,padj_cutoff = padj_cutoff,lfc_threshold = lfc_threshold)
  
  df <- df[!is.na(df$padj) & !is.na(df$log2FC), ]
  
  df$nlp <- -log10(pmax(df$padj, 1e-300))

  lab <- head(sig_hits(df, padj_cutoff = padj_cutoff, lfc_threshold = lfc_threshold), n_label)
  
  lab$nlp <- -log10(pmax(lab$padj, 1e-300))
  
  ggplot(df, aes(log2FC, nlp, colour = status)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_colour_manual(values = status_cols) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    geom_text_repel(data = lab, aes(label = gene), colour = "black",
                    size = 3, max.overlaps = Inf) +
    labs(title = gsub("_", " ", title),
         subtitle = sprintf("Up: %d | Down: %d  (padj < %g, |shrunk log2FC| > %g)",
                            sum(df$status == "Up"), sum(df$status == "Down"),
                            padj_cutoff, lfc_threshold),
         x = "Shrunk log2 fold change",
         y = "-log10 adjusted p-value", colour = NULL) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(hjust = 0.5, face = "bold"))
}

ma_plot <- function(dfs,lfc_threshold,padj_cutoff) {
  
  df <- do.call(rbind, lapply(names(dfs), function(n)
    data.frame(DGE_change(dfs[[n]],lfc_threshold = lfc_threshold,padj_cutoff = padj_cutoff), contrast = n)))
  
  df <- df[!is.na(df$log2FC), ]
  
  ggplot(df, aes(log10(baseMean + 1), log2FC, colour = status)) +
    geom_point(alpha = 0.4, size = 0.7) +
    scale_colour_manual(values = status_cols) +
    geom_hline(yintercept = 0) +
    facet_wrap(~contrast) +
    labs(x = "log10(mean normalised count + 1)", y = "Shrunk log2 fold change",
         colour = NULL) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(hjust = 0.5, face = "bold"))
}

pval_hist <- function(dfs) {
  
  df <- do.call(rbind, lapply(names(dfs), function(n)
    data.frame(pvalue = dfs[[n]]$pvalue, contrast = n)))
  
  ggplot(df[!is.na(df$pvalue), ], aes(pvalue)) +
    geom_histogram(bins = 50, fill = "grey60", colour = "white") +
    facet_wrap(~contrast) +
    labs(title = "Raw p-value histograms",
         subtitle = "Healthy: flat bulk plus a spike near 0",
         x = "p-value", y = "Genes") +
    theme_bw()
}

count_plot <- function(dfs, title,padj_cutoff, lfc_threshold) {
  
  df <- do.call(rbind, lapply(names(dfs), function(n) {
    
    t <- DGE_change(dfs[[n]],padj_cutoff = padj_cutoff,lfc_threshold = lfc_threshold)
    data.frame(contrast = n, direction = c("Up", "Down"),
               n = c(sum(t$status == "Up"), -sum(t$status == "Down")))
  }))
  
  ggplot(df, aes(contrast, n, fill = direction)) +
    geom_col() + geom_hline(yintercept = 0) +
    scale_fill_manual(values = c(Up = "firebrick", Down = "steelblue")) +
    coord_flip() +
    labs(title = title, x = NULL, fill = NULL,
         y = "Significant genes") +
    theme_bw()
}

pca_plot <- function(mat, pheno, colour, shape = NULL, title = "") {
  
  top <- order(apply(mat, 1, var), decreasing = TRUE)[seq_len(min(1000, nrow(mat)))]
  
  pc  <- prcomp(t(mat[top, ]))
  
  pv  <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  
  df  <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                    pheno[rownames(pc$x), c("Histology", "Group", "Batch")])
  
  p <- ggplot(df, aes(PC1, PC2, colour = .data[[colour]])) +
    geom_point(size = 2.5, alpha = 0.85) +
    labs(title = title, x = paste0("PC1 (", pv[1], "%)"),
         y = paste0("PC2 (", pv[2], "%)")) +
    theme_bw()
  if (!is.null(shape)) p <- p + aes(shape = .data[[shape]])
  
  p
}

box_genes <- function(mat, pheno, genes, facet_by_group = FALSE) {
  
  pheno <- droplevels(pheno)
  
  df <- do.call(rbind, lapply(genes, function(g)
    data.frame(gene = g, expr = as.numeric(mat[g, rownames(pheno)]),
               pheno[, c("Histology", "Group", "Batch")], row.names = NULL)))
  
  df$gene <- factor(df$gene, levels = genes)
  
  if (!facet_by_group) {
    
    n <- table(pheno$Histology)
    
    df$Histology <- factor(df$Histology, levels = names(n),
                           labels = paste0(names(n), "\n(n=", n, ")"))
  }
  
  p <- ggplot(df, aes(Histology, expr)) +
    geom_boxplot(aes(fill = Histology), outlier.shape = NA, alpha = 0.4,
                 show.legend = FALSE) +
    geom_jitter(aes(colour = Group), width = 0.15, size = 1.6, alpha = 0.8) +
    labs(y = "VST expression", x = NULL) +
    theme_bw()
  if (facet_by_group) p + facet_grid(gene ~ Group, scales = "free_y")
  else                p + facet_wrap(~gene, scales = "free_y")
}

heatmap_top <- function(mat, pheno, genes, file) {
  
  m <- t(scale(t(mat[genes, rownames(pheno)])))
  
  m[m >  3] <-  3
  
  m[m < -3] <- -3
  
  pheatmap(m, annotation_col = pheno[, c("Histology", "Group", "Batch")],
           show_colnames = FALSE, show_rownames = length(genes) <= 60,
           clustering_method = "ward.D2", filename = file, width = 9, height = 10)
}

concordance_plot <- function(a, b, nm, lab_a = "NMB", lab_b = "non-CPR", padj_cutoff) {
  
  m <- merge(a[, c("gene", "log2FC", "padj")], b[, c("gene", "log2FC", "padj")],
             by = "gene", suffixes = c(".a", ".b"))
  
  m <- m[complete.cases(m), ]
  
  rho <- cor(m$log2FC.a, m$log2FC.b, method = "spearman")
  
  m$sig <- ifelse(m$padj.a < padj_cutoff & m$padj.b < padj_cutoff, "Both",
                  ifelse(m$padj.a < padj_cutoff, paste(lab_a, "only"),
                         ifelse(m$padj.b < padj_cutoff, paste(lab_b, "only"), "Neither")))
  
  ggplot(m, aes(log2FC.a, log2FC.b, colour = sig)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, colour = "grey40") +
    labs(title = paste(gsub("_", " ", nm), "- agreement between cohorts"),
         subtitle = sprintf("Spearman rho = %.2f (%d shared genes)", rho, nrow(m)),
         x = paste(lab_a, "shrunk log2FC"), y = paste(lab_b, "shrunk log2FC"),
         colour = NULL) +
    theme_bw()
}

# ---------------------------- Custom Function for DGE --------------------------

RNA_dir<-"~/Thesis_models/RNA_dir/"

DGE_dir<-"~/Thesis_models/RNA_dir/DGE/"

dir.create(DGE_dir, showWarnings = FALSE, recursive = TRUE)

RNAseq_NMB<-readRDS(file.path(path.expand(RNA_dir), "RNAseq_NMB.rds"))

RNAseq_combined<-readRDS(file.path(path.expand(RNA_dir), "RNAseq_combined.rds"))

contrast_list <- list(
  DN_vs_CLA  = c("Histology", "DN",  "CLA"),
  LCA_vs_CLA = c("Histology", "LCA", "CLA"),
  LCA_vs_DN  = c("Histology", "LCA", "DN")
)

status_cols <- c(Up = "firebrick", Down = "steelblue", NS = "grey75")

GROUP_COLS <- c("WNT-MB"    = "blue",
                "SHH-MB"    = "red3",
                "Group3-MB" = "yellow3",
                "Group4-MB" = "darkgreen")

HISTO_COLS <- c("CLA"    = "purple",
                "LCA"    = "pink3",
                "DN" = "darkgreen")

#---------Helper Functions for DGE--------------------------------

check_cols <- function(cols, lv, what) {
  if (is.null(cols)) return(NULL)
  if (!all(lv %in% names(cols))) {
    warning(what, ": no colour defined for level(s): ",
            paste(setdiff(lv, names(cols)), collapse = ", "),
            " -> using default colours. Names must match levels(pheno$", what, ").")
    return(NULL)
  }
  cols
}

res_table <- function(raw, shrunk) {
  df <- as.data.frame(raw)[, c("baseMean", "log2FoldChange", "pvalue", "padj")]
  
  names(df)[names(df) == "log2FoldChange"] <- "log2FC_raw"
  
  df$log2FC <- shrunk$log2FoldChange
  
  df$gene   <- rownames(df) 
  
  df
}

DGE_change <- function(df, lfc_threshold, padj_cutoff) {
  
  df$status <- "NS"
  
  ok <- !is.na(df$padj) & !is.na(df$log2FC) & df$padj < padj_cutoff
  
  df$status[ok & df$log2FC >  lfc_threshold] <- "Up"
  
  df$status[ok & df$log2FC < -lfc_threshold] <- "Down"
  
  df$status <- factor(df$status, levels = c("Up", "Down", "NS"))
  
  df
}

sig_hits <- function(df,padj_cutoff, lfc_threshold) {
  
  df <- DGE_change(df, padj_cutoff = padj_cutoff, lfc_threshold = lfc_threshold)
  
  df <- df[df$status != "NS", ]
  
  df[order(df$padj), ]
}


volcano_plot <- function(df, title, lfc_threshold, padj_cutoff, n_label = 10) {
  
  df <- DGE_change(df,padj_cutoff = padj_cutoff,lfc_threshold = lfc_threshold)
  
  df <- df[!is.na(df$padj) & !is.na(df$log2FC), ]
  
  df$nlp <- -log10(pmax(df$padj, 1e-300))
  
  lab <- head(sig_hits(df, padj_cutoff = padj_cutoff, lfc_threshold = lfc_threshold), n_label)
  
  lab$nlp <- -log10(pmax(lab$padj, 1e-300))
  
  ggplot(df, aes(log2FC, nlp, colour = status)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_colour_manual(values = status_cols) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    geom_text_repel(data = lab, aes(label = gene), colour = "black",
                    size = 3, max.overlaps = Inf) +
    labs(title = gsub("_", " ", title),
         subtitle = sprintf("Up: %d | Down: %d  (padj < %g, |shrunk log2FC| > %g)",
                            sum(df$status == "Up"), sum(df$status == "Down"),
                            padj_cutoff, lfc_threshold),
         x = "Shrunk log2 fold change",
         y = "-log10 adjusted p-value", colour = NULL) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(hjust = 0.5, face = "bold"))
}

ma_plot <- function(dfs,lfc_threshold,padj_cutoff) {
  
  df <- do.call(rbind, lapply(names(dfs), function(n)
    data.frame(DGE_change(dfs[[n]],lfc_threshold = lfc_threshold,padj_cutoff = padj_cutoff), contrast = n)))
  
  df <- df[!is.na(df$log2FC), ]
  
  ggplot(df, aes(log10(baseMean + 1), log2FC, colour = status)) +
    geom_point(alpha = 0.4, size = 0.7) +
    scale_colour_manual(values = status_cols) +
    geom_hline(yintercept = 0) +
    facet_wrap(~contrast) +
    labs(x = "log10(mean normalised count + 1)", y = "Shrunk log2 fold change",
         colour = NULL) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(hjust = 0.5, face = "bold"))
}

pval_hist <- function(dfs) {
  
  df <- do.call(rbind, lapply(names(dfs), function(n)
    data.frame(pvalue = dfs[[n]]$pvalue, contrast = n)))
  
  ggplot(df[!is.na(df$pvalue), ], aes(pvalue)) +
    geom_histogram(bins = 50, fill = "grey60", colour = "white") +
    facet_wrap(~contrast) +
    labs(title = "Raw p-value histograms",
         subtitle = "Healthy: flat bulk plus a spike near 0",
         x = "p-value", y = "Genes") +
    theme_bw()
}

count_plot <- function(dfs, title,padj_cutoff, lfc_threshold) {
  
  df <- do.call(rbind, lapply(names(dfs), function(n) {
    
    t <- DGE_change(dfs[[n]],padj_cutoff = padj_cutoff,lfc_threshold = lfc_threshold)
    data.frame(contrast = n, direction = c("Up", "Down"),
               n = c(sum(t$status == "Up"), -sum(t$status == "Down")))
  }))
  
  ggplot(df, aes(contrast, n, fill = direction)) +
    geom_col() + geom_hline(yintercept = 0) +
    scale_fill_manual(values = c(Up = "firebrick", Down = "steelblue")) +
    coord_flip() +
    labs(title = title, x = NULL, fill = NULL,
         y = "Significant genes") +
    theme_bw()
}

pca_plot <- function(mat, pheno, colour, shape = NULL, title = "", cols = NULL) {
  
  top <- order(apply(mat, 1, var), decreasing = TRUE)[seq_len(min(1000, nrow(mat)))]
  
  pc  <- prcomp(t(mat[top, ]))
  
  pv  <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  
  df  <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                    pheno[rownames(pc$x), c("Histology", "Group", "Batch")])
  
  p <- ggplot(df, aes(PC1, PC2, colour = .data[[colour]])) +
    geom_point(size = 2.5, alpha = 0.85) +
    labs(title = title, x = paste0("PC1 (", pv[1], "%)"),
         y = paste0("PC2 (", pv[2], "%)")) +
    theme_bw()
  if (!is.null(shape)) p <- p + aes(shape = .data[[shape]])
  
  if (!is.null(cols)) p <- p + scale_colour_manual(values = cols)   # CHANGED
  
  p
}

box_genes <- function(mat, pheno, genes, facet_by_group = FALSE,
                      histo_cols = NULL, group_cols = NULL) {
  
  pheno <- droplevels(pheno)
  
  df <- do.call(rbind, lapply(genes, function(g)
    data.frame(gene = g, expr = as.numeric(mat[g, rownames(pheno)]),
               pheno[, c("Histology", "Group", "Batch")], row.names = NULL)))
  
  df$gene <- factor(df$gene, levels = genes)
  
  if (!facet_by_group) {
    
    n <- table(pheno$Histology)
    
    df$Histology <- factor(df$Histology, levels = names(n),
                           labels = paste0(names(n), "\n(n=", n, ")"))
  }
  
  p <- ggplot(df, aes(Histology, expr)) +
    geom_boxplot(aes(fill = Histology), outlier.shape = NA, alpha = 0.4,
                 show.legend = FALSE) +
    geom_jitter(aes(colour = Group), width = 0.15, size = 1.6, alpha = 0.8) +
    labs(y = "VST expression", x = NULL) +
    theme_bw()
  
  if (!is.null(histo_cols) && all(levels(pheno$Histology) %in% names(histo_cols))) {
    fill_vals <- histo_cols[levels(pheno$Histology)]
    names(fill_vals) <- levels(df$Histology)
    p <- p + scale_fill_manual(values = fill_vals)
  }
  
  if (!is.null(group_cols) && all(levels(pheno$Group) %in% names(group_cols)))
    p <- p + scale_colour_manual(values = group_cols)
  
  if (facet_by_group) p + facet_grid(gene ~ Group, scales = "free_y")
  else                p + facet_wrap(~gene, scales = "free_y")
}

heatmap_top <- function(mat, pheno, genes, file, histo_cols = NULL, group_cols = NULL) {
  
  m <- t(scale(t(mat[genes, rownames(pheno)])))
  
  m[m >  3] <-  3
  
  m[m < -3] <- -3
  
  ann_cols <- Filter(Negate(is.null), list(Histology = histo_cols, Group = group_cols))  
  
  pheatmap(m, annotation_col = pheno[, c("Histology", "Group", "Batch")],
           annotation_colors = if (length(ann_cols) > 0) ann_cols else NA,                
           show_colnames = FALSE, show_rownames = length(genes) <= 60,
           clustering_method = "ward.D2", filename = file, width = 9, height = 10)
}

concordance_plot <- function(a, b, nm, lab_a = "NMB", lab_b = "non-CPR", padj_cutoff) {
  
  m <- merge(a[, c("gene", "log2FC", "padj")], b[, c("gene", "log2FC", "padj")],
             by = "gene", suffixes = c(".a", ".b"))
  
  m <- m[complete.cases(m), ]
  
  rho <- cor(m$log2FC.a, m$log2FC.b, method = "spearman")
  
  m$sig <- ifelse(m$padj.a < padj_cutoff & m$padj.b < padj_cutoff, "Both",
                  ifelse(m$padj.a < padj_cutoff, paste(lab_a, "only"),
                         ifelse(m$padj.b < padj_cutoff, paste(lab_b, "only"), "Neither")))
  
  ggplot(m, aes(log2FC.a, log2FC.b, colour = sig)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, colour = "grey40") +
    labs(title = paste(gsub("_", " ", nm), "- agreement between cohorts"),
         subtitle = sprintf("Spearman rho = %.2f (%d shared genes)", rho, nrow(m)),
         x = paste(lab_a, "shrunk log2FC"), y = paste(lab_b, "shrunk log2FC"),
         colour = NULL) +
    theme_bw()
}

# ---------------------------- Custom Function for DGE --------------------------

run_dge <- function(se, label, lfc_threshold, padj_cutoff, include_batch = TRUE) {
  
  out_dir <- file.path(DGE_dir, label) 
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  required_pgs <- c("dplyr", "limma","ggplot2","multcompView","ggpubr","scales", "ggplot2","tidyr", "ashr", "ggrepel","pheatmap")
  
  required_bioc_pgs<-c("sva","DESeq2","edgeR","BatchQC", "biomaRt","org.Hs.eg.db")
  
  missing_pgs  <- required_pgs[!sapply(required_pgs, requireNamespace, quietly = TRUE)]
  
  missing_bioc_pgs  <- required_bioc_pgs[!sapply(required_bioc_pgs, requireNamespace, quietly = TRUE)]
  
  if (length(missing_pgs) > 0) install.packages(missing_pgs, dependencies = TRUE)
  
  if (length(missing_bioc_pgs) > 0) BiocManager::install(missing_bioc_pgs)
  
  invisible(sapply(required_pgs, library, character.only = TRUE, quietly = TRUE))
  
  invisible(sapply(required_bioc_pgs, library, character.only = TRUE, quietly = TRUE))
  
  
  pheno <- as.data.frame(colData(se))
  
  pheno <- pheno[!is.na(pheno$Group) & !is.na(pheno$Histology) &
                   pheno$Group != "MBNOS", ]
  
  pheno$Group     <- factor(make.names(as.character(pheno$Group)))
  
  pheno$Histology <- relevel(factor(pheno$Histology), ref = "CLA")  # CLA = baseline
  
  pheno$Batch     <- factor(pheno$Batch)
  
  pheno <- droplevels(pheno)
  
  ids <- as.character(pheno$Sample_ID)
  
  rownames(pheno) <- ids
  
  print(table(pheno$Group, pheno$Histology))
  
  cnt <- as.matrix(counts(se))[, ids]
  
  if (include_batch) {design_DGE <- ~ Batch + Group + Histology } else { 
    
    design_DGE <- ~ Group + Histology }
  
  mm <- model.matrix(design_DGE, pheno)
  
  stopifnot(qr(mm)$rank == ncol(mm))          # design should be full rank
  
  dds  <- DESeqDataSetFromMatrix(countData = cnt, colData = pheno, design = design_DGE)
  
  dds  <- DESeq(dds)
  
  pdf(file.path(out_dir, "dispersion.pdf")); plotDispEsts(dds); invisible(dev.off())
  
  raw <- lapply(contrast_list, function(ct) results(dds, contrast = ct, alpha = padj_cutoff))
  
  shrunk <- Map(function(ct, r)
    lfcShrink(dds, contrast = ct, res = r, type = "ashr", quiet = TRUE),
    contrast_list, raw)
  
  for (n in names(raw)) { cat("\n", label, n, "\n"); summary(raw[[n]]) }
  
  tabs <- Map(res_table, raw, shrunk)
  
  hits <- lapply(tabs, function(tab) sig_hits(tab, padj_cutoff = padj_cutoff, lfc_threshold = lfc_threshold))
  
  print(sapply(hits, nrow))
  
  vsd <- vst(dds, blind = TRUE)
  
  mat <- assay(vsd)
  
  mat_adj <- if (nlevels(pheno$Batch) > 1) {
    removeBatchEffect(mat, batch = pheno$Batch,
                      design = model.matrix(~ Group + Histology, pheno))
  } else mat
  
  GROUP_COLS <- c("WNT-MB"    = "blue",
                  "SHH-MB"    = "red3",
                  "Group3-MB" = "yellow3",
                  "Group4-MB" = "darkgreen")
  
  names(GROUP_COLS) <- make.names(names(GROUP_COLS))   
  
  HISTO_COLS <- c(CLA = "purple", LCA = "pink3", DN = "darkgreen")
  
  GROUP_COLS <- check_cols(GROUP_COLS, levels(pheno$Group),     "Group")
  HISTO_COLS <- check_cols(HISTO_COLS, levels(pheno$Histology), "Histology")
  
  for (v in c("Histology", "Group", "Batch")) {
    
    if (v == "Group") {
      COL <- GROUP_COLS
    } else if (v == "Histology") {
      COL <- HISTO_COLS
    } else {
      COL <- NULL      
    }
        ggsave(filename = file.path(out_dir, paste0("pca_by_", v, ".pdf")),
           plot = pca_plot(mat, pheno, v, shape = if (v == "Batch") NULL else "Batch",
                           title = paste(label, "PCA (VST, before batch adjustment)"),
                           cols = COL),
           width = 6.5, height = 5)
    
  }   
  
  ggsave(filename = file.path(out_dir, "pca_batch_adjusted.pdf"),
         plot = pca_plot(mat_adj, pheno, "Histology", shape = "Group",
                         title = paste(label, "PCA (batch-adjusted)"),
                         cols = HISTO_COLS),                                   
         width = 6.5, height = 5)
  
  ggsave(filename = file.path(out_dir, "pvalue_histograms.pdf"),
         plot = pval_hist(tabs), width = 9, height = 4)
  
  ggsave(filename = file.path(out_dir, "ma_plots.pdf"),
         plot = ma_plot(tabs, lfc_threshold = lfc_threshold, padj_cutoff = padj_cutoff),
         width = 10, height = 4)
  
  ggsave(filename = file.path(out_dir, "significant_gene_counts.pdf"),
         plot = count_plot(tabs, paste(label, "- significant genes per contrast"),
                           padj_cutoff = padj_cutoff, lfc_threshold = lfc_threshold),
         width = 6.5, height = 3.5)
  
  
  for (nm in names(tabs)) {
    ggsave(filename = file.path(out_dir, paste0("volcano_", nm, ".pdf")),
           plot = volcano_plot(tabs[[nm]], nm, lfc_threshold = lfc_threshold, padj_cutoff = padj_cutoff),
           width = 6.5, height = 5.5)
    
    top <- head(hits[[nm]]$gene, 6)
    
    if (length(top) > 0) {
      ph <- pheno[pheno$Histology %in% contrast_list[[nm]][2:3], ]
      ggsave(filename = file.path(out_dir, paste0("boxplot_top_", nm, ".pdf")),
             plot = box_genes(mat_adj, ph, top,
                              histo_cols = HISTO_COLS, group_cols = GROUP_COLS),   
             width = 9, height = 6)
    }
  }
  
  hm_genes <- unique(unlist(lapply(hits, function(h) head(h$gene, 25))))
  
  if (length(hm_genes) > 1)
    heatmap_top(mat_adj, pheno, hm_genes, file.path(out_dir, "heatmap_top_hits.pdf"),
                histo_cols = HISTO_COLS, group_cols = GROUP_COLS)                  
  
  ## ------------With interaction term in the design matrix - to see if histology effects could be modelled specific to each molecular group - Likelihood Ratio Test -------------
  
  # Applying interactions isnt full rank as a few Molecular groups do not carry all histology types- so dropping the zero count/less count cells
  
  cell     <- paste(pheno$Group, pheno$Histology, sep = "_")
  
  cell_tab <- table(cell)
  
  pheno_i  <- droplevels(pheno[cell %in% names(cell_tab)[cell_tab >= 5], ])
  
  pheno_i$GH <- factor(paste(pheno_i$Group, pheno_i$Histology, sep = "_"))
  
  print(table(pheno_i$Group, pheno_i$Histology))
  
  full_model    <- ~ Batch + GH
  
  reduced_model <- ~ Batch + Group + Histology
  
  m_full <- model.matrix(full_model, pheno_i)
  
  m_red  <- model.matrix(reduced_model, pheno_i)
  
  lrt_res <- NULL; gh_tabs <- NULL
  
  if (qr(m_full)$rank == ncol(m_full) && qr(m_red)$rank == ncol(m_red)) {
    
    dds_i <- DESeqDataSetFromMatrix(
      countData = cnt[rownames(dds), rownames(pheno_i)],
      colData = pheno_i, design = full_model)
    
    dds_lrt <- DESeq(dds_i, test = "LRT", full = m_full, reduced = m_red)
    
    lrt_res <- results(dds_lrt)
    
    lrt_res <- lrt_res[order(lrt_res$padj), ]
    
    message(label, ": LRT genes with padj < 0.1 = ",
            sum(lrt_res$padj < 0.1, na.rm = TRUE))
    
    ggsave(filename = file.path(out_dir, "lrt_pvalue_histogram.pdf"),
           plot = ggplot(data.frame(pvalue = lrt_res$pvalue[!is.na(lrt_res$pvalue)]), aes(pvalue)) +
             geom_histogram(bins = 50, fill = "grey60", colour = "white") +
             labs(title = paste(label, "- LRT (interaction) p-values"), y = "Genes") +
             theme_bw(),
           width = 6, height = 4)
    
    top_int <- head(rownames(lrt_res)[!is.na(lrt_res$padj) & lrt_res$padj < 0.1], 4)
    
    if (length(top_int) > 0)
      ggsave(filename = file.path(out_dir, "lrt_top_interaction_genes.pdf"),
             plot = box_genes(mat_adj, pheno_i, top_int, facet_by_group = TRUE,
                              histo_cols = HISTO_COLS, group_cols = GROUP_COLS),  
             width = 9, height = 8)
    
    dds_w <- DESeq(dds_i)
    
    tabc  <- table(pheno_i$Group, pheno_i$Histology)
    
    pairs <- list(c("DN", "CLA"), c("LCA", "CLA"), c("LCA", "DN"))
    
    cts <- list()
    
    for (g in rownames(tabc)) for (p in pairs)
      
      if (all(p %in% colnames(tabc)) && all(tabc[g, p] >= 5))
        
        cts[[paste(g, p[1], "vs", p[2], sep = "_")]] <-
      c("GH", paste(g, p[1], sep = "_"), paste(g, p[2], sep = "_"))
    
    if (length(cts) > 0) {
      
      gh_raw <- lapply(cts, function(ct) results(dds_w, contrast = ct, alpha = padj_cutoff))
      
      gh_shr <- Map(function(ct, r)
        lfcShrink(dds_w, contrast = ct, res = r, type = "ashr", quiet = TRUE), cts, gh_raw)
      
      gh_tabs <- Map(res_table, gh_raw, gh_shr)
      
      print(sapply(gh_tabs, function(t) nrow(sig_hits(t,padj_cutoff = padj_cutoff, lfc_threshold = lfc_threshold))))
      
      for (nm in names(gh_tabs))
        
        ggsave(filename = file.path(out_dir, paste0("volcano_groupspecific_", nm, ".pdf")),
               plot = volcano_plot(gh_tabs[[nm]], nm, padj_cutoff = padj_cutoff, lfc_threshold = lfc_threshold),
               width = 6.5, height = 5.5)
      
      ggsave(filename = file.path(out_dir, "groupspecific_gene_counts.pdf"),
             plot = count_plot(gh_tabs, paste(label, "- group-specific contrasts"),
                               padj_cutoff = padj_cutoff, lfc_threshold = lfc_threshold),
             width = 7, height = 4)
    }
  } else {
    warning(label, ": interaction design not full rank after dropping small cells; ",
            "LRT skipped")
  }
  
  list(pheno = pheno, dds = dds, raw = raw, shrunk = shrunk, tabs = tabs,
       hits = hits, mat = mat, mat_adj = mat_adj, lrt = lrt_res, gh_tabs = gh_tabs)
}


NMB <- run_dge(train_RNAseq, "NMB", lfc_threshold = 1, padj_cutoff = 0.05)

Non_CPR <- run_dge(test_RNAseq, "non_CPR", lfc_threshold = 1, padj_cutoff = 0.05)

Non_CPR_bc<-run_dge(test_RNAseq_bc, "non_CPR_bc", lfc_threshold = 1, padj_cutoff = 0.05, include_batch =  FALSE)

saveRDS(NMB$hits,     file.path(path.expand(DGE_dir), "NMB_sig_genes.rds"))

saveRDS(Non_CPR$hits, file.path(path.expand(DGE_dir), "non_CPR_sig_genes.rds"))

saveRDS(list(NMB = NMB$tabs, non_CPR = Non_CPR$tabs),
        file.path(path.expand(DGE_dir), "all_DGE_tables.rds"))

# ------ Cross-cohort agreement -----------

for (nm in names(contrast_list))
  ggsave(filename = file.path(DGE_dir, paste0("concordance_", nm, ".pdf")),
         plot = concordance_plot(NMB$tabs[[nm]], Non_CPR$tabs[[nm]], nm, padj_cutoff = 0.05),
         width = 6, height = 5.5)

lrt_n <- as.data.frame(NMB$lrt)
lrt_c <- as.data.frame(Non_CPR$lrt)

shared <- intersect(rownames(lrt_n)[!is.na(lrt_n$padj)],
                    rownames(lrt_c)[!is.na(lrt_c$padj)])

sig_n <- intersect(shared, rownames(lrt_n)[lrt_n$padj < 0.1])

sig_c <- intersect(shared, rownames(lrt_c)[lrt_c$padj < 0.1])

both  <- intersect(sig_n, sig_c)

length(shared); length(sig_n); length(sig_c); length(both)

length(sig_n) * length(sig_c) / length(shared)

tab <- table(NMB = factor(shared %in% sig_n, c(FALSE, TRUE)),
             nonCPR = factor(shared %in% sig_c, c(FALSE, TRUE)))

fisher.test(tab)

cor(lrt_n[shared, "stat"], lrt_c[shared, "stat"], method = "spearman")

common <- intersect(names(NMB$gh_tabs), names(Non_CPR$gh_tabs))
common

for (nm in common)
  print(concordance_plot(NMB$gh_tabs[[nm]], Non_CPR$gh_tabs[[nm]], nm, padj_cutoff = 0.05))


for (nm in common) {
  a <- NMB$gh_tabs[[nm]];  b <- Non_CPR$gh_tabs[[nm]]
  m <- merge(a[, c("gene", "log2FC", "padj")], b[, c("gene", "log2FC", "padj")],
             by = "gene", suffixes = c(".n", ".c"))
  m <- m[complete.cases(m), ]
  
  disc  <- m$padj.n < 0.05 & abs(m$log2FC.n) > 1      
  agree <- sign(m$log2FC.n[disc]) == sign(m$log2FC.c[disc])
  
  cat(nm, "| NMB hits:", sum(disc),
      "| same direction in non-CPR:", round(100 * mean(agree)), "%\n")
}

# -------- LASSO Feature Selection ----------------------

library(glmnet)

RNA_dir<-"~/Thesis_models/RNA_dir/"

LASSO_RNA_dir<-"~/Thesis_models/RNA_dir/LASSO/"

dir.create(LASSO_RNA_dir, showWarnings = FALSE, recursive = TRUE) 

vst_train<-readRDS(file.path(path.expand(RNA_dir), "Train_RNAseq_vst.rds")) 

train_NMB_mat<- readRDS(file.path(path.expand(RNA_dir), "Train_RNAseq.rds"))

train_pheno<-as.data.frame(colData(vst_train))

set.seed(42)

x_train <- as.matrix(train_NMB_mat[,!colnames(train_NMB_mat) %in% "histo"])

if(all(rownames(x_train) == train_pheno$Sample_ID)) { y_train<- as.factor(train_pheno$Histology) 

} else {
  
  train_pheno<-train_pheno[match(rownames(x_train), train_pheno$Sample_ID),]
  
  y_train<- as.factor(train_pheno$Histology) }

cv.lasso <- cv.glmnet(
  x = x_train,
  y = y_train,
  nfolds = 5,
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

saveRDS(selected_features,file.path(path.expand(LASSO_RNA_dir), "Train_RNAseq_Lasso_features.rds"))

cat("Number of features selected by LASSO:", length(selected_features), "\n")

train_df <- data.frame(x_train[, selected_features, drop = FALSE])

saveRDS(train_df, file = file.path(path.expand(LASSO_RNA_dir), "LASSO_RNAseq.rds") )




