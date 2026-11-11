
#=========================================================================================================#
#------------------------- Similarity Network Fusion------------------------------------------------------
#=========================================================================================================#

#This is a template script,and in this case, Methylation data and CNV data are being integrated 

require(SNFtool) 
require(dplyr)
require(igraph) 


NMB_450k<-read.csv("~/NMB_Final_Cohort.csv") 

NMB_450k<-NMB_450k[NMB_450k$Platform == "450K",] 

NMB_meth<-as.data.frame(readRDS("~/NMB/NMBTC_450k_betas.rds"))

NMB_meth<-as.data.frame(t(NMB_meth)) 

NMB_meth$histo<-as.factor(NMB_450k$CPR_Histology)

NMB_meth$mol.grp<-as.factor(NMB_450k$Sub.group.called.7grp.) 

NMB_meth<-as.data.frame(NMB_meth[!NMB_meth$histo %in% "NOS",])   

NMB_meth$histo<-droplevels(NMB_meth$histo) 

NMB_meth$mol.grp<-gsub("SHH_Inf", "SHH", NMB_meth$mol.grp)

NMB_meth$mol.grp<-gsub("SHH_Old", "SHH", NMB_meth$mol.grp) 

NMB_meth<-as.data.frame(NMB_meth[NMB_meth$mol.grp == "SHH",]) 

NMB_meth<-as.data.frame(NMB_meth[NMB_meth$histo %in% c("CLA","LCA"),])

NMB_meth$histo<-droplevels(NMB_meth$histo) 

NMB_meth<-NMB_meth[,!colnames(NMB_meth) %in% "mol.grp"] 

BINS<-read.table("~/NMB/Conumee2.0/NMB_bins.igv", header = TRUE, sep = "")

bins<-as.data.frame(t(BINS[,!colnames(BINS) %in% c( "Chromosome", "Start", "End")])) 

colnames(bins)<-bins[1,] 

bins<-as.data.frame(bins[-1,])

bins<-as.data.frame(sapply(bins, function(x) as.numeric(x))) 

minMax<-function(x) { x-min(x)/ max(x) - min(x)} 

scaled_bins<-sapply(bins,FUN = minMax)

NMB_bins<-as.data.frame(scaled_bins) 

NMB_bins$histo<-as.factor(NMB_450k$CPR_Histology) 

NMB_bins$mol.grp<-as.factor(NMB_450k$Sub.group.called.7grp.)

NMB_bins$mol.grp<-gsub("SHH_Inf", "SHH", NMB_bins$mol.grp)

NMB_bins$mol.grp<-gsub("SHH_Old", "SHH", NMB_bins$mol.grp) 

NMB_bins<-NMB_bins[!NMB_bins$histo %in% "NOS",]  

NMB_bins$histo<-droplevels(NMB_bins$histo) 

NMB_bins<-as.data.frame(NMB_bins[NMB_bins$mol.grp == "SHH",])

NMB_bins<-as.data.frame(NMB_bins[NMB_bins$histo %in% c("CLA","LCA"),])

NMB_bins$histo<-droplevels(NMB_bins$histo) 

NMB_bins<-NMB_bins[,!colnames(NMB_bins) %in% "mol.grp"] 

rownames(NMB_bins)<-as.factor(rownames(NMB_meth))

NMB_meth<-standardNormalization(NMB_meth[,!colnames(NMB_meth) %in% "histo"])

NMB_meth<-as.data.frame(NMB_meth)

NMB_meth$histo<-as.factor(NMB_bins$histo) 

NMB_bins<-standardNormalization(NMB_bins[,!colnames(NMB_bins) %in% "histo"])

NMB_bins<-as.data.frame(NMB_bins)

NMB_bins$histo<-as.factor(NMB_meth$histo) 

#-------------------Distance matrices--------------------------------

dist1_NMB <- as.matrix(dist(NMB_meth[,!colnames(NMB_meth) %in% "histo"]))

#write.csv(as.data.frame(dist1_NMB),"~/SHH_SNF/SHH_NMB_meth_distmatrix.csv")

dist2_NMB <- as.matrix(dist(NMB_bins[,!colnames(NMB_bins) %in% "histo"]))  

#write.csv(as.data.frame(dist2_NMB),"~/SHH_SNF/SHH_NMB_bins_distmatrix.csv")

#---------------------------------Affinity Matrices(Similarity Matrices)---------------------------------

W1_NMB <- affinityMatrix(dist1_NMB, K = 5, sigma = 0.5) 

#write.csv(as.data.frame(W1_NMB),"~/SHH_SNF/NMB_meth_affinitymatrix.csv") 


W2_NMB <- affinityMatrix(dist2_NMB, K = 5, sigma = 0.5)

#write.csv(as.data.frame(W2_NMB),"~/SHH_SNF/NMB_bins_affinitymatrix.csv") 

#-----------------------Heatmap cluster visualisation----------------------------------------------

#jpeg("~/SHH_SNF/SHH_NMB_bins_clusters.jpeg")  

SHH_NMB_bins_clusters<-displayClustersWithHeatmap(W2_NMB, spectralClustering(W2_NMB, K = 2)) 

#dev.off() 

#-----------------------------Network Fusion----------------------------------------------

SNF_NMB = SNF(list(W1_NMB,W2_NMB), 5, 20) 

#write.csv(as.data.frame(SNF_NMB),"~/SHH_SNF/SHH_NMBMeth+bins.csv") 

#-------------------------- Fused Heatmap cluster visualisation---------------------------------- 

#jpeg("~/SHH_SNF/SHH_NMB_SNF_Meth+bins_clusters.jpeg")  

SHH_NMB_SNF_Methbins<-displayClustersWithHeatmap(SNF_NMB, spectralClustering(SNF_NMB, K = 3)) 

#dev.off() 


#------------------------------Extracting features as graph embeddings-----------------------------

g <- graph_from_adjacency_matrix(SNF_NMB, mode = "undirected", weighted = TRUE)

#plot(g, vertex.size = 5, vertex.label = NA, main = "SHH_NMB Meth+bins Fused Similarity Network")


#-------------------------------------Centrality measures-------------------------------------------- 

degree_features <- degree(g)  
betweenness_features <- betweenness(g)  
closeness_features <- closeness(g) 

node_features_NMB_SHH <- data.frame(degree = degree_features, 
                                betweenness = betweenness_features, 
                                closeness = closeness_features, histo = as.factor(NMB_bins$histo))



#-------------------------Cavalli_SNF------------------------------------

cavalli_meth<-as.data.frame(readRDS("~/450k/Cavalli/Cavalli_450k_betas.rds"))

cavalli_meth<-as.data.frame(t(cavalli_meth))   

cavalli_labels<-read.csv("~/Cavalli-test.csv",header = TRUE) 

cavalli_labels$Histology<-gsub("Classic","CLA",cavalli_labels$Histology)

cavalli_labels$Histology<-gsub("Desmoplastic","DN",cavalli_labels$Histology)

#cavalli_labels$Histology<-gsub("MBEN","DN",cavalli_labels$Histology)

#cavalli_labels$Sub.group<-gsub("Group3","Grp3",cavalli_labels$Sub.group)

cavalli_meth$histo<-as.factor(cavalli_labels$Histology)

cavalli_meth$mol.grp<-as.factor(cavalli_labels$Sub.group) 


cavalli_meth<-as.data.frame(cavalli_meth[cavalli_meth$mol.grp == "WNT",]) 

cavalli_meth<-as.data.frame(cavalli_meth[!is.na(cavalli_meth$histo),])

#cavalli_meth<-cavalli_meth[cavalli_meth$histo != "DN",]

#cavalli_meth$histo<-droplevels(cavalli_meth$histo) 

cavalli_meth<-cavalli_meth[,!colnames(cavalli_meth) %in% "mol.grp"] 


cavalli_BINS<-as.data.frame(read.table("~/450k/Cavalli/Conumee2.0/cavalli_bins.igv",sep = "",header = TRUE)) 

cavalli_bins<-as.data.frame(t(cavalli_BINS[,!colnames(cavalli_BINS) %in% c( "Chromosome", "Start", "End")]))


colnames(cavalli_bins)<-cavalli_bins[1,] 

cavalli_bins<-as.data.frame(cavalli_bins[-1,])


rownames(cavalli_bins)<-as.factor(rownames(cavalli_bins))


colnames(cavalli_bins)<-as.factor(colnames(cavalli_bins)) 


cavalli_bins<-as.data.frame(sapply(cavalli_bins, function(x) as.numeric(x))) 

minMax<-function(x) { x-min(x)/ (max(x) - min(x)) } 

cavallibins<-as.data.frame(sapply(cavalli_bins, function(x) as.numeric (x)))  


cavalli_bins<-as.data.frame(sapply(cavallibins,FUN = minMax))  

pheno<-read.csv("~/Cavalli-test.csv") 

pheno$Histology<-gsub("Classic","CLA",pheno$Histology)

pheno$Histology<-gsub("Desmoplastic","DN",pheno$Histology)

#pheno$Histology<-gsub("MBEN","DN",pheno$Histology)

#pheno$Sub.group<-gsub("Group3","Grp3",pheno$Sub.group)

cavalli_bins$histo<-as.factor(pheno$Histology) 

cavalli_bins$histo<-droplevels(cavalli_bins$histo) 

cavalli_bins$mol.grp<-as.factor(pheno$Sub.group)


cavalli_bins$histo<-as.factor(cavalli_bins$histo) 

cavalli_bins<-as.data.frame(cavalli_bins[cavalli_bins$mol.grp == "WNT",])

cavalli_bins<- as.data.frame(cavalli_bins[!is.na(cavalli_bins$histo),])  

#cavalli_bins$histo<-droplevels(cavalli_bins$histo) 

rownames(cavalli_bins)<-as.factor(rownames(cavalli_meth)) 

cavalli_bins<-cavalli_bins[,!colnames(cavalli_bins) %in% "mol.grp"] 

#cavalli_bins<-cavalli_bins[cavalli_bins$histo != "DN",]

cavalli_bins$histo<-droplevels(cavalli_bins$histo)

cavalli_bins<-standardNormalization(cavalli_bins[,!colnames(cavalli_bins) %in% "histo"])

cavalli_bins<-as.data.frame(cavalli_bins) 

cavalli_bins$histo<-as.factor(cavalli_meth$histo) 

cavalli_meth<-standardNormalization(cavalli_meth[,!colnames(cavalli_meth) %in% "histo"])

cavalli_meth<-as.data.frame(cavalli_meth)  

cavalli_meth$histo<-as.factor(cavalli_bins$histo)  

rownames(cavalli_bins)<-as.factor(rownames(cavalli_meth)) 



#-----------------------------Distance matrices---------------------------------------

dist1_cavalli <- as.matrix(dist(cavalli_meth[,!colnames(cavalli_meth) %in% "histo"]))


#write.csv(as.data.frame(dist1_cavalli),"~/SHH_SNF/SHH_cavalli_meth_distmatrix.csv")


dist2_cavalli <- as.matrix(dist(cavalli_bins[,!colnames(cavalli_bins) %in% "histo"]))

#write.csv(as.data.frame(dist2_cavalli),"~/SHH_SNF/SHH_cavalli_bins_distmatrix.csv")

#-----------------------Affinity Matrices(Similarity Matrices)--------------------------

W1_cavalli <- affinityMatrix(dist1_cavalli, K = 10, sigma = 0.5)

#write.csv(W1_cavalli,"~/SHH_SNF/SHH_cavalli_meth_affinitymatrix.csv") 

W2_cavalli <- affinityMatrix(dist2_cavalli, K = 10, sigma = 0.5)

#write.csv(W2_cavalli,"~/SHH_SNF/SHH_cavalli_bins_affinitymatrix.csv") 

#-------------------------------Heatmap cluster visualisation---------------------------------------

#jpeg("~/SHH_SNF/SHH_cavalli_meth_clusters.jpeg")  

SHH_cavalli_meth_clusters<-displayClustersWithHeatmap(W1_cavalli, spectralClustering(W1_cavalli, K = 4)) 

#jpeg("~/SHH_SNF/SHH_cavalli_bins_clusters.jpeg")  

SHH_cavalli_bins_clusters<-displayClustersWithHeatmap(W2_cavalli, spectralClustering(W2_cavalli, K = 4)) 

#dev.off() 

#-------------------------------------Network Fusion--------------------------------------

SNF_cavalli = SNF(list(W1_cavalli,W2_cavalli), 10, 20)

#write.csv(SNF_cavalli,"~/SHH_SNF/SHH_cavalliMeth+bins.csv") 

#---------------------------Fused Heatmap cluster visualisation-----------------------------

#jpeg("~/SHH_SNF/SHH_cavalli_SNF_Meth+bins_clusters.jpeg")  

SHH_cavalli_SNFMethbins<-displayClustersWithHeatmap(SNF_cavalli, spectralClustering(SNF_cavalli, K = 4)) 

#dev.off() 

#-----------------------------Extracting features as graph embeddings-------------------------

g <- graph_from_adjacency_matrix(SNF_cavalli, mode = "undirected", weighted = TRUE)

#p<-plot(g, vertex.size = 5, vertex.label = NA, main = "SHH_Cavalli Meth+bins Fused Similarity Network")


#-----------------------------------Centrality measures-----------------------------------------

degree_features <- degree(g)  
betweenness_features <- betweenness(g)  
closeness_features <- closeness(g) 

node_features_cavalli <- data.frame(degree = degree_features, 
                                    betweenness = betweenness_features, 
                                    closeness = closeness_features, histo = as.factor(cavalli_meth$histo))




