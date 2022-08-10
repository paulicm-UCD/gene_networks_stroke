###################################################################################################################
#This script is based on and modified from Jeremy Miller's "Meta-analyses of data 
#from two (or more) microarray data sets" tutorial

#Miller's tutorial and related files can be found at: 
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/JMiller/

#Miller's tutorial is based on:
#Miller JA, Horvath S, Geschwind DH. (2010) Divergence of human and mouse brain 
#transcriptome highlights Alzheimer disease pathways. Proc Natl Acad Sci U S A. 2010 Jul 13;107(28):12698-703.
###################################################################################################################

###################################################################################################################
#In addition to Miller's tutorial, the following resources were used to create this script:

#Horvath, S. Weighted Network Analysis. Applications in Genomics and Systems Biology. Book. 2011.
#https://link.springer.com/book/10.1007/978-1-4419-8819-5

#Steve Horvath's Tutorial "Weighted Gene Co-Expression Network Analysis (WGCNA) R Tutorial"
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/ASPMgene/
#based on:
#Horvath S, Zhang B, Carlson M, Lu KV, Zhu S, Felciano RM, Laurance MF, Zhao W, Shu, Q, Lee Y, Scheck AC, 
#Liau LM, Wu H, Geschwind DH, Febbo PG, Kornblum HI, Cloughesy TF, Nelson SF, Mischel PS (2006) "Analysis of 
#Oncogenic Signaling Networks in Glioblastoma Identifies ASPM as a Novel Molecular Target", PNAS 
#November 14, 2006 | vol. 103 | no. 46 | 17402-17407

#Langfelder, P. Signed vs. Unsigned Topological Overlap Matrix. Technical Report. 2013.
#https://www.researchgate.net/file.PostFileLoader.html?id=57bdeaad40485404eb0753d4&assetKey=AS%3A398680193552384%401472064173254

#Langfelder and Horvath's Tutorial "Network analysis of liver expression data from female mice: 
#finding modules related to body weight"
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-auto.pdf
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-blockwise.pdf
#based on:
#Ghazalpour A, Doss S, Zhang B, Wang S, Plaisier C, et al. (2006) Integrating Genetic and Network Analysis to Characterize 
#Genes Related to Mouse Weight. PLOS Genetics 2(8): e130. https://doi.org/10.1371/journal.pgen.0020130

#Additional Resource:
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
###################################################################################################################

###################################################################################################################
#The WGCNA package is described in:
#Langfelder P, Horvath S (2008). "WGCNA: an R package for weighted correlation network analysis." 
#BMC Bioinformatics, 559. https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559.
#and
#Langfelder P, Horvath S (2012). "Fast R Functions for Robust Correlations and Hierarchical Clustering." 
#Journal of Statistical Software, 46(11), 1-17. https://www.jstatsoft.org/v46/i11/.
###################################################################################################################



#The parameters included by default in this script were used to generate the Monocyte (MON) network
#Parameters used for Neutrophil (NEU) and Whole Blood (WB) networks commented in
#Additionally, output files are named for MON Network analysis, and will need to be changed for other networks



##########ONLY INSTALL PACKAGES IF THIS IS THE FIRST TIME USING WGCNA ON THIS MACHINE##########
###############################################################################################
install.packages(c("WGCNA", "dynamicTreeCut", "flashClust", "ggplot2",
"Hmisc", "MASS", "class", "cluster", "survival"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GO.db")

source("https://bioconductor.org/biocLite.R")
biocLite("impute")

source("https://bioconductor.org/biocLite.R")
biocLite("preprocessCore")
###############################################################################################



####LOAD REQUIRED LIBRARIES####

library(dynamicTreeCut)
library(flashClust)
library(ggplot2)
library(Hmisc)
library(WGCNA)
require(WGCNA)
library(MASS)
library(class)
library(cluster)
library(impute)
library(survival)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(4)

#############################################################################
#############################################################################
###############CONDUCT ANALYSIS ON COMBINED GROUPS (EXP+CTRL)################
#############################################################################
#############################################################################

#####LOAD DATA(COMB)#####

datCOMBtemp = read.csv("EXPRESSION_DATA_FILE_NAME.csv", sep=",", row.names=1) ##########INSERT FILE NAME FOR DESIRED DATASET / NETWORK
datCOMB <- as.matrix(datCOMBtemp) #changes .csv into matrix from a data frame

#check for missing values
gsg = goodSamplesGenes(datCOMB, verbose=5)
gsg$allOK #if output=TRUE, continue; otherwise, refer to tutorial I.1.a.

datCOMBp = datCOMB

class(datCOMBp) #confirm in matrix/array form




#####DETERMINE SOFT-THRESHOLDING POWER(COMBINED)#####

datExprdataCOMB = t(datCOMBp)

powers1 = c(seq(1,10,by=1),seq(12,30,by=2))
RpowerTable = pickSoftThreshold(datExprdataCOMB, powerVector=powers1, networkType = "signed", verbose=5)[[2]]
#For above line: Col1=power beta, col2=scale free topology fitting index R^2
#col3=slope of fitting line, col4=fitting index for truncated exponential scale free model

pdf("Soft-thresholding Power_DynamicsOfIS_MON_Network.pdf", height=10, width=9) 
cex1=0.7
par(mfrow=c(1,2))
plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],xlab="Soft Threshold (power)",
ylab="Scale Free Topology Model Fit, signed R^2", main="DynamicsOfIS_MON_Network", type="n")
text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],
labels = powers1 ,cex=cex1, col="red")
abline(h=0.9, col="green")
abline(h=0.8, col="red")
plot(RpowerTable[,1], RpowerTable[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity",
main="DynamicsOfIS_MON_Network", type="n")
text(RpowerTable[,1], RpowerTable[,5], labels=powers1, cex=cex1, col="red")
grid(nx = NULL, ny = 70, col = "gray", lty = "dotted")
abline(h=100, col="red")
abline(h=200, col="green")
dev.off()

beta1= 14 ##beta 14 used for MON network; beta 8 used for NEU network; beta 16 used for WB network

k.dataCOMB = softConnectivity(datExprdataCOMB, power=beta1, type="signed")-1




##########MODULE DETECTION: SELECTING ONLY MOST CONNECTED GENES(COMB)##########

kCut =  14955 ##full dataset 14 used for MON network; full dataset 13921 used for NEU network; full dataset 15360 used for WB network
kRank = rank(-k.dataCOMB)
vardataCOMB = apply(datExprdataCOMB, 2, var)
restk = kRank <= kCut & vardataCOMB>0

sum(restk)

datCOMBg = t(datExprdataCOMB[,restk])

kRankDF = as.matrix(kRank)
geneListCOMBtemp = t(colnames(datExprdataCOMB))
geneListCOMBtemp2 = t(geneListCOMBtemp)
geneListCOMB = cbind(geneListCOMBtemp2, kRankDF)
colnames(geneListCOMB) = c("Gene", "Rank")

write.csv(geneListCOMB, "DynamicsOfIS_MON_Network_Genes ranked by interconnectivity.csv")##########THIS LINE GENERATES A .CSV WHERE EACH ROW IS A GENE NAME AND ITS RANKING FOR INTERCONNECTIVITY
write.csv(datCOMBg, "DynamicsOfIS_MON_Network_ALL_mostInterconn.csv", row.names=TRUE)




##########RUN WGCNA(COMBINED)##########

ADJdataCOMB = adjacency(datExprdataCOMB[,restk], power = beta1, type="signed")
dissTOMdataCOMB = TOMdist(ADJdataCOMB, TOMType="signed")
hierTOMdataCOMB = hclust(as.dist(dissTOMdataCOMB), method = "average");

pdf("Dendrogram_DynamicsOfIS_MON_Network.pdf", height=8, width=12)
par(mfrow = c(1,1))
plot(hierTOMdataCOMB, labels=F, main="Dendrogram: ALL Most Interconnected Genes in DynamicsOfIS_MON_Network") 
dev.off()




##########DETERMINE MODULES(COMBINED)##########

minModuleSize = 50;

# ##...deepSplit = 0 ##
# 
# dynamicModsCOMB = cutreeDynamic(dendro = hierTOMdataCOMB, distM = dissTOMdataCOMB,
# 	deepSplit = FALSE, pamRespectsDendro = FALSE, pamStage = FALSE,
# 	minClusterSize = minModuleSize, method = "tree");
# #USE 'pamStage = TRUE' IF LARGE NUMBER OF GREY GENES [ONLY FOR method = "hybrid"]
# #method = "tree" ONLY ALLOWS pamRespectsDendro = pamStage = FALSE
# #FOR method = "tree" LET deepSplit = TRUE (1) or FALSE (0)
# table(dynamicModsCOMB)
# 
# dynamicColorsCOMB = labels2colors(dynamicModsCOMB)
# modulesCOMB = dynamicColorsCOMB
# table(dynamicColorsCOMB)
# 
# pdf("Module_Choices_DynamicsOfIS_MON_Network_DynamicTreeCut_ds0.pdf", height = 8, width = 12)
# par(mfrow = c(1,1))
# plotDendroAndColors(hierTOMdataCOMB, dynamicColorsCOMB, "Dynamic Tree Cut",
# dendroLabels = FALSE, hang = 0.03,
# addGuide = TRUE, guideHang = 0.05,
# main = "Dendrogram_DynamicsOfIS_MON_Network (ds0)")
# dev.off()

##...deepSplit = 1 ##

dynamicModsCOMB = cutreeDynamic(dendro = hierTOMdataCOMB, distM = dissTOMdataCOMB,
	deepSplit = TRUE, pamRespectsDendro = FALSE, pamStage = FALSE,
	minClusterSize = minModuleSize, method = "tree");
#USE 'pamStage = TRUE' IF LARGE NUMBER OF GREY GENES [ONLY FOR method = "hybrid"]
#method = "tree" ONLY ALLOWS pamRespectsDendro = pamStage = FALSE
#FOR method = "tree" LET deepSplit = TRUE (1) or FALSE (0)
table(dynamicModsCOMB)

dynamicColorsCOMB = labels2colors(dynamicModsCOMB)
modulesCOMB = dynamicColorsCOMB
table(dynamicColorsCOMB)

pdf("Module_Choices_DynamicsOfIS_MON_Network_DynamicTreeCut_ds1.pdf", height = 8, width = 12)
par(mfrow = c(1,1))
plotDendroAndColors(hierTOMdataCOMB, dynamicColorsCOMB, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Dendrogram_DynamicsOfIS_MON_Network (ds1)") 
dev.off()

##USING THE PDF'S GENERATED FOR SOFT-THRESHOLDING POWER, DETERMINE WHICH deepSplit VALUE IS DESIRED AND RUN THAT BLOCK MOST RECENTLY BEFORE PROCEEDING##

##deepSplit 1 was used for MON, NEU, and WB networks##



##DEFINE TERMS FOR FUTURE PROCESSES(COMB)##

adjCOMB = ADJdataCOMB;
diag(adjCOMB)=0
dissTOMACOMB = dissTOMdataCOMB;
geneTreeCOMB = hierTOMdataCOMB;




#####CALCULATE PRINCIPLE COMPONENTS FOR VISUALIZATIONS(COMB)#####

PCsCOMB = moduleEigengenes(t(datCOMBg), colors=modulesCOMB)
ME_COMB = PCsCOMB$eigengenes
distPCCOMB = 1-abs(cor(ME_COMB, use="p"))
distPCCOMB = ifelse(is.na(distPCCOMB), 0, distPCCOMB)
pcTreeCOMB = hclust(as.dist(distPCCOMB), method="a")
MDS_COMB = cmdscale(as.dist(distPCCOMB),2)
colorsCOMB = names(table(modulesCOMB))

ME_COMB_Temp = ME_COMB
rownames(ME_COMB_Temp) = colnames(datCOMBg)
write.csv(ME_COMB_Temp, "DynamicsOfIS_MON_Network_ModuleEigengeneValues.csv")##########THIS FILE IS WHAT WAS USED TO DETERMINE MODULE SIGNIFICANCE FOR Dx, OTHER CLINICAL PARAMETERS




#####MODULE MEMBERSHIP (kME) AND COMPARING NETWORKS(COMBINED)#####

geneModuleMembershipCOMB = signedKME(t(datCOMBg), ME_COMB)
colnames(geneModuleMembershipCOMB) = paste("PC",colorsCOMB, ".cor", sep="");
MMPvalueCOMB = corPvalueStudent(as.matrix(geneModuleMembershipCOMB), dim(datCOMBg)[[2]]);
colnames(MMPvalueCOMB) = paste("PC", colorsCOMB, ".pval", sep="");

Gene = rownames(datCOMBg)
kMEtableCOMB = cbind(Gene, Gene, modulesCOMB)
for(i in 1:length(colorsCOMB))
	kMEtableCOMB = cbind(kMEtableCOMB, geneModuleMembershipCOMB[,i], MMPvalueCOMB[,i])
colnames(kMEtableCOMB) = c("PSID", "Gene", "Module", sort(c(colnames(geneModuleMembershipCOMB),
colnames(MMPvalueCOMB))))

write.csv(kMEtableCOMB, "kMEtable_DynamicsOfIS_MON_Network.csv", row.names=FALSE) 




#####GENERAL VISUALIZATIONS######

pdf("ModuleEigengeneVisualizations_DynamicsOfIS_MON_Network.pdf", height=6, width=6)
par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)

plot(pcTreeCOMB, xlab="", ylab="", main="", sub="")
plot(MDS_COMB, col=colorsCOMB, main="MDS Plot", cex=2, pch=19)

for(which.module in names(table(modulesCOMB))){
	ME = ME_COMB[, paste("ME", which.module, sep="")]
	barplot(ME, col=which.module, main=which.module, cex.main=2,
	ylab="Eigengene Expression DynamicsOfIS_MON_Network", xlab="Sample") 
};

dev.off()




####MODULE NAME COLOR LEGEND####

modCOMB = unique(modulesCOMB)

pdf("DynamicsOfIS_MON_Network_ModuleColorLegends.pdf", height=26, width=8)
par(mfrow=c(1,1))
plot(1,1,legend("center", c(modCOMB), col=c(modCOMB), lwd=10))
dev.off()
#AN ERROR WILL BE DISPLAYED FOR THESE PLOTS, BUT THE OUTPUT PDF WILL STILL PROVIDE LEGENDS FOR THE COLORS OF THE MODULES


#########################################################
#########################################################
##########SCRIPT TO MODULE MEMBERSHIP####################

allModuleColors = colorsCOMB
geneModuleAssignments = modulesCOMB
geneNames = rownames(datCOMBg)

###############################################################################

for (co in allModuleColors){
  
  cIndex = (geneModuleAssignments==co)
  
  genes = geneNames
  
  modGenes = genes[cIndex]
  
  fn = paste(co, "_moduleGeneNames.csv", sep="")
  
  write.table(modGenes, file=fn, row.names=F, col.names=F)
}

#######################################################################################
#######################################################################################
##########SCRIPT TO OUTPUT DATA FOR NETWORK VISUALIZATION IN VisANT####################
#####MUST USE EXISTING DATA STORED IN WORKSPACE TO OUTPUT VisANT FILES

source("tutorialFunctions.R") 

#THE NEXT TWO LINES OUTPUT THE VisANT FILES AND HUB LISTS FOR THE MODULES
for (co in colorsCOMB[colorsCOMB!="grey"])
visantPrepOverall(modulesCOMB, co, t(datCOMBg), rownames(datCOMBg), 10000, beta1, TRUE)

#visantPrepOverall <- function(colorFinal2, moduleColor, datExprrest, genes, numInt, power, signed=FALSE)

## This file returns the overall 10,000 strongest connections within a module for one group of subjects

## USER inputs
# colorFinal2 = vector with module association for each probe
# moduleColor = color of the module... should be a member of colorFinal2
# datExprrest = expression data for the genes corresponding to cIndex
# genes = list of genes that correspond to the probes
# numInt = number of interactions to output to the visant file


#################################################################################
#Testing for module relationships to Diagnosis and other clinical parameters was 
#conducted in another program - Partek Genomics Suite®
#Spearman Correlations and Kruskal-Wallis tests were used to determine continuous and
#categorical paramteres' association with the ModuleEigengeneValues file (noted above)
#################################################################################