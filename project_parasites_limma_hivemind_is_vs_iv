## This file codes for Limma analysis of Pst under attack by 
setwd("C:/Users/u6047942/JWPHD/project_parasites/bioinformatics/R-studio_transcriptomics")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("tximport")

#importing the Salmon files
library(edgeR)
library(limma)
library(tximport)

## Hyperparasite attack conditions have had penx reads removed by aligning to the penx genome "10kb8QCanu". All Pst reads quantified by aligning to Ritas Pst 104E- genome
FRF_insitu1 <-tximport("./salmon_outputs/FRF_salmon_insitu_pst_removed1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
FRF_insitu2 <-tximport("./salmon_outputs/FRF_salmon_insitu_pst_removed2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
FRF_insitu3 <-tximport("./salmon_outputs/FRF_salmon_insitu_pst_removed3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
FRF_invitro1 <-tximport("./salmon_outputs/FRF_salmon_invitro1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
FRF_invitro2 <-tximport("./salmon_outputs/FRF_salmon_invitro2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
FRF_invitro3 <-tximport("./salmon_outputs/FRF_salmon_invitro3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)

PBLUE_insitu1 <-tximport("./salmon_outputs/PBLUE_salmon_insitu_pst_removed1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
PBLUE_insitu2 <-tximport("./salmon_outputs/PBLUE_salmon_insitu_pst_removed2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
PBLUE_insitu3 <-tximport("./salmon_outputs/PBLUE_salmon_insitu_pst_removed3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
PBLUE_invitro1 <-tximport("./salmon_outputs/PBLUE_salmon_invitro1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
PBLUE_invitro2 <-tximport("./salmon_outputs/PBLUE_salmon_invitro2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
PBLUE_invitro3 <-tximport("./salmon_outputs/PBLUE_salmon_invitro3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)

SRM_insitu1 <- tximport("./salmon_outputs/SRM_salmon_insitu_pst_removed1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRM_insitu2 <- tximport("./salmon_outputs/SRM_salmon_insitu_pst_removed2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRM_insitu3 <- tximport("./salmon_outputs/SRM_salmon_insitu_pst_removed3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRM_invitro1 <- tximport("./salmon_outputs/SRM_salmon_invitro1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRM_invitro2 <- tximport("./salmon_outputs/SRM_salmon_invitro2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRM_invitro3 <- tximport("./salmon_outputs/SRM_salmon_invitro3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)

PMUD_insitu1 <- tximport("./salmon_outputs/PMUD_salmon_insitu_pst_removed1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
PMUD_insitu2 <- tximport("./salmon_outputs/PMUD_salmon_insitu_pst_removed2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
PMUD_insitu3 <- tximport("./salmon_outputs/PMUD_salmon_insitu_pst_removed3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
PMUD_invitro1 <- tximport("./salmon_outputs/PMUD_salmon_invitro1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
PMUD_invitro2 <- tximport("./salmon_outputs/PMUD_salmon_invitro2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
PMUD_invitro3 <- tximport("./salmon_outputs/PMUD_salmon_invitro3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)

SRB_insitu1 <- tximport("./salmon_outputs/SRB_salmon_insitu_pst_removed1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRB_insitu2 <- tximport("./salmon_outputs/SRB_salmon_insitu_pst_removed2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRB_insitu3 <- tximport("./salmon_outputs/SRB_salmon_insitu_pst_removed3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRB_invitro1 <- tximport("./salmon_outputs/SRB_salmon_invitro1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRB_invitro2 <- tximport("./salmon_outputs/SRB_salmon_invitro2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRB_invitro3 <- tximport("./salmon_outputs/SRB_salmon_invitro3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)

PCB_insitu1 <- tximport("./salmon_outputs/PCB_salmon_insitu_pst_removed1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
PCB_insitu2 <- tximport("./salmon_outputs/PCB_salmon_insitu_pst_removed2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
PCB_insitu3 <- tximport("./salmon_outputs/PCB_salmon_insitu_pst_removed3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
PCB_invitro1 <- tximport("./salmon_outputs/PCB_salmon_invitro1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
PCB_invitro2 <- tximport("./salmon_outputs/PCB_salmon_invitro2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
PCB_invitro3 <- tximport("./salmon_outputs/PCB_salmon_invitro3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)

SRSh_insitu1 <- tximport("./salmon_outputs/SRSh_salmon_insitu_pst_removed1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRSh_insitu2 <- tximport("./salmon_outputs/SRSh_salmon_insitu_pst_removed2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRSh_insitu3 <- tximport("./salmon_outputs/SRSh_salmon_insitu_pst_removed3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRSh_invitro1 <- tximport("./salmon_outputs/SRSh_salmon_invitro1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRSh_invitro2 <- tximport("./salmon_outputs/SRSh_salmon_invitro2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRSh_invitro3 <- tximport("./salmon_outputs/SRSh_salmon_invitro3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)

SRSi_insitu1 <- tximport("./salmon_outputs/SRSi_salmon_insitu_pst_removed1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRSi_insitu2 <- tximport("./salmon_outputs/SRSi_salmon_insitu_pst_removed2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRSi_insitu3 <- tximport("./salmon_outputs/SRSi_salmon_insitu_pst_removed3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRSi_invitro1 <- tximport("./salmon_outputs/SRSi_salmon_invitro1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRSi_invitro2 <- tximport("./salmon_outputs/SRSi_salmon_invitro2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)
SRSi_invitro3 <- tximport("./salmon_outputs/SRSi_salmon_invitro3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"), dropInfReps = TRUE)

countmatrix_FRF_insitu_against_invitro<-data.frame(FRF_insitu1$counts, FRF_insitu2$counts, FRF_insitu3$counts, FRF_invitro1$counts, FRF_invitro2$counts, FRF_invitro3$counts)
row.names(countmatrix_FRF_insitu_against_invitro)<-row.names(FRF_insitu1$counts)
plotMDS(countmatrix_FRF_insitu_against_invitro)

countmatrix_PBLUE_insitu_against_invitro<-data.frame(PBLUE_insitu1$counts, PBLUE_insitu2$counts, PBLUE_insitu3$counts, PBLUE_invitro1$counts, PBLUE_invitro2$counts, PBLUE_invitro3$counts)
row.names(countmatrix_PBLUE_insitu_against_invitro)<-row.names(PBLUE_insitu1$counts)
plotMDS(countmatrix_PBLUE_insitu_against_invitro)

countmatrix_SRM_insitu_against_invitro<-data.frame(SRM_insitu1$counts, SRM_insitu2$counts, SRM_insitu3$counts, SRM_invitro1$counts, SRM_invitro2$counts, SRM_invitro3$counts)
row.names(countmatrix_SRM_insitu_against_invitro)<-row.names(SRM_insitu1$counts)
plotMDS(countmatrix_SRM_insitu_against_invitro)

countmatrix_PMUD_insitu_against_invitro<-data.frame(PMUD_insitu1$counts, PMUD_insitu2$counts, PMUD_insitu3$counts, PMUD_invitro1$counts, PMUD_invitro2$counts, PMUD_invitro3$counts)
row.names(countmatrix_PMUD_insitu_against_invitro)<-row.names(PMUD_insitu1$counts)
plotMDS(countmatrix_PMUD_insitu_against_invitro)

countmatrix_SRB_insitu_against_invitro<-data.frame(SRB_insitu1$counts, SRB_insitu2$counts, SRB_insitu3$counts, SRB_invitro1$counts, SRB_invitro2$counts, SRB_invitro3$counts)
row.names(countmatrix_SRB_insitu_against_invitro)<-row.names(SRB_insitu1$counts)
plotMDS(countmatrix_SRB_insitu_against_invitro)

countmatrix_PCB_insitu_against_invitro<-data.frame(PCB_insitu1$counts, PCB_insitu2$counts, PCB_insitu3$counts, PCB_invitro1$counts, PCB_invitro2$counts, PCB_invitro3$counts)
row.names(countmatrix_PCB_insitu_against_invitro)<-row.names(PCB_insitu1$counts)
plotMDS(countmatrix_PCB_insitu_against_invitro)

countmatrix_SRSh_insitu_against_invitro<-data.frame(SRSh_insitu1$counts, SRSh_insitu2$counts, SRSh_insitu3$counts, SRSh_invitro1$counts, SRSh_invitro2$counts, SRSh_invitro3$counts)
row.names(countmatrix_SRSh_insitu_against_invitro)<-row.names(SRSh_insitu1$counts)
plotMDS(countmatrix_SRSh_insitu_against_invitro)

countmatrix_SRSi_insitu_against_invitro<-data.frame(SRSi_insitu1$counts, SRSi_insitu2$counts, SRSi_insitu3$counts, SRSi_invitro1$counts, SRSi_invitro2$counts, SRSi_invitro3$counts)
row.names(countmatrix_SRSi_insitu_against_invitro)<-row.names(SRSi_insitu1$counts)
plotMDS(countmatrix_SRSi_insitu_against_invitro)

#import design matrix (generated seperately in excel)
design_insitu_against_invitro<-read.csv("./design_insitu_against_invitro.csv")

#filtering genes by expression
dge_FRF_insitu_against_invitro <- DGEList(counts=countmatrix_FRF_insitu_against_invitro)
dge_PBLUE_insitu_against_invitro <- DGEList(counts=countmatrix_PBLUE_insitu_against_invitro)
dge_SRM_insitu_against_invitro <- DGEList(counts=countmatrix_SRM_insitu_against_invitro)
dge_PMUD_insitu_against_invitro <- DGEList(counts=countmatrix_PMUD_insitu_against_invitro)
dge_SRB_insitu_against_invitro <- DGEList(counts=countmatrix_SRB_insitu_against_invitro)
dge_PCB_insitu_against_invitro <- DGEList(counts=countmatrix_PCB_insitu_against_invitro)
dge_SRSh_insitu_against_invitro <- DGEList(counts=countmatrix_SRSh_insitu_against_invitro)
dge_SRSi_insitu_against_invitro <- DGEList(counts=countmatrix_SRSi_insitu_against_invitro)

group_FRF_insitu_against_invitro <- as.factor(c("insitu","insitu","insitu","invitro","invitro","invitro"))
group_PBLUE_insitu_against_invitro <- as.factor(c("insitu","insitu","insitu","invitro","invitro","invitro"))
group_SRM_insitu_against_invitro <- as.factor(c("insitu","insitu","insitu","invitro","invitro","invitro"))
group_PMUD_insitu_against_invitro <- as.factor(c("insitu","insitu","insitu","invitro","invitro","invitro"))
group_SRB_insitu_against_invitro <- as.factor(c("insitu","insitu","insitu","invitro","invitro","invitro"))
group_PCB_insitu_against_invitro <- as.factor(c("insitu","insitu","insitu","invitro","invitro","invitro"))
group_SRSh_insitu_against_invitro <- as.factor(c("insitu","insitu","insitu","invitro","invitro","invitro"))
group_SRSi_insitu_against_invitro <- as.factor(c("insitu","insitu","insitu","invitro","invitro","invitro"))

keep_FRF_insitu_against_invitro <- filterByExpr(dge_FRF_insitu_against_invitro, group_FRF_insitu_against_invitro=group_FRF_insitu_against_invitro)
keep_PBLUE_insitu_against_invitro <- filterByExpr(dge_PBLUE_insitu_against_invitro, group_PBLUE_insitu_against_invitro=group_PBLUE_insitu_against_invitro)
keep_SRM_insitu_against_invitro <- filterByExpr(dge_SRM_insitu_against_invitro, group_SRM_insitu_against_invitro=group_SRM_insitu_against_invitro)
keep_PMUD_insitu_against_invitro <- filterByExpr(dge_PMUD_insitu_against_invitro, group_PMUD_insitu_against_invitro=group_PMUD_insitu_against_invitro)
keep_SRB_insitu_against_invitro <- filterByExpr(dge_SRB_insitu_against_invitro, group_SRB_insitu_against_invitro=group_SRB_insitu_against_invitro)
keep_PCB_insitu_against_invitro <- filterByExpr(dge_PCB_insitu_against_invitro, group_PCB_insitu_against_invitro=group_PCB_insitu_against_invitro)
keep_SRSh_insitu_against_invitro <- filterByExpr(dge_SRSh_insitu_against_invitro, group_SRSh_insitu_against_invitro=group_SRSh_insitu_against_invitro)
keep_SRSi_insitu_against_invitro <- filterByExpr(dge_SRSi_insitu_against_invitro, group_SRSi_insitu_against_invitro=group_SRSi_insitu_against_invitro)

dge_FRF_insitu_against_invitro <- dge_FRF_insitu_against_invitro[keep_FRF_insitu_against_invitro,,keep.lib.sizes=FALSE]
dge_PBLUE_insitu_against_invitro <- dge_PBLUE_insitu_against_invitro[keep_PBLUE_insitu_against_invitro,,keep.lib.sizes=FALSE]
dge_SRM_insitu_against_invitro <- dge_SRM_insitu_against_invitro[keep_SRM_insitu_against_invitro,,keep.lib.sizes=FALSE]
dge_PMUD_insitu_against_invitro <- dge_PMUD_insitu_against_invitro[keep_PMUD_insitu_against_invitro,,keep.lib.sizes=FALSE]
dge_SRB_insitu_against_invitro <- dge_SRB_insitu_against_invitro[keep_SRB_insitu_against_invitro,,keep.lib.sizes=FALSE]
dge_PCB_insitu_against_invitro <- dge_PCB_insitu_against_invitro[keep_PCB_insitu_against_invitro,,keep.lib.sizes=FALSE]
dge_SRSh_insitu_against_invitro <- dge_SRSh_insitu_against_invitro[keep_SRSh_insitu_against_invitro,,keep.lib.sizes=FALSE]
dge_SRSi_insitu_against_invitro <- dge_SRSi_insitu_against_invitro[keep_SRSi_insitu_against_invitro,,keep.lib.sizes=FALSE]

#checking the results of filtering
dim(dge_FRF_insitu_against_invitro)
dim(dge_PBLUE_insitu_against_invitro)
dim(dge_SRM_insitu_against_invitro)
dim(dge_PMUD_insitu_against_invitro)
dim(dge_SRB_insitu_against_invitro)
dim(dge_PCB_insitu_against_invitro)
dim(dge_SRSh_insitu_against_invitro)
dim(dge_SRSi_insitu_against_invitro)

#Use voom to determine weights for each gene to be passed to limma, voom normalises the data
v_FRF_insitu_against_invitro <- voom(dge_FRF_insitu_against_invitro, design_insitu_against_invitro, plot=TRUE, normalize="quantile")
v_PBLUE_insitu_against_invitro <- voom(dge_PBLUE_insitu_against_invitro, design_insitu_against_invitro, plot = TRUE, normalize = "quantile")
v_SRM_insitu_against_invitro <- voom(dge_SRM_insitu_against_invitro, design_insitu_against_invitro, plot = TRUE, normalize = "quantile")
v_PMUD_insitu_against_invitro <- voom(dge_PMUD_insitu_against_invitro, design_insitu_against_invitro, plot = TRUE, normalize = "quantile")
v_SRB_insitu_against_invitro <- voom(dge_SRB_insitu_against_invitro, design_insitu_against_invitro, plot = TRUE, normalize = "quantile")
v_PCB_insitu_against_invitro <- voom(dge_PCB_insitu_against_invitro, design_insitu_against_invitro, plot = TRUE, normalize = "quantile")
v_SRSh_insitu_against_invitro <- voom(dge_SRSh_insitu_against_invitro, design_insitu_against_invitro, plot = TRUE, normalize = "quantile")
v_SRSi_insitu_against_invitro <- voom(dge_SRSi_insitu_against_invitro, design_insitu_against_invitro, plot = TRUE, normalize = "quantile")

#fit linear model for each gene
fit_FRF_insitu_against_invitro <- lmFit(v_FRF_insitu_against_invitro,design_insitu_against_invitro)
fit_PBLUE_insitu_against_invitro <- lmFit(v_PBLUE_insitu_against_invitro, design_insitu_against_invitro)
fit_SRM_insitu_against_invitro <- lmFit(v_SRM_insitu_against_invitro, design_insitu_against_invitro)
fit_PMUD_insitu_against_invitro <- lmFit(v_PMUD_insitu_against_invitro, design_insitu_against_invitro)
fit_SRB_insitu_against_invitro <- lmFit(v_SRB_insitu_against_invitro, design_insitu_against_invitro)
fit_PCB_insitu_against_invitro <- lmFit(v_PCB_insitu_against_invitro, design_insitu_against_invitro)
fit_SRSh_insitu_against_invitro <- lmFit(v_SRSh_insitu_against_invitro, design_insitu_against_invitro)
fit_SRSi_insitu_against_invitro <- lmFit(v_SRSi_insitu_against_invitro, design_insitu_against_invitro)

#define the comparisons you want to examine
contrast.matrix_FRF_insitu_against_invitro <- makeContrasts(insitu - invitro, levels = design_insitu_against_invitro)
contrast.matrix_PBLUE_insitu_against_invitro <- makeContrasts(insitu - invitro, levels = design_insitu_against_invitro)
contrast.matrix_SRM_insitu_against_invitro <- makeContrasts(insitu - invitro, levels = design_insitu_against_invitro)
contrast.matrix_PMUD_insitu_against_invitro <- makeContrasts(insitu - invitro, levels = design_insitu_against_invitro)
contrast.matrix_SRB_insitu_against_invitro <- makeContrasts(insitu - invitro, levels = design_insitu_against_invitro)
contrast.matrix_PCB_insitu_against_invitro <- makeContrasts(insitu - invitro, levels = design_insitu_against_invitro)
contrast.matrix_SRSh_insitu_against_invitro <- makeContrasts(insitu - invitro, levels = design_insitu_against_invitro)
contrast.matrix_SRSi_insitu_against_invitro <- makeContrasts(insitu - invitro, levels = design_insitu_against_invitro)

#estimate contrasts for each gene
fit2_FRF_insitu_against_invitro <- contrasts.fit(fit_FRF_insitu_against_invitro, contrast.matrix_FRF_insitu_against_invitro)
fit2_PBLUE_insitu_against_invitro <- contrasts.fit(fit_PBLUE_insitu_against_invitro, contrast.matrix_PBLUE_insitu_against_invitro)
fit2_SRM_insitu_against_invitro <- contrasts.fit(fit_SRM_insitu_against_invitro, contrast.matrix_SRM_insitu_against_invitro)
fit2_PMUD_insitu_against_invitro <- contrasts.fit(fit_PMUD_insitu_against_invitro, contrast.matrix_PMUD_insitu_against_invitro)
fit2_SRB_insitu_against_invitro <- contrasts.fit(fit_SRB_insitu_against_invitro, contrast.matrix_SRB_insitu_against_invitro)
fit2_PCB_insitu_against_invitro <- contrasts.fit(fit_PCB_insitu_against_invitro, contrast.matrix_PCB_insitu_against_invitro)
fit2_SRSh_insitu_against_invitro <- contrasts.fit(fit_SRSh_insitu_against_invitro, contrast.matrix_SRSh_insitu_against_invitro)
fit2_SRSi_insitu_against_invitro <- contrasts.fit(fit_SRSi_insitu_against_invitro, contrast.matrix_SRSi_insitu_against_invitro)

#Smoothing of standard error
fit2_FRF_insitu_against_invitro <- eBayes(fit2_FRF_insitu_against_invitro)
fit2_PBLUE_insitu_against_invitro <- eBayes(fit2_PBLUE_insitu_against_invitro)
fit2_SRM_insitu_against_invitro <- eBayes(fit2_SRM_insitu_against_invitro)
fit2_PMUD_insitu_against_invitro <- eBayes(fit2_PMUD_insitu_against_invitro)
fit2_SRB_insitu_against_invitro <- eBayes(fit2_SRB_insitu_against_invitro)
fit2_PCB_insitu_against_invitro <- eBayes(fit2_PCB_insitu_against_invitro)
fit2_SRSh_insitu_against_invitro <- eBayes(fit2_SRSh_insitu_against_invitro)
fit2_SRSi_insitu_against_invitro <- eBayes(fit2_SRSi_insitu_against_invitro)

#Generate results tables
FRF_insitu_against_invitro <- topTable(fit2_FRF_insitu_against_invitro,coef=1, number=Inf)
PBLUE_insitu_against_invitro <- topTable(fit2_PBLUE_insitu_against_invitro, coef = 1, number = Inf)
SRM_insitu_against_invitro <- topTable(fit2_SRM_insitu_against_invitro, coef = 1, number = Inf)
PMUD_insitu_against_invitro <- topTable(fit2_PMUD_insitu_against_invitro, coef = 1, number = Inf)
SRB_insitu_against_invitro <- topTable(fit2_SRB_insitu_against_invitro, coef = 1, number = Inf)
PCB_insitu_against_invitro <- topTable(fit2_PCB_insitu_against_invitro, coef = 1, number = Inf)
SRSh_insitu_against_invitro <- topTable(fit2_SRSh_insitu_against_invitro, coef = 1, number = Inf)
SRSi_insitu_against_invitro <- topTable(fit2_SRSi_insitu_against_invitro, coef = 1, number = Inf)

#Export transcriptome profiles
write.csv(FRF_insitu_against_invitro, file = "./pst_removed_limma_outputs_insitu_vs_invitro/FRF_insitu_against_invitro.csv")
write.csv(PBLUE_insitu_against_invitro, file = "./pst_removed_limma_outputs_insitu_vs_invitro/PBLUE_insitu_against_invitro.csv")
write.csv(SRM_insitu_against_invitro, file = "./pst_removed_limma_outputs_insitu_vs_invitro/SRM_insitu_against_invitro.csv")
write.csv(PMUD_insitu_against_invitro, file = "./pst_removed_limma_outputs_insitu_vs_invitro/PMUD_insitu_against_invitro.csv")
write.csv(SRB_insitu_against_invitro, file = "./pst_removed_limma_outputs_insitu_vs_invitro/SRB_insitu_against_invitro.csv")
write.csv(PCB_insitu_against_invitro, file = "./pst_removed_limma_outputs_insitu_vs_invitro/PCB_insitu_against_invitro.csv")
write.csv(SRSh_insitu_against_invitro, file = "./pst_removed_limma_outputs_insitu_vs_invitro/SRSh_insitu_against_invitro.csv")
write.csv(SRSi_insitu_against_invitro, file = "./pst_removed_limma_outputs_insitu_vs_invitro/SRSi_insitu_against_invitro.csv")

