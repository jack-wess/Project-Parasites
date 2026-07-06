setwd("C:/Users/u6047942/JWPHD/project_parasites/bioinformatics/R-studio_transcriptomics")



library(DOSE)
library(tidyverse)
library(topGO)
library(lintr)
library(lattice)
library(dplyr)
library(Rgraphviz)

# Loading Data
# create GO background mapping
# using dataframe "XXXX6.csv" with two columns, gene name and GO terms http://127.0.0.1:45391/graphics/0d893d81-4f7a-45c6-b176-0d6fad76b074.png
FRF_gene_to_go_mapping <- readMappings(file = "./GO_terms/FRF_GO_terms_comma_separated.csv", sep = ";", IDsep = ",")
PBLUE_gene_to_go_mapping <- readMappings(file = "./GO_terms/PBLUE_GO_terms_comma_separated.csv", sep = ";", IDsep = ",")
SRM_gene_to_go_mapping <- readMappings(file = "./GO_terms/SRM_GO_terms_comma_separated.csv", sep = ";", IDsep = ",")
PMUD_gene_to_go_mapping <- readMappings(file = "./GO_terms/PMUD_GO_terms_comma_separated.csv", sep = ";", IDsep = ",")
SRB_gene_to_go_mapping <- readMappings(file = "./GO_terms/SRB_GO_terms_comma_separated.csv", sep = ";", IDsep = ",")
PCB_gene_to_go_mapping <- readMappings(file = "./GO_terms/PCB_GO_terms_comma_separated.csv", sep = ";", IDsep = ",")
SRSh_gene_to_go_mapping <- readMappings(file = "./GO_terms/SRSh_GO_terms_comma_separated.csv", sep = ";", IDsep = ",")
SRSi_gene_to_go_mapping <- readMappings(file = "./GO_terms/SRSi_GO_terms_comma_separated.csv", sep = ";", IDsep = ",")

# get background gene list
# using dataframe "pst_underattackbypenx_filtered_vs_day12.csv" with two columns, gene name and LFC value

go_FRF_insitu_against_invitro <-  read.csv("./pst_removed_limma_outputs_insitu_vs_invitro/FRF_insitu_against_invitro.csv", header = TRUE)
go_PBLUE_insitu_against_invitro <- read.csv("./pst_removed_limma_outputs_insitu_vs_invitro/PBLUE_insitu_against_invitro.csv", header = TRUE)
go_SRM_insitu_against_invitro <- read.csv("./pst_removed_limma_outputs_insitu_vs_invitro/SRM_insitu_against_invitro.csv", header = TRUE)
go_PMUD_insitu_against_invitro <- read.csv("./pst_removed_limma_outputs_insitu_vs_invitro/PMUD_insitu_against_invitro.csv", header = TRUE)
go_SRB_insitu_against_invitro <- read.csv("./pst_removed_limma_outputs_insitu_vs_invitro/SRB_insitu_against_invitro.csv", header = TRUE)
go_PCB_insitu_against_invitro <- read.csv("./pst_removed_limma_outputs_insitu_vs_invitro/PCB_insitu_against_invitro.csv", header = TRUE)
go_SRSh_insitu_against_invitro <- read.csv("./pst_removed_limma_outputs_insitu_vs_invitro/SRSh_insitu_against_invitro.csv", header = TRUE)
go_SRSi_insitu_against_invitro <- read.csv("./pst_removed_limma_outputs_insitu_vs_invitro/SRSi_insitu_against_invitro.csv", header = TRUE)

geneNames_FRF_insitu_against_invitro <- go_FRF_insitu_against_invitro$Column1
geneNames_PBLUE_insitu_against_invitro <- go_PBLUE_insitu_against_invitro$Column1
geneNames_SRM_insitu_against_invitro <- go_SRM_insitu_against_invitro$Column1
geneNames_PMUD_insitu_against_invitro <- go_PMUD_insitu_against_invitro$Column1
geneNames_SRB_insitu_against_invitro <- go_SRB_insitu_against_invitro$Column1
geneNames_PCB_insitu_against_invitro <- go_PCB_insitu_against_invitro$Column1
geneNames_SRSh_insitu_against_invitro <- go_SRSh_insitu_against_invitro$Column1
geneNames_SRSi_insitu_against_invitro <- go_SRSi_insitu_against_invitro$Column1

# select genes of interest
# Basically this is subsetting down genes based on whether they are significantly upregulated, maybe go higher than 1 LFC but lets see
myInterestingGenes_FRF_insitu_against_invitro <- subset(go_FRF_insitu_against_invitro, adj.P.Val < 0.05 & logFC > 1)$Column1
myInterestingGenes_PBLUE_insitu_against_invitro <- subset(go_PBLUE_insitu_against_invitro, adj.P.Val < 0.05 & logFC > 1)$Column1
myInterestingGenes_SRM_insitu_against_invitro <- subset(go_SRM_insitu_against_invitro, adj.P.Val < 0.05 & logFC > 1)$Column1
myInterestingGenes_PMUD_insitu_against_invitro <- subset(go_PMUD_insitu_against_invitro, adj.P.Val < 0.05 & logFC > 1)$Column1
myInterestingGenes_SRB_insitu_against_invitro <- subset(go_SRB_insitu_against_invitro, adj.P.Val < 0.05 & logFC > 1)$Column1
myInterestingGenes_PCB_insitu_against_invitro <- subset(go_PCB_insitu_against_invitro, adj.P.Val < 0.05 & logFC > 1)$Column1
myInterestingGenes_SRSh_insitu_against_invitro <- subset(go_SRSh_insitu_against_invitro, adj.P.Val < 0.05 & logFC > 1)$Column1
myInterestingGenes_SRSi_insitu_against_invitro <- subset(go_SRSi_insitu_against_invitro, adj.P.Val < 0.05 & logFC > 1)$Column1

# create selection and gene list object
geneList_FRF_insitu_against_invitro <- factor(as.integer(geneNames_FRF_insitu_against_invitro %in% myInterestingGenes_FRF_insitu_against_invitro)) 
names(geneList_FRF_insitu_against_invitro) <- geneNames_FRF_insitu_against_invitro
str(geneList_FRF_insitu_against_invitro)

geneList_PBLUE_insitu_against_invitro <- factor(as.integer(geneNames_PBLUE_insitu_against_invitro %in% myInterestingGenes_PBLUE_insitu_against_invitro))
names(geneList_PBLUE_insitu_against_invitro) <- geneNames_PBLUE_insitu_against_invitro
str(geneList_PBLUE_insitu_against_invitro)

geneList_SRM_insitu_against_invitro <- factor(as.integer(geneNames_SRM_insitu_against_invitro %in% myInterestingGenes_SRM_insitu_against_invitro))
names(geneList_SRM_insitu_against_invitro) <- geneNames_SRM_insitu_against_invitro
str(geneList_SRM_insitu_against_invitro)

geneList_PMUD_insitu_against_invitro <- factor(as.integer(geneNames_PMUD_insitu_against_invitro %in% myInterestingGenes_PMUD_insitu_against_invitro))
names(geneList_PMUD_insitu_against_invitro) <- geneNames_PMUD_insitu_against_invitro
str(geneList_PMUD_insitu_against_invitro)

geneList_SRB_insitu_against_invitro <- factor(as.integer(geneNames_SRB_insitu_against_invitro %in% myInterestingGenes_SRB_insitu_against_invitro))
names(geneList_SRB_insitu_against_invitro) <- geneNames_SRB_insitu_against_invitro
str(geneList_SRB_insitu_against_invitro)

geneList_PCB_insitu_against_invitro <- factor(as.integer(geneNames_PCB_insitu_against_invitro %in% myInterestingGenes_PCB_insitu_against_invitro))
names(geneList_PCB_insitu_against_invitro) <- geneNames_PCB_insitu_against_invitro
str(geneList_PCB_insitu_against_invitro)

geneList_SRSh_insitu_against_invitro <- factor(as.integer(geneNames_SRSh_insitu_against_invitro %in% myInterestingGenes_SRSh_insitu_against_invitro))
names(geneList_SRSh_insitu_against_invitro) <- geneNames_SRSh_insitu_against_invitro
str(geneList_SRSh_insitu_against_invitro)

geneList_SRSi_insitu_against_invitro <- factor(as.integer(geneNames_SRSi_insitu_against_invitro %in% myInterestingGenes_SRSi_insitu_against_invitro))
names(geneList_SRSi_insitu_against_invitro) <- geneNames_SRSi_insitu_against_invitro
str(geneList_SRSi_insitu_against_invitro)

# create topgo object
my_go_data_FRF_insitu_against_invitro <- new("topGOdata",
                                             ontology    = "BP",
                                             allGenes    = geneList_FRF_insitu_against_invitro,
                                             annot       = annFUN.gene2GO,
                                             gene2GO     = FRF_gene_to_go_mapping,
                                             nodeSize    = 5) # nodesize is min number of genes in GO_term

my_go_data_PBLUE_insitu_against_invitro <- new("topGOdata",
                                               ontology    = "BP",
                                               allGenes    = geneList_PBLUE_insitu_against_invitro,
                                               annot       = annFUN.gene2GO,
                                               gene2GO     = PBLUE_gene_to_go_mapping,
                                               nodeSize    = 5)

my_go_data_SRM_insitu_against_invitro <- new("topGOdata",
                                             ontology    = "BP",
                                             allGenes    = geneList_SRM_insitu_against_invitro,
                                             annot       = annFUN.gene2GO,
                                             gene2GO     = SRM_gene_to_go_mapping,
                                             nodeSize    = 5)

my_go_data_PMUD_insitu_against_invitro <- new("topGOdata",
                                              ontology    = "BP",
                                              allGenes    = geneList_PMUD_insitu_against_invitro,
                                              annot       = annFUN.gene2GO,
                                              gene2GO     = PMUD_gene_to_go_mapping,
                                              nodeSize    = 5)

my_go_data_SRB_insitu_against_invitro <- new("topGOdata",
                                             ontology    = "BP",
                                             allGenes    = geneList_SRB_insitu_against_invitro,
                                             annot       = annFUN.gene2GO,
                                             gene2GO     = SRB_gene_to_go_mapping,
                                             nodeSize    = 5)

my_go_data_PCB_insitu_against_invitro <- new("topGOdata",
                                             ontology    = "BP",
                                             allGenes    = geneList_PCB_insitu_against_invitro,
                                             annot       = annFUN.gene2GO,
                                             gene2GO     = PCB_gene_to_go_mapping,
                                             nodeSize    = 5)

my_go_data_SRSh_insitu_against_invitro <- new("topGOdata",
                                              ontology    = "BP",
                                              allGenes    = geneList_SRSh_insitu_against_invitro,
                                              annot       = annFUN.gene2GO,
                                              gene2GO     = SRSh_gene_to_go_mapping,
                                              nodeSize    = 5)

my_go_data_SRSi_insitu_against_invitro <- new("topGOdata",
                                              ontology    = "BP",
                                              allGenes    = geneList_SRSi_insitu_against_invitro,
                                              annot       = annFUN.gene2GO,
                                              gene2GO     = SRSi_gene_to_go_mapping,
                                              nodeSize    = 5)

## Enrichment Tests
# use fisher for gene set number enrichment (not using gene values)
result_weight_fisher_FRF_insitu_against_invitro <- runTest(object = my_go_data_FRF_insitu_against_invitro, algorithm = "weight01", statistic = "fisher") # result object
result_weight_fisher_PBLUE_insitu_against_invitro <- runTest(object = my_go_data_PBLUE_insitu_against_invitro, algorithm = "weight01", statistic = "fisher")
result_weight_fisher_SRM_insitu_against_invitro <- runTest(object = my_go_data_SRM_insitu_against_invitro, algorithm = "weight01", statistic = "fisher")
result_weight_fisher_PMUD_insitu_against_invitro <- runTest(object = my_go_data_PMUD_insitu_against_invitro, algorithm = "weight01", statistic = "fisher")
result_weight_fisher_SRB_insitu_against_invitro <- runTest(object = my_go_data_SRB_insitu_against_invitro, algorithm = "weight01", statistic = "fisher")
result_weight_fisher_PCB_insitu_against_invitro <- runTest(object = my_go_data_PCB_insitu_against_invitro, algorithm = "weight01", statistic = "fisher")
result_weight_fisher_SRSh_insitu_against_invitro <- runTest(object = my_go_data_SRSh_insitu_against_invitro, algorithm = "weight01", statistic = "fisher")
result_weight_fisher_SRSi_insitu_against_invitro <- runTest(object = my_go_data_SRSi_insitu_against_invitro, algorithm = "weight01", statistic = "fisher")

# summarise results
result_weight_output_FRF_insitu_against_invitro <- GenTable(object = my_go_data_FRF_insitu_against_invitro, 
                                                            weight_fisher_FRF_insitu_against_invitro = result_weight_fisher_FRF_insitu_against_invitro,
                                                            orderBy   = "weight_fisher_FRF_insitu_against_invitro",
                                                            topNodes  = length(score(result_weight_fisher_FRF_insitu_against_invitro)))

result_weight_output_PBLUE_insitu_against_invitro <- GenTable(object = my_go_data_PBLUE_insitu_against_invitro, 
                                                              weight_fisher_PBLUE_insitu_against_invitro = result_weight_fisher_PBLUE_insitu_against_invitro,
                                                              orderBy   = "weight_fisher_PBLUE_insitu_against_invitro",
                                                              topNodes  = length(score(result_weight_fisher_PBLUE_insitu_against_invitro)))

result_weight_output_SRM_insitu_against_invitro <- GenTable(object = my_go_data_SRM_insitu_against_invitro, 
                                                            weight_fisher_SRM_insitu_against_invitro = result_weight_fisher_SRM_insitu_against_invitro,
                                                            orderBy   = "weight_fisher_SRM_insitu_against_invitro",
                                                            topNodes  = length(score(result_weight_fisher_SRM_insitu_against_invitro)))

result_weight_output_PMUD_insitu_against_invitro <- GenTable(object = my_go_data_PMUD_insitu_against_invitro, 
                                                             weight_fisher_PMUD_insitu_against_invitro = result_weight_fisher_PMUD_insitu_against_invitro,
                                                             orderBy   = "weight_fisher_PMUD_insitu_against_invitro",
                                                             topNodes  = length(score(result_weight_fisher_PMUD_insitu_against_invitro)))

result_weight_output_SRB_insitu_against_invitro <- GenTable(object = my_go_data_SRB_insitu_against_invitro, 
                                                            weight_fisher_SRB_insitu_against_invitro = result_weight_fisher_SRB_insitu_against_invitro,
                                                            orderBy   = "weight_fisher_SRB_insitu_against_invitro",
                                                            topNodes  = length(score(result_weight_fisher_SRB_insitu_against_invitro)))

result_weight_output_PCB_insitu_against_invitro <- GenTable(object = my_go_data_PCB_insitu_against_invitro, 
                                                            weight_fisher_PCB_insitu_against_invitro = result_weight_fisher_PCB_insitu_against_invitro,
                                                            orderBy   = "weight_fisher_PCB_insitu_against_invitro",
                                                            topNodes  = length(score(result_weight_fisher_PCB_insitu_against_invitro)))

result_weight_output_SRSh_insitu_against_invitro <- GenTable(object = my_go_data_SRSh_insitu_against_invitro, 
                                                             weight_fisher_SRSh_insitu_against_invitro = result_weight_fisher_SRSh_insitu_against_invitro,
                                                             orderBy   = "weight_fisher_SRSh_insitu_against_invitro",
                                                             topNodes  = length(score(result_weight_fisher_SRSh_insitu_against_invitro)))

result_weight_output_SRSi_insitu_against_invitro <- GenTable(object = my_go_data_SRSi_insitu_against_invitro, 
                                                             weight_fisher_SRSi_insitu_against_invitro = result_weight_fisher_SRSi_insitu_against_invitro,
                                                             orderBy   = "weight_fisher_SRSi_insitu_against_invitro",
                                                             topNodes  = length(score(result_weight_fisher_SRSi_insitu_against_invitro)))

# multiple testing correction (may not be necessary)
result_weight_output_FRF_insitu_against_invitro$weight_fisher_adjusted_FRF_insitu_against_invitro <- p.adjust(p = result_weight_output_FRF_insitu_against_invitro$weight_fisher_FRF_insitu_against_invitro, method = c("BH")) 
result_weight_output_PBLUE_insitu_against_invitro$weight_fisher_adjusted_PBLUE_insitu_against_invitro <- p.adjust(p = result_weight_output_PBLUE_insitu_against_invitro$weight_fisher_PBLUE_insitu_against_invitro, method = "BH")
result_weight_output_SRM_insitu_against_invitro$weight_fisher_adjusted_SRM_insitu_against_invitro <- p.adjust(p = result_weight_output_SRM_insitu_against_invitro$weight_fisher_SRM_insitu_against_invitro, method = "BH")
result_weight_output_PMUD_insitu_against_invitro$weight_fisher_adjusted_PMUD_insitu_against_invitro <- p.adjust(p = result_weight_output_PMUD_insitu_against_invitro$weight_fisher_PMUD_insitu_against_invitro, method = "BH")
result_weight_output_SRB_insitu_against_invitro$weight_fisher_adjusted_SRB_insitu_against_invitro <- p.adjust(p = result_weight_output_SRB_insitu_against_invitro$weight_fisher_SRB_insitu_against_invitro, method = "BH")
result_weight_output_PCB_insitu_against_invitro$weight_fisher_adjusted_PCB_insitu_against_invitro <- p.adjust(p = result_weight_output_PCB_insitu_against_invitro$weight_fisher_PCB_insitu_against_invitro, method = "BH")
result_weight_output_SRSh_insitu_against_invitro$weight_fisher_adjusted_SRSh_insitu_against_invitro <- p.adjust(p = result_weight_output_SRSh_insitu_against_invitro$weight_fisher_SRSh_insitu_against_invitro, method = "BH")
result_weight_output_SRSi_insitu_against_invitro$weight_fisher_adjusted_SRSi_insitu_against_invitro <- p.adjust(p = result_weight_output_SRSi_insitu_against_invitro$weight_fisher_SRSi_insitu_against_invitro, method = "BH")

# P-value distribution from Fisher test
ggplot(result_weight_output_FRF_insitu_against_invitro, aes(x = as.numeric(weight_fisher_FRF_insitu_against_invitro))) + geom_histogram()
ggplot(result_weight_output_PBLUE_insitu_against_invitro, aes(x = as.numeric(weight_fisher_PBLUE_insitu_against_invitro))) + geom_histogram()
ggplot(result_weight_output_SRM_insitu_against_invitro, aes(x = as.numeric(weight_fisher_SRM_insitu_against_invitro))) + geom_histogram()
ggplot(result_weight_output_PMUD_insitu_against_invitro, aes(x = as.numeric(weight_fisher_PMUD_insitu_against_invitro))) + geom_histogram()
ggplot(result_weight_output_SRB_insitu_against_invitro, aes(x = as.numeric(weight_fisher_SRB_insitu_against_invitro))) + geom_histogram()
ggplot(result_weight_output_PCB_insitu_against_invitro, aes(x = as.numeric(weight_fisher_PCB_insitu_against_invitro))) + geom_histogram()
ggplot(result_weight_output_SRSh_insitu_against_invitro, aes(x = as.numeric(weight_fisher_SRSh_insitu_against_invitro))) + geom_histogram()
ggplot(result_weight_output_SRSi_insitu_against_invitro, aes(x = as.numeric(weight_fisher_SRSi_insitu_against_invitro))) + geom_histogram()

### Create a GO plot 
prefixes <- c("FRF", "PBLUE", "SRM", "PMUD", "SRB", "PCB", "SRSh", "SRSi")

for (prefix in prefixes) {
  output_dir <- "go_plots"
  output_file <- file.path(output_dir, paste0("topGO_", prefix, "_insitu_against_invitro.png"))
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  png(filename = output_file, width = 1200, height = 1000, res = 300)
  par(cex = 0.22)
  showSigOfNodes(
    get(paste0("my_go_data_", prefix, "_insitu_against_invitro")),
    score(get(paste0("result_weight_fisher_", prefix, "_insitu_against_invitro"))),
    firstSigNodes = 5,
    useInfo = 'def',
    .NO.CHAR = 50
  )
  dev.off()
  
  assign(
    paste0("GO_mapped_", prefix, "_insitu_against_invitro"),
    genesInTerm(get(paste0("my_go_data_", prefix, "_insitu_against_invitro")))
  )
}

## Plotting genes in sig GO terms
# get sig GO terms

sig_terms_FRF_insitu_against_invitro <- result_weight_output_FRF_insitu_against_invitro[as.numeric(result_weight_output_FRF_insitu_against_invitro$weight_fisher_FRF_insitu_against_invitro) < 0.05,]
sig_terms_PBLUE_insitu_against_invitro <- result_weight_output_PBLUE_insitu_against_invitro[as.numeric(result_weight_output_PBLUE_insitu_against_invitro$weight_fisher_PBLUE_insitu_against_invitro) < 0.05,]
sig_terms_SRM_insitu_against_invitro <- result_weight_output_SRM_insitu_against_invitro[as.numeric(result_weight_output_SRM_insitu_against_invitro$weight_fisher_SRM_insitu_against_invitro) < 0.05,]
sig_terms_PMUD_insitu_against_invitro <- result_weight_output_PMUD_insitu_against_invitro[as.numeric(result_weight_output_PMUD_insitu_against_invitro$weight_fisher_PMUD_insitu_against_invitro) < 0.05,]
sig_terms_SRB_insitu_against_invitro <- result_weight_output_SRB_insitu_against_invitro[as.numeric(result_weight_output_SRB_insitu_against_invitro$weight_fisher_SRB_insitu_against_invitro) < 0.05,]
sig_terms_PCB_insitu_against_invitro <- result_weight_output_PCB_insitu_against_invitro[as.numeric(result_weight_output_PCB_insitu_against_invitro$weight_fisher_PCB_insitu_against_invitro) < 0.05,]
sig_terms_SRSh_insitu_against_invitro <- result_weight_output_SRSh_insitu_against_invitro[as.numeric(result_weight_output_SRSh_insitu_against_invitro$weight_fisher_SRSh_insitu_against_invitro) < 0.05,]
sig_terms_SRSi_insitu_against_invitro <- result_weight_output_SRSi_insitu_against_invitro[as.numeric(result_weight_output_SRSi_insitu_against_invitro$weight_fisher_SRSi_insitu_against_invitro) < 0.05,]

# Code to generate a dotplot of enriched GO terms, with Gene Ratio (meaning number of significant genes within an individual GO term) as the x axis and GO term as the y axis. The dot size is relative to the total number of significant genes in the term and its colour
# reflects the significance of the term overall (Fisher test)
theme_update(text = element_text(size=18))
prefixes <- c("FRF","PBLUE", "SRM", "PMUD", "SRB", "PCB", "SRSh", "SRSi")
for (prefix in prefixes) {
  sig_terms <- get(paste0("sig_terms_", prefix, "_insitu_against_invitro"))
  
p <- ggplot(sig_terms, aes(x = Significant / Annotated, y = reorder(Term, Significant / Annotated))) +
    geom_point(aes(size = Significant, color = get(paste0("weight_fisher_adjusted_", prefix, "_insitu_against_invitro"))), binaxis = "y") +
    xlim(0, 1.0) +
    labs(
      y = "GO term",
      x = "Gene ratio",
      size = "# Significant Genes",
      colour = "Significance"
    ) +
    ggtitle(paste(prefix, "insitu vs invitro")) +
    theme(axis.text = element_text(size = 8))
  
  output_dir <- "go_plots"
  output_file <- file.path(output_dir, paste0(prefix, "_insitu_vs_invitro_jacks_ratio_plot.png"))
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  ggsave(filename = output_file, plot = p, width = 10, height = 8, dpi = 300)
}

# Required libraries
library(ggplot2)
library(patchwork)
library(dplyr)

theme_update(text = element_text(size = 18))

prefixes <- c("FRF", "PBLUE", "SRM", "PMUD", "SRB", "PCB", "SRSh", "SRSi")
go_plots <- list()

for (prefix in prefixes) {
  df <- get(paste0("sig_terms_", prefix, "_insitu_against_invitro"))
  
  # Find the correct column for adjusted p-value
  adj_col <- grep(paste0("weight_fisher_adjusted_", prefix, "_insitu_against_invitro"),
                  colnames(df), value = TRUE)
  
  if (length(adj_col) == 0) {
    warning(paste("No adjusted p-value column found for", prefix))
    next
  }
  
  # Filter to significance threshold < 0.05
  filtered_df <- df %>% filter(.data[[adj_col]] < 0.05)
  
  if (nrow(filtered_df) == 0) {
    warning(paste("No significant terms (adj.p < 0.05) for", prefix))
    next
  }
  
  p <- ggplot(filtered_df, aes(x = Significant / Annotated, y = reorder(Term, Significant / Annotated))) +
    geom_point(aes(size = Significant, color = .data[[adj_col]])) +
    xlim(0, 1.0) +
    labs(
      y = "GO term",
      x = "Gene ratio",
      size = "# Significant Genes",
      colour = "Significance"
    ) +
    ggtitle(paste(prefix, "insitu vs invitro")) +
    theme(axis.text = element_text(size = 8))
  
  go_plots[[prefix]] <- p
}

# Combine available plots into a 4x2 composite
if (length(go_plots) > 0) {
  composite_go_plot <- wrap_plots(go_plots, ncol = 2)
  
  ggsave("filtered_go_composite_dotplot.png",
         composite_go_plot,
         width = 20, height = 25, dpi = 300)
} else {
  message("No valid plots after filtering.")
}


# get genes in sig go terms
sig_genes_FRF_insitu_against_invitro <- data.frame()
for(go_term in sig_terms_FRF_insitu_against_invitro$GO.ID){
  sig_gene_FRF_insitu_against_invitro <- data.frame(gene = GO_mapped_FRF_insitu_against_invitro[go_term], GO.ID = go_term)
  names(sig_gene_FRF_insitu_against_invitro)[1] <- "Column1"
  sig_genes_FRF_insitu_against_invitro <- rbind(sig_genes_FRF_insitu_against_invitro, sig_gene_FRF_insitu_against_invitro)
}

sig_genes_PBLUE_insitu_against_invitro <- data.frame()
for (go_term in sig_terms_PBLUE_insitu_against_invitro$GO.ID) {
  sig_gene_df <- data.frame(gene = GO_mapped_PBLUE_insitu_against_invitro[[go_term]], GO.ID = go_term)
  names(sig_gene_df)[1] <- "Column1"
  sig_genes_PBLUE_insitu_against_invitro <- rbind(sig_genes_PBLUE_insitu_against_invitro, sig_gene_df)
}

sig_genes_SRM_insitu_against_invitro <- data.frame()
for (go_term in sig_terms_SRM_insitu_against_invitro$GO.ID) {
  sig_gene_df <- data.frame(gene = GO_mapped_SRM_insitu_against_invitro[[go_term]], GO.ID = go_term)
  names(sig_gene_df)[1] <- "Column1"
  sig_genes_SRM_insitu_against_invitro <- rbind(sig_genes_SRM_insitu_against_invitro, sig_gene_df)
}

sig_genes_PMUD_insitu_against_invitro <- data.frame()
for (go_term in sig_terms_PMUD_insitu_against_invitro$GO.ID) {
  sig_gene_df <- data.frame(gene = GO_mapped_PMUD_insitu_against_invitro[[go_term]], GO.ID = go_term)
  names(sig_gene_df)[1] <- "Column1"
  sig_genes_PMUD_insitu_against_invitro <- rbind(sig_genes_PMUD_insitu_against_invitro, sig_gene_df)
}

sig_genes_SRB_insitu_against_invitro <- data.frame()
for (go_term in sig_terms_SRB_insitu_against_invitro$GO.ID) {
  sig_gene_df <- data.frame(gene = GO_mapped_SRB_insitu_against_invitro[[go_term]], GO.ID = go_term)
  names(sig_gene_df)[1] <- "Column1"
  sig_genes_SRB_insitu_against_invitro <- rbind(sig_genes_SRB_insitu_against_invitro, sig_gene_df)
}

sig_genes_PCB_insitu_against_invitro <- data.frame()
for (go_term in sig_terms_PCB_insitu_against_invitro$GO.ID) {
  sig_gene_df <- data.frame(gene = GO_mapped_PCB_insitu_against_invitro[[go_term]], GO.ID = go_term)
  names(sig_gene_df)[1] <- "Column1"
  sig_genes_PCB_insitu_against_invitro <- rbind(sig_genes_PCB_insitu_against_invitro, sig_gene_df)
}

sig_genes_SRSh_insitu_against_invitro <- data.frame()
for (go_term in sig_terms_SRSh_insitu_against_invitro$GO.ID) {
  sig_gene_df <- data.frame(gene = GO_mapped_SRSh_insitu_against_invitro[[go_term]], GO.ID = go_term)
  names(sig_gene_df)[1] <- "Column1"
  sig_genes_SRSh_insitu_against_invitro <- rbind(sig_genes_SRSh_insitu_against_invitro, sig_gene_df)
}

sig_genes_SRSi_insitu_against_invitro <- data.frame()
for (go_term in sig_terms_SRSi_insitu_against_invitro$GO.ID) {
  sig_gene_df <- data.frame(gene = GO_mapped_SRSi_insitu_against_invitro[[go_term]], GO.ID = go_term)
  names(sig_gene_df)[1] <- "Column1"
  sig_genes_SRSi_insitu_against_invitro <- rbind(sig_genes_SRSi_insitu_against_invitro, sig_gene_df)
}

sig_terms_genes_FRF_insitu_against_invitro <- right_join(sig_terms_FRF_insitu_against_invitro, sig_genes_FRF_insitu_against_invitro, by = "GO.ID")
sig_terms_genes_PBLUE_insitu_against_invitro <- right_join(sig_terms_PBLUE_insitu_against_invitro, sig_genes_PBLUE_insitu_against_invitro, by = "GO.ID")
sig_terms_genes_SRM_insitu_against_invitro <- right_join(sig_terms_SRM_insitu_against_invitro, sig_genes_SRM_insitu_against_invitro, by = "GO.ID")
sig_terms_genes_PMUD_insitu_against_invitro <- right_join(sig_terms_PMUD_insitu_against_invitro, sig_genes_PMUD_insitu_against_invitro, by = "GO.ID")
sig_terms_genes_SRB_insitu_against_invitro <- right_join(sig_terms_SRB_insitu_against_invitro, sig_genes_SRB_insitu_against_invitro, by = "GO.ID")
sig_terms_genes_PCB_insitu_against_invitro <- right_join(sig_terms_PCB_insitu_against_invitro, sig_genes_PCB_insitu_against_invitro, by = "GO.ID")
sig_terms_genes_SRSh_insitu_against_invitro <- right_join(sig_terms_SRSh_insitu_against_invitro, sig_genes_SRSh_insitu_against_invitro, by = "GO.ID")
sig_terms_genes_SRSi_insitu_against_invitro <- right_join(sig_terms_SRSi_insitu_against_invitro, sig_genes_SRSi_insitu_against_invitro, by = "GO.ID")

# associate LFC with sig genes
sig_LFC_FRF_insitu_against_invitro <- right_join(go_FRF_insitu_against_invitro, sig_terms_genes_FRF_insitu_against_invitro, by = "Column1")
sig_LFC_PBLUE_insitu_against_invitro <- right_join(go_PBLUE_insitu_against_invitro, sig_terms_genes_PBLUE_insitu_against_invitro, by = "Column1")
sig_LFC_SRM_insitu_against_invitro <- right_join(go_SRM_insitu_against_invitro, sig_terms_genes_SRM_insitu_against_invitro, by = "Column1")
sig_LFC_PMUD_insitu_against_invitro <- right_join(go_PMUD_insitu_against_invitro, sig_terms_genes_PMUD_insitu_against_invitro, by = "Column1")
sig_LFC_SRB_insitu_against_invitro <- right_join(go_SRB_insitu_against_invitro, sig_terms_genes_SRB_insitu_against_invitro, by = "Column1")
sig_LFC_PCB_insitu_against_invitro <- right_join(go_PCB_insitu_against_invitro, sig_terms_genes_PCB_insitu_against_invitro, by = "Column1")
sig_LFC_SRSh_insitu_against_invitro <- right_join(go_SRSh_insitu_against_invitro, sig_terms_genes_SRSh_insitu_against_invitro, by = "Column1")
sig_LFC_SRSi_insitu_against_invitro <- right_join(go_SRSi_insitu_against_invitro, sig_terms_genes_SRSi_insitu_against_invitro, by = "Column1")

# order by most sig GO term
sig_LFC_FRF_insitu_against_invitro <- sig_LFC_FRF_insitu_against_invitro[order(as.numeric(sig_LFC_FRF_insitu_against_invitro$weight_fisher_FRF_insitu_against_invitro)),]
sig_LFC_FRF_insitu_against_invitro$Term <- factor(sig_LFC_FRF_insitu_against_invitro$Term , levels = rev(unique(sig_LFC_FRF_insitu_against_invitro$Term )))
sig_LFC_PBLUE_insitu_against_invitro <- sig_LFC_PBLUE_insitu_against_invitro[order(as.numeric(sig_LFC_PBLUE_insitu_against_invitro$weight_fisher_PBLUE_insitu_against_invitro)), ]
sig_LFC_PBLUE_insitu_against_invitro$Term <- factor(sig_LFC_PBLUE_insitu_against_invitro$Term, levels = rev(unique(sig_LFC_PBLUE_insitu_against_invitro$Term)))
sig_LFC_SRM_insitu_against_invitro <- sig_LFC_SRM_insitu_against_invitro[order(as.numeric(sig_LFC_SRM_insitu_against_invitro$weight_fisher_SRM_insitu_against_invitro)), ]
sig_LFC_SRM_insitu_against_invitro$Term <- factor(sig_LFC_SRM_insitu_against_invitro$Term, levels = rev(unique(sig_LFC_SRM_insitu_against_invitro$Term)))
sig_LFC_PMUD_insitu_against_invitro <- sig_LFC_PMUD_insitu_against_invitro[order(as.numeric(sig_LFC_PMUD_insitu_against_invitro$weight_fisher_PMUD_insitu_against_invitro)), ]
sig_LFC_PMUD_insitu_against_invitro$Term <- factor(sig_LFC_PMUD_insitu_against_invitro$Term, levels = rev(unique(sig_LFC_PMUD_insitu_against_invitro$Term)))
sig_LFC_SRB_insitu_against_invitro <- sig_LFC_SRB_insitu_against_invitro[order(as.numeric(sig_LFC_SRB_insitu_against_invitro$weight_fisher_SRB_insitu_against_invitro)), ]
sig_LFC_SRB_insitu_against_invitro$Term <- factor(sig_LFC_SRB_insitu_against_invitro$Term, levels = rev(unique(sig_LFC_SRB_insitu_against_invitro$Term)))
sig_LFC_PCB_insitu_against_invitro <- sig_LFC_PCB_insitu_against_invitro[order(as.numeric(sig_LFC_PCB_insitu_against_invitro$weight_fisher_PCB_insitu_against_invitro)), ]
sig_LFC_PCB_insitu_against_invitro$Term <- factor(sig_LFC_PCB_insitu_against_invitro$Term, levels = rev(unique(sig_LFC_PCB_insitu_against_invitro$Term)))
sig_LFC_SRSh_insitu_against_invitro <- sig_LFC_SRSh_insitu_against_invitro[order(as.numeric(sig_LFC_SRSh_insitu_against_invitro$weight_fisher_SRSh_insitu_against_invitro)), ]
sig_LFC_SRSh_insitu_against_invitro$Term <- factor(sig_LFC_SRSh_insitu_against_invitro$Term, levels = rev(unique(sig_LFC_SRSh_insitu_against_invitro$Term)))
sig_LFC_SRSi_insitu_against_invitro <- sig_LFC_SRSi_insitu_against_invitro[order(as.numeric(sig_LFC_SRSi_insitu_against_invitro$weight_fisher_SRSi_insitu_against_invitro)), ]
sig_LFC_SRSi_insitu_against_invitro$Term <- factor(sig_LFC_SRSi_insitu_against_invitro$Term, levels = rev(unique(sig_LFC_SRSi_insitu_against_invitro$Term)))

# plot LFC
# use "GO.ID" if GO "Term" is not unique

# Create the plot
library(ggplot2)

output_dir <- "go_plots"
if (!dir.exists(output_dir)) dir.create(output_dir)

# FRF
sig_LFC_FRF_insitu_against_invitro$Term <- factor(
  sig_LFC_FRF_insitu_against_invitro$Term,
  levels = rev(unique(sig_LFC_FRF_insitu_against_invitro$Term))
)
num_terms <- length(unique(sig_LFC_FRF_insitu_against_invitro$Term))
plot_height <- max(8, num_terms * 0.35)
p_FRF <- ggplot(sig_LFC_FRF_insitu_against_invitro, aes(
  x = logFC, y = Term,
  fill = as.numeric(weight_fisher_FRF_insitu_against_invitro))
) +
  geom_violin(width = 0.6, scale = "width") +
  theme_bw(base_size = 14) +
  geom_vline(xintercept = 0, colour = "grey") +
  scale_fill_gradient(low = "#2E8B57", high = "#AFEEEE") +
  labs(y = "GO term", x = "logFC") +
  theme(axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
ggsave(filename = file.path(output_dir, "FRF_insitu_against_invitro_violin_plot.png"),
       plot = p_FRF, width = 12, height = plot_height, dpi = 300)

# PBLUE
sig_LFC_PBLUE_insitu_against_invitro$Term <- factor(
  sig_LFC_PBLUE_insitu_against_invitro$Term,
  levels = rev(unique(sig_LFC_PBLUE_insitu_against_invitro$Term))
)
num_terms <- length(unique(sig_LFC_PBLUE_insitu_against_invitro$Term))
plot_height <- max(8, num_terms * 0.35)
p_PBLUE <- ggplot(sig_LFC_PBLUE_insitu_against_invitro, aes(
  x = logFC, y = Term,
  fill = as.numeric(weight_fisher_PBLUE_insitu_against_invitro))
) +
  geom_violin(width = 0.6, scale = "width") +
  theme_bw(base_size = 14) +
  geom_vline(xintercept = 0, colour = "grey") +
  scale_fill_gradient(low = "blue", high = "hotpink") +
  labs(y = "GO term", x = "logFC") +
  theme(axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
ggsave(filename = file.path(output_dir, "PBLUE_insitu_against_invitro_violin_plot.png"),
       plot = p_PBLUE, width = 12, height = plot_height, dpi = 300)

# SRM
sig_LFC_SRM_insitu_against_invitro$Term <- factor(
  sig_LFC_SRM_insitu_against_invitro$Term,
  levels = rev(unique(sig_LFC_SRM_insitu_against_invitro$Term))
)
num_terms <- length(unique(sig_LFC_SRM_insitu_against_invitro$Term))
plot_height <- max(8, num_terms * 0.35)
p_SRM <- ggplot(sig_LFC_SRM_insitu_against_invitro, aes(
  x = logFC, y = Term,
  fill = as.numeric(weight_fisher_SRM_insitu_against_invitro))
) +
  geom_violin(width = 0.6, scale = "width") +
  theme_bw(base_size = 14) +
  geom_vline(xintercept = 0, colour = "grey") +
  scale_fill_gradient(low = "red", high = "hotpink") +
  labs(y = "GO term", x = "logFC") +
  theme(axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
ggsave(filename = file.path(output_dir, "SRM_insitu_against_invitro_violin_plot.png"),
       plot = p_SRM, width = 12, height = plot_height, dpi = 300)

# PMUD
sig_LFC_PMUD_insitu_against_invitro$Term <- factor(
  sig_LFC_PMUD_insitu_against_invitro$Term,
  levels = rev(unique(sig_LFC_PMUD_insitu_against_invitro$Term))
)
num_terms <- length(unique(sig_LFC_PMUD_insitu_against_invitro$Term))
plot_height <- max(8, num_terms * 0.35)
p_PMUD <- ggplot(sig_LFC_PMUD_insitu_against_invitro, aes(
  x = logFC, y = Term,
  fill = as.numeric(weight_fisher_PMUD_insitu_against_invitro))
) +
  geom_violin(width = 0.6, scale = "width") +
  theme_bw(base_size = 14) +
  geom_vline(xintercept = 0, colour = "grey") +
  scale_fill_gradient(low = "#5A3E2B", high = "#C4A484") +
  labs(y = "GO term", x = "logFC") +
  theme(axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
ggsave(filename = file.path(output_dir, "PMUD_insitu_against_invitro_violin_plot.png"),
       plot = p_PMUD, width = 12, height = plot_height, dpi = 300)

# SRB
sig_LFC_SRB_insitu_against_invitro$Term <- factor(
  sig_LFC_SRB_insitu_against_invitro$Term,
  levels = rev(unique(sig_LFC_SRB_insitu_against_invitro$Term))
)
num_terms <- length(unique(sig_LFC_SRB_insitu_against_invitro$Term))
plot_height <- max(8, num_terms * 0.35)
p_SRB <- ggplot(sig_LFC_SRB_insitu_against_invitro, aes(
  x = logFC, y = Term,
  fill = as.numeric(weight_fisher_SRB_insitu_against_invitro))
) +
  geom_violin(width = 0.6, scale = "width") +
  theme_bw(base_size = 14) +
  geom_vline(xintercept = 0, colour = "grey") +
  scale_fill_gradient(low = "red", high = "lightblue") +
  labs(y = "GO term", x = "logFC") +
  theme(axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
ggsave(filename = file.path(output_dir, "SRB_insitu_against_invitro_violin_plot.png"),
       plot = p_SRB, width = 12, height = plot_height, dpi = 300)

# PCB
sig_LFC_PCB_insitu_against_invitro$Term <- factor(
  sig_LFC_PCB_insitu_against_invitro$Term,
  levels = rev(unique(sig_LFC_PCB_insitu_against_invitro$Term))
)
num_terms <- length(unique(sig_LFC_PCB_insitu_against_invitro$Term))
plot_height <- max(8, num_terms * 0.35)
p_PCB <- ggplot(sig_LFC_PCB_insitu_against_invitro, aes(
  x = logFC, y = Term,
  fill = as.numeric(weight_fisher_PCB_insitu_against_invitro))
) +
  geom_violin(width = 0.6, scale = "width") +
  theme_bw(base_size = 14) +
  geom_vline(xintercept = 0, colour = "grey") +
  scale_fill_gradient(low = "pink", high = "white") +
  labs(y = "GO term", x = "logFC") +
  theme(axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
ggsave(filename = file.path(output_dir, "PCB_insitu_against_invitro_violin_plot.png"),
       plot = p_PCB, width = 12, height = plot_height, dpi = 300)

# SRSh
sig_LFC_SRSh_insitu_against_invitro$Term <- factor(
  sig_LFC_SRSh_insitu_against_invitro$Term,
  levels = rev(unique(sig_LFC_SRSh_insitu_against_invitro$Term))
)
num_terms <- length(unique(sig_LFC_SRSh_insitu_against_invitro$Term))
plot_height <- max(8, num_terms * 0.35)
p_SRSh <- ggplot(sig_LFC_SRSh_insitu_against_invitro, aes(
  x = logFC, y = Term,
  fill = as.numeric(weight_fisher_SRSh_insitu_against_invitro))
) +
  geom_violin(width = 0.6, scale = "width") +
  theme_bw(base_size = 14) +
  geom_vline(xintercept = 0, colour = "grey") +
  scale_fill_gradient(low = "purple", high = "yellow") +
  labs(y = "GO term", x = "logFC") +
  theme(axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
ggsave(filename = file.path(output_dir, "SRSh_insitu_against_invitro_violin_plot.png"),
       plot = p_SRSh, width = 12, height = plot_height, dpi = 300)

# SRSi
sig_LFC_SRSi_insitu_against_invitro$Term <- factor(
  sig_LFC_SRSi_insitu_against_invitro$Term,
  levels = rev(unique(sig_LFC_SRSi_insitu_against_invitro$Term))
)
num_terms <- length(unique(sig_LFC_SRSi_insitu_against_invitro$Term))
plot_height <- max(8, num_terms * 0.35)
p_SRSi <- ggplot(sig_LFC_SRSi_insitu_against_invitro, aes(
  x = logFC, y = Term,
  fill = as.numeric(weight_fisher_SRSi_insitu_against_invitro))
) +
  geom_violin(width = 0.6, scale = "width") +
  theme_bw(base_size = 14) +
  geom_vline(xintercept = 0, colour = "grey") +
  scale_fill_gradient(low = "black", high = "orange") +
  labs(y = "GO term", x = "logFC") +
  theme(axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
ggsave(filename = file.path(output_dir, "SRSi_insitu_against_invitro_violin_plot.png"),
       plot = p_SRSi, width = 12, height = plot_height, dpi = 300)




























#### TO GENERATE A PCA FOR GO ENRICHMENT 

# Helper function to process each table
process_table <- function(df, go_col, pval_col, species_name) {
  df <- df[, c(go_col, pval_col)]
  colnames(df) <- c("GO.ID", "pval")
  df$pval <- as.numeric(df$pval)
  df$score <- -log10(df$pval)
  df$Species <- species_name
  return(df[, c("GO.ID", "Species", "score")])
}

# Process each species
FRF  <- process_table(sig_terms_FRF_insitu_against_invitro, "GO.ID", "weight_fisher_FRF_insitu_against_invitro", "FRF")
PBLUE <- process_table(sig_terms_PBLUE_insitu_against_invitro, "GO.ID", "weight_fisher_PBLUE_insitu_against_invitro", "PBLUE")
SRM <- process_table(sig_terms_SRM_insitu_against_invitro, "GO.ID", "weight_fisher_SRM_insitu_against_invitro", "SRM")
PMUD <- process_table(sig_terms_PMUD_insitu_against_invitro, "GO.ID", "weight_fisher_PMUD_insitu_against_invitro", "PMUD")
SRB <- process_table(sig_terms_SRB_insitu_against_invitro, "GO.ID", "weight_fisher_SRB_insitu_against_invitro", "SRB")
PCB <- process_table(sig_terms_PCB_insitu_against_invitro, "GO.ID", "weight_fisher_PCB_insitu_against_invitro", "PCB")
SRSh <- process_table(sig_terms_SRSh_insitu_against_invitro, "GO.ID", "weight_fisher_SRSh_insitu_against_invitro", "SRSh")
SRSi <- process_table(sig_terms_SRSi_insitu_against_invitro, "GO.ID", "weight_fisher_SRSi_insitu_against_invitro", "SRSi")

# Combine all
all_data <- rbind(FRF, PBLUE, SRM, PMUD, SRB, PCB, SRSh, SRSi)

# Pivot to GO x Species matrix
library(tidyr)
library(dplyr)

go_matrix <- all_data %>%
  pivot_wider(names_from = Species, values_from = score, values_fill = 0) %>%
  column_to_rownames("GO.ID")

# View matrix
head(go_matrix)

# Optional: Save for later use
write.csv(go_matrix, "GO_term_enrichment_matrix.csv")


library(ggplot2)
library(dplyr)

# Example inputs:
# pca: your PCA object from prcomp()
# scores_df: data.frame with columns PC1, PC2, Species
# loadings_df: data.frame with columns PC1, PC2, GO.ID


go_terms_desc <- bind_rows(
  sig_terms_FRF_insitu_against_invitro,
  sig_terms_PBLUE_insitu_against_invitro,
  sig_terms_SRM_insitu_against_invitro,
  sig_terms_PMUD_insitu_against_invitro,
  sig_terms_SRB_insitu_against_invitro,
  sig_terms_PCB_insitu_against_invitro,
  sig_terms_SRSh_insitu_against_invitro,
  sig_terms_SRSi_insitu_against_invitro
) %>%
  dplyr::select(GO.ID, Term) %>%
  dplyr::distinct()


# Custom color palette
color_palette <- c(
  FRF = "#2c7c71", PCB = "#d97c7c", PBLUE = "#4a66b2",
  PMUD = "#7f4f24", SRB = "red", SRM = "#cc4c99",
  SRSh = "#8b74a3", SRSi = "#c9a800"
)

# Calculate scaling factor for loadings to spread GO terms in PCA plot space
scale_factor_x <- diff(range(scores_df$PC1)) / diff(range(loadings_df$PC1))
scale_factor_y <- diff(range(scores_df$PC2)) / diff(range(loadings_df$PC2))
scale_factor <- min(scale_factor_x, scale_factor_y) * 0.8  # Adjust 0.8 as needed

# Scale GO term loadings
loadings_df <- loadings_df %>%
  mutate(
    PC1_scaled = PC1 * scale_factor,
    PC2_scaled = PC2 * scale_factor
  )

# Plot PCA scores + GO term labels
ggplot() +
  # Species points colored by Species
  geom_point(data = scores_df, aes(x = PC1, y = PC2, color = Species), size = 4) +
  
  # GO term labels at scaled loading positions (smaller text, no arrows)
  geom_text(data = loadings_df, aes(x = PC1_scaled, y = PC2_scaled, label = Term),
            size = 2.5, alpha = 0.6, color = "black", check_overlap = TRUE) +
  
  # Apply custom colors for species
  scale_color_manual(values = color_palette) +
  
  # Clean theme & axis labels with variance %
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA of GO Term Enrichment During Growth on Pst Urediniospores",
    x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "% variance)")
  )

top_pc1 <- pca$rotation %>%
  as.data.frame() %>%
  rownames_to_column("GO.ID") %>%
  arrange(desc(abs(PC1))) %>%
  left_join(go_terms_desc, by = "GO.ID") %>%
  head(10)

# View top contributors to PC2
top_pc2 <- pca$rotation %>%
  as.data.frame() %>%
  rownames_to_column("GO.ID") %>%
  arrange(desc(abs(PC2))) %>%
  left_join(go_terms_desc, by = "GO.ID") %>%
  head(10)


