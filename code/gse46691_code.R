stringsAsFactors = FALSE

library(readr)
library(tidyverse)
library(xCell)

# import data
# data downloaded from : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46691
gse46691 <- read_delim("/Volumes/Lab_Gerke/dataDump/GSE46691_quantile_normalized.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

# annotation files downloaded from http://www.affymetrix.com/support/technical/byproduct.affx?product=huexon-st
anno <- read_csv("/Volumes/Lab_Gerke/dataDump/HuEx-1_0-st-v2.na36.hg19.probeset.csv", 
                 comment = "#")
anno <- anno[anno$probeset_id %in% gse46691$ID_REF, ]

# subset to those with known gene assignment 
anno <- anno[anno$gene_assignment!="---",]
gse <- gse46691[gse46691$ID_REF %in% anno$probeset_id,]

# pull last gene name
anno$gene <- gsub(".*// ", "", anno$gene_assignment)

# use gene names
gse <- merge(gse,anno[,c("probeset_id","gene")], by.x="ID_REF", by.y="probeset_id")

# group by gene names and take median for each person to use as single exp value
gse_single_gene <- gse[,-1] %>%
  group_by(gene) %>%
  summarise_all(median)

# grab names for matrix later
gene_names <- gse_single_gene$gene
sample_names <- colnames(gse_single_gene)[2:ncol(gse_single_gene)]

# transpose and numeric - might be able to skip since need to transpose again
exp <- as.data.frame(t(gse_single_gene))
colnames(exp) <- gene_names
exp <- exp[-1,]
exp[] <- lapply(exp, function(x) as.numeric(as.character(x)))
exp <- as.data.frame(cbind(sample_id = sample_names,exp))

x <- data.matrix(exp[,2:ncol(exp)])
xx <- t(x)

# run xCell
xy <- xCellAnalysis(xx)

mayo_xcell <- t(xy)
mayo_xcell <- as.data.frame(mayo_xcell)
mayo_xcell <- as.data.frame(cbind(sample_id = row.names(mayo_xcell),
                                  mayo_xcell))

write.table(mayo_xcell, file="/Volumes/Lab_Gerke/dataDump/GSE46691_xcell.txt",
            row.names = FALSE, quote=FALSE,sep="\t")

#########################################################################################
# group by gene names and take highest for each person to use as single exp value
gse_single_gene_max <- gse[,-1] %>%
  group_by(gene) %>%
  summarise_all(max)

# grab names for matrix later
gene_names_max <- gse_single_gene_max$gene
sample_names_max <- colnames(gse_single_gene_max)[2:ncol(gse_single_gene_max)]

# transpose and numeric - might be able to skip since need to transpose again
exp_max <- as.data.frame(t(gse_single_gene_max))
colnames(exp_max) <- gene_names_max
exp_max <- exp_max[-1,]
exp_max[] <- lapply(exp_max, function(x) as.numeric(as.character(x)))
exp_max <- as.data.frame(cbind(sample_id = sample_names_max,exp_max))

x_max <- data.matrix(exp_max[,2:ncol(exp_max)])
xx_max <- t(x_max)

# run xCell
xy_max <- xCellAnalysis(xx_max)

mayo_xcell_max <- t(xy_max)
mayo_xcell_max <- as.data.frame(mayo_xcell_max)
mayo_xcell_max <- as.data.frame(cbind(sample_id = row.names(mayo_xcell_max),
                                  mayo_xcell_max))

write.table(mayo_xcell_max, file="/Volumes/Lab_Gerke/dataDump/GSE46691_xcell_max.txt",
            row.names = FALSE, quote=FALSE,sep="\t")
