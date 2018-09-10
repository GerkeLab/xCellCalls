.libPaths(c("~/Documents/Rpackagesv2",.libPaths()))
stringsAsFactors = FALSE

library(readr)
library(xCell)

# import data
# data imported from cbioportal : http://www.cbioportal.org/study?id=prad_tcga_pub#summary
exp <- read_delim("/Volumes/Lab_Gerke/TCGA/prad_tcga_pub/prad_tcga_pub/data_RNA_Seq_v2_expression_median.txt", 
                  "\t", escape_double = FALSE, trim_ws = TRUE)
exp <- exp[,-2]

gene_names <- exp$Hugo_Symbol
sample_names <- colnames(exp)[2:ncol(exp)]

expT <- as.data.frame(t(exp))
colnames(expT) <- gene_names
expT <- expT[-1,]
expT[] <- lapply(expT, function(x) as.numeric(as.character(x)))
expT <- as.data.frame(cbind(sample_id = sample_names,expT))
expT[,2:ncol(expT)] <- expT[,2:ncol(expT)] + 1
expT[,2:ncol(expT)] <- log(expT[,2:ncol(expT)], 2)

x <- data.matrix(expT[,2:ncol(expT)])
xx <- t(x)
xy <- xCellAnalysis(xx)

tcga_xcell <- t(xy)
tcga_xcell <- as.data.frame(tcga_xcell)
tcga_xcell <- as.data.frame(cbind(sample_id = row.names(tcga_xcell),
                                  tcga_xcell))

write.table(tcga_xcell, file="/Volumes/Lab_Gerke/TCGA/prad_tcga_pub/tcga_xcell.txt",
            row.names = FALSE, quote=FALSE,sep="\t")