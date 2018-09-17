# xCellCalls

`tcga_code.R`

Makes xCell calls for the 333 TCGA prostate samples from the 2015 Cell paper. Data can be downloaded from: http://www.cbioportal.org/

`gse46691_code.R`

Makes two sets of xCell calls for GSE46691. One uses th median expression value from multiple probes for each gene(`GSE46691_xcell.txt`) while the other uses the probe with the max expression value for each gene (`GSE46691_xcell_max.txt`). Probe-level data can be downloaded from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46691 and the most up to date annotation files can be downloaded from Affymetrix at: http://www.affymetrix.com/support/technical/byproduct.affx?product=huexon-st. 
