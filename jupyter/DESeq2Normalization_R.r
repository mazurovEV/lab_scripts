
options(warn=-1)

library("DESeq2")

normalize <- function(matrix_path, ctl, exp, norm_matrix_path) {
    #countdata <- read.table("~/workspace/peaks_signal_matrix.csv", header=TRUE, row.names=1)
    countdata <- read.table(matrix_path, header=TRUE, row.names=1)

    countdata <- as.matrix(countdata)

    (condition <- factor(c(rep("ctl", ctl), rep("exp", exp))))

    (coldata <- data.frame(row.names=colnames(countdata), condition))

    dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
    dds <- estimateSizeFactors(dds)

    nz <- counts(dds, normalized = TRUE)
    nz <- round(nz, 2)
    #write.csv(nz, file="peaks_signal_matrix_normalized.csv", sep="\t") 
    write.table(nz, file=norm_matrix_path, quote=FALSE, sep="\t") 
}
