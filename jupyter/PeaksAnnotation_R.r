
options(warn=-1)

library("ChIPpeakAnno")
library("biomaRt")
library("GenomicFeatures")

peaksAnno<-function(filename, output_filename) {
    txdb <- makeTxDbFromGFF('../annotation/gencode.v31.annotation.gff3')
    
    anno <- toGRanges(txdb, format='gene')
    gr1 <- toGRanges(filename, format="BED", skip=1)

    overlaps.anno <- annotatePeakInBatch(gr1, AnnotationData=anno, output="overlapping", maxgap=1000L)
    print("count annotation, write to file")
    write.table(overlaps.anno, file=output_filename, quote=FALSE, sep="\t") 
}

#peaksAnno("../all_marks/H3K27ac/merged_peaks_first_in_biosample.bed", "../all_marks/H3K27ac/peaks_anno.csv")
