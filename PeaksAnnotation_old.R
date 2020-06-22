#biocLite("ChIPpeakAnno")
#biocLite("EnsDb.Hsapiens.v75")

library("ChIPpeakAnno")
library(EnsDb.Hsapiens.v75) ##(hg19)

peaksAnno<-function(filename, output_filename) {
  annoData <- toGRanges(EnsDb.Hsapiens.v75, feature="gene")

  gr1 <- toGRanges(filename, format="BED", skip=1)

  overlaps.anno <- annotatePeakInBatch(gr1, AnnotationData=annoData, output="overlapping", maxgap=1000L)
  overlaps.anno[2, 1:3]

  write.csv2(overlaps.anno, file = output_filename, sep = ";")
}
  