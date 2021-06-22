
options(warn=-1)

library('GenometriCorr')
library("IRanges")
library("GenomicRanges")
library("rtracklayer")
library("ChIPpeakAnno")

human.chrom.length <- c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,
                        138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,
                        83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415)
names(human.chrom.length) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
                               "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                               "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                               "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

#f это будет название без расширения
#
run_genometric <- function(rna, target, result_folder, margi_folder, peaks_folder) {
    possibleError <- tryCatch(
        corr <- genometric(rna, target, margi_folder, peaks_folder),
        error=function(e) e
    )

    if(!inherits(possibleError, "error")) {
        write.table(res.to.text(corr), paste("~/all_marks/", target, "/", result_folder, "/", rna, ".tsv", sep=''), row.names = FALSE, quote = FALSE, sep="\t")
        graphical.report(corr, pdffile = paste("~/all_marks/", target, "/", result_folder, "/", rna, ".pdf", sep=''))
        #visualize(corr, pdffile = paste("~/H3K27me3/genometric_results/vis_", f, ".pdf", sep=''))
    } else {
        print(possibleError)
    }
}

genometric <- function(rna, target, margi_folder, peaks_folder) {
    #print(paste("../MARGI/common_rnas/", unlist(strsplit(rna, "_"))[2], ".bed", sep=''))
    peaks <- toGRanges(paste("~/all_marks/", target, "/", peaks_folder, "/", rna, ".bed", sep=''), format="BED", header=FALSE)
    grid <- toGRanges(paste("~/MARGI/", target, "/", margi_folder, "/", rna, ".bed", sep=''), format="BED", header=FALSE)
        
    corr <- GenometriCorrelation(grid, peaks, chromosomes.length = human.chrom.length,
                                    permut.number = 1000,
                                    keep.distributions = TRUE, showProgressBar = TRUE)
    
    return(corr)
}

atac_genometric <- function(atac_peaks, h3k27c_peaks) {
    peaks <- toGRanges(atac_peaks, format="BED", header=FALSE)
    grid <- toGRanges(h3k27c_peaks, format="BED", header=FALSE)
        
    corr <- GenometriCorrelation(grid, peaks, chromosomes.length = human.chrom.length,
                                    permut.number = 1000,
                                    keep.distributions = TRUE, showProgressBar = TRUE)
    
    return(corr)
}

run_atac_genometric <- function(atac_peaks, h3k27c_peaks) {
    possibleError <- tryCatch(
        corr <- atac_genometric(atac_peaks, h3k27c_peaks),
        error=function(e) e
    )

    if(!inherits(possibleError, "error")) {
        write.table(res.to.text(corr), "~/ATAC_data/genometric_results.tsv", row.names = FALSE, quote = FALSE, sep="\t")
        graphical.report(corr, pdffile = "~/ATAC_data/genometric_results.pdf")
    } else {
        print(possibleError)
    }
}

#query population
#reference population
#relative Ks p-value
#relative ecdf deviation area
#relative ecdf area correlation
#relative ecdf deviation area p-value
#Scaled Absolute min. distance p−value
#Scaled Absolute min. lower tail
#Jaccard Measure p−value
#Jaccard Measure lower tail 
#Projection test p−value
#Projection test lower tail 
#Projection test observed to expected ratio
res.to.text <- function(corr) {
    chr.names = names(corr)
    query.population=c()
    reference.population=c()
    relative.distances.ks.p.value=c()
    relative.distances.ecdf.deviation.area=c()
    relative.distances.ecdf.area.correlation=c()
    relative.distances.ecdf.deviation.area.p.value=c()
    scaled.absolute.min.distance.sum.p.value=c()
    scaled.absolute.min.distance.sum.lower.tail=c()
    jaccard.measure.p.value=c()
    jaccard.measure.lower.tail=c()
    projection.test.p.value=c()
    projection.test.lower.tail=c()
    projection.test.obs.to.exp=c()
    for (name in chr.names) { 
        query.population <- append(query.population, corr[name][[1]][['query.population']])
        reference.population <- append(reference.population, corr[name][[1]][['reference.population']])
        relative.distances.ks.p.value <- append(relative.distances.ks.p.value, corr[name][[1]][['relative.distances.ks.p.value']])
        relative.distances.ecdf.deviation.area <- append(relative.distances.ecdf.deviation.area, corr[name][[1]][['relative.distances.ecdf.deviation.area']])
        relative.distances.ecdf.area.correlation <- append(relative.distances.ecdf.area.correlation, corr[name][[1]][['relative.distances.ecdf.area.correlation']])
        relative.distances.ecdf.deviation.area.p.value <- append(relative.distances.ecdf.deviation.area.p.value, corr[name][[1]][['relative.distances.ecdf.deviation.area.p.value']])
        scaled.absolute.min.distance.sum.p.value <- append(scaled.absolute.min.distance.sum.p.value, corr[name][[1]][['scaled.absolute.min.distance.sum.p.value']])
        scaled.absolute.min.distance.sum.lower.tail <- append(scaled.absolute.min.distance.sum.lower.tail, corr[name][[1]][['scaled.absolute.min.distance.sum.lower.tail']])
        jaccard.measure.p.value <- append(jaccard.measure.p.value, corr[name][[1]][['jaccard.measure.p.value']])
        jaccard.measure.lower.tail <- append(jaccard.measure.lower.tail, corr[name][[1]][['jaccard.measure.lower.tail']])
        projection.test.p.value <- append(projection.test.p.value, corr[name][[1]][['projection.test.p.value']])
        projection.test.lower.tail <- append(projection.test.lower.tail, corr[name][[1]][['projection.test.lower.tail']])
        projection.test.obs.to.exp <- append(projection.test.obs.to.exp, corr[name][[1]][['projection.test.obs.to.exp']])
    }

    data <- data.frame(chr.names, query.population, reference.population, relative.distances.ks.p.value, 
                       relative.distances.ecdf.deviation.area, relative.distances.ecdf.area.correlation, 
                      relative.distances.ecdf.deviation.area.p.value, scaled.absolute.min.distance.sum.p.value,
                      scaled.absolute.min.distance.sum.lower.tail, jaccard.measure.p.value, jaccard.measure.lower.tail,
                      projection.test.p.value, projection.test.lower.tail, projection.test.obs.to.exp, stringsAsFactors = FALSE)
    
    return(data)
}
