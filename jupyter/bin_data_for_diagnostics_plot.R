library(NCIS)
source("ChIP_diagnostics_plot.R")

countFactorAndDiagnosticPlot<-function(path, control_path, file_chip, file_input, plot_name) {
  
  frag.len = 200#get from MACS?
  
  chip=readBEDbyChromosomes(path, file_chip)
  input=readBEDbyChromosomes(control_path, file_input)
  
  print("read files, count factor")
  
  shift.size <- round(frag.len/2)
  ncis <- NCIS:::NCIS.internal(chip, input, shift.size=shift.size, min.binsize=100, max.binsize=20000, 
                               binsize.shift=100, min.stop.binsize=100, chr.vec=NULL, 
                               chr.len.vec=NULL, quant=0.75)
  
  print(paste("Count norm factor:", ncis['est']))
  print("Prepare for diagnostics...")
  
  binsize=200
  shift.size=100
  
  bindata1 <- NCIS:::bin.data(chip.pos=chip, input.pos=input, binsize, shift.size=shift.size, shift.half.size=FALSE, zero.filter=FALSE, by.strand=FALSE, chr.end.max=NULL, by.chr=TRUE)
  
  # reshape the output of ??the bin.data function 
  nchr=length(names(bindata1$chip)) 
  chr_list=names(bindata1$chip)
  chr_vec=NULL
  
  for (ichr in 1:nchr){
    chr=chr_list[ichr]
    windows_chr=length(bindata1$chip[[chr]])
    temp=rep(chr,windows_chr)
    chr_vec=c(chr_vec, temp)
  }
  
  bindata=list(chip_counts=unlist(bindata1$chip, use.names = FALSE), 
               input_counts=unlist(bindata1$input, use.names = FALSE),
               bin_chr=chr_vec)
  
 
  write(c(ncis[['est']], ncis[['binsize.est']], ncis[['r.seq.depth']], ncis[['pi0']]), file = paste("/data/mazurovev/norm_factors/plots/",file_chip,".txt",sep=""))
  print("Plot...")
  chip_diagnostics(bindata, plot_name, K=200 , rhat_vect=c(ncis['est']), rhat_legend=c("NCIS"), rhat_cols=c("navy"))
  return(c(ncis['est'], ncis['binsize.est'], ncis['r.seq.depth'], ncis['pi0']))
}

readBEDbyChromosomes<-function(path, file_prefix) {
  files <- list.files(path=path, file_prefix)
 
  res <- list()
  for(f in files) {
    
    res <- append(res, NCIS:::read.BED(paste(path, f, sep="")))
    
  }
  
  return(res)
}

