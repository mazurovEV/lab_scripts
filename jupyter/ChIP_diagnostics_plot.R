chip_diagnostics  <-  function(bindata, plot_name, K = 200, rhat_vect, rhat_legend, rhat_ltype=rep(1,length(rhat_vect)), rhat_cols=rep("black",length(rhat_vect)),rhat_lwd=1,dens_ltype=c(1,3,4,5,6), dens_cols=c("black","red","green","blue","magenta"),dens_lwd=2,dens_quantiles=c(0,.25,.5,.75,1),xlim_max=0.99){
  
  #data
  tot.counts   <- bindata$chip_counts + bindata$input_counts
  tot.chip   <- sum(bindata$chip_counts)
  tot.input   <- sum(bindata$input_counts)
  tot.sum   <- sum(tot.counts)
  tot.cumsum    <- cumsum(tot.counts)
  ind.end.chr   <- cumsum(table(bindata$bin_chr))
  
  edec.tot.count.num   <- K
  
  #  2. Create EDEC bins
  
  brks.vec          <- sort(unique(c(tot.cumsum[ind.end.chr],seq(1,tot.sum,by=edec.tot.count.num))))
  cut.vec               <- cut(tot.cumsum,brks.vec)
  
  edec.input            <- sapply(split(bindata$input_counts,cut.vec),sum)
  edec.chip         <- sapply(split(bindata$chip_counts,cut.vec),sum)
  edec.length           <- sapply(split(bindata$chip_counts,cut.vec),length)
  edec.tot          <- edec.input + edec.chip
  
  #  2.a Find subset of "good" bin  
  plt.ind   <- (edec.tot.count.num*0.95 < edec.tot)  & (edec.tot  < quantile(edec.tot,xlim_max)) &
    (quantile(edec.length,0.01) < edec.length) &  (edec.length < quantile(edec.length,xlim_max))
  
  edec.input        <- edec.input[plt.ind]
  edec.chip     <- edec.chip[plt.ind]
  edec.length       <- edec.length[plt.ind]
  edec.tot      <- edec.tot[plt.ind]
  
  
  # 2.b Parition "good" bins into acording to bin length into 4 sparsity classes according to 
  # the length of the bin (= number of windows in bin)  short bins are dense and long bins are sparse:
  #  QQ4 is 4th quartile of longest and sparsest bins
  
  cut.vec               <- cut(edec.length,quantile(edec.length,prob=dens_quantiles))
  edec.input.QQ4        <- split(edec.input,cut.vec)[[4]]
  edec.chip.QQ4     <- split(edec.chip,cut.vec)[[4]]
  
  
  #  3. Draw logit densities for all "good" bins
  
  log.RR <- log(edec.chip / edec.input)
  dens.tot  <- density(log.RR)
  dens.QQ   <- lapply(split(log.RR,cut.vec),density)
  
  
  #************************ Diagnostics Plot
  pdf(paste(plot_name, ".pdf", sep=""))
  plot(exp(dens.tot$x), dens.tot$y,type="l",xlim=exp(quantile(log.RR,prob=c(0.001,xlim_max))),main="",ylab = "Density",xlab="Relative Risk ", lty=dens_ltype[1],col = dens_cols[1], log = "x")
  lines(exp(dens.QQ[[1]]$x),dens.QQ[[1]]$y/4,lty=dens_ltype[2], col=dens_cols[2],lwd=dens_lwd)
  lines(exp(dens.QQ[[2]]$x),dens.QQ[[2]]$y/4,lty=dens_ltype[3], col=dens_cols[3],lwd=dens_lwd)
  lines(exp(dens.QQ[[3]]$x),dens.QQ[[3]]$y/4,lty=dens_ltype[4], col=dens_cols[4],lwd=dens_lwd)
  lines(exp(dens.QQ[[4]]$x),dens.QQ[[4]]$y/4,lty=dens_ltype[5], col=dens_cols[5],lwd=dens_lwd)
  
  logit.pi.ch  <- (rhat_vect)
  
  for (i in seq(1:length(rhat_vect))){
    lines(logit.pi.ch[i],max(dens.tot$y),type='h',lty=rhat_ltype[i],col=rhat_cols[i],lwd=rhat_lwd) 
  }
  
  legend("topright",legend=c("Density",paste("bins with length in ",levels(cut.vec),sep=""),
                             paste(rhat_legend," (",rhat_vect,")",sep="")),lty=c(dens_ltype,rhat_ltype), 
         col=c(dens_cols,rhat_cols), 
         lwd=c(rep(dens_lwd,5),rep(1,length(rhat_vect))), cex=0.85) 
  dev.off()
  #************************  
}