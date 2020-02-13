
dir_data=paste0("/users/ywang/10_25_2019/scRNA/data/GSE45719_RAW/")
dir_inpu = paste0("/users/ywang/10_25_2019/scRNA/data/GSE45719_RAW/")
dir_output="/users/ywang/10_25_2019/scRNA/data/"
dir_fig="/users/ywang/10_25_2019/scRNA/fig/"

library(scales)

percent_="0_1"

for(i in 1){
  #nsva=list_nsv[i]
  for(n_PC in c(16)){
    load(paste0(dir_output,percent_,"_percent_",n_PC,"PCs_",i,"_d1.RData"))
    load(paste0(dir_output,percent_,"_percent_",n_PC,"PCs_",i,"_d0.RData"))
    # load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d2.RData"))
    
    yrange=range(d1$y,d0$y)
    xrange=range(d1$x,d0$x)
    
    pdf(paste0(dir_fig,percent_,"_percent_ave_abs_",n_PC,"PCs_",i,"th_tissue_bef.pdf"),width = 2,height=1.5, pointsize = 6)
    #pdf(paste0(dir_fig,i,"th_tissue_center.pdf"),width = 2,height=1.5, pointsize = 6)
    
    par(oma = c(0,0,0,0), mar = c(5,5,1,1))
    
    plot(d1,  xaxt = "n",  bty = "n",yaxt = "n",col=alpha("grey80",1),main="",xlim=xrange,ylim=yrange, cex.lab=1, cex.axis=1.1,xlab="average expression level (log2PPKM)",ylab = "density")
    #"darkorange", "grey80"
    #lines(d2,col="red")
    lines(d0,col=alpha("black",1),lty=2)
    axis(side =1, at = c(0,5,10,15),cex.axis=1)
    axis(side =2, at = c(0.1,0.2),cex.axis=1)
    
    legend("topright",legend=c("background","observed"),col=c(alpha("black",1),"grey80"),lty=c(2,1,1), cex=0.7)
    dev.off()
  }
}

for(i in 1){
  #nsva=list_nsv[i]
  for(n_PC in c(4)){
    load(paste0(dir_output,percent_,"_percent_",n_PC,"PCs_",i,"_d1.RData"))
    load(paste0(dir_output,percent_,"_percent_",n_PC,"PCs_",i,"_d0.RData"))
    # load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d2.RData"))

    yrange=range(d1$y,d0$y)
    xrange=range(d1$x,d0$x)
    
    pdf(paste0(dir_fig,percent_,"_percent_ave_abs_",n_PC,"PCs_",i,"th_tissue_bef.pdf"),width = 2,height=1.5, pointsize = 6)
    #pdf(paste0(dir_fig,i,"th_tissue_center.pdf"),width = 2,height=1.5, pointsize = 6)
    
    par(oma = c(0,0,0,0), mar = c(5,5,1,1))

    plot(d1,  xaxt = "n",  bty = "n",yaxt = "n",col=alpha("grey80",1),main="",xlim=xrange,ylim=yrange, cex.lab=1, cex.axis=1.1,xlab="average expression level (log2PPKM)",ylab = "density")
    #"darkorange", "grey80"
    #lines(d2,col="red")
    lines(d0,col=alpha("black",1),lty=2)
    axis(side =1, at = c(0,5,10,15),cex.axis=1)
    axis(side =2, at = c(0.1,0.2),cex.axis=1)

    legend("topright",legend=c("background","observed"),col=c(alpha("black",1),"grey80"),lty=c(2,1,1), cex=0.7)
    dev.off()
  }
}


for(n_PC in c(16)){
  load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d1.RData"))
  load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d0.RData"))
  load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d2.RData"))
  
  yrange=range(d1$y,d2$y,d0$y)
  xrange=range(d1$x,d2$x,d0$x)
  
  pdf(paste0(dir_fig,percent_,"_percent_ave_abs_",n_PC,"PCs_",i,"th_tissue.pdf"),width = 2,height=1.5, pointsize = 6)
  #pdf(paste0(dir_fig,i,"th_tissue_center.pdf"),width = 2,height=1.5, pointsize = 6)
  
  par(oma = c(0,0,0,0), mar = c(5,5,1,1))
  plot(d1,col=alpha("grey80",1),  xaxt = "n",  bty = "n",yaxt = "n",main="",xlim=xrange,ylim=yrange, cex.lab=1, cex.axis=1.1,xlab="average expression level (log2RPKM)",ylab = "density")
  #"darkorange", "grey80"
  lines(d2,col=alpha("darkorange",1))
  lines(d0,col=alpha("black",1),lty=2)
  axis(side =1, at = c(0,5,10,15),cex.axis=1)
  axis(side =2, at = c(0.1,0.2),cex.axis=1)
  legend("topright",legend=c("background","observed","adjusted"),col=c(alpha("black",1),"grey80",alpha("darkorange",1)),lty=c(2,1,1), cex=0.7)
  dev.off()
}



for(n_PC in c(4)){
  load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d1.RData"))
  load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d0.RData"))
  load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d2.RData"))
  
  yrange=range(d1$y,d2$y,d0$y)
  xrange=range(d1$x,d2$x,d0$x)
  
  pdf(paste0(dir_fig,percent_,"_percent_ave_abs_",n_PC,"PCs_",i,"th_tissue.pdf"),width = 2,height=1.5, pointsize = 6)
  #pdf(paste0(dir_fig,i,"th_tissue_center.pdf"),width = 2,height=1.5, pointsize = 6)
  
  par(oma = c(0,0,0,0), mar = c(5,5,1,1))
  plot(d1,col=alpha("grey80",1),  xaxt = "n",  bty = "n",yaxt = "n",main="",xlim=xrange,ylim=yrange, cex.lab=1, cex.axis=1.1,xlab="average expression level (log2RPKM)",ylab = "density")
  #"darkorange", "grey80"
  lines(d2,col=alpha("darkorange",1))
  lines(d0,col=alpha("black",1),lty=2)
  axis(side =1, at = c(0,5,10,15),cex.axis=1)
  axis(side =2, at = c(0.1,0.2),cex.axis=1)
  legend("topright",legend=c("background","observed","adjusted"),col=c(alpha("black",1),"grey80",alpha("darkorange",1)),lty=c(2,1,1), cex=0.7)
  dev.off()
}
