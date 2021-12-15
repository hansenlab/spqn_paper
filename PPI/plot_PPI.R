library(readr)
library(matrixStats)

dir_input="/PPI/data/"
dir_fig="/PPI/fig/"
dir_input2="/exp_signal_4PC/data/"

percent_="0_1"

###########
########### 4PCs

yrange_all_4PC=xrange_all_4PC=c()
for(i in 1:9){
  for(n_PC in c(4)){
    load(paste0(dir_input,i,"_th_tissue_dist_ave_exp_PPI.RData" ))#dist_ave_exp_PPI
    load(paste0(dir_input2,percent_,"_percent_",n_PC,"PCs_",i,"_d0.RData"))

    yrange=range(dist_ave_exp_PPI$y,d0$y)
    yrange_all_4PC=c(yrange_all_4PC,yrange)
    xrange=range(dist_ave_exp_PPI$x,d0$x)
    xrange_all_4PC=c(xrange_all_4PC,xrange)

  }
}


yrange_all_svaPC=xrange_all_svaPC=c()
for(i in 1:9){
  dir_figure = paste0("/users/ywang/10_25_2019/corplot_adj_4_svapcs/fig/")
  load(paste0(dir_figure,"nsv_",i,"_th_tissue.RData"))
  nsva=nsv
  for(n_PC in nsva){#nsva
    load(paste0(dir_input,i,"_th_tissue_dist_ave_exp_PPI.RData" ))#dist_ave_exp_PPI
    load(paste0(dir_input2,percent_,"_percent_",n_PC,"PCs_",i,"_d0.RData"))

    yrange=range(dist_ave_exp_PPI$y,d0$y)
    yrange_all_svaPC=c(yrange_all_svaPC,yrange)
    xrange=range(dist_ave_exp_PPI$x,d0$x)
    xrange_all_svaPC=c(xrange_all_svaPC,xrange)
  }
}


xrange_all_4PC=range(xrange_all_4PC)
yrange_all_4PC=range(yrange_all_4PC)
xrange_all_svaPC=range(xrange_all_svaPC)
yrange_all_svaPC=range(yrange_all_svaPC)




########### PPI vs. background distribution


for(i in 1:9){
  dir_figure = paste0("/users/ywang/10_25_2019/corplot_adj_4_svapcs/fig/")
  load(paste0(dir_figure,"nsv_",i,"_th_tissue.RData"))
  n_PC=nsv
  
  load(paste0(dir_input,i,"_th_tissue_dist_ave_exp_PPI.RData" ))#dist_ave_exp_PPI
  load(paste0(dir_input2,percent_,"_percent_",n_PC,"PCs_",i,"_d0.RData"))

  pdf(paste0(dir_fig,"PPI_",i,"_th_tissue.pdf"))

  par(oma = c(0,0,0,0), mar = c(5,5,0.6,0.6))

  plot(d0,  xaxt = "n",  bty = "n",yaxt = "n",col="black",main="",xlim=xrange_all_svaPC,ylim=yrange_all_svaPC, cex.lab=1.5, cex.axis=1.1,xlab="average expression level (log2RPKM)",ylab = "density")
  lines(dist_ave_exp_PPI,col="blue")
  axis(side =1, at = c(0,5,10,15),cex.axis=1)
  axis(side =2, at = c(0.1,0.2),cex.axis=1)

  legend("topright",legend=c("background","PPI"),col=c("black","blue"),lty=1, cex=1.3)

  dev.off()
}



