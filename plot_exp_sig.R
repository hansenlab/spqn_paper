dir_input="/users/ywang/10_25_2019/exp_signal_4PC/data/"
#dir_output="/users/ywang/10_25_2019/exp_signal_4PC/data/"
dir_fig="/users/ywang/10_25_2019/exp_signal_4PC/fig/"
dir.create(dir_fig)
library(scales)

percent_="0_1"

yrange_all_4PC=xrange_all_4PC=c()
for(i in 1:9){
  #nsva=list_nsv[i]
  for(n_PC in c(4)){#nsva
    load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d1.RData"))
    load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d0.RData"))
    load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d2.RData"))
    
    yrange=range(d1$y,d2$y,d0$y)
    yrange_all_4PC=c(yrange_all_4PC,yrange)
    xrange=range(d1$x,d2$x,d0$x)
    xrange_all_4PC=c(xrange_all_4PC,xrange)
    
  }
}


yrange_all_svaPC=xrange_all_svaPC=c()
for(i in 1:9){
  #nsva=list_nsv[i]
  dir_figure = paste0("/users/ywang/10_25_2019/corplot_adj_4_svapcs/fig2/")
  load(paste0(dir_figure,"nsv_",i,"_th_tissue.RData"))
  nsva=nsv
  for(n_PC in nsva){#nsva
    load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d1.RData"))
    load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d0.RData"))
    load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d2.RData"))
    
    yrange=range(d1$y,d2$y,d0$y)
    yrange_all_svaPC=c(yrange_all_svaPC,yrange)
    xrange=range(d1$x,d2$x,d0$x)
    xrange_all_svaPC=c(xrange_all_svaPC,xrange)
  }
}

xrange_all_4PC=range(xrange_all_4PC)
yrange_all_4PC=range(yrange_all_4PC)
xrange_all_svaPC=range(xrange_all_svaPC)
yrange_all_svaPC=range(yrange_all_svaPC)

for(i in 1:9){
  #nsva=list_nsv[i]
  dir_figure = paste0("/users/ywang/10_25_2019/corplot_adj_4_svapcs/fig2/")
  load(paste0(dir_figure,"nsv_",i,"_th_tissue.RData"))
  nsva=nsv
  for(n_PC in nsva){
    load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d1.RData"))
    load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d0.RData"))
    load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d2.RData"))
    
    yrange=range(d1$y,d2$y,d0$y)
    pdf(paste0(dir_fig,percent_,"_percent_ave_abs_",n_PC,"PCs_",i,"th_tissue.pdf"),width = 2,height=1.5, pointsize = 6)
    #pdf(paste0(dir_fig,i,"th_tissue_center.pdf"),width = 2,height=1.5, pointsize = 6)
    
    par(oma = c(0,0,0,0), mar = c(5,5,1,1))
    plot(d1,col=alpha("grey80",1),  xaxt = "n",  bty = "n",yaxt = "n",main="",xlim=xrange_all_4PC,ylim=yrange_all_4PC, cex.lab=1, cex.axis=1.1,xlab="average expression level (log2RPKM)",ylab = "density")
    #"darkorange", "grey80"
    lines(d2,col=alpha("darkorange",1))
    lines(d0,col=alpha("black",1),lty=2)
    axis(side =1, at = c(0,5,10,15),cex.axis=1)
    axis(side =2, at = c(0.1,0.2),cex.axis=1)
    legend("topright",legend=c("background","observed","adjusted"),col=c(alpha("black",1),"grey80",alpha("darkorange",1)),lty=c(2,1,1), cex=0.7)
    dev.off()
  }
}




for(i in 1:9){
  #nsva=list_nsv[i]
  for(n_PC in c(4)){
    load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d1.RData"))
    load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d0.RData"))
    load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d2.RData"))
    
    yrange=range(d1$y,d2$y,d0$y)
    pdf(paste0(dir_fig,percent_,"_percent_ave_abs_",n_PC,"PCs_",i,"th_tissue.pdf"),width = 2,height=1.5, pointsize = 6)
    #pdf(paste0(dir_fig,i,"th_tissue_center.pdf"),width = 2,height=1.5, pointsize = 6)
    
    par(oma = c(0,0,0,0), mar = c(5,5,1,1))
    
    plot(d1,  xaxt = "n",  bty = "n",yaxt = "n",col=alpha("grey80",1),main="",xlim=xrange_all_svaPC,ylim=yrange_all_svaPC, cex.lab=1, cex.axis=1.1,xlab="average expression level (log2RPKM)",ylab = "density")
    #"darkorange", "grey80"
    lines(d2,col=alpha("darkorange",1))
    lines(d0,col=alpha("black",1),lty=2)
    axis(side =1, at = c(0,5,10,15),cex.axis=1)
    axis(side =2, at = c(0.1,0.2),cex.axis=1)
    
    legend("topright",legend=c("background","observed","adjusted"),col=c(alpha("black",1),"grey80",alpha("darkorange",1)),lty=c(2,1,1), cex=0.7)
    dev.off()
  }
}
# plot(d1,col=alpha("grey80",0.5),  xaxt = "n",  bty = "n",yaxt = "n",main="",xlim=xrange_all_4PC,ylim=yrange_all_4PC, cex.lab=1.5, cex.axis=1.1,xlab="average expression level (log2RPKM)",ylab = "density")
# #"darkorange", "grey80"
# lines(d2,col=alpha("darkorange",0.5))
# lines(d0,col=alpha("black",0.5))
# axis(side =1, at = c(0,5,10,15),cex.axis=1)
# axis(side =2, at = c(0.1,0.2),cex.axis=1)
# legend("topright",legend=c("background","observed","adjusted"),col=c("black","grey80","darkorange"),lty=1, cex=1.3)
# dev.off()



# 
# for(i in 1){
#   #nsva=list_nsv[i]
#   dir_figure = paste0("/users/ywang/10_25_2019/corplot_adj_4_svapcs/fig/")
#   load(paste0(dir_figure,"nsv_",i,"_th_tissue.RData"))
#   nsva=nsv
#   for(n_PC in nsva){
#     load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d1.RData"))
#     load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d0.RData"))
#     load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d2.RData"))
# 
#     yrange=range(d1$y,d2$y,d0$y)
#     pdf(paste0(dir_fig,percent_,"_percent_ave_abs_",n_PC,"PCs_",i,"th_tissue_bef.pdf"),width = 2,height=1.5, pointsize = 6)
#     #pdf(paste0(dir_fig,i,"th_tissue_center.pdf"),width = 2,height=1.5, pointsize = 6)
#     par(oma = c(0,0,0,0), mar = c(5,5,1,1))
#     plot(d1,col="grey80",  xaxt = "n",  bty = "n",yaxt = "n",main="",xlim=xrange_all_4PC,ylim=yrange_all_4PC, cex.lab=1.5, cex.axis=1.1,xlab="average expression level (log2CPM)",ylab = "density")
#     #"darkorange", "grey80"
#     #lines(d2,col="red")
#     lines(d0)
#     axis(side =1, at = c(0,5,10,15),cex.axis=1)
#     axis(side =2, at = c(0.1,0.2),cex.axis=1)
#     legend("topright",legend=c("background","observed"),col=c("black","blue"),lty=1, cex=1.3)
#     dev.off()
#   }
# }
# 
# 
# 
# 
for(i in 1){
  #nsva=list_nsv[i]
  for(n_PC in c(4)){
    load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d1.RData"))
    load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d0.RData"))
    load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d2.RData"))

    yrange=range(d1$y,d2$y,d0$y)
    pdf(paste0(dir_fig,percent_,"_percent_ave_abs_",n_PC,"PCs_",i,"th_tissue_bef.pdf"),width = 2,height=1.5, pointsize = 6)
    #pdf(paste0(dir_fig,i,"th_tissue_center.pdf"),width = 2,height=1.5, pointsize = 6)
    
    par(oma = c(0,0,0,0), mar = c(5,5,1,1))

    plot(d1,  xaxt = "n",  bty = "n",yaxt = "n",col=alpha("grey80",1),main="",xlim=xrange_all_svaPC,ylim=yrange_all_svaPC, cex.lab=1, cex.axis=1.1,xlab="average expression level (log2PPKM)",ylab = "density")
    #"darkorange", "grey80"
    #lines(d2,col="red")
    lines(d0,col=alpha("black",1),lty=2)
    axis(side =1, at = c(0,5,10,15),cex.axis=1)
    axis(side =2, at = c(0.1,0.2),cex.axis=1)

    legend("topright",legend=c("background","observed"),col=c(alpha("black",1),"grey80"),lty=c(2,1,1), cex=0.7)
    dev.off()
  }
}
# 

# plot(d1,col=alpha("grey80",0.5),  xaxt = "n",  bty = "n",yaxt = "n",main="",xlim=xrange_all_4PC,ylim=yrange_all_4PC, cex.lab=1.5, cex.axis=1.1,xlab="average expression level (log2RPKM)",ylab = "density")
# #"darkorange", "grey80"
# lines(d2,col=alpha("darkorange",0.5))
# lines(d0,col=alpha("black",0.5))
# axis(side =1, at = c(0,5,10,15),cex.axis=1)
# axis(side =2, at = c(0.1,0.2),cex.axis=1)
# legend("topright",legend=c("background","observed","adjusted"),col=c("black","grey80","darkorange"),lty=1, cex=1.3)
# dev.off()



# for(i in 1:9){
#   nsva=list_nsv[i]
#   for(n_PC in c(4,nsva)){
#     load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d1.RData"))
#     load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d0.RData"))
#     load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_d2.RData"))
#   
#     yrange=range(d1$y,d2$y,d0$y)
#     pdf(paste0(dir_fig,percent_,"_percent_ave_abs_",n_PC,"PCs_",i,"th_tissue.pdf"))
#     plot(d1,col="blue",main="",ylim=yrange, cex.lab=1.5, cex.axis=1.5,xlab="average expression level (log2CPM)",ylab = "density")
#     lines(d2,col="red")
#     lines(d0)
#     legend("topright",legend=c("background","observed","adjusted"),col=c("black","blue","red"),lty=1, cex=1.2)
#     dev.off()
#   }
# }

