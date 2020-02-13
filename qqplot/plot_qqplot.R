library(matrixStats)
library(WGCNA)
library(ggplot2)
library(ggridges)

source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")
source("/users/ywang/10_25_2019/functions/functions_evaluation_Aug_5.R")

dir_data="/users/ywang/10_25_2019/data/"
dir_fig="/users/ywang/10_25_2019/qqplot_4pc_est/fig/"

dir_input=paste0("/users/ywang/10_25_2019/corplot_ori_svapcs/data/")
qqplot_density4<-function(cor_matrix,i,j,n,xlim_,ylim_){
  grp_loc=get_grp_loc(cor_matrix)
  cor_tmp=cor_matrix[grp_loc[[i]],grp_loc[[i]]]
  cor_tmp=cor_tmp[upper.tri(cor_tmp)]
  cor_ref=cor_matrix[grp_loc[[j]],grp_loc[[j]]]
  cor_ref_vec=cor_ref[upper.tri(cor_ref)]
  qqplot_pre_=qqplot_pre(as.vector(cor_tmp),cor_ref_vec)
  rasterPlotFun2(qqplot_pre_$x,qqplot_pre_$y,i,i,j,j,n,xlim_,ylim_)
}

qqplot_pre<-function (x, y) 
{
  sx <- sort(x)
  sy <- sort(y)
  lenx <- length(sx)
  leny <- length(sy)
  if (leny < lenx) 
    sx <- approx(1L:lenx, sx, n = leny)$y
  if (leny > lenx) 
    sy <- approx(1L:leny, sy, n = lenx)$y
  
  list(x = sx, y = sy)
}
  
rasterPlotFun2<-function(r2s.1, r2s.2,i1,j1,i2,j2,n,xlim_,ylim_){
  r2s <- cbind(r1=r2s.1, r2=r2s.2)
  r <- raster(nrows = 500, ncols = 500, xmn = r2s.1[1], xmx = r2s.1[length(r2s.1)], ymn = r2s.2[1], ymx = r2s.2[length(r2s.2)],
              resolution=c(0.008,0.008))
  r2s.rast <- rasterize(r2s, r, fun = "count")
  r2s.rast2 <- r2s.rast
  values(r2s.rast2)[which(values(r2s.rast2) > 10)] <- 10
  par(oma = c(0,0,0,0), mar = c(5,5,0.6,0.6))
  #colR <- colorRampPalette(c("#E5E5E5", "black"))(10)
  colR <- "black"
  #qq_=qqplot(cor_mid,cor_edge,ylab=paste0("correlations in (",i1,",",j1,") - all"),xlab=paste0("correlations in (",i1,",",j1,") - edges"),cex.lab=1.5,cex.axis=1.5,cex=0.2)
  plot(0,0, type = "n", xaxt = "n", yaxt = "n", bty = "n",main=paste0(" (",i1,",",j1,") "),
       xlim = xlim_, ylim=ylim_, xlab = "", ylab = "", maxpixels=n,cex=20,pch=20)
  axis(side =1, at = c(-0.5,0,0.5),cex.axis=1.5)
  axis(side =2, at = c(-0.5,0,0.5),cex.axis=1.5)
  title(xlab=paste0("correlations in (",i1,",",j1,") "),cex.lab=2, ylab = paste0("correlations in (",i2,",",j2,") "))
  abline(0,1, col = "blue",alpha=0.5)
  rasterImage(as.raster(r2s.rast2, col = colR, maxpixels = n,cex=20,pch=20),cex=20,pch=20,
              xleft = r2s.1[1], xright = r2s.1[length(r2s.1)], ybottom = r2s.2[1], ytop = r2s.2[length(r2s.2)],
              interpolate = FALSE)
  
}

  
library(SummarizedExperiment)
i=1

list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
               "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")
name_tissue=list_tissues[[i]]
#load and input raw counts, normalize and filter low expression genes(average(log2cpm)>0)
name_var= paste0("gtex.",name_tissue)
load(paste0(dir_data, name_var,".rda" ))
data_raw=eval(as.name(name_var))
counts_raw=assay(data_raw)
counts_raw=counts_raw[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA")),]
list_filt=rpkm_filt(counts_raw,rowData(data_raw)$gene_length[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))])
log2rpkm_kp=list_filt$counts_log2rpkm_keep

detach("package:SummarizedExperiment", unload=TRUE)


#########################################################################################################
#sva PCs
library(raster)
load(paste0("/users/ywang/10_25_2019/corplot_adj_4_svapcs/fig/","nsv_",i,"_th_tissue.RData"))
dir_figure = paste0("/users/ywang/10_25_2019/corplot_adj_4_svapcs/fig/")
load(paste0(dir_figure,"nsv_",i,"_th_tissue.RData"))
n_PC=nsv
# 
# # calculate cor_ori
# if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
# corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
# xlim_ori_est=range(corr_ori[upper.tri(corr_ori)])
# ylim_ori_est=xlim_ori_est
# ylim_ori_est
# 
# for(k in 1:10){
#     for(j in 9){
#       pdf(paste0(dir_fig,j,"_",k,"_",n_PC,"_PCs_ori.pdf"))
#       qqplot_density4(corr_ori,j,k,50000000,xlim_ori_est,ylim_ori_est)
#       dev.off()
#   }
# }
# 
# 
# for(k in 9){
#     for(j in 1:10){
#       pdf(paste0(dir_fig,j,"_",k,"_",n_PC,"_PCs_ori.pdf"))
#       qqplot_density4(corr_ori,j,k,50000000,xlim_ori_est,ylim_ori_est)
#       dev.off()
#   }
# }
# remove(corr_ori)
# 
load(paste0(dir_input,"cor_est_",n_PC,"_PCs_",i,"th_tissue.RData"))
for(j in 1:10){
    for(k in 9){
      pdf(paste0(dir_fig,j,"_",k,"_",n_PC,"_PCs_est.pdf"))
      qqplot_density4(cor_est,j,k,50000000,xlim_ori_est,ylim_ori_est)
      dev.off()
  }
}

for(k in 1:10){
    for(j in 9){
      pdf(paste0(dir_fig,j,"_",k,"_",n_PC,"_PCs_est.pdf"))
      qqplot_density4(cor_est,j,k,50000000,xlim_ori_est,ylim_ori_est)
      dev.off()
  }
}
remove(cor_est)

#########################################################################################################
#4 PCs
n_PC=4
# if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
# corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
# xlim_ori_est=range(corr_ori[upper.tri(corr_ori)])
# ylim_ori_est=xlim_ori_est
# ylim_ori_est
# 
# library(raster)
# for(k in 1:10){
#   for(j in 9){
#     pdf(paste0(dir_fig,j,"_",k,"_",n_PC,"_PCs_ori.pdf"))
#     qqplot_density4(corr_ori,j,k,50000000,xlim_ori_est,ylim_ori_est)
#     dev.off()
#   }
# }


# for(k in 9){
#   for(j in 1:10){
#     pdf(paste0(dir_fig,j,"_",k,"_",n_PC,"_PCs_ori.pdf"))
#     qqplot_density4(corr_ori,j,k,50000000,xlim_ori_est,ylim_ori_est)
#     dev.off()
#   }
# }
# remove(corr_ori)


load(paste0(dir_input,"cor_est_",n_PC,"_PCs_",i,"th_tissue.RData"))

for(j in 1:10){
  for(k in 9){
    pdf(paste0(dir_fig,j,"_",k,"_",n_PC,"_PCs_est.pdf"))
    qqplot_density4(cor_est,j,k,50000000,xlim_ori_est,ylim_ori_est)
    dev.off()
  }
}

for(k in 1:10){
  for(j in 9){
    pdf(paste0(dir_fig,j,"_",k,"_",n_PC,"_PCs_est.pdf"))
    qqplot_density4(cor_est,j,k,50000000,xlim_ori_est,ylim_ori_est)
    dev.off()
  }
}



