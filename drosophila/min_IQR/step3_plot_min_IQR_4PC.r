

# dir_fig="/users/ywang/10_25_2019/corplot_adj_4_svapcs/fig/"
# dir.create(dir_fig)
dir_output="/Users/yiwang/Dropbox/project_Kasper/July_28_2020/drosophila/out/"
dir_fig=dir_output
dir_input=dir_output

for(n_PC in 5:8){
  # load(paste0(dir_output,"min_IQR_4PC_",i,"th_tissue.RData"))
  # save(min_IQR,file=paste0(dir_output,"min_IQR_4PC_",i,"th_tissue.RData"))
  
  load(paste0(dir_output,"min_IQR_dro_rm_",n_PC,"PCs.RData"))
  # save(min_IQR,file=paste0(dir_output,"min_IQR_dro.RData"))
  # save(min_IQR,file=paste0(dir_output,"min_IQR_dro_rm_",n_PC,"PCs.RData"))
  # 
  range_y=range(min_IQR$IQR_vec_ori)
  
  pdf(paste0(dir_fig,"min_IQR_ori_4PC_dro_rm",n_PC,"PCs.pdf"))
  plot(min_IQR$min_mean,min_IQR$IQR_vec_ori,xlab=paste0("min(average(log2RPKM)), remove ",n_PC," PCs"),ylab="IQR",cex.lab=1.5,cex.axis=1.2,col="blue",ylim=range_y)
  dev.off()
  # 
  # pdf(paste0(dir_fig,"min_IQR_ori_4PC_",i,"th_tissue.pdf"))
  # plot(min_IQR$min_mean,min_IQR$IQR_vec_ori,xlab="min(average(log2RPKM))",ylab="IQR",cex.lab=1.5,cex.axis=1.2,col="blue",ylim=range_y)
  # dev.off()
  
  # pdf(paste0(dir_fig,"min_IQR_adj_4PC.pdf"))
  # plot(min_IQR$min_mean,min_IQR$IQR_vec_adj,xlab="min(average(log2CPM))",ylab="IQR",cex.lab=1.5,cex.axis=1.2,col="red",ylim=range_y)
  # dev.off()
}z



