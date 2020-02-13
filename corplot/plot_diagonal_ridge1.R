library(matrixStats)
library(WGCNA)
library(ggplot2)
library(ggridges)
library(SummarizedExperiment)

source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")
source("/users/ywang/10_25_2019/functions/functions_evaluation_Aug_5.R")
dir_data="/users/ywang/10_25_2019/data/"
dir_fig="/users/ywang/10_25_2019/corplot_ori_0pc/"
dir_fig="/users/ywang/10_25_2019/corplot_ori_0pc/diagonal_adj_ridge4/"
dir_fig="/users/ywang/10_25_2019/corplot_ori_0pc/diagonal_adj_ridge4/fig/"
dir_input="/users/ywang/10_25_2019/corplot_ori_svapcs/data/"


plot_diagonal_ridge<-function(cor_matrix,ngrp=10){
  ngene = ncol(cor_matrix)
  grp_label=cut(1:ngene,ngrp)
  grp_loc=split(1:ngene,grp_label) 
  for(i in 1:10){
    cor_tmp=cor_matrix[grp_loc[[i]],grp_loc[[i]]]
    cor_tmp=cor_tmp[upper.tri(cor_tmp)]
    if(i==1){cor_vec_all= cbind(as.numeric(cor_tmp),rep(i,length(as.numeric(cor_tmp))))
    names(cor_vec_all)=c("correlation","group")}else{
      cor_vec_all=rbind(cor_vec_all,cbind(as.numeric(cor_tmp),rep(i,length(as.numeric(cor_tmp)))))}
  }
  cor_vec_all=data.frame(cor_vec_all)
  names(cor_vec_all)=c("correlation","group")
  cor_vec_all$group=as.factor(cor_vec_all$group)
  
  ggplot(cor_vec_all, aes(x = `correlation`, y = `group`)) +
    geom_density_ridges2(fill="blue")+
    theme_ridges(grid = TRUE) + theme(axis.title.x = element_blank())+ 
    geom_vline(xintercept = 0, linetype="dotted", color = "black", size=0.3)
  
}


########### for cor_ori
list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
               "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")
i=1
name_tissue=list_tissues[[i]]
#dir_figure = paste0("/users/ywang/8_25_2019/figures/")

#load and input raw counts, normalize and filter low expression genes(average(log2cpm)>0)
name_var= paste0("gtex.",name_tissue)
load(paste0(dir_data, name_var,".rda" ))
data_raw=eval(as.name(name_var))
counts_raw=assay(data_raw)
counts_raw=counts_raw[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA")),]
list_filt=rpkm_filt(counts_raw,rowData(data_raw)$gene_length[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))])
log2rpkm_kp=list_filt$counts_log2rpkm_keep

# calculate cor_ori

dir_figure = paste0("/users/ywang/10_25_2019/corplot_adj_4_svapcs/fig/")
load(paste0(dir_figure,"nsv_",i,"_th_tissue.RData"))#nsv

for(n_PC in c(4,nsv)){
  
  if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
  corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
  #ave_logcpm_kp=rowMeans(log2cpm_kp)[order(rowMeans(log2cpm_kp))]
  
  pdf(paste0(dir_fig,"diagonal_ori_ridge_",n_PC,"PCs_",i,"th_tissue_ori.pdf"),width=7,height = 7)
  print(plot_diagonal_ridge(corr_ori))
  dev.off()
  
  load(paste0(dir_input,"cor_est_",n_PC,"_PCs_",i,"th_tissue.RData"))
  pdf(paste0(dir_fig,"diagonal_ori_ridge_",n_PC,"PCs_",i,"th_tissue_adj.pdf"),width=7,height = 7)
  print(plot_diagonal_ridge(cor_est))
  dev.off()
  
}








