
# 
# 1.extract the TF data
# 2.find one tissue - 1hr
# extract correlations of adj and ori
# Extract trank of adj and ori  - on cluster
# Scatter plot for ranks before/after spqn  - on cluster



#library("sva")
#library("recount", quietly = T)
#library("WGCNA", quietly = T)
library(matrixStats)
library(matrixStats)
library(WGCNA)
library(ggplot2)
library(ggridges)
library(SummarizedExperiment)
library(readr)

source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")
source("/users/ywang/10_25_2019/functions/functions_2d_quantile_diagnal_only.R")
source("/users/ywang/10_25_2019/functions/functions_evaluation_Aug_5.R")

source("/users/ywang/10_25_2019/WGCNA/code/functions/functions_construct_modules.R")


dir_data="/users/ywang/10_25_2019/data/" # counts
dir_input="/users/ywang/10_25_2019/corplot_ori_svapcs/data/" # cor_est 
dir_output="/users/ywang/10_25_2019/TF/data/"
dir_fig="/users/ywang/10_25_2019/TF/fig/"

names_TF=read.csv(paste0(dir_output,"Barrera_2016_Science_TableS2.csv"),header=T)
# > dim(names_TF)
# [1] 1254    2
# > ls(names_TF)
# [1] "Enseml.Gene.ID" "Gene.Symbol"  

i=1
n_PC=4
list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
               "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")

# load and input raw counts, normalize and filter low expression genes(average(log2cpm)>0)
name_tissue=list_tissues[[i]]

#load and input raw counts, normalize and filter low expression genes(average(log2cpm)>0)
name_var= paste0("gtex.",name_tissue)
load(paste0(dir_data, name_var,".rda" ))
data_raw=eval(as.name(name_var))
counts_raw=assay(data_raw)
counts_raw=counts_raw[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA")),]
list_filt=rpkm_filt(counts_raw,rowData(data_raw)$gene_length[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))])
log2rpkm_kp=list_filt$counts_log2rpkm_keep

# if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
# corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))

#gene_type_filt1=rowData(data_raw)$gene_type[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))]
#gene_type_filt12=gene_type_filt1[list_filt$keep][order(rowMeans(log2rpkm_kp))]
#genes_ensemble_keep2_format=genes_ensemble_keep2_format[which(gene_type_filt12=="protein_coding")]

gene_id_filt1=rownames(data_raw)[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))]
gene_id_filt12=gene_id_filt1[list_filt$keep][order(rowMeans(log2rpkm_kp))]

gene_id_filt12_2=strsplit(as.vector(gene_id_filt12),".",fixed=TRUE)
gene_id_filt12_3=matrix(unlist(gene_id_filt12_2), ncol = 2, byrow = TRUE)
gene_id_filt12_format=gene_id_filt12_3[,1]

###### get exp level of TF
TF_select=which(gene_id_filt12_format %in% names_TF$"Enseml.Gene.ID")

log2rpkm_kp_ordered=log2rpkm_kp[order(rowMeans(log2rpkm_kp)),]
range(log2rpkm_kp_ordered[TF_select])


###### get cor_ori and cor_est of TF-other
if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
cor_ori_TF=corr_ori[TF_select,]
save(cor_ori_TF,file=paste0(dir_output,"cor_ori_TF_others_",i,"th_tissue.RData"))
remove(corr_ori)
# 
load(paste0(dir_input,"cor_est_",n_PC,"_PCs_",i,"th_tissue.RData"))
cor_est_TF=cor_est[TF_select,]
save(cor_est_TF,file=paste0(dir_output,"cor_est_TF_others_",i,"th_tissue.RData"))
remove(cor_est)

# load(paste0(dir_output,"cor_est_TF_others_",i,"th_tissue.RData"))# cor_ori_TF
# load(paste0(dir_output,"cor_ori_TF_others_",i,"th_tissue.RData"))# cor_est_TF
# 
pdf(paste0(dir_fig,"TF_others_cor_d_",i,"th_tissue.pdf"))
plot(density(cor_ori_TF[upper.tri(cor_ori_TF)]))
lines(density(cor_est_TF[upper.tri(cor_est_TF)]),col="blue")
legend("topright",legend=c("ori","adj"),col=c("black","blue"),lty=c(1,1), cex=0.7)
dev.off()


###### get threshold exp levels of cor_ori and cor_est -- 0.1% and 1%
list_threshold_ratio=c(1,c(1:10)*3)/1000
ngene=nrow(log2rpkm_kp_ordered)
n_edges=ngene*(ngene-1)/2
library(dplyr)
threshold_value_ori=threshold_value_est=array(dim=length(list_threshold_ratio))
for(j in 1:length(list_threshold_ratio)){
  threshold_ratio=list_threshold_ratio[j]
  n_signal=round(n_edges*threshold_ratio)
  sort(n_signal)
  if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
  corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
  corr_ori=corr_ori[upper.tri(corr_ori)]
  corr_ori=abs(corr_ori)
  corr_ori=data.frame(corr_ori = corr_ori)
  top_cor=corr_ori %>% top_n(n_signal)
  #remove(corr_ori)
  threshold_value_ori[j]=min(top_cor)
  remove(top_cor)
  load(paste0(dir_input,"cor_est_",n_PC,"_PCs_",i,"th_tissue.RData"))
  cor_est=cor_est[upper.tri(cor_est)]
  cor_est=abs(cor_est)
  cor_est=data.frame(cor_est = cor_est)
  top_cor=cor_est %>% top_n(n_signal)
  #remove(corr_ori)
  threshold_value_est[j]=min(top_cor)

}

save(threshold_value_est,file=paste0(dir_output,"threshold_value_est_",i,"th_tissue.RData"))
save(threshold_value_ori,file=paste0(dir_output,"threshold_value_ori_",i,"th_tissue.RData"))

load(paste0(dir_output,"threshold_value_est_",i,"th_tissue.RData"))
load(paste0(dir_output,"threshold_value_ori_",i,"th_tissue.RData"))

###### check % of TF-other within the signals
# load(paste0(dir_output,"cor_ori_TF_others_",i,"th_tissue.RData"))# cor_ori_TF
# load(paste0(dir_output,"cor_est_TF_others_",i,"th_tissue.RData"))# cor_est_TF
# load(paste0(dir_output,"threshold_value_est_",i,"th_tissue.RData"))#threshold_value_est
# load(paste0(dir_output,"threshold_value_ori_",i,"th_tissue.RData"))#threshold_value_ori
list_threshold_ratio=c(1,c(1:10)*3)/1000

dim(cor_ori_TF)
edges_ori=edges_est=list()
for(j in 1:length(list_threshold_ratio)){
  
  cor_ori_TF_tmp=cor_ori_TF
  cor_ori_TF_tmp[which(! cor_ori_TF_tmp < threshold_value_ori[j])]=1
  cor_ori_TF_tmp[which( cor_ori_TF_tmp < threshold_value_ori[j])]=0
  edges_ori[[j]]=cor_ori_TF_tmp
  
  cor_est_TF_tmp=cor_est_TF
  cor_est_TF_tmp[which(! cor_est_TF_tmp < threshold_value_est[j])]=1
  cor_est_TF_tmp[which( cor_est_TF_tmp < threshold_value_est[j])]=0
  edges_est[[j]]=cor_est_TF_tmp
}
save(edges_ori,file=paste0(dir_output,"edges_ori_others_",i,"th_tissue.RData"))
save(edges_est,file=paste0(dir_output,"edges_est_others_",i,"th_tissue.RData"))
# _",i,"th_tissue
# 
# 
# # > sum(edges_ori[[1]])
# # [1] 5174
# # > sum(edges_ori[[2]])
# # [1] 54454
# # > sum(edges_est[[1]])
# # [1] 13751
# # > sum(edges_est[[2]])
# # [1] 147183
# 
edges_est_abs=edges_ori_abs=list()
for(j in 1:length(list_threshold_ratio)){
  
  cor_ori_TF_tmp=abs(cor_ori_TF)
  cor_ori_TF_tmp[which(! cor_ori_TF_tmp < threshold_value_ori[j])]=1
  cor_ori_TF_tmp[which( cor_ori_TF_tmp < threshold_value_ori[j])]=0
  edges_ori_abs[[j]]=cor_ori_TF_tmp
  
  cor_est_TF_tmp=abs(cor_est_TF)
  cor_est_TF_tmp[which(! cor_est_TF_tmp < threshold_value_est[j])]=1
  cor_est_TF_tmp[which( cor_est_TF_tmp < threshold_value_est[j])]=0
  edges_est_abs[[j]]=cor_est_TF_tmp
  
}

save(edges_ori_abs,file=paste0(dir_output,"edges_ori_abs_others_",i,"th_tissue.RData"))
save(edges_est_abs,file=paste0(dir_output,"edges_est_abs_others_",i,"th_tissue.RData"))


###### compare with PPI - changes of true edges after spqn
### extract a subset of PPI matrix that contains filtered TFs
# 
# dir_output2="/users/ywang/10_25_2019/WGCNA/data/"
# load(paste0(dir_output2,"mat_PPI_HURI.RData")) #mat_PPI_HURI
# TF_select=which(gene_id_filt12_format %in% names_TF$"Enseml.Gene.ID")
# mat_PPI_TF_others_HURI = transfer_ppi(mat_PPI_HURI,gene_id_filt12_format)[TF_select,]
# save(mat_PPI_TF_others_HURI,file=paste0(dir_output,"mat_PPI_TF_others_HURI_",i,"th_tissue.RData"))
# _",i,"th_tissue

# ### get the number of TP of TF in ori and est - cor>threhold
# load(paste0(dir_output,"mat_PPI_TF_others_HURI_",i,"th_tissue.RData"))#mat_PPI_TF_others_HURI
# load(paste0(dir_output,"edges_ori_others_",i,"th_tissue.RData"))#edges_ori
# load(paste0(dir_output,"edges_est_others_",i,"th_tissue.RData"))#edges_est

# > a=mat_PPI_TF_others_HURI+edges_ori[[1]]
# > sum(a==2)-sum(diag(mat_PPI_TF_others_HURI+edges_ori[[1]])==2)
# [1] 24
# > 
#   > b=mat_PPI_TF_others_HURI+edges_est[[1]]
# > sum(b==2)-sum(diag(mat_PPI_TF_others_HURI+edges_est[[1]])==2)
# [1] 30
# > 
#   > a=mat_PPI_TF_others_HURI+edges_ori[[2]]
# > sum(a==2)-sum(diag(mat_PPI_TF_others_HURI+edges_ori[[2]])==2)
# 
# [1] 66
# > b=mat_PPI_TF_others_HURI+edges_est[[2]]
# > sum(b==2)-sum(diag(mat_PPI_TF_others_HURI+edges_est[[2]])==2)
# [1] 131
# # sum(mat_PPI_TF_others_HURI)
# 3957


# ### get the number of TP of TF in ori and est - abs(cor)>threhold
# load(paste0(dir_output,"mat_PPI_TF_others_HURI_",i,"th_tissue.RData"))#mat_PPI_TF_others_HURI
# load(paste0(dir_output,"edges_ori_abs_others_",i,"th_tissue.RData"))#edges_ori
# load(paste0(dir_output,"edges_est_abs_others_",i,"th_tissue.RData"))#edges_est

# 
# > sum(a==2)-sum(diag(mat_PPI_TF_others_HURI+edges_ori[[1]])==2)
# [1] 24
# > b=mat_PPI_TF_others_HURI+edges_est[[1]]
# > sum(b==2)-sum(diag(mat_PPI_TF_others_HURI+edges_est[[1]])==2)
# [1] 30
# > a=mat_PPI_TF_others_HURI+edges_ori[[2]]
# > sum(a==2)-sum(diag(mat_PPI_TF_others_HURI+edges_ori[[2]])==2)
# [1] 66
# > b=mat_PPI_TF_others_HURI+edges_est[[2]]
# > sum(b==2)-sum(diag(mat_PPI_TF_others_HURI+edges_est[[2]])==2)
# [1] 131
