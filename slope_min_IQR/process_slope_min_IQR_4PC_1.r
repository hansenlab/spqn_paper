library(matrixStats)
library(WGCNA)
library(SummarizedExperiment)

source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")
source("/users/ywang/10_25_2019/functions/functions_evaluation_Aug_5.R")

dir_data="/users/ywang/10_25_2019/data/"

dir_input="/users/ywang/10_25_2019/slope_min_QIR/data/"
dir_output="/users/ywang/10_25_2019/slope_min_QIR/data/"

i=1

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
ave_logrpm_kp=rowMeans(log2rpkm_kp)[order(rowMeans(log2rpkm_kp))]


list_IQR_mat=list()

for(nPC in 0:50){
  
  if(nPC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,nPC)}
  
  corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
  
  list_IQR_mat[[nPC+1]] = cal_m_sd_IQR(corr_ori,ave_logrpm_kp)
  print(nPC)
}


save(list_IQR_mat,file=paste0(dir_output,i,"th_tissue_list_IQR_mat.RData"))

