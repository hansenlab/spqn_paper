
########### # get the adjusted correlation for shared expressed genes only, for each tissue of gtex data
########### # by applying spqn on cor_ori that containing shared genes only
library(matrixStats)
library(WGCNA)
library(spqn)
library(SummarizedExperiment)
source("/functions/functions_2d_quantile_no_surf_July_18.R")
source("/functions/plot_signal_condition_exp.R")

source("/functions/functions_evaluation_Aug_5.R")

dir_data="/data/"
dir_output="/pool_tissues/output_4pc/"
dir.create(dir_output)
setwd(dir_output)

i=1


list_tissues_combinations=list()
list_tissues_combinations[[1]]=c(1:3)
list_tissues_combinations[[2]]=c(4:6)
list_tissues_combinations[[3]]=c(7:9)

#### 
for(i in c(1:3)){
  print('------------------------------------------')
  
  print(i)
  
  n_PC=4
  
  ########### for each tissue, format the data 
  list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
                 "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")
 
  
  names_tissue=list_tissues_combinations[[i]]
  
  #load and input raw counts, normalize and filter low expression genes(average(log2cpm)>0)
  for(k in 1:3){
    name_tissue=list_tissues[names_tissue[k]]
    name_var= paste0("gtex.",name_tissue)
    load(paste0(dir_data, name_var,".rda" ))
    data_raw=eval(as.name(name_var))
    counts_raw=assay(data_raw)
    set.seed(111)
    rd_=sample(1:ncol(counts_raw),100)
    counts_raw=counts_raw[,rd_]
    if(k==1){counts_raw_combine=counts_raw}else{
      counts_raw_combine=cbind(counts_raw_combine,counts_raw)
    }
  }

  
  
  counts_raw_combine=counts_raw_combine[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA")),]
  
  
  list_filt=rpkm_filt(counts_raw_combine,rowData(data_raw)$gene_length[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))])
  
  log2rpkm_kp=list_filt$counts_log2rpkm_keep
  ls(list_filt)
  
  ### name the tissue-wise expressed gene exp matrix, as their ensemble id
  
  log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)#remove 4PC from the expressed genes in this tissue
  
  cor_ori=cor(t(log2rpkm_kp_rmPC))
  
  #### plot  - before using vst
  library(spqn)
  library(ggplot2)
  library(matrixStats)
  library(stats)
  
  # pdf("signal_condition_exp_rm4PCs.pdf")
  # pdf(paste0("signal_condition_exp_rm4PCs_poolTissues_",i,"_.pdf"))
  p=(plot_signal_condition_exp(cor_ori, rowMeans(log2rpkm_kp), signal=0.001))
  # dev.off()
  ggsave(p, 
         filename = paste0("signal_condition_exp_rm4PCs_poolTissues_",i,"_.pdf"),
         height = 2, width = 2, units = "in")
  

  
}




