library(readr)
library(matrixStats)
library(SummarizedExperiment)

original_files = "/data/dat_HURI_PPI_full.tsv"
dir_output="/PPI/data/"
source("/functions/functions_2d_quantile_no_surf_July_18.R")
source("/functions/functions_evaluation_Aug_5.R")
dir_data="/data/"

HURI=read_tsv(original_files)
#sum(HURI$Ensembl_gene_id_b=="ENSG00000000005")

#dim(HURI_coding)
#[1] 105002      2


#### map HURI genes to GTEX-including genes
# load gene names of GTEX data
load("/users/ywang/10_25_2019/PPI/data/gene_list_gtex.RData") #gene_list_gtex
list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
               "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")



for(i in 1:9){
  name_tissue=list_tissues[[i]]
  
  # load and input raw counts, normalize and filter low expression genes(average(log2cpm)>0)
  name_var= paste0("gtex.",name_tissue)
  load(paste0(dir_data, name_var,".rda" ))
  data_raw=eval(as.name(name_var))
  counts_raw=assay(data_raw)
  counts_raw=counts_raw[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA")),]
  list_filt=rpkm_filt(counts_raw,rowData(data_raw)$gene_length[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))])
  log2rpkm_kp=list_filt$counts_log2rpkm_keep
  
  genes_ensemble_keep2=gene_list_gtex$Name[list_filt$keep]
  genes_ensemble_keep2_2=strsplit(as.vector(genes_ensemble_keep2),".",fixed=TRUE)
  genes_ensemble_keep2_3=matrix(unlist(genes_ensemble_keep2_2), ncol = 2, byrow = TRUE)
  genes_ensemble_keep2_format=genes_ensemble_keep2_3[,1]


## map gene name to gene expression level for GTEX data

  ave_logrpkm_kp=rowMeans(list_filt$counts_log2rpkm_keep)
  exp_name=cbind(ave_logrpkm_kp,genes_ensemble_keep2_format)
  save(ave_logrpkm_kp,file=paste0(dir_output,i,"_th_tissue_ave_logrpkm_kp.RData" ))
  save(genes_ensemble_keep2_format,file=paste0(dir_output,i,"_th_tissue_genes_ensemble_keep2_format.RData" ))
  save(exp_name,file=paste0(dir_output,i,"_th_tissue_exp_name.RData" ))
}






 
