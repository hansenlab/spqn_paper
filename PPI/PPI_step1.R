library(readr)
library(matrixStats)
library(SummarizedExperiment)

#original_files = "/Users/yiwang/Dropbox/project_Kasper/4_26_2019/threshold_PPI/data/dat_HURI_PPI_full.tsv"
original_files = "/users/ywang/10_25_2019/PPI/data/dat_HURI_PPI_full.tsv"
dir_output="/users/ywang/10_25_2019/PPI/data/"
source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")
source("/users/ywang/10_25_2019/functions/functions_evaluation_Aug_5.R")
dir_data="/users/ywang/10_25_2019/data/"

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







# 
# 
# 
# 
# ## filter out noncoding genes from PPI data
# load("/users/ywang/8_25_2019/overleaf_code_data/PPI/data/humanProteinCodingGenes_ensembl.RData")#humanProteinCodingGenes_ensembl
# PPI_coding=HURI[which(HURI$gene1 %in% humanProteinCodingGenes_ensembl & HURI$gene2%in% humanProteinCodingGenes_ensembl),]
# save(PPI_coding,file=paste0(dir_output,"PPI_coding.RData" ))
# 
# 
# ## keep only genes overlapped with GTEX data in PPI 
# PPI_coding_overlap=PPI_coding[which(PPI_coding$gene1 %in% genes_ensemble_keep2_format & PPI_coding$gene2%in% genes_ensemble_keep2_format),]
# ## match expression level to the pairs of PPI data
# rownames(ave_logcpm_kp)=genes_ensemble_keep2_format
# exp_level_gene1_PPI= ave_logcpm_kp[PPI_coding_overlap$gene1]
# exp_level_gene2_PPI= ave_logcpm_kp[PPI_coding_overlap$gene2]
# ave_exp_PPI=(exp_level_gene1_PPI+exp_level_gene2_PPI)/2
# save(ave_exp_PPI,file=paste0(dir_output,i,"_th_tissue_ave_exp_PPI.RData" ))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
