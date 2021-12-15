library(readr)
library(matrixStats)

original_files = "/PPI/data/dat_HURI_PPI_full.tsv"
dir_output="/PPI/data/"
source("/functions/functions_2d_quantile_no_surf_July_18.R")
source("/functions/functions_evaluation_Aug_5.R")

dir_data="/data/"

HURI=read_tsv(original_files)

## keep only genes overlapped with GTEX data in PPI
for(i in 1:9){
  load(paste0(dir_output,i,"_th_tissue_ave_logrpkm_kp.RData" ))
  load(paste0(dir_output,i,"_th_tissue_genes_ensemble_keep2_format.RData" ))
  PPI_overlap=HURI[which(HURI$gene1 %in% genes_ensemble_keep2_format & HURI$gene2%in% genes_ensemble_keep2_format),]
  ## match expression level to the pairs of PPI data
  names(ave_logrpkm_kp)=genes_ensemble_keep2_format
  exp_level_gene1_PPI= ave_logrpkm_kp[PPI_overlap$gene1]
  exp_level_gene2_PPI= ave_logrpkm_kp[PPI_overlap$gene2]
  ave_exp_PPI=(exp_level_gene1_PPI+exp_level_gene2_PPI)/2
  save(ave_exp_PPI,file=paste0(dir_output,i,"_th_tissue_ave_exp_PPI.RData" ))
  dist_ave_exp_PPI=density(ave_exp_PPI)
  save(dist_ave_exp_PPI,file=paste0(dir_output,i,"_th_tissue_dist_ave_exp_PPI.RData" ))
}

### both
for(i in 1:9){
  load(paste0(dir_output,i,"_th_tissue_ave_logrpkm_kp.RData" ))
  load(paste0(dir_output,i,"_th_tissue_genes_ensemble_keep2_format.RData" ))
  PPI_overlap=HURI[which(HURI$gene1 %in% genes_ensemble_keep2_format & HURI$gene2%in% genes_ensemble_keep2_format),]
  ## match expression level to the pairs of PPI data
  names(ave_logrpkm_kp)=genes_ensemble_keep2_format
  exp_level_gene1_PPI= ave_logrpkm_kp[PPI_overlap$gene1]
  exp_level_gene2_PPI= ave_logrpkm_kp[PPI_overlap$gene2]
  ave_exp_PPI=c(exp_level_gene1_PPI,exp_level_gene2_PPI)
  save(ave_exp_PPI,file=paste0(dir_output,i,"_th_tissue_both_exp_PPI.RData" ))
  dist_ave_exp_PPI=density(ave_exp_PPI)
  save(dist_ave_exp_PPI,file=paste0(dir_output,i,"_th_tissue_dist_both_exp_PPI.RData" ))#dist_ave_exp_PPI
}



# 
# 
# 
# 
# 
# 
