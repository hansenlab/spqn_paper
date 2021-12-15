library(pscl)
library(DESeq)
library(edgeR)
library(WGCNA)

# root="/Users/yiwang/Dropbox/project_Kasper/2_18_2019/"
dir_output="/Users/yiwang/Dropbox/project_Kasper/ongoing/July_28_2020/drosophila/out/"

######### drosophila, large biological variance
#30 whole-animal biological samples. We discarded the larval, pupal and adult stages and kept only the 12 embryonic samples
# data downloaded from: https://github.com/TheCodingCollective/clusterFux/blob/master/data/modencodefly_pooledreps_count_table.txt
droso_counts_ <- read.table("/Users/yiwang/Dropbox/project_Kasper/July_28_2020/drosophila/data//modencodefly_pooledreps_count_table.txt",header=TRUE)

droso_counts=droso_counts_[,-1]

head(droso_counts)

# counts_raw=counts_raw[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA")),]
source('/Users/yiwang/Dropbox/project_Kasper/ongoing/10_25_2019/cluster/functions/functions_2d_quantile_no_surf_July_18.R')
list_filt=rpkm_filt(droso_counts,rep(1000,nrow(droso_counts)))

log2rpkm_kp=list_filt$counts_log2rpkm_keep
ls(list_filt)


dim(log2rpkm_kp)
# [1] 9719   30

### name the tissue-wise expressed gene exp matrix, as their ensemble id
n_PC=4
for(n_PC in c(5:8)){
  log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)#remove 4PC from the expressed genes in this tissue
  
  # cor_ori_noRmPc=cor(t(log2rpkm_kp))
  
  # cor_ori=cor(t(log2rpkm_kp_rmPC))
  
  corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
  ave_logrpkm_kp=rowMeans(log2rpkm_kp)[order(rowMeans(log2rpkm_kp))]
  
  
  
  list_IQR_ori = cal_m_sd_IQR(corr_ori,ave_logrpkm_kp)
  save(list_IQR_ori,file=paste0(dir_output,"list_IQR_ori_combat_dro_rm",n_PC,"PCs.RData"))
  
}






