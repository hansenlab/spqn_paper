library(pscl)
library(DESeq)
library(edgeR)
library(WGCNA)

source("//Users/yiwang/Dropbox/project_Kasper/ongoing/July_28_2020/functions/plot_signal_condition_exp.R")
source("//Users/yiwang/Dropbox/project_Kasper/ongoing/10_25_2019/cluster/functions/functions_evaluation_Aug_5.R")
library(spqn)

root="/Users/yiwang/Dropbox/project_Kasper/2_18_2019/"
# rpkm_filt<-function(counts_raw,gene_length){
#   counts_log2rpkm=log2(rpkm(counts_raw,gene_length))
#   keep=which(rowMedians(data.matrix(counts_log2rpkm))>0)
#   counts_log2rpkm_keep=counts_log2rpkm[keep,]
#   counts_raw_keep=counts_raw[keep,]
#   return(list(counts_raw_keep=counts_raw_keep,counts_log2rpkm_keep=counts_log2rpkm_keep,keep=keep))
# }

######### drosophila, large biological variance
#30 whole-animal biological samples. We discarded the larval, pupal and adult stages and kept only the 12 embryonic samples
# data downloaded from: https://github.com/TheCodingCollective/clusterFux/blob/master/data/modencodefly_pooledreps_count_table.txt
droso_counts_ <- read.table("/Users/yiwang/Dropbox/project_Kasper/ongoing/July_28_2020/drosophila/data//modencodefly_pooledreps_count_table.txt",header=TRUE)
# droso_counts_2 <- read.csv("/Users/yiwang/Dropbox/project_Kasper/July_28_2020/drosophila/data//normalized_count_file.csv",header=TRUE)
# dim(droso_counts_2)
# [1] 14869    31

# > dim(droso_counts_)
# [1] 14869    31


# droso_counts_[1:2,1]

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
n_PC=5
for(n_PC in c(5)){
  log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)#remove 4PC from the expressed genes in this tissue
  
  # cor_ori_noRmPc=cor(t(log2rpkm_kp))
  
  cor_ori=cor(t(log2rpkm_kp_rmPC))
  
  #### plot  - before using vst
  library(spqn)
  library(ggplot2)
  library(matrixStats)
  library(stats)

  # pdf("signal_condition_exp_rm4PCs.pdf")
  name_tissue="drosophila_development"
  # pdf(paste0("/Users/yiwang/Dropbox/project_Kasper/ongoing/July_28_2020//drosophila/out/signal_condition_exp_rm4PCs_,",name_tissue,".pdf"))
  p=(plot_signal_condition_exp(cor_ori, rowMeans(log2rpkm_kp), signal=0.001))
  # dev.off()
  ggsave(p, 
         filename = paste0("/Users/yiwang/Dropbox/project_Kasper/ongoing/July_28_2020//drosophila/out/signal_condition_exp_rm",n_PC,"PCs_,",name_tissue,".pdf"),
         height = 2, width = 2, units = "in")
  
}

# pdf(paste0("/Users/yiwang/Dropbox/project_Kasper/ongoing/July_28_2020/drosophila/out/signal_condition_exp_before_rmPCs_,",name_tissue,".pdf"))
# p=(plot_signal_condition_exp(cor_ori_noRmPc, rowMeans(log2rpkm_kp), signal=0))
# # dev.off()
# ggsave(p, 
#        filename = paste0("/Users/yiwang/Dropbox/project_Kasper/ongoing/July_28_2020/drosophila/out/signal_condition_exp_before_rmPCs_,",name_tissue,".pdf"),
#        height = 2, width = 2, units = "in")
# 

# dim(cor_ori)
# dim(log2rpkm_kp_rmPC)
# dim(log2rpkm_kp)





