library(matrixStats)
library(WGCNA)
library(ggplot2)
library(ggridges)
library(SummarizedExperiment)
library("sva")

source("/functions/functions_2d_quantile_no_surf_July_18.R")
source("/functions/functions_evaluation_Aug_5.R")

dir_data=paste0("/scRNA/data/GSE45719_RAW/")
dir_inpu = paste0("/scRNA/data/GSE45719_RAW/")
dir_output="/scRNA/data/"
dir_figure="/scRNA/fig/"



n_PC=4
######load data and filtration
# load counts data
i=1
load(paste0(dir_output,"rpkm_midblast.RData"))#rpkm_midblast
# dim(rpkm_midblast )
# [1] 22958    60
logrpkm_all=log2(rpkm_midblast+0.5)#not_sure
keep=which(rowMedians(data.matrix(logrpkm_all))>0)
log2rpkm_kp=logrpkm_all[keep,]


mod=matrix(1,nrow=dim(log2rpkm_kp)[2],ncol=1)
colnames(mod)="Intercept"
nsv=num.sv(log2rpkm_kp,mod, method = "be") ## num.sv requires data matrix with features(genes) in the rows and samples in the column
print(paste("Number of PCs estimated to be removed:", nsv))
save(nsv,file=paste0(dir_figure,"nsv_",i,"_th_tissue.RData"))
print(nsv)
