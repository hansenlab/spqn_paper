
########### # get the adjusted correlation for shared expressed genes only, for each tissue of gtex data
########### # by applying spqn on cor_ori that containing shared genes only
library(matrixStats)
library(WGCNA)
library(spqn)
library(SummarizedExperiment)
source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")

dir_data="/users/ywang/10_25_2019/data/"
dir_input="/users/ywang/July_28_2020/regulatome/output/"
dir_output="/users/ywang/July_28_2020/bias_regulatome/output/"
# dir.create(dir_output)

# load(paste0(dir_output,"TF_data_unique_en_filt.RData")) #TF_data_unique_en_filt

load(paste0(dir_input,"mat_reg_filt.RData")) #mat_reg_filt
diag(mat_reg_filt)=0
mat_reg_filt[lower.tri(mat_reg_filt)]=0
# mat_reg_filt[is.na(mat_TF_reg)]=0

n_PC=4
########### for each tissue, format the data 
list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
               "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")




i=1
for( i in 1:9){
  name_tissue=list_tissues[[i]]
  #dir_figure = paste0("/users/ywang/8_25_2019/figures/")
  
  #load and input raw counts, normalize and filter low expression genes(average(log2cpm)>0)
  name_var= paste0("gtex.",name_tissue)
  load(paste0(dir_data, name_var,".rda" ))
  data_raw=eval(as.name(name_var))
  counts_raw=assay(data_raw)
  counts_raw=counts_raw[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA")),]
  list_filt=rpkm_filt(counts_raw,rowData(data_raw)$gene_length[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))])
  format_geneID<-function(gene_id_in){
    gene_id_in=strsplit(as.vector(gene_id_in),".",fixed=TRUE)
    gene_id_in=matrix(unlist(gene_id_in), ncol = 2, byrow = TRUE)
    gene_id_in=gene_id_in[,1]
    gene_id_in
  }
  log2rpkm_kp=list_filt$counts_log2rpkm_keep
  
  gene_id_allgenes=rowData(data_raw)$gene_id
  gene_id_allgenes=format_geneID(gene_id_allgenes)
  gene_id_filt=gene_id_allgenes[list_filt$keep]
  rownames(log2rpkm_kp) = gene_id_filt
  
  
  
  # name the tissue-wise expressed gene exp matrix, as their ensemble id
  
  # log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,4)#remove 4PC from the expressed genes in this tissue
  
  reg_tissue=gene_id_filt[which(gene_id_filt %in% rownames(mat_reg_filt))]
  
  mat_reg_shared_tisse=mat_reg_filt[reg_tissue,reg_tissue]
  log2rpkm_kp=log2rpkm_kp[reg_tissue,]
  
  w=colSums(mat_reg_shared_tisse)
  w=w/sum(w)
  # w_adj=colSums(thisnet_adj)
  d_reg=density(rowMeans(log2rpkm_kp),weights = w)
  
  
  d_back=density(rowMeans(log2rpkm_kp))
  
  save(d_reg,file=paste0(dir_output,"d_reg_",n_PC,"_PCs_",i,"th_tissue.RData"))
  save(d_back,file=paste0(dir_output,"d_back_",n_PC,"_PCs_",i,"th_tissue.RData"))
  
  print(i)
}







