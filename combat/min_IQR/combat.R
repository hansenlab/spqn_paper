
########### # get the adjusted correlation for shared expressed genes only, for each tissue of gtex data
########### # by applying spqn on cor_ori that containing shared genes only
library(matrixStats)
library(WGCNA)
library(spqn)
library(SummarizedExperiment)
source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")
source("/users/ywang/July_28_2020/functions/plot_signal_condition_exp.R")

# source("/users/ywang/10_25_2019/functions/functions_evaluation_Aug_5.R")

# dir_data="/users/ywang/10_25_2019/data/"
dir_data="/dcl01/hansen/data/meanCoexp/"
# dir_input="/users/ywang/July_28_2020/regulatome/output/"
dir_output="/users/ywang/July_28_2020/combat/output_4pc/"
# dir.create(dir_output)
setwd(dir_output)
# load(paste0(dir_output,"TF_data_unique_en_filt.RData")) #TF_data_unique_en_filt

i=1


#### load sample metadata,
# variable description: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/dataset.cgi?study_id=phs000424.v5.p1&pht=2743

metadata=read.delim2("/users/ywang/July_28_2020/combat/data/GTEx_Data_V6_Annotations_SampleAttributesDS.txt")

sampleDate=metadata$SMNABTCHD
sampleDate=as.character(sampleDate)
names(sampleDate)=metadata$SAMPID

# 11983
# length(unique(sampleDate))
# 513

# 
# library(bladderbatch)
# data(bladderdata)
# dat <- bladderEset[1:50,]

# pheno = pData(dat)
# dim(pheno)
# edata = exprs(dat)
# edata[1:4,1:4]
# dim(edata)


  
  
  
  
  
  
  
  
#### 
i=1
for(i in c(1:9)){
  print('------------------------------------------')
  
  print(i)
  
  n_PC=4
  
  ########### for each tissue, format the data 
  list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
                 "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")
  name_tissue=list_tissues[[i]]
  
  #load and input raw counts, normalize and filter low expression genes(average(log2cpm)>0)
  name_var= paste0("gtex.",name_tissue)
  load(paste0(dir_data, name_var,".rda" ))
  data_raw=eval(as.name(name_var))
  counts_raw=assay(data_raw)
  counts_raw=counts_raw[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA")),]
  list_filt=rpkm_filt(counts_raw,rowData(data_raw)$gene_length[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))])
  
  log2rpkm_kp=list_filt$counts_log2rpkm_keep
  ls(list_filt)
  
  ### name the tissue-wise expressed gene exp matrix, as their ensemble id
  
  log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)#remove 4PC from the expressed genes in this tissue
  
  cor_ori=cor(t(log2rpkm_kp_rmPC))
  
  #### plot  - before using vst
  # library(spqn)
  library(ggplot2)
  library(matrixStats)
  library(stats)
  
  # pdf("signal_condition_exp_rm4PCs.pdf")
  # pdf(paste0("signal_condition_exp_rm4PCs_,",name_tissue,",.pdf"))
  p=(plot_signal_condition_exp(cor_ori, rowMeans(log2rpkm_kp), signal=0.001))
  # dev.off()
  # dim(cor_ori)
  # dim(log2rpkm_kp_rmPC)
  # dim(log2rpkm_kp)
  ggsave(p, 
         filename = paste0("signal_condition_exp_rm4PCs_,",name_tissue,".pdf"),
         height = 2, width = 2, units = "in")
  ### using combat 
  library(sva)

  # tb_=table(as.factor(sampleDate[colnames(log2rpkm_kp)]))
  # batch_kp=names(tb_)[which(tb_>1)]
  # which_sample_kp=which(sampleDate[colnames(log2rpkm_kp)]%in%batch_kp)
  # log2rpkm_kp_kp=log2rpkm_kp[,which_sample_kp]
  # dim(log2rpkm_kp_kp)
  # 
  # table(table(as.factor(sampleDate[colnames(log2rpkm_kp_kp)])))
  # 
  # dim(log2rpkm_kp)
  # 0   1   2   3   4   6   7
  # 419  63  19   9   1   1   1
  # class(sampleDate)
  # sampleDate=as.character(sampleDate)
  
  log2rpkm_kp_rmPC_combat=ComBat(log2rpkm_kp[,],batch=as.factor(sampleDate[colnames(log2rpkm_kp)]))
  
  cor_ori_combat=cor(t(log2rpkm_kp_rmPC_combat))
  
  # pdf(paste0("signal_condition_exp_rm4PCs_,",name_tissue,",_combat.pdf"))
  p=print(plot_signal_condition_exp(cor_ori_combat, rowMeans(log2rpkm_kp), signal=0.001))
  # dev.off()
  ggsave(p, 
         filename = paste0("signal_condition_exp_rm4PCs_,",name_tissue,",_combat.pdf"),
         height = 2, width = 2, units = "in")
  
  
  
}

  



