#library("sva")
#library("recount", quietly = T)
library(matrixStats)
library(spqn)

# get the rank of abs(correlation) for shared expressed genes only, for each tissue of gtex data
library(matrixStats)
library(WGCNA)
library(SummarizedExperiment)
source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")

dir_data="/users/ywang/10_25_2019/data/"
dir_output="/users/ywang/July_28_2020/graphical_lasso_morepoints/output_27pc/"
# dir_output="/users/ywang/July_28_2020/graphical_lasso/output_27pc/"

dir_fig="/users/ywang/July_28_2020/graphical_lasso_morepoints/fig_27PC/"
dir.create(dir_fig)

lambda = c(5:20)/100*4

# i=1
# j=1


list_file=list.files(dir_output)
#for i==6:
# setwd("/users/ywang/July_28_2020/graphical_lasso/output_27pc/")
dir_output2="/users/ywang/July_28_2020/graphical_lasso/output_27pc/"

for(i in  6){
  n_PC=27
  
  list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
                 "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")
  
  name_tissue=list_tissues[[i]]
  name_var= paste0("gtex.",name_tissue)
  load(paste0(dir_data, name_var,".rda" ))
  data_raw=eval(as.name(name_var))
  counts_raw=assay(data_raw)
  counts_raw=counts_raw[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA")),]
  list_filt=rpkm_filt(counts_raw,rowData(data_raw)$gene_length[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))])
  
  log2rpkm_kp=list_filt$counts_log2rpkm_keep
  
  set.seed(111)
  kep=sample(1:nrow(log2rpkm_kp),4000)
  
  log2rpkm_kp=log2rpkm_kp[kep,]
  
  
  d_ori=d_adj=list()
  nedge_adj=nedge_ori=mean_adj=mean_ori=mean_background=array(dim=length(lambda))
  
  for(j in 1:length(lambda)){
    
    if(paste0("thisnet_ori_",i,"th_tissue_",j,"_4000genes.RData") %in% list_file){
      # save(thisnet,file=paste0(dir_output,"thisnet_adj_",i,"th_tissue_",j,"_4000genes.RData"))
      
      load(paste0(dir_output2,"thisnet_ori_",i,"th_tissue_",j,"_4000genes.RData")) #thisnet,ori
      thisnet_ori=thisnet
      load(paste0(dir_output2,"thisnet_adj_",i,"th_tissue_",j,"_4000genes.RData")) #thisnet,adj
      thisnet_adj=thisnet
      
      
      thisnet_ori[upper.tri(thisnet_ori)]=0
      diag(thisnet_ori)=0
      thisnet_adj[upper.tri(thisnet_adj)]=0
      diag(thisnet_adj)=0
      
      w_ori=colSums(thisnet_ori)
      w_adj=colSums(thisnet_adj)
      
      nedge_ori[j]=sum(w_ori)
      nedge_adj[j]=sum(w_adj)
      
      
      if(sum(w_ori)>0){
        w_ori=w_ori/sum(w_ori)
        d_ori[[j]]=density(rowMeans(log2rpkm_kp),weights = w_ori)
        mean_ori[j]=weighted.mean(rowMeans(log2rpkm_kp),w = w_ori)
      }
    }
    
    
    # 
    if(sum(w_adj)>0){
      w_adj=w_adj/sum(w_adj)
      d_adj[[j]]=density(rowMeans(log2rpkm_kp),weights = w_adj)
      mean_adj[j]=weighted.mean(rowMeans(log2rpkm_kp),w = w_adj)
    }
    
    print(j)
    
  }
  # j=1
  mean_background=mean(rowMeans(log2rpkm_kp))
  save(mean_background,file=paste0(dir_output,"mean_background_",n_PC,"_PCs_",i,"th_tissue.RData"))
  save(mean_ori,file=paste0(dir_output,"mean_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))
  save(mean_adj,file=paste0(dir_output,"mean_adj_",n_PC,"_PCs_",i,"th_tissue.RData"))
  save(nedge_ori,file=paste0(dir_output,"nedge_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))
  save(nedge_adj,file=paste0(dir_output,"nedge_adj_",n_PC,"_PCs_",i,"th_tissue.RData"))
  save(d_ori,file=paste0(dir_output,"d_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))
  save(d_adj,file=paste0(dir_output,"d_adj_",n_PC,"_PCs_",i,"th_tissue.RData"))
  
  
  for(j in 1:length(lambda)){
    if(j <= length(d_ori)){
      pdf(paste0(dir_fig,"d_ori_",j,"_",i,"th_tissue.pdf"))
      d_ori_tmp=d_ori[[j]]
      d_background=density(rowMeans(log2rpkm_kp))
      xlim_=range(d_background$x,d_ori_tmp$x)
      ylim_=range(d_background$y,d_ori_tmp$y)
      par(mar =c(7, 5, 7, 5))
      main_=paste0("before adjustment, rho=",lambda[j],", nedge=",formatC(nedge_ori[[j]],format="e"))
      plot(d_background,xlim=xlim_,ylim=ylim_,xlab="average log2rpkm",cex.lab=1.5,main=main_,col="black")
      lines(d_ori_tmp,col="blue")
      legend("topright",legend=c("background","signal"),col=c("black","blue"),lty=1,cex=1.5)
      dev.off()
    }
    
    
    if(j <= length(d_adj)){
      pdf(paste0(dir_fig,"/d_adj_",j,"_",i,"th_tissue.pdf"))
      d_adj_tmp=d_adj[[j]]
      # d_background=density(rowMeans(log2rpkm_kp))
      xlim_=range(d_background$x,d_adj_tmp$x)
      ylim_=range(d_background$y,d_adj_tmp$y)
      par(mar =c(7, 5, 7, 5))
      main_=paste0("after adjustment, rho=",lambda[j],", nedge=",formatC(nedge_adj[[j]],format="e"))
      plot(d_background,xlim=xlim_,ylim=ylim_,xlab="average log2rpkm",cex.lab=1.5,main=main_,col="black")
      lines(d_adj_tmp,col="blue")
      legend("topright",legend=c("background","signal"),col=c("black","blue"),lty=1,cex=1.5)
      dev.off()
    }
    
    
  }
  
}



