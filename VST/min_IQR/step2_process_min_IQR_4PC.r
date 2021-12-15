# dir_input="/users/ywang/10_25_2019/boxplot_9tissues/data/"

dir_output="/users/ywang/July_28_2020/VST/output_4pc/"
setwd(dir_output)
dir_input=dir_output

# load(paste0(dir_input,"list_IQR_est_",4,"_PCs_vst_",1,"th_tissue.RData"))
i=1
for(i in 1:9){
  
  
  # IQR - without vst

  
  load(paste0(dir_output,"list_IQR_ori_",4,"_PCs_",i,"th_tissue.RData"))
  # save(list_IQR_ori,file=paste0(dir_output,"list_IQR_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))
  IQR_grp_ori=list_IQR_ori$IQR_cor_mat_m
  

  # IQR - after vst
  load(paste0(dir_input,"list_IQR_ori_",4,"_PCs_vst_",i,"th_tissue.RData"))
  # save(list_IQR_ori_vst,file=paste0(dir_output,"list_IQR_ori_",n_PC,"_PCs_vst_",i,"th_tissue.RData"))
  IQR_grp_ori_vst=list_IQR_ori_vst$IQR_cor_mat_m
  
  
  # get group mean
  grp_mean=list_IQR_ori$grp_mean
  
  mean_min=c()
  IQR_vec_ori=IQR_vec_ori_vst=c()
  for(j1 in 1:10){
    for(j2  in j1:10){
      mean_min=c(mean_min,grp_mean[j1])
      IQR_vec_ori_vst=c(IQR_vec_ori_vst,IQR_grp_ori_vst[j1,j2])
      IQR_vec_ori=c(IQR_vec_ori,IQR_grp_ori[j1,j2])
      
      # IQR_vec_adj=c(IQR_vec_adj,IQR_grp_adj[i,j])
    }
  }
  min_IQR_vst=list(IQR_vec_ori=IQR_vec_ori_vst,min_mean=mean_min)
  min_IQR=list(IQR_vec_ori=IQR_vec_ori,min_mean=mean_min)
  
  save(min_IQR_vst,file=paste0(dir_output,"min_IQR_4PC_vst_",i,"th_tissue.RData"))
  save(min_IQR,file=paste0(dir_output,"min_IQR_4PC_",i,"th_tissue.RData"))
  
}

