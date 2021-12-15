
dir_output="/drosophila/out/"
setwd(dir_output)
dir_input=dir_output

# load(paste0(dir_input,"list_IQR_est_",4,"_PCs_vst_",1,"th_tissue.RData"))
# i=1
for(n_PC in c(5:8)){
  
  
  # IQR - without vst
  # load(paste0(dir_output,"list_IQR_ori_",4,"_PCs_",i,"th_tissue.RData"))
  # save(list_IQR_ori,file=paste0(dir_output,"list_IQR_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))
  # IQR_grp_ori=list_IQR_ori$IQR_cor_mat_m
  

  # IQR - after combat
# save(list_IQR_ori,file=paste0(dir_output,"list_IQR_ori_combat_dro_rm",n_PC,"PCs.RData"))
# save(list_IQR_ori,file=paste0(dir_output,"list_IQR_ori_combat_dro.RData"))
load(paste0(dir_output,"list_IQR_ori_combat_dro_rm",n_PC,"PCs.RData"))
  IQR_grp_ori=list_IQR_ori$IQR_cor_mat_m
  # ls(IQR_grp_ori_combat)
  
  # get group mean
  grp_mean=list_IQR_ori$grp_mean
  
  mean_min=c()
  IQR_vec_ori=c()
  for(j1 in 1:10){
    for(j2  in j1:10){
      mean_min=c(mean_min,grp_mean[j1])
      IQR_vec_ori=c(IQR_vec_ori,IQR_grp_ori[j1,j2])
      # IQR_vec_ori=c(IQR_vec_ori,IQR_grp_ori[j1,j2])
      
      # IQR_vec_adj=c(IQR_vec_adj,IQR_grp_adj[i,j])
    }
  }
  min_IQR=list(IQR_vec_ori=IQR_vec_ori,min_mean=mean_min)
  # min_IQR=list(IQR_vec_ori=IQR_vec_ori,min_mean=mean_min)
  
  save(min_IQR,file=paste0(dir_output,"min_IQR_dro_rm_",n_PC,"PCs.RData"))

}

