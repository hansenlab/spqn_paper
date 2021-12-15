# dir_input="/users/ywang/10_25_2019/boxplot_9tissues/data/"

dir_output="/users/ywang/July_28_2020/pool_tissues/output_4pc/"
setwd(dir_output)
dir_input=dir_output

# load(paste0(dir_input,"list_IQR_est_",4,"_PCs_vst_",1,"th_tissue.RData"))
i=1
for(i in 1:3){
  
  
  # IQR - without vst
  # load(paste0(dir_output,"list_IQR_ori_",4,"_PCs_",i,"th_tissue.RData"))
  # save(list_IQR_ori,file=paste0(dir_output,"list_IQR_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))
  # IQR_grp_ori=list_IQR_ori$IQR_cor_mat_m
  

  # IQR - after combat
  load(paste0(dir_input,"list_IQR_ori_",i,"th_grp.RData"))
  # save(list_IQR_ori,file=paste0(dir_output,"list_IQR_ori_",i,"th_grp.RData"))
  IQR_grp_ori_combat=list_IQR_ori$IQR_cor_mat_m
  
  
  # get group mean
  grp_mean=list_IQR_ori$grp_mean
  
  mean_min=c()
 IQR_vec_ori_combat=c()
  for(j1 in 1:10){
    for(j2  in j1:10){
      mean_min=c(mean_min,grp_mean[j1])
      IQR_vec_ori_combat=c(IQR_vec_ori_combat,IQR_grp_ori_combat[j1,j2])
      # IQR_vec_ori=c(IQR_vec_ori,IQR_grp_ori[j1,j2])
      
      # IQR_vec_adj=c(IQR_vec_adj,IQR_grp_adj[i,j])
    }
  }
  min_IQR_combat=list(IQR_vec_ori_combat=IQR_vec_ori_combat,min_mean=mean_min)
  # min_IQR=list(IQR_vec_ori=IQR_vec_ori,min_mean=mean_min)
  
  save(min_IQR_combat,file=paste0(dir_output,"min_IQR_combat_",i,"th_grp.RData"))

}

