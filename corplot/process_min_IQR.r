dir_input="/users/ywang/10_25_2019/boxplot_9tissues/data/"
dir_output="/users/ywang/10_25_2019/corplot_adj_4_svapcs/data/"

load(paste0(dir_input,"list_IQR_est_",30,"_PCs_",1,"th_tissue.RData"))
load(paste0(dir_input,"list_IQR_ori_",30,"_PCs_",1,"th_tissue.RData"))
IQR_grp_ori=list_IQR_ori$IQR_cor_mat_m
IQR_grp_adj=list_IQR_est$IQR_cor_mat_m
grp_mean=list_IQR_ori$grp_mean

mean_min=c()
IQR_vec_ori=IQR_vec_adj=c()
for(i in 1:10){
  for(j  in i:10){
    mean_min=c(mean_min,grp_mean[i])
    IQR_vec_ori=c(IQR_vec_ori,IQR_grp_ori[i,j])
    IQR_vec_adj=c(IQR_vec_adj,IQR_grp_adj[i,j])
  }
}
min_IQR=list(IQR_vec_ori=IQR_vec_ori,IQR_vec_adj=IQR_vec_adj,min_mean=mean_min)

save(min_IQR,file=paste0(dir_output,"min_IQR.RData"))

