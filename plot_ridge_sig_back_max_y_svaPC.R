source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")
source("/users/ywang/10_25_2019/functions/functions_evaluation_Aug_5.R")
library(matrixStats)
#library(WGCNA)
library(ggplot2)
library(ggridges)


dir_data="/users/ywang/10_25_2019/data/"
dir_input="/users/ywang/10_25_2019/ridge_signal_0pc_ori/data/"
dir_fig="/users/ywang/ridge_signal_0pc_ori/fig/"





maxy_0PC_ori_back=maxy_4PC_ori_back=maxy_svaPC_ori_back=maxy_4PC_est_back=maxy_svaPC_est_back=c()
maxy_0PC_ori_sig=maxy_4PC_ori_sig=maxy_svaPC_ori_sig=maxy_4PC_est_sig=maxy_svaPC_est_sig=c()

get_maxy<-function(df,group = "background"){
  df <- df[df$group == group,]
  cor.split <- split(df$correlation, df$bin)
  dens.split <- lapply(cor.split, density, from = -1, to = 1)
  dens.split <- lapply(names(dens.split), function(nam) {
    dens  <- dens.split[[nam]]
    data.frame(x = dens$x, y = dens$y, bin = nam)
  })
  out <- cbind(do.call(rbind, dens.split), group = group)
  max(out$y)
}





for(i in 1:9){
  #nsva=list_nsv[i] 
  dir_figure = paste0("/users/ywang/10_25_2019/corplot_adj_4_svapcs/fig/")
  load(paste0(dir_figure,"nsv_",i,"_th_tissue.RData"))
  n_PC=nsv
  # load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_df_merge_ori.RData"))
  # maxy_svaPC_ori_back=max(maxy_svaPC_ori_back,get_maxy(df_merge))
  # maxy_svaPC_ori_sig=max(maxy_svaPC_ori_sig,get_maxy(df_merge, group = "signal"))
  
  load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_df_merge_est.RData"))
  maxy_svaPC_est_back=max(maxy_svaPC_est_back,get_maxy(df_merge))
  maxy_svaPC_est_sig=max(maxy_svaPC_est_sig,get_maxy(df_merge, group = "signal"))
}
save(maxy_svaPC_est_back,file=paste0(dir_input,"maxy_svaPC_est_back.RData"))
save(maxy_svaPC_est_sig,file=paste0(dir_input,"maxy_svaPC_est_sig.RData"))
# save(maxy_svaPC_ori_back,file=paste0(dir_input,"maxy_svaPC_ori_back.RData"))
# save(maxy_svaPC_ori_sig,file=paste0(dir_input,"maxy_svaPC_ori_sig.RData"))
