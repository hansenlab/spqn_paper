source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")
source("/users/ywang/10_25_2019/functions/functions_evaluation_Aug_5.R")
library(matrixStats)
#library(WGCNA)
library(ggplot2)
library(ggridges)


dir_input="/users/ywang/10_25_2019/ridge_signal_0pc_ori/data/"
dir_fig="/users/ywang/10_25_2019/ridge_signal_0pc_ori/fig/"


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


percent_="0_1"


for(i in 1:9){
  #nsva=list_nsv[i] 
  n_PC=4
  # load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_df_merge_ori.RData"))
  # maxy_4PC_ori_back=max(maxy_4PC_ori_back,get_maxy(df_merge))
  # maxy_4PC_ori_sig=max(maxy_4PC_ori_sig,get_maxy(df_merge, group = "signal"))
  load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_df_merge_est.RData"))
  
  maxy_4PC_est_back=max(maxy_4PC_est_back,get_maxy(df_merge))
  maxy_4PC_est_sig=max(maxy_4PC_est_sig,get_maxy(df_merge, group = "signal"))
}
save(maxy_4PC_est_back,file=paste0(dir_input,"maxy_4PC_est_back.RData"))
save(maxy_4PC_est_sig,file=paste0(dir_input,"maxy_4PC_est_sig.RData"))
# save(maxy_4PC_ori_back,file=paste0(dir_input,"maxy_4PC_ori_back.RData"))
# save(maxy_4PC_ori_sig,file=paste0(dir_input,"maxy_4PC_ori_sig.RData"))

