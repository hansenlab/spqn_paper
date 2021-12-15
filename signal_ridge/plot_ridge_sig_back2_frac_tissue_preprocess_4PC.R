source("/functions/functions_2d_quantile_no_surf_July_18.R")
source("/functions/functions_evaluation_Aug_5.R")
library(matrixStats)
#library(WGCNA)
library(ggplot2)
library(ggridges)


dir_data="/data/"
dir_input="/ridge_signal_0pc_ori/data/"
dir_fig="/ridge_signal_0pc_ori/fig/"
dir.create(dir_fig)

list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
               "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")
tissue_list=c( "adiposesub", "adrenalgl","arteryti","braince","brainco","breastma",
               "colontr","esophagusmu","heartleft")
#list_nsv=c(27,15,25,16,12,12,9,21,13)


load(paste0(dir_input,"maxy_4PC_ori_back.RData"))
load(paste0(dir_input,"maxy_4PC_ori_sig.RData"))


load(paste0(dir_input,"maxy_4PC_est_back.RData"))
load(paste0(dir_input,"maxy_4PC_est_sig.RData"))
load(paste0(dir_input,"maxy_svaPC_est_back.RData"))
load(paste0(dir_input,"maxy_svaPC_est_sig.RData"))

summarizeDensity2 <- function(df, maxy,group = "background", miny = 0.01) {
  df <- df[df$group == group,]
  cor.split <- split(df$correlation, df$bin)
  dens.split <- lapply(cor.split, density, from = -1, to = 1)
  dens.split <- lapply(names(dens.split), function(nam) {
    dens  <- dens.split[[nam]]
    data.frame(x = dens$x, y = dens$y, bin = nam)
  })
  out <- cbind(do.call(rbind, dens.split), group = group)
  #out$y <- out$y/max(out$y)
  out$y <- out$y/maxy
  
  out <- out[out$y > miny,]
  out
}






for(i in 1:9){
  #nsva=list_nsv[i]

  percent_="0_1"  
  n_PC=4 
  # load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_df_merge_ori.RData"))#df_merge
  # df_plot <- rbind(summarizeDensity2(df_merge, maxy=maxy_4PC_ori_back),
  #                  summarizeDensity2(df_merge, maxy=maxy_4PC_ori_sig,group = "signal"))
  # df_plot$tissue=list_tissues[i]
  # save(df_plot,file=paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_df_plot_ori.RData"))
  # 
  load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_df_merge_est.RData"))#df_merge
  df_plot <- rbind(summarizeDensity2(df_merge, maxy=maxy_4PC_est_back),
                   summarizeDensity2(df_merge, maxy=maxy_4PC_est_sig,group = "signal"))
  df_plot$tissue=list_tissues[i]
  save(df_plot,file=paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_df_plot_est.RData"))
  print(i)
}


