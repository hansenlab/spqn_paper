source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")
source("/users/ywang/10_25_2019/functions/functions_evaluation_Aug_5.R")
library(matrixStats)
library(WGCNA)
library(SummarizedExperiment)
library(ggplot2)
library(ggridges)

dir_data=paste0("/users/ywang/10_25_2019/scRNA/data/GSE45719_RAW/")
dir_input = paste0("/users/ywang/10_25_2019/scRNA/data/")
dir_output="/users/ywang/10_25_2019/scRNA/data/"
dir_fig="/users/ywang/10_25_2019/scRNA/fig/"



#########################################################################################################
# cor_ori scRNA  4 PCs
for(i in 1){
  # dir_sva = paste0("/users/ywang/10_25_2019/corplot_adj_4_svapcs/fig/")
  # load(paste0(dir_sva,"nsv_",i,"_th_tissue.RData"))
  # nsva=nsv

  for(n_PC in c(16)){
  #for(n_PC in c(4,nsva)){
    load(paste0(dir_output,"rpkm_midblast.RData"))#rpkm_midblast
    # dim(rpkm_midblast )
    # [1] 22958    60
    logrpkm_all=log2(rpkm_midblast+0.5)#not_sure
    keep=which(rowMedians(data.matrix(logrpkm_all))>0)
    log2rpkm_kp=logrpkm_all[keep,]

    if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
    corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))


    n_up=length(upper.tri(corr_ori))
    list_nsignal=c(round(n_up/10000),round(n_up/2000),round(n_up/1000),round(n_up/500),round(n_up/100),round(n_up/20))
    percent_list=c(1/10000*100,1/2000*100,1/1000*100,1/500*100,1/100*100,1/20*100)
    percent_list2=c("0_01","0_05", "0_1","0_2", "1", "5")


    grp_locs = get_grp_loc(corr_ori)
    for(k in 3){
      nsignal=list_nsignal[k]
      percent_=percent_list2[k]
      nsignal_grp=nsignal/100

      list_cor_sig=list_cor_back=grp_back=grp_sig=c()

      for(ngrp in 1:10){
        loc_back=grp_locs[[ngrp]]
        cor_back=corr_ori[loc_back,loc_back]
        cor_back=cor_back[upper.tri(cor_back)]
        order_cor_ori_grp=order(abs(cor_back),decreasing=TRUE)
        cor_signal_ori=cor_back[order_cor_ori_grp[1:nsignal_grp]]

        list_cor_sig=c(list_cor_sig,cor_signal_ori)
        list_cor_back=c(list_cor_back,cor_back)
        grp_back=c(grp_back,rep(ngrp,length(cor_back)))
        grp_sig=c(grp_sig,rep(ngrp,length(cor_signal_ori)))
      }

      df_sig=data.frame(correlation=list_cor_sig,bin=grp_sig,group="signal")
      df_back=data.frame(correlation=list_cor_back,bin=grp_back,group="background")

      df_merge=rbind(df_sig,df_back)
      save(df_merge,file=paste0(dir_output,percent_,"_percent_",n_PC,"PCs_",i,"_df_merge_ori.RData"))

    }
  }
  print(i)

}
# 
for(i in 1){
  # dir_sva = paste0("/users/ywang/10_25_2019/corplot_adj_4_svapcs/fig/")
  # load(paste0(dir_sva,"nsv_",i,"_th_tissue.RData"))
  # nsva=nsv
  
  for(n_PC in c(4,16)){
    #for(n_PC in c(4,nsva)){
    load(paste0(dir_output,"rpkm_midblast.RData"))#rpkm_midblast
    # dim(rpkm_midblast )
    # [1] 22958    60
    logrpkm_all=log2(rpkm_midblast+0.5)#not_sure
    keep=which(rowMedians(data.matrix(logrpkm_all))>0)
    log2rpkm_kp=logrpkm_all[keep,]
    
    if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
    #corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
    load(cor_est)
    cor_est
    
    n_up=length(upper.tri(corr_ori))
    list_nsignal=c(round(n_up/10000),round(n_up/2000),round(n_up/1000),round(n_up/500),round(n_up/100),round(n_up/20))
    percent_list=c(1/10000*100,1/2000*100,1/1000*100,1/500*100,1/100*100,1/20*100)
    percent_list2=c("0_01","0_05", "0_1","0_2", "1", "5")
    
    
    grp_locs = get_grp_loc(corr_ori)
    for(k in 3){
      nsignal=list_nsignal[k]
      percent_=percent_list2[k]
      nsignal_grp=nsignal/100
      
      list_cor_sig=list_cor_back=grp_back=grp_sig=c()
      
      for(ngrp in 1:10){
        loc_back=grp_locs[[ngrp]]
        cor_back=corr_ori[loc_back,loc_back]
        cor_back=cor_back[upper.tri(cor_back)]
        order_cor_ori_grp=order(abs(cor_back),decreasing=TRUE)
        cor_signal_ori=cor_back[order_cor_ori_grp[1:nsignal_grp]]
        
        list_cor_sig=c(list_cor_sig,cor_signal_ori)
        list_cor_back=c(list_cor_back,cor_back)
        grp_back=c(grp_back,rep(ngrp,length(cor_back)))
        grp_sig=c(grp_sig,rep(ngrp,length(cor_signal_ori)))
      }
      
      df_sig=data.frame(correlation=list_cor_sig,bin=grp_sig,group="signal")
      df_back=data.frame(correlation=list_cor_back,bin=grp_back,group="background")
      
      df_merge=rbind(df_sig,df_back)
      save(df_merge,file=paste0(dir_output,percent_,"_percent_",n_PC,"PCs_",i,"_df_merge_ori.RData"))
      
    }
  }
  print(i)
  
}


############################ plot for  2PCs
summarizeDensity <- function(df, group = "background", miny = 0.01) {
  df <- df[df$group == group,]
  cor.split <- split(df$correlation, df$bin)
  dens.split <- lapply(cor.split, density, from = -1, to = 1)
  dens.split <- lapply(names(dens.split), function(nam) {
    dens  <- dens.split[[nam]]
    data.frame(x = dens$x, y = dens$y, bin = nam)
  })
  out <- cbind(do.call(rbind, dens.split), group = group)
  out$y <- out$y/max(out$y)
  out <- out[out$y > miny,]
  out
}






for(i in 1){
  for(n_PC in 16){
  
  percent_="0_1"
  
  load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_df_merge_ori.RData"))
  
  pdf(paste0(dir_fig,percent_,"_percent_",i,"th_tissue",n_PC,"PCs_ridge_back_sig_ori.pdf"))
  df.plot <- rbind(summarizeDensity(df_merge),
                   summarizeDensity(df_merge, group = "signal"))
  print(
    
    ggplot(df.plot, aes(x = x, height = y, y = bin, fill = group)) +
      geom_ridgeline() +
      scale_y_discrete(expand = c(0.01, 0)) +
      scale_x_continuous(expand = c(0.01, 0)) +
      scale_fill_cyclical(
        #labels = c("background","signal"),
        values = c("#0000ff", "#ff0000"),
        name = "group", guide = "legend"
      )+  theme_ridges(grid = FALSE)+labs(
        x = "correlation",
        y = "bin" ))
  dev.off()
  }
}



# for(i in 1){
#   
#   n_PC=2
#   percent_="0_1"
#   
#   load(paste0(dir_input,percent_,"_percent_",n_PC,"PCs_",i,"_df_merge_est.RData"))
#   
#   pdf(paste0(dir_fig,percent_,"_percent_",i,"th_tissue",n_PC,"PCs_ridge_back_sig_est.pdf"))
#   df.plot <- rbind(summarizeDensity(df_merge),
#                    summarizeDensity(df_merge, group = "signal"))
#   print(
#     
#     ggplot(df.plot, aes(x = x, height = y, y = bin, fill = group)) +
#       geom_ridgeline() +
#       scale_y_discrete(expand = c(0.01, 0)) +
#       scale_x_continuous(expand = c(0.01, 0)) +
#       scale_fill_cyclical(
#         #labels = c("background","signal"),
#         values = c("#0000ff", "#ff0000"),
#         name = "group", guide = "legend"
#       )+  theme_ridges(grid = FALSE)+labs(
#         x = "correlation",
#         y = "bin" ))
#   dev.off()
#   
# }
