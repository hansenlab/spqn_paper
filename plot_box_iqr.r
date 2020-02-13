library(ggplot2)


dir_input="/users/ywang/10_25_2019/boxplot_9tissues/data/"
dir_fig="/users/ywang/10_25_2019/corplot_adj_4_svapcs/fig/"

i=1
load(paste0("/users/ywang/10_25_2019/corplot_adj_4_svapcs/fig/","nsv_",i,"_th_tissue.RData"))
n_PC=nsv

load(paste0(dir_input,"list_IQR_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))
load(paste0(dir_input,"list_IQR_est_",n_PC,"_PCs_",i,"th_tissue.RData"))

box_plot<-function(group_mat,grp_medians){
  group_mat2=round(group_mat,3)
  group_mat=group_mat/max(group_mat)
  max_sd=max(group_mat)
  sd_grps_offset_y=t(matrix(rep(max_sd,100),nrow=10))
  for(i in 2:10){
    sd_grps_offset_y[i,]=rep(max_sd,10)+sd_grps_offset_y[i-1,]+min(group_mat)/10
  }
  
  sd_grps_offset_x=t(matrix(rep(max_sd,100),nrow=10))
  for(i in 2:10){
    sd_grps_offset_x[,i]=t(rep(max_sd,10)+sd_grps_offset_y[i-1,])+min(group_mat)/10
  }
  sd=as.numeric(group_mat)
  x1=as.numeric((-group_mat/2)+sd_grps_offset_x)#+c(0:9)*min(group_mat)/10
  x2=as.numeric((group_mat/2)+sd_grps_offset_x)
  y1=as.numeric((-group_mat/2)+sd_grps_offset_y)
  y2=as.numeric((group_mat/2)+sd_grps_offset_y)
  
  group_mat=round(group_mat,3)
  d=data.frame(x1,x2,y1,y2,sd)
  ggplot() + 
    geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="black", alpha=0) +
    xlab("average(log2RPKM)")+ylab("average(log2RPKM)")   +
    scale_y_discrete(limits=c(1:10)+c(0:9)*min(group_mat)/10,
                     labels=grp_medians)+
    scale_x_discrete(limits=c(1:10)+c(0:9)*min(group_mat)/10,
                     labels=grp_medians) + theme_bw()+
    theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text=element_text(size=20))
}



pdf(paste0(dir_fig,"box_iqr_ori_svaPC.pdf"))
box_plot(list_IQR_ori$IQR_cor_mat_m,round(list_IQR_ori$grp_mean,1))
dev.off()

pdf(paste0(dir_fig,"box_iqr_adj_svaPC.pdf"))
box_plot(list_IQR_est$IQR_cor_mat_m,round(list_IQR_est$grp_mean,1))
dev.off()




