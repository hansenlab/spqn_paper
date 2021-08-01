
cal_m_sd_IQR <- function(cor_matrix,ave_logcpm,ngrp=10){
  ngene = length(ave_logcpm)
  
  grp_label=cut(1:ngene,ngrp)
  grp_loc=split(1:ngene,grp_label)
  
  sd_cor_mat_m =IQR_cor_mat_m= array(dim=c(10,10))
  grp_mean=c()
  for(i in 1:10){
    for(j in 1:10){
      cor_tmp=cor_matrix[grp_loc[[i]],grp_loc[[j]]]
      if(i==j){
        sd_cor_mat_m[i,j]=sd(cor_tmp[upper.tri(cor_tmp)])#remove diagnal elements
        IQR_cor_mat_m[i,j]=IQR(cor_tmp[upper.tri(cor_tmp)])#remove diagnal elements
      }
      if(i!=j){sd_cor_mat_m[i,j]=sd(cor_tmp);IQR_cor_mat_m[i,j]=IQR(cor_tmp)}
    }
    grp_mean=c(grp_mean,mean(ave_logcpm[grp_loc[[i]]]))
  }
  return(list(sd_cor_mat_m=sd_cor_mat_m,IQR_cor_mat_m=IQR_cor_mat_m,grp_mean=grp_mean))
}

get_grp_loc<-function(cor_matrix){
  ngrp=10
  ngene=nrow(cor_matrix)
  grp_label=cut(1:ngene,ngrp)
  grp_loc=split(1:ngene,grp_label)
  return(grp_loc)
}

cal_plot_diagonal <- function(cor_matrix,ngrp=10){
  ngene = ncol(cor_matrix)
  grp_label=cut(1:ngene,ngrp)
  grp_loc=split(1:ngene,grp_label)
  
  colfunc <- colorRampPalette(c("gray77", "black"))
  color_10=colfunc(10)
  sd_cor_mat_m = array(dim=c(10,10))
  d=list()
  for(i in 1:10){
    for(j in 1:10){
      cor_tmp=cor_matrix[grp_loc[[i]],grp_loc[[j]]]
      if(i==j){
        d[[i]] = density(cor_tmp[upper.tri(cor_tmp)])
        sd_cor_mat_m[i,j]=sd(cor_tmp[upper.tri(cor_tmp)])#remove diagnal elements
      }
      if(i!=j){sd_cor_mat_m[i,j]=sd(cor_tmp)}
    }
  }
  yrange=range(d[[1]]$y, d[[2]]$y,d[[3]]$y,d[[4]]$y,d[[5]]$y,d[[6]]$y,d[[7]]$y,d[[8]]$y,d[[9]]$y,d[[10]]$y)
  plot(d[[1]],col=color_10[1],type="l",xlab="correlation",main="",xlim=c(-1,1),ylim=yrange)
  abline(v=0, col="black",lty=2)
  for(i in 2:10){
    lines(d[[i]],col=color_10[i])
  }
  return(list(sd_cor_mat_m=sd_cor_mat_m,d=d))
}


plot_diagonal_ridge <- function(cor_matrix,ngrp=10){
  ngene = ncol(cor_matrix)
  grp_label=cut(1:ngene,ngrp)
  grp_loc=split(1:ngene,grp_label) 
  for(i in 1:10){
    cor_tmp=cor_matrix[grp_loc[[i]],grp_loc[[i]]]
    cor_tmp=cor_tmp[upper.tri(cor_tmp)]
    if(i==1){cor_vec_all= cbind(as.numeric(cor_tmp),rep(i,length(as.numeric(cor_tmp))))
    names(cor_vec_all)=c("correlation","group")}else{
      cor_vec_all=rbind(cor_vec_all,cbind(as.numeric(cor_tmp),rep(i,length(as.numeric(cor_tmp)))))}
  }
  cor_vec_all=data.frame(cor_vec_all)
  names(cor_vec_all)=c("correlation","group")
  cor_vec_all$group=as.factor(cor_vec_all$group)
  ggplot(cor_vec_all, aes(x = correlation, y = group)) + geom_density_ridges2()+
    theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+ xlim(-1,1)+
    geom_vline(linetype="dotted",xintercept = 0,color = "black")
}


plot_diagonal_ridge2 <- function(cor_matrix){
  ngrp=10
  ngene = ncol(cor_matrix)
  grp_label=cut(1:ngene,ngrp)
  grp_loc=split(1:ngene,grp_label) 
  for(i in 1:10){
    cor_tmp=cor_matrix[grp_loc[[i]],grp_loc[[i]]]
    cor_tmp=cor_tmp[upper.tri(cor_tmp)]
    if(i==1){cor_vec_all= cbind(as.numeric(cor_tmp),rep(i,length(as.numeric(cor_tmp))))
    names(cor_vec_all)=c("correlation","group")}else{
      cor_vec_all=rbind(cor_vec_all,cbind(as.numeric(cor_tmp),rep(i,length(as.numeric(cor_tmp)))))}
  }
  cor_vec_all=data.frame(cor_vec_all)
  names(cor_vec_all)=c("correlation","group")
  cor_vec_all$group=as.factor(cor_vec_all$group)
  ggplot(cor_vec_all, aes(x = `correlation`, y = `group`, fill = ..x..)) +
    geom_density_ridges_gradient(scale = 3, gradient_lwd = 1.)+
    scale_fill_viridis(name = "correlation", option = "C")+
    theme_ridges( grid = TRUE) + theme(axis.title.x = element_blank())+ 
    geom_vline(xintercept = 0, linetype="dotted", color = "black", size=0.3)
}

plot_diagonal_ridge3 <- function(cor_matrix){
  ngrp=10
  ngene = ncol(cor_matrix)
  grp_label=cut(1:ngene,ngrp)
  grp_loc=split(1:ngene,grp_label) 
  for(i in 1:10){
    cor_tmp=cor_matrix[grp_loc[[i]],grp_loc[[i]]]
    cor_tmp=cor_tmp[upper.tri(cor_tmp)]
    if(i==1){cor_vec_all= cbind(as.numeric(cor_tmp),rep(i,length(as.numeric(cor_tmp))))
    names(cor_vec_all)=c("correlation","group")}else{
      cor_vec_all=rbind(cor_vec_all,cbind(as.numeric(cor_tmp),rep(i,length(as.numeric(cor_tmp)))))}
  }
  cor_vec_all=data.frame(cor_vec_all)
  names(cor_vec_all)=c("correlation","group")
  cor_vec_all$group=as.factor(cor_vec_all$group)
  ggplot(cor_vec_all, aes(x = `correlation`, y = `group`, fill = ..x..)) + xlim(-1,1)+
    geom_density_ridges_gradient(scale = 3, gradient_lwd = 1.)+
    scale_fill_viridis(name = "correlation", option = "C")+
    theme_ridges( grid = TRUE) + theme(axis.title.x = element_blank())+ 
    geom_vline(xintercept = 0, linetype="dotted", color = "black", size=0.3)
}

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
    xlab("average(log2(cpm+0.5))")+ylab("average(log2(cpm+0.5))")   +
    scale_y_discrete(limits=c(1:10)+c(0:9)*min(group_mat)/10,
                     labels=grp_medians)+
    scale_x_discrete(limits=c(1:10)+c(0:9)*min(group_mat)/10,
                     labels=grp_medians) + theme_bw()+
    theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text=element_text(size=20))
}


box_plot2<-function(group_mat,grp_medians){
  group_mat2=round(group_mat,3)
  group_mat=group_mat/max(group_mat)
  max_sd=max(group_mat)
  sd_grps_offset_y=t(matrix(rep(max_sd,100),nrow=10))
  for(i in 2:10){
    sd_grps_offset_y[i,]=rep(max_sd,10)+sd_grps_offset_y[i-1,]
  }
  
  sd_grps_offset_x=t(matrix(rep(max_sd,100),nrow=10))
  for(i in 2:10){
    sd_grps_offset_x[,i]=t(rep(max_sd,10)+sd_grps_offset_y[i-1,])
  }
  sd=as.numeric(group_mat)
  x1=as.numeric((-group_mat/2)+sd_grps_offset_x)
  x2=as.numeric((group_mat/2)+sd_grps_offset_x)
  y1=as.numeric((-group_mat/2)+sd_grps_offset_y)
  y2=as.numeric((group_mat/2)+sd_grps_offset_y)
  group_mat=round(group_mat,3)
  r=c(group_mat2[1,1],rep("",10),
      group_mat2[2,2],rep("",10),
      group_mat2[3,3],rep("",10),
      group_mat2[4,4],rep("",10),
      group_mat2[5,5],rep("",10),
      group_mat2[6,6],rep("",10),
      group_mat2[7,7],rep("",10),
      group_mat2[8,8],rep("",10),
      group_mat2[9,9],rep("",10),
      group_mat2[10,10])
  d=data.frame(x1,x2,y1,y2,r,sd)
  ggplot() + 
    geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="black", alpha=0) +
    geom_text(data=d, aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=r), size=4) +
    xlab("average(log2(cpm+0.5))")+ylab("average(log2(cpm+0.5))")   +
    scale_y_discrete(limits=1:10,
                     labels=grp_medians)+
    scale_x_discrete(limits=1:10,
                     labels=grp_medians) + theme_bw()+
    theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text=element_text(size=20))
}

box_plot3<-function(group_mat,grp_medians,title_=""){
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
    xlab("average(log2(cpm+0.5))")+ylab("average(log2(cpm+0.5))")   +
    scale_y_discrete(limits=c(1:10)+c(0:9)*min(group_mat)/10,
                     labels=grp_medians)+
    scale_x_discrete(limits=c(1:10)+c(0:9)*min(group_mat)/10,
                     labels=grp_medians) + theme_bw()+
    theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text=element_text(size=20))+
    ggtitle(title_)

}
