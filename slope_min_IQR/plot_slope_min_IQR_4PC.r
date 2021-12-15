library(matrixStats)

dir_data="/data/"

dir_input="/slope_min_QIR/data/"
dir_output="/slope_min_QIR/data/"

get_slope_min_IQR<-function(list_input){
  IQR_grp=list_input$IQR_cor_mat_m
  grp_mean=list_input$grp_mean
  mean_min=IQR_vec=c()
  for(k in 1:10){
    for(j in k:10){
      mean_min=c(mean_min,grp_mean[k])
      IQR_vec=c(IQR_vec,IQR_grp[k,j])
    }
  }
  lm(IQR_vec~mean_min)$coefficients[2]
}

for(i in 1:8){
  list_slope=c()
  
  load(paste0(dir_output,i,"th_tissue_list_IQR_mat.RData")) #list_IQR_mat
  for(nPC in 0:50){
    list_slope=c(list_slope,get_slope_min_IQR(list_IQR_mat[[nPC+1]]))
  }
  save(list_slope,file=paste0(dir_output,i,"th_tissue_list_slope.RData")) #list_slope
  
  print(i)
}




### on mac
dir_data="/Users/yiwang/Dropbox/project_Kasper/10_25_2019/cluster/slope_min_QIR/data/"
dir_fig="/Users/yiwang/Dropbox/project_Kasper/10_25_2019/cluster/slope_min_QIR/fig/"

range_y=c()
for(i in 1:9){
  load(paste0(dir_data,i,"th_tissue_list_slope.RData"))
  range_y=range(range_y,list_slope)
}
#[1] 0.001084519 0.043152610

for(i in 1:9){
  load(paste0(dir_data,i,"th_tissue_list_slope.RData"))
  pdf(paste0(dir_fig,i,"th_tissue_slope_min_IQR.pdf"))
  plot(c(0:50),list_slope,ylim=range_y,xaxt = "n", yaxt = "n",  bty = "n",xlab="number of PCs removed",ylab="coefficient",cex.lab=1.5,cex.axis=1.5, pch=1)
  axis(side =1, at = c(10,25,40),cex.axis=1.2)
  axis(side =2, at = c(0,0.02,0.04),cex.axis=1.2)
  dev.off()
}


for(i in 1:9){
  load(paste0(dir_data,i,"th_tissue_list_slope.RData"))
  pdf(paste0(dir_fig,i,"th_tissue_slope_min_IQR_cutoff.pdf"))
  list_slope[which(list_slope>0.02)]=0.02
  plot(c(0:50)[which(!list_slope>0.02)],list_slope[which(!list_slope>0.02)],ylim=c(0,0.02),xaxt = "n", yaxt = "n",cex.lab=1.5,cex.axis=1.5, pch=1,  bty = "n",xlab="number of PCs removed",ylab="coefficient")
  points(c(0:50)[which(list_slope==0.02)],list_slope[which(list_slope==0.02)],col="red")
  axis(side =1, at = c(10,25,40),cex.axis=1.2)
  axis(side =2, at = c(0,0.01,0.02),cex.axis=1.2)
  dev.off()
}


