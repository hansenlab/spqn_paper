
library(matrixStats)
library(spqn)

# get the rank of abs(correlation) for shared expressed genes only, for each tissue of gtex data
library(matrixStats)
library(WGCNA)
library(SummarizedExperiment)
source("/functions/functions_2d_quantile_no_surf_July_18.R")

dir_data="/data/"
dir_output="/graphical_lasso_morepoints/output_4pc/"
# dir_output="/graphical_lasso/output_27pc/"

dir_fig="/graphical_lasso_morepoints/fig_4PC/"
dir.create(dir_fig)

lambda = c(5:20)/100*4

# i=1
# j=1

n_PC=4

list_file=list.files(dir_output)

for(i in  1:9){
  list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
                 "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")
  
  # for(i in 1:9){
  name_tissue=list_tissues[[i]]
  # j=1
  # mean_background=mean(rowMeans(log2rpkm_kp))
  load(paste0(dir_output,"mean_background_",n_PC,"_PCs_",i,"th_tissue.RData"))#mean_background
  load(paste0(dir_output,"mean_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))#mean_ori
  load(paste0(dir_output,"mean_adj_",n_PC,"_PCs_",i,"th_tissue.RData"))#mean_adj
  load(paste0(dir_output,"nedge_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))#nedge_ori
  load(paste0(dir_output,"nedge_adj_",n_PC,"_PCs_",i,"th_tissue.RData"))#nedge_adj
  # save(d_ori,file=paste0(dir_output,"d_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))
  # save(d_adj,file=paste0(dir_output,"d_adj_",n_PC,"_PCs_",i,"th_tissue.RData"))
  
  pdf(paste0(dir_fig,i,"th_tissue_4PC_glasso_bias.pdf"))
  plot(log(nedge_ori),mean_ori-mean_background,ylab="blas",main=name_tissue,xlab="log(nedge)",
       ylim=range(mean_ori[!is.na(mean_ori)]-mean_background,mean_adj[!is.na(mean_adj)]-mean_background),
       xlim=range(log(1+nedge_ori[!is.na(nedge_ori)]),log(1+nedge_adj[!is.na(nedge_adj)])))
  points(log(nedge_adj),mean_adj-mean_background,col="blue")
  legend("bottomright",legend = c("before adjustment","after adjustment"),col = c("black","blue"),pch = 1)
  abline(0,0,col="grey")
  dev.off()
  
}

ntotal=4000*(4000-1)/2

i=1
list_tissues_name=c("Adipose Subcutaneous","Adrenal Gland","Artery Tibial","Brain Cerebellum",
               "Brain Cortex", "Breast Mammary","Colon Transverse","Esophagus Mucosa","Heart Left Ventricle")

name_tissue=list_tissues_name[[i]]

load(paste0(dir_output,"mean_background_",n_PC,"_PCs_",i,"th_tissue.RData"))#mean_background
load(paste0(dir_output,"mean_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))#mean_ori
load(paste0(dir_output,"mean_adj_",n_PC,"_PCs_",i,"th_tissue.RData"))#mean_adj
load(paste0(dir_output,"nedge_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))#nedge_ori
load(paste0(dir_output,"nedge_adj_",n_PC,"_PCs_",i,"th_tissue.RData"))#nedge_adj
df1=data.frame(mean_est=mean_ori,mean_background=mean_background,signal_percent=(nedge_ori/ntotal*100),log_nedge=log(nedge_ori+1),type="before adjustment",tissue=name_tissue)
df2=data.frame(mean_est=mean_adj,mean_background=mean_background,signal_percent=(nedge_adj/ntotal*100),log_nedge=log(nedge_adj+1),type="after adjustment",tissue=name_tissue)
df=rbind(df1,df2)

for(i in 2:9){
  name_tissue=list_tissues_name[[i]]
  
  load(paste0(dir_output,"mean_background_",n_PC,"_PCs_",i,"th_tissue.RData"))#mean_background
  load(paste0(dir_output,"mean_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))#mean_ori
  load(paste0(dir_output,"mean_adj_",n_PC,"_PCs_",i,"th_tissue.RData"))#mean_adj
  load(paste0(dir_output,"nedge_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))#nedge_ori
  load(paste0(dir_output,"nedge_adj_",n_PC,"_PCs_",i,"th_tissue.RData"))#nedge_adj
  df1=data.frame(mean_est=mean_ori,mean_background=mean_background,signal_percent=(nedge_ori/ntotal*100),log_nedge=log(nedge_ori+1),type="before adjustment",tissue=name_tissue)
  df2=data.frame(mean_est=mean_adj,mean_background=mean_background,signal_percent=(nedge_adj/ntotal*100),log_nedge=log(nedge_adj+1),type="after adjustment",tissue=name_tissue)
  df=rbind(df,df1)
  df=rbind(df,df2)
}
# df_=df
save(df,file=paste0(dir_output,"df_4pc_summary_tissues.RData"))
df=df[which(df$signal_percent>0.05),]

pdf(paste0(dir_fig,"bias_glasso_4pc.pdf"))
plot(df$log_nedge[which(df$type=="ori")],df$bias[which(df$type=="ori")],main="all tissues together",xlab="log(nedge)",ylab="bias",
     xlim=range(df$log_nedge[!is.na(df$log_nedge)]),ylim=range(df$bias[!is.na(df$bias)]))
points(df$log_nedge[which(df$type=="adj")],df$bias[which(df$type=="adj")],col="blue")
legend("topright",legend = c("before adjustment","after adjustment"),col = c("black","blue"),pch = 1)
dev.off()





####################fracet plot
library(ggplot2)
library(matrixStats)
library(spqn)

# get the rank of abs(correlation) for shared expressed genes only, for each tissue of gtex data
library(matrixStats)
library(WGCNA)
library(SummarizedExperiment)
source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")

dir_data="/users/ywang/10_25_2019/data/"
# dir_output="/users/ywang/July_28_2020/graphical_lasso/output_4pc/"

# dir_fig="/users/ywang/July_28_2020/graphical_lasso/fig_4PC/"

# lambda = seq(0.3,1.0,length.out=40)

n_PC=4

list_file=list.files(dir_output)

load(paste0(dir_output,"df_4pc_summary_tissues.RData")) #df

# > ls(df)
# [1] "bias"      "log_nedge" "tissue"    "type"  

df=df[which(!is.na(df$mean_est)),]

df$type=factor(df$type,levels=c("before adjustment","after adjustment"))
  
df=df[which(df$signal_percent>0.05),]

pdf(paste0(dir_fig,"df_4pc_summary_tissues.pdf"))
mg <- ggplot(df, aes(x = signal_percent, y = mean_est)) + geom_point(aes(colour = type),size=1.5,alpha=0.5)+ 
  geom_hline(data = df,linetype="dashed",  aes(yintercept = mean_background))+
  ylab("average expression level (log2RPKM)")+xlab("percentage of signal")+#xlim(0,1.2)+

#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank())
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_rect(
          color="white", fill="white", size=1.5, linetype="solid"
        ))  +
  scale_x_continuous(breaks=c(0,0.6,1.2))

# theme(
  
mg + facet_wrap( ~ tissue)+ theme(panel.spacing = unit(1.5, "lines"))
dev.off()










