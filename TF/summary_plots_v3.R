

#library("sva")
#library("recount", quietly = T)
#library("WGCNA", quietly = T)
# library(matrixStats)
# library(matrixStats)
# library(WGCNA)
# library(ggplot2)
# library(ggridges)
# library(SummarizedExperiment)
# library(readr)
# 
# source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")
# source("/users/ywang/10_25_2019/functions/functions_2d_quantile_diagnal_only.R")
# source("/users/ywang/10_25_2019/functions/functions_evaluation_Aug_5.R")
# 
# source("/users/ywang/10_25_2019/WGCNA/code/functions/functions_construct_modules.R")
# 
# 
# dir_data="/users/ywang/10_25_2019/data/" # counts
# dir_input="/users/ywang/10_25_2019/corplot_ori_svapcs/data/" # cor_est 
# dir_output="/users/ywang/10_25_2019/TF/data/"
# dir_fig="/users/ywang/10_25_2019/TF/fig/"
# 
# names_TF=read.csv(paste0(dir_output,"Barrera_2016_Science_TableS2.csv"),header=T)
# # > dim(names_TF)
# # [1] 1254    2
# # > ls(names_TF)
# # [1] "Enseml.Gene.ID" "Gene.Symbol"  
# 
# list_threshold_ratio=c(1,c(1:10)*3)/1000
# 
# 
# n_ori=n_est=list()
# n_ori_overlap=n_est_overlap=list()
# 
# load(paste0(dir_output,"n_est_overlap.RData"))
# load(paste0(dir_output,"n_ori_overlap.RData"))
# load(paste0(dir_output,"n_ori.RData"))
# load(paste0(dir_output,"n_est.RData"))
# 
# 
# for(i in 1){
# 
#   n_ori[[i]]=n_est[[i]]=n_ori_overlap[[i]]=n_est_overlap[[i]]=array(dim=length(list_threshold_ratio))
#   
#   n_PC=4
#   list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
#                  "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")
#   
#   
#   name_tissue=list_tissues[[i]]
#   
#   
#   # load and input raw counts, normalize and filter low expression genes(average(log2cpm)>0)
#   name_tissue=list_tissues[[i]]
#   
#   #load and input raw counts, normalize and filter low expression genes(average(log2cpm)>0)
#   name_var= paste0("gtex.",name_tissue)
#   load(paste0(dir_data, name_var,".rda" ))
#   data_raw=eval(as.name(name_var))
#   counts_raw=assay(data_raw)
#   counts_raw=counts_raw[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA")),]
#   list_filt=rpkm_filt(counts_raw,rowData(data_raw)$gene_length[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))])
#   log2rpkm_kp=list_filt$counts_log2rpkm_keep
#   
#   
#   gene_id_filt1=rownames(data_raw)[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))]
#   gene_id_filt12=gene_id_filt1[list_filt$keep][order(rowMeans(log2rpkm_kp))]
#   
#   gene_id_filt12_2=strsplit(as.vector(gene_id_filt12),".",fixed=TRUE)
#   gene_id_filt12_3=matrix(unlist(gene_id_filt12_2), ncol = 2, byrow = TRUE)
#   gene_id_filt12_format=gene_id_filt12_3[,1]
#   
#   ###### get exp level of TF
#   TF_select=which(gene_id_filt12_format %in% names_TF$"Enseml.Gene.ID")
#   # length(TF_select)
#   # [1] 755
#   log2rpkm_kp_ordered=log2rpkm_kp[order(rowMeans(log2rpkm_kp)),]
#   log2rpkm_TF=log2rpkm_kp_ordered[TF_select,]
#   
#   ave_log2RPKM=c(rowMeans(log2rpkm_TF),rowMeans(log2rpkm_kp_ordered))
#   tissue=rep(name_tissue,nrow(log2rpkm_TF)+nrow(log2rpkm_kp_ordered))
#   type=c(rep("TF",nrow(log2rpkm_TF)),rep('background',nrow(log2rpkm_kp_ordered)))
#   
#   # if(i==1){
#   #   df_exp=data.frame(ave_log2RPKM=ave_log2RPKM,tissue=tissue,type=type)
#   # }else{df_exp=rbind(df_exp,data.frame(ave_log2RPKM=ave_log2RPKM,tissue=tissue,type=type))}
#   # 
#   # ### get the number of TP of TF in ori and est - abs(cor)>threhold
#   load(paste0(dir_output,"edges_ori_abs_others_",i,"th_tissue.RData"))#edges_ori_abs
#   load(paste0(dir_output,"edges_est_abs_others_",i,"th_tissue.RData"))#edges_est_abs
#   load(paste0(dir_output,"mat_PPI_TF_others_HURI_",i,"th_tissue.RData"))#mat_PPI_TF_others_HURI
# 
#   nTF=length(TF_select)
#   
#   
#   for(j in 1:length(edges_ori_abs)){
#     edges_ori_abs[[j]]=edges_ori_abs[[j]][,c(TF_select,which(! 1:ncol(edges_ori_abs[[j]]) %in% TF_select))]
#     edges_est_abs[[j]]=edges_est_abs[[j]][,c(TF_select,which(! 1:ncol(edges_est_abs[[j]]) %in% TF_select))]
#     n_ori[[i]][j]=sum(edges_ori_abs[[j]][upper.tri(edges_ori_abs[[j]])])
#     n_est[[i]][j]=sum(edges_est_abs[[j]][upper.tri(edges_est_abs[[j]])])
#     print(j)
#   }
# 
# 
#   ###### compare with PPI - changes of true edges after spqn
#   load(paste0(dir_output,"mat_PPI_TF_others_HURI_",i,"th_tissue.RData"))#mat_PPI_TF_others_HURI
#   load(paste0(dir_output,"edges_ori_abs_others_",i,"th_tissue.RData"))#edges_ori
#   load(paste0(dir_output,"edges_est_abs_others_",i,"th_tissue.RData"))#edges_est_abs
#   
#   mat_PPI_TF_others_HURI=mat_PPI_TF_others_HURI[,c(TF_select,which(! 1:ncol(mat_PPI_TF_others_HURI) %in% TF_select))]
#   
#   for(j in 1:length(edges_ori_abs)){
#     edges_ori_abs[[j]]=edges_ori_abs[[j]][,c(TF_select,which(! 1:ncol(edges_ori_abs[[j]]) %in% TF_select))]
#     edges_est_abs[[j]]=edges_est_abs[[j]][,c(TF_select,which(! 1:ncol(edges_est_abs[[j]]) %in% TF_select))]
#     
#     a=mat_PPI_TF_others_HURI+edges_ori_abs[[j]]
#     n_ori_overlap[[i]][j]=sum(a[upper.tri(a)]==2)
#     
#     b=mat_PPI_TF_others_HURI+edges_est_abs[[j]]
#     n_est_overlap[[i]][j]=sum(b[upper.tri(b)]==2)
#     print(j)
#   }
#   
#  
#   # a=a[,1:length(TF_select)]
#   # sum(a[upper.tri(a)]==2)
#   
# }
# 
# # save(df_exp,file=paste0(dir_output,"df_exp.RData"))#edges_ori
# 
# # pdf(paste0(dir_fig,"TF_exp_box.pdf"))
# # p<-ggplot(df_exp, aes(x=tissue, y=ave_log2RPKM, color=type)) +
# #   geom_boxplot()
# # print(p)
# # dev.off()
# 
# save(n_est_overlap,file=paste0(dir_output,"n_est_overlap.RData"))#n_est_overlap
# save(n_ori_overlap,file=paste0(dir_output,"n_ori_overlap.RData"))#n_ori_overlap
# 
# save(n_ori,file=paste0(dir_output,"n_ori.RData"))#n_ori
# save(n_est,file=paste0(dir_output,"n_est.RData"))#n_est
# 



############ plots
dir_data="/users/ywang/10_25_2019/data/" # counts
dir_input="/users/ywang/10_25_2019/corplot_ori_svapcs/data/" # cor_est 
dir_output="/users/ywang/10_25_2019/TF/data/"
dir_fig="/users/ywang/10_25_2019/TF/fig/"

load(paste0(dir_output,"n_est_overlap.RData"))
load(paste0(dir_output,"n_ori_overlap.RData"))
load(paste0(dir_output,"n_ori.RData"))
load(paste0(dir_output,"n_est.RData"))

list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
               "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")

list_threshold_ratio=c(1,c(1:10)*3)/1000

  
for(i in 1:9){
  # n_est_overlap_tmp=n_est_overlap[[i]]
  # n_ori_overlap_tmp=n_ori_overlap[[i]]
  n_ori_tmp=n_ori[[i]]
  n_est_tmp=n_est[[i]]
  
  #increase_percent=(n_est_tmp-n_ori_tmp)/((n_est_tmp+n_ori_tmp)/2) *100
  increase_percent=((n_est_tmp-n_ori_tmp)/n_ori_tmp) *100
  
  increase_edge=n_est_tmp-n_ori_tmp
  nedge=c(n_ori_tmp,n_est_tmp)
  df_tmp=data.frame(increase_edge=increase_edge,increase_percent=increase_percent,
                    tissue=rep(list_tissues[i],length(list_threshold_ratio)),signal_percent=list_threshold_ratio*100)
  df_tmp2=data.frame(nedge=nedge,type=c(rep("ori",length(list_threshold_ratio)),rep("adj",length(list_threshold_ratio))),
                     tissue=rep(list_tissues[i],length(list_threshold_ratio)*2),
                     signal_percent=c(list_threshold_ratio,list_threshold_ratio)*100)
  if(i==1){df_all=df_tmp;df_all2=df_tmp2}else{df_all=rbind(df_all,df_tmp);df_all2=rbind(df_all2,df_tmp2)}
  
  
  # pdf(paste0(dir_fig,"TF_others_nedges_",i,"th_tissue2.pdf"))
  # yrange=range(n_ori_tmp,n_est_tmp)
  # plot(list_threshold_ratio,n_ori_tmp,ylim=yrange,xlab="ratio of signals",ylab="#edges",type="b")
  # points(list_threshold_ratio,n_est_tmp,col="blue",type="b")
  # dev.off()
  # 
  # pdf(paste0(dir_fig,"TF_others_nedges_overlapped_",i,"th_tissue2.pdf"))
  # yrange=range(n_ori_overlap_tmp,n_est_overlap_tmp)
  # plot(list_threshold_ratio,n_ori_overlap_tmp,col="black",ylim=yrange,xlab="ratio of signals",ylab="#edges overlapped with PPI",type="b")
  # points(list_threshold_ratio,n_est_overlap_tmp,col="blue",type="b")
  # dev.off()
  # 
  # pdf(paste0(dir_fig,"TF_others_nedges_overall_",i,"th_tissue.pdf"))
  # yrange=range(log(n_ori_overlap_tmp+1),log(n_est_overlap_tmp+1),log(n_ori_tmp+1),log(n_est_tmp+1))
  # plot(list_threshold_ratio,log(n_ori_tmp+1),col="black",ylim=yrange,type="l",xlab="ratio of signals",ylab="log (#edges + 1)")
  # lines(list_threshold_ratio,log(n_est_overlap_tmp+1),col="blue",lty=2)
  # lines(list_threshold_ratio,log(n_ori_overlap_tmp+1),col="black" ,lty=2)
  # lines(list_threshold_ratio,log(n_est_tmp+1),col="blue",lty=1)
  # legend("bottomright",legend=c("before","after","before, overlapped with PPI","after,overlapped with PPI"),,col=c("black","blue","black","blue"),lty=c(1,1,2,2))
  # dev.off()
}




library(ggplot2)
df_all$signal_percent=as.character(df_all$signal_percent)
pdf(paste0(dir_fig,"increase_percent_alltissues_v3.pdf"))
# p <- ggplot(data=df_all, aes(x=signal_percent, y=increase_percent,group = signal_percent)) +  geom_boxplot()+
#   geom_point(position="jitter", alpha=0.6)
# print(p)
ggplot(df_all, aes(x=signal_percent, y=increase_percent)) +
  geom_boxplot(position=position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8))+ labs(x = "percentage of signals",y="percentage increase of estimated edges")
dev.off()


df_all$signal_percent=as.character(df_all$signal_percent)
pdf(paste0(dir_fig,"increase_percent_alltissues.pdf"))
# p <- ggplot(data=df_all, aes(x=signal_percent, y=increase_percent,group = signal_percent)) +  geom_boxplot()+
#   geom_point(position="jitter", alpha=0.6)
# print(p)
ggplot(df_all, aes(x=signal_percent, y=increase_percent)) +
  geom_boxplot(position=position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8))
dev.off()








#################### need renew
pdf(paste0(dir_fig,"increase_percent_alltissues2.pdf"))
colfunc <- colorRampPalette(c("green", "red"))
cols=colfunc(9)
n_ori_tmp=n_ori[[1]]
n_est_tmp=n_est[[1]]
increase_percent=(n_est_tmp-n_ori_tmp)/((n_est_tmp+n_ori_tmp)/2) *100
plot(list_threshold_ratio,increase_percent,col=cols[1])
for(i in 2:9){
  n_ori_tmp=n_ori[[i]]
  n_est_tmp=n_est[[i]]
  increase_percent=(n_est_tmp-n_ori_tmp)/((n_est_tmp+n_ori_tmp)/2) *100
  lines(list_threshold_ratio,increase_percent,col=cols[i])
}
#legend()
dev.off()


library(ggplot2)
pdf(paste0(dir_fig,"increase_nedge_alltissues.pdf"))
# p <- ggplot(data=df_all, aes(x=signal_percent, y=log(increase_edge+1),group = signal_percent)) +  geom_boxplot()+
#   geom_point(position="jitter", alpha=0.6)
# print(p)
df_all$signal_percent=as.character(df_all$signal_percent)
ggplot(df_all, aes(x=signal_percent, y=increase_edge)) +
  geom_boxplot(position=position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8))
dev.off()


pdf(paste0(dir_fig,"increase_nedge_alltissues2.pdf"))
colfunc <- colorRampPalette(c("green", "red"))
cols=colfunc(9)
n_ori_tmp=n_ori[[1]]
n_est_tmp=n_est[[1]]
increase_edge=(n_est_tmp-n_ori_tmp)
plot(list_threshold_ratio,increase_edge,col=cols[1],type="l")
for(i in 1:9){
  n_ori_tmp=n_ori[[i]]
  n_est_tmp=n_est[[i]]
  increase_edge=(n_est_tmp-n_ori_tmp)
  lines(list_threshold_ratio,increase_edge,col=cols[i])
}
dev.off()


library(ggplot2)
pdf(paste0(dir_fig,"nedge_compare_alltissues.pdf"))
df_all2$signal_percent=as.character(df_all2$signal_percent)
df_all2$type=factor(df_all2$type,levels=c("ori","adj"))
ggplot(df_all2, aes(x=signal_percent, y=log(nedge+1), color=type)) +
  geom_boxplot(position=position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8))
dev.off()


library(ggplot2)
pdf(paste0(dir_fig,"nedge_compare_alltissues2.pdf"))
df_all2$signal_percent=as.character(df_all2$signal_percent)
df_all2$type=factor(df_all2$type,levels=c("ori","adj"))
ggplot(df_all2, aes(x=signal_percent, y=log(nedge+1))) +
  geom_boxplot(position=position_dodge(0.8), aes(col=type))+
  geom_jitter(position=position_dodge(0.8),aes(col=tissue))
dev.off()

# 
# nedges_estimated=c(n_1_ori,n_1_est)
# type=c(rep("before SpQN",9),rep("after SpQN",9))
# tissue=c(list_tissues,list_tissues)
# df=data.frame(nedges_estimated=nedges_estimated,type=type,tissue=tissue)
# df$type=factor(df$type,levels=c("before SpQN","after SpQN"))
# p <- ggplot(data=df, aes(x=tissue, y=nedges_estimated, fill=type)) +
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   theme_minimal()+theme(axis.text.x = element_text(angle = 45))
# pdf(paste0(dir_fig,"TF_others_nedges_estimated_1percent.pdf"))
# print(p)
# dev.off()
# 
# 
# nedges_estimated=c(n_01_ori,n_01_est)
# type=c(rep("before SpQN",9),rep("after SpQN",9))
# tissue=c(list_tissues,list_tissues)
# df=data.frame(nedges_estimated=nedges_estimated,type=type,tissue=tissue)
# df$type=factor(df$type,levels=c("before SpQN","after SpQN"))
# p <- ggplot(data=df, aes(x=tissue, y=nedges_estimated, fill=type)) +
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   theme_minimal()+theme(axis.text.x = element_text(angle = 45))
# pdf(paste0(dir_fig,"TF_others_nedges_estimated_01percent.pdf"))
# print(p)
# dev.off()
# 
# 
# nedges_estimated=c(n_01_ori_overlap,n_01_est_overlap)
# type=c(rep("before SpQN",9),rep("after SpQN",9))
# tissue=c(list_tissues,list_tissues)
# df=data.frame(nedges_estimated=nedges_estimated,type=type,tissue=tissue)
# df$type=factor(df$type,levels=c("before SpQN","after SpQN"))
# p <- ggplot(data=df, aes(x=tissue, y=nedges_estimated, fill=type)) +
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   theme_minimal()+theme(axis.text.x = element_text(angle = 45))
# pdf(paste0(dir_fig,"TF_others_nedges_estimated_01percent_overlap.pdf"))
# print(p)
# dev.off()
# 
# nedges_estimated=c(n_1_ori_overlap,n_1_est_overlap)
# type=c(rep("before SpQN",9),rep("after SpQN",9))
# tissue=c(list_tissues,list_tissues)
# df=data.frame(nedges_estimated=nedges_estimated,type=type,tissue=tissue)
# df$type=factor(df$type,levels=c("before SpQN","after SpQN"))
# p <- ggplot(data=df, aes(x=tissue, y=nedges_estimated, fill=type)) +
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   theme_minimal()+theme(axis.text.x = element_text(angle = 45))
# pdf(paste0(dir_fig,"TF_others_nedges_estimated_1percent_overlap.pdf"))
# print(p)
# dev.off()
