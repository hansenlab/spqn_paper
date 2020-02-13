source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")
source("/users/ywang/10_25_2019/functions/functions_evaluation_Aug_5.R")
library(matrixStats)
#library(WGCNA)
library(ggplot2)
library(ggridges)
library(stringr)


dir_data="/users/ywang/10_25_2019/data/"
dir_input="/users/ywang/10_25_2019/ridge_signal_0pc_ori/data/"
dir_fig="/users/ywang/10_25_2019/ridge_signal_0pc_ori/fig/"
dir.create(dir_fig)

list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
               "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")
tissue_list=c( "adiposesub", "adrenalgl","arteryti","braince","brainco","breastma",
               "colontr","esophagusmu","heartleft")

###########################

#####

#n_PC=#sva
percent_="0_1"
# for(i in 1:9){
#   n_PC=list_nsv[i]
#   load(paste0(dir_input,percent_,"_percent_","svaPCs_",i,"_df_plot_ori.RData"))#df_plot
#   if(i==1){df_plot_merge=df_plot}else{
#     df_plot_merge=rbind(df_plot_merge,df_plot)
#   }
# }
# #save(df_plot_merge,file=paste0(dir_input,"sva_df_plot_merge_ori_facet.pdf"))
# 
# 
# 
# pdf(paste0(dir_fig,"svaPCs_",i,"_ridge_back_sig_ori_facet.pdf"))
# df_plot_merge$tissue=str_replace(df_plot_merge$tissue,"_"," ")
# 
# print(
#   ggplot(df_plot_merge, aes(x = x, height = y, y = bin, fill = group)) +
#     geom_ridgeline(size=0.1) +
#     scale_y_discrete(expand = c(0.01, 0)) +
#     #scale_x_continuous(expand = c(0.01, 0),limits =range_cor_0PC_ori) +
#     scale_x_continuous(expand = c(0.01, 0)) +
#     scale_fill_cyclical(
#       #labels = c("background","signal"),
#       values = c("#0000ff", "#ff0000"),
#       name = "group", guide = "legend"
#     )+  theme_ridges(grid = FALSE)+labs(
#       x = "correlation",
#       y = "bin" )+ facet_wrap(~tissue, ncol = 3)+  theme(strip.text.x = element_text(size = 7),strip.background =element_rect(fill="white"))#+ facet_grid(rows = vars(tissue))+ facet_wrap(~tissue, ncol = 3)
#   
# )
# dev.off()
# 
# 
# remove(df_plot_merge)
for(i in 1:9){
  n_PC=4
  load(paste0(dir_input,percent_,"_percent_","4PCs_",i,"_df_plot_est.RData"))#df_plot
  if(i==1){df_plot_merge=df_plot}else{
    df_plot_merge=rbind(df_plot_merge,df_plot)
  }
}
#save(df_plot_merge,file=paste0(dir_input,"sva_df_plot_merge_est_facet.pdf"))



pdf(paste0(dir_fig,"4PCs_",i,"_ridge_back_sig_est_facet.pdf"))
df_plot_merge$tissue=str_replace(df_plot_merge$tissue,"_"," ")
print(
  ggplot(df_plot_merge, aes(x = x, height = y, y = bin, fill = group)) +
    geom_ridgeline(size=0.1) +
    scale_y_discrete(expand = c(0.01, 0)) +
    #scale_x_continuous(expand = c(0.01, 0),limits =range_cor_0PC_ori) +
    scale_x_continuous(expand = c(0.01, 0)) +
    scale_fill_cyclical(
      #labels = c("background","signal"),
      values = c("#0000ff", "#ff0000"),
      name = "group", guide = "legend"
    )+  theme_ridges(grid = FALSE)+labs(
      x = "correlation",
      y = "bin" )+ facet_wrap(~tissue, ncol = 3)+  theme(strip.text.x = element_text(size = 7),strip.background =element_rect(fill="white"))#+ facet_grid(rows = vars(tissue))+ facet_wrap(~tissue, ncol = 3)
  
)
dev.off()

