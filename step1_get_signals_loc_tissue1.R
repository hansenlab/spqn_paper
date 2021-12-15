
########### # get the adjusted correlation for shared expressed genes only, for each tissue of gtex data
########### # by applying spqn on cor_ori that containing shared genes only
library(matrixStats)
library(WGCNA)
library(spqn)
library(SummarizedExperiment)
source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")

source("/users/ywang/July_28_2020/functions/plot_signal_condition_exp.R")

# source("/users/ywang/10_25_2019/functions/functions_evaluation_Aug_5.R")

# dir_data="/users/ywang/10_25_2019/data/"
dir_data="/dcl01/hansen/data/meanCoexp/"
# dir_input="/users/ywang/July_28_2020/regulatome/output/"
dir_output="/users/ywang/July_28_2020/VST/output_4pc/"
setwd(dir_output)
# dir.create(dir_output)
# load(paste0(dir_output,"TF_data_unique_en_filt.RData")) #TF_data_unique_en_filt


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



plot_signal_condition_exp<-function (cor_mat, ave_exp, signal){
  idx <- order(ave_exp)
  cor_mat <- cor_mat[idx, idx]
  if (isTRUE(all.equal(signal, 0))) {
    ngrp <- 10
    ngene <- ncol(cor_mat)
    idx <- seq_len(ngene)
    grp_label <- cut(idx, ngrp)
    grp_loc <- split(idx, grp_label)
    length_group <- lapply(grp_loc, length)
    length_group <- unlist(length_group)
    length_group <- length_group * (length_group - 1)/2
    length_group_cumulate <- cumsum(length_group)
    cor_vec_all <- data.frame(matrix(ncol = 2, nrow = length_group_cumulate[10]))
    colnames(cor_vec_all) <- c("correlation", "group")
    for (i in seq_len(10)) {
      cor_tmp <- cor_mat[grp_loc[[i]], grp_loc[[i]]]
      cor_tmp <- cor_tmp[upper.tri(cor_tmp)]
      cor_tmp <- as.numeric(cor_tmp)
      if (i == 1) {
        idx1 <- 1
      }
      else {
        idx1 <- length_group_cumulate[i - 1] + 1
      }
      idx2 <- length_group_cumulate[i]
      idx <- c(idx1:idx2)
      cor_vec_all$correlation[idx] <-cor_tmp
      length_tmp <- length_group[i]
      cor_vec_all$group[idx] <- rep(i, length_tmp)
    }
    cor_vec_all <- data.frame(cor_vec_all)
    names(cor_vec_all) <- c("correlation", "group")
    cor_vec_all$group <- as.factor(cor_vec_all$group)
    ggplot(cor_vec_all, aes_string(x = "correlation", y = "group")) +
      geom_density_ridges2(fill = "blue") + theme_ridges(grid = TRUE) +
      theme(axis.title.x = element_blank()) + geom_vline(xintercept = 0,
                                                         linetype = "dotted", color = "black", size = 0.3)+
      geom_vline(xintercept = c(-.5,0,.5), linetype="dotted", color = "black", size=0.3)+
      theme(legend.position = "none")
    
  }
  else {
    ncor <- length(which(upper.tri(cor_mat)))
    nsignal <- ncor * signal
    nsignal_grp <- round(nsignal/100)
    grp_locs <- get_grp_loc(cor_mat)
    list_cor_sig <- grp_sig <- numeric(10 * nsignal_grp)
    nbackground_group <- lapply(grp_locs, length)
    nbackground_group <- unlist(nbackground_group)
    nbackground_group <- nbackground_group * (nbackground_group -
                                                1)/2
    list_cor_back_cumulate <- cumsum(nbackground_group)
    list_cor_back <- grp_back <- numeric(list_cor_back_cumulate[10])
    for (ngrp in seq_len(10)) {
      loc_back <- grp_locs[[ngrp]]
      cor_back <- cor_mat[loc_back, loc_back]
      cor_back <- cor_back[upper.tri(cor_back)]
      order_cor_ori_grp <- order(abs(cor_back), decreasing = TRUE)
      idx <- seq_len(nsignal_grp)
      id_signal <- order_cor_ori_grp[idx]
      cor_signal_ori <- cor_back[id_signal]
      idx1 <- (ngrp - 1) * nsignal_grp + 1
      idx2 <- ngrp * nsignal_grp
      idx <- c(idx1:idx2)
      list_cor_sig[idx] <- cor_signal_ori
      grp_sig[idx] <- ngrp
      if (ngrp == 1) {
        idx <- seq_len(list_cor_back_cumulate[ngrp])
        list_cor_back[idx] <- cor_back
        grp_back[idx] <- ngrp
      }
      else {
        idx1 <- list_cor_back_cumulate[ngrp - 1] + 1
        idx2 <- list_cor_back_cumulate[ngrp]
        idx <- c(idx1:idx2)
        list_cor_back[idx] <- cor_back
        grp_back[idx] <- ngrp
      }
    }
    df_sig <- data.frame(correlation = list_cor_sig, bin = grp_sig,
                         group = "signal")
    df_back <- data.frame(correlation =list_cor_back, bin = grp_back,
                          group = "background")
    df_merge <- rbind(df_sig, df_back)
    df_merge$bin_group <- paste(df_merge$bin, df_merge$group)
    # ggplot(df_merge, aes_string(y = "bin")) + geom_density_ridges(aes_string(x = "correlation",
    #                                                                          fill = "bin_group"), alpha = 0.8, color = "white") +
    #   labs(x = "correlation", y = "bin") + scale_y_discrete(limits = seq_len(10))+
    #   scale_fill_cyclical(labels = c("background", "signal"),
    #                       values = c("#0000ff", "#ff0000"), name = "group",
    #                       guide = "legend") + theme_ridges(grid = FALSE)+
    #   geom_vline(xintercept = c(-.5,0,.5), linetype="dotted", color = "black", size=0.3)
    maxy_back=get_maxy(df_merge)
    maxy_sig=get_maxy(df_merge, group = "signal")
    df_plot <- rbind(summarizeDensity2(df_merge, maxy=maxy_back),
                     summarizeDensity2(df_merge, maxy=maxy_sig,group = "signal"))
    
    
     p_= ggplot(df_plot, aes(x = x,height=y, y = bin, fill = group)) +
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
        y = "bin" )+   theme(strip.text.x = element_text(size = 7),strip.background =element_rect(fill="white"))+#+ facet_grid(rows = vars(tissue))+ facet_wrap(~tissue, ncol = 3)
      geom_vline(xintercept = c(-.5,0,.5), linetype="dotted", color = "black", size=0.3)+
      theme(legend.position = "none")
    p_
  }
}

# i=1
# 
for(i in c(1:9)){
  print('------------------------------------------')
  
  print(i)
  
  n_PC=4
  
  ########### for each tissue, format the data 
  list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
                 "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")
  name_tissue=list_tissues[[i]]
  
  #load and input raw counts, normalize and filter low expression genes(average(log2cpm)>0)
  name_var= paste0("gtex.",name_tissue)
  load(paste0(dir_data, name_var,".rda" ))
  data_raw=eval(as.name(name_var))
  counts_raw=assay(data_raw)
  counts_raw=counts_raw[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA")),]
  list_filt=rpkm_filt(counts_raw,rowData(data_raw)$gene_length[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))])
  
  log2rpkm_kp=list_filt$counts_log2rpkm_keep
  ls(list_filt)
  
  dim(log2rpkm_kp)
  
  ### name the tissue-wise expressed gene exp matrix, as their ensemble id
  
  log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)#remove 4PC from the expressed genes in this tissue
  
  cor_ori=cor(t(log2rpkm_kp_rmPC))
  
  #### plot  - before using vst
  library(spqn)
  library(ggplot2)
  library(matrixStats)
  library(stats)
  
  # pdf("signal_condition_exp_rm4PCs.pdf")
  # pdf(paste0("signal_condition_exp_rm4PCs_,",name_tissue,",.pdf"))
  p=(plot_signal_condition_exp(cor_ori, rowMeans(log2rpkm_kp), signal=0.001))
  # dev.off()
  
  ggsave(p, 
         filename = paste0("signal_condition_exp_rm4PCs_,",name_tissue,",.pdf"),
         height = 2, width = 2, units = "in")
  
  ### using vst in DESeq2
  library(DESeq2)
  counts_raw_keep_vst=varianceStabilizingTransformation(list_filt$counts_raw_keep)  
  save(counts_raw_keep_vst,file=paste0("counts_raw_keep_vst_",name_tissue,".RData"))
  
  sum(rownames(counts_raw_keep_vst)==rownames(list_filt$counts_raw_keep))
  dim(counts_raw_keep_vst)
  
  # normalize by lib.size
  list_filt_vst=rpkm_filt(counts_raw_keep_vst,rowData(data_raw)$gene_length[rownames(counts_raw_keep_vst)])
  log2rpkm_kp_vst=list_filt_vst$counts_log2rpkm_keep
  # rm 4PCs
  dim(log2rpkm_kp_vst)
  # [1] 12267   350
  dim(counts_raw_keep_vst)
  
  log2rpkm_kp_rmPCt_vst = removePCs(log2rpkm_kp_vst,n_PC)#remove 4PC from the expressed genes in this tissue
  
  # get cor
  cor_ori_vst=cor(t(log2rpkm_kp_rmPCt_vst))
  # var(as.vector(cor_ori_vst[rowMeans(log2rpkm_kp)>7,rowMeans(log2rpkm_kp)>7]))
  
  summary(rowMeans(log2rpkm_kp))
  # -2.303   1.750   2.883   2.982   4.008  13.700
  
  # pdf(paste0("signal_condition_exp_VST_rm4PCs_",name_tissue,".pdf"))
  p=(plot_signal_condition_exp(cor_ori_vst, rowMeans(log2rpkm_kp), signal=0.001))
  # dev.off()
  ggsave(p,
         filename = paste0("signal_condition_exp_VST_rm4PCs_,",name_tissue,".pdf"),
         height = 2, width = 2, units = "in")
  # ggsave(p, 
         # filename = paste0("signal_condition_exp_VST_rm4PCs_,",name_tissue,".pdf"))
  
  # pdf(paste0("signal_condition_exp_VST_rm4PCs_",name_tissue,"_v2.pdf"))
  p=(plot_signal_condition_exp(cor_ori_vst, rowMeans(log2rpkm_kp_vst), signal=0.001))
  # dev.off()
  ggsave(p, 
         filename = paste0("signal_condition_exp_VST_rm4PCs_,",name_tissue,"_v2.pdf"),
         height = 2, width = 2, units = "in")
  
}

  
  



