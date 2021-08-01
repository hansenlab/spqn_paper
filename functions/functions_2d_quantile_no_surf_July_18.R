rpkm<-function(counts,gene_length){
  counts_rpkm=t( (t(counts+0.5)/colSums(counts)) * (1e+9))/gene_length
  counts_rpkm
}

rpkm_filt <- function(counts_raw,gene_length){
  counts_log2rpkm=log2(rpkm(counts_raw,gene_length))
  keep=which(rowMedians(data.matrix(counts_log2rpkm))>0)
  counts_log2rpkm_keep=counts_log2rpkm[keep,]
  counts_raw_keep=counts_raw[keep,]
  return(list(counts_raw_keep=counts_raw_keep,counts_log2rpkm_keep=counts_log2rpkm_keep,keep=keep))
}

removePCs <- function(counts_logrpkm, num_PCs = 4){
  if(num_PCs > 0) {
    counts_removePCs= t(removePrincipalComponents( scale(t(counts_logrpkm)), num_PCs))
    return(counts_removePCs)} else{return(counts_logrpkm)}
}



##### step 1, asssign running bins
#group the sorted genes into equal-size running groups, each group containing 400 genes 
#avelog2cpm=ave_logcpm_kp
#ngrp=20
#size_grp=400
get_grps<-function(avelog2cpm,ngrp=20,size_grp=400){
  ngene=length(avelog2cpm) 
  
  if(size_grp-length(grp_loc0[[1]])<5){
    grp_label=cut(1:ngene,ngrp)
    grp_loc0=split(1:ngene,grp_label)
    for(i in 1:(ngrp)){
      grp_loc[[i]]=grp_loc0[[i]]
    }
  }else{
    grp_label=cut(1:(ngene-size_grp+1),ngrp-1)
    grp_loc0=split(1:(ngene-size_grp+1),grp_label)
    grp_loc=list()
    for(i in 1:(ngrp-1)){
      grp_loc[[i]]=c(grp_loc0[[i]][1]:(grp_loc0[[i]][1]+size_grp-1))
    }
    grp_loc[[ngrp]]=c((ngene-size_grp+1):ngene)
  }
  grp_loc
}


##### step 2, asssign inner bins
# in each running group, get a inner group
get_grps_inner<-function(grp_loc){
  grp_loc_inner=list()
  ngrp=length(grp_loc)
  size_bin=length(grp_loc[[1]])
  ngene=max(grp_loc[[ngrp]])
  
  width_tmp=grp_loc[[2]][1] - grp_loc[[1]][1] 
  grp_loc_inner[[1]]= c(1:round(size_bin/2+width_tmp/2))
  
  for(i in 2:(ngrp-1)){
    width_tmp=grp_loc[[i+1]][1] - grp_loc[[i]][1] 
    grp_loc_inner[[i]]=c( (tail(grp_loc_inner[[i-1]],1)+1) :(tail(grp_loc_inner[[i-1]],1) + width_tmp))
  }  
  
  grp_loc_inner[[ngrp]]=c( (tail(grp_loc_inner[[ngrp-1]],1)+1) : ngene)
  
  grp_loc_inner
}


##### step 3, get rank for each running bin
#rank_bin=get_bin_rank(corr_ori,group_loc,group_loc_adj,cor_ref) #6s for 20*20 bins, 4000 genes ;2min for all genes, 20*20 bins
#cor_obs=corr_ori
#grp_loc=group_loc
#grp_loc_inner=group_loc_adj
#cor_ref=cor_ref
get_bin_rank<-function(cor_obs,grp_loc,grp_loc_inner,cor_ref){
  rank_bin=rank_bin_pre=array(dim=dim(cor_obs))
  ngrp=length(grp_loc)
  size_bin=length(grp_loc[[1]])

  l_cor_tmp_ref=length(cor_ref[upper.tri(cor_ref)])
  
  for(i in 1:ngrp){
    for(j in i:ngrp){
      cor_bin_tmp=cor_obs[grp_loc[[i]],grp_loc[[j]]]
      rank_bin_tmp=array(dim=dim(cor_bin_tmp))
      l_cor_tmp=length(cor_bin_tmp)
      
      rank_bin_tmp[1:l_cor_tmp]=Rank(cor_bin_tmp[1:l_cor_tmp])
      
      #number of diagonals(of full correlation matrix) contained in the bin
      n_diag=sum(grp_loc[[i]] %in% grp_loc[[j]])
      
      #scale the rank of each bin to same scale as rank(cor_ref), 
      #after scaling, rank_bin could contain non-integars 
      rank_bin_pre[grp_loc[[i]],grp_loc[[j]]]=rank_bin_tmp
      rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]]=1 + ((rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]]-1)/(l_cor_tmp-n_diag) *(l_cor_tmp_ref-1))
        
      rank_bin[grp_loc_inner[[i]],grp_loc_inner[[j]]]=rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]]
      
      print(c(i,j))
    }
  }
  #rank_bin[lower.tri(rank_bin)]=t(rank_bin)[lower.tri(rank_bin)]
  #diag(rank_bin)=1
  rank_bin
}

get_bin_rank2<-function(cor_obs,grp_loc,grp_loc_inner,cor_ref){
  rank_bin=rank_bin_pre=array(dim=dim(cor_obs))
  ngrp=length(grp_loc)
  size_bin=length(grp_loc[[1]])
  
  l_cor_tmp_ref=length(cor_ref[upper.tri(cor_ref)])
  
  for(i in 1:ngrp){
    for(j in i:ngrp){
      cor_bin_tmp=cor_obs[grp_loc[[i]],grp_loc[[j]]]
      rank_bin_tmp=array(dim=dim(cor_bin_tmp))
      l_cor_tmp=length(cor_bin_tmp)
      
      rank_bin_tmp[1:l_cor_tmp]=rank(cor_bin_tmp[1:l_cor_tmp])
      
      #number of diagonals(of full correlation matrix) contained in the bin
      n_diag=sum(grp_loc[[i]] %in% grp_loc[[j]])
      
      #scale the rank of each bin to same scale as rank(cor_ref),
      #after scaling, rank_bin could contain non-integars
      rank_bin_pre[grp_loc[[i]],grp_loc[[j]]]=rank_bin_tmp
      #rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]]=1 + ((rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]]-1)/(l_cor_tmp-n_diag) *(l_cor_tmp_ref-1))
      #rank_bin[grp_loc_inner[[i]],grp_loc_inner[[j]]]=rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]]
      
      #rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]]=1 + ((rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]]-1)/(l_cor_tmp-n_diag) *(l_cor_tmp_ref-1))
      rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]]= 1+((rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]]-1)/(l_cor_tmp-n_diag-1) *(l_cor_tmp_ref-1))
      rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]][which(rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]]>l_cor_tmp_ref)]=l_cor_tmp_ref
      
      rank_bin[grp_loc_inner[[i]],grp_loc_inner[[j]]]=rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]]
      
      
      
      print(c(i,j))
    }
  }
  #rank_bin[lower.tri(rank_bin)]=t(rank_bin)[lower.tri(rank_bin)]
  #diag(rank_bin)=1
  rank_bin
}










##### step 4, transform rank to cor_est 

#cor_est=est_cor(rank_bin,cor_ref)#6s for 20*20 bins

est_cor<-function(rank_bin,cor_ref){
  
  cor_adj=array(dim=dim(rank_bin))
  cor_ref_sorted=Sort(cor_ref[upper.tri(cor_ref)])
  
  up_tri=upper.tri(cor_adj)#15s for all genes
  
  #find to nearest integars to rank_bin, since rank_bin could contain non-integars
  rank_bin=rank_bin[up_tri]
  rank_bin1=rank_bin %/% 1
  
  #assign weights to each rank according to the distance to rank_bin
  rank_bin_w2 = (rank_bin - rank_bin1)#15s for all genes
  
  rank_bin2=rank_bin1+1
  rank_bin_w1 = 1-rank_bin_w2
  
  #find the correlations in the cor_ref corresponding to the two nearest ranks to rank_bin
  #estimate the correlation using weighted average based on the distance to rank_bin
  cor_adj[up_tri] =  rank_bin_w1*cor_ref_sorted[rank_bin1]+ rank_bin_w2*cor_ref_sorted[rank_bin2]#1.3min for all genes
  
  remove(rank_bin1)
  remove(rank_bin2)
  remove(rank_bin_w1)
  remove(rank_bin_w2)
  remove(up_tri)
  remove(rank_bin)
  
  low_tri=lower.tri(cor_adj)
  cor_adj[low_tri]=t(cor_adj)[low_tri]
  diag(cor_adj)=1
  cor_adj
}

est_cor2<-function(rank_bin,cor_ref){
  
  cor_adj=array(dim=dim(rank_bin))
  cor_ref_sorted=sort(cor_ref[upper.tri(cor_ref)])
  
  up_tri=upper.tri(cor_adj)#15s for all genes
  
  #find to nearest integars to rank_bin, since rank_bin could contain non-integars
  rank_bin=rank_bin[up_tri]
  rank_bin1=rank_bin %/% 1
  
  #assign weights to each rank according to the distance to rank_bin
  rank_bin_w2 = (rank_bin - rank_bin1)#15s for all genes
  
  rank_bin2=rank_bin1+1
  rank_bin_w1 = 1-rank_bin_w2
  
  #find the correlations in the cor_ref corresponding to the two nearest ranks to rank_bin
  #estimate the correlation using weighted average based on the distance to rank_bin
  rank_bin2[which(rank_bin2>length(cor_ref_sorted))]=length(cor_ref_sorted)
  cor_adj[up_tri] =  rank_bin_w1*cor_ref_sorted[rank_bin1]+ rank_bin_w2*cor_ref_sorted[rank_bin2]#1.3min for all genes
  
  remove(rank_bin1)
  remove(rank_bin2)
  remove(rank_bin_w1)
  remove(rank_bin_w2)
  remove(up_tri)
  remove(rank_bin)
  
  low_tri=lower.tri(cor_adj)
  cor_adj[low_tri]=t(cor_adj)[low_tri]
  diag(cor_adj)=1
  cor_adj
}
