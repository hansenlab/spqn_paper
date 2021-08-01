


##### step 1, asssign running bins
#group the sorted genes into equal-size running groups, each group containing 400 genes 
#avelog2cpm=ave_logcpm_kp
#ngrp=20
#size_grp=400
get_grps<-function(avelog2cpm,ngrp=20,size_grp=400){
  ngene=length(avelog2cpm) 
  
  grp_label=cut(1:(ngene-size_grp+1),ngrp-1)
  grp_loc0=split(1:(ngene-size_grp+1),grp_label)
  
  grp_loc=list()
  for(i in 1:(ngrp-1)){
    grp_loc[[i]]=c(grp_loc0[[i]][1]:(grp_loc0[[i]][1]+size_grp-1))
  }
  grp_loc[[ngrp]]=c((ngene-size_grp+1):ngene)
  
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


##### step 3, get rank for each running bin-diag only
#extract locations of diagonal bins in the matrix
get_diag_loc<-function(bin_list){
  n_grp=length(bin_list)
  list_bin=combn(bin_list[[1]],2)
  for(k in 2:n_grp){
    list_bin=cbind(list_bin,combn(bin_list[[k]],2))
  }
  paste0(list_bin[1,],"_",list_bin[2,])
}

#=
#group_loc_eva=split(1:n_gene,cut(1:n_gene,n_eva))

# get_diag_loc2<-function(group_loc_eva,avelog2cpm){
#   n_gene=length(avelog2cpm)
#   binlist=split(1:n_gene,group_loc_eva)
#   get_diag_loc(binlist)
# }


#create a list of adjustment inner bins that overlapped with evaluation bins
overlapped_pairs<-function(group_loc_eva,avelog2cpm,grp_inner){
  list_overlapped_pairs1=list_overlapped_pairs2=c()
  #create evaluation bins list
  list_diag_eva=get_diag_loc(group_loc_eva)
  for(j in 1:length(grp_inner)){
    for(k in j:length(grp_inner)){
      #create adjustment inner bins list
      #list_bin_tmp=expand.grid(grp_inner[[j]], grp_inner[[k]])
      # find overlap between the two lists
      loc_overlap=paste0(list_bin_tmp[,1],"_",list_bin_tmp[,2]) %in% list_diag_eva
      if(sum(loc_overlap)>0){
        list_overlapped_pairs1=c(list_overlapped_pairs1,j)
        list_overlapped_pairs2=c(list_overlapped_pairs2,k)
      }
    }
  }
  cbind(list_overlapped_pairs1,list_overlapped_pairs2)
}

overlapped_pairs2<-function(group_loc_eva,grp_inner,size_neval,size_inner){
  
  list_overlapped_pairs1=list_overlapped_pairs2=c()
  for(j in 1:length(grp_inner)){
    for(k in j:length(grp_inner)){
      
      list_choose=c()
      for(l in 1:length(group_loc_eva)){
        ifin=grp_inner[[j]][length(grp_inner[[j]])] %in% group_loc_eva[[l]] &&  grp_inner[[k]][1] %in% group_loc_eva[[l]]
        list_choose=c(list_choose,ifin)
        ifin=grp_inner[[j]][length(grp_inner[[j]])] %in% group_loc_eva[[l]] &&  grp_inner[[k]][length(grp_inner[[k]])] %in% group_loc_eva[[l]]
        list_choose=c(list_choose,ifin)
        ifin=grp_inner[[j]][1] %in% group_loc_eva[[l]] &&  grp_inner[[k]][1] %in% group_loc_eva[[l]]
        list_choose=c(list_choose,ifin)
        ifin=grp_inner[[j]][1] %in% group_loc_eva[[l]] &&  grp_inner[[k]][length(grp_inner[[k]])] %in% group_loc_eva[[l]]
        list_choose=c(list_choose,ifin)
      }
      
      for(l in 1:length(group_loc_eva)){
        ifin=group_loc_eva[[l]][length(group_loc_eva[[l]])] %in% grp_inner[[j]] &&  group_loc_eva[[l]][1] %in% grp_inner[[k]]
        list_choose=c(list_choose,ifin)
        ifin=group_loc_eva[[l]][length(group_loc_eva[[l]])] %in% grp_inner[[j]] &&  group_loc_eva[[l]][length(group_loc_eva[[l]])] %in% grp_inner[[k]]
        list_choose=c(list_choose,ifin)
        ifin=group_loc_eva[[l]][1] %in% grp_inner[[j]] &&  group_loc_eva[[l]][1] %in% grp_inner[[k]]
        list_choose=c(list_choose,ifin)
        ifin=group_loc_eva[[l]][1] %in% grp_inner[[j]] &&  group_loc_eva[[l]][length(group_loc_eva[[l]])] %in% grp_inner[[k]]
        list_choose=c(list_choose,ifin)
      }
      
      
      if(sum(list_choose)>0){
        list_overlapped_pairs1=c(list_overlapped_pairs1,j)
        list_overlapped_pairs2=c(list_overlapped_pairs2,k)
      }
      
    }
  }
  
  cbind(list_overlapped_pairs1,list_overlapped_pairs2)
}

#not using rfast
#rank_bin=get_bin_rank2_diag_only(corr_ori,group_loc,group_loc_adj,cor_ref,list_overlapped_pairs) #6s for 20*20 bins, 4000 genes; 2min for all genes,ngrp=20,size_grp=1000

#cor_obs=corr_ori;grp_loc=group_loc;grp_loc_inner=group_loc_adj;
get_bin_rank2_diag_only<-function(cor_obs,grp_loc,grp_loc_inner,cor_ref,list_overlapped_pairs){
  rank_bin=rank_bin_pre=array(dim=dim(cor_obs))
  ngrp=length(grp_loc)
  size_bin=length(grp_loc[[1]])
  
  l_cor_tmp_ref=length(cor_ref[upper.tri(cor_ref)])
  
  for(i in 1:ngrp){
    list_j=list_overlapped_pairs[,2][which(list_overlapped_pairs[,1]==i)]
    list_j_upper=intersect(c(i:ngrp),list_j)
    for(j in list_j_upper){
      cor_bin_tmp=cor_obs[grp_loc[[i]],grp_loc[[j]]]
      rank_bin_tmp=array(dim=dim(cor_bin_tmp))
      l_cor_tmp=length(cor_bin_tmp)
      
      #loc_diag=grp_loc[[i]][which(grp_loc[[i]] %in% grp_loc[[j]])]
      rank_bin_tmp[1:l_cor_tmp]=rank(cor_bin_tmp[1:l_cor_tmp])
      #number of diagonals(of full correlation matrix) contained in the bin
      n_diag=sum(grp_loc[[i]] %in% grp_loc[[j]])
      rank_bin_tmp[which(rank_bin_tmp> l_cor_tmp-n_diag)]=l_cor_tmp-n_diag
      
      #scale the rank of each bin to same scale as rank(cor_ref),
      #after scaling, rank_bin could contain non-integars
      rank_bin_pre[grp_loc[[i]],grp_loc[[j]]]=rank_bin_tmp
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
#cor_est=est_cor2(rank_bin,cor_ref)#6s for 20*20 bins, 4000 genes; 3.5min for all genes, 20*20 bins,siez_grp=1000
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


