get_grps_disjoint<-function(avelog2cpm,ngrp=20){
  grp_loc_label=cut(1:length(avelog2cpm),ngrp)
  grp_loc=split(1:length(avelog2cpm),grp_loc_label)
  grp_loc
}

get_grps_disjoint_size<-function(avelog2cpm,size_grp=400){
  ngrp=round(length(avelog2cpm)/size_grp)
  grp_loc_label=cut(1:length(avelog2cpm),ngrp)
  grp_loc=split(1:length(avelog2cpm),grp_loc_label)
  print(ngrp)
  grp_loc
}
