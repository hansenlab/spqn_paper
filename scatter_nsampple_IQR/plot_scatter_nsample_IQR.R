
dir_input="=/scatter_nsampple_IQR/data/"
dir_fig="/scatter_nsampple_IQR/fig/"
dir.create(dir_fig)

load(paste0(dir_input,"list_nsample.RData"))


###### plot
# 0PCs
n_PC=0
# load(paste0(dir_input,"IQR_ori_all_",n_PC,"_PCs.RData"))
# 
# pdf(paste0(dir_fig,"scatter_nsample_IQR_0PC.pdf"))
# plot(list_nsample,IQR_ori_all,col="blue",cex=0.2,cex.lab=1.4,pch=1,xlab="number of samples",ylab="IQR")
# #points(list_nsample,IQR_adj_all,col="red",cex=0.2,cex.lab=1.4,pch=1)
# #legend("topright",legend=c("observed","adjusted"),col=c("blue","red"),pch=1,cex=1.4,pt.cex=0.5)
# dev.off()
# # > for(i in 1:9){print(list_nsample[1+100*(i-1)])}
# # [1] 238
# # [1] 100
# # [1] 227
# # [1] 100
# # [1] 87
# # [1] 130
# # [1] 130
# # [1] 193
# # [1] 157
# #c(2,4)
# #c(6,7)
# 
# 
# 
# # 4PCs
# n_PC=4
# load(paste0(dir_input,"IQR_ori_all_",n_PC,"_PCs.RData"))
# load(paste0(dir_input,"IQR_adj_all_",n_PC,"_PCs.RData"))
# 
# pdf(paste0(dir_fig,"scatter_nsample_IQR_4PC.pdf"))
# plot(list_nsample,IQR_ori_all,col="blue",cex=0.2,cex.lab=1.4,pch=1,xlab="number of samples",ylab="IQR")
# points(list_nsample,IQR_adj_all,col="red",cex=0.2,cex.lab=1.4,pch=1)
# legend("topright",legend=c("observed","adjusted"),col=c("blue","red"),pch=1,cex=1.4,pt.cex=0.5)
# dev.off()


# sva_PCs
load(paste0(dir_input,"IQR_ori_all_","sva","_PCs.RData"))
load(paste0(dir_input,"IQR_adj_all_","sva","_PCs.RData"))

pdf(paste0(dir_fig,"scatter_nsample_IQR_svaPC.pdf"))
plot(list_nsample,IQR_ori_all,col="blue",cex=0.2,cex.lab=1.4,pch=1,xlab="number of samples",ylab="IQR")
points(list_nsample,IQR_adj_all,col="red",cex=0.2,cex.lab=1.4,pch=1)
legend("topright",legend=c("observed","adjusted"),col=c("blue","red"),pch=1,cex=1.4,pt.cex=0.5)
dev.off()


