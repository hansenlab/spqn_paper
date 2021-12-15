dir_fig="/corplot_adj_4_svapcs/fig/"
dir.create(dir_fig)

dir_input="/corplot_adj_4_svapcs/data/"
load(paste0(dir_input,"min_IQR.RData"))

range_y=range(min_IQR$IQR_vec_ori,min_IQR$IQR_vec_adj)

pdf(paste0(dir_fig,"min_IQR_ori.pdf"))
plot(min_IQR$min_mean,min_IQR$IQR_vec_ori,xlab="min(average(log2CPM))",ylab="IQR",cex.lab=1.5,cex.axis=1.2,col="blue",ylim=range_y)
dev.off()


pdf(paste0(dir_fig,"min_IQR_adj.pdf"))
plot(min_IQR$min_mean,min_IQR$IQR_vec_adj,xlab="min(average(log2CPM))",ylab="IQR",cex.lab=1.5,cex.axis=1.2,col="red",ylim=range_y)
dev.off()


