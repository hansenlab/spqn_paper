

library(matrixStats)
library(matrixStats)
library(WGCNA)
library(ggplot2)
library(ggridges)
library(SummarizedExperiment)
library(readr)
library(stringr)


dir_data="/data/"
dir_input="/corplot_ori_svapcs/data/" # cor_est
dir_output="/TF/data/"
dir_fig="/TP_updatedPartition/fig/"
setwd(dir_fig)

### P
# conditional on exp
load(paste0(dir_fig,"df_all_P_middle.RData"))
df_all$expression_level="middle"
df_all_=df_all
load(paste0(dir_fig,"df_all_P_high.RData"))
df_all$expression_level="high"
df_all_=rbind(df_all_,df_all)
load(paste0(dir_fig,"df_all_P_low.RData"))
df_all$expression_level="low"
df_all_=rbind(df_all_,df_all)

ls(df_all_)


# overall, unconditional on exp
dir_fig_old="/users/ywang/10_25_2019/TF/fig/"

load(paste0("/users/ywang/10_25_2019/TP_Nov1_2021/fig/df_all_P_overall.RData"))
df_all$expression_level="overall"
# df_all_=df_all
# ls(df_all_)
df_all_=rbind(df_all_,df_all)
ls(df_all)
ls(df_all_)
df_all_$tissue=str_replace((df_all_$tissue),"_"," ")
df_all_$tissue=str_replace((df_all_$tissue),"_"," ")

gg_color <- function(n) {
  
  hues = seq(15, 375, length = n + 1)
  
  hcl(h = hues, l = 65, c = 100)[1:n]
  
}
# n = 3
# cols = gg_color(n)

pdf(paste0(dir_fig,"scatter_P.pdf"))
df_all_$signal_percent=as.numeric(df_all_$signal_percent)
df_all_$expression_level=factor(df_all_$expression_level,levels=c("low","middle","high","overall"))
ggplot(df_all_[as.numeric(df_all_$signal_percent)>0.3,], 
       aes(x=signal_percent, y=increase_percent))+ 
  geom_point(size=1,
             aes(colour = expression_level),alpha=0.3)+
  scale_colour_manual(name="expression level",
                    values=c(gg_color(3),"#000000"))+ 
  facet_wrap(~ tissue, ncol=3)+ scale_x_continuous(breaks=c(.6, 1.4, 2.2, 3))+
  # theme(strip.text.x = element_text(size = 7),strip.background =element_rect(fill="white"))+
  # theme_classic()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(
          color="white", fill="white", size=1.5, linetype="solid"
        ))+
  geom_hline(linetype="dashed",  aes(yintercept = 0))+
  xlab("signal threshold (in percent)") +ylab("percent increase of edges")
dev.off()

# pdf("1.pdf")
# p <- ggplot(mtcars, aes(mpg, wt)) +
#   geom_point(aes(colour = factor(cyl)))
# p + scale_colour_manual(values = c("red", "blue", "green"))
# dev.off()
pdf(paste0(dir_fig,"scatter_P_AdiposeSubcutaneous.pdf"),width=5,height=3)
# for(tissue_tmp in unique(df_all_$tissue)){
tissue_tmp="Adipose Subcutaneous"
# plot(df_all_[which(df_all_$tissue == tissue_tmp ),])
ggplot(df_all_[which(df_all_$tissue == tissue_tmp & as.numeric(df_all_$signal_percent)>0.3 ),], 
       aes(x=signal_percent, y=increase_percent))+ 
  geom_point(aes(colour = expression_level),alpha=0.5)+ ggtitle("percentage increase of positives") + 
  # scale_colour_discrete(name="expression level") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_rect(
          color="white", fill="white", size=1.5, linetype="solid"
        )) +
  xlab("signal threshold (in percent)") +ylab("percent increase of edges")+
  geom_hline(linetype="dashed",  aes(yintercept = 0))+
  scale_colour_manual(name="expression level",
                      values=c(gg_color(3),"#000000"))
# guides(fill=guide_legend(title="expression level")))
# }
dev.off()


pdf(paste0(dir_fig,"nedge_compare_alltissues_P_cond_exp_v2.pdf"),width=10,height=5)
# df_all2=df_all
df_all_$signal_percent=as.character(df_all_$signal_percent)
# df_all_$increase_percent[df_all_$increase_percent>40]=40
df_all_$expression_level=factor(df_all_$expression_level,level=c("low","middle","high","overall"))
ggplot(df_all_[as.numeric(df_all_$signal_percent)>0.3,], 
       aes(x=signal_percent, y=increase_percent,fill=expression_level,)) +
  geom_boxplot(position=position_dodge(0.8),outlier.shape = NA,alpha=0.5)+
  ggtitle("conditional on expression")+
  geom_point(aes(col=expression_level),position=position_jitterdodge(),size=.7,alpha=0.9)+
  # geom_jitter(size=.5,position=position_dodge(0.8))+ 
  scale_colour_manual(name="expression level",
                    values=c(gg_color(3),"grey"))+
  scale_fill_manual(name="expression level",
                      values=c(gg_color(3),"grey"))+
  labs(x = "percentage of signals",y="percentage increase of coexpression signals")#+ 
  # scale_colour_discrete(name="expression level")
  # guides(fill=guide_legend(title="expression level")) 
dev.off()

############ TP
# conditional on exp

load(paste0(dir_fig,"df_all_TP_middle.RData"))
df_all$expression_level="middle"
df_all_=df_all
load(paste0(dir_fig,"df_all_TP_high.RData"))
df_all$expression_level="high"
df_all_=rbind(df_all_,df_all)
load(paste0(dir_fig,"df_all_TP_low.RData"))
df_all$expression_level="low"
df_all_=rbind(df_all_,df_all)

# library(ggplot2)


# overall, unconditional on exp
load(paste0("/users/ywang/10_25_2019/TP_Nov1_2021/fig/df_all_TP_overall.RData"))
df_all$expression_level="overall"
# df_all_=df_all
# ls(df_all_)
df_all_=rbind(df_all_,df_all)

df_all_$tissue=str_replace((df_all_$tissue),"_"," ")
df_all_$tissue=str_replace((df_all_$tissue),"_"," ")


# range(df_all_$increase_percent[which(abs(df_all_$increase_percent)<1000)])
# -100  200


pdf(paste0(dir_fig,"scatter_TP.pdf"))
# for(tissue_tmp in unique(df_all_$tissue)){
# plot(df_all_[which(df_all_$tissue == tissue_tmp &),])
# }
# for(tissue_tmp in unique(df_all_$tissue)){
# unique(df_all_$signal_percent)
df_all_$signal_percent=as.numeric(df_all_$signal_percent)
df_all_$expression_level=factor(df_all_$expression_level,levels=c("low","middle","high","overall"))
ggplot(df_all_[as.numeric(df_all_$signal_percent)>0.3,], 
       aes(x=signal_percent, y=increase_percent))+ 
  geom_point(size=1,
             aes(colour = expression_level),alpha=0.3)+
  scale_colour_manual(name="expression level",
    values=c(gg_color(3),"#000000"))+ 
  facet_wrap(~ tissue, ncol=3)+ scale_x_continuous(breaks=c(.6, 1.4, 2.2, 3))+
  # theme(strip.background =element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"),
        
        strip.background = element_rect(
          color="white", fill="white", size=1.5, linetype="solid"
        ))+
  # theme(strip.text.x = element_text(size = 7),strip.background =element_rect(fill="white"))+
  # theme_classic()+
  geom_hline(linetype="dashed",  aes(yintercept = 0))+
  xlab("signal threshold (in percent)") +ylab("percent increase of true positives")

  # guides(fill=guide_legend(title="expression level"))
  # labs(fill='expression level') 
  # scale_fill_manual(legend_title,values=c("orange","red"),name="expression level")

# }
dev.off()


pdf(paste0(dir_fig,"scatter_TP_AdiposeSubcutaneous.pdf"),width=5,height=3)
# for(tissue_tmp in unique(df_all_$tissue)){
tissue_tmp="Adipose Subcutaneous"
  # plot(df_all_[which(df_all_$tissue == tissue_tmp ),])
  ggplot(df_all_[which(df_all_$tissue == tissue_tmp&as.numeric(df_all_$signal_percent)>0.3 ),], 
               aes(x=signal_percent, y=increase_percent))+ 
          geom_point(
                     aes(colour = expression_level),alpha=0.3)+ ggtitle("percentage increase of true positives") + 
          # scale_colour_discrete(name="expression level") +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                                                                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                                strip.background = element_rect(
                                                                  color="white", fill="white", size=1.5, linetype="solid"
                                                                )) +
  xlab("signal threshold (in percent)") +ylab("percent increase of true positives")+
  geom_hline(linetype="dashed",  aes(yintercept = 0))+
  scale_colour_manual(name="expression level",
                        values=c(gg_color(3),"#000000"))
  # guides(fill=guide_legend(title="expression level")))
# }
dev.off()



pdf(paste0(dir_fig,"nedge_compare_alltissues_TP_cond_exp_v2.pdf"),width=10,height=5)
# df_all2=df_all
df_all_$signal_percent=as.character(df_all_$signal_percent)
# df_all_$increase_percent[df_all_$increase_percent>40]=40
df_all_$expression_level=factor(df_all_$expression_level,level=c("low","middle","high","overall"))
# ggplot(df_all_[as.numeric(df_all_$signal_percent)>0.3,], aes(x=signal_percent, y=increase_percent,fill=expression_level)) +
#   geom_boxplot(position=position_dodge(0.8),outlier.shape = NA)+ggtitle("conditional on expression")+
#   geom_jitter(position=position_dodge(0.8))+ labs(x = "percentage of signals",y="percentage increase of true positives")+ 
#   df_all_$signal_percent=as.character(df_all_$signal_percent)
# df_all_$increase_percent[df_all_$increase_percent>40]=40
# df_all_$expression_level=factor(df_all_$expression_level,level=c("low","middle","high","overall"))
ggplot(df_all_[as.numeric(df_all_$signal_percent)>0.3,], 
       aes(x=signal_percent, y=increase_percent,fill=expression_level,)) +
  geom_boxplot(position=position_dodge(0.8),outlier.shape = NA,alpha=0.5)+
  ggtitle("conditional on expression")+
  geom_point(aes(col=expression_level),position=position_jitterdodge(),size=.7,alpha=0.9)+
  # geom_jitter(size=.5,position=position_dodge(0.8))+ 
  scale_colour_manual(name="expression level",
                      values=c(gg_color(3),"grey"))+
  scale_fill_manual(name="expression level",
                    values=c(gg_color(3),"grey"))+
  labs(x = "percentage of signals",y="percentage increase of true positives")#+ 
# scale_colour_discrete(name="expression level")
# guides(fill=guide_legend(title="expression level"))   scale_colour_discrete(name="expression level")+ guides(fill=guide_legend(title="expression level")) 
dev.off()
# scale_colour_discrete(name="expression level")  )

# strip.background = element_blank()





