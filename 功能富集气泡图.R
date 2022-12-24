library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr)
library(stringr)
data_918 <- read.csv("C:/Users/ttt/Desktop/ttt/data918.csv",header=T,as.is=T)
#qipao+boxplot###################################################################################
all_mrna_target <- list()
list_name <- c("ATO_induce_all","ATO_induce_up","ATO_induce_down",
               "ATO_inhibit_all","ATO_inhibit_up","ATO_inhibit_down",
               "syn_induce_all","syn_induce_up","syn_induce_down",
               "syn_inhibit_all","syn_inhibit_up","syn_inhibit_down",
               "deto_induce_all","deto_induce_up","deto_induce_down",
               "deto_inhibit_all","deto_inhibit_up","deto_inhibit_down")
up_regulate <- c("Dephosphorylation,Upregulation","Methylation,Upregulation",
                 "Nuclear translocation,Upregulation","Upregulation",
                 "Upregulation,Downregulation","Upregulation,Phosphorylation")
down_regulate <- c("Downregulation","Downregulation,Dephosphorylation",
                   "Phosphorylation,Downregulation","Upregulation,Downregulation")

#去掉target为“-”的记录
new_data_918 <- subset(data_918,Target!="-")
#591-344-273,有些靶点是”Upregulation,Downregulation“，有些在不同记录中有”up..“和”down..“
all_mrna_target[[1]] <- unique(subset(new_data_918,Drug.Combination=="-"&Induce.Inhibit=="Induce")[,7]) 
all_mrna_target[[2]] <- unique(subset(new_data_918,Drug.Combination=="-"&Induce.Inhibit=="Induce"&Regulation%in%up_regulate)[,7])
all_mrna_target[[3]] <- unique(subset(new_data_918,Drug.Combination=="-"&Induce.Inhibit=="Induce"&Regulation%in%down_regulate)[,7])

#1674 1137  698
all_mrna_target[[4]] <- unique(subset(new_data_918,Drug.Combination=="-"&Induce.Inhibit=="Inhibit")[,7]) 
all_mrna_target[[5]] <- unique(subset(new_data_918,Drug.Combination=="-"&Induce.Inhibit=="Inhibit"&Regulation%in%up_regulate)[,7])
all_mrna_target[[6]] <- unique(subset(new_data_918,Drug.Combination=="-"&Induce.Inhibit=="Inhibit"&Regulation%in%down_regulate)[,7])

#61   42   19
all_mrna_target[[7]] <- unique(subset(new_data_918,Effect=="synergy"&Induce.Inhibit=="Induce")[,7]) 
all_mrna_target[[8]] <- unique(subset(new_data_918,Effect=="synergy"&Induce.Inhibit=="Induce"&Regulation%in%up_regulate)[,7])
all_mrna_target[[9]] <- unique(subset(new_data_918,Effect=="synergy"&Induce.Inhibit=="Induce"&Regulation%in%down_regulate)[,7])

#358  183  223
all_mrna_target[[10]] <- unique(subset(new_data_918,Effect=="synergy"&Induce.Inhibit=="Inhibit")[,7]) 
all_mrna_target[[11]] <- unique(subset(new_data_918,Effect=="synergy"&Induce.Inhibit=="Inhibit"&Regulation%in%up_regulate)[,7])
all_mrna_target[[12]] <- unique(subset(new_data_918,Effect=="synergy"&Induce.Inhibit=="Inhibit"&Regulation%in%down_regulate)[,7])

#146    89   63
all_mrna_target[[13]] <- unique(subset(new_data_918,Effect=="detoxicity"&Induce.Inhibit=="Induce")[,7]) 
all_mrna_target[[14]] <- unique(subset(new_data_918,Effect=="detoxicity"&Induce.Inhibit=="Induce"&Regulation%in%up_regulate)[,7])
all_mrna_target[[15]] <- unique(subset(new_data_918,Effect=="detoxicity"&Induce.Inhibit=="Induce"&Regulation%in%down_regulate)[,7])

#5    3    2
all_mrna_target[[16]] <- unique(subset(new_data_918,Effect=="detoxicity"&Induce.Inhibit=="Inhibit")[,7]) 
all_mrna_target[[17]] <- unique(subset(new_data_918,Effect=="detoxicity"&Induce.Inhibit=="Inhibit"&Regulation%in%up_regulate)[,7])
all_mrna_target[[18]] <- unique(subset(new_data_918,Effect=="detoxicity"&Induce.Inhibit=="Inhibit"&Regulation%in%down_regulate)[,7])

sapply(all_mrna_target,length)
names(all_mrna_target) <- list_name

#可以sapply生成所有文件
mtar <- all_mrna_target[[1]]

all_mtarget <- read.table("totalmRNA_ENSG.txt",header=T,as.is=T)
true_mrna_target <- all_mrna_target
for(i in 1:18){
  index <- which(true_mrna_target[[i]] %in% all_mtarget[,1])
  true_mrna_target[[i]] <-true_mrna_target[[i]][index]
}

#qipqo###############################
setwd("C:/Users/ttt/Desktop/qipao2")
qipao_name <- c("ATO_induce","ATO_inhibit","syn_induce","syn_inhibit","deto_induce","deto_inhibit")
sapply(true_mrna_target,length) #deto_inhibit最少，syn_induce也少

draw_mf <- function(name,i){
  index <- c(1,2,3)+3*(i-1)
  mf_all <- enrichGO(gene=true_mrna_target[[index[1]]],OrgDb="org.Hs.eg.db",keyType = "SYMBOL",ont="MF",pvalueCutoff=0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 3,readable=FALSE)
  mf_up <- enrichGO(gene=true_mrna_target[[index[2]]],OrgDb="org.Hs.eg.db",keyType = "SYMBOL",ont="MF",pvalueCutoff=0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 3,readable=FALSE)
  mf_down <- enrichGO(gene=true_mrna_target[[index[3]]],OrgDb="org.Hs.eg.db",keyType = "SYMBOL",ont="MF",pvalueCutoff=0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 3,readable=FALSE)
  
  if(nrow(mf_all)>10){
    mf_all2 <- as.data.frame(mf_all)[1:10,c(2,6,9)]
  }else{
    mf_all2 <- as.data.frame(mf_all)[,c(2,6,9)]
  }
  if(nrow(mf_up)>10){
    mf_up2 <- as.data.frame(mf_up)[1:10,c(2,6,9)]
  }else{
    mf_up2 <- as.data.frame(mf_up)[,c(2,6,9)]
  }
  if(nrow(mf_down)>10){
    mf_down2 <- as.data.frame(mf_down)[1:10,c(2,6,9)]
  }else{
    mf_down2 <- as.data.frame(mf_down)[,c(2,6,9)]
  }
  
  mf_up2$Type <- "Upregulation"
  mf_all2$Type <- "All"
  mf_down2$Type <- "Downregulation"
  #重复的rownames自动加上12的后缀
  
  mf_all_up_down <- rbind(mf_all2,mf_up2,mf_down2)
  names(mf_all_up_down) <- c("Term","FDR","Count","Type")
  
  p = ggplot(mf_all_up_down,aes(x=Type,y=reorder(Term,Count)))
  p=p + geom_point()  
  # 修稿点的大小
  p=p + geom_point(aes(size=Count))
  # 展示三维数据
  pbubble = p+ geom_point(aes(size=Count,color=-1*log10(FDR)))
  # 设置渐变色
  pr = pbubble+scale_color_gradient(low="#FBEA2F",high = "#14A48F")
  # 绘制p气泡图
  pr = pr+labs(color=expression(-log[10](FDR)),size="Count",  
               x="Type",y="Pathway name",title="")+
    scale_y_discrete(labels=function(x) str_wrap(x, width=50))+
    theme_bw()+
    theme(axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12), #, face =  "bold"
          legend.title = element_text(size = 13), 
          legend.text = element_text(size = 13),
          axis.title = element_text(size=15)
    )
  ggsave(paste(name,"_mf.pdf"),width=12,height=8)
}
draw_mf(qipao_name[1],1)
draw_mf(qipao_name[2],2)
draw_mf(qipao_name[3],3)
draw_mf(qipao_name[4],4)
draw_mf(qipao_name[5],5)
draw_mf(qipao_name[6],6)



draw_cc <- function(name,i){
  index <- c(1,2,3)+3*(i-1)
  cc_all <- enrichGO(gene=true_mrna_target[[index[1]]],OrgDb="org.Hs.eg.db",keyType = "SYMBOL",ont="CC",pvalueCutoff=0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 3,readable=FALSE)
  cc_up <- enrichGO(gene=true_mrna_target[[index[2]]],OrgDb="org.Hs.eg.db",keyType = "SYMBOL",ont="CC",pvalueCutoff=0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 3,readable=FALSE)
  cc_down <- enrichGO(gene=true_mrna_target[[index[3]]],OrgDb="org.Hs.eg.db",keyType = "SYMBOL",ont="CC",pvalueCutoff=0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 3,readable=FALSE)
  
  if(nrow(cc_all)>10){
    cc_all2 <- as.data.frame(cc_all)[1:10,c(2,6,9)]
  }else{
    cc_all2 <- as.data.frame(cc_all)[,c(2,6,9)]
  }
  if(nrow(cc_up)>10){
    cc_up2 <- as.data.frame(cc_up)[1:10,c(2,6,9)]
  }else{
    cc_up2 <- as.data.frame(cc_up)[,c(2,6,9)]
  }
  if(nrow(cc_down)>10){
    cc_down2 <- as.data.frame(cc_down)[1:10,c(2,6,9)]
  }else{
    cc_down2 <- as.data.frame(cc_down)[,c(2,6,9)]
  }
  
  cc_up2$Type <- "Upregulation"
  cc_all2$Type <- "All"
  cc_down2$Type <- "Downregulation"
  #重复的rownames自动加上12的后缀
  
  cc_all_up_down <- rbind(cc_all2,cc_up2,cc_down2)
  names(cc_all_up_down) <- c("Term","FDR","Count","Type")
  
  p = ggplot(cc_all_up_down,aes(x=Type,y=reorder(Term,Count)))
  p=p + geom_point()  
  # 修稿点的大小
  p=p + geom_point(aes(size=Count))
  # 展示三维数据
  pbubble = p+ geom_point(aes(size=Count,color=-1*log10(FDR)))
  # 设置渐变色
  pr = pbubble+scale_color_gradient(low="#FBEA2F",high = "#14A48F")
  # 绘制p气泡图
  pr = pr+labs(color=expression(-log[10](FDR)),size="Count",  
               x="Type",y="Pathway name",title="")+
    scale_y_discrete(labels=function(x) str_wrap(x, width=60))+
    theme_bw()+
    theme(axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12), #, face =  "bold"
          legend.title = element_text(size = 13), 
          legend.text = element_text(size = 13),
          axis.title = element_text(size=15)
    )
  ggsave(paste(name,"_cc.pdf"),width=12,height=8)
}

draw_cc(qipao_name[1],1)
draw_cc(qipao_name[2],2)
draw_cc(qipao_name[3],3)
draw_cc(qipao_name[4],4)
draw_cc(qipao_name[5],5)
draw_cc(qipao_name[6],6)

draw_bp <- function(name,i){
  index <- c(1,2,3)+3*(i-1)
  bp_all <- enrichGO(gene=true_mrna_target[[index[1]]],OrgDb="org.Hs.eg.db",keyType = "SYMBOL",ont="BP",pvalueCutoff=0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 3,readable=FALSE)
  bp_up <- enrichGO(gene=true_mrna_target[[index[2]]],OrgDb="org.Hs.eg.db",keyType = "SYMBOL",ont="BP",pvalueCutoff=0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 3,readable=FALSE)
  bp_down <- enrichGO(gene=true_mrna_target[[index[3]]],OrgDb="org.Hs.eg.db",keyType = "SYMBOL",ont="BP",pvalueCutoff=0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 3,readable=FALSE)
  
  if(nrow(bp_all)>10){
    bp_all2 <- as.data.frame(bp_all)[1:10,c(2,6,9)]
  }else{
    bp_all2 <- as.data.frame(bp_all)[,c(2,6,9)]
  }
  if(nrow(bp_up)>10){
    bp_up2 <- as.data.frame(bp_up)[1:10,c(2,6,9)]
  }else{
    bp_up2 <- as.data.frame(bp_up)[,c(2,6,9)]
  }
  if(nrow(bp_down)>10){
    bp_down2 <- as.data.frame(bp_down)[1:10,c(2,6,9)]
  }else{
    bp_down2 <- as.data.frame(bp_down)[,c(2,6,9)]
  }
  
  bp_up2$Type <- "Upregulation"
  bp_all2$Type <- "All"
  bp_down2$Type <- "Downregulation"
  #重复的rownames自动加上12的后缀
  
  bp_all_up_down <- rbind(bp_all2,bp_up2,bp_down2)
  names(bp_all_up_down) <- c("Term","FDR","Count","Type")
  
  p = ggplot(bp_all_up_down,aes(x=Type,y=reorder(Term,Count)))
  p=p + geom_point()  
  # 修稿点的大小
  p=p + geom_point(aes(size=Count))
  # 展示三维数据
  pbubble = p+ geom_point(aes(size=Count,color=-1*log10(FDR)))
  # 设置渐变色
  pr = pbubble+scale_color_gradient(low="#FBEA2F",high = "#14A48F")
  # 绘制p气泡图
  pr = pr+labs(color=expression(-log[10](FDR)),size="Count",  
               x="Type",y="Pathway name",title="")+
    scale_y_discrete(labels=function(x) str_wrap(x, width=60))+
    theme_bw()+
    theme(axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12), #, face =  "bold"
          legend.title = element_text(size = 13), 
          legend.text = element_text(size = 13),
          axis.title = element_text(size=15)
    )
  ggsave(paste(name,"_bp.pdf"),width=12,height=8)
}
draw_bp(qipao_name[1],1)
draw_bp(qipao_name[2],2)
draw_bp(qipao_name[3],3)
draw_bp(qipao_name[4],4)
draw_bp(qipao_name[5],5)
draw_bp(qipao_name[6],6)

draw_kegg <- function(name,i){
  index <- c(1,2,3)+3*(i-1)
  all = bitr(true_mrna_target[[index[1]]], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") 
  up = bitr(true_mrna_target[[index[2]]], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") 
  down = bitr(true_mrna_target[[index[3]]],fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") 
  kegg_all<- enrichKEGG(as.character(all[,2]),keyType="kegg",organism="hsa",pvalueCutoff=0.05)
  kegg_up<- enrichKEGG(as.character(up[,2]),keyType="kegg",organism="hsa",pvalueCutoff=0.05)
  kegg_down<- enrichKEGG(as.character(down[,2]),keyType="kegg",organism="hsa",pvalueCutoff=0.05)
  
  if(nrow(kegg_all)>10){
    kegg_all2 <- as.data.frame(kegg_all)[1:10,c(2,6,9)]
  }else{
    kegg_all2 <- as.data.frame(kegg_all)[,c(2,6,9)]
  }
  if(nrow(kegg_up)>10){
    kegg_up2 <- as.data.frame(kegg_up)[1:10,c(2,6,9)]
  }else{
    kegg_up2 <- as.data.frame(kegg_up)[,c(2,6,9)]
  }
  if(nrow(kegg_down)>10){
    kegg_down2 <- as.data.frame(kegg_down)[1:10,c(2,6,9)]
  }else{
    kegg_down2 <- as.data.frame(kegg_down)[,c(2,6,9)]
  }
  
  kegg_up2$Type <- "Upregulation"
  kegg_all2$Type <- "All"
  kegg_down2$Type <- "Downregulation"
  
  #重复的rownames自动加上12的后缀
  
  kegg_all_up_down <- rbind(kegg_all2,kegg_up2,kegg_down2)
  names(kegg_all_up_down) <- c("Term","FDR","Count","Type")
  
  p = ggplot(kegg_all_up_down,aes(x=Type,y=reorder(Term,Count)))
  p=p + geom_point()  
  # 修稿点的大小
  p=p + geom_point(aes(size=Count))
  # 展示三维数据
  pbubble = p+ geom_point(aes(size=Count,color=-1*log10(FDR)))
  # 设置渐变色
  pr = pbubble+scale_color_gradient(low="#FBEA2F",high = "#14A48F")
  # 绘制p气泡图
  pr = pr+labs(color=expression(-log[10](FDR)),size="Count",  
               x="Type",y="Pathway name",title="")+
    scale_y_discrete(labels=function(x) str_wrap(x, width=45))+
    theme_bw()+
    theme(axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12), #, face =  "bold"
          legend.title = element_text(size = 13), 
          legend.text = element_text(size = 13),
          axis.title = element_text(size=15)
    )
  ggsave(paste(name,"_kegg.pdf"),width=12,height=8)
}
draw_kegg(qipao_name[1],1)
draw_kegg(qipao_name[2],2)
draw_kegg(qipao_name[3],3)
draw_kegg(qipao_name[4],4)
draw_kegg(qipao_name[5],5)
draw_kegg(qipao_name[6],6)