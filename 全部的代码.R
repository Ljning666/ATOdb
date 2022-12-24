library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr)
library(stringr)
#drug-drug############################################################
data_918 <- read.csv("C:/Users/ttt/Desktop/ttt/data918.csv",header=T,as.is=T)
drug <- unique(toupper(data_918[,2])) #有些对不上，转大写
drug <- drug[-which(drug=="-")] #280个药
drug <- gsub(" ","",drug)
write.table(drug,"all-drug.txt",sep="\t",quote=F,row.names = F) #2640条

drug_tar <- read.table("C:/Users/ttt/Desktop/ttt/drug_drug/tar_drugin.txt",sep="\t",header = T,as.is=T,quote="",fill=NA)#29000
drug_tar[,2] <- toupper(drug_tar[,2])
drug_index <- which(drug_tar[,2] %in% drug)
drug_tar2 <- drug_tar[drug_index,c(1,2)] #全转大写为1952条
uni_drug <- unique(drug_tar2[,2]) #76个药
#length(which(drug%in%drug_tar[,2])) 76
#write.table(drug_tar2,"target-drug.txt",sep="\t",quote=F,row.names = F)
names(drug_tar2) <- c("Target","Drug")

drug_tar_HD <- read.table("C:/Users/ttt/Desktop/ttt/drug_drug/human_drug_target_TTD.txt",sep="\t",header = T,as.is=T,quote="",fill=NA)#40009
drug_tar_HD <- unique(drug_tar_HD[,c(2,4)]) #39804条
drug_tar_HD[,2] <- toupper(drug_tar_HD[,2]) 
drug_index_HD <- which(drug_tar_HD[,2] %in% drug) 
drug_tar_HD_2 <- drug_tar_HD[drug_index_HD,c(1,2)] #133条关系
uni_drug_HD <- unique(drug_tar_HD_2[,2]) #65个药，76+65=94
names(drug_tar_HD_2) <- c("Target","Drug")
drug_12 <- rbind(drug_tar2,drug_tar_HD_2) #unique 2023条

drug_tar_IN <- read.table("C:/Users/ttt/Desktop/ttt/drug_drug/interactions.tsv",sep="\t",header = T,as.is=T,fill=NA,quote="")#100273
drug_tar_IN <- unique(drug_tar_IN[,c(1,8)]) #59010
drug_index_IN <- which(drug_tar_IN[,2] %in% drug) 
drug_tar_IN_2 <- drug_tar_IN[drug_index_IN,c(1,2)] #2649条关系,85个药
uni_drug_IN <- unique(drug_tar_IN_2[,2])
#union(union(uni_drug,uni_drug_HD),uni_drug_IN) #一共107个药
names(drug_tar_IN_2) <- c("Target","Drug")
drug_123 <- rbind(drug_12,drug_tar_IN_2) 
drug_123 <- unique(drug_123) #unique 3188条
drug_123 <- drug_123[order(drug_123[,2]),]
drug_123 <- drug_123[-which(drug_123[,1]==""),] #有些基因名为空
write.table(drug_123,"target-drug.txt",sep="\t",quote=F,row.names = F) #3188条

uni_drug_all <- unique(drug_123[,2])
#生成list
drug_tar_list <- list()
for(i in 1:length(uni_drug_all)){
  index <- which(drug_123[,2]==uni_drug_all[i])
  drug_tar_list[[i]] <- drug_123[index,1]
}
names(drug_tar_list) <- uni_drug_all

#算共享靶点数
index_mat <- combn(107,2)
drug_drug_target_num <- matrix(NA,1,3)
colnames(drug_drug_target_num) <- c("Drug1","Drug2","number")

for(i in 1:length(index_mat[1,])){
  a <- index_mat[1,i]
  b <- index_mat[2,i]
  intersect_num <- length(intersect(drug_tar_list[[a]],drug_tar_list[[b]]))
  if(intersect_num>0){
    drug_drug_target_num <- rbind(drug_drug_target_num,c(names(drug_tar_list)[a],names(drug_tar_list)[b],intersect_num))
  }
}
drug_drug_target_num <- drug_drug_target_num[-1,] #去掉第一行NA,2209条关系
write.table(drug_drug_target_num,file="drug_drug_target_num.txt",sep="\t",row.names=F,quote=F)

#生成结点属性文件，只包括63个匹配到靶点的药物
#ind <- which(data_918[,2]==" verteporfin")
#data_918[ind,2] <- "VERTEPORFIN" 
drug_f_index <- which(toupper(data_918[,2]) %in% uni_drug_all)
drug_function <- unique(data_918[drug_f_index,c(2,6)])
ind2 <- which(drug_function[,1]=="cisplatin")
drug_function[ind2,1] <- "Cisplatin"
#sort((toupper(unique(drug_function[,1])))) 105个药，CISPLATIN重复了,VERTEPORFIN 原本有空格没对上,修改后为106个
#dim(drug_function) #196条关系
#length(unique(drug_function[,1])) #107
drug_fun_sort <- drug_function[order(drug_function[,1]),]
drug_fun_sort <- unique(drug_fun_sort)
write.table(drug_fun_sort,file="drug_att.txt",sep="\t",row.names=F,quote=F)
#write.table(sort(uni_drug),file="node_att.txt",sep="\t",row.names=F,quote=F)

function_system <- read.table("C:/Users/ttt/Desktop/function_system.txt",sep="\t",header = T,as.is=T)
node_att <- matrix(0,107,11)
row.names(node_att) <- uni_drug_all

for(i in 1:length(drug_fun_sort[,1])){
  drug_index <- which(uni_drug_all==toupper(drug_fun_sort[i,1]))
  function_index <- which(function_system[,1]==drug_fun_sort[i,2])
  system_index <- function_system[function_index,2]
  node_att[drug_index,system_index] <- 1
}

write.table(node_att,file="node_att2.txt",sep="\t",row.names=T,quote=F)

#根据drug_fun_sort改一下drug_drug_target_num和node_att里面的药物名称
drug_drug_target_num2 <- drug_drug_target_num
drug_name <- unique(drug_function[,1])
for(i in 1:length(drug_drug_target_num2[,1])){
  ind1 <- which(toupper(drug_name)==drug_drug_target_num2[i,1])
  ind2 <- which(toupper(drug_name)==drug_drug_target_num2[i,2])
  drug_drug_target_num2[i,1] <- drug_name[ind1]
  drug_drug_target_num2[i,2] <- drug_name[ind2]
}
#length(union(drug_drug_target_num[,1],drug_drug_target_num[,2])) 
write.table(drug_drug_target_num2,file="drug_drug_target_num2.txt",sep="\t",row.names=F,quote=F)

node_att3 <- node_att
for(i in 1:length(node_att3[,1])){
  ind <- which(toupper(drug_name)==row.names(node_att3)[i])
  row.names(node_att3)[i] <- drug_name[ind]
}
write.table(node_att3,file="node_att3.txt",sep="\t",row.names=T,quote=F)

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


#boxplot-fantom################################################################
fantom_data<-read.table("/Users/ttt/Desktop/ttt/GTEx_FANTOM/FANTOM_mRNA_cpm_34.txt",sep = "\t",header = T,stringsAsFactors = F) #23615
fantom_data[,35]<-rownames(fantom_data)
colnames(fantom_data)[35]<-'gene'
gene.df <- bitr(fantom_data[,35], fromType = "ENSEMBL",toType = c("SYMBOL"), OrgDb = org.Hs.eg.db) #22%没map上,剩18502个，有CATG开头的
new_fantom_data <- merge(fantom_data,gene.df,by.x="gene",by.y="ENSEMBL") #merge后剩18502行
#sum(duplicated(new_fantom_data[,36])) 有22个重复的gene name

#先把重复的gene 去掉
duplicate_index <- which(duplicated(new_fantom_data[,36]))
new_fantom_data2 <- new_fantom_data[-duplicate_index,] #剩18480
row.names(new_fantom_data2) <- new_fantom_data2[,36]
new_fantom_data2 <- new_fantom_data2[,c(2:35)]

frameto3frame <- function(a){
  row <- nrow(a)
  col <- ncol(a)
  row_name <- row.names(a)
  col_name <- colnames(a)
  V1 <- rep(row_name,time=col)
  V2 <- rep(col_name,each=row)
  V3 <- matrix(as.matrix(a),ncol=1)
  result <- as.data.frame(cbind(V1,V2,V3),stringsAsFactors = F)
  result[,3]<-as.numeric(result[,3])
  result
}


draw_boxplot <- function(x){
  #将数据集分成 mrna-target和非mtar两部分
  mtar_fantom_data <- new_fantom_data2[(which(row.names(new_fantom_data2)%in%x)),] 
  rest_fantom_data <- new_fantom_data2[-(which(row.names(new_fantom_data2)%in%x)),]
  print(dim(mtar_fantom_data))
  
  #在非tar随机抽取与mrna-target数量相等的
  #set.seed(1) 
  random_index <- sample(1:length(rest_fantom_data[,1]),length(mtar_fantom_data[,1]),replace = F)
  random_fantom_data <- rest_fantom_data[random_index,]
  
  #利用函数，整理为画图的格式
  new_mtar_fantom_data <- frameto3frame(mtar_fantom_data)
  new_mtar_fantom_data$type=rep("target",length(mtar_fantom_data[,1]))
  colnames(new_mtar_fantom_data)<-c("gene","tissue","expression","type")
  
  new_random_fantom_data <- frameto3frame(random_fantom_data)
  new_random_fantom_data$type=rep("random",length(random_fantom_data[,1]))
  colnames(new_random_fantom_data)<-c("gene","tissue","expression","type")
  all_data <- rbind(new_mtar_fantom_data,new_random_fantom_data)
  all_data$expression<-log2(all_data$expression) 
  
  #去掉log2后产生的-Inf值(log2(0)=-inf)
  aa<-ggplot(all_data[!is.infinite(all_data[,3]),],aes(x =tissue,y=expression,fill =type)) + 
    geom_boxplot(outlier.shape = 21,outlier.alpha =0) + theme_bw() + 
    labs(x = "Tissue", y = "Expression") +
    theme(axis.ticks.x = element_blank()) +
    scale_fill_manual(values = c("#94B300","#FFD36B")) +
    theme_classic()+
    theme(axis.text.x = element_text(size=11,angle=45,hjust=0.8, face =  "bold"),
          axis.text.y = element_text(size=11, face =  "bold"),
          axis.line = element_line(size=0.5),
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 14),
          axis.title = element_text(size=16)
    )+
    stat_compare_means(aes(group =type),label = "p.signif",cex=4)
  aa
}

#用循环画总是空图，手动画1~18张
setwd("C:/Users/ttt/Desktop/fantom")
pdf(file = paste0(names(all_mrna_target[1]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot(all_mrna_target[[1]])
dev.off()
pdf(file = paste0(names(all_mrna_target[2]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot(all_mrna_target[[2]])
dev.off()
pdf(file = paste0(names(all_mrna_target[3]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot(all_mrna_target[[3]])
dev.off()
pdf(file = paste0(names(all_mrna_target[4]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot(all_mrna_target[[4]])
dev.off()
pdf(file = paste0(names(all_mrna_target[5]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot(all_mrna_target[[5]])
dev.off()
pdf(file = paste0(names(all_mrna_target[6]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot(all_mrna_target[[6]])
dev.off()
pdf(file = paste0(names(all_mrna_target[7]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot(all_mrna_target[[7]])
dev.off()
pdf(file = paste0(names(all_mrna_target[8]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot(all_mrna_target[[8]])
dev.off()
pdf(file = paste0(names(all_mrna_target[9]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot(all_mrna_target[[9]])
dev.off()
pdf(file = paste0(names(all_mrna_target[10]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot(all_mrna_target[[10]])
dev.off()
pdf(file = paste0(names(all_mrna_target[11]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot(all_mrna_target[[11]])
dev.off()
pdf(file = paste0(names(all_mrna_target[12]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot(all_mrna_target[[12]])
dev.off()
pdf(file = paste0(names(all_mrna_target[13]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot(all_mrna_target[[13]])
dev.off()
pdf(file = paste0(names(all_mrna_target[14]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot(all_mrna_target[[14]])
dev.off()
pdf(file = paste0(names(all_mrna_target[15]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot(all_mrna_target[[15]])
dev.off()
pdf(file = paste0(names(all_mrna_target[16]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot(all_mrna_target[[16]])
dev.off()
pdf(file = paste0(names(all_mrna_target[17]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot(all_mrna_target[[17]])
dev.off()
pdf(file = paste0(names(all_mrna_target[18]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot(all_mrna_target[[18]])
dev.off()

#gtex############################################################################
gtex_data<-read.table("/Users/ttt/Desktop/ttt/GTEx_FANTOM/GTEx_GeneExp.txt",sep = "\t",header = T,stringsAsFactors = F) #56238
gene.df2 <- bitr(gtex_data[,1], fromType = "ENSEMBL",toType = c("SYMBOL"), OrgDb = org.Hs.eg.db) #52%没map上,剩27097个
new_gtex_data <- merge(gtex_data,gene.df2,by.x="Gene",by.y="ENSEMBL") #merge后剩27097行
#sum(duplicated(new_gtex_data[,32])) 有54个重复的gene name
#先把重复的gene 去掉
duplicate_index2 <- which(duplicated(new_gtex_data[,32]))
new_gtex_data2 <- new_gtex_data[-duplicate_index2,] #剩27043
row.names(new_gtex_data2) <- new_gtex_data2[,32]
new_gtex_data2 <- new_gtex_data2[,c(2:31)]

draw_boxplot_gtex <- function(x){
  #将数据集分成 mrna-target和非mtar两部分
  mtar_gtex_data <- new_gtex_data2[(which(row.names(new_gtex_data2)%in%x)),] 
  rest_gtex_data <- new_gtex_data2[-(which(row.names(new_gtex_data2)%in%x)),]
  print(dim(mtar_gtex_data))
  
  #在非tar随机抽取与mrna-target数量相等的
  #set.seed(1) 
  random_index <- sample(1:length(rest_gtex_data[,1]),length(mtar_gtex_data[,1]),replace = F)
  random_gtex_data <- rest_gtex_data[random_index,]
  
  #利用函数，整理为画图的格式
  new_mtar_gtex_data <- frameto3frame(mtar_gtex_data)
  new_mtar_gtex_data$type=rep("target",length(mtar_gtex_data[,1]))
  colnames(new_mtar_gtex_data)<-c("gene","tissue","expression","type")
  
  new_random_gtex_data <- frameto3frame(random_gtex_data)
  new_random_gtex_data$type=rep("random",length(random_gtex_data[,1]))
  colnames(new_random_gtex_data)<-c("gene","tissue","expression","type")
  all_data <- rbind(new_mtar_gtex_data,new_random_gtex_data)
  all_data$expression<-log2(all_data$expression) 
  
  #去掉log2后产生的-Inf值(log2(0)=-inf)
  aa<-ggplot(all_data[!is.infinite(all_data[,3]),],aes(x =tissue,y=expression,fill =type)) + 
    geom_boxplot(outlier.shape = 21,outlier.alpha =0) + theme_bw() + 
    labs(x = "Tissue", y = "Expression") +
    theme(axis.ticks.x = element_blank()) +
    scale_fill_manual(values = c("#94B300","#FFD36B")) +
    theme_classic()+
    theme(axis.text.x = element_text(size=11,angle=45,hjust=0.8, face =  "bold"),
          axis.text.y = element_text(size=11, face =  "bold"),
          axis.line = element_line(size=0.5),
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 14),
          axis.title = element_text(size=16)
    )+
    stat_compare_means(aes(group =type),label = "p.signif",cex=4)
  aa
}

#用循环画总是空图，手动画1~18张
setwd("C:/Users/ttt/Desktop/boxplot2")
pdf(file = paste0(names(all_mrna_target[1]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot_gtex(all_mrna_target[[1]])
dev.off()
pdf(file = paste0(names(all_mrna_target[2]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot_gtex(all_mrna_target[[2]])
dev.off()
pdf(file = paste0(names(all_mrna_target[3]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot_gtex(all_mrna_target[[3]])
dev.off()
pdf(file = paste0(names(all_mrna_target[4]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot_gtex(all_mrna_target[[4]])
dev.off()
pdf(file = paste0(names(all_mrna_target[5]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot_gtex(all_mrna_target[[5]])
dev.off()
pdf(file = paste0(names(all_mrna_target[6]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot_gtex(all_mrna_target[[6]])
dev.off()
pdf(file = paste0(names(all_mrna_target[7]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot_gtex(all_mrna_target[[7]])
dev.off()
pdf(file = paste0(names(all_mrna_target[8]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot_gtex(all_mrna_target[[8]])
dev.off()
pdf(file = paste0(names(all_mrna_target[9]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot_gtex(all_mrna_target[[9]])
dev.off()
pdf(file = paste0(names(all_mrna_target[10]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot_gtex(all_mrna_target[[10]])
dev.off()
pdf(file = paste0(names(all_mrna_target[11]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot_gtex(all_mrna_target[[11]])
dev.off()
pdf(file = paste0(names(all_mrna_target[12]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot_gtex(all_mrna_target[[12]])
dev.off()
pdf(file = paste0(names(all_mrna_target[13]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot_gtex(all_mrna_target[[13]])
dev.off()
pdf(file = paste0(names(all_mrna_target[14]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot_gtex(all_mrna_target[[14]])
dev.off()
pdf(file = paste0(names(all_mrna_target[15]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot_gtex(all_mrna_target[[15]])
dev.off()
pdf(file = paste0(names(all_mrna_target[16]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot_gtex(all_mrna_target[[16]])
dev.off()
pdf(file = paste0(names(all_mrna_target[17]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot_gtex(all_mrna_target[[17]])
dev.off()
pdf(file = paste0(names(all_mrna_target[18]), ".pdf"),width = 15,height=8,onefile=F)
draw_boxplot_gtex(all_mrna_target[[18]])
dev.off()

####功能堆叠图###################################################################
fun_in_drug <- data_918[,5:7]
fun_in_drug <- unique(fun_in_drug)
function_10 <- sort(table(fun_in_drug[,2]),decreasing = T)[1:10]
#Breast-cancer,Leukemia,Brain-cancer,Liver-cancer,Vascular-toxicity
#Heart-toxicity,Immune-toxicity,Lung-cancer,Cervical-cancer,Skin-toxicity 
fun_in_drug_10 <- subset(fun_in_drug,Function %in% names(function_10))

#unique(fun_in_drug_10[,c(1,2)]) 只有liver-cancer有induce和inhibit
#fun_in_drug_10[which(fun_in_drug_10[,2]=="Liver-cancer"),],induce只有一个，图上没看出来
data <- fun_in_drug_10
data$Function<-factor(data$Function,levels = names(function_10))
up<-ggplot(data, aes(Function, fill = Induce.Inhibit)) +
  geom_bar(alpha=0.5) +  theme_classic()+
  theme(axis.text.x = element_text(size = 13,angle=45,vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 13,vjust = 0.5, hjust = 0.5),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        axis.title = element_text(size=15)
  )+
  scale_fill_manual(values=c("#A4CA93","#E86E13"))
up

######协同条形###################################################
syn_drug_fun <- subset(data_918,Effect=="synergy")[,c(1,2,3,6)] #加上PMID，unique
syn_drug_fun <- unique(syn_drug_fun)
drug_10 <- sort(table(syn_drug_fun[,2]),decreasing = T)[1:10]

syn_drug_fun_10 <- subset(syn_drug_fun,Drug.Combination %in% names(drug_10))
data <- syn_drug_fun_10
#调整顺序，从高到低
data$Drug.Combination<-factor(data$Drug.Combination,levels = names(drug_10))

#function匹配系统
drug_system <- matrix(0,length(unique(data$Function)),2)
for(i in 1:length(drug_system[,2])){
  ind <- which(function_system[,1]==unique(data$Function)[i])
  drug_system[i,] <- c(unique(data$Function)[i],function_system[ind,2])
}
drug_system <- drug_system[order(as.numeric(drug_system[,2])),]
#调整图例顺序，相同系统色系同
data$Function<-factor(data$Function,levels = drug_system[,1])

up<-ggplot(data, aes(Drug.Combination, fill = Function)) +
  geom_bar(alpha=0.5) +  theme_classic()+
  theme(axis.text.x = element_text(size = 13,angle=45,vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 13,vjust = 0.5, hjust = 0.5),
        legend.title = element_text(size = 13), 
        legend.text = element_text(size = 10),
        axis.title = element_text(size=15))+
  scale_fill_manual(values=c("#7BB358","#86B768","#9CD47A","#AEEB89","#C2F8A1","#BFD1A9","#CCDFB5",#1-green
                             "#E08E4D","#E39A61","#EAA46C",#2-orange
                             "#838C9C","#98A2B3", #3-grey
                             "#E26A57", #5-red
                             "#EFCDA7","#F0DDC8",#7
                             "#94C6DC","#7AC5F3","#A7DCFC","#DEEEF5",#8-blue
                             "#ECAA43","#F6CC8b","#F2D78D", #9-yellow
                             "#BAA38F",#10
                             "#D39E5D"))
up

####减毒条形#####################################################
deto_drug_fun <- subset(data_918,Effect=="detoxicity")[,c(1,2,3,6)] #加上PMID，unique？
deto_drug_fun <- unique(deto_drug_fun)
drug_10 <- sort(table(deto_drug_fun[,2]),decreasing = T)[1:10]

deto_drug_fun_10 <- subset(deto_drug_fun,Drug.Combination %in% names(drug_10))
data <- deto_drug_fun_10

data$Drug.Combination<-factor(data$Drug.Combination,levels = names(drug_10))

#function匹配系统
drug_system <- matrix(0,length(unique(data$Function)),2)
for(i in 1:length(drug_system[,2])){
  ind <- which(function_system[,1]==unique(data$Function)[i])
  drug_system[i,] <- c(unique(data$Function)[i],function_system[ind,2])
}
drug_system[11,2] <- "8"
drug_system[10,2] <- "5"
drug_system <- drug_system[order(drug_system[,2]),]
#调整图例顺序，相同系统色系同
data$Function<-factor(data$Function,levels = drug_system[,1])


up<-ggplot(data, aes(Drug.Combination, fill = Function)) +
  geom_bar(alpha=0.5) +  theme_classic()+
  theme(axis.text.x = element_text(size = 12,angle=45,vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 12,vjust = 0.5, hjust = 0.5),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10),
        axis.title = element_text(size=14))+
  scale_fill_manual(values=c("#7BB358","#AEEB89","#C2F8A1",#1-green
                             "#E08E4D",#2-orange
                             "#9EA5B2", #3-grey
                             "#EDEF8B",#4-yellow
                             "#E26A57", #5
                             "#EFCDA7",#7
                             "#94C6DC","#A7DCFC","#D3ECF6","#DEEEF5",#8-blue
                             "#EBC889","#ECAA43","#F2D78D" #9-yellow
  )) 
up



