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
#length(union(drug_drug_target_num[,1],drug_drug_target_num[,2])) 只有99个结点，1866条关系
write.table(drug_drug_target_num2,file="drug_drug_target_num2.txt",sep="\t",row.names=F,quote=F)

node_att3 <- node_att
for(i in 1:length(node_att3[,1])){
  ind <- which(toupper(drug_name)==row.names(node_att3)[i])
  row.names(node_att3)[i] <- drug_name[ind]
}
write.table(node_att3,file="node_att3.txt",sep="\t",row.names=T,quote=F)