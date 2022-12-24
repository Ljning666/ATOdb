BiocManager::install("org.Hs.eg.db")
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)

mtar<-read.table("mRNA_target.txt",sep="\t",header=F,as.is=T)

###GO

ego_mf <- enrichGO(gene=as.character(mtar[,1]),OrgDb="org.Hs.eg.db",keyType = "SYMBOL",ont="MF",pvalueCutoff=0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 3,readable=FALSE)
ego_cc <- enrichGO(gene=as.character(mtar[,1]),OrgDb="org.Hs.eg.db",keyType = "SYMBOL",ont="CC",pvalueCutoff=0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 3,readable=FALSE)
ego_bp <- enrichGO(gene=as.character(mtar[,1]),OrgDb="org.Hs.eg.db",keyType = "SYMBOL",ont="BP",pvalueCutoff=0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 3,readable=FALSE)

write.table(as.data.frame(ego_mf),file="ego_mf.txt",row.names =F,quote = F,sep="\t")
write.table(as.data.frame(ego_cc),file="ego_cc.txt",row.names =F,quote = F,sep="\t")
write.table(as.data.frame(ego_bp),file="ego_bp.txt",row.names =F,quote = F,sep="\t")

##KEGG，匹配不上
#x <- marker1
test = bitr(mtar[,1], #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") 
kegg<- enrichKEGG(as.character(test[,2]),keyType="kegg",organism="hsa",pvalueCutoff=0.05)
kegg<- enrichKEGG(as.character(test[,2]),keyType="kegg",pvalueCutoff=0.05)
write.table(as.data.frame(kegg),file = "kegg.txt",sep = "\t",quote = F,row.names = F)



#直接整理，第八列是基因名
mf_mat <- matrix(NA,1,4)
colnames(mf_mat) <- c("geneID",	"ID",	"Description",	"Type")
gene_list <-  strsplit(ego_mf[,8],"/")

for(i in 1:length(ego_mf[,8])){
  for(j in 1:length(gene_list[[i]])){
    new_mat <- c(gene_list[[i]][j],ego_mf[i,1],ego_mf[i,2],"MF")
    #names(mf_mat) <- c("geneID",	"ID",	"Description",	"Type")
    mf_mat <- rbind(mf_mat,new_mat)
  }
}
mf_mat <- mf_mat[-1,]


cc_mat <- matrix(NA,1,4)
colnames(cc_mat) <- c("geneID",	"ID",	"Description",	"Type")
gene_list <-  strsplit(ego_cc[,8],"/")

for(i in 1:length(ego_cc[,8])){
  for(j in 1:length(gene_list[[i]])){
    new_mat <- c(gene_list[[i]][j],ego_cc[i,1],ego_cc[i,2],"CC")
    #names(mf_mat) <- c("geneID",	"ID",	"Description",	"Type")
    cc_mat <- rbind(cc_mat,new_mat)
  }
}
cc_mat <- cc_mat[-1,]


bp_mat <- matrix(NA,1,4)
colnames(bp_mat) <- c("geneID",	"ID",	"Description",	"Type")
gene_list <-  strsplit(ego_bp[,8],"/")

for(i in 1:length(ego_bp[,8])){
  for(j in 1:length(gene_list[[i]])){
    new_mat <- c(gene_list[[i]][j],ego_bp[i,1],ego_bp[i,2],"BP")
    #names(mf_mat) <- c("geneID",	"ID",	"Description",	"Type")
    bp_mat <- rbind(bp_mat,new_mat)
  }
}
bp_mat <- bp_mat[-1,]

kegg_mat <- matrix(NA,1,3)
colnames(kegg_mat) <- c("geneID",	"ID",	"Description")
gene_list <-  strsplit(kegg[,8],"/")

for(i in 1:length(kegg[,8])){
  for(j in 1:length(gene_list[[i]])){
    new_mat <- c(gene_list[[i]][j],kegg[i,1],kegg[i,2])
    #names(mf_mat) <- c("geneID",	"ID",	"Description",	"Type")
    kegg_mat <- rbind(kegg_mat,new_mat)
  }
}
kegg_mat <- kegg_mat[-1,]


#排序
go_result <- rbind(bp_mat,cc_mat,mf_mat)
new_go_result <-go_result[order(go_result[,1]),] 
write.table(new_go_result,"go_result_final.txt",sep="\t",row.names = F,col.names = T,quote=F)

#整理kegg结果
gene.df <- bitr(kegg_mat[,1], fromType = "ENTREZID", #fromType是指你的数据ID类型是属于哪一类的
                toType ="SYMBOL", #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = org.Hs.eg.db)#Orgdb是指对应的注释包是哪个

indall <- NULL
for(i in 1:length(kegg_mat[,1])){
  ind <- which(gene.df[,1]==kegg_mat[i,1])
  indall <- append(indall,ind)
}

newdata <- cbind(gene.df[indall,2],kegg_mat[,c(2,3)])
newdata <- newdata[order(newdata[,1]),] 
names(newdata) <- c("geneID","ID","Description")
write.table(newdata,"kegg_result_final.txt",sep="\t",row.names = F,col.names = T,quote=F)

