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
