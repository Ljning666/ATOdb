data_929 <- read.csv("C:/Users/ttt/Desktop/data929.csv",header=T,as.is=T)

#靶点数目统计##################################################################
ATO_all_target <- unique(subset(data_929, Species=="Homo sapiens ")[,17])
ATO_all_target <- ATO_all_target[-1] #1955(不筛人类2325)

#匹配为mrna靶点（所有药物联用的）
prcod_mrna_target <- read.table("totalmRNA_ENSG.txt",sep="\t",header=T,as.is=T) #原本是mRNATYPE.txt
ATO_all_mrna_target <- ATO_all_target[which(ATO_all_target %in% prcod_mrna_target[,1])] #1259个
#write.table(ATO_all_mrna_target,"ATO_all_mrna_target.txt",sep="\t",col.names = F,row.names = F,quote = F)

#ATO单独处理对应的mrna靶点
ATO_alone_target <- unique(subset(data_929,Drug.Combination=="-" & Species=="Homo sapiens ")[,17]) #人类的
ATO_alone_target <- ATO_alone_target[-1] #1807（含mrna外）
ATO_alone_mrna_target <- ATO_alone_target[which(ATO_alone_target %in% prcod_mrna_target[,1])] #1174个
#write.table(ATO_alone_mrna_target,"ATO_alone_mrna_target.txt",sep="\t",col.names = F,row.names = F,quote = F)


#组合药物数目统计##############################################################
ATO_comb_drug <- unique(data_929[,6])
ATO_comb_drug <- ATO_comb_drug[-1] #284种与ATO联用的药物

#miRNA靶点统计####################################################
count_mirna_target <- ATO_all_target[grep("hsa",ATO_all_target)] #剩104个
yuanben_mirna_target <- read.table("miRNA_target.txt",header=T,sep="\t")
#sort(count_mirna_target)
#sort(yuanben_mirna_target[,1]) #原本是105个miRNA，存在不是人类的物种有hsa开头的靶点（hsa-miR-133）



#靶点分类############################################################
ato_alone_mrna_target <- list()
list_name <- c("ATO_induce_all","ATO_induce_up","ATO_induce_down",
               "ATO_inhibit_all","ATO_inhibit_up","ATO_inhibit_down")
up_regulate <- c("Dephosphorylation,Upregulation","Methylation,Upregulation",
                 "Nuclear translocation,Upregulation","Upregulation",
                 "Upregulation,Downregulation","Upregulation,Phosphorylation")
down_regulate <- c("Downregulation","Downregulation,Dephosphorylation",
                   "Phosphorylation,Downregulation","Upregulation,Downregulation")

#去掉target为“-”的记录
new_data_929 <- subset(data_929,Target!="-")
#dim(unique(subset(new_data_929,Drug.Combination=="-"& Species=="Homo sapiens ")[17])) #1806与上面数目1807一致（少“-”）
#基于此进一步筛选；筛选后再过滤非mrna靶点

#299-187-120，有些靶点是”Upregulation,Downregulation“，有些在不同记录中有”up..“和”down..“
ato_alone_mrna_target[[1]] <- unique(subset(new_data_929,Drug.Combination=="-"& Species=="Homo sapiens "&Induce.Inhibit=="Induce")[,17]) 
ato_alone_mrna_target[[2]] <- unique(subset(new_data_929,Drug.Combination=="-"& Species=="Homo sapiens "&Induce.Inhibit=="Induce"&Regulation%in%up_regulate)[,17])
ato_alone_mrna_target[[3]] <- unique(subset(new_data_929,Drug.Combination=="-"& Species=="Homo sapiens "&Induce.Inhibit=="Induce"&Regulation%in%down_regulate)[,17])
#unique(subset(new_data_929,Drug.Combination=="-"& Species=="Homo sapiens "&Induce.Inhibit=="Induce")[,c(1,6,7,9,10,17,18)])

#1592-1104-635 
ato_alone_mrna_target[[4]] <- unique(subset(new_data_929,Drug.Combination=="-"& Species=="Homo sapiens "&Induce.Inhibit=="Inhibit")[,17]) 
ato_alone_mrna_target[[5]] <- unique(subset(new_data_929,Drug.Combination=="-"& Species=="Homo sapiens "&Induce.Inhibit=="Inhibit"&Regulation%in%up_regulate)[,17])
ato_alone_mrna_target[[6]] <- unique(subset(new_data_929,Drug.Combination=="-"& Species=="Homo sapiens "&Induce.Inhibit=="Inhibit"&Regulation%in%down_regulate)[,17])

sapply(ato_alone_mrna_target,length)
names(ato_alone_mrna_target) <- list_name
true_ato_alone_mrna_target <- ato_alone_mrna_target

for(i in 1:6){
  index <- which(ato_alone_mrna_target[[i]] %in% prcod_mrna_target[,1])
  true_ato_alone_mrna_target[[i]] <- true_ato_alone_mrna_target[[i]][index]
}
sapply(true_ato_alone_mrna_target,length)

#ATO_induce_all    ATO_induce_up  ATO_induce_down  ATO_inhibit_all   ATO_inhibit_up ATO_inhibit_down
#218              136               94             1013              711              377
#没差太多

#原本1179(没筛选人类)的结果
#ATO_induce_all    ATO_induce_up  ATO_induce_down  ATO_inhibit_all   ATO_inhibit_up ATO_inhibit_down
#231              145              102             1015              711            380 

#数量统计
#length(which(unique(subset(new_data_929,Drug.Combination=="-" & Species=="Homo sapiens ")[,17]) %in% prcod_mrna_target[,1]) ) 
#ATO直接靶向的mrna target有1174个
#induce共218个，136上调，94下调 
#inhibit 1013个, 711上调，377下调
#（有些记录可能induce和inhibit同时出现了）


#计算拓扑属性#####################################################################
#计算整个大网的拓扑属性然后再提取是ATO靶点的出来
library(igraph) 
ppi_network <- read.table("ppis_duplicated.txt",sep="\t",header = T,as.is=T)
ppi_network <- graph.data.frame(ppi_network[,c(1,2)], directed=F) #7937个结点

#度
all_degree <- as.matrix(degree(ppi_network))
#聚类系数分布情况
all_cluster <- as.matrix(transitivity(ppi_network,type = "local",isolates = "zero"))
row.names(all_cluster) <- names(V(ppi_network))
#中心性betweenness
all_between <- as.matrix(betweenness(ppi_network,directed = F,weights = NA))
#closeness接近中心性
all_closeness<-as.matrix(closeness(ppi_network))

all_tuopu <- cbind(all_degree,all_cluster,all_between,all_closeness)
colnames(all_tuopu) <- c("degree","cluster_coefficient","betweenness","closeness")
#和原来计算的是一致的

#重新提取
ato_alone_tuopu <- all_tuopu[which(row.names(all_tuopu) %in% ATO_alone_mrna_target),] #855个（原本所有靶点1268计算是929个）
rest_tuopu <- all_tuopu[which(!(row.names(all_tuopu) %in% ATO_alone_mrna_target)),] #7082个

#write.table(ato_alone_tuopu,file = "ato_alone_tuopu.txt",sep = "\t",quote = F,row.names = T)
#write.table(rest_tuopu,file = "rest_tuopu.txt",sep = "\t",quote = F,row.names = T)


#根据induce和inhibit 分组画拓扑属性图########################
target_data <- as.data.frame(ato_alone_tuopu)
target_data$type=rep("target",855)

rest_data <- as.data.frame(rest_tuopu)
rest_data$type=rep("others",7082)
rest_data$gene <- row.names(rest_data)
row.names(rest_data) <- 1:7082

#induce###############################
target_data_induce <- target_data[row.names(target_data) %in% true_ato_alone_mrna_target[[1]],] #166个induce
target_data_induce_up <- target_data_induce[row.names(target_data_induce) %in% true_ato_alone_mrna_target[[2]],] #107个up
target_data_induce_down <- target_data[row.names(target_data) %in% true_ato_alone_mrna_target[[3]],] #70个down

target_data_induce_up$type <- "up"
target_data_induce_down$type <- "down"

target_data_induce$gene <- row.names(target_data_induce)
target_data_induce_up$gene <- row.names(target_data_induce_up)
target_data_induce_down$gene <- row.names(target_data_induce_down)
row.names(target_data_induce) <- 1:166 
row.names(target_data_induce_up) <- 1:107
row.names(target_data_induce_down) <- 1:70


###四个靶基因集
# ato直接的induce靶点，up /down ，等量随机
set.seed(10) #567都有不显著的
random_index_induce <- sample(1:7082,166,replace = F) 
random_data_induce <- rest_data[random_index_induce,]

all_induce_tuopu_data <- rbind(target_data_induce,random_data_induce)
all_induce_tuopu_data$gene <- row.names(all_induce_tuopu_data)

all_induce_tuopu_data <- rbind(all_induce_tuopu_data,target_data_induce_up)
all_induce_tuopu_data <- rbind(all_induce_tuopu_data,target_data_induce_down)
all_induce_tuopu_data$type <- factor(all_induce_tuopu_data$type,levels = c("target", "up", "down","others"))

all_tuopu_data_test_induce <- all_induce_tuopu_data

library(ggplot2)
library(ggpubr)

my_comparisons <- list( c("target","others"),c("up","others"),c("down","others"))

#度#### 1~24，改刻度0~30,三组都显著，可以p上去
summary(all_tuopu_data_test_induce[,1]) 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   2.000   4.000   6.045   8.000  81.000 
head(sort(all_tuopu_data_test_induce[,1],decreasing = T),30)

#开始
plot11<-ggplot(all_tuopu_data_test_induce[,c(1,5)], aes(x=type, y=degree,fill=type)) +
  geom_boxplot(outlier.shape = 21,outlier.alpha =0)+
  labs(x = "Type", y = "Degree") +
  scale_fill_manual(values=c("#E94849","#FF9933","#8FC9E8","#5CB356"))+#c("#F5E5BA","#9ECB55"))+
  theme_classic()+
  theme(axis.text.x = element_text(size=14),#,face="bold"),
        axis.text.y = element_text(size=14),#,face="bold"),
        axis.line = element_line(size=0.7),
        legend.text = element_text(size = 13),
        legend.title  = element_text(size = 13),
        axis.title = element_text(size=16)#,face="bold")
  )+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",cex=3)+#,label.y = c(18,20))+
  stat_compare_means(label.y = 29.5)+
  coord_cartesian(ylim = c(0,30))
plot11


#看一下是否显著(都是四星显著)
ggplot(all_tuopu_data_test_induce[,c(1,5)], aes(x=type, y=degree,fill=type)) +
  geom_boxplot(outlier.shape = 21,outlier.alpha =0)+
  labs(x = "Type", y = "Degree") +
  scale_fill_manual(values=c("#E94849","#FF9933","#8FC9E8","#5CB356"))+#c("#F5E5BA","#9ECB55"))+
  theme_classic()+
  theme(axis.text.x = element_text(size=14),#,face="bold"),
        axis.text.y = element_text(size=14),#,face="bold"),
        axis.line = element_line(size=0.7),
        legend.text = element_text(size = 13),
        legend.title  = element_text(size = 13),
        axis.title = element_text(size=16)#,face="bold")
  )+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",cex=3)


#聚类系数，0~1，大部分小
summary(all_tuopu_data_test_induce[,2]) 
#> summary(all_tuopu_data_test_induce[,2]) 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000000 0.000000 0.043230 0.003333 1.000000

plot2 <- ggplot(all_tuopu_data_test_induce[,c(2,5)], aes(x=type, y=log10(cluster_coefficient+0.001),fill=type)) +
  geom_boxplot(outlier.shape = 21,outlier.alpha =0)+
  labs(x = "Type", y = "Cluster") +
  scale_fill_manual(values=c("#E94849","#FF9933","#8FC9E8","#5CB356"))+
  theme_classic()+
  theme(axis.text.x = element_text(size=14),#,face="bold"),
        axis.text.y = element_text(size=14),#,face="bold"),
        axis.line = element_line(size=0.7),
        legend.text = element_text(size = 13),
        legend.title  = element_text(size = 13),
        axis.title = element_text(size=16)#,face="bold")
  )+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",cex=3)+#,label.y = c(0.045,0.05,0.055))+
  stat_compare_means(label.y = 1.3)+
  coord_cartesian(ylim = c(-3.1,1.5))
plot2  


#bet
all_tuopu_data_test_induce$betweenness<-log10(all_tuopu_data_test_induce$betweenness+0.01) #需要改刻度

plot3<-ggplot(all_tuopu_data_test_induce[,c(3,5)], aes(x=type, y=betweenness,fill=type)) +
  geom_boxplot(outlier.shape = 21,outlier.alpha =0)+
  labs(x = "Type", y = "Betweenness") +
  scale_fill_manual(values=c("#E94849","#FF9933","#8FC9E8","#5CB356"))+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",cex=3)+
  stat_compare_means(label.y =9.5)+
  theme_classic()+
  theme(axis.text.x = element_text(size=14),#,face="bold"),
        axis.text.y = element_text(size=14),#,face="bold"),
        axis.line = element_line(size=0.7),
        legend.text = element_text(size = 15),
        legend.title  = element_text(size = 15),
        axis.title = element_text(size=16)#,face="bold")
  )
# coord_cartesian(ylim = c(0,0.08)) 
plot3

#中心
all_tuopu_data_test_induce$closeness<-log10(all_tuopu_data_test_induce$closeness)
#plot4<-ggplot(all_tuopu_data_test[which(all_tuopu_data_test$closeness>(-6.123927)),c(4,5)], aes(x=type, y=closeness,fill=type)) +
plot4<-ggplot(all_tuopu_data_test_induce[,c(4,5)], aes(x=type, y=closeness,fill=type)) +
  geom_boxplot(outlier.shape = 21,outlier.alpha =0)+
  labs(x = "Type", y = "Closeness") +
  scale_fill_manual(values=c("#E94849","#FF9933","#8FC9E8","#5CB356"))+
  theme_classic()+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.line = element_line(size=0.7),
        legend.text = element_text(size = 15),
        legend.title  = element_text(size = 15),
        axis.title = element_text(size=16)
  )+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",cex=3)+#,label.y = c(7.5,7.510,7.520))+
  stat_compare_means(label.y = -6.105)+ #全局p值
  coord_cartesian(ylim = c(-6.120,-6.104))
plot4

#也都是四星显著
ggplot(all_tuopu_data_test_induce[,c(4,5)], aes(x=type, y=closeness,fill=type)) +
  geom_boxplot(outlier.shape = 21,outlier.alpha =0)+
  labs(x = "Type", y = "Closeness") +
  scale_fill_manual(values=c("#E94849","#FF9933","#8FC9E8","#5CB356"))+
  theme_classic()+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.line = element_line(size=0.7),
        legend.text = element_text(size = 15),
        legend.title  = element_text(size = 15),
        axis.title = element_text(size=16)
  )+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",cex=3)+#,label.y = c(7.5,7.510,7.520))+
  stat_compare_means(label.y = -6.105)

ggarrange(plot11, plot2, plot3, plot4,  
          ncol = 2, nrow = 2,
          common.legend = T,legend = "top")
#第一个第四个都是四星显著


#inhibit###############################
target_data_inhibit <- target_data[row.names(target_data) %in% true_ato_alone_mrna_target[[4]],] #741个inhibit
target_data_inhibit_up <- target_data_inhibit[row.names(target_data_inhibit) %in% true_ato_alone_mrna_target[[5]],] #511个up
target_data_inhibit_down <- target_data[row.names(target_data) %in% true_ato_alone_mrna_target[[6]],] #299个down

target_data_inhibit_up$type <- "up"
target_data_inhibit_down$type <- "down"

target_data_inhibit$gene <- row.names(target_data_inhibit)
target_data_inhibit_up$gene <- row.names(target_data_inhibit_up)
target_data_inhibit_down$gene <- row.names(target_data_inhibit_down)
row.names(target_data_inhibit) <- 1:741 
row.names(target_data_inhibit_up) <- 1:511
row.names(target_data_inhibit_down) <- 1:299


###四个靶基因集
# ato直接的inhibit靶点，up /down ，等量随机
set.seed(10) 
random_index_inhibit <- sample(1:7082,741,replace = F) 
random_data_inhibit <- rest_data[random_index_inhibit,]

all_inhibit_tuopu_data <- rbind(target_data_inhibit,random_data_inhibit)
all_inhibit_tuopu_data$gene <- row.names(all_inhibit_tuopu_data)

all_inhibit_tuopu_data <- rbind(all_inhibit_tuopu_data,target_data_inhibit_up)
all_inhibit_tuopu_data <- rbind(all_inhibit_tuopu_data,target_data_inhibit_down)
all_inhibit_tuopu_data$type <- factor(all_inhibit_tuopu_data$type,levels = c("target", "up", "down","others"))

all_tuopu_data_test_inhibit <- all_inhibit_tuopu_data


#度#### 1~24，改刻度0~30,三组都显著，可以p上去
summary(all_tuopu_data_test_inhibit[,1]) 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   2.000   4.000   6.493   8.000  81.000 
head(sort(all_tuopu_data_test_inhibit[,1],decreasing = T),30)

#开始
plot111<-ggplot(all_tuopu_data_test_inhibit[,c(1,5)], aes(x=type, y=degree,fill=type)) +
  geom_boxplot(outlier.shape = 21,outlier.alpha =0)+
  labs(x = "Type", y = "Degree") +
  scale_fill_manual(values=c("#E94849","#FF9933","#8FC9E8","#5CB356"))+#c("#F5E5BA","#9ECB55"))+
  theme_classic()+
  theme(axis.text.x = element_text(size=14),#,face="bold"),
        axis.text.y = element_text(size=14),#,face="bold"),
        axis.line = element_line(size=0.7),
        legend.text = element_text(size = 13),
        legend.title  = element_text(size = 13),
        axis.title = element_text(size=16)#,face="bold")
  )+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",cex=3)+#,label.y = c(18,20))+
  stat_compare_means(label.y = 29.5)+
  coord_cartesian(ylim = c(0,30))
plot111


#都是四星显著
ggplot(all_tuopu_data_test_inhibit[,c(1,5)], aes(x=type, y=degree,fill=type)) +
  geom_boxplot(outlier.shape = 21,outlier.alpha =0)+
  labs(x = "Type", y = "Degree") +
  scale_fill_manual(values=c("#E94849","#FF9933","#8FC9E8","#5CB356"))+#c("#F5E5BA","#9ECB55"))+
  theme_classic()+
  theme(axis.text.x = element_text(size=14),#,face="bold"),
        axis.text.y = element_text(size=14),#,face="bold"),
        axis.line = element_line(size=0.7),
        legend.text = element_text(size = 13),
        legend.title  = element_text(size = 13),
        axis.title = element_text(size=16)#,face="bold")
  )+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",cex=3)


#聚类系数，0~1，大部分小
summary(all_tuopu_data_test_inhibit[,2]) 
#> summary(all_tuopu_data_test_inhibit[,2]) 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.00000 0.00000 0.00000 0.03127 0.01515 1.00000 

#
plot22 <- ggplot(all_tuopu_data_test_inhibit[,c(2,5)], aes(x=type, y=log10(cluster_coefficient+0.001),fill=type)) +
  geom_boxplot(outlier.shape = 21,outlier.alpha =0)+
  labs(x = "Type", y = "Cluster") +
  scale_fill_manual(values=c("#E94849","#FF9933","#8FC9E8","#5CB356"))+
  theme_classic()+
  theme(axis.text.x = element_text(size=14),#,face="bold"),
        axis.text.y = element_text(size=14),#,face="bold"),
        axis.line = element_line(size=0.7),
        legend.text = element_text(size = 13),
        legend.title  = element_text(size = 13),
        axis.title = element_text(size=16)#,face="bold")
  )+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",cex=3)+#,label.y = c(0.045,0.05,0.055))+
  stat_compare_means(label.y = 1.3)+
  coord_cartesian(ylim = c(-3.1,1.5))
plot22  


#bet
all_tuopu_data_test_inhibit$betweenness<-log10(all_tuopu_data_test_inhibit$betweenness+0.01) #需要改刻度

plot33 <-ggplot(all_tuopu_data_test_inhibit[,c(3,5)], aes(x=type, y=betweenness,fill=type)) +
  geom_boxplot(outlier.shape = 21,outlier.alpha =0)+
  labs(x = "Type", y = "Betweenness") +
  scale_fill_manual(values=c("#E94849","#FF9933","#8FC9E8","#5CB356"))+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",cex=3)+
  stat_compare_means(label.y =9.5)+
  theme_classic()+
  theme(axis.text.x = element_text(size=14),#,face="bold"),
        axis.text.y = element_text(size=14),#,face="bold"),
        axis.line = element_line(size=0.7),
        legend.text = element_text(size = 15),
        legend.title  = element_text(size = 15),
        axis.title = element_text(size=16)#,face="bold")
  )
# coord_cartesian(ylim = c(0,0.08)) 
plot33

#中心
all_tuopu_data_test_inhibit$closeness<-log10(all_tuopu_data_test_inhibit$closeness)
#plot4<-ggplot(all_tuopu_data_test[which(all_tuopu_data_test$closeness>(-6.123927)),c(4,5)], aes(x=type, y=closeness,fill=type)) +
plot44 <- ggplot(all_tuopu_data_test_inhibit[,c(4,5)], aes(x=type, y=closeness,fill=type)) +
  geom_boxplot(outlier.shape = 21,outlier.alpha =0)+
  labs(x = "Type", y = "Closeness") +
  scale_fill_manual(values=c("#E94849","#FF9933","#8FC9E8","#5CB356"))+
  theme_classic()+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.line = element_line(size=0.7),
        legend.text = element_text(size = 15),
        legend.title  = element_text(size = 15),
        axis.title = element_text(size=16)
  )+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",cex=3)+#,label.y = c(7.5,7.510,7.520))+
  stat_compare_means(label.y = -6.105)+ #全局p值
  coord_cartesian(ylim = c(-6.120,-6.104))
plot44

#都是四星显著
ggplot(all_tuopu_data_test_inhibit[,c(4,5)], aes(x=type, y=closeness,fill=type)) +
  geom_boxplot(outlier.shape = 21,outlier.alpha =0)+
  labs(x = "Type", y = "Closeness") +
  scale_fill_manual(values=c("#E94849","#FF9933","#8FC9E8","#5CB356"))+
  theme_classic()+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.line = element_line(size=0.7),
        legend.text = element_text(size = 15),
        legend.title  = element_text(size = 15),
        axis.title = element_text(size=16)
  )+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",cex=3)+#,label.y = c(7.5,7.510,7.520))+
  stat_compare_means(label.y = -6.105)

ggarrange(plot111, plot22, plot33, plot44,  # 要排版的图形
          ncol = 2, nrow = 2,
          common.legend = T,legend = "top")#
#第一、四个都四星显著



#药物-靶点图###########################################
drug_target <- unique(data_929[,c(6,12,17)])
drug_target <- drug_target[-which(drug_target$Drug.Combination=="-"),]
drug_target <- drug_target[-which(drug_target$Target=="-"),]
drug_target[,1] <- gsub(" ","",drug_target[,1]) #有些药前面有空格，需要处理
drug_target <- unique(drug_target)

 # drug_target[,1] <- tolower(drug_target[,1])
 # for(i in 1:length(drug_target[,1])){
 #   x <- drug_target[i,1]
 #   substr(x,1,1) <- toupper(substring(x,1,1)) #替换为首字母大写
 #   drug_target[i,1] <- x
 # }
 drug_188 <- unique(drug_target[,1]) #188种组合药物（有作用靶点的，没筛选mrna靶点）#也没有筛选人类

#和原本的文件比对即可
drug_att <- read.csv("drug_node.csv",header = T)
sum(drug_att[,12] %in% drug_188) 
#根据104看是否在188里面，只有83个匹配到了(其他都没有靶点)
drug_att[which(!drug_att[,12] %in% drug_188),12]

#筛选网络
index <- which(drug_target[,1] %in% drug_att[,12])
drug_target_net <- unique(drug_target[index,c(1,3)]) #575条边,83个药物，280个靶点
write.table(drug_target_net,"drug_target_net.txt",sep="\t",row.names = F,col.names = T,quote=F)

new_att <- drug_att[which(drug_att[,12] %in% drug_188),]

target_uni <- unique(drug_target_net[,2])
target_att <- matrix(0,length(target_uni),11)
target_att <- as.data.frame(target_att)
target_att$name <- target_uni
names(target_att)[1:11] <- names(new_att)[1:11]
target_att$type <- 1
new_att$type <- 0  
new_att <- rbind(new_att,target_att)

write.table(new_att,"drug_target_att.txt",sep="\t",row.names = F,col.names = T,quote=F)

#必要基因的ppi网络提取
ppi_network2 <- read.table("ppis_duplicated.txt",sep="\t",header = T,as.is=T)
biyao_gene <- read.table("biyaotarget.txt",sep="\t",header=F,as.is=T)
index1 <- which(ppi_network2[,1] %in% biyao_gene[,1])
index2 <- which(ppi_network2[,2] %in% biyao_gene[,1])
index_all <- union(index1,index2) #2270条关系
#index_all <- intersect(index1,index2) #150条
biyao_net <- ppi_network2[index_all,c(1,2)]
biyao_nodes <- union(biyao_net[,1],biyao_net[,2]) #1753个结点

#结点属性
all_nodes_att <- cbind(biyao_nodes,rep(0,1753))  
index_biyao <- which(all_nodes_att[,1] %in% biyao_gene[,1])  #229个必要基因在ppi里
#union(ppi_network2[index1,1],ppi_network2[index2,2])
all_nodes_att[index_biyao,2] <- 1 
all_nodes_att <- as.data.frame(all_nodes_att)
names(all_nodes_att) <- c("gene","type")
write.table(all_nodes_att,"biyao_att.txt",sep="\t",row.names = F,col.names = T,quote=F)
write.table(biyao_net,"biyao_net.txt",sep="\t",row.names = F,col.names = T,quote=F)

  
  
  