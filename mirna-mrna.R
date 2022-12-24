mrna <- read.table("mRNA_target.txt",header = F,as.is=T) 
ppi <- read.table("ppis_duplicated.txt",header = T,as.is=T)
mimrna <- read.table("miRNA-mRNA.txt",header = T,as.is=T,skip=3)

#一个个来，apply，对每个mrna，找其所在的1 or 2，是1 就用2作为蛋白

allrowto <- NULL
getindex <- function(x){
  aindex <- which(ppi[,1] %in% x)
  arow <- ppi[aindex,]
  ppi_totala <- rowSums(arow[,3:9])
  allrowto <- rbind(allrowto,cbind(arow[,1:2],ppi_totala))
  colnames(allrowto)<-c("V1","V2","V3")
  
  bindex <- which(ppi[,2] %in% x)
  brow <- ppi[bindex,]
  ppi_totalb <- rowSums(brow[,3:9])
  brow2 <- as.data.frame(cbind(brow[,2],brow[,1],ppi_totalb))
  colnames(brow2)<-c("V1","V2","V3")
  allrowto <- rbind(allrowto,brow2)
  return(allrowto) 
}

art<-getindex(mrna[,1])
names(art) <- c("mRNA","Protein","PPI_total")
art<-art[order(art[,1]),]
write.table(art,"mRNAppi.txt",sep="\t",col.names = T,row.names = F,quote=F)


newindex<-which(mimrna[,4] %in% mrna[,1])
mim<-mimrna[newindex,c(4,2)]
names(mim)<-c("mRNA","miRNA")
mim <- unique(mim)
mim<-mim[order(mim[,1]),]
write.table(mim,"mRNA-miRNA.txt",sep="\t",col.names = T,row.names = F,quote=F)
