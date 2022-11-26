library(htmltools)
library( "DESeq2" )
library(ggplot2)




liverNormalFile="C:/Users/dip/Documents/Sourcecodes/Pairwise/data/STAD_Normal.gct"
LIHC_Normal= read.delim(file=liverNormalFile, skip=2)

liverCancerFile="C:/Users/dip/Documents/Sourcecodes/Pairwise/data/STAD_Cancer.cct"
LIHC_Cancer= read.delim(file=liverCancerFile)

LIHC_Normal<-LIHC_Normal[!duplicated(LIHC_Normal$Description), ]
ins<-intersect(LIHC_Normal$Description,LIHC_Cancer$attrib_name)
ins<-unique(as.data.frame(ins))
LIHC_Normal<-merge(x=LIHC_Normal,y=ins,by.x ="Description" , by.y = "ins")


LIHC_Cancer<-LIHC_Cancer[!duplicated(LIHC_Cancer$attrib_name), ]
ins<-intersect(LIHC_Cancer$attrib_name,LIHC_Normal$Description)
ins<-unique(as.data.frame(ins))
LIHC_Cancer<-merge(x=LIHC_Cancer,y=ins,by.x ="attrib_name" , by.y = "ins")


rownames(LIHC_Normal)<-LIHC_Normal$Description
LIHC_Normal$id<-NULL
LIHC_Normal$Description<-NULL
LIHC_Normal$Name<-NULL
LIHC_Normal<-LIHC_Normal[order(row.names(LIHC_Normal)), ]

LIHC_Normal= apply(LIHC_Normal[1:nrow(LIHC_Normal),1:ncol(LIHC_Normal)],c(1,2),function(x)  log2(1+as.numeric (x)))

rownames(LIHC_Cancer)<-LIHC_Cancer$attrib_name
LIHC_Cancer$attrib_name<-NULL
LIHC_Cancer<-LIHC_Cancer[order(row.names(LIHC_Cancer)), ]



LIHC_Normal<-t(LIHC_Normal)
new <- c()
for (val in 1: nrow(LIHC_Normal))
{
  # statement
  new[val]=0
}

LIHC_Normal <- cbind(new, LIHC_Normal)
LIHC_Normal<-as.data.frame(LIHC_Normal)
colnames(LIHC_Normal)[1] <- "labels"




LIHC_Cancer<-t(LIHC_Cancer)
new <- c()
for (val in 1: nrow(LIHC_Cancer))
{
  # statement
  new[val]=1
}


LIHC_Cancer <- cbind(new, LIHC_Cancer)
LIHC_Cancer<-as.data.frame(LIHC_Cancer)
colnames(LIHC_Cancer)[1] <- "labels"



dataset<-rbind(LIHC_Normal,LIHC_Cancer)

#write.csv(dataset,  file="C:/Users/dip/Documents/Sourcecodes/Pairwise/data/dataset_LIHC.csv")


View(dataset[1:10,1:10])

newdata <- dataset["labels"]
meta<-rownames(dataset)
meta = data.frame(newdata, meta)

View(meta)


dataset$labels<-NULL
dataset<-t(dataset)


dds <-  DESeqDataSetFromMatrix(countData=round(dataset),   colData=meta,  design=~labels, tidy = FALSE)
dds<-  DESeq(dds)





res <-  results(dds)
head(results(dds, tidy=TRUE))

summary(res)


res <-  res[order(res$padj),]
View(head(res))



up <- subset(res, res$log2FoldChange > 0,)
head(up)

down <- subset(res, res$log2FoldChange < 0,)
head(down)

write.csv(down,  file="C:/Users/dip/Documents/Sourcecodes/Pairwise/data/LIHC_5_down.csv")

write.csv(up,  file="C:/Users/dip/Documents/Sourcecodes/Pairwise/data/LIHC_5_up.csv")

View(up)
View(down)

plotMA(res, ylim=c(-2,2))




dds$labels <- factor(dds$labels, levels = c("0","1"))



-------------------STAD-----------------------------
plotCounts
par(mfrow=c(4,4))

plotCounts(dds, gene="NXF3",cex.axis=.8, intgroup="labels")
plotCounts(dds, gene="GPR18", intgroup="labels")
plotCounts(dds, gene="ILF2", intgroup="labels")
plotCounts(dds, gene="CD84", intgroup="labels")

plotCounts(dds, gene="EMP1", intgroup="labels")
plotCounts(dds, gene="ILF2", intgroup="labels")
plotCounts(dds, gene="SLC30A2", intgroup="labels")
plotCounts(dds, gene="POU3F1", intgroup="labels")

plotCounts(dds, gene="DDTL", intgroup="labels")
plotCounts(dds, gene="CLLU", intgroup="labels")
plotCounts(dds, gene="KRTAP4-12", intgroup="labels")
plotCounts(dds, gene="VRK1", intgroup="labels")


plotCounts(dds, gene="ARL5A", intgroup="labels")
plotCounts(dds, gene="DDTL", intgroup="labels")
plotCounts(dds, gene="GPS2", intgroup="labels")
plotCounts(dds, gene="SUPT3H", intgroup="labels")


plotCounts
par(mfrow=c(4,4))

plotCounts(dds, gene="GATA6",cex.axis=.8, intgroup="labels")
plotCounts(dds, gene="ZNF492", intgroup="labels")
plotCounts(dds, gene="KCNK18", intgroup="labels")
plotCounts(dds, gene="PRDM15", intgroup="labels")

plotCounts(dds, gene="FAM47B", intgroup="labels")
plotCounts(dds, gene="PRSS16", intgroup="labels")
plotCounts(dds, gene="BST1", intgroup="labels")
plotCounts(dds, gene="PCDHB1", intgroup="labels")


plotCounts(dds, gene="NXF3", intgroup="labels")
plotCounts(dds, gene="PLCXD2", intgroup="labels")
plotCounts(dds, gene="CBX5", intgroup="labels")
plotCounts(dds, gene="WISP1", intgroup="labels")


plotCounts(dds, gene="CHD7", intgroup="labels")
plotCounts(dds, gene="ZNF215", intgroup="labels")





---------------------------------PAAD
plotCounts
par(mfrow=c(4,4))

plotCounts(dds, gene="MT2A",cex.axis=.8, intgroup="labels")
plotCounts(dds, gene="TMEM115", intgroup="labels")
plotCounts(dds, gene="INS", intgroup="labels")
plotCounts(dds, gene="PAX6", intgroup="labels")

plotCounts(dds, gene="CLTA", intgroup="labels")
plotCounts(dds, gene="VCAM1", intgroup="labels")
plotCounts(dds, gene="OLIG1", intgroup="labels")
plotCounts(dds, gene="POU3F1", intgroup="labels")

plotCounts(dds, gene="CD46", intgroup="labels")
plotCounts(dds, gene="ZNF598", intgroup="labels")
plotCounts(dds, gene="KRTAP4-12", intgroup="labels")
plotCounts(dds, gene="VRK1", intgroup="labels")

plotCounts(dds, gene="INTS5", intgroup="labels")
plotCounts(dds, gene="PTPN18", intgroup="labels")
plotCounts(dds, gene="PRKX", intgroup="labels")
plotCounts(dds, gene="RAMP3", intgroup="labels")

plotCounts
par(mfrow=c(3,3))

plotCounts(dds, gene="REPS2", intgroup="labels")
plotCounts(dds, gene="VRK1", intgroup="labels")
plotCounts(dds, gene="ABCB11", intgroup="labels")

plotCounts(dds, gene="CCDC136", intgroup="labels")
plotCounts(dds, gene="CNGA4", intgroup="labels")
plotCounts(dds, gene="MT2A", intgroup="labels")

plotCounts(dds, gene="PAX6", intgroup="labels")
plotCounts(dds, gene="POU3F2", intgroup="labels")
plotCounts(dds, gene="TTC7A", intgroup="labels")



-------------------LIHC-------------------------
  
plotCounts
par(mfrow=c(4,4))

plotCounts(dds, gene="AFF3",cex.axis=.8, intgroup="labels")
plotCounts(dds, gene="BEST2", intgroup="labels")
plotCounts(dds, gene="ENO1", intgroup="labels")
plotCounts(dds, gene="TREM1", intgroup="labels")

plotCounts(dds, gene="UBAF", intgroup="labels")
plotCounts(dds, gene="TBX21", intgroup="labels")
plotCounts(dds, gene="PKD2L1", intgroup="labels")
plotCounts(dds, gene="ATXN3L", intgroup="labels")

plotCounts(dds, gene="DPP6", intgroup="labels")
plotCounts(dds, gene="HBA1", intgroup="labels")
plotCounts(dds, gene="KCNK4", intgroup="labels")
plotCounts(dds, gene="LBH", intgroup="labels")


plotCounts(dds, gene="KLK1", intgroup="labels")
plotCounts(dds, gene="SEMA6A", intgroup="labels")
plotCounts(dds, gene="TMEM9B", intgroup="labels")
plotCounts(dds, gene="ACOT1", intgroup="labels")


plotCounts
par(mfrow=c(3,2))

plotCounts(dds, gene="PPAN", intgroup="labels")
plotCounts(dds, gene="ENDOG", intgroup="labels")
plotCounts(dds, gene="ZNF174", intgroup="labels")

plotCounts(dds, gene="SOX9", intgroup="labels")
plotCounts(dds, gene="TBX21", intgroup="labels")


