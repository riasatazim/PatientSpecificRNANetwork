
library("igraph")

PAAD_mutationFile="C:/Users/dip/Documents/Sourcecodes/Pairwise/data/PAAD_mutation.csv"
PAAD_mutation= read.csv(file=PAAD_mutationFile)


liverNormalFile="C:/Users/dip/Documents/Sourcecodes/Pairwise/data/PAAD_Normal.gct"
PAAD_Normal= read.delim(file=liverNormalFile, skip=2)

liverCancerFile="C:/Users/dip/Documents/Sourcecodes/Pairwise/data/PAAD_Cancer.cct"
PAAD_Cancer= read.delim(file=liverCancerFile)

PAAD_Normal<-PAAD_Normal[!duplicated(PAAD_Normal$Description), ]
ins<-intersect(PAAD_Normal$Description,PAAD_Cancer$attrib_name)
ins<-unique(as.data.frame(ins))
PAAD_Normal<-merge(x=PAAD_Normal,y=ins,by.x ="Description" , by.y = "ins")


PAAD_Cancer<-PAAD_Cancer[!duplicated(PAAD_Cancer$attrib_name), ]
ins<-intersect(PAAD_Cancer$attrib_name,PAAD_Normal$Description)
ins<-unique(as.data.frame(ins))
PAAD_Cancer<-merge(x=PAAD_Cancer,y=ins,by.x ="attrib_name" , by.y = "ins")


rownames(PAAD_Normal)<-PAAD_Normal$Description
PAAD_Normal$id<-NULL
PAAD_Normal$Description<-NULL
PAAD_Normal$Name<-NULL
PAAD_Normal<-PAAD_Normal[order(row.names(PAAD_Normal)), ]

PAAD_Normal= apply(PAAD_Normal[1:nrow(PAAD_Normal),1:ncol(PAAD_Normal)],c(1,2),function(x)  log2(1+as.numeric (x)))

rownames(PAAD_Cancer)<-PAAD_Cancer$attrib_name
PAAD_Cancer$attrib_name<-NULL
PAAD_Cancer<-PAAD_Cancer[order(row.names(PAAD_Cancer)), ]

PAAD_Normal<-t(PAAD_Normal)
PAAD_Normal<-as.data.frame(PAAD_Normal)
PAAD_Cancer<-t(PAAD_Cancer)
PAAD_Cancer<-as.data.frame(PAAD_Cancer)

dataset<-rbind(PAAD_Cancer[5,],PAAD_Normal)
tdataset<-t(dataset)
tdataset<-as.data.frame(tdataset)

n<-ncol(tdataset)-1

#Clean
dataset<-NULL
ins<-NULL
PAAD_Cancer<-NULL
PAAD_Normal<-NULL
gc()


#Network Construction

tdataset<-t(tdataset)
cellNetwork<-cor(tdataset)

cellNetwork[is.na(cellNetwork)] <- 0

gc()
perturbed=cor(tdataset[-c(1), ])

gc()
perturbed[is.na(perturbed)] <- 0
perturbed=cellNetwork-perturbed
perturbed<-abs(perturbed)


cellNetwork=cellNetwork*cellNetwork
cellNetwork=1-cellNetwork



perturbed=perturbed/cellNetwork
perturbed=perturbed*(n)
tdataset<-NULL
cellNetwork<-NULL
#perturbed=abs(perturbed)


perturbed=abs(perturbed)
save(perturbed, file = "D:/Rdat/PAAD_5.RData")



perturbed["ABCA11P","PALMD"]
#View(perturbed[1:100,1:100])
perturbed[perturbed >0.0000001]<-0
perturbed[is.na(perturbed)] <- 0




network=graph_from_adjacency_matrix( perturbed, weighted=T, mode="undirected", diag=F)
edg=get.edgelist(network)
nrow(edg)
View(edg)


write.csv(edg,'C:/Users/dip/Desktop/LHC.csv')



PAAD_mutationFile="C:/Users/dip/Documents/Sourcecodes/Pairwise/data/Result/path/STAD/Result_PAAD_DesCom__global_3.csv"
PAAD_mutation= read.csv(file=PAAD_mutationFile)
X<-unique(PAAD_mutation$V1)
X<-X[1:1150]
X<-as.data.frame(X)
Y<-unique(PAAD_mutation$V2)
Y<-Y[1:1150]
Y<-as.data.frame(Y)
colnames(Y)<-"X"
X<-rbind(X,Y)

Z<-unique(X)
View(Z)

write.csv(Z,'C:/Users/dip/Documents/Sourcecodes/Pairwise/data/Result/path/temp.csv')





