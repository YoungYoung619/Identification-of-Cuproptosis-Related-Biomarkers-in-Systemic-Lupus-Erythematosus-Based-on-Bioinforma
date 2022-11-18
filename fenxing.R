

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")


library(ConsensusClusterPlus)        #���ð�
cellFile="37tong.txt"     #����ϸ�������ļ�
workDir="D:\\Maney\\SLE\\GSE126307\\fenxing"     #����Ŀ¼
setwd(workDir)       #���ù���Ŀ¼



data=read.table(cellFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)



maxK=9
results=ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="km",
              distance="euclidean",
              seed=123456,
              plot="png")




clusterNum=3        #�ּ��࣬�����жϱ�׼�ж�
cluster=results[[clusterNum]][["consensusClass"]]
write.table(cluster,file="ICIcluster.txt",sep="\t",quote=F,col.names=F)

