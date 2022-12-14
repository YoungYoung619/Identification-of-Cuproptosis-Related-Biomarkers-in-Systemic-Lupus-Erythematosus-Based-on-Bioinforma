
library(randomForest)
library(limma)
library(ggpubr)
set.seed(123)

inputFile="基因表达谱.txt"      
C="C"                        

#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)
data=data[,read.table("disease.txt", header=F, sep="\t", check.names=F)[,1]]
sample=read.table("sample.txt",sep="\t",header=F,check.names=F,row.names = 1)
data=data[rownames(sample),]
colnames(data)=gsub("-", "afaf", colnames(data))


afcon=sum(sample[,1]==C)
group=c(rep("con",afcon),rep("treat",nrow(data)-afcon))



rf=randomForest(as.factor(group)~., data=data, ntree=500)
pdf(file="forest.pdf", width=6, height=6)
plot(rf, main="Random forest", lwd=2)
dev.off()


optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)



importance=importance(x=rf2)
importance=as.data.frame(importance)
importance$size=gsub("-", "afaf", importance$size)
importance$size=rownames(importance)
importance=importance[,c(2,1)]
names(importance)=c("Gene","importance")


af=importance[order(importance$importance,decreasing = T),]
af=af[1:20,]
p=ggdotchart(af, x = "Gene", y = "importance",
             color = "importance", # Custom color palette
             sorting = "descending",                       # Sort value in descending order
             add = "segments",                             # Add segments from y = 0 to dots
             add.params = list(color = "lightgray", size = 2), # Change segment color and size
             dot.size = 6,                        # Add mpg values as dot labels
             font.label = list(color = "white", size = 9,
                               vjust = 0.5),               # Adjust label parameters
             ggtheme = theme_bw()         ,               # ggplot2 theme
             rotate=TRUE                                       )#翻转坐标轴 
p1=p+ geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  gradient_color(palette =c(ggsci::pal_npg()(2)[2],ggsci::pal_npg()(2)[1])      ) +#颜色
  grids()   
#保存图片
pdf(file="importance.pdf", width=6, height=6)
print(p1)
dev.off()


rfGenes=importance[order(importance[,"importance"], decreasing = TRUE),]
write.table(rfGenes, file="rfGenes.xls", sep="\t", quote=F, col.names=T, row.names=F)
