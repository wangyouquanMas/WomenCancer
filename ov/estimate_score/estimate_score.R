#计算免疫评分


#导入avereps依赖包
library(limma)
library(estimate)
setwd("/Users/admin/Desktop/ov/3 estimate_score")          
inputFile="symbol.txt"                                                  #?????ļ?????


rt=read.table(inputFile,sep="\t",header=T,check.names=F)

#查看对象类型
class(rt)
#data.frame转为matrix
rt_matrix = as.matrix(rt)
#将rt_matrix 第一列设置为matrix行名称
rownames(rt_matrix)=rt_matrix[,1]
exp=rt_matrix[,2:ncol(rt_matrix)]
dimnames=list(rownames(exp),colnames(exp))
#因为 data 中的类型是 character需要转换为 numeric形式才可以进行下面rowMeans运算
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
#将相同行数据取均值
res = avereps(data,ID=rownames(data))
#要求行均值大于0
out=res[rowMeans(res)>0,]
out=rbind(ID=colnames(out),out)

write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=T,row.names = T)

filterCommonGenes(input.f="uniq.symbol.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct", 
              platform="illumina")

scores=read.table("estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(ID=colnames(scores),scores)
write.table(out,file="scores.txt",sep="\t",quote=F,col.names=F)

