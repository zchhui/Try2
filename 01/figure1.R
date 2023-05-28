rm(list = ls())
options(stringsAsFactors = F)
setwd("D:/R数据")
getwd()
exp <- read.table("CO_prog_fdr05_adjust_5per(1).txt",header=T,as.is=T,sep = "\t")
C=exp[c(1:2)]
EXP<-read.table("co_edit_sample_matrix-1(1).txt",header=T,as.is=T,sep = "\t")
#write.table(EXP,"D:/R数据/EXPtumor.exp.txt",sep="\t",quote=F)#保存tumor数据
axp<-read.table("LIHC.htseq_fpkm_symbol(1).txt",header=T,as.is=T,sep = "\t")
colnames(axp)<- as.factor(substring(colnames(axp),1,15))#保留名字15位
rownames(axp)=axp[,1]
axp=axp[,-1]
####整和5个位点
colnames(C) <- c("site1","site2")
D<-merge(C,EXP)#提取5个位点
D=D[,-2]
D=D[,-1]
D=as.data.frame(D)
# 提取协同（>=1)和非协同（=0)的数据
xt<-D[,colSums(D^2) !=0]
##1.1 去除那些每个值都是零的行
##转换思维，每个数据都为0的行，则该列的和也为0
bxt<-D[,colSums(D) == 0]
bxt=bxt[1,]
bxt=-t(bxt)

xt=xt[1,]
xt[xt == "0"] = 1
xt<-t (xt)
bxt=as.data.frame(bxt)
xt=as.data.frame(xt)
###xt与fxt分别与表达数据结合
co1<-axp[,rownames(xt)]
nco1<-axp[,rownames(bxt)]
#寻找差异表达基因（5对）
data1<- cbind(co1,nco1);#合并数据
t_p.value <- function(x){wilcox.test(x[1:length(co1)],x[(length(co1)+1):ncol(data1)])$p.value};#t.test独立双样本t检验；ncol()返回exp.pro中的列数
fc <- function(x){mean(x[1:length(co1)]) / mean(x[(length(co1)+1):ncol(data1)]) };#求fc
p.value <- apply(data1,1,t_p.value);#apply对t_p.value这个指定运算函数求行平均
fc <-  apply(2^data1,1,fc);#2^data对fc函数求行平均，数字1为行，2为列
p.value <- p.value;
deg1<- cbind(p.value,p.adjust(p.value,"fdr"),fc);#行名 p.value ，fdt矫正后的p值 ，le.fc
deg1<-data.frame(deg1)
colnames(deg1)=c("p.value","fdr","fc")#给deg1改名
rownames(deg1) <- rownames(data1);
deg1.fdr <- deg1[which(deg1[,2]<0.05 & (abs(log2(deg1[,3]+10^(-10)))  >= log2(2) ) ),];#筛选p<0.05,logfc<1.5
deg1.fdr<-data.frame(deg1.fdr)#转换为数据
deg1.fdr$log2fc<-log2(abs(deg1.fdr$fc)) #取logfc值
deg1.up<- deg1.fdr[which(deg1.fdr[,2]<0.05 & (log2(deg1.fdr[,3]+10^(-10)))  >= log2(2) ) ,]
deg1.down<- deg1.fdr[which(deg1.fdr[,2]<0.05 & (log2(deg1.fdr[,3]+10^(-10)))  <= -log2(2) ) ,]
write.table(deg1.up,"D:/R数据/deg1.up.txt",sep="\t",quote=F)
write.table(deg1.down,"D:/R数据/deg1.down.txt",sep="\t",quote=F)

deg1.fdr<-rbind(deg1.up,deg1.down)
#热图.
deg1.fdr=as.data.frame(deg1.fdr)
x=deg1.fdr$log2fc#单独筛选logfc
names(x)=row.names(deg1.fdr)
x[1:4]
x[1:length(x)]
cg=c(names(head(sort(x),)),
     names(tail(sort(x),)))#取前100与后100个显著差异基因
cg=rownames(deg1.fdr)
data1=as.matrix(data1)
library(pheatmap)
pheatmap(data1[cg,],show_colnames = F,show_rownames = F)    #画图
n=t(scale(t(data1[cg,])))
n[n>4]=4
n[n< -4]=4
n[1:4,1:4]
pheatmap(n,show_colnames = F,show_rownames = F)#上下调程度

group_list1 <- factor(c(rep("synergism",83),rep("Non-synergistic",290)))#给协同非协同分组
group_list2 <- factor(c(rep("up",4),rep("down",32)))#给基因上调下调分组

ac=data.frame(g=group_list1)
cc=data.frame(g=group_list2)
rownames(cc)=rownames(n)
rownames(ac)=colnames(n) #将ac的行名也就分组信息（是‘Non_Tumor’还是‘Tumor’）给到n的列名，即热图中位于上方的分组信息
pheatmap(n,show_colnames =F,
         show_rownames = F,
         cluster_cols = F, 
         annotation_col=ac)




fz <- read.csv("clinical.csv",header = T,row.names = 1)
#fz<-fz[,c(2,5,6)]
#fz<-t(fz)
data3<- rbind(xt,bxt)
data3$group<-group_list1
data3=as.data.frame(data3)

data3$gender <- fz[as.factor(substring(rownames(data3),1,12)),5]
data3$grade <- fz[as.factor(substring(rownames(data3),1,12)),6]
data3$stage.T <- fz[as.factor(substring(rownames(data3),1,12)),7]
data3=data3[, -1]
colnames(data3)=c("group","gender","grade","stage.T")#?????ļ?????????

pheatmap(n,cluster_rows = F,cluster_cols = F,
         show_colnames = F,border_color = NA,show_rownames =F,
         annotation_col = data3, )

my <- read.table("CIBERSORT-Results.csv",header=T,as.is=T,sep = ",")
my<- t(my)
colnames(my)=my[1,]
my=my[-1,]
my<-data.frame(my)
colnames(my)=substring(colnames(my),1,15)


options(stringsAsFactors = F)
co <- my[,rownames(xt)]
nco <- my[,rownames(bxt)]
data2 <- cbind(co,nco);#合并数据
data2 = data2[c(1:22),]
data2 = data.matrix(data2)
t_p.value <- function(x){wilcox.test(x[1:length(co)],x[(length(co)+1):ncol(data2)])$p.value};#t.test独立双样本t检验；ncol()返回exp.pro中的列数
fc <- function(x){mean(x[1:length(co)]) / mean(x[(length(co)+1):ncol(data2)]) };#求fc
p.value <- apply(data2,1,t_p.value);#apply对t_p.value这个指定运算函数求行平均
fc <-  apply(2^data2,1,fc);#2^data对fc函数求行平均，数字1为行，2为列

deg1<- cbind(p.value,p.adjust(p.value,"fdr"),fc);#行名 p.value ，fdt矫正后的p值 ，le.fc
deg1<-data.frame(deg1)
colnames(deg1)=c("p.value","fdr","fc")#给deg1改名

rownames(deg1) <- rownames(data2);
deg1.fdr <- deg1[which(deg1[,1]<0.05 & (abs(log2(deg1[,3]+10^(-10)))  >= log2(1.2) ) ),];#筛选p<0.05,logfc<1.5
deg1.fdr<-data.frame(deg1.fdr)#转换为数据
deg1.fdr$log2fc<-log2(abs(deg1.fdr$fc)) #取logfc值



my <- read.table("CIBERSORT-Results.txt",header=T,as.is=T,sep = "\t")
data3$B.cells.naive <- my[as.factor(substring(rownames(data3),1,15)),2]

data3$T.cells.CD8 <- my[as.factor(substring(rownames(data3),1,15)),5]
data3$T.cells.CD4.memory.activated <- my[as.factor(substring(rownames(data3),1,15)),6]

colnames(data3)=c("group","gender","grade","stage.T","B.cells.naive","T.cells.CD8","T.cells.CD4")
colnames(cc)=c("Regulatory genes")
pheatmap(n,cluster_rows = F,cluster_cols = F,
         show_colnames = F,border_color = NA,show_rownames =F,
         annotation_col = data3,
         annotation_row = cc
         )
library(ComplexHeatmap)





write.table(data3,"D:/R数据/data3.txt",sep="\t",quote=F)#保存



