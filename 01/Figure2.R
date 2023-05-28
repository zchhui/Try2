rm(list = ls())
library(circlize)
library(openxlsx)
setwd("D:/R图")
getwd()

circos.clear()
#####读取UCSC数据库的染色体区段信息
human_cytoband1 = read.cytoband(species = "hg38")$df
head(human_cytoband1) 

human_cytoband2 = read.cytoband(species = "hg38")$df
head(human_cytoband2)
######将第一列染色体分为两组 ：human1和human2，方便后面分两组
human_cytoband1[ ,1] = paste0("human1_", human_cytoband1[, 1])
head(human_cytoband1)

human_cytoband2[ ,1] = paste0("human2_", human_cytoband2[, 1])
head(human_cytoband2)
#####合并数据
cytoband = rbind(human_cytoband1, human_cytoband2)
head(cytoband)
dim(cytoband)
table(cytoband$V1) 
#####添加24条染色体
chromosome.index = c(paste0("human1_chr", c(1:22, "X", "Y")), 
                     rev(paste0("human2_chr", c(1:22, "X", "Y"))))
#初始化看一下是不是两组
circos.initializeWithIdeogram(cytoband, chromosome.index = chromosome.index)

circos.clear()
#####设置自圈图布局
circos.par(gap.after = c(rep(1, 23), 5, rep(1, 23), 5))
circos.initializeWithIdeogram(cytoband, plotType = NULL, 
                              chromosome.index = chromosome.index)
####创建绘图区域
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), 
              gsub(".*chr", "", CELL_META$sector.index), cex = 0.7, niceFacing = TRUE)
}, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)
#######用颜色凸显两组数据
highlight.chromosome(paste0("human1_chr", c(1:22, "X", "Y")), 
                     col = "#FF6666", track.index = 1)
highlight.chromosome(paste0("human2_chr", c(1:22, "X", "Y")), 
                     col = "#0099CC", track.index = 1)
######添加hg38的染色体信息
circos.genomicIdeogram(cytoband)
circos.clear()

######确定染色体范围，read.chrominfo直接读取ucsc数据库的染色体信息
human_chromInfo1 = read.chromInfo(species = "hg38")$df
head(human_chromInfo1)
human_chromInfo2 = read.chromInfo(species = "hg18")$df
head(human_chromInfo2)

human_chromInfo1[ ,1] = paste0("human1_", human_chromInfo1[, 1])
head(human_chromInfo1)
human_chromInfo2[ ,1] = paste0("human2_", human_chromInfo2[, 1])
head(human_chromInfo2)

chromInfo = rbind(human_chromInfo1, human_chromInfo2)
head(chromInfo);dim(chromInfo)
length(chromosome.index)
#######染色体号按因子排序
chromInfo[, 1] = factor(chromInfo[ ,1], levels = chromosome.index)
head(chromInfo)
table(chromInfo$chr)
####加上染色体框
circos.par(gap.after = c(rep(1, 23), 5, rep(1, 23), 5))
circos.genomicInitialize(chromInfo, plotType = NULL)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), 
              gsub(".*chr", "", CELL_META$sector.index), cex = 0.6, niceFacing = TRUE)
}, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)
highlight.chromosome(paste0("human1_chr", c(1:22, "X", "Y")), 
                     col = "#FF6666", track.index = 1)
highlight.chromosome(paste0("human2_chr", c(1:22, "X", "Y")), 
                     col = "#0099CC", track.index = 1)

circos.track(ylim = c(0, 1))
circos.clear()

circos.par(gap.after = c(rep(1, 23), 5, rep(1, 23), 5))
circos.genomicInitialize(chromInfo, plotType = NULL)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), 
              gsub(".*chr", "", CELL_META$sector.index), cex = 0.6, niceFacing = TRUE)
}, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)
highlight.chromosome(paste0("human1_chr", c(1:22, "X", "Y")), 
                     col = "#FF6666", track.index = 1)
highlight.chromosome(paste0("human2_chr", c(1:22, "X", "Y")), 
                     col = "#0099CC", track.index = 1)

circos.genomicIdeogram(cytoband)
######添加点信息
human_df1 = as.data.frame(read.xlsx("data1.xlsx",sheet = 6))
head(human_df1)
human_df2 = as.data.frame(read.xlsx("data1.xlsx",sheet = 7))
head(human_df2)


human_df1[ ,1] = paste0("human1_", human_df1[, 1])
human_df2[ ,1] = paste0("human2_", human_df2[, 1])
df = rbind(human_df1, human_df2)
head(df)


circos.genomicTrack(df, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, col = rand_color(1), cex = 0.5, ...)
})

circos.genomicLink(human_df1, human_df2,
                   lwd=1, 
                   lty="11", 
                   col="grey90") 




# 从新读取数据，这里没有分组了
human_df1 = as.data.frame(read.xlsx("data1.xlsx",sheet = 4))
head(human_df1)
dim(human_df1)
human_df2 = as.data.frame(read.xlsx("data1.xlsx",sheet = 5))
head(human_df2)
dim(human_df2)

circos.clear()
# 调用hg38染色体
pdf("circos3.pdf",width=10, height=10) # 设置PDF保存图像
circos.initializeWithIdeogram(species="hg38")

circos.genomicLink(human_df1, human_df2,  # 传入bed格式文件,chr,start,end,value
                   col="grey90", # 设置颜色
                                   )

dev.off()

