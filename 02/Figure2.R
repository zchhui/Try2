rm(list = ls())
options(stringsAsFactors = F)
setwd("D:/R数据")
getwd()
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggalluvial)
data3<-read.table("data3.txt")
Ratio1 <- data3 %>%
  group_by(group, grade) %>%#分组
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))
Ratio1=na.omit(Ratio1)

Ratio1<-Ratio1 %>%
  group_by(group) %>%
  mutate(Percentage=paste0(round(n/sum(n)*100,3)))
library(ggplot2)

library(RColorBrewer)

library(ggplot2)
ggplot( Ratio1, aes( x = group, weight = as.numeric(Percentage), fill = grade))+
  geom_bar( position = "stack")+
  labs(x="", y = "Percentage") +
  scale_fill_manual(values=c("G1"="#FDC6D5","G2" = "#E1B1D5","G3" = "#A19BCB","G4" = "#6D749B"))
# fill为修改图例标题
