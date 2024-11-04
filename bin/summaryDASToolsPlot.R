#!/usr/bin/env Rscript

#library ------ 
library("ggpubr")
library("ggplot2")
library("ggprism")
library("rstatix")
library("gtools")
library("tidyverse")
library("dplyr") 
library("reshape2")


Args <- commandArgs(TRUE)


dasToolsResult <- read.table(Args[1],header=T, sep = "\t", check.name = F)


d3 <- dasToolsResult %>% dplyr::mutate(BinScore = case_when(
  bin_score > 0.9  ~ ">0.9",
  bin_score > 0.8  ~ ">0.8",
  bin_score > 0.7  ~ ">0.7",
  bin_score > 0.6  ~ ">0.6",
  bin_score >= 0.5  ~ ">=0.5",
  any(bin_score < 0.5) ~ "<0.5"
))


df <- d3 %>%
  group_by(binner, BinScore,sampleID) %>%
  summarise(counts = n())%>%
  group_by(binner, BinScore) %>% 
  summarise(avg_counts = mean(counts))


type_col <- rev(c("#f1eef6","#bdc9e1","#74a9cf","#2b8cbe","#045a8d","#062439"))
F1_list <- c(">0.9",">0.8",">0.7",">0.6",'>=0.5',"<0.5")

df$BinScore <- factor(df$BinScore,levels=F1_list)

df_sorted <- df %>%
  filter(BinScore %in% c( ">0.9",">0.8",">0.7",">0.6")) %>%
  group_by(binner) %>%
  summarise(counts = sum(avg_counts))%>%
  arrange(desc(counts))


df$binner <- factor(df$binner,levels=df_sorted$binner)

print(df_sorted$binner)

dd <- ggbarplot(df, x = "binner", y = "avg_counts", fill = "BinScore",size = 0.2,
            ggtheme = theme_bw() +
              theme(axis.text.x = element_text(color="black",size=10,angle=-90,hjust= 0 ,vjust = 0.5 ,face="bold"),
                    axis.text.y = element_text(color="black",size=10,face="bold"),
                    axis.title.y=element_text(color="black",size=10,face="bold"),
                    axis.title.x=element_text(color="black",size=10,face="bold"),
                    axis.line = element_line(color="black"),
                    axis.ticks = element_line(color="black"),
                    strip.text.x = element_text(size = 8, colour = "black", face="bold"), 
                    strip.background  = element_blank(),
                    legend.position = "top",
                    legend.text = element_text(size = 10, colour = "black"),
                    legend.title = element_text(size = 10),
                    legend.key.width = unit(0.3, 'cm'),
                    legend.key.size = unit(0, 'lines'),
                    panel.grid = element_blank(),
                    plot.margin = margin(l =5 , r = 20),
                    panel.background = element_blank()),
            legend = "top",title = "",xlab = 'Binner Combination', ylab = 'Avg (#Bins) Per Sample',width = 0.7)+ 
    scale_fill_manual(values=type_col,labels=F1_list,limits = F1_list, breaks = F1_list) 

drawwidth <-  length(levels(df$binner))*0.2 + 1
drawheight <- max(nchar(levels(df$binner)))*0.1 + 5

png(paste(Args[2],"_binner_DASToolsCombine_preSample_avgBinScore_counts.png",sep = ""),width = drawwidth,height = drawheight,units='in',res=600)
dd
dev.off()

write.table(df,file=paste(Args[2],"_binner_DASToolsCombine_preSample_avgBinScore_counts.xls",sep = ""),quote = F, sep = "\t", row.names= T, col.names = NA)



#bin percent
d3 <- dasToolsResult %>% dplyr::mutate(BinType = case_when(
  bin_score >= 0.5  ~ "bins",
  any(bin_score < 0.5) ~ 'unbins'
))

df <- d3 %>%
  group_by(binner,BinType,sampleID) %>%
  summarise(counts = n())%>%
  group_by(sampleID,binner)%>%
  summarize(
            counts = counts,
            BinType = BinType,
            Percentage = (counts /  sum(counts))*100)%>% 
  group_by(binner, BinType) %>% 
  summarise(avg_Prevalence = mean(Percentage))


bin_class <- c("#fb9a99","#80b1d3")
bin_class_order <- c('unbins',"bins")



df$BinType <- factor(df$BinType,levels=bin_class_order)

df_sorted <- df %>%
  filter(BinType == "bins") %>%
  group_by(binner) %>%
  summarise(counts = sum(avg_Prevalence))%>%
  arrange(desc(counts))

df$binner <- factor(df$binner,levels=df_sorted$binner)

dBin <- ggbarplot(df, x = "binner", y = "avg_Prevalence", fill = "BinType",size = 0.2,
            ggtheme = theme_bw() +
              theme(axis.text.x = element_text(color="black",size=10,angle=-90,hjust= 0 ,vjust = 0.5 ,face="bold"),
                    axis.text.y = element_text(color="black",size=10,face="bold"),
                    axis.title.y=element_text(color="black",size=10,face="bold"),
                    axis.title.x=element_text(color="black",size=10,face="bold"),
                    axis.line = element_line(color="black"),
                    axis.ticks = element_line(color="black"),
                    strip.text.x = element_text(size = 8, colour = "black", face="bold"), 
                    strip.background  = element_blank(),
                    legend.position = "top",
                    legend.text = element_text(size = 10, colour = "black"),
                    legend.title = element_text(size = 10),
                    legend.key.width = unit(0.3, 'cm'),
                    legend.key.size = unit(0, 'lines'),
                    panel.grid = element_blank(),
                    plot.margin = margin(l =5 , r = 20),
                    
                    panel.background = element_blank()),
            legend = "top",title = "",xlab = 'Binner Combination', ylab = 'Bining Percentage %',width = 0.7)+ 
    scale_fill_manual(values=bin_class,labels=bin_class_order,limits = bin_class_order, breaks = bin_class_order)



png(paste(Args[2],"_binner_DASToolsCombine_Bining_Percentage.png",sep = ""),width = drawwidth,height = drawheight,units='in',res=600)
dBin 
dev.off()

write.table(df,file=paste(Args[2],"_binner_DASToolsCombine_Bining_Percentage.xls",sep = ""),quote = F, sep = "\t", row.names= T, col.names = NA)
