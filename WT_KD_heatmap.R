#!/usr/bin/Rscript

options(scipen=10)
library(Matrix)
library(reshape2)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)

######################load data###############################
setwd("~/wang_lab/hichip2")
WT<-read.table("./A549_DDX5_WT/dump_10kb/chr1_NONE_10kb.txt",header = F, as.is = T)
#KD<-read.table("./A549_DDX5_KD/dump_10kb/chr1_NONE_10kb.txt",header = F, as.is = T)
#Hi-C
KD<-read.table("~/wang_lab/HiC/A549_SRR5129660/dump_10kb/chr1_NONE_10kb.txt",header = F, as.is = T)


setwd("~/wang_lab/DDX5_ChIPseq/A549_merge/macs2_res")
peak.data1<-read.table("DDX5_A549_peaks.narrowPeak",header=F,as.is=T)
setwd("~/wang_lab/hichip2/chipseq")
peak.data6<-read.table("A549_H3K27ac.bed",header=F,as.is=T)


##########plot parameter####################
your_chr<-"chr1"
dump.resolution<-10000
your.plot.start<-0
your.plot.range<-4000000

#heatmap.title<-c("A549 WT-KD HiChIP")
heatmap.title<-c("A549 WT_HiChIP diff Hi-C")
peak1.ylab<-c("DDX5")
peak6.ylab<-c("H3K27ac")
WT.ylab<-c("WT HiChIP")
#KD.ylab<-c("KD HiChIP")
KD.ylab<-c("A549 Hi-C")

chr.length<-max(WT[,2])
plot.area<-seq(your.plot.start,chr.length,by=your.plot.range)
#test
#plot.area<-seq(your.plot.start,10000000,by=your.plot.range)
plot.range.start<-plot.area[-length(plot.area)]

###############################plot in PDF#####################################

setwd("~/wang_lab/hichip2")
#pdf(file="A549_WT-KD_HiChIP_chipseq.pdf",width=8, height=10)
pdf(file="A549_WT_HiChIP-HiC.pdf",width=8, height=10)

for(i in plot.range.start){
  your_start<-i
  your_end<-your_start+your.plot.range
  
  #your_start<-20000000
  #your_end<-24000000
  
  
  ###########WT matrix#############
  colnames(WT)<-c("start","end","counts")
  dump.your.range<-WT[WT$start>your_start & WT$start<your_end & WT$end>your_start & WT$end<your_end,]
  colnames(dump.your.range)<-c("seq1","seq2","count")
  
  seq1<-seq(your_start,your_end,by=dump.resolution)
  seq2<-seq(your_start,your_end,by=dump.resolution)
  hic.dataframe<-data.frame(seq1,seq2)
  hic.dataframe$count<-0
  
  inter.triple<-rbind(hic.dataframe, dump.your.range)
  
  #generate the lower matrix data(exchange seq1 and seq2)
  lower.seq1<-inter.triple$seq2
  lower.seq2<-inter.triple$seq1
  lower.count<-inter.triple$count
  lower.ma<-data.frame(lower.seq1,lower.seq2,lower.count)
  colnames(lower.ma)<-c("seq1","seq2","count")
  
  WT.triple.all<-rbind(inter.triple, lower.ma)
  
  WT.data<-xtabs(count~seq1+seq2,WT.triple.all)
  WT.ma<-as.matrix(WT.data)
  
  #################KD matrix####################
  colnames(KD)<-c("start","end","counts")
  dump.your.range<-KD[KD$start>your_start & KD$start<your_end & KD$end>your_start & KD$end<your_end,]
  colnames(dump.your.range)<-c("seq1","seq2","count")
  
  seq1<-seq(your_start,your_end,by=dump.resolution)
  seq2<-seq(your_start,your_end,by=dump.resolution)
  hic.dataframe<-data.frame(seq1,seq2)
  hic.dataframe$count<-0
  
  inter.triple<-rbind(hic.dataframe, dump.your.range)
  
  #generate the lower matrix data(exchange seq1 and seq2)
  lower.seq1<-inter.triple$seq2
  lower.seq2<-inter.triple$seq1
  lower.count<-inter.triple$count
  lower.ma<-data.frame(lower.seq1,lower.seq2,lower.count)
  colnames(lower.ma)<-c("seq1","seq2","count")
  
  KD.triple.all<-rbind(inter.triple, lower.ma)
  
  KD.data<-xtabs(count~seq1+seq2,KD.triple.all)
  KD.ma<-as.matrix(KD.data)
  
  
  #############normalization############
  proportion<-sum(WT.ma)/sum(KD.ma)
  WT.ma<-WT.ma/proportion
  
  
  ################WT-KD##################
  diff.ma<-WT.ma-KD.ma
  diff.triple<-as.data.frame(as.table(diff.ma))
  colnames(diff.triple)<-c("seq1", "seq2", "count")
  diff.triple$seq1<-as.numeric(as.character(diff.triple$seq1))
  diff.triple$seq2<-as.numeric(as.character(diff.triple$seq2))
  
  #############################heatmap####################################
  #weaken the max section(95%)
  diff.triple[diff.triple$count>quantile(diff.triple$count, 0.95),3]<-quantile(diff.triple$count,0.95)
  #weaken the min section(5%)
  diff.triple[diff.triple$count<quantile(diff.triple$count, 0.05),3]<-quantile(diff.triple$count,0.05)
  #diff.triple[diff.triple$count<quantile(diff.triple$count, 0.2),3]<-quantile(diff.triple$count,0.2)
  
  heatmap.gg<-
    ggplot(data = diff.triple)+
    geom_tile(aes(x=seq1,y=seq2,fill=count))+
    theme_classic()+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    theme(axis.title.x = element_blank())+
    labs(y=your_chr,title=heatmap.title)+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_fill_gradient2(low="blue",high="red")+
    guides(fill=F)
  
  #print(heatmap.gg)
  
  ###########################plot peak################################
  
  peak1<-peak.data1[peak.data1[,1] %in% your_chr & peak.data1[,2]>your_start & peak.data1[,2]<your_end, c(2,5)]
  colnames(peak1)<-c("position","score")
  peak1.gg<-
    ggplot(data = peak1,aes(x=peak1$position, y=peak1$score))+
    geom_segment(xend=peak1$position,yend=0,color = "darkred")+
    scale_y_continuous(expand = c(0,0))+
    theme_classic()+
    ylab(peak1.ylab)+
    scale_x_continuous(limits=c(your_start,your_end), expand = c(0,0))+
    theme(axis.line.x = element_line(color = "darkred"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_text(size = 10,vjust = 0.5, hjust = 0.5,angle = 0))
  
  
  peak6<-peak.data6[peak.data6[,1] %in% your_chr & peak.data6[,2]>your_start & peak.data6[,2]<your_end, c(2,5)]
  colnames(peak6)<-c("position","score")
  peak6.gg<-
    ggplot(data = peak6,aes(x=peak6$position, y=peak6$score))+
    geom_segment(xend=peak6$position,yend=0,color = "darkorchid4")+
    scale_y_continuous(expand = c(0,0))+
    theme_classic()+
    ylab(peak6.ylab)+
    scale_x_continuous(limits=c(your_start,your_end), expand = c(0,0))+
    theme(axis.line.x = element_line(color = "darkorchid4"),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_text(size = 10,vjust = 0.5, hjust = 0.5,angle = 0))
  
  #########################plot reads enrichment####################
  WT.reads<-rowSums(WT.ma)
  WT.reads.pos<-names(WT.reads)
  WT.contact<-data.frame(as.numeric(WT.reads.pos), WT.reads)
  colnames(WT.contact)<-c("position","count")
  WT.gg<-
    ggplot(data = WT.contact,aes(x=WT.contact$position, y=WT.contact$count))+
    geom_segment(xend=WT.contact$position,yend=0,color = "red")+
    scale_y_continuous(expand = c(0,0))+
    theme_classic()+
    ylab(WT.ylab)+
    scale_x_continuous(limits=c(your_start,your_end), expand = c(0,0))+
    theme(axis.line.x = element_line(color = "red"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          #axis.ticks.y = element_blank(),
          #axis.text.y = element_blank(),
          #axis.line.y = element_blank(),
          axis.title.y = element_text(size = 10,vjust = 0.5, hjust = 0.5,angle = 0))
  
  KD.reads<-rowSums(KD.ma)
  KD.reads.pos<-names(KD.reads)
  KD.contact<-data.frame(as.numeric(KD.reads.pos), KD.reads)
  colnames(KD.contact)<-c("position","count")
  KD.gg<-
    ggplot(data = KD.contact,aes(x=KD.contact$position, y=KD.contact$count))+
    geom_segment(xend=KD.contact$position,yend=0,color = "blue")+
    scale_y_continuous(expand = c(0,0))+
    theme_classic()+
    ylab(KD.ylab)+
    scale_x_continuous(limits=c(your_start,your_end), expand = c(0,0))+
    theme(axis.line.x = element_line(color = "blue"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          #axis.ticks.y = element_blank(),
          #axis.text.y = element_blank(),
          #axis.line.y = element_blank(),
          axis.title.y = element_text(size = 10,vjust = 0.5, hjust = 0.5,angle = 0))
  
  #######################plot arrange#################################
  hg<-ggplotGrob(heatmap.gg)
  pg1<-ggplotGrob(peak1.gg)
  pg6<-ggplotGrob(peak6.gg)
  WT.pg<-ggplotGrob(WT.gg)
  KD.pg<-ggplotGrob(KD.gg)
  
  maxWidth<-grid::unit.pmax(hg$widths[2:3], pg1$widths[2:3], pg6$widths[2:3], WT.pg$widths[2:3], KD.pg$widths[2:3])
  hg$widths[2:3] <- as.list(maxWidth)
  pg1$widths[2:3] <- as.list(maxWidth)
  pg6$widths[2:3] <- as.list(maxWidth)
  WT.pg$widths[2:3] <- as.list(maxWidth)
  KD.pg$widths[2:3] <- as.list(maxWidth)
  
  grid.arrange(hg, pg1, pg6, WT.pg, KD.pg, 
               ncol=1,nrow=5,
               heights=c(6,1,1,1,1),widths=6)
}

dev.off()


#############normalization############
#proportion<-sum(na.omit(WT[,3]))/sum(na.omit(KD[,3]))
#WT[,3]<-WT[,3]/mean(WT[,3])
#KD[,3]<-KD[,3]/mean(KD[,3])

