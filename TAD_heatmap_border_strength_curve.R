
options(scipen=10)
library(Matrix)
library(reshape2)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)

setwd("~/others/xenopus_hic/blood/dump_50kb_KR/")
dump.txt<-read.table("chr01_50kb_KR.txt", header=F, as.is=T)

setwd("~/others/xenopus_hic/blood/TAD_res")
border_strength<-read.table("chr01_border_strength_50kb.txt", header=T, as.is=T)
TAD_result<-read.table("chr01_TAD_50kb.txt",header=T, as.is=T)

setwd("~/others/xenopus_hic/stage11/TAD_res/")
border_strength11<-read.table("chr01_border_strength_50kb.txt", header=T, as.is=T)


#############plot parameter#############
your_chr<-"chr01"
dump.resolution<-50000
your.plot.start<-0
your.plot.range<-10000000

heatmap.title<-c("blood chr01")

chr.length<-max(dump.txt[,2])
##test:
#plot.area<-seq(your.plot.start,20000000,by=your.plot.range)
plot.area<-seq(your.plot.start,chr.length,by=your.plot.range)
plot.range.start<-plot.area[-length(plot.area)]

###########################plot in PDF#####################################

setwd("~/others/xenopus_hic/blood")
pdf(file="TAD_heatmap_2_border_strength_chr01.pdf",width=6, height=9)

for(i in plot.range.start){
  your_start<-i
  your_end<-your_start+your.plot.range
  

#############################heatmap####################################
colnames(dump.txt)<-c("start","end","counts")
dump.your.range<-dump.txt[dump.txt$start>your_start & dump.txt$start<your_end & dump.txt$end>your_start & dump.txt$end<your_end,]
colnames(dump.your.range)<-c("seq1","seq2","count")

seq1<-seq(your_start,your_end,by=dump.resolution)
seq2<-seq(your_start,your_end,by=dump.resolution)
hic.dataframe<-data.frame(seq1,seq2)
hic.dataframe$count<-0

inter.triple<-rbind(hic.dataframe, dump.your.range)

#weaken the max section(90%)
inter.triple[inter.triple$count>quantile(inter.triple$count, 0.95),3]<-quantile(inter.triple$count,0.95)

#generate the lower matrix data(exchange seq1 and seq2)
lower.seq1<-inter.triple$seq2
lower.seq2<-inter.triple$seq1
lower.count<-inter.triple$count
lower.ma<-data.frame(lower.seq1,lower.seq2,lower.count)
colnames(lower.ma)<-c("seq1","seq2","count")

inter.triple.all<-rbind(inter.triple, lower.ma)

#######
my_TAD_result=TAD_result[TAD_result[,2]>your_start & TAD_result[,1]<your_end,1:2]

###########plot heatmap
heatmap.gg<-
  ggplot(data = inter.triple.all)+
  geom_tile(aes(x=seq1,y=seq2,fill=count))+
  theme_classic()+
  scale_x_continuous(limits = c(your_start,your_end), 
                     breaks = seq(your_start, your_end, by=2000000),
                     labels = paste(seq(your_start, your_end, by=2000000)/1000000, "M", sep=""),
                     expand = c(0,0))+
  scale_y_continuous(limits=c(your_start,your_end),
                     breaks = seq(your_start, your_end, by=2000000),
                     labels = paste(seq(your_start, your_end, by=2000000)/1000000, "M", sep=""),
                     expand = c(0,0))+
  theme(axis.title.x = element_blank())+
  labs(y=your_chr, title=heatmap.title)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_gradient2(low="white",high="red")+
  guides(fill=F)+
  annotate("segment",x=my_TAD_result[,1],y=my_TAD_result[,1],
           xend=my_TAD_result[,1],yend=my_TAD_result[,2])+
  annotate("segment",x=my_TAD_result[,1],y=my_TAD_result[,2],
           xend=my_TAD_result[,2],yend=my_TAD_result[,2])


####################plot border_strength###################################
sample1<-rep("blood", nrow(border_strength))
bs1<-cbind(border_strength, sample1)
colnames(bs1)<-c("chr_num","border.position","border.strength","sample")

sample2<-rep("stage11", nrow(border_strength11))
bs2<-cbind(border_strength11, sample2)
colnames(bs2)<-c("chr_num","border.position","border.strength","sample")

bs_merge<-rbind(bs1, bs2)
my_border_strength<-bs_merge[bs_merge[,2]>your_start & bs_merge[,2]<your_end,]


#my_border_strength<-border_strength[border_strength[,2]>your_start & border_strength[,2]<your_end,]

bs.gg<-
  ggplot(data = my_border_strength, aes(x=my_border_strength[,2], y=my_border_strength[,3], color=sample))+
  #geom_line(size=0.6,color="darkblue")+
  geom_line(size=0.6)+
  theme_classic()+
  scale_x_continuous(limits = c(your_start,your_end), 
                     breaks = seq(your_start, your_end, by=2000000),
                     labels = paste(seq(your_start, your_end, by=2000000)/1000000, "M", sep=""),
                     expand = c(0,0))+
  xlab(your_chr)+
  ylab("border strength")+
  theme(legend.position = "bottom")

#plot.margin = margin(0.5,1,0.5,0.5,"cm"),
#legend.margin = unit(-2,"inches")

###########################
hg<-ggplotGrob(heatmap.gg)
bs<-ggplotGrob(bs.gg)
maxWidth<-grid::unit.pmax(hg$widths[2:3], bs$widths[2:3])
hg$widths[2:3] <- as.list(maxWidth)
bs$widths[2:3] <- as.list(maxWidth)

grid.arrange(hg,bs,ncol=1,nrow=2,heights=c(4,2),widths=4)

}

dev.off()
