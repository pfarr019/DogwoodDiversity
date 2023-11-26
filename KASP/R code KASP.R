# R script for producing KASP result graphic
setwd("")

library(RColorBrewer)
library(ggplot2)

### KASP graphic
 
KASP.raw=(read.table("KASP data for R.txt", header=TRUE))
KASP.matrix <- as.matrix(KASP.raw)
#sort by overall
#KASP.matrix.reordered <-KASP.matrix[order(KASP.matrix[,8],decreasing=FALSE),]
#KASP.matrix.reordered
#Remove "X" from the beginning of column names
for ( col in 1:ncol(KASP.matrix)){colnames(KASP.matrix)[col] <- sub("X","", colnames(KASP.matrix)[col])}
KASP.matrix
library(reshape2)
KASP.matrix.2<-melt(KASP.matrix)
KASPplot<-ggplot(KASP.matrix.2, aes(as.factor(Var2), Var1, fill= value)) + 
  geom_tile() +
  ggtitle("KASP data")+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")+
  #scale_fill_gradient(low="red",high="green",middle="yellow", na.value="black")
  scale_fill_gradient2(low = "#1B7837", mid = "lightsteelblue1", high = "#762A83", midpoint = 0.5, na.value = "black")+
  ylab("Cultivar")+
  xlab("Polymorphism")
KASPplot

#now do ddRadSeq data
ddRadSeq.raw=(read.table("ddRadSeq data for R.txt", header=TRUE))
ddRadSeq.matrix <- as.matrix(ddRadSeq.raw)
#sort by overall
#ddRadSeq.matrix.reordered <-ddRadSeq.matrix[order(ddRadSeq.matrix[,8],decreasing=FALSE),]
#ddRadSeq.matrix.reordered
#Remove "X" from the beginning of column names
for ( col in 1:ncol(ddRadSeq.matrix)){colnames(ddRadSeq.matrix)[col] <- sub("X","", colnames(ddRadSeq.matrix)[col])}
ddRadSeq.matrix
library(reshape2)
ddRadSeq.matrix.2<-melt(ddRadSeq.matrix)
ddRadSeqplot<-ggplot(ddRadSeq.matrix.2, aes(as.factor(Var2), Var1, fill= value)) + 
  geom_tile() +
  ggtitle("ddRadSeq data")+
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position = "none")+
  #scale_fill_gradient(low="red",high="green",middle="yellow", na.value="black")
  scale_fill_gradient2(low = "#1B7837", mid = "lightsteelblue1", high = "#762A83", midpoint = 0.5, na.value = "black")+
  ylab("Cultivar")+
  xlab("Polymorphism")
ddRadSeqplot

#now do overall data
color1=rgb(.2, .6, .4, 1)
overall.raw=(read.table("Overall data for R.txt", header=TRUE))
overall.matrix <- as.matrix(overall.raw)
#Remove "X" from the beginning of column names
library(reshape2)
overall.matrix.2<-melt(overall.matrix)
overallplot<-ggplot(overall.matrix.2, aes(Var2, Var1, fill= value)) + 
  geom_tile() +
  ggtitle("Overall data")+
  labs(fill = "Proportion
ssp. kousa")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=45, hjust=1),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank())+
  scale_x_discrete(labels = c("KASP", "ddRadSeq", "Structure\n1250\nMarkers"))+
  scale_fill_gradient2(low = "#1B7837", mid = "lightsteelblue1", high = "#762A83", midpoint = 0.5, na.value = "black")+
  #scale_fill_gradientn(colors = c(muted(rgb(.2, .6, .4, 1)),"lightsteelblue1","lightsteelblue1","lightsteelblue1", muted(rgb(.6, 0, .8, 1))), values= scales::rescale(c(0, 0.25, 0.5, 0.75, 1)))+
  ylab("Cultivar")
overallplot

#put the two plots together
library(ggpubr)
arrangedplots <-ggarrange(KASPplot, ddRadSeqplot, overallplot,
                          labels = c("a", "b", "c"),
                          nrow = 1,
                          align= "h",
                          widths = c(1.5, 0.95, 0.70))
arrangedplots
ggexport(arrangedplots, filename="arranged heatmap for KASP.tiff", width=6800, height=2700, res=600)

########## create scatterplot of KASP results for powerpoint presentation (recreating the plot from stepone in ggplot with better colors and legend)##########################
library(ggplot2)
Kasp645 <- read.csv("KASP 645 data for R.csv")
###for some reason, there are extra empty lines after the data if the csv file is resaved, must be deleted before plotting
ggplot(data=Kasp645, aes(x=Allele1, y=Allele2, color=Call, fill=Call, shape=Status))+
  geom_point(size=4)+
  guides(color = guide_legend(reverse=TRUE), shape=guide_legend(reverse=TRUE), fill="none")+
  scale_fill_manual(values=c(rgb(0, 0, 0, 0.5), "#1B7837B3", "lightsteelblue1", "#762A83B3"))+
  scale_color_manual(values=c("black", "#1B7837",  '#789ac7', "#762A83"))+
  scale_shape_manual(values=c(8,24))+
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5)) + #centers title
  ggtitle("KASP Allelic Discrimination Plot\nMarker 645") +
  xlab("ssp. kousa allele (VIC)") +
  ylab("ssp. chinensis allele (FAM)")
ggsave("KASP allelic discrimination plot 645.tiff", width=15, height=12, units="cm", dpi=300)
