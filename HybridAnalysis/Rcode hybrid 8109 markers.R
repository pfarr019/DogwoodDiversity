setwd("C:/Users/pfarr/Documents/Dogwood research/Full Diversity study/Species evaluation 12-15-21/R")


# PCO analysis ------------------------------------------------------------

chybrid.raw8109=read.table("8109markersforR.txt", row = 1, skipNul = TRUE)
library(vegan)
chybrid.dist8109 <- vegdist(chybrid.raw8109, method = "gower", na.rm = TRUE)

#PCO
library(ape)
chybrid.pco8109ape <- pcoa(chybrid.dist8109)
#do not use one of the pre-defined corrections for negative eigenvalues, instead, treat them as zeros by only adding up the positive eigenvalues
pco1.8109ape <- chybrid.pco8109ape$values$Eigenvalues[1]/sum(chybrid.pco8109ape$values$Eigenvalues[1:sum(chybrid.pco8109ape$values$Eigenvalues>0)]) #only specify the positive values
pco2.8109ape <- chybrid.pco8109ape$values$Eigenvalues[2]/sum(chybrid.pco8109ape$values$Eigenvalues[1:sum(chybrid.pco8109ape$values$Eigenvalues>0)]) #only specify the positive values
pco1.8109ape
pco2.8109ape

plot(chybrid.pco8109ape$vectors[,1], chybrid.pco8109ape$vectors[,2], type = "p", xlab = "PCO1", ylab = "PCO2", axes = TRUE, main = "PCO All Samples")
#add small labels that make the PCO too busy in order to see trends
plot(chybrid.pco8109ape$vectors[,1], chybrid.pco8109ape$vectors[,2], type = "n", xlab = "PCO1", ylab = "PCO2", axes = TRUE, main = "PCO All Samples")
text(chybrid.pco8109ape$vectors[,1], chybrid.pco8109ape$vectors[,2], labels(chybrid.pco8109ape$vectors[,1]), cex = 0.1, xpd = TRUE)

library(ggplot2)
chybrid.pcogroups=scan("hybridgroupsPCO.txt")
chybrid.pcogroups2 <- as.factor(chybrid.pcogroups) #convert groups of bract color into factors so can be read by ggplot
chybrid.pco8109dataframe <- data.frame(chybrid.pco8109ape$vectors[,1], chybrid.pco8109ape$vectors[,2], chybrid.pcogroups) #make the eigenvectors and PCO groupsinto a data frame
ggplot(chybrid.pco8109dataframe, aes(x=chybrid.pco8109ape$vectors[,1], y=chybrid.pco8109ape$vectors[,2], color=chybrid.pcogroups2, fill=chybrid.pcogroups2)) +
  geom_point(size=3, shape=21) + #makes dots larger and gives them shape with outline (outline looks good when you print it out)
  scale_fill_manual(values=c(rgb(0, 0, 1, .1), rgb(0.4, 0, 0.6, .3),rgb(0.7, 0, 0.3, .2),rgb(1, 0, 0, 0.1),rgb(1, 0.6, 0, 0.4),  rgb(1, 1, 0, 0.3), rgb(0, 0.6, 0, .3)),
                    labels=c("C. nuttallii","C. x elwinortonii","C. x elwinortonii \nbackcrossed\nto C. kousa","C. kousa","C. x rutgersensis", "C. florida", "C. nuttallii x C. florida")) + #color scale for fill
  scale_color_manual(values=c(rgb(0, 0, 1, 1), rgb(0.4, 0, 0.6, 1),rgb(0.7, 0, 0.3, 1),rgb(1, 0, 0, 0.8),rgb(1, 0.6, 0, 1),   rgb(0, 0, 0, 0.5),  rgb(0, 0.6, 0, 1)),
                     labels=c("C. nuttallii", "C. x elwinortonii","C. x elwinortonii \nbackcrossed\nto C. kousa","C. kousa","C. x rutgersensis", "C. florida", "C. nuttallii x C. florida")) + #color scale for outline  theme_bw() + #removes gray background
  theme_bw()+
  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), #gets rid of the minor grid lines
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), #gets rid of the major grid lines
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5)) + #centers title
  guides(color = guide_legend(override.aes = list(size = 3))) + #changes the legend to make the dots larger
  coord_fixed(ratio=1)+ #makes sure the axes have the same scale
  ggtitle("PCO Cornus hybrid analysis") +
  xlab("PCO1 62.6%") +
  ylab("PCO2 32.2%")
ggsave("PCO Cornus hybrid.tiff", scale=1, width=7, units="in", dpi=600)

#larger plot for powerpoint
chybrid.pcogroups=scan("hybridgroupsPCO.txt")
chybrid.pcogroups2 <- as.factor(chybrid.pcogroups) #convert groups of bract color into factors so can be read by ggplot
chybrid.pco8109dataframe <- data.frame(chybrid.pco8109ape$vectors[,1], chybrid.pco8109ape$vectors[,2], chybrid.pcogroups) #make the eigenvectors and PCO groupsinto a data frame
ggplot(chybrid.pco8109dataframe, aes(x=chybrid.pco8109ape$vectors[,1], y=chybrid.pco8109ape$vectors[,2], color=chybrid.pcogroups2, fill=chybrid.pcogroups2)) +
  geom_point(size=5, shape=21) + #makes dots larger and gives them shape with outline (outline looks good when you print it out)
  scale_fill_manual(values=c(rgb(0, 0, 1, .1), rgb(0.4, 0, 0.6, .3),rgb(0.7, 0, 0.3, .2),rgb(1, 0, 0, 0.1),rgb(1, 0.6, 0, 0.4),  rgb(1, 1, 0, 0.3), rgb(0, 0.6, 0, .3)),
                    labels=c("C. nuttallii","C. x elwinortonii","C. x elwinortonii \nbackcrossed\nto C. kousa","C. kousa","C. x rutgersensis", "C. florida", "C. nuttallii x C. florida")) + #color scale for fill
  scale_color_manual(values=c(rgb(0, 0, 1, 1), rgb(0.4, 0, 0.6, 1),rgb(0.7, 0, 0.3, 1),rgb(1, 0, 0, 0.8),rgb(1, 0.6, 0, 1),   rgb(0, 0, 0, 0.5),  rgb(0, 0.6, 0, 1)),
                     labels=c("C. nuttallii", "C. x elwinortonii","C. x elwinortonii \nbackcrossed\nto C. kousa","C. kousa","C. x rutgersensis", "C. florida", "C. nuttallii x C. florida")) + #color scale for outline  theme_bw() + #removes gray background
  theme_bw()+
  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), #gets rid of the minor grid lines
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), #gets rid of the major grid lines
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5)) + #centers title
  guides(color = guide_legend(override.aes = list(size = 5))) + #changes the legend to make the dots larger
  coord_fixed(ratio=1)+ #makes sure the axes have the same scale
  ggtitle("PCO Cornus hybrid analysis") +
  xlab("PCO1 62.6%") +
  ylab("PCO2 32.2%")
ggsave("PCO Cornus hybrid larger.tiff", scale=1, width=12, units="in", dpi=600)

# visualizing and comparing the structure results, observed species composition, and expected species composition --------

library("pophelper")
library(gridExtra)

setwd("C:/Users/Pfarr/Documents/Dogwood research/Full Diversity study/Species evaluation 12-15-21/Structure")

sfiles <- list.files(path="./pophelper", full.names=TRUE) #read all of the files in a folder
sfiles

sfiles1 <- list.files(path="./pophelper3-24-22", full.names=TRUE) #read all of the files in a folder
sfiles1

# basic usage
slist <- readQ(files=sfiles)
slist

p1 <- plotQ(qlist=slist[c(4,5,11)],imgoutput="join",returnplot=TRUE,exportplot=FALSE,basesize=11,exportpath =getwd())
grid.arrange(p1$plot[[1]])


#comparing structure results after two different SNP data sets (filtered differently)
slist2 <- alignK(slist[c(1,2)])
slist2
p1 <- plotQ(qlist=slist2,imgoutput="join",returnplot=TRUE,exportplot=FALSE,basesize=11, showyaxis=TRUE, showticks=TRUE)
grid.arrange(p1$plot[[1]])
#result, the structure results are the same, even if I filter the SNPs slightly differently

#import names of individuals
nam=read.csv("species names for pophelper.csv",header=FALSE,sep=",")
nam
if(length(unique(sapply(slist1,nrow)))==1) slist1 <- lapply(slist1,"row.names<-",nam[,1])


inds <- read.delim("species names for pophelper.txt",header=FALSE,stringsAsFactors=F)
rownames(slist1[[2]]) <- inds$V1

rownames(slist[[15]]) <- inds$V1
#careful because the lab legend isn't necessarily in order
plotQ(slist1[2],returnplot=TRUE,#exportplot=FALSE,
      showtitle=T, titlelab="Species Analysis of Rutgers Dogwood Breeding Material", titlesize=5,
      #sortind="all",
      clustercol=c("blue", "red","yellow"), height=5, width=10, barsize=1, barbordersize = 0, barbordercolour = "black",
      showlegend=T, legendlab=c("nuttallii   ","kousa   ","florida"), legendkeysize=3, legendtextsize=3, 
      showsp=TRUE, splab="K=3", splabsize=3, 
      showindlab=TRUE, useindlab=TRUE, indlabsize=2,
      outputfilename="plotq",imgtype="png", dpi=600, exportpath=getwd()
)


#######final data, in the right order
#using the calculated species composition instead of the structure results because the structure results miss small, common introgressions
setwd("C:/Users/pfarr/Documents/Dogwood research/Full Diversity study/Species evaluation 12-15-21/Structure/Final Final data for paper graphic")
sfiles <- "calculated results.txt"
sfiles

slist <- readQ(files=sfiles)
slist

#import names of individuals
names=read.csv("species names for pophelper.csv",header=FALSE,sep=",")
names
if(length(unique(sapply(slist,nrow)))==1) slist <- lapply(slist,"row.names<-",names[,1])
slist

labelgroups <- read.delim("species groups for pophelper.txt",header=FALSE,stringsAsFactors=FALSE) #read in file with groups
labelgroups1 <- labelgroups[,2:3]
labelgroups2 <- labelgroups[,2, drop=FALSE]
colnames(labelgroups2) <- c('type')
plotQ(slist,returnplot=TRUE,#exportplot=FALSE,
      showtitle=T, titlelab="Species Analysis of Rutgers Dogwood Breeding Material", titlesize=5,
      clustercol=c("yellow", "red","blue"), height=5, width=10, barsize=1, barbordersize = 0, barbordercolour = "black",
      showlegend=T, legendlab=c("florida   ","kousa   ","nuttallii"), legendkeysize=3, legendtextsize=3,  
      showsp=TRUE, splab="K=3", splabsize=3, grplabheight=0,
      showindlab=TRUE, useindlab=TRUE, indlabsize=2,indlabangle=90,
      outputfilename="plotq1",imgtype="png", dpi=600, exportpath=getwd())

#plot observed and expected results right next to one another
sfiles1 <- "expected results.txt" #read in the expected results based on pedigree calculations
sfiles1

slist1 <- readQ(files=sfiles1)
slist1
if(length(unique(sapply(slist1,nrow)))==1) slist1 <- lapply(slist1,"row.names<-",names[,1])
slist1

slistcombined <-joinQ(slist1,slist)
slistcombined

###make a combined plot that can be rotated on its side for easier viewing of the names###
#will have to cut and paste the title to the top
plotQ(slistcombined, imgoutput= 'Join', returnplot=TRUE,#exportplot=FALSE,
      showtitle=T, titlelab="Species Analysis of Rutgers Dogwood Breeding Material ", titlesize=15,
      clustercol=c("yellow", "red","blue"), height=9, width=25, barsize=1, barbordersize = 0, barbordercolour = "black", titlehjust=1, 
      showlegend=T, legendlab=c("florida   ","kousa   ","nuttallii  "), legendkeysize=11, legendtextsize=11, legendpos="right",  
      showsp=TRUE, splab=c("Expected\nspecies composition", "Observed\nspecies composition"), splabsize=12, splabangle=90, sppos="left",
      showindlab=TRUE, useindlab=TRUE, indlabsize=9,indlabangle=90,
      outputfilename="plotqObsvsExp to turn2",imgtype="tiff", dpi=600, exportpath=getwd()
) #keep at dpi=600, otherwise not clear in picture


# scatterplots showing hybrid DNA genomic locations -----------------------

####try importing the data as is, then melting it to get it to "lay it flat"
dataraw <- read.csv("8109markersToMelt.csv", check.names=FALSE)
dataraw
#note: it seems to be ok that there are both NA and <NA> in the data frame, they both show up as TRUE when is.na(dataraw) is run

class(dataraw) #should return a data.frame
library(reshape2)
molten.data <- melt(dataraw, id= "row", na.rm = TRUE) #melt the data
molten.data

molten.data$value<-gsub("lcl[|]Contig","",as.character(molten.data$value))
molten.data
library(dplyr)
library(tidyr)
molten.data1 <-molten.data %>% separate(value, c("contig", "bp_contig", "het_or_homo"), sep= " ") #separate value column based on spaces
molten.data1

molten.data1$contig<-gsub("\\<0\\>","11",molten.data1$contig) #replace an exact match of 0 with 11

categories <- unique(molten.data1$contig) 
categories

numberOfCategories <- length(categories)
numberOfCategories

molten.data2 <-molten.data1 %>% separate(variable, c("distinguishing_locus", "individual"), sep= " ") #separate value column based on spaces
molten.data2

write.csv(molten.data2,"8109toCSV.csv", row.names = TRUE)
#### in the csv, searched and replaced "Contreras_Garden" to "2021_RC", D-038 to KF137-62, and Ck_K149-7 to Ck-K194-7 in order to update names
##Also changed three of the unlabeled accessions to shorter labeled names
#Rutgers_unlabeled_1 = Unlabeled_behind_pin_oak_Orton's_House
#Rutgers_unlabeled_2 = unlabeled_between_F70R06P26_and_KN161-119
#Rutgers_unlabeled_3 = unlabeled_hybrid_in_trio_with_R2P21_and_R1P17

alldata <- read.csv("8109fromCSVtoR.csv", check.names=FALSE)
alldata
alldata$individual<-gsub("\"","",alldata$individual)

###for some reason, there are extra empty lines after the data if the csv file is resaved, must be deleted before plotting

###try the same plot, but with layers in order to get the different points to not overlap
specified_order2 <-c("KN144-12", "KF83-1")
level_order <- factor(dataforgraphing2$individual, level = c("KN144-19", "KF83-1"))
level_order

library(magrittr)
subsetd <- alldata[alldata$individual %in% c("KN144-19", "D-038","KF83-1"), ] #create a subset with only these three individual's data

dataforgraphing1 <- alldata[alldata$order_in_graph<=21,] #data for first graph
dataforgraphing1$new_name = with(dataforgraphing1, reorder(new_name, order_in_graph)) #reorder based on the number in order_in_graph. Need to remember that 
#the lower number is at the bottom of the graph
dataforgraphing1

dataforgraphing2 <- alldata[alldata$order_in_graph>21 & alldata$order_in_graph<=42,]
dataforgraphing2$new_name = with(dataforgraphing2, reorder(new_name, order_in_graph)) #reorder based on the number in order_in_graph. Need to remember that
#the lower number is at the bottom of the graph

dataforgraphing3 <- alldata[alldata$order_in_graph>42,] #data for third graph
dataforgraphing3$new_name = with(dataforgraphing3, reorder(new_name, order_in_graph)) #reorder based on the number in order_in_graph. Need to remember that
#the lower number is at the bottom of the graph

#plotting the data in three separate plots so can see all individuals' data clearly
ggplot()+
  geom_point(data=subset(dataforgraphing3,distinguishing_locus %in% "K"), 
             aes(x=overall_bp, y=new_name, color="red"),
             size=1, alpha=0.2)+ #position=position_jitter(h=0.15))+
  scale_x_continuous(breaks=c(0,129942925,252856830,359034702,479829781,607726880,730331874,829888132,935043946,1040321751,1165293952,1318504356), 
                     labels=c("Contig 1", "Contig 2", "Contig 3", "Contig 4", "Contig 5", "Contig 6", "Contig 7","Contig 8",
                              "Contig 9", "Contig 10", "Contig 0 (11)", ""))+
  #scale_y_discrete(limits=specified_order2)+
  theme(panel.grid.minor.x = element_blank(), #axis.text.x=element_blank(), #gets rid of the minor grid lines
        panel.grid.major.x=element_line(color=rgb(0,0,0,0.4)), #turns the x axis gridlines gray
        axis.text.x= element_text(hjust=0, size=8),
        #legend.title=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5))+
  guides(color = guide_legend(override.aes = list(size =3, alpha= 0.7)))+ #changes the legend to make the dots larger and not transparent
  ggtitle("Species specific loci genome distribution") +
  labs(x=expression(paste("Contig in ",italic("C. florida"), " reference genome"))) +
  ylab("Breeding Selection")+
  geom_point(data=subset(dataforgraphing3,distinguishing_locus %in% "N"), 
             aes(x=overall_bp, y=new_name, color="blue"),
             size=1, alpha=0.2, position=position_nudge(y=0.17))+
  geom_point(data=subset(dataforgraphing3,distinguishing_locus %in% "F"), #gives the yellow dots a black outline to make them easier to see
             aes(x=overall_bp, y=new_name),
             size=1.1, alpha=0.9, color="black", position=position_nudge(y=-0.165))+
  geom_point(data=subset(dataforgraphing3,distinguishing_locus %in% "F"), 
             aes(x=overall_bp, y=new_name, color="yellow"),
             size=1, alpha=1, position=position_nudge(y=-0.17))+
  scale_color_identity(guide = "legend", #forces the aes to recognize the color names when they are inside aes
                       name= "Species Specific Loci",
                       labels=c("C. nuttallii   ", "C. kousa   ", "C. florida")) 
ggsave("Species specific loci1.tiff", width=22.23, height=19.05, units="cm", dpi=300)
ggsave("Fig14.tiff", width=22.23, height=19.05, units="cm", dpi=300)

ggplot()+
  geom_point(data=subset(dataforgraphing2,distinguishing_locus %in% "K"), 
             aes(x=overall_bp, y=new_name, color="red"),
             size=1, alpha=0.2)+ #position=position_jitter(h=0.15))+
  scale_x_continuous(breaks=c(0,129942925,252856830,359034702,479829781,607726880,730331874,829888132,935043946,1040321751,1165293952,1318504356), 
                     labels=c("Contig 1", "Contig 2", "Contig 3", "Contig 4", "Contig 5", "Contig 6", "Contig 7","Contig 8",
                              "Contig 9", "Contig 10", "Contig 0 (11)", ""))+
  #scale_y_discrete(limits=specified_order2)+
  theme(panel.grid.minor.x = element_blank(), #axis.text.x=element_blank(), #gets rid of the minor grid lines
        panel.grid.major.x=element_line(color=rgb(0,0,0,0.4)), #turns the x axis gridlines gray
        axis.text.x= element_text(hjust=0, size=8),
        #legend.title=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5))+
  guides(color = guide_legend(override.aes = list(size =3, alpha= 0.7)))+ #changes the legend to make the dots larger and not transparent
  ggtitle("Species specific loci genome distribution") +
  labs(x=expression(paste("Contig in ",italic("C. florida"), " reference genome"))) +
  ylab("Breeding Selection")+
  geom_point(data=subset(dataforgraphing2,distinguishing_locus %in% "N"), 
             aes(x=overall_bp, y=new_name, color="blue"),
             size=1, alpha=0.2, position=position_nudge(y=0.17))+
  geom_point(data=subset(dataforgraphing2,distinguishing_locus %in% "F"), #gives the yellow dots a black outline to make them easier to see
             aes(x=overall_bp, y=new_name),
             size=1.1, alpha=0.9, color="black", position=position_nudge(y=-0.165))+
  geom_point(data=subset(dataforgraphing2,distinguishing_locus %in% "F"), 
             aes(x=overall_bp, y=new_name, color="yellow"),
             size=1, alpha=1, position=position_nudge(y=-0.17))+
  scale_color_identity(guide = "legend", #forces the aes to recognize the color names when they are inside aes
                       name= "Species Specific Loci",
                       labels=c("C. nuttallii   ", "C. kousa   ", "C. florida")) 
ggsave("Species specific loci2.tiff", width=22.23, height=19.05, units="cm", dpi=300)
ggsave("Fig15.tiff", width=22.23, height=19.05, units="cm", dpi=300)


ggplot()+
  geom_point(data=subset(dataforgraphing1,distinguishing_locus %in% "K"), 
             aes(x=overall_bp, y=new_name, color="red"),
             size=1, alpha=0.2)+ #position=position_jitter(h=0.15))+
  scale_x_continuous(breaks=c(0,129942925,252856830,359034702,479829781,607726880,730331874,829888132,935043946,1040321751,1165293952,1318504356), 
                     labels=c("Contig 1", "Contig 2", "Contig 3", "Contig 4", "Contig 5", "Contig 6", "Contig 7","Contig 8",
                              "Contig 9", "Contig 10", "Contig 0 (11)", ""))+
  #scale_y_discrete(limits=specified_order2)+
  theme(panel.grid.minor.x = element_blank(), #axis.text.x=element_blank(), #gets rid of the minor grid lines
        panel.grid.major.x=element_line(color=rgb(0,0,0,0.4)), #turns the x axis gridlines gray
        axis.text.x= element_text(hjust=0, size=8),
        #legend.title=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5))+
  guides(color = guide_legend(override.aes = list(size =3, alpha= 0.7)))+ #changes the legend to make the dots larger and not transparent
  ggtitle("Species specific loci genome distribution") +
  labs(x=expression(paste("Contig in ",italic("C. florida"), " reference genome"))) +
  ylab("Breeding Selection")+
  geom_point(data=subset(dataforgraphing1,distinguishing_locus %in% "N"), 
             aes(x=overall_bp, y=new_name, color="blue"),
             size=1, alpha=0.2, position=position_nudge(y=0.17))+
  geom_point(data=subset(dataforgraphing1,distinguishing_locus %in% "F"), #gives the yellow dots a black outline to make them easier to see
             aes(x=overall_bp, y=new_name),
             size=1.1, alpha=0.9, color="black", position=position_nudge(y=-0.165))+
  geom_point(data=subset(dataforgraphing1,distinguishing_locus %in% "F"), 
             aes(x=overall_bp, y=new_name, color="yellow"),
             size=1, alpha=1, position=position_nudge(y=-0.17))+
  scale_color_identity(guide = "legend", #forces the aes to recognize the color names when they are inside aes
                       name= "Species Specific Loci",
                       labels=c("C. nuttallii   ", "C. kousa   ", "C. florida")) 
ggsave("Species specific loci3.tiff", width=22.23, height=19.05, units="cm", dpi=300)
ggsave("Fig16.tiff", width=22.23, height=19.05, units="cm", dpi=300)

#making dots smaller so that gaps are easier to see
ggplot()+
  geom_point(data=subset(dataforgraphing3,distinguishing_locus %in% "K"), 
             aes(x=overall_bp, y=new_name, color="red"),
             size=0.2, alpha=0.2)+ #position=position_jitter(h=0.15))+
  scale_x_continuous(breaks=c(0,129942925,252856830,359034702,479829781,607726880,730331874,829888132,935043946,1040321751,1165293952,1318504356), 
                     labels=c("Contig 1", "Contig 2", "Contig 3", "Contig 4", "Contig 5", "Contig 6", "Contig 7","Contig 8",
                              "Contig 9", "Contig 10", "Contig 0 (11)", ""))+
  #scale_y_discrete(limits=specified_order2)+
  theme(panel.grid.minor.x = element_blank(), #axis.text.x=element_blank(), #gets rid of the minor grid lines
        panel.grid.major.x=element_line(color=rgb(0,0,0,0.4)), #turns the x axis gridlines gray
        axis.text.x= element_text(hjust=0),
        #legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5))+
  guides(color = guide_legend(override.aes = list(size =3, alpha= 0.7)))+ #changes the legend to make the dots larger and not transparent
  ggtitle("Species specific loci genome distribution") +
  labs(x=expression(paste("Contig in ",italic("C. florida"), " reference genome"))) +
  ylab("Breeding Selection")+
  geom_point(data=subset(dataforgraphing3,distinguishing_locus %in% "N"), 
             aes(x=overall_bp, y=new_name, color="blue"),
             size=0.2, alpha=0.2, position=position_nudge(y=0.17))+
  geom_point(data=subset(dataforgraphing3,distinguishing_locus %in% "F"), #gives the yellow dots a black outline to make them easier to see
             aes(x=overall_bp, y=new_name),
             size=0.3, alpha=0.2, color="black", position=position_nudge(y=-0.165))+
  geom_point(data=subset(dataforgraphing3,distinguishing_locus %in% "F"), 
             aes(x=overall_bp, y=new_name, color="yellow"),
             size=0.2, alpha=0.2, position=position_nudge(y=-0.17))+
  scale_color_identity(guide = "legend", #forces the aes to recognize the color names when they are inside aes
                       name= "Species Specific\nLoci",
                       labels=c("C. nuttallii", "C. kousa", "C. florida")) 
ggsave("Species specific loci smaller.tiff", width=30, height=15, units="cm", dpi=300)
