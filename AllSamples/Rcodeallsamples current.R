setwd("")

cornus.raw3009=read.table("3009SNPsR.txt", row = 1, skipNul = TRUE)
library(vegan)
cAll.dist3009 <- vegdist(cornus.raw3009, method = "gower", na.rm = TRUE)


#Doing PCO with ape package
library(ape)
cAll.pco3009ape <- pcoa(cAll.dist3009)
#do not use one of the pre-defined corrections for negative eigenvalues, instead, treat them as zeros by only adding up the positive eigenvalues
pco1.3009ape <- cAll.pco3009ape$values$Eigenvalues[1]/sum(cAll.pco3009ape$values$Eigenvalues[1:sum(cAll.pco3009ape$values$Eigenvalues>0)]) #only specify the positive values
pco2.3009ape <- cAll.pco3009ape$values$Eigenvalues[2]/sum(cAll.pco3009ape$values$Eigenvalues[1:sum(cAll.pco3009ape$values$Eigenvalues>0)]) #only specify the positive values
pco1.3009ape
pco2.3009ape

plot(cAll.pco3009ape$vectors[,1], cAll.pco3009ape$vectors[,2], type = "p", xlab = "PCO1", ylab = "PCO2", axes = TRUE, main = "PCO All Samples")
#add small labels that make the PCO too busy in order to see trends
plot(cAll.pco3009ape$vectors[,1], cAll.pco3009ape$vectors[,2], type = "n", xlab = "PCO1", ylab = "PCO2", axes = TRUE, main = "PCO All Samples")
text(cAll.pco3009ape$vectors[,1], cAll.pco3009ape$vectors[,2], labels(cAll.pco3009ape$vectors[,1]), cex = 0.1, xpd = TRUE)


#printing out high quality 
library(ggplot2)
cAll.pcogroups=scan("GroupsforPCO.txt")
cAll.pcogroups2 <- as.factor(cAll.pcogroups) #convert groups of bract color into factors so can be read by ggplot
cAll.pco3009dataframe <- data.frame(cAll.pco3009ape$vectors[,1], cAll.pco3009ape$vectors[,2], cAll.pcogroups) #make the eigenvectors and PCO groupsinto a data frame
ggplot(cAll.pco3009dataframe, aes(x=cAll.pco3009ape$vectors[,1], y=cAll.pco3009ape$vectors[,2], color=cAll.pcogroups2, fill=cAll.pcogroups2)) +
geom_point(size=4, shape=21)+
  scale_fill_manual(values=c(rgb(1, 1, 0, 0.3), rgb(1, 0, 0, 0.3), rgb(0, 0, 1, 0.1), rgb(0, 0, 0, 0.3), rgb(0.4, 0, 0.6, .3)),
                    labels=c("C. florida", "C. kousa", "C. nuttallii", "other outgroups", "hybrids")) + #color scale for fill
  scale_color_manual(values=c(rgb(0, 0, 0, 0.02), rgb(1, 0, 0, 1), rgb(0, 0, 1, 1), rgb(0, 0, 0, 1), rgb(0.4, 0, 0.6, 1)),
                     labels=c("C. florida", "C. kousa", "C. nuttallii", "other outgroups", "hybrids")) + #color scale for outline
  theme_bw() + #removes gray background
  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), #gets rid of the minor grid lines
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), #gets rid of the major grid lines
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5)) + #centers title
  guides(fill = guide_legend(override.aes = list(size = 4))) + #changes the legend to make the dots larger
  coord_fixed(ratio=1)+ #makes sure the axes have the same scale
  ggtitle("PCO All Accessions") +
  xlab("PCO1 95.1%") +
  ylab("PCO2 1.8%")
ggsave("PCO Cornus All.tiff", scale=1, width=10, units="in", dpi=600)


#####Now visualize the structure data here 
remotes::install_github('royfrancis/pophelper') #pophelper install is slightly different because not in CRAN
library("pophelper")
library(gridExtra)

sfiles <- list.files(path="pophelper", full.names=TRUE) #read all of the files in a folder
sfiles

# basic usage
slist <- readQ(files=sfiles,indlabfromfile =T) #the previous can be used if names are not cut off in the output file
slist

p1 <- plotQ(qlist=slist[c(3,4,10)],imgoutput="join",returnplot=TRUE,exportplot=FALSE,basesize=11,exportpath =getwd())
grid.arrange(p1$plot[[1]])

#careful because the lab legend isn't necessarily in order
plotQ(slist[1],returnplot=TRUE,exportplot=TRUE,
      showtitle=T, titlelab="Structure all accessions", titlesize=5,
      sortind="all",
      clustercol=c("yellow","red"), height=5, width=10, barsize=1, barbordersize = 0,
      showlegend=T, legendlab=c("florida   ","kousa"), legendkeysize=3, legendtextsize=3, 
      showsp=TRUE, splab="K=2", splabsize=3, 
      showindlab=TRUE, useindlab=TRUE, indlabsize=2,
      outputfilename="plotq",imgtype="png", dpi=600, exportpath=getwd()
)

#plot files with individuals in the same order as in the PCO (first coordinate)
sfile <- "BASIC file for pophelper in PC1 order.txt"
sfile

inorder <- readQ(files=sfile)
inorder

plotQ(inorder,returnplot=TRUE,exportplot=TRUE,
      showtitle=T, titlelab="Structure all accessions", titlesize=7,
      clustercol=c("yellow","red"), height=4, width=10, barsize=1, barbordersize = 0,
      showlegend=T, legendlab=c("C. florida   ","C. kousa"), legendkeysize=5, legendtextsize=5, 
      showsp=TRUE, splab="K=2", splabsize=6,
      showindlab=FALSE, useindlab=TRUE, indlabsize=2,
      outputfilename="StructureplotinPCOorder",imgtype="tiff", dpi=600, exportpath=getwd()
)

#########################Combine the two plots into one graphic for printing#####################################

Structureall <-plotQ(inorder,returnplot=TRUE,exportplot=FALSE,
      showtitle=T, titlelab="Structure All Accessions", titlesize=16, titlecol= "black",
      clustercol=c("yellow","red"), height=4, width=10, barsize=1, barbordersize = 0,
      showlegend=T, legendlab=c("C. florida   ","C. kousa"), legendkeysize=13, legendtextsize=13, 
      showsp=TRUE, splab="K=2", splabsize=14, titlehjust=0.45, splabcol="black",
      showindlab=FALSE, useindlab=TRUE, indlabsize=11, indlabcol="black",
      showyaxis=TRUE
)
Structureall

PCOall <-ggplot(cAll.pco3009dataframe, aes(x=cAll.pco3009ape$vectors[,1], y=cAll.pco3009ape$vectors[,2], color=cAll.pcogroups2, fill=cAll.pcogroups2)) +
  geom_point(size=4, shape=21)+
  scale_fill_manual(values=c(rgb(1, 1, 0, 0.8), rgb(1, 0, 0, 0.3), rgb(0, 0, 1, 0.1), rgb(0, 0, 0, 0.3), rgb(0.4, 0, 0.6, .3)),
                    labels=c("C. florida", "C. kousa", "C. nuttallii", "other outgroups", "hybrids")) + #color scale for fill
  scale_color_manual(values=c(rgb(0, 0, 0, 0.3), rgb(1, 0, 0, 1), rgb(0, 0, 1, 1), rgb(0, 0, 0, 1), rgb(0.4, 0, 0.6, 1)),
                     labels=c("C. florida", "C. kousa", "C. nuttallii", "other outgroups", "hybrids")) + #color scale for outline
  theme_bw() + #removes gray background
  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), #gets rid of the minor grid lines
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), #gets rid of the major grid lines
        legend.title=element_blank(),
        legend.text = element_text(size=12),
        axis.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5, size=16)) + #centers title
  guides(fill = guide_legend(override.aes = list(size = 4))) + #changes the legend to make the dots larger
  coord_fixed(ratio=1)+ #makes sure the axes have the same scale
  ggtitle("PCO All Accessions") +
  xlab("PCO1 95.1%") +
  ylab("PCO2 1.8%")
PCOall

library(ggpubr)
library(gridExtra)
arrangedplotsall <-ggarrange(PCOall, NULL, Structureall$plot[[1]],
                          labels = c("a","","b"), heights = c(1, 0.05, 1), #created an extra row in order to add space between the graphs
                          nrow=3)
arrangedplotsall
ggexport(arrangedplotsall, filename="PCO+structure all samplesdarkenedyellow.tiff", width=6800, height=6000, res=600)
ggexport(arrangedplotsall, filename="PCO+structure all samples res300.tiff", width=3400, height=3000, res=300)
