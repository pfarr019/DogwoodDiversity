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
ggsave("Fig2.tiff", scale=1, width=10, units="in", dpi=600)

#####plotting a general tree with all accessions################
library(ape)
library(vegan)
cAll.nj3009 <-nj(cAll.dist3009)
plot(cAll.nj3009, cex = 0.7, label.offset = 0.001)
ggsave("cAll NJ Tree.tiff", dpi=600)

cAll.nj3009outgroup <- root(cAll.nj3009, "Cc_Mountain_Moon", edgelabel=TRUE) #set outgroup
plot(cAll.nj3009outgroup, cex = 0.3, label.offset = 0.001) #plot the new diagram with the outgroup
#note, in order to present a rooted tree with bootstrap values, one needs to reroot it before doing bootstrapping, That way the bootstraps are associated with the right node
f <- function(x) nj(vegdist(x, method = "gower", na.rm = TRUE)) #first define the function used to make the distance matrix and phylogram
bpcAll3009 <- boot.phylo(cAll.nj3009outgroup, cornus.raw3009, f, B=1000, rooted = is.rooted(cAll.nj3009outgroup))
bpcAll3009

bpcAll3009.1 <- replace(bpcAll3009, bpcAll3009<600, "") #gets rid of bp values that are less than 600

cAll.nj3009outgroup.ladder <- ladderize(cAll.nj3009outgroup) #how to reorder the tree so it is a little easier to see
#printing out ladderized figure
tiff("NJall with bootstrapping.tiff", units="px", width=10000, height=10000, res=600)
plot(cAll.nj3009outgroup.ladder, cex = 0.3)
nodelabels(bpcAll3009.1, adj = c(1.1, -0.2), frame = "none", cex = 0.3, col = "green4")
add.scale.bar(cex=.5) #add a scale bar
dev.off()

#Dendrogram- Add cirles of the same color as the PCO chart

cAll.pcogroups=scan("GroupsforPCO.txt")
cAll.pcogroups2 <- as.factor(cAll.pcogroups) #convert groups of bract color into factors so can be read by ggplot


cflorida.pcogroups=scan("PCOgroups.txt")
cflorida.pcogroupswild=scan("PCOgroupswild.txt")
cflorida.kmeans=read.table("kmeansclusters.txt", header=FALSE, row.names=1, stringsAsFactors = FALSE)
StructureResults=read.table("StructureResultsforR.txt", header=FALSE, row.names=1, stringsAsFactors = FALSE)
StructureResults2 <- matrix(unlist(StructureResults), ncol=4, nrow=94) #converting from a list to a matrix
typeof(StructureResults2) #should return "double"
StructureResults2
rownames(StructureResults2) <- rownames(StructureResults) #re-assigning row names, need to make sure the number of accessions exactly matches
cflorida.kmeans2 <- matrix(unlist(cflorida.kmeans), ncol=1, nrow=94)
rownames(cflorida.kmeans2) <- rownames(cflorida.kmeans) #re-assigning row names, need to make sure the number of accessions exactly matches

#this will print directly to a tiff file and not show in the plots window
# For publication
tiff("Phyloplot.tif", units="px", width=10000, height=10000, res=600)
plot(cAll.nj3009outgroup.ladder, "ph", cex = 0.3, label.offset = .005, x.lim = 1)
nodelabels(bpcAll3009.1, adj = c(1.1, -0.2), frame = "none", cex = 0.3, col = "green4")
tiplabels(type = "p", pch=21,
          bg = c(rgb(1, 1, 0, 0.7), rgb(1, 0, 0, 0.3), rgb(0, 0, 1, 0.1), rgb(0, 0, 0, 0.3), rgb(0.4, 0, 0.6, .3))[cAll.pcogroups2],
          col = c(rgb(0, 0, 0, 0.03), rgb(1, 0, 0, 1), rgb(0, 0, 1, 1), rgb(0, 0, 0, 1), rgb(0.4, 0, 0.6, 1))[cAll.pcogroups2], offset= 0.003, cex = .65)
legend(0.03, 25, legend = c("C. florida", "C. kousa", "C. nuttallii", "other outgroups", "hybrids"),
       bty='n', cex = 1, pch = c(21, 21, 21, 21, 21), 
       pt.bg = c(rgb(1, 1, 0, 0.7), rgb(1, 0, 0, 0.3), rgb(0, 0, 1, 0.1), rgb(0, 0, 0, 0.3), rgb(0.4, 0, 0.6, .3)), 
       col = c(rgb(0, 0, 0, 0.2), rgb(1, 0, 0, 1), rgb(0, 0, 1, 1), rgb(0, 0, 0, 1), rgb(0.4, 0, 0.6, 1)), pt.cex = 2)
add.scale.bar(x=0.17,y=1.1, cex=1.1)
dev.off()
