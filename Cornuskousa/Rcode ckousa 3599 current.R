setwd("")
ckousa.raw3599=read.table("3599markerskousawoutdoublesforNJ.txt", row = 1, skipNul = TRUE)

library(ape)
library(vegan)
ckousa.dist3599 <- vegdist(ckousa.raw3599, method = "gower", na.rm = TRUE)
ckousa.nj3599 <-nj(ckousa.dist3599)
plot(ckousa.nj3599, cex = 0.7, label.offset = 0.001)



#note, in order to present a rooted tree with bootstrap values, one needs to reroot it before doing bootstrapping, 
# That way the bootstraps are associated with the right node
#in this case I did not designate an outgroup
f <- function(x) nj(vegdist(x, method = "gower", na.rm = TRUE)) #first define the function used to make the distance matrix and phylogram
bpckousa3599 <- boot.phylo(ckousa.nj3599, ckousa.raw3599, f, B=1000)

bpckousa3599.1 <- replace(bpckousa3599, bpckousa3599<600, "") #gets rid of bp values that are less than 600
nodelabels(bpckousa3599.1, adj = c(1.1, -0.4), frame = "none", cex = 0.2, col = "green4")
ckousa.nj3599outgroup.ladder <- ladderize(ckousa.nj3599) #how to reorder the tree so it is a little easier to see
plot(ckousa.nj3599outgroup.ladder, cex = 0.2, label.offset = 0.001, lwd=0.1) #lwd is the line thickness, but it doesn't seem to work, figure out later
add.scale.bar(cex=.5) #add a scale bar

#Dendrogram- Add circles of showing pink vs white, and wild vs garden origin
#also, plot alongside a barplot showing the structure results
#additionally, add another plot showing the k=means results
ckousa.pcogroups=scan("NJgroupsforcolorWOduplicates.txt")
ckousa.pcogroupswild=scan("NJgroupswildWOduplicates.txt")
ckousa.kmeans=read.table("kmeansclusterskousa.txt", header=FALSE, row.names=1, stringsAsFactors = FALSE)
StructureResults=read.table("StructureResultsforRkousa.txt", header=FALSE, row.names=1, stringsAsFactors = FALSE)
StructureResults2 <- matrix(unlist(StructureResults), ncol=3, nrow=157) #converting from a list to a matrix
typeof(StructureResults2) #should return "double"
StructureResults2
rownames(StructureResults2) <- rownames(StructureResults) #re-assigning row names, need to make sure the number of accessions exactly matches
ckousa.kmeans2 <- matrix(unlist(ckousa.kmeans), ncol=1, nrow=157)
rownames(ckousa.kmeans2) <- rownames(ckousa.kmeans) #re-assigning row names, need to make sure the number of accessions exactly matches

ckousa.kmeans2
ckousa.nj3599outgroup.ladder

#this will print directly to a tiff file and not show in the plots window
#may have to tweak sizes on final graph in order to make everything legible
tiff("Phyloplot+structure+Kmeans.tif", width = 3600, height = 6000, units = "px", res = 300)
plot(ckousa.nj3599outgroup.ladder, "ph", cex = .7, label.offset = .013, x.lim = 1)
phydataplot(StructureResults2, ckousa.nj3599outgroup.ladder, offset=0.173, scaling=0.2, col=c("#33a02c", "#6a3d9a", rgb(0,0,0,0)), border=NA, main="Structure Results")
f <- function(n) c("#6a3d9a", "#1f78b4", "#33a02c", "#fb9a99", "#e31a1c", "#b2df8a")
#c("#33a02c", "#1f78b4", "#6a3d9a")
phydataplot(ckousa.kmeans2, ckousa.nj3599outgroup.ladder, "mosaic", offset=0.4, width=0.05, funcol=f, border="white", legend="none")
title(main = "NJ Tree + Structure + K means clustering", line= -.1) #play around with this more to put the titles in the right place, or maybe add later with a different program.
nodelabels(bpckousa3599.1, adj = c(1.1, -0.2), frame = "none", cex = .7, col = "green4")
tiplabels(type = "p", pch = c(21,24,22)[ckousa.pcogroupswild],
          bg = c(rgb(0,0,0,0.3), "#fb9a99")[ckousa.pcogroups], 
          col = c(rgb(0,0,0,1), "#e31a1c")[ckousa.pcogroups], offset= 0.007, cex = 1)
legend("bottomleft", legend = c("White-bracted", "Pink-bracted", "", "Cultivars", "Rutgers Material", "Wild-collected"), 
       bty='n', cex = 1, pch = c(21,21,21,21,24,22), 
       pt.bg = c(rgb(0,0,0,0.3), "#fb9a99", rgb(0, 0, 0, 0), rgb(0, 0, 0, 1), rgb(0, 0, 0, 1), rgb(0, 0, 0, 1)), 
       col = c(rgb(0,0,0,1), "#e31a1c", rgb(0, 0, 0, 0), rgb(0, 0, 0, 1), rgb(0, 0, 0, 0), rgb(0, 0, 0, 1)), pt.cex = 2, y.intersp = 1)
add.scale.bar(x=0.2,y=1, cex=1)
dev.off()

plot(ckousa.nj3599outgroup.ladder, "ph", cex = .7, label.offset = .013, x.lim = 1)

########Do the PCO with non-duplicate accessions###########
ckousa.raw3599WOduplicates=read.table("3599markerskousawoutduplicatesforPCO.txt", row = 1, skipNul = TRUE)

library(vegan)
ckousa.distwodup3599 <- vegdist(ckousa.raw3599WOduplicates, method = "gower", na.rm = TRUE)

#making pco with ape
library(ape)
ckousa.pco3599ape <- pcoa(ckousa.distwodup3599)
#do not use one of the pre-defined corrections for negative eigenvalues, instead, treat them as zeros by only adding up the positive eigenvalues
pco1.3599ape <- ckousa.pco3599ape$values$Eigenvalues[1]/sum(ckousa.pco3599ape$values$Eigenvalues[1:sum(ckousa.pco3599ape$values$Eigenvalues>0)]) #only specify the positive values
pco2.3599ape <- ckousa.pco3599ape$values$Eigenvalues[2]/sum(ckousa.pco3599ape$values$Eigenvalues[1:sum(ckousa.pco3599ape$values$Eigenvalues>0)]) #only specify the positive values
pco3.3599ape <- ckousa.pco3599ape$values$Eigenvalues[3]/sum(ckousa.pco3599ape$values$Eigenvalues[1:sum(ckousa.pco3599ape$values$Eigenvalues>0)]) #only specify the positive values
pco1.3599ape
pco2.3599ape
pco3.3599ape

plot(ckousa.pco3599ape$vectors[,1], ckousa.pco3599ape$vectors[,2], type = "p", xlab = "PCO1", ylab = "PCO2", axes = TRUE, main = "PCO All Samples")
#add small labels that make the PCO too busy in order to see trends
plot(ckousa.pco3599ape$vectors[,1], ckousa.pco3599ape$vectors[,2], type = "n", xlab = "PCO1", ylab = "PCO3", axes = TRUE, main = "PCO All Samples")
text(ckousa.pco3599ape$vectors[,1], ckousa.pco3599ape$vectors[,2], labels(ckousa.pco3599ape$vectors[,1]), cex = 0.1, xpd = TRUE)

#Use PCO plot with ggplot to make it easier to print as high quality TIFF
# Adding pink astericks for the pink bracted dogwoods
ckousa.pcogroups2=scan("PCOgroupsforcolorWOduplicates.txt")
ckousa.pcogroupswild2=scan("PCOgroupswildWOduplicates.txt")
ckousa.pcogroupspink=scan("PCOgroupspinkorwhiteWOduplicates.txt")
ckousa.pcogroups3 <- as.factor(ckousa.pcogroups2) #convert groups of bract color into factors so can be read by ggplot
ckousa.pcogroups4 <- as.factor(ckousa.pcogroupswild2)
ckousa.pco3599dataframe <- data.frame(ckousa.pco3599ape$vectors[,1], ckousa.pco3599ape$vectors[,2],ckousa.pco3599ape$vectors[,3], ckousa.pcogroups3, ckousa.pcogroups4, ckousa.pcogroupspink) #make the eigenvectors and PCO groupsinto a data frame
library(ggplot2)
ggplot()+
  #overlay pink data to the graph
  geom_point(data=ckousa.pco3599dataframe, 
             aes(x=ckousa.pco3599ape.vectors...1., y=ckousa.pco3599ape.vectors...2., color=as.factor(ckousa.pcogroups3), fill=as.factor(ckousa.pcogroups3), shape=as.factor(ckousa.pcogroups4)),
             size=3)+
  scale_fill_manual(values=c("#33a02c33", "#1f78b433", "#6a3d9a33"),
                    labels=c("ssp. chinensis", "ssp. hybrid", "ssp. kousa")) + #color scale for fill
  scale_color_manual(values=c("#33a02c", "#1f78b4", "#6a3d9a"),
                     labels=c("ssp. chinensis", "ssp. hybrid", "ssp. kousa")) + #color scale for outline
  scale_shape_manual(values=c(21,24,22), labels=c("Cultivars", "Rutgers material", "Wild-collected"))+ #three different shapes depending on ckousa.pcogroups4
  theme_bw() + #removes gray background
  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), #gets rid of the minor grid lines
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), #gets rid of the major grid lines
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5)) + #centers title
  coord_fixed(ratio=1)+ #makes sure the axes have the same scale
  #guides(fill = guide_legend(override.aes = list(size = 4, order=2)))+ #order=1)) + #changes the legend to make the dots larger and put the colors before shapes
  #fill=FALSE)+ #Takes out the fill menu from the side
  scale_y_continuous(breaks=c(-0.1, 0.0, 0.1))+
  ggtitle("PCO Cornus kousa") +
  xlab("PCO1 49.1%") +
  ylab("PCO2 8.4%")+
  geom_point(data=subset(ckousa.pco3599dataframe,as.factor(ckousa.pcogroupspink)=="1"),
             aes(x=ckousa.pco3599ape.vectors...1., y=ckousa.pco3599ape.vectors...2.), size=3, color="deeppink", shape = "*")
guides(color=guide_legend(override.aes= list(fill=1)))
ggsave("PCO Cornus kousa shapeswpink.tiff", scale=1, width=8.5, units="in", dpi=600)


# unfortunately have to make a separate plot just with the legend, use it to copy and paste onto other graph
ggplot()+
  #overlay pink data to the graph
  geom_point(data=ckousa.pco3599dataframe,
             aes(x=ckousa.pco3599ape.vectors...1., y=ckousa.pco3599ape.vectors...2., color=as.factor(ckousa.pcogroupspink)), size=3, shape="*")+
  scale_color_manual(values=c("deeppink", "black"),
                   labels=c("Pink-bracted", "White-bracted")) + #color scale for outline
  scale_shape_manual(values=c("*",""), labels=c("Pink-bracted", "White-bracted"))+ #two different shapes
  theme_bw() + #removes gray background
  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), #gets rid of the minor grid lines
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), #gets rid of the major grid lines
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5)) + #centers title
  coord_fixed(ratio=1)+ #makes sure the axes have the same scale
  #guides(fill = guide_legend(override.aes = list(size = 4, order=2)))+ #order=1)) + #changes the legend to make the dots larger and put the colors before shapes
  #fill=FALSE)+ #Takes out the fill menu from the side
  scale_y_continuous(breaks=c(-0.1, 0.0, 0.1))+
  ggtitle("PCO Cornus kousa") +
  xlab("PCO1 49.1%") +
  ylab("PCO2 8.4%")
ggsave("PCO Cornus kousa for legend.tiff", scale=1, width=8.5, units="in", dpi=600)


#how to extract the names of a ladderized tree, in order
#First step is to filter out internal nodes from the the second column of the edge matrix:
is_tip <- ckousa.nj3599outgroup.ladder$edge[,2] <= length(ckousa.nj3599outgroup.ladder$tip.label)
ordered_tips <- ckousa.nj3599outgroup.ladder$edge[is_tip, 2]
#Then you can use this vector to extract the tips in the right order:
ckousa.nj3599outgroup.ladder$tip.label[ordered_tips]
#source https://stackoverflow.com/questions/34364660/how-to-get-correct-order-of-tip-labels-in-ape-after-calling-ladderize-function
#write the vector to a file
write(ckousa.nj3599outgroup.ladder$tip.label[ordered_tips], "ckousaNamesLadderized.txt")

#Doing an AMOVA with poppr
#It uses adegenet and genind objects
#first of all, how to import the strata correctly, from https://github.com/thibautjombart/adegenet/blob/master/tutorials/tutorial-strata.pdf
library(poppr)
ckousa.raw3599poppr=read.table("3599markerskousapopper.txt", skipNul = TRUE)
kousa.strata=read.table("KousaStrata.txt", skipNul = TRUE)
#convert into genind object
ckousa.3599genindpoppr <- df2genind(ckousa.raw3599poppr[,-1], ploidy=2, sep="", pop=ckousa.raw3599poppr$Pop, strata=kousa.strata) # [,-1]ignore the column with the pop data in the data frame initially
ckousa.3599genindpoppr

#change genind into genclone object
ckousa.3599genclone <- as.genclone(ckousa.3599genindpoppr)
ckousa.3599genclone

#perform AMOVA
#used the default parameters, so only used markers with 5 percent or less missing data
ckousa.3599amova <- poppr.amova(ckousa.3599genclone, ~Ssp) #note this either ignores or changes loci with greater than 0.05 missing data
ckousa.3599amova
ckousa.3599amova0.02 <- poppr.amova(ckousa.3599genclone, ~Ssp, cutoff = 0.02) #note this either ignores or changes loci with greater than 0.02 missing data
ckousa.3599amova0.02
ckousa.3599amova0.1 <- poppr.amova(ckousa.3599genclone, ~Ssp, cutoff = 0.1) #note this either ignores or changes loci with greater than 0.10 missing data
ckousa.3599amova0.1
ckousa.3599amova0.15 <- poppr.amova(ckousa.3599genclone, ~Ssp, cutoff = 0.15) #note this either ignores or changes loci with greater than 0.15 missing data
ckousa.3599amova0.15

#both cutoff = 0.05 (default), 0.02, 0.1, and 0.15 give similar answers
#use the default

#significance testing by randomly permuting the sample matrix
set.seed(1999)
ckousa.3599amovasignif   <- randtest(ckousa.3599amova, nrepet = 999)
plot(ckousa.3599amovasignif)
ckousa.3599amovasignif

##################next analysis, parentage
#using the apparent package from Iago and Hale to see if it can successfully identify parents of individuals, with confidence assignments.

setwd("./apparent")

library(outliers)
# Load the input file
InputFile <- read.table(file="3599 markers for apparent.txt",sep="\t",h=FALSE)
#make sure to run "apparent_script_Pavel for printing out triad plot" before this step
CkapparentOUT <- apparent(InputFile, MaxIdent=0.10, nloci=300, self=FALSE, plot=TRUE, Dyad=TRUE)
#there is an error with Dyad=TRUE, which is unfortunate, because that is my main interest with this r package
#actually, the error is fixed with Pavel's script

#to check out the results
CkapparentOUT

InputFile

#printing out part of the results to a CSV in order to see all of them in excel
write.csv(CkapparentOUT$Triad_all,"apparentOUT_Triad_all.csv", row.names = FALSE)
write.csv(CkapparentOUT$Dyad_sig,"apparentOUT_Dyad_sig2.csv", row.names = FALSE)

##################
#do dyad analysis in another way, using a similar strategy as the triad analysis in apparent, using proportion of loci in opposite homozygous states
#import distance matrix data to lay it flat

#importing the data as is, then melting it to get it to "lay it flat"
distancedataraw <- read.csv("3599 markers distance matrix for dyad analysis to melt.csv", check.names=FALSE)
distancedataraw
#note: it seems to be ok that there are both NA and <NA> in the data frame, they both show up as TRUE when is.na(dataraw) is run

class(distancedataraw) #should return a data.frame
library(reshape2)
distancemolten.data <- melt(distancedataraw, id= "row", na.rm = TRUE) #melt the data
distancemolten.data
class(distancemolten.data) #should return a data.frame

#Sort by value (the distance between the two, calculated as the number of homozygous mismatches
# divided by the total number of loci in common between the two)
attach(distancemolten.data)
sorted_distance_data <- distancemolten.data[order(value),]
sorted_distance_data
detach(distancemolten.data)

write.csv(sorted_distance_data,"dyad distance data sorted.csv", row.names = FALSE)

#add a column to spell out order of the Pohl values for plotting
sorted_distance_data_1 <- sorted_distance_data
#create a vector from 1 to the number of rows in the data frame
vec <- 1:nrow(sorted_distance_data)
sorted_distance_data_1$Pohlvalueorder <-vec

sorted_distance_data_1
ggplot(sorted_distance_data_1, aes(x=Pohlvalueorder, y=value)) + geom_point(shape=1)

#add Off, Pa, and All values to the data frame, then delete pairs with two Pa or Two offspring to decrease the search space (as with triads)

lookup <- read.csv("Categories for dyad analysis.csv", check.names=FALSE)
library(plyr)
names(lookup) <- c("row", 'Category1', 'Species1')
sorted_distance_data_2 <- join(sorted_distance_data_1, lookup, by = 'row')
names(lookup) <- c("variable", 'Category2', 'Species2')
sorted_distance_data_2 <- join(sorted_distance_data_2, lookup, by = 'variable')

#remove rows with Pa and Pa and also rows with Off and Off
sorted_distance_data_2 <- sorted_distance_data_2[!(sorted_distance_data_2[,5]=="Pa" & sorted_distance_data_2[,7]=="Pa"),]
sorted_distance_data_2 <- sorted_distance_data_2[!(sorted_distance_data_2[,5]=="Off" & sorted_distance_data_2[,7]=="Off"),]

#remove all the hybrids
sorted_distance_data_2 <-sorted_distance_data_2[sorted_distance_data_2$Species1 != "H" & sorted_distance_data_2$Species2 != "H", ]
#redo the Pohlvalueorder (fourth column)
sorted_distance_data_2 <-sorted_distance_data_2[ -c(4) ]
vec <- 1:nrow(sorted_distance_data_2)
sorted_distance_data_2$Pohlvalueorder <-vec

ggplot(sorted_distance_data_2, aes(x=Pohlvalueorder, y=value)) + geom_point(shape=1)

#do in ggplot to combine with the triad analysis plot
CkDyadPlot <- ggplot(sorted_distance_data_2, aes(x=Pohlvalueorder, y=value))+
  geom_point(size=0.5, shape=1)+
  labs(title="Cornus kousa Dyad Analysis Plot",
       x="Test dyads, ordered by POHL", y = "Proportion of loci with opposite \nheterozygous states (POHL)")+
  theme_bw() + 
  theme(panel.border = element_blank(), legend.position="none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks=seq(0,0.35,0.05))+
  scale_x_continuous(breaks=c(1, nrow(sorted_distance_data_2)))

#now the triad plot
SortGD1Ck <-CkapparentOUT$SortGD

vec <- 1:nrow(CkapparentOUT$SortGD)
CkapparentOUT$SortGD$Order <-vec

#replicate triad plot in ggplot in order to print out with the dyad analysis figure
CkTriadPlot <- ggplot(CkapparentOUT$SortGD, aes(x=Order, y=GD, color=Colour))+
  geom_point(size=0.5, shape=1)+
  scale_color_manual(values=c('Black','Red'))+
  labs(title="Cornus kousa Triad Analysis Plot",
       x="Test triads, ordered by GDij|POk", y = "Gower Genetic Dissimilarity (GDij|k)")+
  theme_bw() + 
  theme(panel.border = element_blank(), legend.position="none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks=seq(0,0.60,0.05))+
  scale_x_continuous(breaks=c(1, nrow(CkapparentOUT$SortGD)))+
  geom_hline(yintercept=CkapparentOUT$ThresholdT, linetype="dashed", color = "red")

#put the triad and the dyad plots together
library(ggpubr)
arranged <-ggarrange(CkTriadPlot, CkDyadPlot,
                     labels = c("a", "b"),
                     nrow = 1)
arranged
ggexport(arranged, filename="Triad and dyad plots.tiff", width=5000, height=4000, res=600)


ggplot(sorted_distance_data_2, aes(x=value)) +geom_histogram(binwidth = 0.01)

Tdiff <- vector(mode="numeric",length=0)  
for (i in 1:nrow(sorted_distance_data_2)) {
  Tdiff[i] <- sorted_distance_data_2$value[i+1] - sorted_distance_data_2$value[i]
}
Tdiff
sort(Tdiff, decreasing=TRUE)
length(Tdiff)
which.max(Tdiff)
TMax <-0.003219615
#it is returning the gap at the top


#Searching only the first half of the data
halflength<-nrow(sorted_distance_data_2)/2

Tdiff <- vector(mode="numeric",length=0)  
for (i in 1:halflength) {
  Tdiff[i] <- sorted_distance_data_2$value[i+1] - sorted_distance_data_2$value[i]
}
Tdiff
sort(Tdiff, decreasing=TRUE)
length(Tdiff)
which.max(Tdiff)
TMax <-max(Tdiff)


#set TIndex to the number of the largest gap in Tdiff (on the low end if there are larger gaps at the high end)
TIndex <- which.max(Tdiff)
Tvect1 <- c(sample(na.omit(Tdiff[-TIndex]),29,replace=T),TMax)
Tvect1
TDtGap <- dixon.test(Tvect1)
TDtGap
#if this is an outlier, than can proceed to the next test


TCutoff <- sorted_distance_data_2$value[TIndex]
TCutoff
L <- sorted_distance_data_2[which(sorted_distance_data_2$value <= TCutoff),]
H <- sorted_distance_data_2[which(sorted_distance_data_2$value > TCutoff),]
S <- H$value[1:29]

Tpv <- vector(mode="numeric",length=0)
for (i in 1:nrow(L)) { 
  Tvect2 <- c(S, L$value[i])
  TDt <- dixon.test(Tvect2)
  Tpv <- append (Tpv,TDt$p.value)
}
Tpv
L
class(Tpv)
format(Tpv, scientific = TRUE)
TDt


###testing examples to see if I can collapse half of the tree
library(ape)
## Not run:
tr <- rtree(20)
f <- function(col) {
  o <- identify(tr)
  nodelabels(node=o$nodes, pch = 19, col = col)
}
plot(tr)
f("red") # click close to a node
f("green")


data(bird.families)
tip <- c(
  "Eopsaltriidae", "Acanthisittidae", "Pittidae", "Eurylaimidae",
  "Philepittidae", "Tyrannidae", "Thamnophilidae", "Furnariidae",
  "Formicariidae", "Conopophagidae", "Rhinocryptidae", "Climacteridae",
  "Menuridae", "Ptilonorhynchidae", "Maluridae", "Meliphagidae",
  "Pardalotidae", "Petroicidae", "Irenidae", "Orthonychidae",
  "Pomatostomidae", "Laniidae", "Vireonidae", "Corvidae",
  "Callaeatidae", "Picathartidae", "Bombycillidae", "Cinclidae",
  "Muscicapidae", "Sturnidae", "Sittidae", "Certhiidae",
  "Paridae", "Aegithalidae", "Hirundinidae", "Regulidae",
  "Pycnonotidae", "Hypocoliidae", "Cisticolidae", "Zosteropidae",
  "Sylviidae", "Alaudidae", "Nectariniidae", "Melanocharitidae",
  "Paramythiidae","Passeridae", "Fringillidae")
plot(drop.tip(bird.families, tip))
plot(drop.tip(bird.families, tip, trim.internal = FALSE))
data(bird.orders)
plot(drop.tip(bird.orders, 6:23, subtree = TRUE))
plot(drop.tip(bird.orders, c(1:5, 20:23), subtree = TRUE))
plot(drop.tip(bird.orders, c(1:20, 23), subtree = TRUE))
plot(drop.tip(bird.orders, c(1:20, 23), subtree = TRUE, rooted = FALSE))
### Examples of the use of `root.edge'
tr <- read.tree(text = "(A:1,(B:1,(C:1,(D:1,E:1):1):1):1):1;")
drop.tip(tr, c("A", "B"), root.edge = 0) # = (C:1,(D:1,E:1):1);
drop.tip(tr, c("A", "B"), root.edge = 1) # = (C:1,(D:1,E:1):1):1;
drop.tip(tr, c("A", "B"), root.edge = 2) # = (C:1,(D:1,E:1):1):2;
drop.tip(tr, c("A", "B"), root.edge = 3) # = (C:1,(D:1,E:1):1):3;

