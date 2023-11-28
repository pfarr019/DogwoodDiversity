setwd("")
cflorida.raw7622=read.table("7622markersFloridaWOduplicates.txt", row = 1, skipNul = TRUE)

library(ape)
library(vegan)
cflorida.dist7622 <- vegdist(cflorida.raw7622, method = "gower", na.rm = TRUE)
cflorida.nj7622 <-nj(cflorida.dist7622)
plot(cflorida.nj7622, cex = 0.7, label.offset = 0.001)

cflorida.nj7622outgroup <- root(cflorida.nj7622, "Cf_ssp._urbiniana_USNA", edgelabel=TRUE) #set outgroup
plot(cflorida.nj7622outgroup, cex = 0.7, label.offset = 0.001) #plot the new diagram with the outgroup
#note, in order to present a rooted tree with bootstrap values, one needs to reroot it before doing bootstrapping, That way the bootstraps are associated with the right node
f <- function(x) nj(vegdist(x, method = "gower", na.rm = TRUE)) #first define the function used to make the distance matrix and phylogram
bpcflorida7622 <- boot.phylo(cflorida.nj7622outgroup, cflorida.raw7622, f, B=1000, rooted = is.rooted(cflorida.nj7622outgroup))
bpcflorida7622

bpcflorida7622.1 <- replace(bpcflorida7622, bpcflorida7622<600, "") #gets rid of bp values that are less than 600
nodelabels(bpcflorida7622.1, adj = c(1.1, -0.2), frame = "none", cex = 0.7, col = "green4")
cflorida.nj7622outgroup.ladder <- ladderize(cflorida.nj7622outgroup) #how to reorder the tree so it is a little easier to see
plot(cflorida.nj7622outgroup.ladder, cex = 0.7, label.offset = 0.001)
add.scale.bar(cex=.5) #add a scale bar

#Dendrogram- Add cirles of the same color as the PCO chart
#also, plot alongside a barplot showing the structure results
#additionally, add another plot showing the k=means results
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
tiff("Phyloplot+structure+Kmeans.tif", width = 3600, height = 6000, units = "px", res = 300)
plot(cflorida.nj7622outgroup.ladder, "ph", cex = 1.05, label.offset = .013, x.lim = 1)
phydataplot(StructureResults2, cflorida.nj7622outgroup.ladder, offset=0.225, scaling=0.2, col=c(rgb(0, 1, 0, 1), rgb(1, 0, 0, 1), rgb(0, 0, 1, 1), rgb(0, 0, 0, 0)), border=NA, main="Structure Results")
f <- function(n) c(rgb(0, 0, 1, 0.7), rgb(1, 0, 0, 0.7))
phydataplot(cflorida.kmeans2, cflorida.nj7622outgroup.ladder, "mosaic", offset=0.5, width=0.05, funcol=f, border="white", legend="none")
# title(main = "NJ Tree                                                  Structure                 K means                      \n   clustering                     ", line= -1, adj=1)
title(main = "A.                           NJ Tree                                                  B.   Structure                 C. K means                      \n   clustering                     ", line= -1, adj=1)
nodelabels(bpcflorida7622.1, adj = c(1.1, -0.2), frame = "none", cex = .9, col = "green4")
tiplabels(type = "p", pch = c(21, 22)[cflorida.pcogroupswild],
          bg = c(rgb(0, 0, 1, 0.3), rgb(1, 0, 0, 0.5), rgb(1, 0, 0, 0.05), rgb(0, 1, 0, 0.1), rgb(0, 0, 0, .3))[cflorida.pcogroups], 
          col = c(rgb(0, 0, 1, 1), rgb(1, 0, 0, 1), rgb(1, 0, 0, 1), rgb(0, 1, 0, 1), rgb(0, 0, 0, 1))[cflorida.pcogroups], offset= 0.007, cex = 1.5)
legend(0.01, 10, legend = c("White-bracted", "Pink-bracted", "Blush pink-bracted", "ssp. urbiniana", "Unknown bract color", "", "Garden origin", "Wild-collected"), 
       bty='n', cex = 1, pch = c(21, 21, 21, 21, 21, 21, 21, 22), 
       pt.bg = c(rgb(0, 0, 1, 0.3), rgb(1, 0, 0, 0.5), rgb(1, 0, 0, 0.05), rgb(0, 1, 0, 0.1), rgb(0, 0, 0, .3), rgb(0,0,0,0), rgb(0, 0, 0, 1), rgb(0, 0, 0, 1)), 
       col = c(rgb(0, 0, 1, 1), rgb(1, 0, 0, 1), rgb(1, 0, 0, 1), rgb(0, 1, 0, 1), rgb(0, 0, 0, 1), rgb(0, 0, 0, 0), rgb(0, 0, 0, 1), rgb(0, 0, 0, 1)), pt.cex = 2, y.intersp = 1)
add.scale.bar(x=0.2,y=1, cex=1)
dev.off()

########Do the PCO with non-duplicate accessions###########
cflorida.raw7622WOduplicates=read.table("7622SNPsforRWOduplicatesforPCO.txt", row = 1, skipNul = TRUE)

library(vegan)
cflorida.dist7622 <- vegdist(cflorida.raw7622WOduplicates, method = "gower", na.rm = TRUE)

#making pco with ape
library(ape)
cflorida.pco7622ape <- pcoa(cflorida.dist7622)
#do not use one of the pre-defined corrections for negative eigenvalues, instead, treat them as zeros by only adding up the positive eigenvalues
pco1.7622ape <- cflorida.pco7622ape$values$Eigenvalues[1]/sum(cflorida.pco7622ape$values$Eigenvalues[1:sum(cflorida.pco7622ape$values$Eigenvalues>0)]) #only specify the positive values
pco2.7622ape <- cflorida.pco7622ape$values$Eigenvalues[2]/sum(cflorida.pco7622ape$values$Eigenvalues[1:sum(cflorida.pco7622ape$values$Eigenvalues>0)]) #only specify the positive values
pco1.7622ape
pco2.7622ape

plot(cflorida.pco7622ape$vectors[,1], cflorida.pco7622ape$vectors[,2], type = "p", xlab = "PCO1", ylab = "PCO2", axes = TRUE, main = "PCO All Samples")
#add small labels that make the PCO too busy in order to see trends
plot(cflorida.pco7622ape$vectors[,1], cflorida.pco7622ape$vectors[,2], type = "n", xlab = "PCO1", ylab = "PCO2", axes = TRUE, main = "PCO All Samples")
text(cflorida.pco7622ape$vectors[,1], cflorida.pco7622ape$vectors[,2], labels(cflorida.pco7622ape$vectors[,1]), cex = 0.5, xpd = TRUE)


#Use PCO plot with ggplot to make it easier to print as high quality TIFF
library(ggplot2)
cflorida.pcogroups2=scan("PCOgroupsforcolorWOduplicates.txt")
cflorida.pcogroupswild2=scan("PCOgroupswildWOduplicates.txt")
cflorida.pco7622dataframe <- data.frame(cflorida.pco7622ape$vectors[,1], cflorida.pco7622ape$vectors[,2], cflorida.pcogroups2, cflorida.pcogroupswild2) #make the eigenvectors and PCO groupsinto a data frame
cflorida.pcogroups3 <- as.factor(cflorida.pcogroups2) #convert groups of bract color into factors so can be read by ggplot
cflorida.pcogroups4 <- as.factor(cflorida.pcogroupswild2)

ggplot(cflorida.pco7622dataframe, aes(x=cflorida.pco7622ape.vectors...1., y=cflorida.pco7622ape.vectors...2., color=cflorida.pcogroups3, fill=cflorida.pcogroups3)) +
  geom_point(size=2.5, shape=21)+
  scale_fill_manual(values=c(rgb(0, 0, 1, 0.3), rgb(1, 0, 0, 0.5), rgb(1, 0, 0, 0.05), rgb(0, 1, 0, 0.1), rgb(0, 0, 0, .3)),
                    labels=c("White-bracted", "Pink-bracted", "Blush pink-bracted", "ssp. urbiniana", "Unknown bract color")) + #color scale for fill
  scale_color_manual(values=c(rgb(0, 0, 1, 1), rgb(1, 0, 0, 1), rgb(1, 0, 0, 1), rgb(0, 1, 0, 1), rgb(0, 0, 0, 1)),
                     labels=c("White-bracted", "Pink-bracted", "Blush pink-bracted", "ssp. urbiniana", "Unknown bract color")) + #color scale for outline
  theme_bw() + #removes gray background
  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), #gets rid of the minor grid lines
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), #gets rid of the major grid lines
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5)) + #centers title
  guides(fill = guide_legend(override.aes = list(size = 2.5))) + #changes the legend to make the dots larger
  coord_fixed(ratio=1)+ #makes sure the axes have the same scale
  scale_y_continuous(breaks=c(-0.1, 0.0, 0.1))+
  ggtitle("PCO Cornus florida") +
  xlab("PCO1 11.9%") +
  ylab("PCO2 4.3%")
ggsave("PCO Cornus florida2.tiff", height=3.5, width=5.2, units="in", dpi=600) #setting as the column width for publication


#how to extract the names of a ladderized tree, in order
#First step is to filter out internal nodes from the the second column of the edge matrix:
is_tip <- cflorida.nj7622outgroup.ladder$edge[,2] <= length(cflorida.nj7622outgroup.ladder$tip.label)
ordered_tips <- cflorida.nj7622outgroup.ladder$edge[is_tip, 2]
#Then use this vector to extract the tips in the right order:
cflorida.nj7622outgroup.ladder$tip.label[ordered_tips]
#source https://stackoverflow.com/questions/34364660/how-to-get-correct-order-of-tip-labels-in-ape-after-calling-ladderize-function
#write the vector to a file
write(cflorida.nj7622outgroup.ladder$tip.label[ordered_tips], "CfloridaNamesLadderized.txt")

##################next analysis, parentage
#using the apparent package from Iago and Hale to see if it can successfully identify parents of individuals, with confidence assignments.

library(outliers)
# Load the input file
InputFile <- read.table(file="apparent/7622 markers florida apparent.txt",sep="\t",h=FALSE) #note that the names are included in the excel file of the same name
#try with the test data
# InputFile <- read.table(file="apparent/Apparent_testdata.txt",sep="\t",h=FALSE)

library(outliers)
#make sure to run "apparent_script_Pavel for printing out triad plot" before this step
apparentOUT <- apparent(InputFile, MaxIdent=0.10, nloci=300, self=FALSE, plot=TRUE, Dyad=TRUE)
#there is an error with Dyad=TRUE, which is unfortunate, because that is my main interest with this r package
#actually, the error is fixed with Pavel's script
#but the triad analysis works if I remove Dyad=TRUE, see below
apparentOUT0.01 <- apparent(InputFile, MaxIdent=0.10, nloci=300, self=FALSE, plot=TRUE, Dyad=FALSE)
apparentOUT0.05 <- apparent(InputFile, MaxIdent=0.10, alpha=0.05, nloci=300, self=FALSE, plot=TRUE, Dyad=TRUE)

#to check out the results
apparentOUT
apparentOUT0.01
apparentOUT0.05

InputFile

#how to write out part of the results to a CSV in order to see all of them
write.csv(apparentOUT$Triad_all,"apparent/apparentOUT_Triad_all.csv", row.names = FALSE)
write.csv(apparentOUT$Dyad_sig,"apparent/apparentOUT_Dyad_sig.csv", row.names = FALSE)

#printing out triad plot 
tiff("apparent/triad analysis plot2.tiff", width=800, height=800)
par(mar=c(5, 6, 4, 2) + 0.1)
plot (apparentOUT$SortGD$GD,xlab="Test triads, ordered by GDij|POk",ylab="Gower Genetic Dissimilarity (GDij|k)",
      main="Cornus florida Triad analysis plot",pch=1,cex=1.5, cex.axis=2,col=apparentOUT$SortGD$Colour,xaxt="n", cex.lab=2, cex.main=3)
axis(1,at=c(1,nrow(apparentOUT$SortGD)),labels=c("1",nrow(apparentOUT$SortGD)),cex.axis=2)
abline(h = apparentOUT$ThresholdT, lty = 2, col =  "tomato", lwd = 2)
dev.off()

SortGD1 <-apparentOUT$SortGD

vec <- 1:nrow(apparentOUT$SortGD)
apparentOUT$SortGD$Order <-vec

#replicate plot in ggplot in order to print out with the dyad analysis figure
CfTriadPlot <- ggplot(apparentOUT$SortGD, aes(x=Order, y=GD, color=Colour))+
  geom_point(size=2, shape=1)+
  scale_color_manual(values=c('Black','Red'))+
  labs(title="Cornus florida Triad Analysis Plot",
       x="Test triads, ordered by GDij|POk", y = "Gower Genetic Dissimilarity (GDij|k)")+
  theme_bw() + 
  theme(panel.border = element_blank(), legend.position="none", panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks=seq(0,0.35,0.05))+
  scale_x_continuous(breaks=c(1, nrow(apparentOUT$SortGD)))+
  geom_hline(yintercept=apparentOUT$ThresholdT, linetype="dashed", color = "red")

##################
#do dyad analysis in another way, using a similar strategy as the triad analysis in apparent, using proportion of loci in opposite homozygous states
#import distance matrix data to lay it flat
Cfdistancedataraw <- read.csv("apparent/7622 markers distance matrix for dyad analysis to melt.csv", check.names=FALSE)
Cfdistancedataraw
#note: it seems to be ok that there are both NA and <NA> in the data frame, they both show up as TRUE when is.na(dataraw) is run

class(Cfdistancedataraw) #should return a data.frame
library(reshape2)
Cfdistancemolten.data <- melt(Cfdistancedataraw, id= "row", na.rm = TRUE) #melt the data
Cfdistancemolten.data
class(Cfdistancemolten.data) #should return a data.frame

#Sort by value (the distance between the two, calculated as the number of homozygous mismatches
# divided by the total number of loci in common between the two)
attach(Cfdistancemolten.data)
Cfsorted_distance_data <- Cfdistancemolten.data[order(value),]
Cfsorted_distance_data
detach(Cfdistancemolten.data)

write.csv(Cfsorted_distance_data,"apparent/dyad distance data sorted.csv", row.names = FALSE)

#add a column to spell out order of the Pohl values for plotting
Cfsorted_distance_data_1 <- Cfsorted_distance_data
#create a vector from 1 to the number of rows in the data frame
vec <- 1:nrow(Cfsorted_distance_data)
Cfsorted_distance_data_1$Pohlvalueorder <-vec

#add Off, Pa, and All values to the data frame, then delete pairs with two Pa or Two offspring to decrease the search space (as with triads)

lookup <- read.csv("apparent/pa off all categories for dyad analysis.csv", check.names=FALSE)
library(plyr)
Cfsorted_distance_data_2 <- join(Cfsorted_distance_data_1, lookup, by = 'row')
names(lookup) <- c("variable", 'variable_type')
Cfsorted_distance_data_2 <- join(Cfsorted_distance_data_2, lookup, by = 'variable')

#remove rows with Pa and Pa and also rows with Off and Off
Cfsorted_distance_data_2 <- Cfsorted_distance_data_2[!(Cfsorted_distance_data_2[,5]=="Pa" & Cfsorted_distance_data_2[,6]=="Pa"),]
Cfsorted_distance_data_2 <- Cfsorted_distance_data_2[!(Cfsorted_distance_data_2[,5]=="Off" & Cfsorted_distance_data_2[,6]=="Off"),]

Cfsorted_distance_data_1
ggplot(Cfsorted_distance_data_1, aes(x=Pohlvalueorder, y=value)) + geom_point(shape=1)

Cfsorted_distance_data_2
ggplot(Cfsorted_distance_data_2, aes(x=Pohlvalueorder, y=value)) + geom_point(shape=1)

ggplot(Cfsorted_distance_data_2, aes(x=value)) +geom_histogram(binwidth = 0.01)

#using code repurposed from apparent triad analysis
Tdiff <- vector(mode="numeric",length=0)  
for (i in 1:nrow(Cfsorted_distance_data_2)) {
  Tdiff[i] <- Cfsorted_distance_data_2$value[i+1] - Cfsorted_distance_data_2$value[i]
}
Tdiff
sort(Tdiff, decreasing=TRUE)

#set TIndex to the number of the largest gap in Tdiff (on the low end if there are larger gaps at the high end)
TIndex <- 15
TMax=0.012839579 #set TMax to the size of the largest gap (on the low end if there are larger gaps at the high end)
Tvect1 <- c(sample(na.omit(Tdiff[-TIndex]),29,replace=T),TMax)
Tvect1
TDtGap <- dixon.test(Tvect1)
TDtGap
#if this is an outlier, than can proceed to the next test


TCutoff <- Cfsorted_distance_data_2$value[TIndex]
TCutoff
L <- Cfsorted_distance_data_2[which(Cfsorted_distance_data_2$value <= TCutoff),]
H <- Cfsorted_distance_data_2[which(Cfsorted_distance_data_2$value > TCutoff),]
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

CfThresholdT <- mean(c(Cfsorted_distance_data_2$value[TIndex],Cfsorted_distance_data_2$value[TIndex+1]))

#graphing
Cfsorted_distance_data_3 <- Cfsorted_distance_data_2

vec <- 1:nrow(Cfsorted_distance_data_2)
Cfsorted_distance_data_3$Colour <-vec

Cfsorted_distance_data_3 <- transform(
  Cfsorted_distance_data_3, Colour= ifelse(value<=TCutoff, "Red", "Black"))

#do in ggplot so that it can be combine with the triad analysis plot
CfDyadPlot <- ggplot(Cfsorted_distance_data_3, aes(x=as.numeric(Pohlvalueorder), y=value, color=Colour))+
  geom_point(size=2, shape=1)+
  scale_color_manual(values=c('Black','Red'))+
  labs(title="Cornus florida Dyad Analysis Plot",
       x="Test dyads, ordered by POHL", y = "Proportion of loci with opposite \nheterozygous states (POHL)")+
  theme_bw() + 
  theme(panel.border = element_blank(), legend.position="none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks=seq(0,0.35,0.05))+
  scale_x_continuous(breaks=c(1, nrow(Cfsorted_distance_data_2))) +
  geom_hline(yintercept=CfThresholdT, linetype="dashed", color = "red")
  

#put the triad and the dyad plots together
library(ggpubr)
Cfarranged <-ggarrange(CfTriadPlot, CfDyadPlot,
                          labels = c("a", "b"),
                          nrow = 1)
Cfarranged
ggexport(Cfarranged, filename="apparent/Triad and dyad plots.tiff", width=5000, height=2700, res=600)

  

#Taking out the ssp. urbiniana from spring meadow nursery to compare bootstrapping values with and without it

cflorida.raw7622WO=read.table("7622markersFloridaWOduplicatesWOurbinianaSM.txt", row = 1, skipNul = TRUE)
cflorida.dist7622WO <- vegdist(cflorida.raw7622WO, method = "gower", na.rm = TRUE)
cflorida.nj7622WO <-nj(cflorida.dist7622WO)
plot(cflorida.nj7622WO, cex = 0.7, label.offset = 0.001)

cflorida.nj7622WOoutgroup <- root(cflorida.nj7622WO, "Cf_ssp._urbiniana_USNA", edgelabel=TRUE) #set outgroup
plot(cflorida.nj7622WOoutgroup, cex = 0.7, label.offset = 0.001) #plot the new diagram with the outgroup
#note, if you want to present a rooted tree with bootstrap values, you need to reroot it before you do bootstrapping, That way the bootstraps are associated with the right node
f <- function(x) nj(vegdist(x, method = "gower", na.rm = TRUE)) #first define the function used to make the distance matrix and phylogram
bpcflorida7622WO <- boot.phylo(cflorida.nj7622WOoutgroup, cflorida.raw7622WO, f, B=1000, rooted = is.rooted(cflorida.nj7622WOoutgroup))

bpcflorida7622WO.1 <- replace(bpcflorida7622WO, bpcflorida7622WO<600, "") #gets rid of bp values that are less than 600
nodelabels(bpcflorida7622WO.1, adj = c(1.1, -0.2), frame = "none", cex = 0.7, col = "green4")
cflorida.nj7622WOoutgroup.ladder <- ladderize(cflorida.nj7622WOoutgroup) #how to reorder the tree so it is a little easier to see
plot(cflorida.nj7622WOoutgroup.ladder, cex = 0.7, label.offset = 0.001)
add.scale.bar(cex=.5) #add a scale bar

#writing dist object to a csv file
mat2 <- as.matrix(cflorida.dist7622)
mat2[upper.tri(mat2, diag = FALSE)] <- ""
write.csv(mat2, "mat2.csv")

#print out dendrogram without SMN urbiniana for inclusion in paper
cflorida.pcogroups=scan("PCOgroupsWOSMN.txt")
cflorida.pcogroupswild=scan("PCOgroupswildWOSMN.txt")
tiff("PhyloplotwoutSMN.tif", width = 3600, height = 6000, units = "px", res = 300)
plot(cflorida.nj7622WOoutgroup.ladder, "ph", cex = 1.05, label.offset = .013, x.lim = 1)
title(main = "NJ Tree + Structure + K means clustering", line= -.1) #play around with this more to put the titles in the right place, or maybe add later with a different program.
nodelabels(bpcflorida7622WO.1, adj = c(1.1, -0.2), frame = "none", cex = .9, col = "green4")
tiplabels(type = "p", pch = c(21, 22)[cflorida.pcogroupswild],
          bg = c(rgb(0, 0, 1, 0.3), rgb(1, 0, 0, 0.5), rgb(1, 0, 0, 0.05), rgb(0, 1, 0, 0.1), rgb(0, 0, 0, .3))[cflorida.pcogroups], 
          col = c(rgb(0, 0, 1, 1), rgb(1, 0, 0, 1), rgb(1, 0, 0, 1), rgb(0, 1, 0, 1), rgb(0, 0, 0, 1))[cflorida.pcogroups], offset= 0.007, cex = 1.5)
legend("bottomleft", legend = c("White-bracted", "Pink-bracted", "Blush pink-bracted", "ssp. urbiniana", "Unknown bract color", "", "Garden origin", "Wild-collected"), 
       bty='n', cex = 1, pch = c(21, 21, 21, 21, 21, 21, 21, 22), 
       pt.bg = c(rgb(0, 0, 1, 0.3), rgb(1, 0, 0, 0.5), rgb(1, 0, 0, 0.05), rgb(0, 1, 0, 0.1), rgb(0, 0, 0, .3), rgb(0,0,0,0), rgb(0, 0, 0, 1), rgb(0, 0, 0, 1)), 
       col = c(rgb(0, 0, 1, 1), rgb(1, 0, 0, 1), rgb(1, 0, 0, 1), rgb(0, 1, 0, 1), rgb(0, 0, 0, 1), rgb(0, 0, 0, 0), rgb(0, 0, 0, 1), rgb(0, 0, 0, 1)), pt.cex = 2, y.intersp = 1)
add.scale.bar(x=0.2,y=0.1, cex=1)
dev.off()
