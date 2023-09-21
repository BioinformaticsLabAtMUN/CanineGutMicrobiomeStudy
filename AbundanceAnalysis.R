######
##### Code based on the code provided by Shuangbin Xu and  Guangchuang Yu in the Workshop of microbiome dataset analysis using MicrobiotaProcess
##### available at https://yulab-smu.top/MicrobiotaProcessWorkshop/articles/MicrobiotaProcessWorkshop.html and
##### https://github.com/YuLab-SMU/MicrobiotaProcessWorkshop/blob/master/vignettes/MicrobiotaProcessWorkshop.Rmd
#####

suppressPackageStartupMessages({
  library(MicrobiotaProcess) # an R package for analysis, visualization and biomarker discovery of Microbiome.
  library(phyloseq) # Handling and analysis of high-throughput microbiome census data.
  library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics.
  library(tidyverse) # Easily Install and Load the 'Tidyverse'.
  library(vegan) # Community Ecology Package.
  library(coin) # Conditional Inference Procedures in a Permutation Test Framework.
  library(reshape2) # Flexibly Reshape Data: A Reboot of the Reshape Package.
  library(ggnewscale) # Multiple Fill and Colour Scales in 'ggplot2'.
})


otuda <- read.table("SeqTab_NoChim_SamplesInColumns.tsv", header = TRUE, sep="\t")
otuda <- data.frame(t(otuda), check.names = F)
sampleda <- read.csv("samplevariables.csv", header=TRUE, sep=",", row.names=1, comment.char = "")
taxda <- read.table("taxa_dada2.tsv", header=TRUE, sep="\t", row.names=1, check.names=F, comment.char="")
taxda <- taxda[match(colnames(otuda), rownames(taxda)),]
psraw <- import_dada2(seqtab=otuda, taxatab=taxda, sampleda=sampleda)
#### Not of our samples was removed
psraw <- prune_samples(sample_sums(psraw)>=sort(rowSums(otu_table(psraw)))[3], psraw)
set.seed(1024)
ps <- rarefy_even_depth(psraw)
ps

# We also can use phyloseq to build phyloseq object.
#library(phyloseq)
#ps2 <- phyloseq(otu_table(otuda, taxa_are_rows=FALSE), sample_data(sampleda), tax_table(as.matrix(taxda)))
#ps2

######
colnames(sampleda)
##### MODIFY VAR2PLOT TO GET THE ANALYSIS EITHER FOR "Anxiety" OR "Aggression"
Var2Plot <- "Anxiety"

### Alpha diversity rarefication curves ### 
library(patchwork)
set.seed(1024)
rareres <- get_rarecurve(obj=ps, chunks=400)

p_rareanx <- ggrarecurve(obj=rareres,
                      indexNames=c("Observe","Chao1","ACE"),
) +
  theme(legend.spacing.y=unit(0.01,"cm"),
        legend.text=element_text(size=4))
     
prare1anx <- ggrarecurve(obj=rareres, factorNames=Var2Plot,
                      indexNames=c("Observe", "Chao1", "ACE")
) +
  scale_fill_manual(values=c("#00AED7", "#FD9347"))+
  scale_color_manual(values=c("#00AED7", "#FD9347"))+
  theme_bw()+
  theme(axis.text=element_text(size=8), panel.grid=element_blank(),
        strip.background = element_rect(colour=NA,fill="grey"),
        strip.text.x = element_text(face="bold"))
       
prare2anx <- ggrarecurve(obj=rareres,
                      factorNames=Var2Plot,
                      shadow=FALSE,
                      indexNames=c("Observe", "Chao1", "ACE")
) +
  scale_color_manual(values=c("#00AED7", "#FD9347"))+
  theme_bw()+
  theme(axis.text=element_text(size=8), panel.grid=element_blank(),
        strip.background = element_rect(colour=NA,fill="grey"),
        strip.text.x = element_text(face="bold"))
        

## Alpha Index ##
alphaobj <- get_alphaindex(ps)
head(as.data.frame(alphaobj))
p_alphaanx <- ggbox(alphaobj, geom="violin", factorNames=Var2Plot) +
  scale_fill_manual(values=c("#00AED7", "#FD9347"))+
  theme(strip.background = element_rect(colour=NA, fill="grey"))


################# New code added June 16, 2023 by LPC
### See https://github.com/joey711/phyloseq/issues/1521
### To get the relative abundance per OTU into a data frame and write the data frame into a file
library(dplyr)
library(tidyr)

write.table(ps %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% 
	arrange(OTU) %>% rename(ASV = OTU) %>% 
     select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species, Sample, Abundance) %>%
     spread(Sample, Abundance), 
	file = "ps.relative_abundance.all.tsv", sep = "\t", quote = F, row.names = F, col.names = T)


### To get the relative abundance for a certain level into a data frame and write the data frame into a file
level = "Genus" #Change level to Phylum, Class, Order, Family, Genus as desired
write.table(ps %>% tax_glom(taxrank = level) %>% 
        transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>%
        select(Genus, Sample, Abundance) %>% spread(Sample, Abundance), #CHANGE THE LEVEL IN THIS LINE
	file = paste("ps.relative_abundance.", level,".tsv", sep =""), sep = "\t", quote = F, row.names = F, col.names = T)

### To get the relative abundance summarized per group
mergedAbundances <- merge_samples(ps, Var2Plot)

# Use psmelt to obtain a data.frame
mergedAbundances_DF <- mergedAbundances %>% tax_glom(taxrank = level) %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()
head(mergedAbundances_DF)


###################################
#### Taxonomy composition analysis
### CHOOSE THE LEVEL FOR THE ANALYSIS
classtaxa <- get_taxadf(obj=ps, taxlevel=3) #level 1 is kingdom, 5 family, 6 genus

### To get the raw abundance into a data frame and then one can write the data frame in a file
abundanceDataTable <- psmelt(classtaxa)
dim(abundanceDataTable)
head(abundanceDataTable)


# The 30 most abundant taxonomy will be visualized by default (parameter `topn=30`). 
pclass <- ggbartax(obj=classtaxa, facetNames=Var2Plot ) +
          xlab(NULL) +
          ylab("relative abundance (%)") +
          scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
          guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))


# Show the abundance in different groups.
fclass <- ggbartax(obj=classtaxa, facetNames=Var2Plot , plotgroup=TRUE, topn=10) +
          xlab(NULL) +
          ylab("relative abundance (%)") +
          scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
          guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5, ncol=2))


pclass2 <- ggbartax(obj=classtaxa, count=TRUE, facetNames=Var2Plot ) +
          xlab(NULL) +
          ylab("count reads") +
          scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
          guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))


##### Step 4.2 of the workshop
### Venn diagram
vennlist <- get_vennlist(obj=ps, factorNames=Var2Plot)
vennlist <- list(High = vennlist[,1], Low = vennlist[,2])
library(VennDiagram)
vennp <- venn.diagram(vennlist,
                      height=5,
                      width=5, 
                      filename=NULL, 
                      fill=c("#00AED7", "#FD9347"),
                      cat.col=c("#00AED7", "#FD9347"),
                      alpha = 0.85, 
                      fontfamily = "serif",
                      fontface = "bold",
                      cex = 1.2,
                      cat.cex = 1.3,
                      cat.default.pos = "outer",
                      cat.dist=0.1,
                      margin = 0.1, 
                      lwd = 3,
                      lty ='dotted',
                      imagetype = "svg")
                      
### All OTUs are in both groups!
grid::grid.draw(vennp)

### Beta Analysis

#Step 5.1 PCA Analysis
# If the input was normalized, the method parameter should be setted NULL.
pcares <- get_pca(obj=ps, method="hellinger")
# Visulizing the result
pcaplot1 <- ggordpoint(obj=pcares, biplot=TRUE, speciesannot=TRUE,
                      factorNames=c(Var2Plot), ellipse=TRUE) +
            scale_color_manual(values=c("#00AED7", "#FD9347")) +
            scale_fill_manual(values=c("#00AED7", "#FD9347"))
# pc = c(1, 3) to show the first and third principal components.
pcaplot2 <- ggordpoint(obj=pcares, pc=c(1, 3), biplot=TRUE, speciesannot=TRUE,
                       factorNames=c(Var2Plot), ellipse=TRUE) +
            scale_color_manual(values=c("#00AED7", "#FD9347")) +
            scale_fill_manual(values=c("#00AED7", "#FD9347"))

#Step 5.2 PCoA analysis

# distmethod
# "unifrac",  "wunifrac", "manhattan", "euclidean", "canberra", "bray", "kulczynski" ...(vegdist, dist)
pcoares <- get_pcoa(obj=ps, distmethod="bray", method="hellinger")
# Visualizing the result
pcoaplot1 <- ggordpoint(obj=pcoares, biplot=TRUE, speciesannot=TRUE,
                       factorNames=c(Var2Plot), ellipse=TRUE) +
            scale_color_manual(values=c("#00AED7", "#FD9347")) +
            scale_fill_manual(values=c("#00AED7", "#FD9347"))
# first and third principal co-ordinates
pcoaplot2 <- ggordpoint(obj=pcoares, pc=c(1, 3), biplot=TRUE, speciesannot=TRUE,
                        factorNames=c(Var2Plot), ellipse=TRUE) +
             scale_color_manual(values=c("#00AED7", "#FD9347")) +
             scale_fill_manual(values=c("#00AED7", "#FD9347"))

# Step 5.3 Permutational Multivariate Analysis of Variance
distme <- get_dist(ps, distmethod ="bray", method="hellinger")
sampleda <- data.frame(sample_data(ps), check.names=FALSE)
sampleda <- sampleda[match(colnames(as.matrix(distme)),rownames(sampleda)),,drop=FALSE]
sampleda$Group <- factor(sampleda[,Var2Plot])
set.seed(1024)
adores <- adonis(distme ~ Group, data=sampleda, permutation=9999)
data.frame(adores$aov.tab)



## Step 6 Biomarker discovery ##

library(coin)
set.seed(1024)
deres <- diff_analysis(obj = ps, classgroup = Var2Plot,
                       mlfun = "lda",
                       filtermod = "pvalue",
                       firstcomfun = "kruskal.test",
                       firstalpha = 0.1,
                       strictmod = TRUE,
                       secondcomfun = "wilcox.test",
                       subclmin = 3,
                       subclwilc = TRUE,
                       secondalpha = 0.1,
                       ldascore=3)
deres

## Anxiety
# The original data: 52 features and 48 samples
# The sample data: 1 variables and 48 samples
# The taxda contained 539 by 7 rank
# after first test (kruskal.test) number of feature (pvalue<=0.1):9
# after second test (wilcox.test and generalizedFC) number of significantly discriminative feature:9
# after lda, Number of discriminative features: 9 (certain taxonomy classification:6; uncertain taxonomy classication: 3)

### Visualization of differential analysis results

diffbox <- ggdiffbox(obj=deres, box_notch=FALSE, 
             colorlist=c("#00AED7", "#FD9347"), l_xlabtext="relative abundance")


ggdifftaxbar(obj=deres, xtextsize=1.5, 
             output=paste(Var2Plot,"_biomarker_barplot", sep =""),
             coloslist=c("#00AED7", "#FD9347"))

#### Save all the plots to a PDF
pdf(paste(Var2Plot, "_abundanceAnalysis.pdf", sep =""))
p_rareanx / prare1anx / prare2anx
p_alphaanx
pclass
fclass
pclass2
pcaplot1 | pcaplot2
pcoaplot1 | pcoaplot2
diffbox
dev.off()
