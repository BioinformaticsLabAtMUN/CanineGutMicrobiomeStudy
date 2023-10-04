

## First construct phylogenetic tree as per https://f1000research.com/articles/5-1492

## Create phyloseq object including tree using https://vaulot.github.io/tutorials/Phyloseq_tutorial.html

## Continue with Philr analysis using https://rdrr.io/bioc/philr/f/vignettes/philr-intro.Rmd

## Intro to Philr: https://bioconductor.org/packages/release/bioc/vignettes/philr/inst/doc/philr-intro.html


.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")

sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
set.seed(100)

##### CHANGE PATH TO FILES AS APPROPRIATE
otuda <- read.table("Data/SeqTab_NoChim_SamplesInColumns.tsv", header = TRUE, sep="\t")
otuda <- data.frame(t(otuda), check.names = F)
sampleda <- read.csv("Data/samplevariables.csv", header=TRUE, sep=",", row.names=1, comment.char = "")
taxda <- read.table("Data/taxa_dada2.tsv", header=TRUE, sep="\t", row.names=1, check.names=F, comment.char="")
taxda <- taxda[match(colnames(otuda), rownames(taxda)),]
otuda <- as.matrix(otuda)
taxda <- as.matrix(taxda)
OTU = otu_table(otuda, taxa_are_rows = TRUE)
TAX = tax_table(taxda)
samples = sample_data(sampleda)

cbarq <- read.table("Data/CBARQ_ContinuousData.csv", header = T, sep = ",")
dim(cbarq)
cbarq[,1] <- paste("Sample", cbarq[,1],sep ="")
rownames(cbarq) <- cbarq[,1]

## construct phylogenetic tree

seqs <- getSequences(otuda)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)


phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) 
fit = pml(treeNJ, data=phang.align)


fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

## create phyloseq object

ps <- phyloseq(tax_table(taxda), sample_data(sampleda), otu_table(otuda, taxa_are_rows = FALSE), phy_tree(fitGTR$tree))

## continue PhILR analysis

library(philr); packageVersion("philr")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(ape); packageVersion("ape")


GP <-  filter_taxa(ps, function(x) sum(x > 3) > (0.1*length(x)), TRUE)
GP <-  filter_taxa(GP, function(x) sd(x)/mean(x) > 3.0, TRUE)
GP <- transform_sample_counts(GP, function(x) x+1)

is.rooted(phy_tree(GP)) # Is the tree Rooted? 
is.binary(phy_tree(GP)) # All multichotomies resolved?

phy_tree(GP) <- makeNodeLabel(phy_tree(GP), method="number", prefix='n')

name.balance(phy_tree(GP), tax_table(GP), 'n1')

## Check  variables.

otu.table <- otu_table(GP) #Do not transpose as it is done in the tutorial
tree <- phy_tree(GP)
metadata <- sample_data(GP)
tax <- tax_table(GP)

dim(otu.table)
otu.table[1:2,1:3] # OTU Table
tree # Phylogenetic Tree
head(metadata,2) # Metadata
head(tax,2) # taxonomy table


sbp <- phylo2sbp(tree)
dim(sbp)

gp.philr <- philr(otu.table, tree, sbp,
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')
gp.philr[1:5,1:5]

#gp.philr is the transformed data
write.table(gp.philr, file = "ILR.tsv", sep = "\t")


gp.dist <- dist(gp.philr, method="euclidean")
gp.pcoa <- ordinate(GP, 'PCoA', distance=gp.dist)
pdf(file ="PCOA.pdf")
plot_ordination(GP, gp.pcoa, color='Aggression') + geom_point(size=4)
plot_ordination(GP, gp.pcoa, color='Anxiety') + geom_point(size=4)
dev.off()


## fit a sparse logistic regression model (logistic regression with $l_1$ penalty)
### CHANGE var2Plot to do the analysis for Anxiety or Aggression
var2Plot <- "Aggression"

library(glmnet); packageVersion('glmnet')
glmmod <- glmnet(gp.philr, sample_data(GP)[[var2Plot]], alpha=1, family="binomial")

top.coords <- as.matrix(coefficients(glmmod, s=0.14)) #top5
top.coords <- rownames(top.coords)[which(top.coords != 0)]
(top.coords <- top.coords[2:length(top.coords)])

tc.names <- sapply(top.coords, function(x) name.balance(tree, tax, x))
tc.names

## Visualize results - boxplots
gp.philr.long <- convert_to_long(gp.philr, sample_data(GP)[[var2Plot]]) 
gp.philr.long <- gp.philr.long[gp.philr.long$coord %in% top.coords, ] 

pdf(paste("Top5_PhILR_Balances_", var2Plot, ".pdf", sep = ""))
ggplot(gp.philr.long, aes(x=labels, y=value)) +
  geom_boxplot(fill='lightgrey') +
  facet_grid(.~coord, scales='free_x') +
  xlab(var2Plot) + ylab('Balance Value') +
  theme_bw()
dev.off()  
  
## scatterplots
library(tidyr); packageVersion('tidyr')

#### NOTE: NEED TO CHANGE VARIABLE NAME IN LINE DPLYR::RENAME
createScatterplot <- function(i, j, var2Plot){
		  p <- gp.philr.long %>%
  		dplyr::rename(Aggression = labels) %>%
  		dplyr::filter(coord %in% c(i, j)) %>%
  		tidyr::spread(coord, value) %>%
  		ggplot(aes_string(x=i, y=j, color=var2Plot)) +
  		geom_point(size=4) +
  		labs(x = tc.names[i], y = tc.names[j]) +
  		theme_bw()
  		p
}

pdf(paste("Scatterplots", var2Plot, ".pdf", sep = ""))
for (k in 1:length(top.coords)){
	for (i in setdiff(1:length(top.coords),k)) {
 		print(createScatterplot(top.coords[k],  top.coords[i], var2Plot))
	}
}
dev.off()
#save.image("PhILR_analysis.RDATA")
#### Selbal
#### Software: https://github.com/malucalle/selbal
#### Tutorial: https://htmlpreview.github.io/?https://github.com/malucalle/selbal/blob/master/vignettes/vignette.html
#### Article: https://journals.asm.org/doi/epub/10.1128/mSystems.00053-18

# Create a count matrix by taxa - add read counts of all OTUs identified as belonging to the same taxa - NA is removed
getCountSums <- function(s, colname){
	print(paste(s, colname, sep =" "))
	seqs <- grepl(s, taxda[, colname], fixed = T)
	if (sum(seqs) > 1){
		rowSums(otuda[,row.names(taxda)[seqs]], na.rm = T)
	} else {
		otuda[,row.names(taxda)[seqs]]
	}
}

foo <- sapply(colnames(taxda)[-1], function(i){
	uniqueSp <- unique(taxda[,i])
	uniqueSp <- uniqueSp[!is.na(uniqueSp)]
	print(length(uniqueSp))
	m <- sapply(uniqueSp, getCountSums, i)
	colnames(m) <- paste(substr(i,1,2), colnames(m), sep ="_")
	m
})

nCols <- sum(unlist(sapply(foo, ncol)))

CountMatrix <- matrix(unlist(foo), ncol = nCols, dimnames= list(row.names(foo[[1]]), unlist(sapply(foo, colnames)) ))

dim(CountMatrix)
head(CountMatrix)

#remove columns with more than 80% zeros to avoid sebel warning
CountMatrix <- CountMatrix[,apply(CountMatrix, 2, function(c) { sum(c > 0) > floor(nrow(CountMatrix) *.2)})]

anxietyM <- cbind(CountMatrix, cbarq[row.names(CountMatrix), "Anxiety"])

aggressionM <- cbind(CountMatrix, cbarq[row.names(CountMatrix), "Aggression"])

runSelbal <- function(x2, y2, nfold = 10, niter=10, seed = 1718){
	require(selbal)
	bal <- selbal.cv(x = x2, y = y2, n.fold = nfold, n.iter = niter, covar = NULL, seed = seed)
     print(bal$opt.nvar) #Optimal number of variables to use
     print(summary(bal$cv.accuracy)) #MSE or AUC
	print(bal$glm) #GLM
	return(bal)
}

# Run selbal.cv function for Regression
BAL.anxiety <- runSelbal(anxietyM[,1:ncol(CountMatrix)], anxietyM[,ncol(CountMatrix)+1]) 

BAL.aggression <- runSelbal(aggressionM[,1:ncol(CountMatrix)], aggressionM[,ncol(CountMatrix)+1]) 

#### SDA, DDA, NSF
BAL <- list()
for (i in c("SDA", "ODA", "DDF", "SDF", "NSF", "Sep")){
	m <- cbind(CountMatrix, cbarq[row.names(CountMatrix), i])
    BAL[[i]]<- runSelbal(m[,1:ncol(CountMatrix)], m[,ncol(CountMatrix)+1], niter = 5) 
}

for (i in c("SDA", "ODA", "DDF", "SDF", "NSF", "Sep")){
	pdf(paste("Selbal_plots_class_", i, ".pdf", sep = ""))
	print(BAL[[i]]$accuracy.nvar)
	print(BAL[[i]]$var.barplot)
	plot.new()
	grid.draw(BAL[[i]]$global.plot)
	plot.new()
	plot.tab(BAL[[i]]$cv.tab)
	dev.off()
}

# Run selbal.cv function for Classification
m <- as.data.frame(CountMatrix)
m$Class <- as.factor(sampleda[row.names(CountMatrix), "Anxiety"])
BAL.anxiety_clas <- runSelbal(m[,1:ncol(CountMatrix)], m[,ncol(CountMatrix)+1])

m$Class <- as.factor(sampleda[row.names(CountMatrix), "Aggression"])
BAL.aggression_clas <- runSelbal(m[,1:ncol(CountMatrix)], m[,ncol(CountMatrix)+1]) 


### CHANGE to plot for Anxiety.
## Plot selbal for classification
pdf("Selbal_plots_class_Aggresion.pdf")
	BAL.aggression_clas$ROC.plot #only available for classification
	BAL.aggression_clas$accuracy.nvar
	 BAL.aggression_clas$var.barplot
	 plot.new()
	 grid.draw(BAL.aggression_clas$global.plot)
	plot.new()
	plot.tab(BAL.aggression_clas$cv.tab)
dev.off()


###  BAL.aggression_clas$global.plot$grobs has the components from the plot. 1 is the labels, 2 is the ROC, 3 boxplot, 4 density plot.
# grid.draw(BAL.anxiety_clas$global.plot$grobs[[1]])
library(ggplotify)
library(gridExtra)
gt <- arrangeGrob(BAL.anxiety_clas$global.plot$grobs[[3]], BAL.anxiety_clas$global.plot$grobs[[2]], ncol = 2, nrow = 1)
p <- as.ggplot(gt) + draw_plot_label(label = "Numerator: Oscillospiraceae, Negativicutes\nDenominator: Blautia", size = 9, x = 0.06, y = 1, hjust = 0, vjust = 1.4, font = "italic") + draw_plot_label(label = c("A","B"), size = 11, x=c(0, 0.5), y=c(1,1) )
p









