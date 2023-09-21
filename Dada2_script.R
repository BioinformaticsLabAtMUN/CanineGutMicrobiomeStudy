##########################################
## This code is based on code provided by
## Richa Bharti, Dominik G Grimm, Current challenges and best-practice protocols for microbiome analysis, Briefings in Bioinformatics, Volume 22, Issue 1, January 2021, Pages 178–193, https://doi.org/10.1093/bib/bbz155 at 
## https://github.com/grimmlab/MicrobiomeBestPracticeReview/blob/master/Amplicon_analysis/dada2_workflow.R
## and by Benjamin Callahan in DADA2 Pipeline Tutorial at  https://benjjneb.github.io/dada2/tutorial_1_8.html
############################################

library(dada2); packageVersion("dada2")
sessionInfo()

### SET THE PATH TO THE FASTQ FILES
path <- "~/projects/def-lpenacas/CanineGut/Fastq_files/2_trimmed_fastq_files"
### SET THE NAME OF THE OUTPUT DIRECTORY
outdir <- "output_dada2_taxonomy"

### MAKE SURE THE EXTENSION OF THE FILES IS AS GIVEN TO getFiles function and OR SET AS APPROPRIATE
getFiles <- function(pattern) {
    sort(list.files(path, pattern, full.names = TRUE))
}
fnFs <- getFiles("_R1_paired_001.fastq.gz")
fnRs <- getFiles("_R2_paired_001.fastq.gz")


sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
png("QualityProfiles_forward.png")
plotQualityProfile(fnFs[1:2])
dev.off()
png("QualityProfiles_reverse.png")
plotQualityProfile(fnRs[1:2])
dev.off()
filtFs <- file.path(path, "dada2filter", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "dada2filter", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Filter and trim
### SET truncLen as appropriate for data
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,190), maxN=0, truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
head(out)
write.csv(out, file = paste0(outdir,"/out_dada2_filterAndTrim.csv"))

png("QualityProfiles_filterAndTrim_forward.png")
plotQualityProfile(filtFs[1:2])
dev.off()

png("QualityProfiles_filterAndTrim_reverse.png")
plotQualityProfile(filtRs[1:2)
dev.off()

#Learn the error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
png("dada2_errors.png")
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

################
# Dereplication
################
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


###################
# Sample Inference
###################
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


####################
# Merge paired reads
####################
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

#Print mergers for all samples
sapply(1:length(mergers), function(s){
	filename <- paste(name(mergers)[s], "mergers_dada2.csv", sep="_")
	write.csv(mergers[[s]], file = paste(outdir, filename, sep ="/"))
}

##########################
# Construct sequence table
###########################
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
write.table(seqtab, file = paste0(outdir, "/seqtab_dada2.tab"), sep = "\t", col.names = TRUE, row.names = TRUE)

##################
# Remove chimeras
##################
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab)))
sum(seqtab.nochim)/sum(seqtab)

#a data.frame of merged, error-free, non-chimeric, amplicon sequence variants
write.table(seqtab.nochim, file = paste0(outdir, "/seqtab_nochim_dada2.tab"), sep = "\t", col.names = TRUE, row.names = TRUE)


#################
#Save into Rds
#################
saveRDS(seqtab, paste0(outdir, "/seqtab.Rds"))
saveRDS(seqtab.nochim, paste0(outdir, "/seqtab.nochim.Rds"))


#################################
# rack reads through the pipeline
#################################
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, file = paste(outdir, "trackReads_dada2.csv", sep = "/"))

rm(list= ls())

seqtab.nochim <- readRDS(paste0(outdir, "/seqtab.nochim.Rds"))

##############
# Assign taxonomy
#################
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa", minBoot = 80, outputBootstraps = TRUE, multithread=TRUE)
taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa", allowMultiple = TRUE)
write.table(taxa, file = paste0(outdir,"/taxa_dada2.tsv"), quote=FALSE, sep="\t")
