# CanineGutMicrobiomeStudy
Companion scripts for the manuscript [Pellowe, S.D., Zhang, A., Bignell, D.R.D. et al. Gut microbiota composition is related to anxiety and aggression scores in companion dogs. Sci Rep 15, 24336 (2025)](https://doi.org/10.1038/s41598-025-06178-4)

The first step to reproduce our analysis is to run Dada2. To do this, before running the Dada2_script.R one needs to download the fastq files available at the [NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) under BioProject PRJNA1020865. Note that we performed quality trimmering/filtering of the reads before using Dada2.

For convenience we have uploaded the files obtained from the Dada2 Step in the [Data directory](https://github.com/BioinformaticsLabAtMUN/CanineGutMicrobiomeStudy/tree/main/Data). All the files needed to run the other analyses (Abundance Analysis, Compositional Balance Analysis with PhILR and Selbal, and Random Forest Analysis) are available in this Data directory. All these three R scripts (PhiLR_selbal.R, AbundanceAnalysis.R and RandomForest.R) can be run independently of each other. Note that there are comments in the scripts indicating whether something has to be changed depending on whether the analysis is for Aggression or Anxiety.

