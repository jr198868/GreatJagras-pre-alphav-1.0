### 16S-rRNA-amplicon-sequencing-characterization-of-biosolids-from-a-wastewater-treatment-plant-and-hum

library(dada2); 
packageVersion("dada2")

path<-"C:/Users/Administrator.000/Desktop/USA Stool Microbial 16S rRNA dataset/Fecal 16S rRNA analysis/SRR22" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#### #Forward and reverse fastq filenames have format: SAMPLENAME_R1_01.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

#### #Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:21])
plotQualityProfile(fnRs[1:21])

filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0,truncLen = c(240,140), maxEE=c(2,2), truncQ=2, trimLeft=c(35,21), rm.phix=TRUE, compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
#This is the filtering parameters of American intestinal microbiomes

#Filtering parameters of Chinese intestinal microbiomes
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0,truncLen = c(275,250), maxEE=c(2,2), truncQ=2, trimLeft=c(19,20), rm.phix=TRUE, compress=TRUE, multithread=FALSE) 
#### #On Windows set multithread=FALSE

####Filtering parameters of biosolid microbiomes
####out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0,truncLen = c(290,256), maxEE=c(2,2), truncQ=2, trimLeft=33, rm.phix=TRUE, compress=TRUE, multithread=FALSE) 
#### #On Windows set multithread=FALSE


plotQualityProfile(filtFs[1:3])
plotQualityProfile(filtRs[1:3])

head(out)

errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE, MAX_CONSIST = 13)

plotErrors(errF, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

#### #Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap=25, maxMismatch=0, verbose=TRUE)

#### #Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)

dim(seqtab)

table(nchar(getSequences(seqtab)))

#### #Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

#### #If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#### #Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "~/tax/silva_species_assignment_v128.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#### #Evaluate accuracy
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

