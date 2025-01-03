library(dada2)
path<-"~/Interact/seed_phytate"
list.files(path)
library(ShortRead)
library(Biostrings)
fnFs <- sort(list.files(path, pattern="1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="2.fastq.gz", full.names = TRUE))
sample.names <- lapply(lapply(strsplit(basename(fnFs), "_"), `[`, 2:6), function(x) paste(x,collapse = "_"))
plotQualityProfile(fnFs[6]) # looks okay until 280 nt
plotQualityProfile(fnRs[6]) #looks okay until 200 nt
table(nchar(getSequences(fnFs[6])))
#Searching the reads for primers (V5-V7)
FWD<-"AACMGGATTAGATACCCKG"
REV<-"ACGTCATCCCCACCTTCC"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
#create paths for filtered reads
filtFsNS <- file.path(path, "filtN", paste0(sample.names, "_F_filtN.fastq.gz"))
filtRsNS <- file.path(path, "filtN", paste0(sample.names, "_R_filtN.fastq.gz"))

names(filtFsNS) <- sample.names
names(filtRsNS) <- sample.names
#Trimming the reads to remove N's
filterAndTrim(fnFs, filtFsNS, fnRs, filtRsNS, maxN = 0, multithread = TRUE)

#Check length distribution
table(nchar(getSequences(filtFsNS[6])))
#Function for search of primers
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = filtFsNS[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = filtRsNS[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = filtFsNS[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = filtRsNS[[1]]))
#Removing primers
cutadapt<-"/usr/bin/cutadapt"
system2(cutadapt, args = "--version")
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

R1.flags <- paste("-g", FWD, "-a", REV.RC) 
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             filtFsNS[i], filtRsNS[i])) # input files
} #It gives a warning about incomplete adapter sequence but that sequences is okay
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
cutFs <- sort(list.files(path.cut, pattern = "1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "2.fastq.gz", full.names = TRUE))

plotQualityProfile(cutRs[2])
table(nchar(getSequences(cutFs[2]))) # The majority of the reads are 382 bp - 
table(nchar(getSequences(cutRs[2]))) #The majority of the reads are 280 bp 
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

#The overlap should be (793+282)-(1193-280)=168 bp
#Hence, if I truncate 10 bp from  forward and 70 from reverse reads there should still be lots of overlap. 
#This truncationis suggested by the authors of dada2
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 3), 
                     truncQ = 2, minLen = 200, truncLen = c(272,220), rm.phix = TRUE, compress = TRUE, multithread = TRUE)  
#See dada2_few_samples.R for more information
out
save.image("first_part.RData")
errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ = TRUE) #error rates look okay
errR <- learnErrors(filtRs, multithread=TRUE) 
plotErrors(errR, nominalQ =TRUE)

dadaFs <-dada(filtFs, err = errF, multithread = 20, pool = "pseudo")
dadaRs<-dada(filtRs, err = errR, multithread = 20, pool = "pseudo")
save.image("second.part.RData")
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab,  multithread=TRUE, verbose=TRUE, method = "pooled")
sum(seqtab.nochim)/sum(seqtab) 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim)/out[,1])
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim","ratio")
track

save.image(file = "syncom_all.RData")
#Assign taxonomy
taxa<-assignTaxonomy(seqtab.nochim,tryRC = TRUE,"~/Interact/Syncom_amp/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread = 40) 
saveRDS(taxa,"taxa.rds")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#Convert to a phyloseq object (easy to work with)
samples.out <- rownames(seqtab.nochim)
library(phyloseq)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))
ps
dna <- Biostrings::DNAStringSet(taxa_names(ps)) #Here the DNA sequences are stored in a separate object, which is then put into the phyloseq object. It might be usefull later on in an anlysis.
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps))) #Change taxa names to ASV1, ASV2, etc.
ps
sample_names(ps)<-sample.names

saveRDS(ps, "Syncom_phyloseq.rds") #Save phyloseq object for later use.

