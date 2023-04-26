#!Rscript
# Usage: Rscript myscript.R
##################################################################################################################
#Taking the intersection of lncRNA identification results, to obtain high confidence lncRNAs.
##################################################################################################################
	
args<-commandArgs(TRUE)

#################################################
## R package
myPaths <- .libPaths()
library(tidyverse)

#################################################
#setwd(args[1])

#################################################
##### 
candidate_lncRNA <- read.table(args[1], header = FALSE, col.names = "genes") %>% pull(genes)
#candidate_lncRNA <- read.table("candidate_transcript.txt",  header = FALSE, col.names = "genes") %>% pull(genes)


##### shorter than 200 bp and overlapping with known mRNAs.
length <- read_delim(args[2], delim = "\t", col_names = FALSE) %>% pull(X1)
#length <- read_delim("candidate_lncRNA.txt", delim = "\t", col_names = FALSE) %>% pull(X1)

 
 
##### CPAT-plant
CPAT <- read_delim(args[3], delim = "\t") %>% filter(coding_prob < 0.46) %>% pull(mRNA_size)
#CPAT <- read_delim("CPAT_plant.output", delim = "\t") %>% filter(coding_prob < 0.46) %>% pull(mRNA_size) 



##### LncFinder-plant
LncFinder <- read_delim(args[4], delim = "\t") %>% filter(Coding.Potential == "NonCoding") %>% pull(Pred)
#LncFinder <- read_delim("plant-lncFinder.txt", delim = "\t") %>% filter(Coding.Potential == "NonCoding") %>% pull(Pred)




##### Swissprot
uniprot <- read_delim(args[5], delim = "\t", col_names = FALSE) %>% filter(X3>80, X11<1e-5) %>% pull(X1) %>% unique()
uniprot <- dplyr::setdiff(candidate_lncRNA, uniprot)

#uniprot <- read_delim("uniprotoutput.txt", delim = "\t", col_names = FALSE) %>% filter(X3>80, X11<1e-5) %>% pull(X1) %>% unique()
#uniprot <- dplyr::setdiff(candidate_lncRNA, uniprot)


##### Take intersection
lncRNA = Reduce(intersect, list(length, CPAT, LncFinder, uniprot))
write.table(lncRNA, file = "lncRNA_results.txt",
            row.names = F,
            col.names = F,
            quote = F)














