# **Plant-LncPip**





### **A pipeline for identifying and characterizing lncRNAs in plants.**

#### **Email:** xuechan.tian@bjfu.edu.cn








## **1. Dependencies.** 



    HISAT2
    StringTie
    LncFinder
    CPAT
    Diamond





## **2.Installation**



### **2.1. Install HISAT2:**
    git clone https://github.com/DaehwanKimLab/hisat2.git
    cd hisat2
    make



### **2.2. Install StringTie:**
    wget https://github.com/gpertea/stringtie/releases/download/v2.1.4/stringtie-2.1.4.Linux_x86_64.tar.gz
    tar xzf stringtie-2.1.4.Linux_x86_64.tar.gz

	
	
### **2.3. Install LncFinder-plant:**
 R Package
 
    install.packages("LncFinder")
    install.package("seqinr")



### **2.4. Install CPAT-plant:**
 CPAT (version 1.2.4)
 
    conda create -n py27 python=2.7 -y
    source activate py27
    pip2 install CPAT


	
### **2.5. Install Diamond:**
     conda install -c bioconda diamond
	
	
	
	
## **3. Run HISAT2 to map RNA-seq reads to the reference genome:**



### **3.1. Construct reference genome**

    hisat2-build -p 8 genome.fasta genome.index 


### **3.2. Genome alignment with hisat2**
If the RNA-seq library is strand-specific

    for i in `cat sample.list`; do hisat2 --new-summary --rna-strandness RF -p 10 -x genome.index -1 ${i}_1.fastq -2 ${i}_2.fastq -S ${i}sam; done

	
	
If the RNA-seq library is not strand-specific

    for i in `cat sample.list`; do hisat2 --new-summary -p 10 -x genome.index -1 ${i}_1.fastq -2 ${i}_2.fastq -S ${i}.sam; done


### **3.3. Sort and compress sam files with samtools**
    for i in `cat sample.list`; do samtools sort -o ${i}.bam ${i}.sam; done
	
	

## **4. Assemble transcripts using StringTie:**	
	
	
	
### **4.1. Transcription reconstruction of a single sample**
    for i in `cat sample.list`; do stringtie -p 10 --rf -G genome.gtf -o ${i}.gtf ${i}.bam; done

	
	
	
### **4.2. Merge transcripts of multiple samples**	
    stringtie --merge -o candidate_transcript.gtf -G genome.gtf SRR*.gtf
    gffread -w candidate_transcript.fasta -g genome.fasta candidate_transcript.gtf
    grep '>' candidate_transcript.fasta | awk '{print $1}' | sed 's/>//g' | sort -u > candidate_transcript.txt
	
	
	
	
	
## **5. LncRNA identification:**	



### **5.1. Get the sequence id of length greater than 200bp.**

    awk 'BEGIN{RS=">";FS="\n"}NR>1{seq="";for (i=2;i<=NF;i++) seq=seq""$i; print ">"$1"\n"seq}' candidate_transcript.fasta | awk '{if($0 ~ /^>/) {id=$0;} else {len=length($0); if(len > 200) {print id}}}' | cut -d " " -f 1 |sed 's/>//g' >  filtered_transcript.txt



### **5.2.  Identification of lncRNA with LncFinder-plant.**	

R Package

    R
    library(LncFinder)
    library(seqinr)
	
	
import training data

    mRNA <- seqinr::read.fasta(file ="./data/training/mRNA.fasta")
    lncRNA <- seqinr::read.fasta(file ="./data/training/lncRNA.fasta")
    
    
Use "make_frequencies" function to generate the feature file.

    frequencies <- make_frequencies(cds.seq = mRNA, lncRNA.seq = lncRNA, SS.features = FALSE, cds.format = "DNA", lnc.format = "DNA", check.cds = TRUE, ignore.illegal = TRUE)	
	
	
loading the model

    plant = readRDS("Plant_model.rda")
	
	
Identification of lncRNA 

    Seqs <- seqinr::read.fasta(file ="candidate_transcript.fasta")
    Plant_results <- LncFinder::lnc_finder(Seqs, SS.features = FALSE, format = "DNA", frequencies.file = frequencies, svm.model = plant, parallel.cores = 2)
	
	
Export results

    write.table(Plant_results, file ="plant-lncFinder.txt", sep ="\t",row.names =TRUE, col.names =TRUE,quote =FALSE)




### **5.3. Identification of lncRNA with CPAT-plant.**	
The coding probability (CP) cutoff: 0.46 (CP >=0.46 indicates coding sequence, CP < 0.46 indicates noncoding sequence).

    cpat.py -x Plant_Hexamer.tsv -d Plant.logit.RData -g candidate_transcript.fasta -o CPAT_plant.output





### **5.4. Alignment of sequences to the UniProt protein database with diamond.**
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    gunzip uniprot_sprot.fasta.gz
    diamond makedb --in uniprot_sprot.fasta -d uniprot_out
    diamond blastx -d uniprot_out -q candidate_transcript.fasta -o uniprotoutput.txt




### **5.5. By intersecting the results obtained from the aforementioned steps, a set of high-confidence lncRNAs were obtained.**
	
	

	
	
## **6. Classify the final set of lncRNAs based on their genomic locations and sequence features.**		
	FEELnc_classifier.pl -i lncRNA.gtf -a mRNA.gtf > lncRNA_classes.txt
	
	
	
	
## **7. TE-derived lncRNAs.**		
	bedtools intersect -a lncRNA.bed -b TE.bed -wo | sort -u > TE_lncRNA_intersect.txt 
	
	
	
	
	
	
	
	
