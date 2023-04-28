![LncRNA pipline](https://github.com/xuechantian/LncRNA-pipline/blob/master/Plant-LncPip_workflow.png)





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
     
     
### **2.6. Install FEELnc:**
    git clone https://github.com/tderrien/FEELnc.git
    export FEELNCPATH=/path/to/FEELnc/bin/
    export PERL5LIB=$PERL5LIB:/path/to/FEELnc/lib/
    export PATH=$PATH:/path/to/FEELnc/scripts/

	
	
## **3. Run HISAT2 to map RNA-seq reads to the reference genome:**



### **3.1. Construct reference genome**

    hisat2-build -p 8 genome.fasta genome.index 


### **3.2. Genome alignment with hisat2**

cat sample.list

    SRR3087771
    SRR3087772
    SRR3087773
    SRR3087774
    ...

If the RNA-seq library is strand-specific

    for i in `cat sample.list`; do hisat2 --new-summary --rna-strandness RF -p 10 -x genome.index -1 ${i}_1.fastq -2 ${i}_2.fastq -S ${i}.sam; done

	
	
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



### **5.1. Remove transcripts shorter than 200 bp and overlapping with known mRNAs.**

    FEELnc_filter.pl -i candidate_transcript.gtf -a genome.gtf --monoex=-1 -s 200 -p 20 > candidate_lncRNA.gtf
    cut -d ";" -f 2 candidate_lncRNA.gtf |sed 's/ transcript_id //g' | sed 's/"//g' | sort -u > candidate_lncRNA.txt


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

    plant = readRDS("./Model/Plant_model.rda")
	
	
Identification of lncRNA 

    Seqs <- seqinr::read.fasta(file ="candidate_transcript.fasta")
    Plant_results <- LncFinder::lnc_finder(Seqs, SS.features = FALSE, format = "DNA", frequencies.file = frequencies, svm.model = plant, parallel.cores = 2)
	
	
Export results

    write.table(Plant_results, file ="plant-lncFinder.txt", sep ="\t",row.names =TRUE, col.names =TRUE,quote =FALSE)




### **5.3. Identification of lncRNA with CPAT-plant.**	
The coding probability (CP) cutoff: 0.46 (CP >=0.46 indicates coding sequence, CP < 0.46 indicates noncoding sequence).

    source activate py27
    cpat.py -x ./Model/Plant_Hexamer.tsv -d ./Model/Plant.logit.RData -g candidate_transcript.fasta -o CPAT_plant.output





### **5.4. Alignment of sequences to the UniProt protein database with diamond.**
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    gunzip uniprot_sprot.fasta.gz
    diamond makedb --in uniprot_sprot.fasta -d uniprot_out
    diamond blastx -d uniprot_out -q candidate_transcript.fasta -o uniprotoutput.txt




### **5.5. By intersecting the results obtained from the aforementioned steps, a set of high-confidence lncRNAs were obtained.**

The id of the lncRNA

    Rscript insersection.sh candidate_transcript.txt candidate_lncRNA.txt CPAT_plant.output plant-lncFinder.txt uniprotoutput.txt
    
The lncRNA gtf file

    grep -Fwf final_lncRNA_results.txt candidate_transcript.gtf > lncRNA.gtf

	
	
	
## **6. Classify the final set of lncRNAs based on their genomic locations and sequence features.**
Classification result file

	FEELnc_classifier.pl -i lncRNA.gtf -a genome.gtf > lncRNA_classes.txt
	
	
Antisense_exonic-lncRNA

	awk -F '\t' '{if($6 == "antisense" && $10 == "exonic") {print $0}}' lncRNA_classes.txt > LncRNA_antisense_exonic.gtf


Intronic-lncRNA

	awk -F '\t' '{if($10 == "intronic") {print $0}}' lncRNA_classes.txt > LncRNA_intronic.gtf


Upstream-lncRNA

	awk -F '\t' '{if($7 == "intergenic" && $8 <= 2000 && $10 == "upstream") {print $0}}' lncRNA_classes.txt | awk -F '\t' '{if($8 <= 2000 && $10 == "upstream") print $0}' > LncRNA_upstream.gtf


Intergenic-lncRNA

	awk -F '\t' '{if( $7 == "intergenic" && $8 > 2000) {print $0}}' lncRNA_classes.txt > LncRNA_intergenic.gtf


Bidirectional-lncRNA

	awk -F '\t' '{if( $6 == "antisense" && $7 == "intergenic" && $8 <= 2000 && $9 == "divergent") {print $0}}' lncRNA_classes.txt > LncRNA_Bidirectional.gtf
	
	
	
	
## **7. TE-derived lncRNAs.**
cat TE.bed

	Chr1    15827287        15838845        LTR_Gypsy
	Chr1    13085455        13085593        LTR_Copia
	Chr1    11181821        11181959        LTR_Copia
	Chr1    20699111        20699248        LTR_Copia
	...
	
	
cat LncRNA.bed

	Chr1    11171031        11171031        MSTRG.2781.1
	Chr1    12199350        12199350        MSTRG.2973.1
	Chr1    13466928        13466928        MSTRG.3115.1
	Chr1    13838536        13838536        MSTRG.3127.1
	...
	
	
TE-lncRNA

	bedtools intersect -a lncRNA.bed -b TE.bed -wo | sort -u | awk '{print $1,$2,$3,$4,$6,$7,$8,$9}' | sed 's/ /\t/g' | sed '1iChr\tLncRNA_start\tLncRNA_end\tLncRNA_ID\tTE_start\tTE_end\tTE_ID\tOverlap' > TE_lncRNA_intersect.txt 

	
	
	
	
	
	
	
	
