#!/usr/bin/env bash
module load fastqc
mkdir /home/nattohz/Fun_RNASeq/Qc_Raw_reads_all
fastqc /opt/data/nattohz/Raw_reads/*.fastq.gz \
-o /home/nattohz/Fun_RNASeq/Qc_Raw_reads_all
##
##i want to make sure the module/program used for trimming is present by calling it###########
module load trimmomatic/0.39
## I also want to create a folder where all my trimmed reads will be stored after trimming ##############
mkdir /home/nattohz/Fun_RNASeq/Trimmed_reads_all_samples
for r1 in *_R1.fastq.gz
do
echo $r1
sample=$(basename $r1)
sample=${sample%_R1.fastq.gz}
echo "Processing sample: "$sample
trimmomatic PE ${sample}_R1.fastq.gz ${sample}_R1.fastq.gz \
/home/nattohz/Fun_RNASeq/Trimmed_reads_all_samples/${sample}_trimmed_forward_paired.fq.gz \
/home/nattohz/Fun_RNASeq/Trimmed_reads_all_samples/${sample}_trimmed_forward_unpaired.fq.gz \
/home/nattohz/Fun_RNASeq/Trimmed_reads_all_samples/${sample}_trimmed_reverse_paired.fq.gz \
/home/nattohz/Fun_RNASeq/Trimmed_reads_all_samples/${sample}_trimmed_reverse_unpaired.fq.gz \
ILLUMINACLIP:/home/nattohz/Fun_RNASeq/Adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:30
##
##
module load fastqc
mkdir /home/nattohz/Fun_RNASeq/Qc_Trimmed_all_reads
fastqc \
/home/nattohz/Fun_RNASeq/Trimmed_reads_all_samples/*_paired.fq.gz \
-o /home/nattohz/Fun_RNASeq/Qc_Trimmed_all_reads
##
##
##
##
############### Alignment with the hisat2 ###########################################################
###first make a directory where my index files will be located
mkdir /home/nattohz/Fun_RNASeq/Refseq/index_outputs
## download An. funestus genome
wget http://ftp.ensemblgenomes.org/pub/metazoa/release-53/fasta/anopheles_funestus/dna/Anopheles_funestus.AfunF3.dna_sm.toplevel.fa.gz -o /home/nattohz/Fun_RNASeq/Refseq
##download Ano. funestus GFF file
wget http://ftp.ensemblgenomes.org/pub/metazoa/release-53/gff3/anopheles_funestus/Anopheles_funestus.AfunF3.53.gff3.gz -o /home/nattohz/Fun_RNASeq/Refseq
##download Ano. funestus cds file
wget http://ftp.ensemblgenomes.org/pub/metazoa/release-53/fasta/anopheles_funestus/cds/Anopheles_funestus.AfunF3.cds.all.fa.gz -o /home/nattohz/Fun_RNASeq/Refseq

### Build a HISAT2 index for latest Anopheles funestus reference genome using hisat2-build
hisat2-build /home/nattohz/Fun_RNASeq/Refseq/GCA_003951495.1_AfunF3_genomic.fna \
/home/nattohz/Fun_RNASeq/Refseq/index_outputs/GCA_003951495.1_AfunF3_genomic_hisat2.idx
##
##### Map the reads for all samples sample using HISAT2 #################################################
mkdir /home/nattohz/Fun_RNASeq/Sam_files_all
hisat2 \
-x /home/nattohz/Fun_RNASeq/Refseq/index_outputs/GCA_003951495.1_AfunF3_genomic_hisat2.idx \
-1 /home/nattohz/Fun_RNASeq/Trimmed_reads_all_samples/${sample}_trimmed_forward_paired.fq.gz \
-2 /home/nattohz/Fun_RNASeq/Trimmed_reads_all_samples/${sample}_trimmed_reverse_paired.fq.gz \
-S /home/nattohz/Fun_RNASeq/Sam_files_all/${sample}.sam
##
###### Convert the SAM file to a BAM file###################################3###########################
mkdir /home/nattohz/Fun_RNASeq/Bam_files_all
samtools view -S \
-o /home/nattohz/Fun_RNASeq/Bam_files_all/${sample}.bam \
-b /home/nattohz/Fun_RNASeq/Sam_files_all/${sample}.sam
##
##
##### Sort the BAM file##############################
samtools sort \
-o /home/nattohz/Fun_RNASeq/Sam_files_all/${sample}_sorted.bam /home/nattohz/Fun_RNASeq/Bam_files_all/${sample}.bam
##
##
#####Index the BAM file so that it can be read efficiently by IGV#####################
samtools index \
/home/nattohz/Fun_RNASeq/Bam_files_all/${sample}_sorted.bam
done







