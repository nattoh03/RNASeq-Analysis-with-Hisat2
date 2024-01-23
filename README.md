# RNASeq-Analysis-with-Hisat2



source 
# http://daehwankimlab.github.io/hisat2/
# (https://genome.cshlp.org/content/31/7/1290.full.pdf+html)



# RNA-seq analysis with with Subjuct pipeline

#1 Download the reference genome of An arabiensis from Vectorbase
/
wget https://vectorbase.org/common/downloads/Current_Release/AarabiensisDONGOLA2021/fasta/data/VectorBase-57_AarabiensisDONGOLA2021_Genome.fasta

#2 indexing the reference genome with subhead-build
subread-buildindex VectorBase-57_AarabiensisDONGOLA2021_Genome.fasta -o VectorBase-57_AarabiensisDONGOLA2021_Genome_index

#3  Download the gff file
wget https://vectorbase.org/common/downloads/Current_Release/AarabiensisDONGOLA2021/gff/data/VectorBase-57_AarabiensisDONGOLA2021.gff

#4 convert the gff file to gtf using gffread
gffread -T VectorBase-57_AarabiensisDONGOLA2021.gff -o VectorBase-57_AarabiensisDONGOLA2021.gtf

#5 quality control with fastqc. Several samples one or several samples  could be provided as input to fastqc 
fastqc *.fastq.gz 

#6 quality control filtering using 'fastp'. This command line below with trim bad quality reads and adapters
# the steps 6-10 could be run in one bash script
fastp \
        -i *_R1_001.fastq.gz -I *_R2_001.fastq.gz \
        -o *_R1_001_fp.fastq.gz -O *_R2_001_fp.fastq.gz \
        -h *_fp.html -l 25

#7 Alignment 
subjunc -T 4 -M 5 -d 50 -D 1500 -r Read_R1_001_fp.fastq.gz -R Read_R2_001_fp.fastq.gz -i VectorBase-57_AarabiensisDONGOLA2021_Genome_index -o Read_L001.bam

#8 To get the summary statistics about the alignment
samtools flagstat *.bam > *.bam.summary.txt


#9 filter bam to remove MAPQ<10 (-q 10) and unmapped reads (-F 4)
samtools view -b -q 10 -F 4 *.bam -o *_mq10.bam -@ 16

#10 sorting the filtered *bam file
samtools sort *_mq10.bam -o *_mq10.srt.bam -@16


#11 FeatureCount (part of the subread package)
featureCounts \
-T 16 \
-s 2 \
-p -B -C -Q 10 \
--largestOverlap --minOverlap 1 \
-t exon -g transcript_id \
-a VectorBase-57_AarabiensisDONGOLA2021.gtf -o Aarab.antisense.count.txt \
*_mq10.srt.bam
