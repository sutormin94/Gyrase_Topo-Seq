#!bin/bash

##############
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#Shell script that makes QC of the reads before and after the trimming procedure. 
#Than script maps trimmed and only paired reads to the reference genome, and prepares
#sorted and indexed BAM-files suitable for visualization with IGV

#Requirements: factqc, trimmomatic, bwa mem, samtools 
#This variables should be in the path (or replace them with the path to the particular program)
##############


#Set path to different files
#Path to the working directory, contains /Data folder with raw data
Sample_path=''
print $Sample_path
cd $Sample_path
#Path to the file containing sequencing adapters sequences for trimmomatic uses. Typically in the Trimmomatic-0.36/adapters/All_TruSeq.fa
Adapters=''
#Path to the reference genome
Ref_genome=''


#Initial quality control
mkdir Fastqc_analysis/
fastqc -t 20 -o $Sample_path/Fastqc_analysis/ $Sample_path/Data/*

#Reads trimming
mkdir $Sample_path/Data/Trimmed_gentely/
for i in `ls -a  $Sample_path/Data/ | sed -r "s/(.+)_R[1,2]_00.*/\1/g" | uniq | sort -d`; do
print $i
java -jar trimmomatic PE -threads 10 -phred33 $Sample_path/Data/${i}_R1_001.fastq.gz $Sample_path/Data/${i}_R2_001.fastq.gz $Sample_path/Data/Trimmed_gentely/${i}_paired_R1_001.fastq.gz $Sample_path/Data/Trimmed_gentely/${i}_unpaired_R1_001.fastq.gz $Sample_path/Data/Trimmed_gentely/${i}_paired_R2_001.fastq.gz $Sample_path/Data/Trimmed_gentely/${i}_unpaired_R2_001.fastq.gz ILLUMINACLIP:$Adapters:2:30:10 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:0 MINLEN:30 ; done

#Quality control after the trimming procedure
mkdir Fastqc_analysis/Trimmed_gentely/
fastqc -t 20 -o $Sample_path/Fastqc_analysis/Trimmed_gentely/ $Sample_path/Data/Trimmed_gentely/*

#Reads mapping to the reference genome: make SAM-files
mkdir SAM/
for i in `ls -a $Sample_path/Data/Trimmed_gentely/ | sed -r "s/(.+_paired)_R[1,2]_00.*/\1/g" | uniq | sort -d`; do 
bwa mem -t 20 Ref_genome $Sample_path/Data/Trimmed_gentely/${i}_R1_001.fastq.gz $Sample_path/Data/Trimmed_gentely/${i}_R2_001.fastq.gz > $Sample_path/SAM/$i.sam; done


#Prepares tracks for IGV: makes BAM-files, sorts them, makes index-files
mkdir SAM_sorted/
mkdir BAM/
#Makes BAM-files
for i in `ls -a $Sample_path/SAM/ | sed -r "s/(.+_paired).*/\1/g"`; do 
samtools view -S -b $Sample_path/SAM/${i}.sam > $Sample_path/BAM/${i}.bam ; done
#Sorts BAM-files
for i in `ls -a $Sample_path/BAM/ | sed -r "s/(.+_paired).*/\1/g"`; do 
samtools sort $Sample_path/BAM/${i}.bam $Sample_path/BAM_sorted/${i}_sorted.bam ; done
#Makes index files
for i in `ls -a $Sample_path/BAM_sorted/ | sed -r "s/(.+_paired_sorted).*/\1/g"`; do 
samtools index $Sample_path/BAM_sorted/${i}.bam ; done


