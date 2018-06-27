#!bin/bash

##############
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#Shell script that makes QC of the reads before and after the trimming procedure. 
#Than script maps trimmed only paired reads to the reference genome, prepares
#sorted and indexed BAM-files suitable for visualization with IGV

#Requirements: factqc, trimmomatic, bwa mem, samtools 
#This variables should be in the path (or replace them with the path to the particular program)
##############


#######
#Variables to be defined.
#######

#Path to the working directory, contains /Raw_data folder with raw reads files.
PWD='/data/Gyrase/Data_preparation/Cfx_10mkM'
echo $PWD
cd $PWD

#Path to the file containing sequencing adapters sequences for trimmomatic uses. Typically in the Trimmomatic-0.36/adapters/XXX.fa
Adapters='/home/microcin/Programms/Trimmomatic-0.36/adapters/All_TruSeq.fa'
trimmomatic='/home/microcin/Programms/Trimmomatic-0.36/trimmomatic-0.36.jar'
#Path to the reference genome.
Ref_genome='/data/Gyrase/Genomes_tracks/E_coli_w3110_G_Mu.fasta'

#######
#Quality control and sequencing data preparation.
#######

#Initial quality control
echo '
#######################
Initial quality control is in progress...
#######################
'	
mkdir $PWD/Fastqc_analysis/
fastqc -t 20 -o $PWD/Fastqc_analysis/ $PWD/Raw_data/*

#######
#Reads mapping, alignment conversion to IGV-compatible format (sorted indexed BAM).
#######

#Reads mapping to the reference genome: make SAM-files
echo '
#######################
Reads mapping, SAM files generation...
#######################
'
mkdir $PWD/SAM/
for i in `ls -a $PWD/Raw_data/ | grep 'fastq.gz' | sed -r "s/(.+)_R[1,2]\.fastq\.gz/\1/g" | uniq | sort -d`; do 
bwa mem -t 20 $Ref_genome $PWD/Raw_data/${i}_R1.fastq.gz $PWD/Raw_data/${i}_R2.fastq.gz > $PWD/SAM/$i.sam; done

#Prepares tracks for IGV: makes BAM-files, sorts them, makes index-files
mkdir $PWD/BAM_sorted/
mkdir $PWD/BAM/
#Makes BAM-files
echo '
#######################
BAM files preparation...
#######################
'
for i in `ls -a $PWD/SAM/ | grep '.sam' | sed -r "s/(.+).sam/\1/g"`; do 
samtools view -S -b $PWD/SAM/${i}.sam > $PWD/BAM/${i}.bam ; done
#Sorts BAM-files
echo '
#######################
BAM files sorting...
#######################
'
for i in `ls -a $PWD/BAM/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do 
samtools sort $PWD/BAM/${i}.bam $PWD/BAM_sorted/${i}_sorted.bam ; done
#Makes index files
echo '
#######################
BAM files indexing...
#######################
'
for i in `ls -a $PWD/BAM_sorted/`; do 
samtools index $PWD/BAM_sorted/${i} ; done

echo '
Script ended its work succesfully!
'

"""
#Reads trimming
echo '
#######################
Reads trimming...
#######################
'
mkdir $PWD/Trimmed_gentely/
for i in `ls -a  $PWD/Raw_data/ | grep 'fastq.gz' | sed -r "s/(.+)_R[1,2]\.fastq\.gz/\1/g" | uniq | sort -d`; do
echo $i
java -jar $trimmomatic PE -threads 10 -phred33 $PWD/Raw_data/${i}_R1.fastq.gz $PWD/Raw_data/${i}_R2.fastq.gz $PWD/Trimmed_gentely/${i}_paired_R1.fastq.gz $PWD/Trimmed_gentely/${i}_unpaired_R1.fastq.gz $PWD/Trimmed_gentely/${i}_paired_R2.fastq.gz $PWD/Trimmed_gentely/${i}_unpaired_R2.fastq.gz ILLUMINACLIP:$Adapters:2:30:10 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:0 MINLEN:30 ; done

#Quality control after the trimming procedure
echo '
#######################
Quality control after trimming...
#######################
'
mkdir $PWD/Fastqc_analysis/Trimmed_gentely/
fastqc -t 20 -o $PWD/Fastqc_analysis/Trimmed_gentely/ $PWD/Trimmed_gentely/*

#######
#Reads mapping, alignment conversion to IGV-compatible format (sorted indexed BAM).
#######

#Reads mapping to the reference genome: make SAM-files
echo '
#######################
Reads mapping, SAM files generation...
#######################
'
mkdir $PWD/SAM/
for i in `ls -a $PWD/Trimmed_gentely/ | grep '_paired_' | sed -r "s/(.+)_paired_R[1,2]\.fastq\.gz/\1/g" | uniq | sort -d`; do 
bwa mem -t 20 $Ref_genome $PWD/Trimmed_gentely/${i}_paired_R1.fastq.gz $PWD/Trimmed_gentely/${i}_paired_R2.fastq.gz > $PWD/SAM/$i.sam; done

#Prepares tracks for IGV: makes BAM-files, sorts them, makes index-files
mkdir $PWD/BAM_sorted/
mkdir $PWD/BAM/
#Makes BAM-files
echo '
#######################
BAM files preparation...
#######################
'
for i in `ls -a $PWD/SAM/ | grep '.sam' | sed -r "s/(.+).sam/\1/g"`; do 
samtools view -S -b $PWD/SAM/${i}.sam > $PWD/BAM/${i}.bam ; done
#Sorts BAM-files
echo '
#######################
BAM files sorting...
#######################
'
for i in `ls -a $PWD/BAM/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do 
samtools sort $PWD/BAM/${i}.bam $PWD/BAM_sorted/${i}_sorted.bam ; done
#Makes index files
echo '
#######################
BAM files indexing...
#######################
'
for i in `ls -a $PWD/BAM_sorted/`; do 
samtools index $PWD/BAM_sorted/${i} ; done

echo '
Script ended its work succesfully!
'
"""
