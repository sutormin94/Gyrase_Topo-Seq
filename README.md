# Gyrase_Topo-seq
Exploring gyrase cleavage sites across *E. coli W3110* genome

DNA-gyrase is a type II topoisomerase that introduces negative superhelicity into topologically closed DNA molecules. It operates with two DNA segments - so called G and T. 
During catalysis the enzyme introduces temporal double-stranded break into G-segment, transfers T-segment through it and religates the gap. 
5'-ends of the DNA break are stabilized by formation of an intermediate covalent complex between DNA and gyrase.
Topo-Seq is a ChIP-Seq-like approach that exploits formation of these intermediates to map the gyrase cleavage sites (GCSs) with a single-base precision.

This repository contains a set of bash, python and R scripts which were used for Topo-Seq data analysis and visualization. 
Raw sequencing data and some processed files (coverage depth WIG, N3E WIG, GCSs lists) can be retrieved from GEO datasets with accession GSE117186.


# Main pipeline

![alt text](https://github.com/sutormin94/Gyrase_Topo-seq/blob/master/Pipeline_overview/Main_pipeline.png)

## Raw_reads_processing.sh

Shell script that makes QC of the reads before and after the trimming procedure. 
Than script maps trimmed and paired reads to the reference genome, prepares sorted and 
indexed BAM-files suitable for visualization with IGV.

**Requirements:** factqc, trimmomatic, bwa mem, samtools, shell

**Input:** Raw reads files (FASTQ), Genome file (FASTA)

**Output:** FastQC reports, SAM files, sorted and indexed BAM files

######################

## SAM_to_coverage_and_N5E_N3E.py

Script takes SAM files as an input, performs QC filtering of reads relying on the alignment quality and a presence of the partner: 
only reads pairs that have a score<256 are stored. Than the script computes coverage depth for DNA chains separately and for both. 
Additionally it calculates N5E (number of DNA fragments starts) and N3E (number of DNA fragments ends) values for every genome position. 
Coverage depth, N3E and N5E info returns as WIG files.

**Requirements:** python 3

**Input:** SAM files, chromosome identificator (for output WIG)

**Output:** SAM files contain proper aligned reads only, TAB files with start coordinates of DNA fragments alignments and alignment lengths, a range of WIG files:
coverage depth for forward strand, reverse strand, and for both strands, N3E, N5E, and N3E+N5E data 

######################

## GCSs_calling.py

The script takes WIG files tetrade that contain N3E or N5E values: A+IP+, A+IP-, A-IP+, A-IP-. It smooths A+IP- and A-IP- tracks and 
divides A+IP+ and A-IP+ by them. Once A+IP+_div and A-IP+_div are obtained the script performs Audic-Clavery statistic test (Audic & Claverie, 1997) 
and returns regions of A+IP+_div where i and i+5 positions are significantly higher than corresponding in A-IP+. 
These regions from now are called GCSs. GCSs are stored in the output TXT file. Also two plots are generated: 1) 
signal coverage over the genome for treated and untreated samples; 2) Motif expected to be under the GCSs.

**Requirements:** python 3

**Input:** TAB file with deletions coordinates, N3E or N5E values containing WIG files forms a tetrade (A+IP+, A+IP-, A-IP+, A-IP-), FASTA genome file

**Output:** Plot shows coverage depth normalization before GCSs calling, TAB file with raw GCSs coordinates and N3E values, Plot with raw motif

######################

## GCSs_filtering_and_overlapping.py

The script takes raw GCSs data, returns only trusted GCSs, computes GCSs shared between different conditions, 
draws Venn diagrams of the sets overlappings, writes GCSs sets.

**Requirements:** python 3

**Input:** TAB files with raw GCSs (triplicates)

**Output:** TAB files with trusted GCSs, TAB files with shared GCSs (Cfx_Micro, Cfx_Oxo, Micro_Oxo, Cfx_MIcro_Oxo, Cfx_Rif_Cfx), 
Venn diagrams (Cfx vs Micro vs Oxo, Cfx vs RifCfx, 3 replicas for of the each experimental conditions)

######################

## Motifs_visualization_sequences_extraction.py

The script takes sets of trusted GCSs as input and plots motifs using the sequences under the GCSs.
Also it writes sequences and motif to files.

**Requirements:** python 3

**Input:** TAB files with trusted GCSs, FASTA genome file

**Output:** FASTA file with sequences under GCSs, TAB file with motif PFM (GC degenerate for further Fourier analysis), Gyrase motif plots constructed for trusted GCSs

######################

## Combined_motif_construction_scanning_plotting_writing.py

The script takes sets of trusted GCSs as input, filters GCSs with highest N3E, makes a combined set consists of these GCSs, returns sequences under them and constructs
PSSM matrix by the way getting rid of antibiotic-specific bias at positions forming the cleavage site. Than the script scans a sequence of interest with the PSSM, 
returns the results of scanning, plots combined motif and writes it in a GC% degenerate and in a non-degenerate forms.

**Requirements:** python 3

**Input:** TAB files with trusted GCSs, FASTA genome file, FASTA file with a sequence to scan using the motif constructed, 
dataset name and chromosome identificator (both for output WIG)

**Output:** WIG files with the scanning procedure results (score track), FASTA file with sequences under GCSs that were selected for combined motif construction,
Plot of the combined gyrase motif, TAB file with motif PFM (GC degenerate for further Fourier analysis) vertically oriented, 
TAB file with motif non-degenerate and corrected PFM horizontally oriented.

######################

## Return_GCSs_score_height_correlation.py

The script takes results of scanning procedure (WIG file) and 
returns score for GCSs (writes TAB files contain coordinate\tN3E\tScore info), 
computes Pearson correlation between N3E and score, plots (Score, N3E) scatter plots.

**Requirements:** python 3

**Input:** TAB files with trusted GCSs, WIG file with the scanning procedure results (score track).

**Output:** TAB files with trusted GCSs (score info added),
plot with distributions of N3E and score values for different GCSs sets and overall genome (the last is for score),
(score, N3E) scatter plot with 1d trend line, correlation (score, N3E)

######################

# Additional scripts

Further scripts are not included into the main straightforward pipeline and require additional data:
1. Genome annotation
2. Transcription data
3. GC% track for genome of interest
4. Binding sites of NAPs (Fis, IHF, H-NS, MatP, etc.)
5. Annotation of special regions (BIME-1s, BIME-2s)
6. Annotation of TAD borders
7. VCF file with mutations
8. Processed data for GCSs enrichment in different genome regions (GCSs_association_with_TUs_USUS_USGB_GBDS_DSDS_barplot.R script)

For *E. coli DY330* MuSGS all the additional files required are stored in **Additional_genome_features** folder.


## GCSs_transcription_score_GC_distributions_throughout_genome.py

The script takes sets of trusted GCSs and analyzes the distribution of GCSs throughout the genome. 
Also it plots the distribution of other values such as score, GC% and transcription.

**Requirements:** python 3

**Input:** TAB files with trusted GCSs, WIG genome score file, WIG genome GC file, TAB file with transcription data

**Output:** Plot with all the values examined distributed throughout the genome: GCSs number, score, GC, transcription level

######################

## Genome_intervals_analysis.py

The script analyzes sets of genome intervals (transcription units - TUs, BIME-1s, BIME-2s, IHF sites, Fis sites, H-NS sites, MatP sites, etc.)
for the enrichment of GCSs (binomial test), compares their N3E and score with mean GCSs N3E and score (t-test), 
compares intervals mean score with genome mean score (t-test).

**Note: ** Returns some warning messages due to the ommitting statistics (t-test) for too short sets of values.

**Requirements:** python 3

**Input:** TAB files with trusted GCSs and score info, WIG genome score file, TAB transcription units data files, TAB intervals data files

**Output:** TAB file with numbers of GCSs are associated with TUs compartments (USUS, USGB, GBDS, DSDS), 
TAB file with GCSs-TUs association analysis (GCSs number, GCSs N3E, GCSs score, TUs compartments score), TAB file with normalized numbers of GCSs are associated with TUs,
TAB file with the number of GCSs are associated with particular intervals set and statistics (GCSs number, GCSs N3E, GCSs score),
TAB file with the number of GCSs are associated with particular intervals (BIMEs-1, BIMEs-2), TAB file with intervals score statistics

######################

## Cfx_RifCfx_data_comparison.py

Script compares data from Cfx and RifCfx (conditions with transcription inhibited with rifampicin) experiments. It identifies GCSs
shared between datasets, computes whether the signal (N3E) goes up or down as a response for transcription inhibition. Also it 
analyzes shared GCSs that fall into BIMEs or DS regions of rRNA operons to be associated with signal increase or decrease.

**Requirements:** python 3

**Input:** TAB files with trusted GCSs and score info (Cfx and RifCfx), TAB transcription units data files (rRNA operons), TAB intervals data files (BIMEs-1, BIMEs-2)

**Output:** Plot representing N3E ratio (RifCfx/Cfx), STOUTPUT with statistical information about preferential increase or decrease of 
N3E as a response for transcription inhibition for GCSs subsets

######################

## TAD_Mut_assoc_GCSs_analysis.py

Script analyzes colocalization of mutations (from Foster, 2015) and TADs borders (from Lioy, 2018) with GCSs.

**Requirements:** python 3

**Input:** TAB files with trusted GCSs, BED file with TAD intervals, VCF file with mutations

**Output:** Statistics of GCSs associations

######################

## qPCR_statistics.R

Script that visualizes qPCR data: Cfx vs RifCfx compared for 4 locus - ccmH, MuSGS, rRNA A DS, rRNA A US.
Fold enrichment data is stored in the script body, raw data could be found in Summary_table.xlsx supplementary table, sheet DS1.
Script makes barplot and computes t-test statistic.

**Requirements:** R (tested on 3.4.3 "Kite-Eating Tree")

**Input:** Data is already stored in the script

**Output:** Plot represents qPCR data and t statistics

######################

## GCSs_association_with_TUs_USUS_USGB_GBDS_DSDS_barplot.R

Script visualizes GCSs association with different sets of transcription units (TUs):
AG - all genes, GLE - genes with low transcription level, GHE - genes with high level of transcription,
AO - all operons, OLE - operons characterized by low level of transcription, OHE - operons with high transcription level.
Simultaneously the script is revealing association of GCSs with a range of TUs regions:
Upstream, TU start, TU end, Downstream that are correspond to USUS, USGB, GBDS, DSDS.
The same analysis is performing with score of GCSs and average score of regions working with. Additionally 
GCSs association with rRNA operons is visualizing.

**Requirements:** R (tested on 3.4.3 "Kite-Eating Tree")

**Input:** Excel file with GCSs data and score data 

**Output:** Plots represent GCSs association data and score data

######################