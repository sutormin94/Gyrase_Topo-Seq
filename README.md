# Gyrase_Topo-seq
Exploring gyrase cleavage sites accross *E. coli W3110* genome

DNA-gyrase is a type II topoisomerase that introduces negative superhelicity into topologically closed DNA molecules. It operates with two DNA segments - so called G and T. During catalysis the enzyme introduces temporal double-stranded break into the G-segment, transferres the T-segment through it and religates the gap. The 5'-ends of the DNA break are stabilized by the formation of the intermediate covalent complex between DNA and gyrase.
Topo-Seq is a ChIP-Seq-like approach that exploits the formation the complexes to map the gyrase cleavage sites with a single-base precision.

This repository contains a set of bash, python and R scripts that were used for Topo-Seq data analysis and visualization.

######################

## Raw_reads_processing.sh

Shell script that makes QC of the reads before and after the trimming procedure. 
Than script maps trimmed and paired reads to the reference genome, prepares sorted and 
indexed BAM-files suitable for visualization with IGV.

**Requirements:** factqc, trimmomatic, bwa mem, samtools, shell

**Input:** Raw reads files (FASTQ), Genome file (FASTA)

**Output:** FastQC reports, SAM files, sorted indexed BAM files

######################

## SAM_to_coverage_and_N5E_N3E.py

Script takes SAM files as input, performs QC filtering of reads relying on the alignment quality and a presence of the partner: 
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

## GCSs_transcription_score_GC_distributions_throughout_genome.py

The script takes sets of trusted GCSs and analysis the distribution of GCSs throughout the genome. 
Also it plots the distribution of other values such as score, GC% and transcription.

**Requirements:** python 3

**Input:** TAB files with trusted GCSs, WIG genome score file, WIG genome GC file, TAB file with transcription data

**Output:** Plot with all the values examined distributed throughout the genome: GCSs number, score, GC, transcription level

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

**Input:** TAB files with trusted GCSs, FASTA genome file, FASTA file with a sequence to scan using the motif constructed

**Output:** A pair of PreWIG files with the scanning procedure results (for forward and reverse strands), FASTA file with sequences under GCSs were selected for combined motif construction,
Plot of the combined gyrase motif, TAB file with motif PFM (GC degenerate for further Fourier analysis) vertically oriented, 
TAB file with motif non-degenerate and corrected PFM horizontally oriented.

######################

## Prepare_score_track_GCSs_score_height_score_correlation.py

The script takes results of scanning procedure for forward and reverse strands.
Returns score for every genome position (writes into WIG file), for GCSs (writes TAB files contain coordinate\tN3E\tScore info), 
computes Pearson correlation between N3E and score, plots (Score, N3E) scatter plots.

**Requirements:** python 3

**Input:** TAB files with trusted GCSs, a pair of PreWIG files with the scanning procedure results (for forward and reverse strands), 
dataset name and chromosome identificator (both for output WIG)

**Output:** TAB files with trusted GCSs (score info added), WIG file with score defined for every genome position, 
(score, N3E) scatter plot with 1d trend line, correlation (score, N3E)

######################

## Genome_intervals_analysis.py

The script analysis sets of genome intervals (transcription units - TUs, BIMEs-1, BIMEs-2, IHF sites, Fis sites, H-NS sites, MatP sites, etc.)
for the enrichment of GCSs (binomial test), compares their N3E and score with mean GCSs N3E and score (t-test), 
compares intervals mean score with genome mean score (t-test).

**Requirements:** python 3

**Input:** TAB files with trusted GCSs and score info, WIG genome score file, TAB transcription units data files, TAB intervals data files

**Output:** TAB file with numbers of GCSs are associated with TUs compartments (USUS, USGB, GBDS, DSDS), 
TAB file with GCSs-TUs association analysis (GCSs number, GCSs N3E, GCSs score, TUs compartments score), TAB file with normalized numbers of GCSs are associated with TUs,
TAB file with the number of GCSs are associated with particular intervals set and statistics (GCSs number, GCSs N3E, GCSs score),
TAB file with the number of GCSs are associated with particular intervals (BIMEs-1, BIMEs-2), TAB file with intervals score statistics