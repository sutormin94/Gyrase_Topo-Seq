# Gyrase_Topo-seq
Exploring gyrase cleavage sites accross E. coli W3110 genome

DNA-gyrase is a type II topoisomerase that introduces negative superhelicity into topologically closed DNA molecules. It operates with two DNA segments - so called G and T. During catalysis the enzyme introduces temporal double-stranded break into the G-segment, transferres the T-segment through it and religates the gap. The 5'-ends of the DNA break are stabilized by the formation of the intermediate covalent complex between DNA and gyrase.
Topo-Seq is a ChIP-Seq-like approach that exploits the formation the complexes to map the gyrase cleavage sites with a single-base precision.

This repositorium contains a set of bash, python and R scripts that were used for Topo-Seq data analysis and visualization.

######################

## Raw_reads_processing.sh

Shell script that makes QC of the reads before and after the trimming procedure. Than script maps trimmed and paired reads to the reference genome, prepares sorted and indexed BAM-files suitable for visualization with IGV
Requirements: factqc, trimmomatic, bwa mem, samtools 

######################

## SAM_to_coverage_and_N5E_N3E.py

Script takes SAM files as input, performs QC filtering of reads relying on the alignment quality and a presence of the partner: only reads pairs that have a score<256 are stored. Than the script computes coverage depth for DNA chains separately and for both. Additionally it calculates N5E (number of DNA fragments starts) and N3E (number of DNA fragments ends) values for every genome position. Coverage depth, N3E and N5E info returns as WIG files.

######################

## GCSs_calling.py

The script takes WIG files tetrade that contain N3E or N5E values: A+IP+, A+IP-, A-IP+, A-IP-. It smooths A+IP- and A-IP- tracks and divides A+IP+ and A-IP+ by them. Once A+IP+_div and A-IP+_div are obtained the script performs Audic-Clavery statistic test (Audic & Claverie, 1997) and returns regions of A+IP+_div where i and i+5 positions are significantly higher than corresponding in A-IP+. These regions from now are called GCSs. GCSs are stored in the output TXT file. Also two plots are generated: 1) signal coverage over the genome for treated and untreated samples; 2) Motif expected to be under the GCSs.
Requirements: TAB file with deletions coordinates.

## GCSs_filtering_and_overlapping.py

The script takes raw GCSs data, returns only trusted GCSs, computes GCSs shared between different conditions, draws Venn diagrams of the sets overlappings, writes GCSs sets.

######################
