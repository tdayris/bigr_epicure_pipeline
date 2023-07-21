# 0.12.2 – Fixes

## Fix

* Samtools stats regions

# 0.12.1 – No input

## Features

* Design without input file is now handled.

## Fix

* Expand not allowed in params, replaced by list
* Bigr Launcher now finds Conda


# 0.12.0 – Quality Controls

## Features

* New QC in CSaw (cout, filters, normalize)
* New QC in EdgeR (disp, qldisp, tables)
* New QC before CSaw import
* Explicit fragment length estimation with corresponding graph
* Rename and concat file, automatically retrieve files from iRODS in Gustave Roussy

## Fix

* CSaw normalize now returns GRanges in any case

# 0.11.2 – Fixes

## Documentation

* Mat&Met updates

## Fix

* Gene body coverage missing tag matrix fixed
* Fixed OOM occurring on gene body coverage analysis
* Deeptools PCA arguments fixed
* ROI filter now considers only regions defined in the reference genomic fasta file

# 0.11.1 – Flamingo-launcher

## Documentation

* Job composer in open-on-demand

## Fix

* Launcher version

# 0.11.0 – ROI

## Features

* Select regions of interest
* Merge sample sharing condition for heatmaps
* Xenome for BiGR flamingo cluster only

## Changes

* Log2(threshold) changed from 3 to 1.1 in csaw filters

# 0.10.0 – Motifs

## Features

* Homer find motifs installation enhanced

# 0.9.0 – Cut&Run, Cut&Tag

## Features

* Enable SEACR
* Deal with deduplication in Cut&Tag/Cut&Run
* Enable schema validation for both config and design file
* Rulegrpahs available in readme

## Changes

* DeepTools bins have been reduced from 50 to 5.

# 0.8.0 – Motifs

## Features

* Enable motif analysis with Homer
* Fixes signle-end designs

# 0.7.0 – Differential binding

## Features

* Enable csaw
* Ensble MEDIPS
* deeptools Fingerprint
* deeptools Correlation
* deeptools PCA

## Fix

* Wildcard error fixed in Medipseq DB

# 0.6.0 – Annotation

## Featues

* ChIPSeeker annotation
* ChIPSeeker distance to tss graph
* DeepTools fingerprint metrics for MultiQC
* Upset-plot in ChIPSeeker
* Cumulative histogram in ChIPSeeker
* New gene body and genome coverage plots

## Fix

* Sambamba index bam key error fixed
* Sambamba filter output format
* Proper closing syntax in R scripts
* Unsorted bam before markduplicates
* Merged blacklist bed files in order to fix spurious normalization in DeepTools
* MultiQC mapping output
* Fix Macs2 IO with input files
* Fix annotation syntax error
* Fix KeyError in chipseeker plot
* Fix naming scheme
* Fix missing ggupset library

# 0.5.0 – Single-sample Peak Calling

## Fix

* Macs2: fixed ambiguous rule resulution

## Features

* Single-sample calling
* csaw initial release

# 0.4.1 – Mapping

## Fix

* Snakemake logging does not go to stdout anymore
* Bam indexed before Samtools sequence retrieval 

# 0.4.0 – Bins

## Features

* Rose binsize in DeepTools from 25 to 50 to reduce matrix size
* MultiQC now sticks to the input file list
* General bam index rule that works for all `{tool}/{subcommand}/{sample}.bam` files
* `pipeline_lint.sh` now also performs a dry-run


## Fix

* MultiQC now detects all log files
* Picard missing input files


## Misc

* Blacklist download protocol changed


# 0.3.2 – BiGR

## Features

* CSAW scripts available
* Downgrade Python to 3.9

## Fix

* Fastp SE/PE detection fixed
* Bowtie2 expected input files now includes complete list of index files
* MultiQC now looks for correct fastp repository


## Known bug

* FastScreen index download creates sub-directories leading to Snakemake error


# 0.3.1 – Fixes

## Features

* General usage documentation
* Report available
* Medip-seq DGE + Annotation
* Linter script

## Fix

* User-defined fastq-screen configuration is not used when provided.
* Sambamba filters now consider bam index
* Missing modules in MultiQC
* Multiple linter suggestions accepted


# 0.3 – Coverage

## Features

* DeepTools coverage and Heatmap
* MeDIPS initial commit

# 0.2 – Mapping

## Features

* Fastp trimming
* Bowtie2 mapping for ChIP, MeDIP, Atac.
* Sambamba/Samtools/DeepTools quality filters
* Picard/Samtools/DeepTools QC
* MultiQC report

# 0.1 – Indexes

## Features

* Download genome sequences
* Index genome sequences for Bowtie2, FastqScreen, Picard, and Samtools