# 0.6.0 – SPC Fixes

## Featues

* ChIPSeeker annotation
* ChIPSeeker distance to tss graph

## Fix

* Sambamba index bam key error fixed
* Sambamba filter output format
* Proper closing syntax in R scripts

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