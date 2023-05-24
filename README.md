# bigr_epicure_pipeline
Snakemake pipeline for Epicure analyses: Chip-Seq, Atac-Seq, Cut&Tag, Cut&Run, MeDIP-Seq, 8-OxoG-Seq


# Summary

1. [Usage]()
    1. [Installation](https://github.com/tdayris/bigr_epicure_pipeline#installation-following-snakemake-workflows-guidelines)
    1. [Deployment](https://github.com/tdayris/bigr_epicure_pipeline#deployment-following-snakemake-workflows-guidelines)
    1. [Configuration](https://github.com/tdayris/bigr_epicure_pipeline#configure-workflow-following-snakemake-workflows-guidelines)
    1. [Run this workflow](https://github.com/tdayris/bigr_epicure_pipeline#run-workflow-following-snakemake-workflows-guidelines)
    1. [Report](https://github.com/tdayris/bigr_epicure_pipeline#generate-report-following-snakemake-workflows-guidelines)
1. [Pipeline components](https://github.com/tdayris/bigr_epicure_pipeline#pipeline-description)
1. [Material and methods](https://github.com/tdayris/bigr_epicure_pipeline#material-and-methods)
    1. [ChIP-Seq](https://github.com/tdayris/bigr_epicure_pipeline#chip-seq)
    1. [Atac-Seq](https://github.com/tdayris/bigr_epicure_pipeline#atac-seq)
    1. [Cut&Tag](https://github.com/tdayris/bigr_epicure_pipeline#cuttag)
    1. [Cut&Run](https://github.com/tdayris/bigr_epicure_pipeline#cutrun)
    1. [MeDIP-Seq](https://github.com/tdayris/bigr_epicure_pipeline#medip-seq)
    1. [OG-Seq](https://github.com/tdayris/bigr_epicure_pipeline#og-seq)
1. [Up-comming features](https://github.com/tdayris/bigr_epicure_pipeline#roadmap)


# Usage

## Installation (following Snakemake-Workflows guidelines)

_note: This action has already been done for you if you work at Gustave Roussy. See at the end of this section_

1. Install [snakemake](https://snakemake.readthedocs.io) and [snakedeploy](https://snakedeploy.readthedocs.io/en/latest/) with [mamba](https://github.com/mamba-org/mamba) package manager. Given that Mamba is installed, run:

`mamba create -c conda-forge -c bioconda --name bigr_epicure_pipeline snakemake snakedeploy pandas`

2. Ensure your conda environment is activated in your bash terminal:

`conda shell.bash activate bigr_epicure_pipeline`

Alternatively, if you work at Gustave Roussy, you can use our shared environment:

`conda shell.bash activate /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/shared_conda/bigr_epicure_pipeline`


## Deployment (following Snakemake-Workflows guidelines)

_note: This action has been made easier for you if you work at Gustave Roussy. See at the end of this section_

Given that Snakemake and Snakedeploy are installed and available (see [Installation]()), the workflow can be deployed as follows.

1. Go to your working directory:

`cd path/to/my/project`

2. Deploy workflow:

`snakedeploy deploy-workflow https://github.com/tdayris/bigr_epicure_pipeline . --tag v0.3.0`

Snakedeploy will create two folders `workflow` and `config`. The former contains the deployment of the chosen workflow as a [Snakemake module](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-and-combining-pre-exising-workflows), the latter contains configuration files which will be modified in the next step in order to configure the workflow to your needs. Later, when executing the workflow, Snakemake will automatically find the main `Snakefile` in the `workflow` subfolder.

3. Consider to put the exact version of the pipeline and all modifications you might want perform under version control. _e.g._ by [managing it via a (private) Github repository](https://docs.github.com/en/github/importing-your-projects-to-github/adding-an-existing-project-to-github-using-the-command-line)


## Configure workflow (following Snakemake-Workflows guidelines)

See dedicated [`config/README.md`](https://github.com/tdayris/bigr_epicure_pipeline/blob/main/config/README.md) file for dedicated explanations of all options and consequences.

## Run workflow (following Snakemake-Workflows guidelines)

_note: This action has been made easier for you if you work at Gustave Roussy. See at the end of this section_

Given that the workflow has been properly deployed and configured, it can be executed as follows.

Fow running the workflow while deploying any necessary software via conda (using the Mamba package manager by default), run Snakemake with

`snakemake --cores all --use-conda `

Alternatively, for users at Gustave Roussy, you may use:

`bash workflow/scripts/misc/bigr_launcher.sh`

Snakemake will automatically detect the main `Snakefile` in the `workflow` subfolder and execute the workflow module that has been defined by the deployment in step 2.

For further options, e.g. for cluster and cloud execution, see [the docs](https://snakemake.readthedocs.io/).

## Generate report (following Snakemake-Workflows guidelines)

After finalizing your data analysis, you can automatically generate an interactive visual HTML report for inspection of results together with parameters and code inside of the browser via 

`snakemake --report report.zip`

The resulting `report.zip` file can be passed on to collaborators, provided as a supplementary file in publications, or uploaded to a service like [Zenodo](https://zenodo.org/) in order to obtain a citable [DOI](https://en.wikipedia.org/wiki/Digital_object_identifier). 

# Pipeline description

_note: The following steps may not be perform in that exact order._

## Pre-pocessing

| Step                        | Tool             | Documentation                                                                                                                             | Reason                                                                                         |
| --------------------------- | ---------------- | ----------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------- |
| Download genome sequence    | curl             | [Snakemake-Wrapper: download-sequence](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/reference/ensembl-sequence.html)     | Ensure genome sequence are consistent in Epicure analyses                                      |
| Download genome annotation  | curl             | [Snakemake-Wrapper: download-annotation](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/reference/ensembl-annotation.html) | Ensure genome annotation are consistent in Epicure analyses                                    |
| Download blacklised regions | manual shell FTP |                                                                                                                                           | Ensure blacklist regions are consistent in Epicure analyses                                    |
| Trimming + QC               | Fastp            | [Snakemake-Wrapper: fastp](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/fastp.html)                                      | Perform read quality check and corrections, UMI, adapter removal, QC before and after trimming |
| Quality Control             | FastqScreen      | [Snakemake-Wrapper: fastq-screen](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/fastq_screen.html)                               | Perform library quality check |


## Read mapping

| Step            | Tool      | Documentation                                                                                                                                          | Reason                                                                                                                                                                                                                             |
| --------------- | --------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Indexation      | Bowtie2   | [Snakemake-Wrapper: bowtie2-build](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/bowtie2/build.html)                                   | Index genome for up-coming read mapping                                                                                                                                                                                            |
| Mapping         | Bowtie2   | [Snakemake-Wrapper: bowtie2-align](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/bowtie2/align.html)                                   | Align sequenced reads on the genome                                                                                                                                                                                                |
| Filtering       | Sambamba  | [Snakemake-Wrapper: sambamba-sort](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/sambamba/sort.html)                                   | Sort alignment over chromosome position, this reduce up-coming required computing resources, and reduce alignment-file size.                                                                                                       |
| Filtering       | Sambamba  | [Snakemake-Wrapper: sambamba-view](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/sambamba/view.html)                                   | Remove non-canonical chromosomes and mapping belonging to mitochondrial chromosome.                                                                                                                                                |
| Filtering       | Sambamba  | [Snakemake-Wrapper: sambamba-markdup](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/sambamba/markdup.html)                             | Remove sequencing duplicates.                                                                                                                                                                                                      |
| Filtering       | DeepTools | [Snakemake-Wrapper: sambamba-sort](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/deeptools/alignmentsieve.html)                        | For Atac-seq only. Reads on the positive strand should be shifted 4 bp to the right and reads on the negative strand should be shifted 5 bp to the left as in [Buenrostro et al. 2013](https://pubmed.ncbi.nlm.nih.gov/24097267/). |
| Archive         | Sambamba  | [Snakemake-Wrapper: sambamba-view](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/sambamba/view.html)                                   | Compress alignment fil in CRAM format in order to reduce archive size.                                                                                                                                                             |
| Quality Control | Picard    | [Snakemake-Wrapper: picard-collect-multiple-metrics](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/picard/collectmultiplemetrics.html) | Summarize alignments, GC bias, insert size metrics, and quality score distribution.                                                                                                                                                |
| Quality Control | Samtools  | [Snakemake-Wrapper: samtools-stats](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/samtools/stats.html)                                 | Summarize alignment metrics. Performed before and after mapping-post-processing in order to highlight possible bias.                                                                                                               |
| Quality Control | DeepTools | [Snakemake-Wrapper: deeptools-fingerprint](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/deeptools/plotfingerprint.html)               | Control imuno precipitation signal specificity. |


## Coverage


| Step            | Tool      | Documentation                                                                                                                            | Reason                                                                              |
| --------------- | --------- | ---------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| Coverage        | DeepTools | [Snakemake-Wrapper: deeptools-bamcoverage](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/deeptools/bamcoverage.html)     | Compute genome coverage, normalized to 1M reads                                     |
| Coverage        | MEDIPS    | Incoming                                                                                                                                 | Compute genome coverage with CpG density correction using MEDIPS (MeDIP-Seq only)   |
| Scaled-Coverage | DeepTools | [Snakemake-Wrapper: deeptools-computematrix](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/deeptools/computematrix.html) | Calculate scores per genomic regions. Used for heatmaps and profile coverage plots. |
| Heatmap         | DeepTools | [Snakemake-Wrapper: deeptools-plotheatmap](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/deeptools/plotheatmap.html)     | Plot heatmap and peak coverage among samples                                        |
| Depth        | DeepTools | [Snakemake-Wrapper: deeptools-plotheatmap](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/deeptools/plotcoverage.html)    | Assess the sequencing depth of given samples |


## Peak-Calling

| Step         | Tool | Documentation                                                                                                          | Reason |
| ------------ | ---- | ---------------------------------------------------------------------------------------------------------------------- | ------ |
| Peak-Calling | Mac2 | [Snakemake-Wrapper: macs2-callpeak](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/macs2/callpeak.html) | Search for significant peaks |


## Differential Peak Calling

| Step         | Tool | Documentation                                                                                                          | Reason |
| ------------ | ---- | ---------------------------------------------------------------------------------------------------------------------- | ------ |
| Peak-Calling | MEDIPS | Incoming | Search for significant variation in peak coverage with EdgeR (MeDIP-Seq only) |

# Material and Methods

## ChIP-Seq

Reads are trimmed using [fastp](https://github.com/OpenGene/fastp) version [0.23.2](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/fastp.html#software-dependencies). Trimmed reads are mapped over the genome of interest defined in the configuration file, using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) version [2.5.1](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/bowtie2/align.html#software-dependencies). 

Mapped reads are filtered, deduplicated and archived using [Sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html) version [1.0](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/sambamba/view.html#software-dependencies). Initial filters keep only mapped reads, and reads with a mapped mate. Mapping qualities must be equal or above 30, or are filtered out. Only canonical chromosomes are kept at this step. Deduplication removes duplicated reads and their mates.

A quality control is made _both before and after_ these filters to ensure no valuable data is lost. These quality controls are made using [Picard](https://broadinstitute.github.io/picard/command-line-overview.html#CollectMultipleMetrics) version [3.0.0](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/picard/collectmultiplemetrics.html#software-dependencies), and [Samtools](http://www.htslib.org/doc/samtools-stats.html) version [1.17](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/samtools/stats.html#software-dependencies).

Genome coverage was assessed with [DeepTools](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html?highlight=bamcoverage) version [3.5.1](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html?highlight=bamcoverage). The same tool was used to normalize tracks, and plot both quality controls and peak coverage heatmaps.

## Atac-Seq

Reads are trimmed using [fastp](https://github.com/OpenGene/fastp) version [0.23.2](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/fastp.html#software-dependencies). Trimmed reads are mapped over the genome of interest defined in the configuration file, using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) version [2.5.1](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/bowtie2/align.html#software-dependencies). 

Mapped reads are filtered, deduplicated and archived using [Sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html) version [1.0](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/sambamba/view.html#software-dependencies). Initial filters keep only mapped reads, and reads with a mapped mate. Mapping qualities must be equal or above 30, or are filtered out. Only canonical chromosomes are kept at this step. Deduplication removes duplicated reads and their mates. Read shifting was performed with [DeepTools](https://deeptools.readthedocs.io/en/develop/content/tools/alignmentSieve.html) version [3.5.1](https://deeptools.readthedocs.io/en/develop/content/tools/alignmentSieve.html), as described in [Buenrostro et al. 2013](https://pubmed.ncbi.nlm.nih.gov/24097267/).

A quality control is made _both before and after_ these filters to ensure no valuable data is lost. These quality controls are made using [Picard](https://broadinstitute.github.io/picard/command-line-overview.html#CollectMultipleMetrics) version [3.0.0](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/picard/collectmultiplemetrics.html#software-dependencies), and [Samtools](http://www.htslib.org/doc/samtools-stats.html) version [1.17](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/samtools/stats.html#software-dependencies).

Genome coverage was assessed with [DeepTools](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html?highlight=bamcoverage) version [3.5.1](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html?highlight=bamcoverage). The same tool was used to normalize tracks, and plot both quality controls and peak coverage heatmaps.

## Cut&Tag

Reads are trimmed using [fastp](https://github.com/OpenGene/fastp) version [0.23.2](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/fastp.html#software-dependencies). Trimmed reads are mapped over the genome of interest defined in the configuration file, using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) version [2.5.1](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/bowtie2/align.html#software-dependencies). 

Mapped reads are filtered, deduplicated and archived using [Sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html) version [1.0](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/sambamba/view.html#software-dependencies). Initial filters keep only mapped reads, and reads with a mapped mate. Mapping qualities must be equal or above 30, or are filtered out. Only canonical chromosomes are kept at this step. Deduplication removes duplicated reads and their mates. Read shifting was performed with [DeepTools](https://deeptools.readthedocs.io/en/develop/content/tools/alignmentSieve.html) version [3.5.1](https://deeptools.readthedocs.io/en/develop/content/tools/alignmentSieve.html), as described in [Buenrostro et al. 2013](https://pubmed.ncbi.nlm.nih.gov/24097267/).

A quality control is made _both before and after_ these filters to ensure no valuable data is lost. These quality controls are made using [Picard](https://broadinstitute.github.io/picard/command-line-overview.html#CollectMultipleMetrics) version [3.0.0](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/picard/collectmultiplemetrics.html#software-dependencies), and [Samtools](http://www.htslib.org/doc/samtools-stats.html) version [1.17](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/samtools/stats.html#software-dependencies).

Genome coverage was assessed with [DeepTools](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html?highlight=bamcoverage) version [3.5.1](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html?highlight=bamcoverage). The same tool was used to normalize tracks, and plot both quality controls and peak coverage heatmaps.

## Cut&Run

Reads are trimmed using [fastp](https://github.com/OpenGene/fastp) version [0.23.2](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/fastp.html#software-dependencies). Trimmed reads are mapped over the genome of interest defined in the configuration file, using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) version [2.5.1](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/bowtie2/align.html#software-dependencies). 

Mapped reads are filtered, deduplicated and archived using [Sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html) version [1.0](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/sambamba/view.html#software-dependencies). Initial filters keep only mapped reads, and reads with a mapped mate. Mapping qualities must be equal or above 30, or are filtered out. Only canonical chromosomes are kept at this step. Deduplication removes duplicated reads and their mates. Read shifting was performed with [DeepTools](https://deeptools.readthedocs.io/en/develop/content/tools/alignmentSieve.html) version [3.5.1](https://deeptools.readthedocs.io/en/develop/content/tools/alignmentSieve.html), as described in [Buenrostro et al. 2013](https://pubmed.ncbi.nlm.nih.gov/24097267/).

A quality control is made _both before and after_ these filters to ensure no valuable data is lost. These quality controls are made using [Picard](https://broadinstitute.github.io/picard/command-line-overview.html#CollectMultipleMetrics) version [3.0.0](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/picard/collectmultiplemetrics.html#software-dependencies), and [Samtools](http://www.htslib.org/doc/samtools-stats.html) version [1.17](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/samtools/stats.html#software-dependencies).

Genome coverage was assessed with [DeepTools](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html?highlight=bamcoverage) version [3.5.1](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html?highlight=bamcoverage). The same tool was used to normalize tracks, and plot both quality controls and peak coverage heatmaps.

## MeDIP-Seq

Reads are trimmed using [fastp](https://github.com/OpenGene/fastp) version [0.23.2](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/fastp.html#software-dependencies). Trimmed reads are mapped over the genome of interest defined in the configuration file, using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) version [2.5.1](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/bowtie2/align.html#software-dependencies). 

Mapped reads are filtered, deduplicated and archived using [Sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html) version [1.0](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/sambamba/view.html#software-dependencies). Initial filters keep only mapped reads, and reads with a mapped mate. Mapping qualities must be equal or above 30, or are filtered out. Only canonical chromosomes are kept at this step. Deduplication removes duplicated reads and their mates.

A quality control is made _both before and after_ these filters to ensure no valuable data is lost. These quality controls are made using [Picard](https://broadinstitute.github.io/picard/command-line-overview.html#CollectMultipleMetrics) version [3.0.0](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/picard/collectmultiplemetrics.html#software-dependencies), and [Samtools](http://www.htslib.org/doc/samtools-stats.html) version [1.17](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/wrappers/samtools/stats.html#software-dependencies).

Genome coverage was assessed with [MeDIPS](https://bioconductor.org/packages/release/bioc/vignettes/MEDIPS/inst/doc/MEDIPS.pdf) version [1.50.0](https://github.com/tdayris/bigr_epicure_pipeline/blob/25802a43cc212fa96d39d16e7cee6b98a356a71e/workflow/envs/medips.yaml#L8).

## OG-Seq

Not yes available.

## Citations

1. Fastp:`Chen, Shifu, et al. "fastp: an ultra-fast all-in-one FASTQ preprocessor." Bioinformatics 34.17 (2018): i884-i890.`
1. Bowtie2: `Langmead, Ben, and Steven L. Salzberg. "Fast gapped-read alignment with Bowtie 2." Nature methods 9.4 (2012): 357-359.`
1. Sambamba: `Tarasov, Artem, et al. "Sambamba: fast processing of NGS alignment formats." Bioinformatics 31.12 (2015): 2032-2034.`
1. Ensembl: `Martin, Fergal J., et al. "Ensembl 2023." Nucleic Acids Research 51.D1 (2023): D933-D941.`
1. Snakemake: `Köster, Johannes, and Sven Rahmann. "Snakemake—a scalable bioinformatics workflow engine." Bioinformatics 28.19 (2012): 2520-2522.`
1. Atac-Seq: `Buenrostro, Jason D., et al. "Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position." Nature methods 10.12 (2013): 1213-1218.`
1. MeDIPS: `Lienhard, Matthias, et al. "MEDIPS: genome-wide differential coverage analysis of sequencing data derived from DNA enrichment experiments." Bioinformatics 30.2 (2014): 284-286.`
1. [Snakemake-Wrappers](https://snakemake-wrappers.readthedocs.io/en/v1.31.1/index.html)
1. [Snakemake-Workflows](https://snakemake.github.io/snakemake-workflow-catalog/?rules=true)
1. [BiGR Epicure](https://github.com/tdayris/bigr_epicure_pipeline)

# Roadmap

* Coverage: Fingerprint, PCA, PBC
* Peak-calling: Seacr, FRiP, FDR
* Peak-annotation: Homer, CentriMo
* Differential Peak Calling: DiffBind, EdgeR, csaw
* IGV: screen-shot, igv-reports
* Big freaking multiqc at the end!

Based on Snakemake-Wrappers version 1.31.1