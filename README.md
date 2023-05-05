# bigr_epicure_pipeline
Snakemake pipeline for Epicure analyses: Chip-seq, Atac-seq, Cut&Tag, Cut&Run


# Pipeline description

_note: The following steps may not be perform in that exact order._

## Pre-pocessing

| Step                        | Tool             | Documentation                                                                                                                                          | Reason                                                                                                                                                                                                                             |
| --------------------------- | ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Download genome sequence    | curl             | [Snakemake-Wrapper: download-sequence](https://snakemake-wrappers.readthedocs.io/en/v1.28.0/wrappers/reference/ensembl-sequence.html)                  | Ensure genome sequence are consistent in Epicure analyses                                                                                                                                                                          |
| Download genome annotation  | curl             | [Snakemake-Wrapper: download-annotation](https://snakemake-wrappers.readthedocs.io/en/v1.28.0/wrappers/reference/ensembl-annotation.html)              | Ensure genome annotation are consistent in Epicure analyses                                                                                                                                                                        |
| Download blacklised regions | manual shell FTP |                                                                                                                                                        | Ensure blacklist regions are consistent in Epicure analyses                                                                                                                                                                        |
| Trimming + QC               | Fastp            | [Snakemake-Wrapper: fastp](https://snakemake-wrappers.readthedocs.io/en/v1.28.0/wrappers/fastp.html)                                                   | Perform read quality check and corrections, UMI, adapter removal, QC before and after trimming                                                                                                                                     |


## Read mapping

| Step                        | Tool             | Documentation                                                                                                                                          | Reason                                                                                                                                                                                                                             |
| --------------------------- | ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Indexation                     | Bowtie2          | [Snakemake-Wrapper: bowtie2-build](https://snakemake-wrappers.readthedocs.io/en/v1.28.0/wrappers/bowtie2/build.html)                                   | Index genome for up-coming read mapping                                                                                                                                                                                            |
| Mapping                     | Bowtie2          | [Snakemake-Wrapper: bowtie2-align](https://snakemake-wrappers.readthedocs.io/en/v1.28.0/wrappers/bowtie2/align.html)                                   | Align sequenced reads on the genome                                                                                                                                                                                                |
| Filtering                     | Sambamba         | [Snakemake-Wrapper: sambamba-sort](https://snakemake-wrappers.readthedocs.io/en/v1.28.0/wrappers/sambamba/sort.html)                                   | Sort alignment over chromosome position, this reduce up-coming required computing resources, and reduce alignment-file size.                                                                                                       |
| Filtering                     | Sambamba         | [Snakemake-Wrapper: sambamba-view](https://snakemake-wrappers.readthedocs.io/en/v1.28.0/wrappers/sambamba/view.html)                                   | Remove non-canonical chromosomes and mapping belonging to mitochondrial chromosome.                                                                                                                                                |
| Filtering                     | Sambamba         | [Snakemake-Wrapper: sambamba-markdup](https://snakemake-wrappers.readthedocs.io/en/v1.28.0/wrappers/sambamba/markdup.html)                             | Remove sequencing duplicates.                                                                                                                                                                                                      |
| Filtering                     | DeepTools        | [Snakemake-Wrapper: sambamba-sort](https://snakemake-wrappers.readthedocs.io/en/v1.28.0/wrappers/deeptools/alignmentsieve.html)                        | For Atac-seq only. Reads on the positive strand should be shifted 4 bp to the right and reads on the negative strand should be shifted 5 bp to the left as in [Buenrostro et al. 2013](https://pubmed.ncbi.nlm.nih.gov/24097267/). |
| Archive                     | Sambamba         | [Snakemake-Wrapper: sambamba-view](https://snakemake-wrappers.readthedocs.io/en/v1.28.0/wrappers/sambamba/view.html)                                   | Compress alignment fil in CRAM format in order to reduce archive size.                                                                                                                                                             |
| Quality Control                  | Picard           | [Snakemake-Wrapper: picard-collect-multiple-metrics](https://snakemake-wrappers.readthedocs.io/en/v1.28.0/wrappers/picard/collectmultiplemetrics.html) | Summarize alignments, GC bias, insert size metrics, and quality score distribution. Performed before and after mapping post processing in order to highlight possible bias.                                                        |


## Coverage


| Step                        | Tool             | Documentation                                                                                                                                          | Reason                                                                                                                                                                                                                             |
| --------------------------- | ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Coverage                    | DeepTools        | [Snakemake-Wrapper: picard-collect-multiple-metrics](https://snakemake-wrappers.readthedocs.io/en/v1.28.0/wrappers/deeptools/bamcoverage.html)         | Compute genome coverage, normalized to 1M reads |


# Roadmap

* Preprocessing: FastqScreen
* Coverage: Fingerprint, Heatmaps (matrices), PCA, PBC
* Peak-calling: Macs2, Seacr, FRiP, FDR
* Peak-annotation: Homer, CentriMo
* Differential Peak Calling: DiffBind, DESeq2
* IGV: screen-shot, igv-reports
* Big freaking multiqc at the end!

Awaiting for SNv1.29.1 to be out.