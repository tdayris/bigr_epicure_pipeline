# General configuration

To configure this workflow, modify `config/config.yaml` according to your needs, following the explanations provided in the file.

Users working at Gustave Roussy may want to use pre-indexed files.

# Sample description, _aka._ design file

## Content

The design file is a [TSV](https://en.wikipedia.org/wiki/Tab-separated_values)-formatted text file. It contains the following columns:

1. `Sample_id`: Required. A unique name for each sample.
1. `Upstream_file`: Required. A unique path to the upstream reads. (Fastq formatted, can be gzipped)
1. `Downstream_file`: Optional. Missing value(s) suggest single-ended libraries. A unique path to the downstream reads. A design file can contain samples from both single-ended and pair-ended libraries. (Fastq formatted, can be gzipped)
1. `Fragment_size`: Required if and only if sample is single-ended. The expected size of the fragment.
1. `Input`: Optional Sample_id of the corresponding input file. If the current sample is an input file, leave this cell empty.
1. `Protocol`: Optional. Ignored: not yet implemented. All samples must belong to a same sequencing protocol defined in `config/config.yaml`.
1. `Organism`: Optional. Ignored: not yet implemented. All samples must belong to a same organism defined in `config/config.yaml`.
1. `Build`: Optional. Ignored: not yet implemented. All samples must belong to a same organism build defined in `config/config.yaml`.
1. `Release`: Optional. Ignored: not yet implemented. All samples must belong to a same organism release defined in `config/config.yaml`.
1. `...`: Optional. Ignored: not yet implemented. Any character factor name as column name, and any level character name as cell values.

## Pair-ended example

| Sample_id | Upstream_file | Downstream_file | Input | Status |
| --------- | ------------- | --------------- | ----- | ------ |
| S1        | S1.R1.fq.gz   | S1.R2.fq.gz     | I1    | Mut    |
| S2        | S2.R1.fq.gz   | S2.R2.fq.gz     | I2    | WT     |
| I1        | I1.R1.fq.gz   | I1.R2.fq.gz     |       |        |
| S3        | S3.R1.fq.gz   | S3.R2.fq.gz     | I1    | Mut    |
| S4        | S4.R1.fq.gz   | S4.R2.fq.gz     | I2    | WT     |
| I2        | I2.R1.fq.gz   | I2.R2.fq.gz     |       |        |

## Single-ended example

| Sample_id | Upstream_file | Input | Status | Fragment_size |
| --------- | ------------- | ----- | ------ | ------------- |
| S1        | S1.R1.fq.gz   | I1    | Mut    | 300           |
| S2        | S2.R1.fq.gz   | I2    | WT     | 300           |
| I1        | I1.R1.fq.gz   |       |        | 300           |
| S3        | S3.R1.fq.gz   | I1    | Mut    | 300           |
| S4        | S4.R1.fq.gz   | I2    | WT     | 300           |
| I2        | I2.R1.fq.gz   |       |        | 300           |

## Nested Single/Pair -ended example

| Sample_id | Upstream_file | Downstream_file | Input | Status | Fragment_size |
| --------- | ------------- | --------------- | ----- | ------ | ------------- |
| S1        | S1.R1.fq.gz   | S1.R2.fq.gz     | I1    | Mut    |               |
| S2        | S2.R1.fq.gz   | S2.R2.fq.gz     | I2    | WT     |               |
| I1        | I1.R1.fq.gz   | I1.R2.fq.gz     |       |        |               |
| S3        | S3.R1.fq.gz   |                 | I1    | Mut    | 300           |
| S4        | S4.R1.fq.gz   |                 | I2    | WT     | 300           |
| I2        | I2.R1.fq.gz   |                 |       |        | 300           |

