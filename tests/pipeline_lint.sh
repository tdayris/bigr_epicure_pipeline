set -e

snakefmt ../workflow/rules/{{coverage,differential_peak_calling,indexing,mapping,mapping_qc,medips,multiqc,peak-calling,reference,trimming}/*.smk,*.smk}
black ../workflow/scripts/jobscripts/*.py

export SNAKEMAKE_OUTPUT_CACHE="."

touch {I,S}{1,2}.R{1,2}.fq.gz
echo -e "Sample_id\tUpstream_file\tDownstream_file\tCondition\tInput\tFragment_size" > ../config/design.tsv
echo -e "S1\S1.R1.fq.gz\tS1.R2.fq.gz\tWT\tI1\t" >> ../config/design.tsv
echo -e "S2\S2.R1.fq.gz\tS2.R2.fq.gz\tWT\tI1\t" >> ../config/design.tsv
echo -e "I1\I1.R1.fq.gz\tI1.R2.fq.gz\t\t\t" >> ../config/design.tsv

echo -e "S3\S3.R1.fq.gz\tS3.R2.fq.gz\tMut\tI2\t" >> ../config/design.tsv
echo -e "S4\S4.R1.fq.gz\tS4.R2.fq.gz\tMut\tI2\t" >> ../config/design.tsv
echo -e "I2\I2.R1.fq.gz\tI2.R2.fq.gz\t\t\t" >> ../config/design.tsv

snakemake -c 1 --use-conda -s ../workflow/Snakefile --cache --lint

unset SNAKEMAKE_OUTPUT_CACHE