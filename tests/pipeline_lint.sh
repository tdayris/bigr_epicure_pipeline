set -e

find ../workflow -type f -name "*.smk" | while read SNAKEFILE; do
    echo ${SNAKEFILE}
    snakefmt ${SNAKEFILE}
done
black ../workflow/scripts/*/*.py

export SNAKEMAKE_OUTPUT_CACHE="."

if [ ! -f "config/design.tsv" ]; then
    touch {I,S}{1,2}.R{1,2}.fq.gz S{3,4}.R{1,2}.fq.gz
    echo -e "Sample_id\tUpstream_file\tDownstream_file\tCondition\tInput\tFragment_size" > ../config/design.tsv
    echo -e "S1\tS1.R1.fq.gz\tS1.R2.fq.gz\tWT\tI1\t" >> ../config/design.tsv
    echo -e "S2\tS2.R1.fq.gz\tS2.R2.fq.gz\tWT\tI1\t" >> ../config/design.tsv
    echo -e "I1\tI1.R1.fq.gz\tI1.R2.fq.gz\t\t\t" >> ../config/design.tsv

    echo -e "S3\tS3.R1.fq.gz\tS3.R2.fq.gz\tMut\tI2\t" >> ../config/design.tsv
    echo -e "S4\tS4.R1.fq.gz\tS4.R2.fq.gz\tMut\tI2\t" >> ../config/design.tsv
    echo -e "I2\tI2.R1.fq.gz\tI2.R2.fq.gz\t\t\t" >> ../config/design.tsv
else
    echo "Design already exists"
fi

snakemake -c 1 --use-conda -s ../workflow/Snakefile --cache --lint


snakemake -c 1 --use-conda -s ../workflow/Snakefile --cache -n -p --configfile ../config/config.yaml

unset SNAKEMAKE_OUTPUT_CACHE