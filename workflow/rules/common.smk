import logging
import os
import pandas

from snakemake.remote import FTP
from typing import Any, Dict


########################
### Load design file ###
########################


design = pandas.read_csv(
    config.get("design", "../config/design.tsv"),
    sep="\t",
    header=0,
    index_col=0,
)


########################
### Global variables ###
########################

# FTP handler for possible downloads
FTP = FTP.RemoteProvider()


# This let all temporary file end-up in the same directory
tmp = os.environ.get("TMP", "temp")


# Mouse has less chromosomes, but this will work anyway
# Does not work for organisms with more than 22 autosomal chromosomes
# Does not work for organisms with no X/Y system for sexual chromosomes
canonical_chromosomes += list(map(str, range(1, 23))) + ["X", "Y"]
canonical_chromosomes += [f"chr{chrom}" for chrom in canonical_chromosomes]


# Protocol. It can be either : `chip-seq`, `atac-seq`, `cut&run`, or `cut&tag` 
protocol = config.get("protocol", "chip-seq")

################################
### Paths to reference files ###
################################

# Main genome informations
species = config.get("reference", {}).get("species", "homo_sapiens")
build = config.get("reference", {}).get("species", "GRCh38")
release = config.get("reference", {}).get("species", "109")

# Genome sequence (FASTA, not gzipped)
genome_fasta_path = config.get("reference", {}).get("genome_fasta")
if not genome_fasta_path:
    logging.info(
        "Missing genome FASTA path in the file `config.yaml`. "
        "A new one will be downloaded from Ensembl FTP."
    )
    genome_fasta_path = f"reference/{species}.{build}.{release}.fasta"


# Genome sequences indexes
genome_fai_path = genome_fasta_path + ".fai"
genome_dict_path = ".".join(genome_fasta_path.split(".")[:-1]) + ".dict"


# Genome annotation (GFF/GTF, not gzipped)
genome_annotation_path = config.get("reference", {}).get("genome_gtf")
if not genome_annotation_path:
    logging.info(
        "Missing GTF path in the file `config.yaml`. "
        "A new one will be downloaded from Ensembl FTP."
    )
    genome_annotation_path = f"reference/{species}.{build}.{release}.gtf"


# Bowtie2 genome index
bowtie2_index_path = config.get("reference", {}).get("bowtie2_index")
if not bowtie2_index_path:
    logging.info(
        "Missing Bowtie2 index path in the file `config.yaml`. "
        "A new one will be created."
    )
    bowtie2_index = multiext(
        f"reference/bowtie2_index/{species}.{build}.{release}",
        ".1.bt2",
        ".2.bt2",
        ".3.bt2",
        ".4.bt2",
        ".rev.1.bt2",
        ".rev.2.bt2",
    )


# Genome blacklisted regions
blacklist_path = config.get("reference", {}).get("blacklist")
if not blacklist_path:
    blacklist_path = f"reference/blacklist/{species}.{build}.{release}.bed.gz"

###################
### IO function ###
###################


def get_fastp_input(wildcards) -> Dict[str, List[str]]:
    """
    Return the list of Fastp input files
    """
    fastp_input = {
        "sample": [design["Upstream_file"].loc[wildcards.sample]],
    }
    if "Downstream_file" in design.keys():
        down = design["Downstream_file"].loc[wildcards.sample]
        if down is not None and down != "":
            fastp_input["sample"].append(down)
    return fastp_input


def get_bowtie2_align_input(wildcards, bowtie2_index: str = bowtie2_index) -> Dict[str, List[str]]:
    """
    Return the list of bowtie2 align input
    """
    bowtie2_align_input = {
        "idx": bowtie2_index_path
    }
    if "Downstream_file" in design.keys():
        down = design["Downstream_file"].loc[wildcards.sample]
        if down is not None and down != "":
            bowtie2_align_input["sample"] = expand(
                "fastp/trimmed/pe/{sample}.{stream}.fastq",
                stream=["1", "2"],
                sample=[wildcards.sample]
            )
        else:
            bowtie2_align_input["sample"] = expand(
                "fastp/trimmed/se/{sample}.fastq",
                sample=[wildcards.sample]
            )

    return bowtie2_align_input
    


def targets(config: Dict[str, Any]):
    """
    Return the list of expected output files, depending on the
    choices made by user in configuration file at: `<root>/config/config.yaml`
    """
    expected_targets = {}

    if config.get("steps", {}).get("install", False):
        expected_targets["gtf"] = genome_annotation_path
        expected_targets["genome_fasta"] = genome_fasta_path
        expected_targets["bowtie2_index"] = bowtie2_index_path

    if config.get("step", {}).get("mapping", False):
        if protocol == "atac-seq":
            expected_targets["mapping"] = "data_output/CRAM/{sample}.shift.cram"
        else:
            expected_targets["mapping"] = "data_output/CRAM/{sample}.cram"

    