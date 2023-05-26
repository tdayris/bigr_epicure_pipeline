import logging
import os
import pandas
import snakemake

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from typing import Any, Dict, List, Optional, Tuple, Union


########################
### Load design file ###
########################


design: pandas.DataFrame = pandas.read_csv(
    config.get("designs", "config/design.tsv"),
    sep="\t",
    header=0,
    index_col=0,
)


########################
### Global variables ###
########################

# FTP handler for possible downloads
HTTPS = HTTPRemoteProvider()


# This let all temporary file end-up in the same directory
tmp: str = os.environ.get("TMP", "temp")


# Mouse has less chromosomes, but this will work anyway
# Does not work for organisms with more than 22 autosomal chromosomes
# Does not work for organisms with no X/Y system for sexual chromosomes
canonical_chromosomes = list(map(str, range(1, 23))) + ["X", "Y"]
canonical_chromosomes += [f"chr{chrom}" for chrom in canonical_chromosomes]


# Protocol. It can be either : `chip-seq`, `atac-seq`, `cut&run`, `cut&tag`, or `medip-seq`
protocol: str = config.get("protocol", "chip-seq").lower().replace(" ", "")

################################
### Paths to reference files ###
################################

# Main genome informations
species: str = config.get("reference", {}).get("species", "homo_sapiens")
build: str = config.get("reference", {}).get("build", "GRCh38")
release: str = config.get("reference", {}).get("release", "109")

# Genome sequence (FASTA, not gzipped)
genome_fasta_path: str = config.get("reference", {}).get("genome_fasta")
if not genome_fasta_path:
    logging.info(
        "Missing genome FASTA path in the file `config.yaml`. "
        "A new one will be downloaded from Ensembl FTP."
    )
    genome_fasta_path = f"reference/{species}.{build}.{release}.fasta"


# Genome sequences indexes
genome_fai_path: str = f"{genome_fasta_path}.fai"
genome_dict_path: str = ".".join(genome_fasta_path.split(".")[:-1]) + ".dict"
genome_twobit_path: str = ".".join(genome_fasta_path.split(".")[:-1]) + ".2bit"


# Genome annotation (GFF/GTF, not gzipped)
genome_annotation_path: str = config.get("reference", {}).get("genome_gtf")
if not genome_annotation_path:
    logging.info(
        "Missing GTF path in the file `config.yaml`. "
        "A new one will be downloaded from Ensembl FTP."
    )
    genome_annotation_path = f"reference/{species}.{build}.{release}.gtf"


# Bowtie2 genome index
bowtie2_index_path: List[str] = config.get("reference", {}).get("bowtie2_index")
if not bowtie2_index_path:
    logging.info(
        "Missing Bowtie2 index path in the file `config.yaml`. "
        "A new one will be created."
    )
    bowtie2_index_path = multiext(
        f"reference/bowtie2_index/{species}.{build}.{release}",
        ".1.bt2",
        ".2.bt2",
        ".3.bt2",
        ".4.bt2",
        ".rev.1.bt2",
        ".rev.2.bt2",
    )


# Genome blacklisted regions
blacklist_path: str = config.get("reference", {}).get("blacklist")
if not blacklist_path:
    blacklist_path = f"reference/blacklist/{species}.{build}.{release}.bed.gz"


# See: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
effective_genome_size: int = config.get("library", {}).get("effective_genome_size")
if not effective_genome_size:
    default_effective_genome_size = {
        "GRCz10": {
            "50": 1_195_445_591,
            "75": 1_251_132_686,
            "100": 1_280_189_044,
            "150": 1_312_207_169,
            "200": 1_321_355_241,
        },
        "WBcel235": {
            "50": 95_159_452,
            "75": 96_945_445,
            "100": 98_259_998,
            "150": 98_721_253,
            "200": 98_672_758,
        },
        "dm3": {
            "50": 130_428_560,
            "75": 135_004_462,
            "100": 139_647_232,
            "150": 144_307_808,
            "200": 148_524_010,
        },
        "dm6": {
            "50": 125_464_728,
            "75": 127_324_632,
            "100": 129_789_873,
            "150": 129_941_135,
            "200": 132_509_163,
        },
        "GRCh37": {
            "50": 2_685_511_504,
            "75": 2_736_124_973,
            "100": 2_776_919_808,
            "150": 2_827_437_033,
            "200": 2_855_464_000,
        },
        "GRCh38": {
            "50": 2_701_495_761,
            "75": 2_747_877_777,
            "100": 2_805_636_331,
            "150": 2_862_010_578,
            "200": 2_887_553_303,
        },
        "GRCm37": {
            "50": 2_304_947_926,
            "75": 2_404_646_224,
            "100": 2_462_481_010,
            "150": 2_489_384_235,
            "200": 2_513_019_276,
        },
        "GRCm38": {
            "50": 2_308_125_349,
            "75": 2_407_883_318,
            "100": 2_467_481_108,
            "150": 2_494_787_188,
            "200": 2_520_869_189,
        },
    }
    read_length: int = config.get("library", {}).get("read_length", 100)
    effective_genome_size = default_effective_genome_size.get(build).get(
        read_length, 2_805_636_331
    )


default_fastq_screen_genomes: List[str] = [
    "Adapters/",
    "Arabidopsis/",
    "Drosophila/",
    "E_coli/",
    "Human/",
    "Lambda/",
    "Mitochondria/",
    "Mouse/",
    "PhiX/",
    "Rat/",
    "Vectors/",
    "Worm/",
    "Yeast/",
    "rRNA/",
]


##########################
### Protocol functions ###
##########################


def protocol_is_atac(protocol: str = protocol) -> bool:
    """
    Return `True` if protocol is Atac-seq
    """
    return protocol.startswith("atac")


def protocol_is_medip(protocol: str = protocol) -> bool:
    """
    Return `True` if protocol is MeDIP-seq
    """
    return protocol.startswith("medip")


def protocol_is_chip(protocol: str = protocol) -> bool:
    """
    Return `True` if protocol is ChIP-seq
    """
    return protocol.startswith("chip")


def protocol_is_cutntag(protocol: str = protocol) -> bool:
    """
    Return `True` if protocol is Cut & Tag
    """
    return protocol == "cut&tag"


def protocol_is_cutnrun(protocol: str = protocol) -> bool:
    """
    Return `True` if protocol is Cut & Run
    """
    return protocol == "cut&run"


def protocol_is_ogseq(protocol: str = protocol) -> bool:
    """
    Return `True` if protocol is 8-OxoG-seq
    """
    return ("og" in protocol) or ("ox" in protocol)


#############################################
### Experimental design related functions ###
#############################################


def is_paired(sample: str, design: pandas.DataFrame = design) -> Optional[str]:
    """
    Return path to downstream read file if a sample is pair-ended, else return `None`
    """
    if "Downstream_file" in design.keys():
        down: str = design["Downstream_file"].loc[sample]
        if (down is not None) and (down != ""):
            return down


def has_fragment_size(
    sample: str, design: pandas.DataFrame = design, sample_is_paired: bool = False
) -> Union[str, int]:
    """
    Return the expected fragement size of a single-ended sample.

    Raise warning if this information is missing
    """
    if "Fragment_size" in design.keys():
        fragment_size: Union[str, int] = design["Fragment_size"].loc[sample]
        if (fragment_size is not None) and (fragment_size != ""):
            return fragment_size
        elif not sample_is_paired:
            raise Warning(f"Single-ended sample `{sample}` has no fragment size")


def has_input(sample: str, design: pandas.DataFrame = design) -> Optional[str]:
    """
    Return input sample id if there is an input. Else return `None`
    """
    if "Input_id" in design.keys():
        input_id: str = design["Input_id"].loc[sample]
        if (input_id is not None) and (input_id != ""):
            return input_id


############################################
### Differential Peak Coverage functions ###
############################################


def get_samples_per_condition(
    wildcards, design: pandas.DataFrame = design
) -> List[str]:
    """
    For a `wildcards.condition` return the list of samples id
    belonging to that condition
    """
    return (
        design[design.eq(wildcards.condition).any(axis=1)]
        .index.drop_duplicates()
        .tolist()
    )


def get_input_per_condition(wildcards, design: pandas.DataFrame = design) -> List[str]:
    """
    For a `wildcards.condition` return the list of input id
    corresponding to the samples belonging to that condition
    """
    input_samples: List[str] = []
    for sample in get_samples_per_condition(wildcards, design):
        input_id: str = has_input(sample)
        if input_id:
            input_samples.append(input_id)

    return list(set(input_samples))


def get_sample_genome(
    wildcards,
    design: pandas.DataFrame = design,
    species: str = species,
    build: str = build,
    release: Union[str, int] = release,
) -> Tuple[Union[str, int]]:
    """
    For a `wildcards.sample` return the corresponding dict:
    {genome, build, release}
    """
    if "Species" in design.columns:
        species = design["Species"].loc[wildcards.sample] or species

    if "Release" in design.columns:
        release = design["Release"].loc[wildcards.sample] or release

    if "Build" in design.columns:
        build = design["Build"].loc[wildcards.sample] or build

    return species, build, release


###################
### IO function ###
###################

################
### Trimming ###
################


def get_fastp_input(
    wildcards, design: pandas.DataFrame = design
) -> Dict[str, List[str]]:
    """
    Return the list of Fastp input files
    """
    fastp_input = {
        "sample": [design["Upstream_file"].loc[wildcards.sample]],
    }
    down = is_paired(wildcards.sample, design=design)
    if down:
        fastp_input["sample"].append(down)

    return fastp_input


def get_fastp_output_html(design: pandas.DataFrame = design) -> List[str]:
    """
    Return the list of fastp reports
    """
    html_list = []
    for sample in design.index:
        down = is_paired(sample)
        if down:
            html_list.append(f"data_output/QC/fastp/{sample}.pe.html")
        else:
            html_list.append(f"data_output/QC/fastp/{sample}.se.html")

    return html_list


def get_multiqc_trimming_input(
    wildcards, protocol: str = protocol, design: pandas.DataFrame = design
) -> List[str]:
    """
    Return the expected list of input files for multiqc right after trimming
    """
    multiqc_trimming_input = []
    for sample in design.index:
        down = is_paired(sample)
        if down:
            multiqc_trimming_input.append(f"data_output/QC/fastp/{sample}.pe.html")

            multiqc_trimming_input.append(f"fastq_screen/{sample}.1.fastq_screen.txt")

            multiqc_trimming_input.append(f"fastq_screen/{sample}.2.fastq_screen.txt")

            multiqc_trimming_input.append(f"fastp/report/pe/{sample}.fastp.json")

        else:
            multiqc_trimming_input.append(f"data_output/qc/fastp/{sample}.se.html")

            multiqc_trimming_input.append(f"fastq_screen/{sample}.fastq_screen.txt")

            multiqc_trimming_input.append(f"fastp/report/se/{sample}.fastp.json")

    return multiqc_trimming_input


###############
### Mapping ###
###############


def get_bowtie2_align_input(
    wildcards,
    design: pandas.DataFrame = design,
    bowtie2_index: str = bowtie2_index_path,
) -> Dict[str, List[str]]:
    """
    Return the list of bowtie2 align input
    """
    bowtie2_align_input = {"idx": bowtie2_index_path}
    down = is_paired(wildcards.sample, design=design)
    if down:
        bowtie2_align_input["sample"] = expand(
            "fastp/trimmed/pe/{sample}.{stream}.fastq",
            stream=["1", "2"],
            sample=[wildcards.sample],
        )
    else:
        bowtie2_align_input["sample"] = expand(
            "fastp/trimmed/se/{sample}.fastq", sample=[wildcards.sample]
        )

    return bowtie2_align_input


def get_fastq_screen_input(
    wildcards,
    config: Dict[str, Any] = config,
    design: pandas.DataFrame = design,
    default_fastq_screen_genomes: List[str] = default_fastq_screen_genomes,
) -> List[str]:
    """
    Return the list of expected input files for fastq_screen

    First file in the input file list MUST be the fastq file of interest.
    Second file in the input file list MUST be the configuration file.
    """
    fastq_screen_input = []

    # Let the fastq file be the first file in the input file list
    down = is_paired(wildcards.sample, design=design)
    if down:
        fastq_screen_input.append(
            f"fastp/trimmed/pe/{wildcards.sample}.{wildcards.stream}.fastq"
        )
    else:
        fastq_screen_input.append(f"fastp/trimmed/se/{wildcards.sample}.fastq")

    # Let the config file be the second input
    fq_conf = config.get("reference", {}).get(
        "fastq_screen_config", "../../../config/fastq_screen.conf"
    )
    fastq_screen_input.append(fq_conf)

    # Optional other input files
    if config.get("reference", {}).get("download_fastq_screen_indexes", False):
        for fq_genome in default_fastq_screen_genomes:
            fastq_screen_input.append(
                f"reference/fastq_screen/index/{fq_genome}",
            )

    return fastq_screen_input


def get_samtools_stats_input(wildcards, protocol: str = protocol) -> Dict[str, str]:
    """
    Return expected input files for Samtools stats
    """
    if str(wildcards.step) == "raw":
        return {
            "bam": f"bowtie2/align/{wildcards.sample}.bam",
            "bai": f"bowtie2/align/{wildcards.sample}.bam.bai",
        }

    if protocol_is_ogseq(protocol):
        return {
            "bam": f"deeptools/corrected/{wildcards.sample}.bam",
            "bai": f"deeptools/corrected/{wildcards.sample}.bam.bai",
        }

    # if protocol_is_atac(protocol):
    #     return {
    #         "bam": f"deeptools/alignment_sieve/{wildcards.sample}.bam",
    #         "bai": f"deeptools/alignment_sieve/{wildcards.sample}.bam.bai",
    #     }

    return {
        "bam": f"sambamba/markdup/{wildcards.sample}.bam",
        "bai": f"sambamba/markdup/{wildcards.sample}.bam.bai",
    }


def get_multiqc_mapping_input(
    wildcards, protocol: str = protocol, design: pandas.DataFrame = design
) -> List[str]:
    """
    Return the expected list of input files for multiqc right after mapping
    """
    multiqc_mapping_input = get_multiqc_trimming_input(wildcards, protocol, design)

    picard_files = [
        ".alignment_summary_metrics",
        ".insert_size_metrics",
        ".insert_size_histogram.pdf",
        ".quality_distribution_metrics",
        ".quality_distribution.pdf",
        ".gc_bias.detail_metrics",
        ".gc_bias.summary_metrics",
        ".gc_bias.pdf",
    ]

    for sample in design.index:
        multiqc_mapping_input.append(f"sambamba/markdup/{sample}.bam")

        for step in [".raw", ""]:
            multiqc_mapping_input.append(
                f"samtools/stats/{sample}{step}.txt",
            )

        for picard_file in picard_files:
            multiqc_mapping_input.append(
                f"picard/collectmultiplemetrics/stats/{sample}{picard_file}"
            )

    multiqc_mapping_input.append("deeptools/plot_fingerprint/raw_counts.tab")

    return multiqc_mapping_input


################
### Coverage ###
################


def get_deeptools_bamcoverage_input(
    wildcards, protocol: str = protocol, blacklist_path: str = blacklist_path
) -> Dict[str, str]:
    """
    Return the expected list of DeepTools input file list
    """
    deeptools_bamcoverage_input = {"blacklist": blacklist_path}
    if protocol_is_atac(protocol):
        deeptools_bamcoverage_input["bam"] = "deeptools/alignment_sieve/{sample}.bam"
        deeptools_bamcoverage_input[
            "bai"
        ] = "deeptools/alignment_sieve/{sample}.bam.bai"
    else:
        deeptools_bamcoverage_input["bam"] = "sambamba/markdup/{sample}.bam"
        deeptools_bamcoverage_input["bai"] = "sambamba/markdup/{sample}.bam.bai"

    return deeptools_bamcoverage_input


def get_deeptools_plotfingerprint_input(
    wildcards, protocol: str = protocol, design: pandas.DataFrame = design
) -> Dict[str, List[str]]:
    """
    Return the list of expected input files for deeptools plot fingerprint
    """
    bam_prefix = (
        "deeptools/alignment_sieve"
        if protocol_is_atac(protocol)
        else "sambamba/markdup"
    )

    deeptools_plotfingerprint_input = {"bam_files": [], "bam_idx": []}
    for sample in design.index:
        deeptools_plotfingerprint_input["bam_files"].append(
            f"{bam_prefix}/{sample}.bam"
        )

        deeptools_plotfingerprint_input["bam_idx"].append(
            f"{bam_prefix}/{sample}.bam.bai"
        )

    return deeptools_plotfingerprint_input


def get_medips_params_extra(
    wildcards, design: pandas.DataFrame = design, build: str = build
) -> List[str]:
    """
    Return expected optional parameters for MEDIPS::MEDIPS.createSet(...)
    """
    medips_params_extra = []
    medips_base = config.get("medips", {}).get("medips_createset_extra", "")

    if build.lower() == "grch38":
        medips_base += ", BSgenome = BSgenome.Hsapiens.UCSC.hg38"
    elif build.lower() == "grcm38":
        medips_base += ", BSgenome = BSgenome.Mmusculus.UCSC.mm10"
    else:
        raise NotImplementedError(f"{build} genome not implemented for MeDIP-Seq")

    for sample in get_samples_per_condition(wildcards, design):
        down = is_paired(wilcards.sample)
        fragment_size = has_fragment_size(wilcards.sample, sample_is_paired=bool(down))
        if down:
            medips_params_extra.append(f"{medips_base}, paired = TRUE")
        elif fragment_size:
            # Case the bam is single-ended and properly annotated
            medips_params_extra.append(f"{medips_base}, extend = {fragment_size}")

    return medips_params_extra


def get_medips_meth_coverage_input(
    wildcards: snakemake.utils.Wildcards,
    design: pandas.DataFrame = design,
    config: Dict[str, Any] = config,
) -> Dict[str, Union[str, List[str]]]:
    """
    Return list of input expected by rule medips_meth_coverage
    """
    # Gathering comparison information
    comparisons = config.get("differential_peak_coverage")
    if not comparison:
        raise ValueError(
            "No differential peak coverage information provided in "
            "configuration file. Please fill configuration file."
        )

    # Gathering reference/tested samples
    reference = None
    tested = None
    for models in comparisons:
        if models["model_name"] == wilcards.comparison:
            reference = models["reference"]
            tested = models["tested"]
            break
    else:
        raise ValueError("Could not find a comparison " f"named: {wilcards.comparison}")

    reference_samples = get_samples_per_condition(
        snakemake.io.Wildcards(fromdict={"condition": reference})
    )
    tested_samples = get_samples_per_condition(
        snakemake.io.Wildcards(fromdict={"condition": tested})
    )

    # Building output dictionary
    medips_meth_coverage_input = {
        "mset1": expand("sambamba/markdup/{sample}.bam", sample=tested_samples),
        "mset2": expand("sambamba/markdup/{sample}.bam", sample=reference_samples),
        "cset": f"medips/coupling/{tested}.RDS",
    }

    # Looking for input samples, if any
    # ... in reference condition,
    reference_input_samples = get_input_per_condition(
        snakemake.io.Wildcards(fromdict={"condition": reference})
    )
    if len(reference_input_samples) > 0:
        medips_meth_coverage_input["iset1"] = expand(
            "sambamba/markdup/{sample}.bam", sample=reference_input_samples
        )

    # ... and in tested condition.
    tested_input_samples = get_input_per_condition(
        snakemake.io.Wildcards(fromdict={"condition": tested})
    )
    if len(tested_input_samples) > 0:
        medips_meth_coverage_input["iset1"] = expand(
            "sambamba/markdup/{sample}.bam", sample=tested_input_samples
        )

    return medips_meth_coverage_input


####################
### Peak Calling ###
####################


def get_macs2_callpeak_input(wildcards, protocol: str = "chip-seq") -> Dict[str, str]:
    """
    Return expected list of input files for Macs2 callpeak
    """
    macs2_callpeak_input = {}
    if protocol_is_atac(protocol):
        macs2_callpeak_input[
            "teatment"
        ] = f"deeptools/alignment_sieve/{wildcards.sample}.bam"
    else:
        macs2_callpeak_input["teatment"] = f"sambamba/markdup/{wildcards.sample}.bam"

    input_id = has_input(wilcards.sample, design=design)
    if (input_id is not None) and (input_id != ""):
        if protocol_is_atac(wilcards.sample):
            macs2_callpeak_input[
                "control"
            ] = f"deeptools/alignment_sieve/{input_id}.bam"
        else:
            macs2_callpeak_input["control"] = f"sambamba/markdup/{input_id}.bam"

    return macs2_callpeak_input


def get_macs2_params(
    wildcards,
    effective_genome_size: int = effective_genome_size,
    design: pandas.DataFrame = design,
) -> str:
    """
    Return expected parameters for Macs2 callpeak
    """
    extra = f" --gsize {effective_genome_size} "
    down = is_paired(wildcards.sample, design=design)
    if down:
        extra += " --format BAMPE "
    else:
        extra += " --format BAM "

        fragment_size = has_fragment_size(wilcards.sample, sample_is_paired=False)
        if fragment_size:
            extra += f" --nomodel --extsize {fs} "
        else:
            raise ValueError(
                "Single-ended reads should have a "
                "`Fragment_size` associated in the design file."
                f"{wilcards.sample} has none."
            )

    return extra


#############################
### Wildcards constraints ###
#############################


wildcard_constraints:
    sample=r"|".join(design.index),
    protocol=protocol,
    release=str(release),
    build=str(build),
    species=str(species),
    peaktype=r"|".join(["narrow", "broad", "gapped"]),
    step=r"|".join([".raw", ""]),
    command=r"|".join(["scale-region", "reference-point"]),
    tool=r"|".join(["bowtie2", "sambamba", "deeptools"]),
    subcommand=r"|".join(["align", "markdup", "view", "alignment_sieve"]),


########################
### Main Target rule ###
########################


def targets(
    config: Dict[str, Any] = config,
    design: pandas.DataFrame = design,
    protocol: str = protocol,
    genome_fasta_path: str = genome_fasta_path,
    genome_annotation_path: str = genome_annotation_path,
    bowtie2_index_path: str = bowtie2_index_path,
    default_fastq_screen_genomes: List[str] = default_fastq_screen_genomes,
):
    """
    Return the list of expected output files, depending on the
    choices made by user in configuration file at: `<root>/config/config.yaml`
    """
    expected_targets = {}
    steps = config.get("steps", {})

    if steps.get("install", False):
        expected_targets["gtf"] = genome_annotation_path
        expected_targets["genome_fasta"] = genome_fasta_path
        expected_targets["bowtie2_index"] = bowtie2_index_path
        if steps.get("download_fastq_screen_indexes", False):
            expected_targets["fastq_screen_indexes"] = expand(
                "reference/fastq_screen/index/{fq_genome}",
                fq_genome=default_fastq_screen_genomes,
            )

    if steps.get("trimming", False):
        expected_targets["fastp"] = get_fastp_output_html(design=design)
        expected_targets["multiqc_trim"] = "data_output/QC/Trimming.QC.html"

    if steps.get("mapping", False):
        expected_targets["mapping"] = expand(
            "data_output/CRAM/{sample}.cram", sample=design.index
        )
        expected_targets["multiqc_map"] = ("data_output/QC/Mapping.QC.html",)

    if steps.get("coverage", False):
        expected_targets["bam_coverage"] = expand(
            "data_output/Coverage/{sample}.bw", sample=design.index
        )

    if steps.get("calling", False):
        if config.get("macs2", {}).get("broad", False):
            expected_targets["macs2_broad"] = expand(
                "data_output/Peak_Calling/macs2/{sample}_broad_peaks.xls",
                sample=design.index,
            )

        if config.get("macs2", {}).get("narrow", False):
            expected_targets["macs2_broad"] = expand(
                "data_output/Peak_Calling/macs2/{sample}_narrow_peaks.xls",
                sample=design.index,
            )

    if steps.get("diff_cov", False):
        raise NotImplementedError("Differential coverage analysis not yet implemented")

    if steps.get("motives", False):
        raise NotImplementedError("Mitives analysis not yet implemented")

    print(expected_targets)
    return expected_targets
