# Used for faster/easier iterations in python
import itertools

# Used to define logging behavior
import logging

# Used to check OS environment
import os

# Used to deal with (large) tables
import pandas

# Main pipeline handler
import snakemake

# Display non-blocking warnings
import warnings

# Validate user IO
from snakemake.utils import validate

# Type hinting for future devs
from typing import Any, Dict, List, Optional, Tuple, Union


##########################
### Load configuration ###
##########################

# Load config if and only if user did not provide any
if config == {}:

    configfile: "config/config.yaml"


# Check config required fields, and variable types
validate(config, schema="../schemas/config.schema.yaml")


#########################
### Logging behaviour ###
#########################

logging.basicConfig(filename="epicure_pipeline.log", filemode="w", level=logging.DEBUG)


########################
### Load design file ###
########################

# Load design
design_path: str = config.get("designs", "config/design.tsv")
design: pandas.DataFrame = pandas.read_csv(
    filepath_or_buffer=design_path,
    sep="\t",
    header=0,
    index_col=0,
)

# Check design required fields, and variable types
validate(design, schema="../schemas/design.schema.yaml")


################################
### Main Snakemake variables ###
################################


# Main report content
report: "../report/workflow.rst"


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"


# This is the main version of all snakemake-wrappers
snakemake_wrappers_version: str = "v2.0.0"


########################
### Global variables ###
########################


# This let all temporary file end-up in the same directory
tmp: str = os.environ.get("TMP", "temp")


# Mouse has less chromosomes, but this will work anyway
# Does not work for organisms with more than 22 autosomal chromosomes
# Does not work for organisms with no X/Y system for sexual chromosomes
canonical_chromosomes = list(map(str, range(1, 23))) + ["X", "Y"]
canonical_chromosomes += [f"chr{chrom}" for chrom in canonical_chromosomes]


# Protocol. It can be either : `chip-seq`, `atac-seq`, `cut&run`, `cut&tag`, or `medip-seq`
# Protocol is converted to lower case here. No need to do it later.
protocol: str = config.get("protocol", "chip-seq").lower().replace(" ", "")


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
    return protocol in ("cut&tag", "cutntag")


def protocol_is_cutnrun(protocol: str = protocol) -> bool:
    """
    Return `True` if protocol is Cut & Run
    """
    return protocol in ("cut&run", "cutnrun")


def protocol_is_ogseq(protocol: str = protocol) -> bool:
    """
    Return `True` if protocol is 8-OxoG-seq
    """
    return ("og" in protocol) or ("ox" in protocol)


def protocol_is_mnaseseq(protocol: str = protocol) -> bool:
    """
    Return `True` if protocol is MNase-seq
    """
    return protocol.startswith("mn")


def protocol_is_riboseq(protocol: str = protocol) -> bool:
    """
    Return `True` if protocol is Ribo-seq
    """
    return protocol.startswith("ribo")


def protocol_is_groseq(protocol: str = protocol) -> bool:
    """
    Return `True` if protocol is Ribo-seq
    """
    return protocol.startswith("gro")


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


# Xenome genomes index with host = mouse, and target = human
xenome_index: List[str] = config.get("reference", {}).get("xenome_index")
if not xenome_index:
    logging.info(
        "Missing Xenome index path in the file `config.yaml`. "
        "A new one will be created."
    )
    xenome_index = multiext(
        "reference/xenome/index/pdx-both",
        ".header",
        ".kmers-d0",
        ".kmers-d1",
        ".kmers.header",
        ".kmers.high-bits",
        ".kmers.low-bits.lwr",
        ".kmers.low-bits.upr",
        ".lhs-bits",
        ".rhs-bits",
    )


# Genome blacklisted regions
blacklist_path: str = config.get("reference", {}).get("blacklist")
if not blacklist_path:
    blacklist_path = f"reference/blacklist/{species}.{build}.{release}.merged.bed.gz"


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


# FastQ-Screen indexes
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

# List of expected ChipSeeker plots
chipseeker_plot_list: List[str] = [
    "UpsetVenn",
    "Feature_Distribution",
    "Distance_to_TSS",
    "Genome_Coverage",
    "Gene_Body_Coverage",
]

# List of expected deeptools plots
deeptools_plot_type: List[str] = ["heatmap", "scatterplot"]

# List of expected peaks with Macs2
peaktype_list: List[str] = []
macs2_params = config.get("macs2", {})
if macs2_params is None:
    macs2_params = {}
if macs2_params.get("broad", False):
    peaktype_list.append("broad")
if macs2_params.get("narrow", False):
    peaktype_list.append("narrow")

# List of expected peaks with Seacr
seacr_mode_list: List[str] = []
seacr_params = config.get("seacr", {})
if (seacr_params is None) and (
    not (protocol_is_cutnrun(protocol) or protocol_is_cutntag(protocol))
):
    seacr_params = {}
if seacr_params.get("relaxed", False):
    seacr_mode_list.append("relaxed")
if seacr_params.get("stringent", False):
    seacr_mode_list.append("stringent")

# List of expected peaks with ALL peak-callers
peak_type_and_mode_list: List[str] = peaktype_list + seacr_mode_list

# Optional expected list of motifs to plot
motifs_list: Optional[List[str]] = config.get("finger_prints", {}).get("motifs")

# List of deeptools transformations to do
command_list: List[str] = ["reference-point"]  # ["scale-region", "reference-point"]


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
        if (
            (fragment_size is not None)
            and (fragment_size != "")
            and (not pandas.isna(fragment_size))
        ):
            return fragment_size
        elif not sample_is_paired:
            warnings.warn(
                f"Single-ended sample `{sample}` has no fragment size: {fragment_size}"
            )


def has_input(sample: str, design: pandas.DataFrame = design) -> Optional[str]:
    """
    Return input sample id if there is an input. Else return `None`
    """
    if "Input" in design.keys():
        input_id: str = design["Input"].loc[sample]
        if (input_id is not None) and (input_id != "") and not (pandas.isna(input_id)):
            return input_id


############################################
### Differential Peak Coverage functions ###
############################################


def get_paired_samples(design: pandas.DataFrame = design) -> Optional[List[str]]:
    """
    Return the list of paired sample, if there are any. Else return None.
    """
    result = []
    for sample in design.index:
        if is_paired(sample=sample, design=design):
            result.append(sample)

    if len(result) < 1:
        return None

    return result


def get_samples_per_level(
    wildcards: snakemake.io.Wildcards, design: pandas.DataFrame = design
) -> List[str]:
    """
    For a `wildcards.level` return the list of samples id
    belonging to that level
    """
    if "level" in wildcards.keys():
        return (
            design[design.eq(wildcards.level).any(axis=1)]
            .index.drop_duplicates()
            .tolist()
        )
    if "factor_level" in wildcards.keys():
        return (
            design[design.eq(wildcards.factor_level).any(axis=1)]
            .index.drop_duplicates()
            .tolist()
        )
    return (
        design[design.eq(wildcards.model_name).any(axis=1)]
        .index.drop_duplicates()
        .tolist()
    )


def get_input_per_level(
    wildcards: snakemake.io.Wildcards, design: pandas.DataFrame = design
) -> List[str]:
    """
    For a `wildcards.level` return the list of input id
    corresponding to the samples belonging to that level
    """
    input_samples: List[str] = []
    for sample in get_samples_per_level(wildcards=wildcards, design=design):
        input_id: str = has_input(sample=sample, design=design)
        if input_id:
            input_samples.append(input_id)

    return list(set(input_samples))


def get_level_from_model_name(
    model_name: str, config: Dict[str, Any] = config
) -> Dict[str, Optional[str]]:
    """
    Return a dictionary containing {ref: reference-level-name, test: tested-level-name}
    from a comparison name.
    """
    if "differential_peak_coverage" in config.keys():
        for idx, comparison in enumerate(config["differential_peak_coverage"]):
            if comparison["model_name"] == model_name:
                return {
                    "ref": comparison["reference"],
                    "test": comparison["tested"],
                }
        raise KeyError(
            f"{model_name} is not a comparison name available in "
            "configfile. Available names are: "
            f"{[c['model_name'] for c in config['differential_peak_coverage']]}"
        )
    else:
        return {"ref": None, "test": None}


def get_sample_list_from_model_name(
    model_name: str,
    signal: str,
    design: pandas.DataFrame = design,
    config: Dict[str, Any] = config,
) -> List[str]:
    """
    Return list of samples belonging to a given comparison
    """
    level_dict: Dict[str, Optional[str]] = get_level_from_model_name(
        model_name=model_name, config=config
    )
    sample_list: List[str] = get_samples_per_level(
        wildcards=snakemake.io.Wildcards(fromdict={"level": level_dict["ref"]}),
        design=design,
    )
    sample_list += get_samples_per_level(
        wildcards=snakemake.io.Wildcards(fromdict={"level": level_dict["test"]}),
        design=design,
    )

    if len(sample_list) == 0:
        raise ValueError(f"No sample found in comparison: {model_name}")

    if str(signal).lower() == "input":
        if any(has_input(sample=sample, design=design) for sample in sample_list):
            sample_list = get_input_per_level(
                wildcards=snakemake.io.Wildcards(fromdict={"level": level_dict["ref"]}),
                design=design,
            )
            sample_list += get_input_per_level(
                wildcards=snakemake.io.Wildcards(
                    fromdict={"level": level_dict["test"]}
                ),
                design=design,
            )
        elif not protocol_is_atac(protocol):
            logging.warning(
                "No sample have corresponding input. "
                "Filtering and normalisation cannot be done using Input files."
            )

    return sample_list


def get_model_names(config: Dict[str, Any] = config) -> List[str]:
    """
    Return complete list of comparisons
    """
    model_names: List[str] = []
    if "differential_peak_coverage" in config.keys():
        for comparison in config["differential_peak_coverage"]:
            model_names.append(comparison["model_name"])

    return model_names


def get_sample_genome(
    wildcards: snakemake.io.Wildcards,
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
        species: str = design["Species"].loc[wildcards.sample] or species

    if "Release" in design.columns:
        release: Union[str, int] = design["Release"].loc[wildcards.sample] or release

    if "Build" in design.columns:
        build: str = design["Build"].loc[wildcards.sample] or build

    return species, build, release


def get_tested_sample_list(design: pandas.DataFrame = design) -> List[str]:
    """
    Return list of non-input sample ids
    """
    sample_list: List[str] = design.index.tolist()
    if "Input" in design.columns:
        input_list = [i for i in design.Input if i != ""]
        sample_list = list(set(sample_list) - set(input_list))

    return sample_list


def get_bam_prefix(wildcards: snakemake.io.Wildcards, protocol: str = protocol) -> str:
    """
    Return the best refined file after all corrections
    """
    shift: str = str(
        config.get("mapping", {}).get("deeptools", {}).get("alignment_sieve_extra", "")
    )
    bam_prefix: str = "sambamba/markdup"
    # For raw-bam quality control
    if "step" in wildcards:
        if str(wildcards.step) == "raw":
            bam_prefix = "bowtie2/align"
    # Atac-seq samples should be shifted
    elif protocol_is_atac(protocol) or shift != "":
        bam_prefix = "deeptools/sorted_sieve"
    # OxiDIP-seq requires GC correction and filters
    elif protocol_is_ogseq(protocol):
        bam_prefix = "deeptools/corrected"

    return bam_prefix


def get_peak_file(
    wildcards: snakemake.io.Wildcards,
    peaktype_list: List[str] = peaktype_list,
    seacr_mode_list: List[str] = seacr_mode_list,
) -> str:
    """
    Return the correct bed file after peak-calling, given a peak-calling mode
    """
    if str(wildcards.peaktype) in seacr_mode_list:
        return "seacr/{peaktype}/{sample}.bed"

    if str(wildcards.peaktype) in peaktype_list:
        return "macs2/callpeak_{peaktype}/{sample}_peaks.{peaktype}Peak.bed"

    raise ValueError(
        f"Wildcards.peaktype should be in {peaktype_list} "
        f"or in {seacr_mode_list}. Got: `{wildcards.peaktype}`"
    )


###################
### IO function ###
###################

################
### Trimming ###
################


def get_rename_input(
    wildcards: snakemake.io.Wildcards, design: pandas.DataFrame = design
) -> Dict[str, List[str]]:
    """
    Return the list of Fastp input files
    """
    if "stream" in wildcards.keys():
        if int(wildcards.stream) == 1:
            return design["Upstream_file"].loc[wildcards.sample]
        pair = is_paired(sample=wildcards.sample, design=design)
        return pair
    return design["Upstream_file"].loc[wildcards.sample]


def get_fastp_input(
    wildcards: snakemake.io.Wildcards, design: pandas.DataFrame = design
) -> Dict[str, List[str]]:
    """
    Return the list of Fastp input files
    """
    if is_paired(sample=wildcards.sample, design=design):
        return {
            "sample": [
                f"data_input/{wildcards.sample}.{stream}.fq.gz" for stream in ["1", "2"]
            ]
        }
    return {"sample": [f"data_input/{wildcards.sample}.fq.gz"]}


def get_fastp_params(
    wildcards: snakemake.io.Wildcards, protocol: str = protocol
) -> str:
    """
    Return fastp parameters
    """
    extra: str = f"--report_title 'Fastp report for sample {wildcards.sample}' "
    user_defined: Optional[str] = config.get("trimming", {})
    if user_defined:
        extra += user_defined.get("fastp_extra")

    if protocol_is_ogseq(protocol):
        extra += " --poly_g_min_len 25 "

    return extra


def get_fastp_output_html(design: pandas.DataFrame = design) -> List[str]:
    """
    Return the list of fastp reports
    """
    html_list: List[str] = []
    for sample in design.index:
        down = is_paired(sample=sample, design=design)
        if down:
            html_list.append(f"data_output/QC/fastp/{sample}.pe.html")
        else:
            html_list.append(f"data_output/QC/fastp/{sample}.se.html")

    return html_list


def get_multiqc_trimming_input(
    wildcards: snakemake.io.Wildcards,
    protocol: str = protocol,
    design: pandas.DataFrame = design,
) -> List[str]:
    """
    Return the expected list of input files for multiqc right after trimming
    """
    multiqc_trimming_input: List[str] = []
    for sample in design.index:
        down: Optional[str] = is_paired(sample=sample, design=design)
        if down:
            multiqc_trimming_input.append(f"data_output/QC/fastp/{sample}.pe.html")

            multiqc_trimming_input.append(f"fastq_screen/{sample}.1.fastq_screen.txt")

            multiqc_trimming_input.append(f"fastq_screen/{sample}.2.fastq_screen.txt")

            multiqc_trimming_input.append(f"fastp/report/pe/{sample}.fastp.json")

        else:
            multiqc_trimming_input.append(f"data_output/QC/fastp/{sample}.se.html")

            multiqc_trimming_input.append(f"fastq_screen/{sample}.fastq_screen.txt")

            multiqc_trimming_input.append(f"fastp/report/se/{sample}.fastp.json")

    return multiqc_trimming_input


###########
### PDX ###
###########


def get_xenome_index_input(
    wildcards: snakemake.io.Wildcards, config: Dict[str, Any] = config
) -> Dict[str, str]:
    """
    Return best input files for xenome input
    """
    if species.lower() != "homo_sapiens":
        raise ValueError("Xenographft from human (graft) into mouse (host) only")

    mouse_fasta = config.get("reference", {}).get(
        "mouse_fasta", "reference/mus_musculus.GRCm38.99.fasta"
    )

    return {
        "human": genome_fasta_path,
        "mouse": mouse_fasta,
    }


def get_xenome_classify_input(
    wildcards: snakemake.io.Wildcards,
    design: pandas.DataFrame = design,
    xenome_index: List[str] = xenome_index,
) -> Dict[str, Union[str, List[str]]]:
    """
    Return best input files for xenome classify
    """
    xenome_classify_input: Dict[str, Union[str, List[str]]] = {
        "index": xenome_index,
    }

    pair = is_paired(sample=wildcards.sample, design=design)
    if pair:
        xenome_classify_input["r1"] = f"fastp/trimmed/pe/{wildcards.sample}.1.fastq"
        xenome_classify_input["r2"] = f"fastp/trimmed/pe/{wildcards.sample}.2.fastq"
    else:
        xenome_classify_input["reads"] = f"fastp/trimmed/se/{wildcards.sample}.fastq"

    return xenome_classify_input


###############
### Mapping ###
###############


def get_sambamba_quality_filter_params(
    wildcards: snakemake.io.Wildcards, design: pandas.DataFrame = design
) -> str:
    """
    Return best parameters for Sambamba quality filter
    """
    extra: str = " --with-header --format 'bam' "
    if is_paired(sample=wildcards.sample, design=design):
        extra += (
            " --filter 'mapping_quality >= 30 and not (unmapped or mate_is_unmapped)' "
        )
    else:
        extra += " --filter 'mapping_quality >= 30' "

    return extra


def get_sambamba_merge_per_factors_level_input_input(
    wildcards: snakemake.io.Wildcards,
    design: pandas.DataFrame = design,
    protocol: str = protocol,
) -> List[str]:
    """
    Return correct list of bam file to concat
    """
    samples: List[str] = get_samples_per_level(wildcards=wildcards, design=design)
    bam_prefix: str = get_bam_prefix(wildcards=wildcards, protocol=protocol)
    return expand("{prefix}/{sample}.bam", sample=samples, prefix=[bam_prefix])


def get_bowtie2_align_input(
    wildcards: snakemake.io.Wildcards,
    design: pandas.DataFrame = design,
    bowtie2_index: str = bowtie2_index_path,
    config: Dict[str, Any] = config,
) -> Dict[str, List[str]]:
    """
    Return the list of bowtie2 align input
    """
    bowtie2_align_input: Dict[str, List[str]] = {"idx": bowtie2_index_path}
    down: Optional[str] = is_paired(sample=wildcards.sample, design=design)
    pdx: bool = config.get("reference", {}).get("pdx", False)
    if down:
        if pdx:
            bowtie2_align_input["sample"] = expand(
                "xenome/classify/{sample}_graft_{stream}.fastq",
                stream=["1", "2"],
                sample=[wildcards.sample],
            )
        else:
            bowtie2_align_input["sample"] = expand(
                "fastp/trimmed/pe/{sample}.{stream}.fastq",
                stream=["1", "2"],
                sample=[wildcards.sample],
            )
    elif pdx:
        bowtie2_align_input["sample"] = expand(
            "xenome/classify/{sample}_graft.fastq", sample=wildcards.sample
        )
    else:
        bowtie2_align_input["sample"] = expand(
            "fastp/trimmed/se/{sample}.fastq", sample=wildcards.sample
        )

    return bowtie2_align_input


def get_fastq_screen_input(
    wildcards: snakemake.io.Wildcards,
    config: Dict[str, Any] = config,
    design: pandas.DataFrame = design,
    default_fastq_screen_genomes: List[str] = default_fastq_screen_genomes,
) -> List[str]:
    """
    Return the list of expected input files for fastq_screen

    First file in the input file list MUST be the fastq file of interest.
    Second file in the input file list MUST be the configuration file.
    """
    fastq_screen_input: List[str] = []

    # Let the fastq file be the first file in the input file list
    down: Optional[str] = is_paired(sample=wildcards.sample, design=design)
    if down:
        fastq_screen_input.append(
            f"fastp/trimmed/pe/{wildcards.sample}.{wildcards.stream}.fastq"
        )
    else:
        fastq_screen_input.append(f"fastp/trimmed/se/{wildcards.sample}.fastq")

    # Let the config file be the second input
    fq_conf: str = config.get("reference", {}).get(
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


def get_samtools_stats_input(
    wildcards: snakemake.io.Wildcards,
    protocol: str = protocol,
    config: Dict[str, Any] = config,
) -> Dict[str, str]:
    """
    Return expected input files for Samtools stats
    """

    bam_prefix: str = get_bam_prefix(wildcards=wildcards, protocol=protocol)

    samtools_stats_input: Dict[str, str] = {
        "bam": f"{bam_prefix}/{wildcards.sample}.bam",
        "bai": f"{bam_prefix}/{wildcards.sample}.bam.bai",
    }
    roi_bed: str = config.get("reference", {}).get("roi_bed")
    if roi_bed:
        samtools_stats_input["bed"] = roi_bed

    return samtools_stats_input


def get_sambamba_markdup_params(
    wildcards: snakemake.io.Wildcards, protocol: str = protocol
) -> str:
    """
    Return best sambamba markdup params for a given protocol
    """
    extra = " --overflow-list-size 600000 "

    # Duplicates should not be removed in Cut&Tag sequencings
    if not protocol_is_cutntag(protocol) and not protocol_is_cutnrun(protocol):
        extra += " --remove-duplicates "

    return extra


def get_multiqc_mapping_input(
    wildcards: snakemake.io.Wildcards,
    protocol: str = protocol,
    design: pandas.DataFrame = design,
) -> List[str]:
    """
    Return the expected list of input files for multiqc right after mapping
    """
    multiqc_mapping_input: List[str] = get_multiqc_trimming_input(
        wildcards=wildcards, protocol=protocol, design=design
    )

    picard_files: List[str] = [
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


def get_medips_import_sample_bam_input(
    wildcards: snakemake.io.Wildcards, design: pandas.DataFrame = design
) -> Dict[str, List[str]]:
    """
    Return expected lists of input for medips
    """
    samples: List[str] = get_samples_per_level(wildcards=wildcards, design=design)

    return {
        "bam": expand("sambamba/markdup/{sample}.bam", sample=samples),
        "bai": expand("sambamba/markdup/{sample}.bam.bai", sample=samples),
    }


def get_deeptools_estimate_gc_bias_params(
    wildcards: snakemake.io.Wildcards, design: pandas.DataFrame = design
) -> str:
    """
    Return the best parameters for the provided sample to deeptools computeGCBias
    """
    extra: str = " --plotFileFormat png "

    # Single-ended samples requires fragment size (fs)
    if not is_paired(sample=wildcards.sample, design=design):
        fs: Union[str, int] = has_fragment_size(
            sample=wildcards.sample, design=design, sample_is_paired=False
        )
        extra += f" --fragmentLength {fs} "

    return extra


def get_deeptools_alignment_sieve_params(
    wildcards: snakemake.io.Wildcards,
    protocol: str = protocol,
    config: Dict[str, Any] = config,
) -> str:
    """
    Return best parameters for deeptools alignment sieve
    """
    extra: str = (
        config.get("mapping", {}).get("deeptools", {}).get("alignment_sieve_extra", "")
    )
    if extra == "" and protocol_is_atac(protocol):
        extra = " --ATACshift "

    return extra


def get_bedtools_filter_roi_input(
    wildcards: snakemake.io.Wildcards, config: Dict[str, Any] = config
) -> Dict[str, str]:
    """
    If exists, return path to the regions of interest
    """
    bedtools_filter_roi_input: Dict[str, str] = {
        "left": f"samtools/view/{wildcards.sample}.bam",
        "genome": "reference/genome_contigs.txt",
    }

    roi_bed: Optional[str] = config.get("reference", {}).get("roi_bed")
    if roi_bed:
        bedtools_filter_roi_input["right"] = roi_bed

    return bedtools_filter_roi_input


def get_sambamba_sort_filtered(
    wildcards: snakemake.io.Wildcards, config: Dict[str, Any] = config
) -> List[str]:
    """
    Return best input file for sambamba sort filtered
    """
    if config.get("reference", {}).get("roi_bed"):
        return ["bedtools/filter_roi/{sample}.bam"]
    return ["samtools/view/{sample}.bam"]


def get_rsamtools_bam_file(
    wildcards: snakemake.io.Wildcards, protocol: str = protocol
) -> Dict[str, str]:
    """
    Return correct bam file for RSamtools
    """
    bam_prefix = get_bam_prefix(wildcards=wildcards, protocol=protocol)

    return {
        "bam": f"{bam_prefix}/{wildcards.sample}.bam",
        "bai": f"{bam_prefix}/{wildcards.sample}.bam.bai",
    }


################
### Coverage ###
################


def get_sambamba_merge_input(
    wildcards: snakemake.io.Wildcards,
    protocol: str = protocol,
    design: pandas.DataFrame = design,
) -> List[str]:
    """
    Return required bam files for sambamba merge per factors in each models
    """
    samples: List[str] = get_samples_per_level(wildcards=wildcards, design=design)
    bam_prefix: str = get_bam_prefix(wildcards=wildcards, protocol=protocol)
    return [f"{bam_prefix}/{sample}.bam" for sample in samples]


def get_bedtools_merge_factor_input(
    wildcards: snakemake.io.Wildcards,
    protocol: str = protocol,
    design: pandas.DataFrame = design,
) -> List[str]:
    """
    Return required bam files for bedtools merge per factors in each models
    """
    samples: List[str] = get_samples_per_level(wildcards=wildcards, design=design)
    peak_file: str = get_peak_file(
        wildcards=wildcards,
        peaktype_list=peaktype_list,
        seacr_mode_list=seacr_mode_list,
    )

    return [
        peak_file.format(peaktype=wildcards.peaktype, sample=sample)
        for sample in samples
    ]


def get_deeptools_compute_matrix_per_factor_input(
    wildcards: snakemake.io.Wildcards, config: Dict[str, Any] = config
) -> Dict[str, Union[str, List[str]]]:
    """
    Return list of expected input files in deeptools_compute_matrix_per_factor
    """
    factors: List[str] = list(
        set(
            get_level_from_model_name(
                model_name=wildcards.model_name, config=config
            ).values()
        )
    )
    return {
        "blacklist": blacklist_path,
        "bed": expand(
            "bedtools/merge/{factor_level}.{peaktype}.bed",
            factor_level=factors,
            allow_missing=True,
        ),
        "bigwig": expand(
            "data_output/Coverage/{factor_level}.bw", factor_level=factors
        ),
    }


def get_multiqc_coverage_input(
    wildcards: snakemake.io.Wildcards,
    protocol: str = protocol,
    design: pandas.DataFrame = design,
) -> List[str]:
    """
    Return expected list of input files for MultiQC after coverage analysis
    """
    multiqc_coverage_input: List[str] = get_multiqc_mapping_input(
        wildcards=wildcards, protocol=protocol, design=design
    )

    # Add bampe if any sample is paired
    paired_samples: Optional[List[str]] = get_paired_samples(design=design)
    if paired_samples:
        for sample in paired_samples:
            multiqc_coverage_input.append(f"deeptools/bampe_fs/{sample}_metrics.txt")

    # Add fingerprints
    multiqc_coverage_input.append("deeptools/plot_fingerprint/raw_counts.tab")
    multiqc_coverage_input.append("deeptools/plot_fingerprint/qc_metrics.txt")

    # Add plot correlation
    for plot_type in deeptools_plot_type:
        multiqc_coverage_input.append(f"deeptools/plot_corr/{plot_type}.stats")

    # Add plot pca
    multiqc_coverage_input.append("deeptools/plot_pca.stats")

    # Add plot profile
    for peaktype in ["broad", "narrow"]:
        for command in ["reference-point"]:
            multiqc_coverage_input.append(
                f"deeptools/plot_profile/{peaktype}/{command}.tab"
            )

    return multiqc_coverage_input


def get_deeptools_bamcoverage_input(
    wildcards: snakemake.io.Wildcards,
    protocol: str = protocol,
    blacklist_path: str = blacklist_path,
) -> Dict[str, str]:
    """
    Return the expected list of DeepTools input file list
    """
    bam_prefix = get_bam_prefix(wildcards=wildcards, protocol=protocol)

    return {
        "blacklist": blacklist_path,
        "bam": f"{bam_prefix}/{wildcards.sample}.bam",
        "bam": f"{bam_prefix}/{wildcards.sample}.bam.bai",
    }


def get_deeptools_bamcoverage_params(
    wildcards: snakemake.io.Wildcards, protocol: str = protocol
) -> str:
    """
    Return best parameters for deeptools bamcoverage for a given protocol
    """
    extra = (
        " --normalizeUsing RPKM --binSize 5 "
        "--skipNonCoveredRegions --ignoreForNormalization chrX chrM chrY"
    )
    if not (protocol_is_cutnrun(protocol) or protocol_is_cutntag(protocol)):
        extra += " --ignoreDuplicates "

    if protocol_is_mnaseseq(protocol):
        extra += " --MNase "

    elif protocol_is_riboseq(protocol):
        extra += " --Offset 12 "

    elif protocol_is_groseq(protocol):
        extra += " --Offset 15 "

    return extra


def get_deeptools_plotfingerprint_input(
    wildcards: snakemake.io.Wildcards,
    protocol: str = protocol,
    design: pandas.DataFrame = design,
) -> Dict[str, List[str]]:
    """
    Return the list of expected input files for deeptools plot fingerprint
    """
    bam_prefix = get_bam_prefix(wildcards=wildcards, protocol=protocol)

    deeptools_plotfingerprint_input: Dict[str, List[str]] = {
        "bam_files": [],
        "bam_idx": [],
    }
    for sample in design.index:
        deeptools_plotfingerprint_input["bam_files"].append(
            f"{bam_prefix}/{sample}.bam"
        )

        deeptools_plotfingerprint_input["bam_idx"].append(
            f"{bam_prefix}/{sample}.bam.bai"
        )

    return deeptools_plotfingerprint_input


def get_deeptools_fingerprint_params(
    wildcards: snakemake.io.Wildcards,
    design: pandas.DataFrame = design,
    protocol: str = protocol,
) -> str:
    """
    Return the best parameters for deeptools fingerprint
    """
    labels: str = " ".join(design.index.tolist())
    extra = f" --skipZeros --labels {labels} --plotFileFormat png "
    if not (protocol_is_cutntag(protocol) or protocol_is_cutnrun(protocol)):
        extra += " --ignoreDuplicates "

    return extra


def get_deeptools_bampefragmentsize_params(
    wildcards: snakemake.io.Wildcards, design: pandas.DataFrame = design
) -> str:
    """
    Return the best parameters for deeptools bampe fragment size
    """
    labels: str = " ".join(design.index.tolist())
    return f" --plotFileFormat png --samplesLabel {labels} "


def get_deeptools_plotpca_params(
    wildcards: snakemake.io.Wildcards, design: pandas.DataFrame = design
) -> str:
    """
    Return the best parameters for deeptools plot PCA
    """
    labels: str = " ".join(design.index.tolist())
    return f" --plotFileFormat png --labels {labels} --PCs 1 2"


def get_deeptools_compute_matrix_input(
    wildcards: snakemake.io.Wildcards,
    design: pandas.DataFrame = design,
    peaktype_list: List[str] = peaktype_list,
    seacr_mode_list: List[str] = seacr_mode_list,
) -> str:
    """
    Return the expected list of input files for deeptools compute matrix
    """
    deeptools_compute_matrix_input = {
        "bigwig": expand(
            "data_output/Coverage/{sample}.bw",
            sample=get_tested_sample_list(design=design),
        ),
        "blacklist": blacklist_path,
        "bed": expand(
            get_peak_file(
                wildcards=wildcards,
                peaktype_list=peaktype_list,
                seacr_mode_list=seacr_mode_list,
            ),
            sample=get_tested_sample_list(design=design),
            peaktype=[wildcards.peaktype],
        ),
    }

    return deeptools_compute_matrix_input


def get_deeptools_compute_matrix_params(
    wildcards: snakemake.io.Wildcards,
    design: pandas.DataFrame = design,
    config: Dict[str, Any] = config,
) -> str:
    """
    Return best parameters for deeptools compute matrix
    """
    labels: str = " ".join(get_tested_sample_list(design=design))
    if "model_name" in wildcards.keys():
        labels = " --smartLabels "
    else:
        labels = f" --samplesLabel {labels} "

    extra: str = f" --skipZeros --missingDataAsZero {labels} --binSize 50 "
    if str(wildcards.command) == "reference-point":
        extra += " --referencePoint TSS "
    else:
        extra += " --startLabel TSS "

    return extra


def get_deeptools_multibigwigsummary_params(
    wildcards: snakemake.io.Wildcards, design: pandas.DataFrame = design
) -> str:
    """
    Return best parameters for deeptools multiBigwigSummary
    """
    labels: str = " ".join(get_tested_sample_list(design=design))
    return f" --labels {labels} --chromosomesToSkip chrX chrY chrM"


def get_deeptools_plotcoverage_params(
    wildcards: snakemake.io.Wildcards,
    design: pandas.DataFrame = design,
    protocol: str = protocol,
) -> str:
    """
    Return best parameters for deeptools plot coverage
    """
    labels: str = " ".join(get_tested_sample_list(design=design))
    extra: str = f"--labels {labels} --skipZeros --minMappingQuality 30"
    if not (protocol_is_cutnrun(protocol) or protocol_is_cutntag(protocol)):
        extra += " --ignoreDuplicates "

    return extra


def get_medips_params_extra(
    wildcards: snakemake.io.Wildcards,
    design: pandas.DataFrame = design,
    build: str = build,
) -> List[str]:
    """
    Return expected optional parameters for MEDIPS::MEDIPS.createSet(...)
    """
    medips_params_extra: List[str] = []
    medips_base: str = config.get("medips", {}).get("medips_createset_extra", "")

    if build.lower() == "grch38":
        medips_base += ", BSgenome = BSgenome.Hsapiens.UCSC.hg38"
    elif build.lower() == "grcm38":
        medips_base += ", BSgenome = BSgenome.Mmusculus.UCSC.mm10"
    else:
        raise NotImplementedError(f"{build} genome not implemented for MeDIP-Seq")

    for sample in get_samples_per_level(wildcards=wildcards, design=design):
        down: Optional[str] = is_paired(sample=wildcards.sample, design=design)
        fragment_size = has_fragment_size(
            sample=wildcards.sample, design=design, sample_is_paired=bool(down)
        )
        if down:
            medips_params_extra.append(f"{medips_base}, paired = TRUE")
        elif fragment_size:
            # Case the bam is single-ended and properly annotated
            medips_params_extra.append(f"{medips_base}, extend = {fragment_size}")

    return medips_params_extra


def get_medips_meth_coverage_input(
    wildcards: snakemake.io.Wildcards,
    design: pandas.DataFrame = design,
    config: Dict[str, Any] = config,
) -> Dict[str, Union[str, List[str]]]:
    """
    Return list of input expected by rule medips_meth_coverage
    """
    # Gathering comparison information
    comparisons: Optional[List[Dict[str, str]]] = config.get(
        "differential_peak_coverage"
    )
    if not comparisons:
        raise ValueError(
            "No differential peak coverage information provided in "
            "configuration file. Please fill configuration file."
        )

    # Gathering reference/tested samples
    reference: Optional[str] = None
    tested: Optional[str] = None
    for models in comparisons:
        if models["model_name"] == wildcards.model_name:
            reference = models["reference"]
            tested = models["tested"]
            break
    else:
        raise ValueError(
            "Could not find a comparison " f"named: {wildcards.model_name}"
        )

    reference_samples: List[str] = get_samples_per_level(
        wildcards=snakemake.io.Wildcards(fromdict={"level": reference}), design=design
    )
    tested_samples: List[str] = get_samples_per_level(
        wildcards=snakemake.io.Wildcards(fromdict={"level": tested}), design=design
    )

    # Building output dictionary
    medips_meth_coverage_input: Dict[str, Union[List[str], str]] = {
        "mset1": expand("sambamba/markdup/{sample}.bam", sample=tested_samples),
        "mset2": expand("sambamba/markdup/{sample}.bam", sample=reference_samples),
        "cset": f"medips/coupling/{wildcards.model_name}.RDS",
    }

    # Looking for input samples, if any
    # ... in reference level,
    reference_input_samples: List[str] = get_input_per_level(
        wildcards=snakemake.io.Wildcards(fromdict={"level": reference}), design=design
    )
    if len(reference_input_samples) > 0:
        medips_meth_coverage_input["iset1"] = expand(
            "sambamba/markdup/{sample}.bam", sample=reference_input_samples
        )

    # ... and in tested level.
    tested_input_samples: List[str] = get_input_per_level(
        wildcards=snakemake.io.Wildcards(fromdict={"level": tested}), design=design
    )
    if len(tested_input_samples) > 0:
        medips_meth_coverage_input["iset1"] = expand(
            "sambamba/markdup/{sample}.bam", sample=tested_input_samples
        )

    return medips_meth_coverage_input


def get_csaw_count_input(
    wildcards: snakemake.io.Wildcards,
    design: pandas.DataFrame = design,
    protocol: str = protocol,
    config: Dict[str, Any] = config,
) -> Dict[str, Union[str, List[str]]]:
    """
    Return expected list of input files for csaw-count
    """
    sample_list: List[str] = get_sample_list_from_model_name(
        model_name=wildcards.model_name,
        signal=wildcards.signal,
        design=design,
        config=config,
    )

    bam_prefix = get_bam_prefix(wildcards=wildcards, protocol=protocol)

    # Handle both single-ended and pair-ended sequencing libraries
    library: str = "se"
    if all(is_paired(sample=sample, design=design) for sample in sample_list):
        library = "pe"

    return {
        "bams": [f"{bam_prefix}/{sample}.bam" for sample in sample_list],
        "bais": [f"{bam_prefix}/{sample}.bam.bai" for sample in sample_list],
        "read_params": f"csaw/readparam.{library}.RDS",
        "design": "rsamtools/design.RDS",
        "fragment_length": "csaw/fragment_length.RDS",
    }


def get_csaw_count_params(
    wildcards: snakemake.io.Wildcards,
    design: pandas.DataFrame = design,
    protocol: str = protocol,
    config: Dict[str, Any] = config,
) -> str:
    """
    Return best parameters considering IO files list
    """
    extra: str = "filter = 10 "
    if str(wildcards.signal) in ["binned", "input_binned"]:
        extra += ", bin = TRUE, width=10000 "
    else:
        extra += ",  width=100 "

    # Atac-seq reads should be shifted to account for transposase size
    # if protocol_is_atac(protocol):
    # This step is automatically performed over the BAM
    # file through deeptools alignment sieve
    # extra += ", shift = 4"

    sample_list: List[str] = get_sample_list_from_model_name(
        model_name=wildcards.model_name,
        signal=wildcards.signal,
        design=design,
        config=config,
    )
    if any(is_paired(sample=sample, design=design) for sample in sample_list):
        if not all(is_paired(sample=sample, design=design) for sample in sample_list):
            raise ValueError(
                "Analysis of mixed Single-end / Pair-end libraries are not yet available"
            )
    else:
        fragment_lengths: str = ", ".join(
            map(
                str,
                [
                    has_fragment_size(
                        sample=sample, design=design, sample_is_paired=False
                    )
                    for sample in sample_list
                ],
            )
        )
        extra += f", ext=c({fragment_lengths})"

    return extra


def get_csaw_read_param(
    wildcards: snakemake.io.Wildcards,
    design: pandas.DataFrame = design,
    protocol: str = protocol,
) -> str:
    """
    Return best parameters for csaw::readParam()
    """
    extra: str = "dedup = TRUE, minq = 30"

    if str(wildcards.library) == "pe":
        # Pair-ended case
        extra += ", pe = 'both'"
    else:
        # Single-ended case
        extra += ", pe = 'none'"

    # Strand filtering for Medip-seq
    if protocol_is_medip(protocol):
        extra += ", forward=logical(0)"

    return extra


def get_csaw_filter_input(
    wildcards: snakemake.io.Wildcards,
    design: pandas.DataFrame = design,
    config: Dict[str, Any] = config,
) -> Dict[str, str]:
    """
    Return the list of input files expected by csaw_filter.R
    """
    csaw_filter_input: Dict[str, str] = {
        "counts": f"csaw/count/{wildcards.model_name}.tested.RDS",
    }
    input_list: List[str] = get_sample_list_from_model_name(
        model_name=wildcards.model_name, signal="input", design=design, config=config
    )
    if len(input_list) > 0:
        csaw_filter_input[
            "input_counts"
        ] = f"csaw/count/{wildcards.model_name}.input.RDS"
        csaw_filter_input[
            "input_binned"
        ] = f"csaw/count/{wildcards.model_name}.input_binned.RDS"
        csaw_filter_input["binned"] = f"csaw/count/{wildcards.model_name}.binned.RDS"
    else:
        csaw_filter_input["binned"] = f"csaw/count/{wildcards.model_name}.binned.RDS"

    return csaw_filter_input


####################
### Peak Calling ###
####################


def get_macs2_callpeak_input(
    wildcards: snakemake.io.Wildcards,
    design: pandas.DataFrame = design,
    protocol: str = "chip-seq",
) -> Dict[str, str]:
    """
    Return expected list of input files for Macs2 callpeak
    """
    macs2_callpeak_input = {}

    bam_prefix = get_bam_prefix(wildcards=wildcards, protocol=protocol)

    macs2_callpeak_input["treatment"] = f"{bam_prefix}/{wildcards.sample}.bam"

    # Use input files if available
    input_id = has_input(sample=wildcards.sample, design=design)
    if (input_id is not None) and (input_id != ""):
        macs2_callpeak_input["control"] = f"{bam_prefix}/{input_id}.bam"

    return macs2_callpeak_input


def get_macs2_params(
    wildcards: snakemake.io.Wildcards,
    effective_genome_size: int = effective_genome_size,
    design: pandas.DataFrame = design,
) -> str:
    """
    Return expected parameters for Macs2 callpeak
    """
    extra = f" --gsize {effective_genome_size} "
    down = is_paired(sample=wildcards.sample, design=design)
    if down:
        # Case pair-ended library
        extra += " --format BAMPE "
    else:
        # Case single-ended library, then a fragment size is required
        extra += " --format BAM "

        fragment_size = has_fragment_size(
            sample=wildcards.sample, design=design, sample_is_paired=False
        )
        if fragment_size:
            extra += f" --extsize {fragment_size} "
        else:
            raise ValueError(
                "Single-ended reads should have a "
                "`Fragment_size` associated in the design file."
                f"{wildcards.sample} has none."
            )

    return extra


def get_bedtools_intersect_macs2_input(
    wildcards: snakemake.io.Wildcards,
    protocol: str = protocol,
    peaktype_list: List[str] = peaktype_list,
    seacr_mode_list: List[str] = seacr_mode_list,
):
    """
    Return expected input files for bedtools intersect after macs2 calling.
    This is used to compute FRiP score.
    """
    peak_file: str = get_peak_file(
        wildcards=wildcards,
        peaktype_list=peaktype_list,
        seacr_mode_list=seacr_mode_list,
    )
    bam_prefix = get_bam_prefix(wildcards=wildcards, protocol=protocol)

    return {
        "right": peak_file.format(peaktype=wildcards.peaktype, sample=wildcards.sample),
        "left": f"{bam_prefix}/{wildcards.sample}.bam",
    }


def get_seacr_callpeak_input(
    wildcards: snakemake.io.Wildcards, design: pandas.DataFrame = design
) -> Dict[str, str]:
    """
    Return expected input files for SEACR shell interface
    """
    seacr_callpeak_input: Dict[str, str] = {
        "exp_bg": f"bedtools/genomecov/{wildcards.sample}.bg"
    }

    # If input are available, then use it
    input_signal = has_input(sample=wildcards.sample, design=design)
    if input_signal:
        seacr_callpeak_input["input_bg"] = f"bedtools/genomecov/{input_signal}.bg"

    return seacr_callpeak_input


def get_seacr_callpeak_params(
    wildcards: snakemake.io.Wildcards, design: pandas.DataFrame = design
) -> str:
    """
    Return best parameters for SEACR shell interface
    """
    seacr_callpeak_params = ""
    input_signal = has_input(sample=wildcards.sample, design=design)

    # If input are available, then use it
    if input_signal:
        seacr_callpeak_params += f" bedtools/genomecov/{input_signal}.bg "
    else:
        seacr_callpeak_params += " 0.01 "

    seacr_callpeak_params += (
        f" norm {wildcards.seacr_mode} seacr/raw/{wildcards.sample}"
    )

    return seacr_callpeak_params


def get_chipseeker_annotate_peak_single_sample_input(
    wildcards: snakemake.io.Wildcards,
    peaktype_list: List[str] = peaktype_list,
    seacr_mode_list: List[str] = seacr_mode_list,
) -> Dict[str, str]:
    """
    Return expected input file list for chipseeker, with correct peak-caller
    """
    bed_file: str = get_peak_file(
        wildcards=wildcards,
        peaktype_list=peaktype_list,
        seacr_mode_list=seacr_mode_list,
    )
    return {
        "bed": bed_file.format(peaktype=wildcards.peaktype, sample=wildcards.sample)
    }


############################
### Differential Binding ###
############################
def get_chipseeker_annotate_peak_from_ranges_input(
    wildcards: snakemake.io.Wildcards, protocol: str = protocol
) -> Dict[str, Any]:
    """
    Return expected list of input files for chipseeker annotate
    """
    if protocol_is_medip(protocol):
        return {"ranges": "medips/edger/{model_name}.RDS"}
    return {"ranges": "csaw/results/{model_name}.RDS"}


def get_chipseeker_genome_cov_single_sample_input(
    wildcards: snakemake.io.Wildcards,
    peaktype_list: List[str] = peaktype_list,
    seacr_mode_list: List[str] = seacr_mode_list,
) -> Dict[str, str]:
    """
    Return expected list of input files for chipseeker genome coverage
    """
    bed_file: str = get_peak_file(
        wildcards=wildcards,
        peaktype_list=peaktype_list,
        seacr_mode_list=seacr_mode_list,
    )
    return {
        "bed": bed_file.format(peaktype=wildcards.peaktype, sample=wildcards.sample)
    }


def get_edger_formula(
    wildcards: snakemake.io.Wildcards, config: Dict[str, Any] = config
) -> str:
    """
    Based on level wildcard, return edgeR stats formula
    """
    if "differential_peak_coverage" not in config.keys():
        raise KeyError("No differential coverage design was provided")

    for comparison in config["differential_peak_coverage"]:
        if comparison["model_name"] == wildcards.model_name:
            return comparison["formula"]

    raise ValueError("Could not find a statistical formula for: {wildcards.model_name}")


def get_deeptools_plot_correlation_params(wildcards: snakemake.io.Wildcards) -> str:
    """
    Based on plot-type (in wildcards) return correct parameters
    """
    labels: str = " ".join(design.index.tolist())
    deeptools_plot_correlation_params: str = (
        f" --plotFileFormat png --skipZeros --corMethod spearman --labels {labels} "
    )
    if str(wildcards.plot_type) == "heatmap":
        deeptools_plot_correlation_params += str(
            " --whatToPlot heatmap --colorMap Blues --plotNumbers --plotNumbers "
        )
    else:
        deeptools_plot_correlation_params += " --whatToPlot scatterplot "

    return deeptools_plot_correlation_params


def get_deeptools_plotheatmap_params(wildcards: snakemake.io.Wildcards) -> str:
    """
    Retrun correct parameters for deeptools plot heatmap
    """
    labels: str = " ".join(design.index.tolist())
    if "model_name" in wildcards.keys():
        labels = ""
    else:
        labels = f" --samplesLabel {labels} "
    return str(
        f"{labels} --plotFileFormat 'png' "
        "--colorMap 'Blues' --missingDataColor 1 "
        "--whatToShow 'heatmap and colorbar' "
    )


##############
### Motifs ###
##############


def get_homer_annotatepeaks_params(
    wildcards: snakemake.io.Wildcards,
    protocol: str = protocol,
    design: pandas.DataFrame = design,
) -> str:
    """
    Return best parameters for Homer annotate peaks
    """
    extra: str = " -homer2 "
    fs: str = has_fragment_size(
        sample=wildcards.sample, design=design, sample_is_paired=False
    )
    if fs:
        extra += f" -fragLength {fs} "

    if protocol_is_ogseq(protocol):
        extra += " -CpG "

    return extra


def get_peak_file_list(
    wildcards: snakemake.io.Wildcards,
    peaktype_list: List[str] = peaktype_list,
    seacr_mode_list: List[str] = seacr_mode_list,
):
    return [
        get_peak_file(
            wildcards=wildcards,
            peaktype_list=peaktype_list,
            seacr_mode_list=seacr_mode_list,
        ).format(peaktype=wildcards.peaktype, sample=wildcards.sample)
    ]


def get_homer_annotatepeaks_input(
    wildcards: snakemake.io.Wildcards,
    peaktype_list: List[str] = peaktype_list,
    seacr_mode_list: List[str] = seacr_mode_list,
) -> Dict[str, str]:
    """
    Return correct list of annotation input files for Homer
    """
    return {
        "genome": genome_fasta_path,
        "gtf": genome_annotation_path,
        "wig": f"data_output/Coverage/{wildcards.sample}.bw",
        "peak": f"homer/peaks/{wildcards.sample}.{wildcards.peaktype}.bed",
    }


def get_plot_footprints_input(
    wildcards: snakemake.io.Wildcards,
    protocol: str = protocol,
) -> Dict[str, str]:
    """
    Return required input files for ATACSeqQC
    """
    # bam_prefix: str = get_bam_prefix(wildcards=wildcards, protocol=protocol)
    bam_prefix: str = "sambamba/markdup"
    if "sample" in wildcards.keys():
        return {
            "bam": f"{bam_prefix}/{wildcards.sample}.bam",
            "bai": f"{bam_prefix}/{wildcards.sample}.bam.bai",
        }
    elif "factor_level" in wildcards.keys():
        return {
            "bam": f"{bam_prefix}/{wildcards.factor_level}.bam",
            "bai": f"{bam_prefix}/{wildcards.factor_level}.bam.bai",
        }
    raise ValueError("Could not find correct bam to return to FootPrints")


#############################
### Wildcards constraints ###
#############################

factor_level_list: List[str] = list(
    set(
        itertools.chain.from_iterable(
            [
                get_level_from_model_name(model, config).values()
                for model in get_model_names(config)
            ]
        )
    )
)


wildcard_constraints:
    sample=r"|".join(design.index),
    protocol=protocol,
    release=r"|".join([release, "99"]),
    build=r"|".join([build, "GRCm38"]),
    species=r"|".join([species, "mus_musculus"]),
    step=r"|".join([".raw", ""]),
    command=r"|".join(command_list),
    tool=r"|".join(["bowtie2", "sambamba", "deeptools"]),
    subcommand=r"|".join(["align", "markdup", "view", "alignment_sieve", "corrected"]),
    signal=r"|".join(["tested", "input", "binned", "input_binned"]),
    library=r"|".join(["se", "pe"]),
    chipseeker_plot=r"|".join(chipseeker_plot_list),
    model_name=r"|".join(get_model_names(config)),
    factor_level=r"|".join(factor_level_list),
    plot_type=r"|".join(deeptools_plot_type),
    peaktype=r"|".join(peaktype_list + ["gapped"] + seacr_mode_list),
    seacr_mode=r"|".join(seacr_mode_list),
    macs2_mode=r"|".join(peaktype_list + ["gapped"]),


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
    macs2_params: Dict[str, Any] = macs2_params,
    peak_type_and_mode_list: List[str] = peak_type_and_mode_list,
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
            "data_output/Coverage/{sample}.bw",
            sample=get_tested_sample_list(design=design),
        )

    if steps.get("calling", False):
        if macs2_params.get("broad", False) or macs2_params.get("narrow", False):
            expected_targets["macs2_broad"] = expand(
                "data_output/Peak_Calling/{peaktype}/Macs2/{sample}_peaks.xls",
                sample=get_tested_sample_list(design=design),
                peaktype=peaktype_list,
            )

        if len(peak_type_and_mode_list) > 0:
            expected_targets["dit_tss_broad"] = expand(
                "data_output/Peak_Calling/{peaktype}/{chipseeker_plot}/{sample}.png",
                sample=get_tested_sample_list(design=design),
                chipseeker_plot=chipseeker_plot_list,
                peaktype=peak_type_and_mode_list,
            )
            expected_targets["deeptools_heatmap"] = expand(
                "data_output/Heatmaps/{peaktype}/{command}.png",
                command=command_list,
                peaktype=peak_type_and_mode_list,
            )

        expected_targets["multiqc_coverage_report"] = "data_output/QC/Coverage.QC.html"

    if steps.get("diff_cov", False):
        comparison_list: List[str] = get_model_names(config)
        if len(comparison_list) > 0:
            expected_targets["annotated_csaw_tsv"] = expand(
                "data_output/Differential_Binding/{model_name}.tsv",
                model_name=comparison_list,
            )

            expected_targets["distance_to_tss"] = expand(
                "data_output/Differential_Binding/{model_name}/{chipseeker_plot}.png",
                model_name=comparison_list,
                chipseeker_plot=chipseeker_plot_list,
            )

            expected_targets["heatmap_per_model"] = expand(
                "data_output/Heatmaps/{peaktype}/{command}.{model_name}.png",
                model_name=comparison_list,
                peaktype=peak_type_and_mode_list,
                command=command_list,
            )

    if steps.get("motives", False):
        expected_targets["homer_annotate"] = expand(
            "data_output/Motifs/{peaktype}/{sample}_homer_annot.txt",
            peaktype=peak_type_and_mode_list,
            sample=get_tested_sample_list(design=design),
        )

        expected_targets["homer_denovo"] = expand(
            "data_output/Motifs/{peaktype}/{sample}_homer_annot.txt",
            peaktype=peak_type_and_mode_list,
            sample=get_tested_sample_list(design=design),
        )

        if (motifs_list != None) and (len(motifs_list) > 0):
            expected_targets["motifs_fingerprints"] = expand(
                "data_output/Motifs/Fingerprints/{target}/{target}.{motif}.png",
                target=get_tested_sample_list(design=design) + factor_level_list,
                motif=motifs_list,
            )

    return expected_targets
