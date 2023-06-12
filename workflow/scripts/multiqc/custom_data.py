#! coding: utf-8


import logging
import yaml

from pathlib import Path
from typing import Any, Dict
from snakemake.shell import shell

# Python loggings
logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

# Global variables
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
ln_extra: str = snakemake.params.get(
    "ln_extra", "--symbolic --force --relative --verbose"
)

protocol: str = snakemake.params.get("protocol", "Unknown")
step: str = snakemake.params.get("step", "Unknown")


def link_image(source: str, destination: str, log: str, params: str = ln_extra) -> None:
    """
    Link an image to MultiQC expected destination

    Args:
        source (str): Path to image/data source
        destination (str): Path to image/data destination
        params (str): Parameters for `ln` (bash)
        log (str): Snakemake logging behavior
    """
    shell(f"ln {params} {source} {destination} {log}")


# Main program
url_dict: Dict[str, str] = {
    "chip-seq": "https://github.com/tdayris/bigr_epicure_pipeline#chip-seq",
    "atac-seq": "https://github.com/tdayris/bigr_epicure_pipeline#atac-seq",
    "cut&tag": "https://github.com/tdayris/bigr_epicure_pipeline#cuttag",
    "cut&run": "https://github.com/tdayris/bigr_epicure_pipeline#cutrun",
    "medip-seq": "https://github.com/tdayris/bigr_epicure_pipeline#medip-seq",
    "og-seq": "https://github.com/tdayris/bigr_epicure_pipeline#og-seq",
}
url: str = url_dict.get(protocol, "chip-seq")

# Define header
mqc_config: Dict[str, Any] = {
    "title": f"{protocol.capitalize()} {step.lower()} analysis report",
    "subtitle": "Quality control graphs aggregation",
    "intro_text": f"See material and methods on <a href='{url}'>official pipeline github page</a>.",
    "show_analysis_paths": False,
    "show_analysis_time": False,
    "custom_logo": snakemake.input["logo"],
    "custom_logo_url": "https://bioinfo_gustaveroussy.gitlab.io/bigr/webpage/",
    "custom_logo_title": "BiGR @ Gustave Roussy",
    "report_header_info": [
        {"Contact E-mail": "bigr@gustaveroussy.fr"},
        {"Application Type": protocol.lower()},
    ],
    "custom_data": {},
    "custom_content": [],
    "sp": {},
    "ignore_images": False,
}

# Add peak-calling genomic annotation bar-plot
mqc_config["custom_data"]["peak_annotation_barplot"]: Dict[str, str] = {
    "parent_id": "chipseeker",
    "parent_name": "CHiPSeeker",
    "parent_description": str(
        "ChIPseeker is an R package for annotating peak data analysis."
    ),
    "section_name": "Peak annotation",
    "description": str(
        "This barplot displays the number of peaks "
        "covering a given genomic annotation."
    ),
    "plot_name": snakemake.output["annobar"],
}
mqc_config["custom_content"].append("peak_annotation_barplot")
link_image(
    source=snakemake.input["annobar"],
    destination=snakemake.output["annobar"],
    log=log,
    params=ln_extra,
)

# Add Peak-calling distance to TSS bar-plot
mqc_config["custom_data"]["dist_tss_barplot"]: Dict[str, str] = {
    "parent_id": "chipseeker",
    "section_name": "Peak annotation",
    "description": str(
        "This barplot displays the number of peaks "
        "covering a given genomic annotation."
    ),
    "plot_name": "multiqc/config/peak_annotation_barplot_mqc.png",
}

with open(snakemake.output["mqc"], "w") as mqc_out:
    yaml.dump(data=mqc_config, stream=mqc_out, default_flow_style=False)
