#! coding: utf-8

import yaml

from typing import Dict

protocol: str = snakemake.params.get("protocol", "Unknown")
step: str = snakemake.params.get("step", "Unknown")
url_dict: Dict[str, str] = {
    "chip-seq": "https://github.com/tdayris/bigr_epicure_pipeline#chip-seq",
    "atac-seq": "https://github.com/tdayris/bigr_epicure_pipeline#atac-seq",
    "cut&tag": "https://github.com/tdayris/bigr_epicure_pipeline#cuttag",
    "cut&run": "https://github.com/tdayris/bigr_epicure_pipeline#cutrun",
    "medip-seq": "https://github.com/tdayris/bigr_epicure_pipeline#medip-seq",
    "og-seq": "https://github.com/tdayris/bigr_epicure_pipeline#og-seq",
}
url = url_dict.get(protocol, "chip-seq")

result = {
    "title": f"{step} analysis report",
    "subtitle": "Quality control graphs aggregation",
    "intro_text": f"See material and methods on <a href='{url}'>official pipeline github page</a>.",
    "show_analysis_paths": False,
    "show_analysis_time": False,
    "custom_logo": snakemake.input["logo"],
    "custom_logo_url": "https://bioinfo_gustaveroussy.gitlab.io/bigr/webpage/",
    "custom_logo_title": "BiGR @ Gustave Roussy",
    "report_header_info": [
        {"Contact E-mail": "bigr@gustaveroussy.fr"},
        {"Application Type": protocol},
    ],
}

with open(snakemake.output["yaml"], "w") as yaml_out_stream:
    yaml.dump(data=result, stream=yaml_out_stream, default_flow_style=False)
