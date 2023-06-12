#! coding: utf-8

import yaml

results = {
    "custom_data": {
        "frip_score": {
            "id": "frip_score",
            "section_name": "FRiP Score",
            "description": (
                "FRiP (Fragment of Read in Peaks) is a common "
                "quality control performed right after peak-calling. "
                "The higher the FRiP score is, the more specific "
                "was the signal over noise."
            ),
            "plot_type": "bargraph",
            "pconfig": {
                "id": "frip_score_barplot",
                "title": "Fragment of Read in Peaks",
                "ylab": "Percent of reads in Peaks",
            },
            "data": {}
        }
    }
}

with open(snakemake.input["frip_tsv"], "r") as frip_stream:
    for idx, line in enumerate(frip_stream):
        if idx == 0:
            continue
        chomp = line[:-1].split("\t")
        results["custom_data"]["data"][chomp[0]] = chomp[-1]

with open(snakemake.output["mqc"], 'w')  as mqc_out:
    yaml.dump(data=results, stream=mqc_out, default_flow_style=False)