#! coding: utf-8

import yaml

from snakemake.shell import shell

results = {
    "custom_data": {
        "annobar_list": {
            "id": "annobar_list",
            "section_name": "Genomic Annotation",
            "description": (
                "This barplot displays the number of peaks "
                "covering a given genomic annotation."
            ),
            "plot_name": "multiqc/config/genomic_annotation_barplot_mqc.png",
        }
    }
}

shell(
    "ln -sfrv {input.png} "
    "multiqc/config/genomic_annotation_barplot_mqc.png "
    "> {log} 2>&1"
)

with open(snakemake.output["cfg"], "w") as mqc_out:
    yaml.dump(data=results, stream=mqc_out, default_flow_style=False)
