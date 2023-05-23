import re
import subprocess
import sys

from snakemake.utils import read_job_properties
from typing import Union

"""
From: https://github.com/Snakemake-Profiles/slurm/blob/master/%7B%7Bcookiecutter.profile_name%7D%7D/slurm-submit.py
"""


def get_best_partition(time: int) -> str:
    """
    Get partition name from time reservation

    Args:
        time (str): Time in minutes

    Returns:
        str: Partition name
    """
    if time < 360:
        return "shortq"
    elif time < 1140:
        return "mediumq"
    elif time < 10080:
        return "longq"
    elif time < 86400:
        return "verylongq"

    raise ValueError(f"Too much time requested: {time}, max is 86400 minutes")


jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript) or {}

time = job_properties.get("runtime", 359)
mem = job_properties.get("mem_mb", 1024 * 6)
partition = job_properties.get("partition", get_best_partition(time=time))
threads = job_properties.get("threads", 1)
error = job_properties.get("cluster", {}).get("error", "logs/slurm/slurm-%x-%j-%N.err")
output = job_properties.get("cluster", {}).get(
    "output", "logs/slurm/slurm-%x-%j-%N.out"
)


command = (
    f"sbatch --parsable --mem={mem}M --time={time} "
    f"--cpus-per-task={threads} --partition={partition} "
    f"--error='{error}' --output='{output}' --nodes=1 "
    "--job-name '{name}.{jobid}.slurm.snakejob.sh' "
    f"{jobscript} "
)

try:
    res = subprocess.check_output(args=command)
except subprocess.CalledProcessError as e:
    raise e

res = res.decode()
try:
    jobid = re.search(r"(\d+)", res).group(1)
except Exception as e:
    raise e

print(jobid)
