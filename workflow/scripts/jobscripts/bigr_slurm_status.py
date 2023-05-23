import subprocess
import logging
import re
import shlex
import sys
import time

"""
From: https://github.com/Snakemake-Profiles/slurm/blob/master/%7B%7Bcookiecutter.profile_name%7D%7D/slurm-status.py
"""

logger = logging.getLogger(__name__)
jobid = sys.argv[1]
status_attempt = 10

for i in range(status_attempt):
    res = {}
    try:
        sacct_res = subprocess.check_output(shlex.split(f"sacct -P -b -j {jobid} -n"))
        res = {
            x.split("|")[0]: x.split("|")[1]
            for x in sacct_res.decode().strip().split("\n")
        }
        break

    except subprocess.CalledProcessError as e:
        logger.error("sacct process error")
        logger.error(e)

    except IndexError as e:
        logger.error(e)
        pass

    try:
        sctrl_res = subprocess.check_output(
            shlex.split(f"scontrol -o show job {jobid}")
        )
        m = re.search(r"JobState=(\w+)", sctrl_res.decode())
        res = {jobid: m.group(1)}
        break
    except subprocess.CalledProcessError as e:
        logger.error("scontrol process error")
        logger.error(e)
        if i >= status_attempt - 1:
            print("failed")
            exit(0)
        else:
            time.sleep(1)


failed_status = (
    "BOOT_FAIL",
    "OUT_OF_MEMORY",
    "DEADLINE",
    "NODE_FAIL",
    "FAILED",
    "PREEMPTED",
    "TIMEOUT",
    "CANCELLED",
)

success_status = "COMPLETED"
status = res.get(jobid, "")

if status.startswith(failed_status):
    print("failed")
elif status.startswith(success_status):
    print("success")
else:
    print("running")
