import configparser
import yaml

from util_functions import *

PROJECT_DIR = config["all"]["PROJECT_DIR"]
MAPPING_FP = PROJECT_DIR + "/" + config["all"]["MAPPING"]
SAMPLE_IDS = get_sample(MAPPING_FP)

include: "rules/targets/targets.rules"
include: "rules/demux/dnabc.rules"
include: "rules/import/qiime_import.rules"

workdir: PROJECT_DIR

rule all:
    input: TARGET_ALL

onsuccess:
        print("Workflow finished, no error")
        shell("mail -s 'Workflow finished successfully' " + config["all"]["ADMIN_EMAIL"] + " < {log}")

onerror:
        print("An error occurred")
        shell("mail -s 'An error occurred' " + config["all"]["ADMIN_EMAIL"] + " < {log}")
