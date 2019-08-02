import configparser
import yaml

from util_functions import get_sample

PROJECT_DIR = config["all"]["PROJECT_DIR"]
QIIME_OUTPUT_DIR = PROJECT_DIR + "/" + config["import"]["QIIME"]
MAPPING_FP = PROJECT_DIR + "/" + config["all"]["MAPPING"]
SAMPLE_IDS = get_sample(MAPPING_FP)

DADA2_trim_left_f = config["denoise"]["DADA2"]["trim_left_f"]
DADA2_trunc_len_f = config["denoise"]["DADA2"]["trunc_len_f"]
DADA2_trim_left_r = config["denoise"]["DADA2"]["trim_left_r"]
DADA2_trunc_len_r = config["denoise"]["DADA2"]["trunc_len_r"]
DADA2_denoise_dir = (QIIME_OUTPUT_DIR + "/" + config["denoise"]["DENOISE"] +
                        "_fwd_" + str(DADA2_trim_left_f) + "-" + str(DADA2_trunc_len_f) +
                        "_rev_" + str(DADA2_trim_left_r) + "-" + str(DADA2_trunc_len_r)
                    )
CORE_METRIC_DIR = (DADA2_denoise_dir + "/" + config["diversity"]["core_metrics"] + 
                    "_sampling_depth_" + str(config["diversity"]["sampling_depth"]))

include: "rules/targets/targets.rules"
include: "rules/demux/dnabc.rules"
include: "rules/import/qiime_import.rules"
include: "rules/import/qiime_demux_stat.rules"
include: "rules/denoise/dada2.rules"
include: "rules/taxonomy/taxonomy.rules"
include: "rules/tree/tree.rules"
include: "rules/diversity/diversity.rules"
include: "rules/unassign/unassign.rules"
include: "rules/report/report.rules"

workdir: PROJECT_DIR

rule all:
    input: TARGET_ALL

onsuccess:
        print("Workflow finished, no error")
        shell("mail -s 'Workflow finished successfully' " + config["all"]["ADMIN_EMAIL"] + " < {log}")

onerror:
        print("An error occurred")
        shell("mail -s 'An error occurred' " + config["all"]["ADMIN_EMAIL"] + " < {log}")
