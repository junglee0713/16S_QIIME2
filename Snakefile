import configparser
import yaml

from scripts import util_functions

PROJECT_DIR = config["all"]["project_dir"]
QIIME_OUTPUT_DIR = PROJECT_DIR + "/QIIME_output"
MAPPING_FP = PROJECT_DIR + "/" + config["all"]["mapping"]
SAMPLE_IDS = util_functions.get_sample(MAPPING_FP)

trim_left_f = config["denoise"]["trim_left_f"]
trunc_len_f = config["denoise"]["trunc_len_f"]
trim_left_r = config["denoise"]["trim_left_r"]
trunc_len_r = config["denoise"]["trunc_len_r"]

DENOISE_DIR = (QIIME_OUTPUT_DIR + "/denoise" +
                "_fwd_" + str(trim_left_f) + "-" + str(trunc_len_f) +
                "_rev_" + str(trim_left_r) + "-" + str(trunc_len_r)
                    )
CORE_METRIC_DIR = (DENOISE_DIR + "/core-metrics" +
		  "_sampling_depth_" + str(config["diversity"]["sampling_depth"]))

CORE_METRIC_UNRAREFIED_DIR = (DENOISE_DIR + "/core-metrics-unrarefied")

include: "rules/targets/targets.rules"
include: "rules/demux/dnabc.rules"
include: "rules/import/qiime_import.rules"
include: "rules/import/qiime_demux_stat.rules"
include: "rules/denoise/denoise.rules"
include: "rules/taxonomy/taxonomy.rules"
include: "rules/tree/tree.rules"
include: "rules/diversity/diversity.rules"
include: "rules/unassign/unassign.rules"
include: "rules/report/report.rules"
include: "rules/dada2_species/dada2.rules"
include: "rules/vsearch/vsearch.rules"

workdir: PROJECT_DIR

rule all:
    input: TARGET_ALL

onsuccess:
        print("Workflow finished, no error")
        shell("mail -s 'Workflow finished successfully' " + config["all"]["admin_email"] + " < {log}")

onerror:
        print("An error occurred")
        shell("mail -s 'An error occurred' " + config["all"]["admin_email"] + " < {log}")
