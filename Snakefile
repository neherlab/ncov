from os import environ
from socket import getfqdn
from getpass import getuser
from snakemake.utils import validate
import glob
import copy

configfile: "defaults/parameters.yaml"
validate(config, schema="workflow/schemas/config.schema.yaml")

# default build if none specified in config
if "builds" not in config:
    config["builds"] = {
        "global": {
            "subsampling_scheme": "region_global",
        }
    }

# if want to run builds for clusters
if os.path.isdir("swiss_profile/clusters/") and "cluster" in config['builds'] and "cluster_sampling" in config["subsampling"]:
    cluster_names = [w.replace("swiss_profile/clusters/cluster_","").replace(".txt", "") for w in glob.glob("swiss_profile/clusters/cluster_*.txt")]
    for new_clus in cluster_names:
        new_sample_scheme = "cluster_sampling_{}".format(new_clus)
        # use cluster build as 'template' for each individual cluster build
        config["builds"][new_clus] = copy.deepcopy(config["builds"]["cluster"])
        config["builds"][new_clus]["subsampling_scheme"] = new_sample_scheme
        config["builds"][new_clus]["title"] = config["builds"]["cluster"]["title"]+" - cluster {}".format(new_clus)        

        #make a new subsample scheme for each cluster - excluding that cluster from the non-focal set
        config["subsampling"][new_sample_scheme] = copy.deepcopy(config["subsampling"]["cluster_sampling"])
        config["subsampling"][new_sample_scheme]["global"]["exclude"] = "--exclude swiss_profile/clusters/cluster_{}.txt".format(new_clus)
    config["builds"].pop("cluster")  # get rid of the 'template' build

BUILD_NAMES = list(config["builds"].keys())

# Define patterns we expect for wildcards.
wildcard_constraints:
    # Allow build names to contain alpha characters, underscores, and hyphens
    # but not special strings used for Nextstrain builds.
    build_name = r'(?:[_a-zA-Z-](?!(tip-frequencies|gisaid|zh)))+',
    date = r"[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]"

localrules: download

# Create a standard ncov build for auspice, by default.
rule all:
    input:
        auspice_json = expand("auspice/ncov_{build_name}.json", build_name=BUILD_NAMES),
        tip_frequency_json = expand("auspice/ncov_{build_name}_tip-frequencies.json", build_name=BUILD_NAMES)

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"

# Include small, shared functions that help build inputs and parameters.
include: "workflow/snakemake_rules/common.smk"

# Include rules to handle primary build logic from multiple sequence alignment
# to output of auspice JSONs for a default build.
include: "workflow/snakemake_rules/main_workflow.smk"

# Include a custom Snakefile that specifies `localrules` required by the user's
# workflow environment.
if "localrules" in config:
    include: config["localrules"]

# Include custom rules defined in the config.
if "custom_rules" in config:
    for rule_file in config["custom_rules"]:
        include: rule_file
