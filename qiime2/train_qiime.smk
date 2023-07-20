###########################################
# ---- Qiime2 DB trainning workflow ----- #
###########################################
import sys
sys.path.append(".")


###############################################
# ------- Libraries/Functions set up -------- #
###############################################

import glob
import os
import shutil
import pandas as pd
import numpy as np
import re
from subprocess import call, run
import pathlib
from colorama import Fore, Style
from snakemake.logging import logger

###############################################
# ------------ Config variables ------------- #
###############################################

# Output folder for pipeline
BASE_FOLDER = Path("/home/usuario/Proyectos/Results/db_install_test/qiime_db_v3v4")

# Environments
CONDA_ENVS = {
    "qiime2": "qiime2-2022.8",
}

###############################################
# -------------- Static config -------------- #
###############################################

# Others
SILVA_EXPORTS_URL = "https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports"
SILVA_URLS = {
    "taxonomy": {
        "tax_ranks": f"{SILVA_EXPORTS_URL}/taxonomy/tax_slv_ssu_138.1.txt.gz",
        "tax_tree": f"{SILVA_EXPORTS_URL}/taxonomy/tax_slv_ssu_138.1.tre.gz",
        "tax_map": f"{SILVA_EXPORTS_URL}/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz",
    },
    "data": f"{SILVA_EXPORTS_URL}/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz",
}

PRETRAINED_CLASSIFIERS = {
    "silva-138-99-nb": {
        "data": "https://data.qiime2.org/2022.8/common/silva-138-99-nb-classifier.qza",
        "md5": "b8609f23e9b17bd4a1321a8971303310",
    },
}

# Available primer setups with NO adapters, only the primers that bind to the database
PRIMERS_SETUPS = {
    "KellyV3V4": {
        "fwd": {
            "primer": "CCTACGGGNGGCWGCAG",
            "adapter": "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
            "name": "U_Illu_Bakt_341_F",
        },
        "rev": {
            "primer": "GACTACHVGGGTATCTAATCC",
            "adapter": "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
            "name": "U_Illu_Bakt_805_R",
        }
    }
}
SELECTED_PRIMER_SETUP = "KellyV3V4"
PRIMERS = PRIMERS_SETUPS[SELECTED_PRIMER_SETUP]

###############################################
# -------------- Folder set up -------------- #
###############################################


# ---- Qiime2 DB trainning ----
QIIME_DB_CREATION_OUTPUT = f"{BASE_FOLDER}/qiime2_db"

# Silva download
SILVA_DB_OUTPUT = f"{QIIME_DB_CREATION_OUTPUT}/silva_download"

# qiime imports
QIIME_IMPORT_OUTPUT = f"{QIIME_DB_CREATION_OUTPUT}/import"
QIIME_IMPORT_TAX_TREE = f"{QIIME_IMPORT_OUTPUT}/import_tax_tree"
QIIME_IMPORT_TAX_MAP = f"{QIIME_IMPORT_OUTPUT}/import_tax_map"
QIIME_IMPORT_TAX_RANK = f"{QIIME_IMPORT_OUTPUT}/import_tax_rank"
QIIME_IMPORT_DATA= f"{QIIME_IMPORT_OUTPUT}/import_data"

# rescript
RESCRIPT_OUTPUT = f"{QIIME_DB_CREATION_OUTPUT}/rescript"
RESCRIPT_CULL = f"{RESCRIPT_OUTPUT}/cull_seqs"
RESCRIPT_PARSE_TAX = f"{RESCRIPT_OUTPUT}/parse_tax"
RESCRIPT_REVTRANS = f"{RESCRIPT_OUTPUT}/revtrans_data"
RESCRIPT_DEREPLICATE = f"{RESCRIPT_OUTPUT}/dereplicate"
RESCRIPT_EXTRACT = f"{RESCRIPT_OUTPUT}/extract"
RESCRIPT_DEREP_EXTRACTED = f"{RESCRIPT_OUTPUT}/dereplicate_extracted"

# Trainning
QIIME_DB_TRAINNING_OUTPUT = f"{QIIME_DB_CREATION_OUTPUT}/db_trainning"


###############################################
# ---------------- Output ------------------- #
###############################################
rule all:
    input:
        # ---- Qiime DB Trainning ----
        QIIME_DB_TRAINNING_OUTPUT,


###############################################
# ----------- Set up Silva data ------------- #
###############################################

rule download_silva:
    """
    Download SILVA files
    Output files MUST be uncompressed
    """
    output:
        folder = directory(Path(SILVA_DB_OUTPUT)),
        tax_map = Path(SILVA_DB_OUTPUT, Path(SILVA_URLS['taxonomy']['tax_map']).stem), # Stem converts a name like file.txt.gz in file.txt (removes last suffix)
        tax_ranks = Path(SILVA_DB_OUTPUT, Path(SILVA_URLS['taxonomy']['tax_ranks']).stem),
        tax_tree = Path(SILVA_DB_OUTPUT, Path(SILVA_URLS['taxonomy']['tax_tree']).stem),
        data = Path(SILVA_DB_OUTPUT, Path(SILVA_URLS['data']).stem),
    params:
        tax_map = SILVA_URLS["taxonomy"]["tax_map"],
        tax_ranks = SILVA_URLS["taxonomy"]["tax_ranks"],
        tax_tree = SILVA_URLS["taxonomy"]["tax_tree"],
        data = SILVA_URLS["data"],
        tries = 3,
    shell:
        """
        mkdir -p {output.folder}
        wget --show-progress --tries {params.tries} --output-document {output.tax_map}.gz {params.tax_map}
        wget --show-progress --tries {params.tries} --output-document {output.tax_ranks}.gz {params.tax_ranks}
        wget --show-progress --tries {params.tries} --output-document {output.tax_tree}.gz {params.tax_tree}
        wget --show-progress --tries {params.tries} --output-document {output.data}.gz {params.data}

        gzip --decompress {output.folder}/*.gz
        """

rule qiime_import_tax_ranks:
    input:
        tax_ranks = rules.download_silva.output.tax_ranks,
    output:
        folder = directory(Path(QIIME_IMPORT_TAX_RANK)),
        qza = Path(QIIME_IMPORT_TAX_RANK, "taxranks.qza"),
    params:
        feature_type = "'FeatureData[SILVATaxonomy]'"
    conda: CONDA_ENVS["qiime2"]
    shell:
        """
        qiime tools import \
            --type {params.feature_type} \
            --input-path {input.tax_ranks} \
            --output-path {output.qza}
        """

rule qiime_import_tax_map:
    input:
        tax_map = rules.download_silva.output.tax_map,
    output:
        folder = directory(Path(QIIME_IMPORT_TAX_MAP)),
        qza = Path(QIIME_IMPORT_TAX_MAP, "taxmap.qza"),
    params:
        feature_type = "'FeatureData[SILVATaxidMap]'"
    conda: CONDA_ENVS["qiime2"]
    shell:
        """
        qiime tools import \
            --type {params.feature_type} \
            --input-path {input.tax_map} \
            --output-path {output.qza}
        """


rule qiime_import_tax_tree:
    input:
        tax_tree = rules.download_silva.output.tax_tree,
    output:
        folder = directory(Path(QIIME_IMPORT_TAX_TREE)),
        qza = Path(QIIME_IMPORT_TAX_TREE, "taxtree.qza"),
    params:
        feature_type = "'Phylogeny[Rooted]'"
    conda: CONDA_ENVS["qiime2"]
    shell:
        """
        qiime tools import \
            --type {params.feature_type} \
            --input-path {input.tax_tree} \
            --output-path {output.qza}
        """

rule qiime_import_data:
    input:
        data = rules.download_silva.output.data,
    output:
        folder = directory(QIIME_IMPORT_DATA),
        qza = Path(QIIME_IMPORT_DATA, "silva-seqs.qza"),
    params:
        feature_type = "'FeatureData[RNASequence]'"
    conda: CONDA_ENVS["qiime2"]
    shell:
        """
        qiime tools import --type {params.feature_type} --input-path {input.data} --output-path {output.qza}
        """


###############################################
# ----------- Process Silva data ------------ #
###############################################

rule qiime_rescript_cull_seqs:
    input:
        data = rules.qiime_import_data.output.qza,
    output:
        folder = directory(RESCRIPT_CULL),
        data = Path(RESCRIPT_CULL, "silva-seqs.qza"),
    conda: CONDA_ENVS["qiime2"]
    threads: 60 # Low usage of threads but 40Gb of ram
    shell:
        """
        qiime rescript cull-seqs \
            --i-sequences {input.data} \
            --o-clean-sequences {output.data} \
            --p-n-jobs {threads}
        """


rule qiime_rescript_taxonomy:
    input:
        data = rules.qiime_rescript_cull_seqs.output.data,
        tax_ranks = rules.qiime_import_tax_ranks.output.qza,
        tax_tree = rules.qiime_import_tax_tree.output.qza,
        tax_map = rules.qiime_import_tax_map.output.qza,
    output:
        folder = directory(RESCRIPT_PARSE_TAX),
        tax = Path(RESCRIPT_PARSE_TAX, "silva-tax.qza"),
    conda: CONDA_ENVS["qiime2"]
    shell:
        """
        qiime rescript parse-silva-taxonomy \
            --i-taxonomy-tree {input.tax_tree} \
            --i-taxonomy-map {input.tax_map} \
            --i-taxonomy-ranks {input.tax_ranks} \
            --p-include-species-labels \
            --o-taxonomy {output.tax}
        """

rule qiime_rescript_dereplicate:
    input:
        data = rules.qiime_rescript_cull_seqs.output.data,
        tax = rules.qiime_rescript_taxonomy.output.tax,
    output:
        folder = directory(RESCRIPT_DEREPLICATE),
        data = Path(RESCRIPT_DEREPLICATE, "silva-seqs.qza"),
        tax = Path(RESCRIPT_DEREPLICATE, "silva-tax.qza"),
    conda: CONDA_ENVS["qiime2"]
    threads: 60
    shell:
        """
        qiime rescript dereplicate \
            --i-sequences {input.data} \
            --i-taxa {input.tax} \
            --o-dereplicated-sequences {output.data} \
            --o-dereplicated-taxa  {output.tax} \
            --p-threads {threads}
        """

###############################################
# ------------- Train database -------------- #
###############################################
rule qiime_db_trainning:
    input:
        data = rules.qiime_rescript_cull_seqs.output.data,
        tax = rules.qiime_rescript_taxonomy.output.tax,
    output:
        folder = directory(QIIME_DB_TRAINNING_OUTPUT),
        classifier = Path(QIIME_DB_TRAINNING_OUTPUT, "silva-classifier-full.qza"),
        observed_tax = Path(QIIME_DB_TRAINNING_OUTPUT, "oberserved-taxonomy.qza"),
        classifier_evaluation = Path(QIIME_DB_TRAINNING_OUTPUT, "classifier-evaluation.qzv"),
    conda: CONDA_ENVS["qiime2"]
    threads: 60 # TODO: It seems it is not parallel so this is not needed (?)
    shell:
        """
        # Train classifier and evaluate
        qiime rescript evaluate-fit-classifier \
            --i-sequences {input.data} \
            --i-taxonomy {input.tax}  \
            --o-classifier {output.classifier} \
            --o-observed-taxonomy {output.observed_tax} \
            --o-evaluation {output.classifier_evaluation} \
            --p-n-jobs {threads} \
            --verbose

        # # Only train classifier (faster)
        # qiime feature-classifier fit-classifier-naive-bayes \
        #     --i-reference-reads {input.data} \
        #     --i-reference-taxonomy {input.tax} \
        #     --o-classifier {output.classifier} \
        #     --verbose
        """


