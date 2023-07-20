#################################
# ---- Qiime2 workflow ----- #
#################################
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
from pathlib import Path
from colorama import Fore, Style
from snakemake.logging import logger

def replace_proxy_rule_output(wildcards, expandable_structure, wildcards_dict, zip_expansion="yes"):
    """
    Takes path structure (fakepath/{barcode}/{barcode}_assembly.fasta) and expands the ProxyRule
    Inputs:
        wildcards = snakemake wildcard object produced in rule (it is not used in the function)
        expandable_structure =  "fakepath/{wildcard1}/{wildcard2}_assembly.fasta"
        wildcards_dict = {'wildcard1':['A','B','C], 'wildcard2':['D','E','F']}
    """
    if zip_expansion == "yes":
        l = expand(expandable_structure, zip, **wildcards_dict)
    else:
        l = expand(expandable_structure, **wildcards_dict)

    if l:
        return l
    else:
        raise Exception(f"Error: Structure {expandable_structure} does not contain any known wildcards.")

def get_threads(recommended: int, samples: int, available: int, regulator: float = 1):
    return int(max(min(available, regulator*available/samples), recommended))

###############################################
# ------------ Config variables ------------- #
###############################################

# Input
CLASSIFIER = Path("/home/usuario/Proyectos/Results/db_install_test/qiime2_db/db_trainning/silva-classifier.qza")
READS_FOLDER = Path("/home/usuario/Proyectos/Results/Kelly/KellyCCR/data/reads")
METADATA = Path("/home/usuario/Proyectos/Results/Kelly/KellyCCR/data/metadata/paper_metadata_pablo.tsv")
CONTAMINANTS = Path("/home/usuario/Proyectos/Results/Kelly/KellyBiome/parvimonas89/bacterias_contaminantes_limpias_KC.tsv")
VARIABLES_OF_INTEREST = ["origin","sample-type","sex","age","localization","tnm"] # ["origin","sample-type","sex","age","localization","tnm"]

# Output folder for pipeline
BASE_FOLDER = Path("/home/usuario/Proyectos/Results/Kelly/KellyCCR/output")

# Environments
CONDA_ENVS = {
    "qiime2": "qiime2-2022.8",
    "phyloseq": "base"
}

# Available primer setups with NO adapters, only the primers that bind to the database
PRIMERS_SETUPS = {
    "V3V4_kelly": {
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
SELECTED_PRIMER_SETUP = "V3V4_kelly"
PRIMERS = PRIMERS_SETUPS[SELECTED_PRIMER_SETUP]


###############################################
# -------------- Check input ---------------- #
###############################################

# Wildcards
sequencing_runs = list(pd.read_csv(METADATA, sep='\t', comment='#')["sequencing-run"].unique())
reads_wcs = glob_wildcards(Path(READS_FOLDER, '{barcode}_S{lane}_L001_R{read_orientation}_001.fastq.gz'))


def ensure_no_barcode_repetition(barcodes: list) -> bool:
    """
    Raise exception if glob_wildcards did not find exactly two copies of each barcode
    """
    rep_ids = {
        "repeated_ids": {
            "ids": [],
            "comment": "Following barcodes are repeated. Remember sample-ids, both in metadata and in reads must be unique."
        },
        "less_than_two": {
            "ids": [],
            "comment": "Following barcodes have less than two sets of reads. Illumina PE reads must have both R1 and R2"
        },
    }
    for i in set(barcodes):
        c = barcodes.count(i)            
        if c > 2:
            rep_ids["repeated_ids"]["ids"] += [i]
        elif c < 2:
            rep_ids["less_than_two"]["ids"] += [i]

    raise_error = False
    for k,v in rep_ids.items():
        if len(v["ids"]) > 0:
            logger.error("\nThere is a problem in your reads!")
            logger.error(v["comment"])
            logger.error(v["ids"])
            logger.error("")
            raise_error = True
            
    if raise_error:
        raise Exception
    else:
        return True

_ = ensure_no_barcode_repetition(reads_wcs.barcode)


wildcards_dict = {
    "barcode": list(set(reads_wcs.barcode)),
    "sequencing_run": sequencing_runs,
}
logger.info(f"Unique barcode count: {len(wildcards_dict['barcode'])}")
logger.info(f"Sequencing runs count: {len(wildcards_dict['sequencing_run'])}")
logger.info(f"First 5 barcodes: {'; '.join(wildcards_dict['barcode'][:5])}...")


# Check that metada follows expected structure
def ensure_metadata_structure(metadata_file: Path):
    metadata_df = pd.read_csv(metadata_file, sep='\t', comment='#')
    if metadata_df.columns[0] != "sample-id": 
        logger.error(f"First column in metadata is not 'sample-id'")
        raise Exception
    # elif any([True for i in metadata_df.columns if i not in ["origin","sample-type","subject","sequencing-run","control-applies-to"]])
    #     logger.error(f'A column is missing in metadata, expected columns (not in order) are: ["sample-id","origin","sample-type","subject","sequencing-run","control-applies-to"]')
    #     raise Exception

    return True


# Check that sample-ids are not repeated
def ensure_unique_barcodes(metadata_file: Path):
    metadata_df = pd.read_csv(metadata_file, sep='\t', comment='#')
    dups = metadata_df[metadata_df[["sample-id","sequencing-run"]].duplicated()]
    if len(dups) != 0: 
        logger.error(f"There are duplicated sample-ids in same sequencing-run: {';'.join(dups['sample-id'].to_list())}")
        raise Exception

    dups = metadata_df[metadata_df[["sample-id"]].duplicated()]
    if len(dups) != 0: 
        logger.error(f"There are duplicated sample-ids in different sequencing-runs: {';'.join(dups['sample-id'].to_list())}")
        raise Exception
    
    return True

# Check all barcodes are in the metadata first (DADA2 fails at the end of the procedure if not)
def ensure_barcode_presence(metadata_file: Path, wc_dict: dict = wildcards_dict):
    metadata_df = pd.read_csv(metadata_file, sep='\t', comment='#')

    # Metadata sample-ids not present in reads
    tmp = metadata_df[~metadata_df["sample-id"].isin(wc_dict["barcode"])]
    if len(tmp) != 0:
        logger.error(f"There are metadata sample-ids not present in reads:\n{tmp['sample-id']}")
        raise Exception
    
    # Reads not present in metadata
    tmp = [i for i in wc_dict["barcode"] if i not in metadata_df["sample-id"].to_list()]
    if len(tmp) != 0:
        logger.error(f"Some reads are not present in the metadata:\n{tmp}")
        raise Exception
    
    logger.info(f"All barcodes are present in metadata's sample-ids and vice versa")
    return True

# Check all variables of interest are in metadata
def ensure_voi_presence(metadata_file: Path, variables_of_interest: list):
    metadata_df = pd.read_csv(metadata_file, sep='\t', comment='#')
    for i in variables_of_interest:
        if i not in metadata_df.columns:
            logger.error(f"Variable of interest {i} is not a metadata column")
            raise Exception
    
    logger.info(f"All variables of interest in metadata")
    return True

# Check that all samples belong to a sequencing run (not empty or not None)
def ensure_belonging_seq_run(metadata_file: Path, wc_dict: dict = wildcards_dict, seq_run_col: str ="sequencing-run"):
    metadata_df = pd.read_csv(metadata_file, sep='\t', comment='#')

    # Check barcodes have a value in sequencing-run column
    no_belonging_df = metadata_df[
        (metadata_df[seq_run_col].isna()) 
        & (metadata_df["sample-id"].isin(wc_dict["barcode"]))
    ]
    if len(no_belonging_df) != 0:
        logger.error(f"There are sample that do not belong to any sequencing run. Please fix these: {';'.join(no_belonging_df['sample-id'])}")
        raise Exception    

    logger.info(f"All samples belong to a defined sequencing run")
    return True


_ = ensure_unique_barcodes(METADATA)
_ = ensure_barcode_presence(METADATA)
_ = ensure_metadata_structure(METADATA)
_ = ensure_voi_presence(METADATA, VARIABLES_OF_INTEREST)
_ = ensure_belonging_seq_run(METADATA)

logger.info(f"Input passed checks!\n")

###############################################
# -------------- Folder set up -------------- #
###############################################

# ---- Qiime2 main ----
QIIME_MAIN_OUTPUT = f"{BASE_FOLDER}/qiime2"

# Input preparation
SEPARATED_READS = f"{QIIME_MAIN_OUTPUT}/separate_reads_per_run"
QIIME_INPUT_PREP = f"{QIIME_MAIN_OUTPUT}/input_preparation"
QIIME_METADATA_TO_QZV = f"{QIIME_INPUT_PREP}/metadata"
QIIME_READS_TO_QZA = f"{QIIME_INPUT_PREP}/demux"
QIIME_TRIMM = f"{QIIME_INPUT_PREP}/trimmed_reads"

# ASVs
QIIME_ASV_OUTPUT = f"{QIIME_MAIN_OUTPUT}/asvs"
QIIME_DADA2 = f"{QIIME_ASV_OUTPUT}/dada2"
QIIME_MERGE_DADA2 = f"{QIIME_ASV_OUTPUT}/dada2_merge"
QIIME_CLASSIFY_ASVS = f"{QIIME_ASV_OUTPUT}/classification"

# View ASVs
QIIME_VIEW_ASV_OUTPUT = f"{QIIME_MAIN_OUTPUT}/view_asvs"
QIIME_BARPLOTS = f"{QIIME_VIEW_ASV_OUTPUT}/barplots"
QIIME_COLLAPSE_TO_TAXA = f"{QIIME_VIEW_ASV_OUTPUT}/collapse_to_taxa"

# Rarefaction
QIIME_RAREFACTION_OUTPUT = f"{QIIME_MAIN_OUTPUT}/rarefaction"
QIIME_PHYLOGENY = f"{QIIME_RAREFACTION_OUTPUT}/phylogeny"
QIIME_ALPHA_RAREFACTION = f"{QIIME_RAREFACTION_OUTPUT}/alpha_rarefaction"
QIIME_ALPHA_DIVERSITY = f"{QIIME_RAREFACTION_OUTPUT}/alpha_diversity"

# Core metrics
QIIME_METRICS_OUTPUT = f"{QIIME_MAIN_OUTPUT}/core_metrics"
QIIME_CORE_METRICS = f"{QIIME_METRICS_OUTPUT}/core_metrics"
QIIME_ALPHA_GROUP_SIGNIFICANCE = f"{QIIME_METRICS_OUTPUT}/alpha_group_significance"
QIIME_BETA_GROUP_SIGNIFICANCE = f"{QIIME_METRICS_OUTPUT}/beta_group_significance"

# Export
QIIME_EXPORT_RESULTS = f"{QIIME_MAIN_OUTPUT}/export_results"


# ---- Phyloseq ----
PHYLOSEQ_MAIN_OUTPUT = f"{BASE_FOLDER}/phyloseq"
PHYLOSEQ_SPLIT_TAX = f"{PHYLOSEQ_MAIN_OUTPUT}/split_taxonomy"
PHYLOSEQ_PLOT = f"{PHYLOSEQ_MAIN_OUTPUT}/plots"


###############################################
# ---------------- Output ------------------- #
###############################################
rule all:
    input:
        # ---- Qiime Main ----
        # Input
        # expand(Path(QIIME_METADATA_TO_QZV, "{sequencing_run}"), zip, sequencing_run=wildcards_dict['sequencing_run']),
        expand(Path(QIIME_READS_TO_QZA, "{sequencing_run}", "demux-pe.qza"), zip, sequencing_run=wildcards_dict['sequencing_run']),

        # ASVs
        expand(Path(QIIME_DADA2, "{sequencing_run}", "pet-table.qza"), zip, sequencing_run=wildcards_dict['sequencing_run']),
        Path(QIIME_CLASSIFY_ASVS, "taxonomy.qza"),

        # View ASVs
        Path(QIIME_BARPLOTS, "barplot.qzv"),
        Path(QIIME_COLLAPSE_TO_TAXA, "collapsed-table.qza"),

        # Rarefaction
        Path(QIIME_PHYLOGENY, "rooted-tree.qza"),
        Path(QIIME_ALPHA_RAREFACTION, "alpha-rarefaction.qzv"),
        Path(QIIME_ALPHA_DIVERSITY, "observed_features.qza"),

        # Core
        Path(QIIME_CORE_METRICS),
        Path(QIIME_ALPHA_GROUP_SIGNIFICANCE),
        expand(Path(QIIME_BETA_GROUP_SIGNIFICANCE, "{voi}"), zip, voi=VARIABLES_OF_INTEREST),

        # Export
        Path(QIIME_EXPORT_RESULTS),
        Path(QIIME_EXPORT_RESULTS, "feature-table.txt"),

        # ---- Phyloseq ----
        Path(PHYLOSEQ_SPLIT_TAX),



#############################################################
# -- Separate reads into folders based on sequencing run -- #
#############################################################

def get_belongings(wc, reads_folder: Path, metadata_file: Path):
    metadata_df = pd.read_csv(metadata_file, sep='\t', comment='#')
    bcs = metadata_df[metadata_df["sequencing-run"]==wc.sequencing_run]["sample-id"].to_list()
    belonging_reads = []
    for i in bcs:
        belonging_reads += glob.glob(f'{reads_folder}/{i}_S*_L001_R*_001.fastq.gz')

    # logger.info(f"Sequencing run {wc} has {len(bcs)} samples")
    return belonging_reads

rule separate_reads:
    input:
        reads = lambda wc: get_belongings(wc, READS_FOLDER, METADATA),
    output:
        folder = temp(directory(Path(SEPARATED_READS, "{sequencing_run}"))),
    shell:
        """
        mkdir -p {output.folder}
        cp {input.reads} {output.folder}
        """

###############################################
# ---------- Qiime2 Preparation ------------- #
###############################################

rule qiime_reads_to_qza:
    input:
        reads_folder = rules.separate_reads.output.folder,
    output:
        folder = directory(Path(QIIME_READS_TO_QZA, "{sequencing_run}")),
        reads_qza = Path(QIIME_READS_TO_QZA, "{sequencing_run}", "demux-pe.qza"),
    params:
        sample_type = "'SampleData[PairedEndSequencesWithQuality]'",
        input_format = "CasavaOneEightSingleLanePerSampleDirFmt",
    threads: 4
    conda: CONDA_ENVS["qiime2"]
    shell:
        """
        echo "Transformimg reads to .qza..."
        qiime tools import --type {params.sample_type} \
            --input-path {input.reads_folder} \
            --input-format {params.input_format} \
            --output-path {output.reads_qza}

        echo "Generating .qzv..."
        qiime demux summarize --i-data {output.reads_qza} --o-visualization {output.reads_qza}.qzv
        """

###############################################
# ------------- Qiime2 ASVs ----------------- #
###############################################
"""
overlap = fwd_len + rev_len - amplicon_len - (fwd_len - trunc_len_f) - (rev_len - trunc_len_r)
overlap_v2 = 250 + 250 - 464 - (250-251) - (250-251)
"""

rule qiime_dada2:
    input:
        reads_qza = rules.qiime_reads_to_qza.output.reads_qza,
        metadata = METADATA,
    output:
        folder = directory(Path(QIIME_DADA2, "{sequencing_run}")),
        table = Path(QIIME_DADA2, "{sequencing_run}", "pet-table.qza"),
        representative_sequences = Path(QIIME_DADA2, "{sequencing_run}", "rep-seqs-dada2.qza"),
        denoising_stats = Path(QIIME_DADA2, "{sequencing_run}", "denoising_stats.qza"),
    params:
        dada2_input_type = "denoise-paired",
        trim_left_f = len(PRIMERS["fwd"]["primer"]),
        trim_left_r = len(PRIMERS["rev"]["primer"]),   
        trunc_len_f = 251,
        trunc_len_r = 251,
    threads: get_threads(recommended=20, samples=len(wildcards_dict["sequencing_run"]), available=workflow.cores)
    conda: CONDA_ENVS["qiime2"]
    shell:
        """
        qiime dada2 {params.dada2_input_type} --i-demultiplexed-seqs {input.reads_qza} \
            --p-trim-left-f {params.trim_left_f} --p-trim-left-r {params.trim_left_r} \
            --p-trunc-len-f {params.trunc_len_f} --p-trunc-len-r {params.trunc_len_r} \
            --o-representative-sequences {output.representative_sequences} \
            --o-denoising-stats {output.denoising_stats} \
            --o-table {output.table} \
            --p-n-threads {threads} \
            --verbose

        # View table
        qiime feature-table summarize --i-table {output.table} --o-visualization {output.table}.qzv --m-sample-metadata-file {input.metadata}

        # View representative sequences and its stats
        qiime feature-table tabulate-seqs --i-data {output.representative_sequences} --o-visualization {output.representative_sequences}.qzv

        qiime metadata tabulate --m-input-file {output.denoising_stats} --o-visualization {output.denoising_stats}.qzv
        """

rule qiime_dada2_merge:
    """
    If reads come from different sequencing runs they must be analyzed separatedly by dada2
    One of the reasons is that dada2 removes background and cross-contamination, which is inherent to a sequencing run,
    and will be different in each one.
    Other reason can be that each run might have different length reads (150bp, 250bp, 350bp)
    """
    input:
        tables = lambda wc: replace_proxy_rule_output(wc, rules.qiime_dada2.output.table, wildcards_dict),
        representative_sequences = lambda wc: replace_proxy_rule_output(wc, rules.qiime_dada2.output.representative_sequences, wildcards_dict),
    output:
        folder = directory(Path(QIIME_MERGE_DADA2)),
        table = Path(QIIME_MERGE_DADA2, "pet-table.qza"),
        representative_sequences = Path(QIIME_MERGE_DADA2, "rep-seqs-dada2.qza"),
    conda: CONDA_ENVS["qiime2"]
    shell:
        """
        qiime feature-table merge \
            --i-tables {input.tables} \
            --o-merged-table {output.table}

        qiime feature-table merge-seqs \
            --i-data {input.representative_sequences} \
            --o-merged-data {output.representative_sequences}
        """


rule qiime_classify_asvs:
    """
    Benchmark in 64 thread CPU and 128Gb RAM (no virtual memory):
    - Fails with:
        - threads=20, reads_per_batch=default -> Above +128Gb RAM
    - Succeess with:
        - threads=40, reads_per_batch=10000 -> Around 40Gb RAM
        - threads=60, reads_per_batch=40000 -> Around 40Gb RAM
    """
    input:
        classifier = CLASSIFIER, 
        representative_sequences = rules.qiime_dada2_merge.output.representative_sequences,
    output:
        folder = directory(QIIME_CLASSIFY_ASVS),
        taxonomy = Path(QIIME_CLASSIFY_ASVS, "taxonomy.qza"),
    params:
        reads_per_batch = 40000,
    threads: 60
    conda: CONDA_ENVS["qiime2"]
    shell:
        """
        echo "Classifying ASVs with {input.classifier}..."
        qiime feature-classifier classify-sklearn \
            --i-reads {input.representative_sequences} \
            --i-classifier {input.classifier} \
            --o-classification {output.taxonomy} \
            --p-n-jobs {threads} --p-reads-per-batch {params.reads_per_batch} \
            --verbose

        qiime metadata tabulate --m-input-file {output.taxonomy} --o-visualization {output.taxonomy}.qzv
        """




###############################################
# ---------- View classification ------------ #
###############################################

rule qiime_barplots:
    input:
        table = rules.qiime_dada2_merge.output.table,
        taxonomy = rules.qiime_classify_asvs.output.taxonomy,
        metadata = METADATA,
    output:
        folder = directory(QIIME_BARPLOTS),
        barplot = Path(QIIME_BARPLOTS, "barplot.qzv"),
    threads: 60
    conda: CONDA_ENVS["qiime2"]
    shell:
        """
        qiime taxa barplot --i-table {input.table} \
            --i-taxonomy {input.taxonomy} \
            --m-metadata-file {input.metadata} \
            --o-visualization {output.barplot}
        """

rule qiime_collapse_to_taxa:
    input:
        table = rules.qiime_dada2_merge.output.table,
        taxonomy = rules.qiime_classify_asvs.output.taxonomy,        
    output:
        folder = directory(QIIME_COLLAPSE_TO_TAXA),
        table = Path(QIIME_COLLAPSE_TO_TAXA, "collapsed-table.qza"),
    params:
        taxa_level = 7
    threads: 60
    conda: CONDA_ENVS["qiime2"]
    shell:
        """
        qiime taxa collapse \
            --i-table {input.table} \
            --i-taxonomy {input.taxonomy} \
            --p-level {params.taxa_level} \
            --o-collapsed-table {output.table}

        qiime feature-table summarize \
            --i-table {output.table} \
            --o-visualization {output.table}.summary.qzv

        """


###############################################
# -------------- Rarefaction ---------------- #
###############################################

rule qiime_phylogeny:
    input:
        representative_sequences = rules.qiime_dada2_merge.output.representative_sequences,
    output:
        folder = directory(QIIME_PHYLOGENY),
        aln = Path(QIIME_PHYLOGENY, "aligned-rep-seqs-dada2.qza"),
        masked_aln = Path(QIIME_PHYLOGENY, "masked-aligned-rep-seqs-dada2.qza"),
        unrooted_tree = Path(QIIME_PHYLOGENY, "unrooted-tree.qza"),
        rooted_tree = Path(QIIME_PHYLOGENY, "rooted-tree.qza"),
    params:
        max_depth = 49754,
    threads: 60
    conda: CONDA_ENVS["qiime2"]
    shell:
        """
        qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences {input.representative_sequences} \
            --o-alignment {output.aln} \
            --o-masked-alignment {output.masked_aln} \
            --o-tree {output.unrooted_tree} \
            --o-rooted-tree {output.rooted_tree} \
            --p-n-threads {threads}
        """

rule qiime_alpha_rarefaction:
    input:
        rooted_tree = rules.qiime_phylogeny.output.rooted_tree,
        table = rules.qiime_dada2_merge.output.table,
    output:
        folder = directory(QIIME_ALPHA_RAREFACTION),
        alpha_rarefaction = Path(QIIME_ALPHA_RAREFACTION, "alpha-rarefaction.qzv"),
    params:
        max_depth = 300000,
    threads: 60
    conda: CONDA_ENVS["qiime2"]
    shell:
        """
        qiime diversity alpha-rarefaction \
            --i-table {input.table} \
            --i-phylogeny {input.rooted_tree} \
            --o-visualization {output.alpha_rarefaction}

        # --p-max-depth {params.max_depth}
        """

rule qiime_alpha_diversity:
    input:
        table = rules.qiime_collapse_to_taxa.output.table,
    output:
        folder = directory(QIIME_ALPHA_DIVERSITY),
        features = Path(QIIME_ALPHA_DIVERSITY, "observed_features.qza"),
    params:
        taxa_level = rules.qiime_collapse_to_taxa.params.taxa_level,
        metric = "'observed_features'",
    threads: 60
    conda: CONDA_ENVS["qiime2"]
    shell:
        """
        qiime diversity alpha --i-table {input.table} --p-metric {params.metric} --o-alpha-diversity {output.features}
        
        # View
        qiime metadata tabulate --m-input-file {output.features} --o-visualization {output.features}.qzv
        """

###############################################
# -------------- Core metrics --------------- #
###############################################


rule qiime_core_metrics:
    input:
        metadata = METADATA,
        rooted_tree = rules.qiime_phylogeny.output.rooted_tree,
        table = rules.qiime_dada2_merge.output.table,
    output:
        folder = directory(Path(QIIME_CORE_METRICS)),
        unweighted_unifrac_emperor = Path(QIIME_CORE_METRICS, "unweighted_unifrac_emperor.qzv"),
        weighted_unifrac_emperor = Path(QIIME_CORE_METRICS, "weighted_unifrac_emperor.qzv"),
        jaccard_emperor = Path(QIIME_CORE_METRICS, "jaccard_emperor.qzv"),
        bray_curtis_emperor = Path(QIIME_CORE_METRICS, "bray_curtis_emperor.qzv"),
        faith_pd_vector = Path(QIIME_CORE_METRICS, "faith_pd_vector.qza"),
        evenness_vector = Path(QIIME_CORE_METRICS, "evenness_vector.qza"),
        shannon_vector = Path(QIIME_CORE_METRICS, "shannon_vector.qza"),
        unweighted_unifrac_distance_matrix = Path(QIIME_CORE_METRICS, "unweighted_unifrac_distance_matrix.qza"),
    params:
        sampling_depth = 10000,
    threads: 60
    conda: CONDA_ENVS["qiime2"]
    shell:
        """
        rm -r {output.folder}

        qiime diversity core-metrics-phylogenetic \
            --i-phylogeny {input.rooted_tree} \
            --i-table {input.table} \
            --output-dir {output.folder} \
            --m-metadata-file {input.metadata} \
            --p-sampling-depth 10000 \
            --p-n-jobs-or-threads {threads}
        """

rule qiime_alpha_group_significance:
    input:
        metadata = METADATA,
        faith_pd_vector = rules.qiime_core_metrics.output.faith_pd_vector,
        evenness_vector = rules.qiime_core_metrics.output.evenness_vector,
        shannon_vector = rules.qiime_core_metrics.output.shannon_vector,
    output:
        folder = directory(Path(QIIME_ALPHA_GROUP_SIGNIFICANCE)),
        faith_pd_significance = Path(QIIME_ALPHA_GROUP_SIGNIFICANCE, "faith-pd-group-significance.qzv"),
        evenness_significance = Path(QIIME_ALPHA_GROUP_SIGNIFICANCE, "evenness-group-significance.qzv"),
        shannon_significance = Path(QIIME_ALPHA_GROUP_SIGNIFICANCE, "shannon_group-significance.qzv"),
    threads: 60
    conda: CONDA_ENVS["qiime2"]
    shell:
        """
        qiime diversity alpha-group-significance \
            --i-alpha-diversity {input.faith_pd_vector} \
            --m-metadata-file {input.metadata} \
            --o-visualization {output.faith_pd_significance}

        qiime diversity alpha-group-significance \
            --i-alpha-diversity {input.evenness_vector} \
            --m-metadata-file {input.metadata}  \
            --o-visualization {output.evenness_significance}

        qiime diversity alpha-group-significance \
            --i-alpha-diversity {input.shannon_vector} \
            --m-metadata-file {input.metadata} \
            --o-visualization {output.shannon_significance}
        """

rule qiime_beta_group_significance:
    input:
        unweighted_unifrac_distance_matrix = rules.qiime_core_metrics.output.unweighted_unifrac_distance_matrix,
        metadata = METADATA,
    output:
        folder = directory(Path(QIIME_BETA_GROUP_SIGNIFICANCE, "{voi}")),
        significance_viz = Path(QIIME_BETA_GROUP_SIGNIFICANCE, "{voi}", "unweighted-unifrac-significance.qzv"),
    threads: 60
    conda: CONDA_ENVS["qiime2"]
    shell:
        """
        qiime diversity beta-group-significance \
            --i-distance-matrix {input.unweighted_unifrac_distance_matrix} \
            --m-metadata-file {input.metadata} \
            --m-metadata-column {wildcards.voi} \
            --o-visualization {output.significance_viz} \
            --p-pairwise

        """



###############################################
# ---------- Export qiime results ----------- #
###############################################
rule qiime_export_results:
    input:
        table = rules.qiime_dada2_merge.output.table,
        taxonomy = rules.qiime_classify_asvs.output.taxonomy,
        rooted_tree = rules.qiime_phylogeny.output.rooted_tree,
    output:
        folder = directory(Path(QIIME_EXPORT_RESULTS)),
        biom_table = Path(QIIME_EXPORT_RESULTS, "feature-table.biom"),
        txt_biom_table = Path(QIIME_EXPORT_RESULTS, "feature-table.txt"),
        json_biom_table = Path(QIIME_EXPORT_RESULTS, "feature-table.json"),
        taxonomy = Path(QIIME_EXPORT_RESULTS, "taxonomy.tsv"),
        rooted_tree = Path(QIIME_EXPORT_RESULTS, "tree.nwk"),
    threads: 60
    conda: CONDA_ENVS["qiime2"]
    shell:
        """
        # Export to .biom
        qiime tools export --input-path {input.table} --output-path {output.folder}

        # Export to .json.biom
        biom convert -i {output.biom_table} -o {output.json_biom_table} --table-type="OTU table" --to-json

        # Export taxonomy
        qiime tools export --input-path {input.taxonomy} --output-path {output.folder}

        # Export tree
        qiime tools export --input-path {input.rooted_tree} --output-path {output.folder}

        # Transform to readable format
        biom summarize-table -i {output.biom_table} -o {output.biom_table}.summary.txt
        biom convert -i {output.biom_table} -o {output.txt_biom_table} \
            --to-tsv --header-key taxonomy --table-type "OTU table"
        """


###############################################
# -------------- Phyloseq ------------------- #
###############################################


rule split_taxonomy:
    input:
        taxonomy =  rules.qiime_export_results.output.taxonomy,
    output:
        folder = directory(Path(PHYLOSEQ_SPLIT_TAX)),
        taxonomy = Path(PHYLOSEQ_SPLIT_TAX, "separated_taxonomy.tsv"),
    params:
        taxa_levels = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    run:
        try: 
            df = pd.read_csv(input.taxonomy, sep='\t')
            df[params.taxa_levels] = df['Taxon'].str.split('; ',expand=True)

            df = df[params.taxa_levels]
            df.to_csv(output.taxonomy, sep='\t', index=False)       
        except Exception as e:
            print(e)
            raise Exception(e) 
