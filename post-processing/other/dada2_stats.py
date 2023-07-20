import pandas as pd
import glob
from pathlib import Path
from zipfile import ZipFile
import plotly.express as px

# ---- DADA2 filtering stats ----
df_list = []
zip_files = glob.glob(
    "/home/usuario/Proyectos/Results/Kelly/KellyCCR/output/qiime2/asvs/dada2/*/denoising_stats.qza"
)
for zip_file in zip_files:
    # Variables
    zip_file = Path(zip_file)
    run_id = zip_file.parent.name

    # Read stats.tsv file
    with ZipFile(zip_file, "r") as handle:
        for f in handle.namelist()[1:]:
            if "data/stats.tsv" in f:
                with handle.open(f) as file:
                    df = pd.read_csv(file, comment="#", sep="\t")
                    df["run_id"] = run_id
                    df_list.append(df)

# Show stats
df = pd.concat(df_list)
input_reads = df["input"].median()
final_kept_pct = round((100 * df["non-chimeric"] / df["input"]).median(), 2)
final_median_reads = df["non-chimeric"].median()
print(f"Total input reads: {df['input'].sum()}")
print(f"Input median reads: {input_reads}")
print(f"Final median reads: {final_median_reads}")
print(f"Final pct of reads: {final_kept_pct}%")
print(f"Total final reads: {df['non-chimeric'].sum()}")


# ---- DADA2 asv length ----
asv_list = []
zip_files = glob.glob(
    "/home/usuario/Proyectos/Results/Kelly/KellyCCR/output/qiime2/asvs/dada2/*/rep-seqs-dada2.qza.qzv"
)
for zip_file in zip_files:
    # Variables
    zip_file = Path(zip_file)
    run_id = zip_file.parent.name

    # Read stats.tsv file
    with ZipFile(zip_file, "r") as handle:
        for f in handle.namelist()[1:]:
            if "data/sequences.fasta" in f:
                with handle.open(f) as file:
                    for line in file:
                        line = line.decode("utf-8").strip()
                        if line.startswith(">"):
                            feature_id = line.replace(">", "")
                        else:
                            feature_length = len(line)
                            feature_seq = line
                            asv_list.append(
                                {
                                    "feature_id": feature_id,
                                    "feature_length": feature_length,
                                    "feature_seq": feature_seq,
                                }
                            )

df = pd.DataFrame(asv_list).drop_duplicates()
df["feature_length"].median()
px.histogram(df, x="feature_length")


# ---- DADA2 merge stats ----
zip_file = "/home/usuario/Proyectos/Results/Kelly/KellyCCR/output/qiime2/asvs/dada2_merge/pet-table.qza.qzv"
with ZipFile(zip_file, "r") as handle:
    for f in handle.namelist()[1:]:
        if (
            f
            == "b01269b9-a29a-4801-aa08-31b32a2953e6/data/feature-frequency-detail.csv"
        ):
            with handle.open(f) as file:
                df = pd.read_csv(file, comment="#", sep=",")

df.columns = ["feature", "frequency"]

(df["frequency"] >= df["frequency"].median()).sum()


df[df["frequency"] > 5]["frequency"].median()
