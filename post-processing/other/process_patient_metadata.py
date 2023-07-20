import pandas as pd
import numpy as np
from unidecode import unidecode

f = "/home/usuario/Descargas/metadata_cuest.xlsx"

# ---- Input ----

# Read excel
xls = pd.ExcelFile(f)
xls.sheet_names

# Read sheets
md_ccr = pd.read_excel(xls, "metadata_cuestionarios_ccr")
md_healthy = pd.read_excel(xls, "metadata_cuestionarios_sanos")

# Add sample.type
md_ccr["sample.type"] = "crc"
md_healthy["sample.type"] = "non-crc"

# Concatenate
df = pd.concat([md_ccr, md_healthy])

# Rename
column_corr = {
    "ID": "subject",
    "sample.type": "sample.type",
    "Edad": "age",
    "Sexo": "sex",
    "Altura": "height",
    "Peso": "weight",
    "Heces": "bristol_scale",
    "Frec-dep": "defecation_freq",
    "Enf-oral": "has_oral_disease",
    # "Tipo-enf-oral": "oral_disease",
    # "Dieta-esp": "has_diet",
    # "Tipo-dieta-esp": "diet",
    # "Intolerancias": "has_alimentary_intolerance",
    # "Tipo-intolerancia": "alimentary_intolerance",
    "Otra-enf": "has_other_disease",
    "Tipo-otra-enf": "other_disease",
    "Deporte": "exercises",
    # "Medicacion-cronica": "has_cronic_medication",
    # "Tipo-medicacion": "medication",
    "Trat-antb-largo": "had_lengthy_antibiotic_treatment",
    "Alcohol": "consumes_alcohol",
    "Tabaco": "consumes_tobacco",
    "Cafeina": "consumes_caffeine",
    "Trast-sueño": "sleep_disorder",
    "Ambiente": "environment",
}
df = df.rename(columns=column_corr)

# Filter
df = df[column_corr.values()]

# See values
for i in df.columns:
    print()
    print(i)
    print(df[i].unique())


# ---- Format values ----
# Change in bristol_scale
df.loc[df["bristol_scale"] == "ind", "bristol_scale"] = np.nan

# Change in environment
df.loc[df["environment"] == 2, "environment"] = "both"


# ---- Expand ";" separated columns ----
def test_fun(row, col):
    if row[col] is not np.nan:
        elements = row[col].split(";")
        for i in elements:
            i = unidecode(i)
            row[f"{col}.{i}"] = True
    return row


col = "other_disease"
expanded_col_df = df[["subject", "sample.type", col]].apply(
    lambda x: test_fun(x, col), axis=1
)
expanded_col_df = expanded_col_df.fillna(False)

# Merge
df = pd.merge(df, expanded_col_df[["subject", "other_disease.diabetes"]], on="subject")


# ---- See counts of each variable per sample.type ----
for i in df.columns:
    temp = df.value_counts(["sample.type", i], dropna=False).sort_index()
    temp = pd.concat(
        [
            (100 * temp[temp.index.map(lambda x: x[0] == "crc")] / 93).round(),
            (100 * temp[temp.index.map(lambda x: x[0] == "non-crc")] / 30).round(),
        ]
    )
    print()
    print(f"{i}".center(30, "-"))
    print(temp)

df.to_csv("/home/usuario/Proyectos/Results/patients.tsv", sep="\t")
