import pandas as pd
metadata_file = "/home/usuario/Proyectos/CRC-16S-study/data/non-ffpe_metadata.tsv"
df = pd.read_csv(metadata_file, sep='\t')

# Extract tnm
x = df["tnm"].str.extract("p(T.*)(N.*)(M.*)|p(T.*)(N.*)")
x.columns = ["T1","N1","M1","T2","N2"]
x["tnm"] = df["tnm"]
x["T1"] = x["T1"].fillna(x["T2"])
x["N1"] = x["N1"].fillna(x["N2"])
x["M1"] = x["M1"].fillna("")
x = x.rename(columns={"T1": "tnm_t", "N1": "tnm_n", "M1": "tnm_m"})[["tnm","tnm_t","tnm_n","tnm_m"]]

# Check
temp = x[(~x["tnm"].isna()) & (x["tnm"]!="nd")]
cond = "p" + temp["tnm_t"] + temp["tnm_n"] + temp["tnm_m"] == temp["tnm"]
assert len(temp[~cond])==0

# Remove small subclassifications
x = x.replace("[a-z]","", regex=True)

# Add to metadata
df[["tnm_t","tnm_n","tnm_m"]] = x[["tnm_t","tnm_n","tnm_m"]]

df.to_csv(f"{metadata_file}_tnmsep.tsv", sep='\t', index=None)

