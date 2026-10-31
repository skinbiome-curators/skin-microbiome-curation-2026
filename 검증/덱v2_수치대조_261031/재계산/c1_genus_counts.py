import pandas as pd, sys
R = r"C:\Users\109so\OneDrive\Desktop\학교\micro3\새주제_v11"
its = pd.read_csv(R + r"\data\03_sensitive_ITS_genus_pct.tsv", sep="\t", index_col=0)
print("ITS cols", its.shape)
det = (its > 0).sum()
print("ITS detected (any>0)", (det > 0).sum(), "excluding Unassigned", ((det > 0) & (det.index != "Unassigned")).sum())
b = pd.read_csv(R + r"\results\595_genus_relabund.tsv", sep="\t", index_col=0)
print("595 shape", b.shape)
b16 = pd.read_csv(R + r"\data\03_sensitive_16S_genus_pct.tsv", sep="\t", index_col=0)
print("42 16S shape", b16.shape)
print("16S 42 samples", b16.shape[0], "ITS samples", its.shape[0])
