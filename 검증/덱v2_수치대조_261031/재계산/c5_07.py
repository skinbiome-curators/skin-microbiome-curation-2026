import pandas as pd
R = r"C:\Users\109so\OneDrive\Desktop\학교\micro3\새주제_v11\data"
g = pd.read_csv(R + r"\07_skincare_16S_genus_pct.tsv", sep="\t", index_col=0)
s = pd.read_csv(R + r"\07_skincare_16S_sample_info.tsv", sep="\t", index_col=0)
print(s.columns.tolist()[:30])
print(g.shape, s.shape)
six = ["Wolbachia", "Sphingomonas", "Akkermansia", "Acetobacter", "Bacteroides", "Lachnospiraceae NK4A136 group"]
cols = [c for c in g.columns if any(c.endswith(x) or c == x for x in six)]
print(cols)
tp = [c for c in s.columns if "time" in c.lower() or "week" in c.lower()]
print(tp)
x = g[cols].sum(axis=1)
d = pd.DataFrame({"sum": x, "tp": s.loc[g.index, "timepoint"]})
print(d.groupby("tp")["sum"].agg(["count", "median", "min", "max"]))
print("subjects", s["subject"].nunique(), s["sex"].value_counts().to_dict(), "age", s["age"].min(), s["age"].max(), s["body_site"].value_counts().to_dict())
print(s.groupby("timepoint").size().to_dict())
