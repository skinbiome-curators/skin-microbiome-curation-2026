import pandas as pd
R = r"C:\Users\109so\OneDrive\Desktop\학교\micro3\새주제_v11\data"
for f in ["03_sensitive_ITS_sample_info.tsv", "03_sensitive_16S_sample_info.tsv"]:
    s = pd.read_csv(R + "\\" + f, sep="\t")
    g = s.groupby("group")["age"].agg(["count", "mean", "median", "min", "max"])
    print(f); print(g)
    print("diff mean", g.loc["Non-sensitive","mean"] - g.loc["Sensitive","mean"], "diff median", g.loc["Non-sensitive","median"] - g.loc["Sensitive","median"])
    print(s["sex"].value_counts().to_dict(), s["body_site"].value_counts().to_dict())
