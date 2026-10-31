import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu, spearmanr
R = r"C:\Users\109so\OneDrive\Desktop\학교\micro3\새주제_v11"
its = pd.read_csv(R + r"\data\03_sensitive_ITS_genus_pct.tsv", sep="\t", index_col=0)
inf = pd.read_csv(R + r"\data\03_sensitive_ITS_sample_info.tsv", sep="\t", index_col=0)
g = inf.loc[its.index, "group"]
m = its["g__Malassezia"]
S = m[g == "Sensitive"]; N = m[g == "Non-sensitive"]
u = mannwhitneyu(N, S).statistic
print("Malassezia AUC(N>S)", u / (len(N) * len(S)), "medians S/N", S.median(), N.median())
# bootstrap person-level
rng = np.random.default_rng(1)
b = []
for _ in range(4000):
    s = rng.choice(S.values, len(S)); n = rng.choice(N.values, len(N))
    b.append(mannwhitneyu(n, s).statistic / (len(n) * len(s)))
print("boot CI", np.percentile(b, [2.5, 97.5]))
# Lactobacillus
b595 = pd.read_csv(R + r"\results\595_genus_relabund.tsv", sep="\t", index_col=0)
b42 = pd.read_csv(R + r"\data\03_sensitive_16S_genus_pct.tsv", sep="\t", index_col=0)
print([c for c in b595.columns if "Lactobac" in c][:5], [c for c in b42.columns if "Lactobac" in c][:5])
for gen in ["Lactobacillus", "Streptococcus"]:
    c595 = [c for c in b595.columns if c.endswith(gen) or c == gen]
    c42 = [c for c in b42.columns if c.endswith(gen) or c == gen]
    x = b595[c595[0]]; y = b42[c42[0]]
    scale = 100 if x.max() <= 1.0 else 1
    print(gen, c595, c42, "595 detect", (x > 0).mean() * 100, "median595(%)", x.median() * scale, "42 median", y.median(), "fold", (x.median() * scale) / y.median() if y.median() else None)
print(b595.iloc[:2, :3])
