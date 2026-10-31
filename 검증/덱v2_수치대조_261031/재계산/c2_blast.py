import json, collections
R = r"C:\Users\109so\OneDrive\Desktop\학교\micro3\새주제_v11\results"
d = json.load(open(R + r"\포도상구균_BLAST_원본_core_nt.json", encoding="utf-8"))
for rep in d["reports"]:
    s = rep["results"]["search"]
    if s["query_title"] != "ASV0006": continue
    hits = s["hits"]
    print("query_len", s["query_len"], "n_hits", len(hits))
    sp = collections.Counter()
    allfull = True
    for h in hits:
        hs = h["hsps"][0]
        full = hs["identity"] == 375 and hs["align_len"] == 375 and hs["gaps"] == 0 and hs["query_from"] == 1 and hs["query_to"] == 375
        allfull &= full
        sp[h["description"][0]["sciname"]] += 1
    print("all 50 full 375/375 gap0:", allfull)
    print("species:", len(sp)); print(sp.most_common())
    hs = hits[0]["hsps"][0]
    print({k: hs[k] for k in ["identity","align_len","gaps","query_from","query_to","bit_score"]})
for rep in d["reports"]:
    s = rep["results"]["search"]
    if s["query_title"] == "ASV0006":
        L = sorted(h["len"] for h in s["hits"]); print("hit len min/max", L[:5], L[-3:])
