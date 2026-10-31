# -*- coding: utf-8 -*-
"""r1b — 15계절: 두 시점 사이 간격이 길수록 변화량 분산이 커지나(사람별 표류 → 독립 잡음 가정 점검).
07: 0W→2W · 2W→4W · 0W→4W 분산(0W 는 6속 오염 — 참고만). 평균 이동은 뺀 분산(2원과 같은 뜻)."""
import csv, json, math, os, sys
import numpy as np
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
HERE = os.path.dirname(os.path.abspath(__file__))
RAW = r'C:\Users\109so\OneDrive\Desktop\학교\micro3\새주제_v11\data'
CAND = ['Cutibacterium', 'Staphylococcus', 'Corynebacterium']
EPS = 1e-5

def load(tag, site, tpkey):
    si = {r['run']: r for r in csv.DictReader(open(os.path.join(RAW, tag + '_sample_info.tsv'), encoding='utf-8-sig'), delimiter='\t')}
    rd = csv.reader(open(os.path.join(RAW, tag + '_genus_pct.tsv'), encoding='utf-8'), delimiter='\t')
    hdr = next(rd); G = hdr[1:]
    rows = [r for r in rd if r[0] in si and si[r[0]].get('body_site') == site and si[r[0]].get('group') != '음성대조']
    X = np.array([[float(x) if x not in ('', 'NA') else 0.0 for x in r[1:]] for r in rows])
    meta = [(si[r[0]]['subject'], si[r[0]][tpkey]) for r in rows]
    rest = X[:, [g not in CAND for g in G]].sum(1)
    L = {c: np.log((X[:, G.index(c)] + EPS) / (rest + EPS)) for c in CAND}
    return meta, L

out = {}
for tag, site, key, tps in (('15_season_16S', 'Forehead', 'collection_date', None), ('07_skincare_16S', 'Cheek', 'timepoint', ['0W', '2W', '4W'])):
    meta, L = load(tag, site, key)
    tps = tps or sorted({t for _, t in meta})
    res = {}
    for c in CAND:
        d = {}
        for (s, t), v in zip(meta, L[c]):
            d.setdefault(s, {})[t] = v
        subj = [s for s in d if all(t in d[s] for t in tps)]
        M = np.array([[d[s][t] for t in tps] for s in subj])
        rows = []
        for i in range(len(tps)):
            for j in range(i + 1, len(tps)):
                dd = M[:, j] - M[:, i]
                rows.append({'pair': f'{tps[i]}→{tps[j]}', 'var_centered': float(np.var(dd, ddof=1)),
                             'implied_Sw(=sqrt(var/2))': float(math.sqrt(np.var(dd, ddof=1) / 2)), 'n': len(subj)})
        res[c] = rows
    out[tag] = res
json.dump(out, open(os.path.join(HERE, 'r1b_간격별분산.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
for tag, res in out.items():
    print('■', tag)
    for c, rows in res.items():
        print('  ', c, ' · '.join(f"{r['pair']} S_w≈{r['implied_Sw(=sqrt(var/2))']:.3f}" for r in rows))
