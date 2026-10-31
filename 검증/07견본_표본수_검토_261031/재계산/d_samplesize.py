"""표본수계산 틀 3·5·6·7장 숫자 재계산 + S_w 를 data/ 원표에서 다시 구함(코드 24 를 import 하지 않음)."""
import csv, math, json, sys
from statistics import NormalDist
sys.stdout.reconfigure(encoding='utf-8')
BASE = sys.argv[1]
RAW = BASE + '/data/'
CAND = ['Cutibacterium', 'Staphylococcus', 'Corynebacterium']
EPS = 1e-5
z = NormalDist().inv_cdf

# ---- S_w from raw
def load(tag, site, tpkey):
    si = list(csv.DictReader(open(RAW + tag + '_sample_info.tsv', encoding='utf-8-sig'), delimiter='\t'))
    sid = {r['run']: r for r in si if r.get('body_site') == site and r.get('group') != '음성대조'}
    with open(RAW + tag + '_genus_pct.tsv', encoding='utf-8', newline='') as fh:
        rd = csv.reader(fh, delimiter='\t'); hdr = next(rd); rows = [r for r in rd if r[0] in sid]
    gen = hdr[1:]
    X = {g: [float(r[j]) if r[j] not in ('', 'NA') else 0.0 for r in rows] for j, g in enumerate(gen, 1)}
    rest = [sum(X[g][i] for g in gen if g not in CAND) for i in range(len(rows))]
    lr = {c: [math.log((X[c][i] + EPS) / (rest[i] + EPS)) for i in range(len(rows))] for c in CAND}
    subj = [sid[r[0]]['subject'] for r in rows]; tp = [sid[r[0]][tpkey] for r in rows]
    return lr, subj, tp, si

def groups(lr, subj, tp, tps):
    d = {}
    for s, t, v in zip(subj, tp, lr):
        d.setdefault(s, {})[t] = v
    return [[d[k][t] for t in tps] for k in sorted(d) if all(t in d[k] for t in tps)]

def sw1(g):
    n, k = len(g), len(g[0])
    return math.sqrt(sum((x - sum(r) / k) ** 2 for r in g for x in r) / (n * (k - 1)))

def sw2(g):
    n, k = len(g), len(g[0]); gm = sum(x for r in g for x in r) / (n * k)
    rm = [sum(r) / k for r in g]; cm = [sum(r[j] for r in g) / n for j in range(k)]
    sse = sum((g[i][j] - rm[i] - cm[j] + gm) ** 2 for i in range(n) for j in range(k))
    return math.sqrt(sse / ((n - 1) * (k - 1)))

lr15, s15, t15, si15 = load('15_season_16S', 'Forehead', 'collection_date')
tps15 = sorted(set(t15))
print('15 timepoints', tps15, 'n samples', len(s15))
eth = {}
for r in si15:
    if r.get('body_site') == 'Forehead' and r.get('group') != '음성대조':
        eth[r.get('ethnicity')] = eth.get(r.get('ethnicity'), 0) + 1
print('15 ethnicity', eth, '| region', {r.get('region') for r in si15})
lr07, s07, t07, _ = load('07_skincare_16S', 'Cheek', 'timepoint')
SW = {}
for c in CAND:
    g = groups(lr15[c], s15, t15, tps15)
    g7 = groups(lr07[c], s07, t07, ['2W', '4W'])
    SW[c] = dict(s15_1=sw1(g), s15_2=sw2(g), s07_1=sw1(g7), s07_2=sw2(g7), n15=len(g), n07=len(g7))
    print(c, {k: round(v, 4) for k, v in SW[c].items()})

# ---- constants
za = z(1 - 0.05 / 2); zb = z(0.80)
c1 = 2 * (za + zb) ** 2
za3 = z(1 - 0.05 / 6)
c3 = 2 * (za3 + zb) ** 2
ctost = 2 * (z(0.95) + z(1 - 0.20 / 2)) ** 2
print('const m=1 %.4f  m=3 %.4f  TOST %.4f' % (c1, c3, ctost))
for d in (0.5, 0.8, 1.0):
    print(' d', d, 'n0', c1 / d ** 2, math.ceil(c1 / d ** 2))

deltas = {'x1.5': math.log(1.5), 'x2': math.log(2), 'x3': math.log(3)}
deltas_r = {'x1.5': 0.41, 'x2': 0.69, 'x3': 1.10}
print('\nTable 5 (m=3, r=2) using ln and rounded deltas; S_w from raw (rounded 3dp in parens)')
for src, key in (('15계절 2원', 's15_2'), ('07 2원', 's07_2')):
    for c in CAND:
        sw = SW[c][key]
        row = []
        for k in deltas:
            n_ln = c3 * 1.5 * sw ** 2 / deltas[k] ** 2
            n_r = c3 * 1.5 * sw ** 2 / deltas_r[k] ** 2
            row.append('%s ln:%d(%.2f) rnd:%d' % (k, math.ceil(n_ln), n_ln, math.ceil(n_r)))
        print(src, c, 'Sw %.4f' % sw, ' | '.join(row))

print('\nSection 6')
def plan(sw, m_const, s, T=3, q=0.10, are=0.864):
    n0 = m_const * 1.5 * sw ** 2 / math.log(2) ** 2
    n1 = math.ceil(n0 / are)
    n1_mult = math.ceil(n0 * are)
    rec = math.ceil(n1 / (s ** T * (1 - q)))
    return n0, n1, n1_mult, rec, 2 * rec, 2 * rec * T
for src, key in (('07', 's07_2'), ('15', 's15_2')):
    sw = SW['Cutibacterium'][key]
    for s in (0.95, 0.90, 0.775):
        print(src, 'm=3 s', s, plan(sw, c3, s))
    print(src, 'm=1 s 0.95', plan(sw, c1, 0.95))

# ratios
print('\nSw ratio 15/07 (2원):', {c: round(SW[c]['s15_2'] / SW[c]['s07_2'], 3) for c in CAND})
print('n ratio:', {c: round((SW[c]['s15_2'] / SW[c]['s07_2']) ** 2, 3) for c in CAND})
print('delta x1.5 vs x3 n ratio', (math.log(3) / math.log(1.5)) ** 2)
print('cost 75*63', 75 * 63, '75*278', 75 * 278)
print('93/120', 93 / 120)
