"""07 견본 A-5 · C-1 · C-5 의 숫자를 data/ 원표에서 직접 다시 구한다(코드 24 를 import 하지 않음)."""
import csv, math, statistics as st, sys, json
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
BASE = sys.argv[1]  # 새주제_v11
RAW = BASE + '/data/'
CAND = ['Cutibacterium', 'Staphylococcus', 'Corynebacterium']
EPS = 1e-5

si = list(csv.DictReader(open(RAW + '07_skincare_16S_sample_info.tsv', encoding='utf-8-sig'), delimiter='\t'))
print('n sample_info rows', len(si))
print('columns:', list(si[0].keys())[:30])
for col in ('body_site', 'group', 'sex', 'age', 'age_group', 'region', 'timepoint', 'project', 'panel'):
    vals = {}
    for r in si:
        vals[r.get(col)] = vals.get(r.get(col), 0) + 1
    print(col, vals)
subj = sorted({r['subject'] for r in si})
print('n subjects', len(subj))

# observed per timepoint
for col in ('Observed', 'observed_otus_paper', 'read_count', 'reads_paper', 'n_genus_detected'):
    for t in ('0W', '2W', '4W'):
        v = [float(r[col]) for r in si if r['timepoint'] == t and r[col] not in ('', 'NA')]
        print(col, t, 'n', len(v), 'median', st.median(v) if v else None)

# genus table
with open(RAW + '07_skincare_16S_genus_pct.tsv', encoding='utf-8', newline='') as fh:
    rd = csv.reader(fh, delimiter='\t')
    hdr = next(rd)
    rows = list(rd)
sidx = {r['run']: r for r in si}
rows = [r for r in rows if r[0] in sidx and sidx[r[0]]['body_site'] == 'Cheek']
print('genus rows', len(rows), 'n genera', len(hdr) - 1)
genera = hdr[1:]
X = {g: [float(r[j]) if r[j] not in ('', 'NA') else 0.0 for r in rows] for j, g in enumerate(genera, 1)}
tps = [sidx[r[0]]['timepoint'] for r in rows]
subs = [sidx[r[0]]['subject'] for r in rows]
tot = [sum(X[g][i] for g in genera) for i in range(len(rows))]
print('row sums range', min(tot), max(tot))

# confined at 1%
def share(g, t, thr):
    v = [x for x, tt in zip(X[g], tps) if tt == t]
    return sum((x >= thr) if thr > 0 else (x > 0) for x in v) / len(v)
conf = []
for g in genera:
    sh = {t: share(g, t, 1.0) for t in ('0W', '2W', '4W')}
    for t in sh:
        if sh[t] >= 0.8 and all(sh[o] <= 0.2 for o in sh if o != t):
            conf.append((g, t, sh))
print('confined 1%:', [(g, t, {k: round(v, 3) for k, v in sh.items()}) for g, t, sh in conf])
six = [g for g, t, sh in conf]
for t in ('0W', '2W', '4W'):
    s = [sum(X[g][i] for g in six) for i in range(len(rows)) if tps[i] == t]
    print('6-genus sum', t, 'median %.3f' % st.median(s), 'min %.2f max %.2f' % (min(s), max(s)), '>50%:', sum(x > 50 for x in s), '/', len(s))
# each genus median at 0W
for g in six:
    v = [X[g][i] for i in range(len(rows)) if tps[i] == '0W']
    print('   ', g, '0W median %.2f' % st.median(v), 'max %.2f' % max(v))

# logratio
def logratio(excl=()):
    rest = [sum(X[g][i] for g in genera if g not in CAND and g not in excl) for i in range(len(rows))]
    return {c: [math.log((X[c][i] + EPS) / (rest[i] + EPS)) for i in range(len(rows))] for c in CAND}

def paired(vals, t1, t2):
    d = {}
    for s, t, v in zip(subs, tps, vals):
        d.setdefault(s, {})[t] = v
    ks = sorted(k for k in d if t1 in d[k] and t2 in d[k])
    return [d[k][t2] - d[k][t1] for k in ks]

def hl_ci(d):
    n = len(d)
    walsh = sorted((d[i] + d[j]) / 2 for i in range(n) for j in range(i, n))
    # exact signed-rank null distribution
    M = n * (n + 1) // 2
    c = [0] * (M + 1); c[0] = 1
    for k in range(1, n + 1):
        for w in range(M, k - 1, -1):
            c[w] += c[w - k]
    acc = 0; cdf = []
    for x in c:
        acc += x; cdf.append(acc / 2 ** n)
    k = max(w for w in range(M + 1) if cdf[w] <= 0.025)
    return st.median(walsh), walsh[k], walsh[len(walsh) - 1 - k]

def holm(ps):
    m = len(ps); order = sorted(range(m), key=lambda i: ps[i]); adj = [0] * m; run = 0
    for r, i in enumerate(order):
        run = max(run, min(1, (m - r) * ps[i])); adj[i] = run
    return adj

alt = {'Cutibacterium': 'greater', 'Staphylococcus': 'greater', 'Corynebacterium': 'two-sided'}
for name, excl in (('raw', ()), ('excl_6', tuple(six))):
    lr = logratio(excl)
    res = {}
    for c in CAND:
        d = paired(lr[c], '0W', '4W')
        hl, lo, hi = hl_ci(d)
        p = stats.wilcoxon(d, alternative=alt[c], method='exact').pvalue
        p2 = stats.wilcoxon(d, alternative='two-sided', method='exact').pvalue
        res[c] = (len(d), st.median(d), hl, lo, hi, p, p2)
    ph = holm([res[c][5] for c in CAND])
    p2h = holm([res[c][6] for c in CAND])
    for c, a, b in zip(CAND, ph, p2h):
        n, med, hl, lo, hi, p, p2 = res[c]
        print(f'{name} 0W→4W {c}: n {n} median {med:+.4f} HL {hl:+.4f} CI [{lo:+.4f}, {hi:+.4f}] p({alt[c]}) {p:.4f} Holm {a:.4f} | p2 {p2:.4f} Holm2 {b:.4f}')

# percent-scale Cutibacterium medians per timepoint (for curiosity: original report style)
for t in ('0W', '2W', '4W'):
    v = [X['Cutibacterium'][i] for i in range(len(rows)) if tps[i] == t]
    s = [X['Staphylococcus'][i] for i in range(len(rows)) if tps[i] == t]
    print('pct median', t, 'Cuti %.2f' % st.median(v), 'Staph %.2f' % st.median(s))
