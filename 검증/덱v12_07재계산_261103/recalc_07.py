# 독립 재계산 — 07 공개 화장품 시험 0W→4W, 6속 제외 로그비, 정확 Wilcoxon 부호순위, Holm(3속)
# code/24 를 import 하지 않고 명세(덱 8장 꼬리말 · 효능판정기준 초안 #4)만 보고 새로 짬. 표준 라이브러리만.
import csv, math, sys
sys.stdout.reconfigure(encoding='utf-8')
root = sys.argv[1]
info = {r['run']: r for r in csv.DictReader(open(f'{root}/data/07_skincare_16S_sample_info.tsv', encoding='utf-8'), delimiter='\t')}
rows = list(csv.DictReader(open(f'{root}/data/07_skincare_16S_genus_pct.tsv', encoding='utf-8'), delimiter='\t'))
genera = [g for g in rows[0].keys() if g != 'run']
CAND = ['Cutibacterium', 'Staphylococcus', 'Corynebacterium']
EXCL = ['Wolbachia', 'Sphingomonas', 'Akkermansia', 'Acetobacter', 'Bacteroides', 'Lachnospiraceae NK4A136 group']
miss = [g for g in CAND + EXCL if g not in genera]
print('없는 속:', miss)
EPS = 1e-5
lr = {}
for r in rows:
    tp, subj = info[r['run']]['timepoint'], info[r['run']]['subject']
    v = {g: float(r[g] or 0) for g in genera}
    rest = sum(v[g] for g in genera if g not in CAND and g not in EXCL)
    for c in CAND:
        lr[(subj, tp, c)] = math.log((v[c] + EPS) / (rest + EPS))
subs = sorted({k[0] for k in lr})
def exact_counts(n):
    m = n * (n + 1) // 2
    cnt = [0] * (m + 1); cnt[0] = 1
    for k in range(1, n + 1):
        for s in range(m, k - 1, -1):
            cnt[s] += cnt[s - k]
    return cnt
def wilcoxon(d, alt):
    d = [x for x in d if x != 0]
    n = len(d)
    order = sorted(range(n), key=lambda i: abs(d[i]))
    absd = [abs(d[i]) for i in order]
    assert len(set(absd)) == n, 'ties'
    wplus = sum(r + 1 for r, i in enumerate(order) if d[i] > 0)
    cnt = exact_counts(n); tot = 2 ** n
    ge = sum(cnt[wplus:]) / tot; le = sum(cnt[:wplus + 1]) / tot
    if alt == 'greater': return ge, min(1, 2 * min(ge, le)), n, wplus
    return min(1, 2 * min(ge, le)), min(1, 2 * min(ge, le)), n, wplus
def holm(ps):
    idx = sorted(range(len(ps)), key=lambda i: ps[i]); out = [0] * len(ps); run = 0
    for j, i in enumerate(idx):
        run = max(run, min(1, (len(ps) - j) * ps[i])); out[i] = run
    return out
alt = {'Cutibacterium': 'greater', 'Staphylococcus': 'greater', 'Corynebacterium': 'two-sided'}
res = {}
for c in CAND:
    d = [lr[(s, '4W', c)] - lr[(s, '0W', c)] for s in subs if (s, '4W', c) in lr and (s, '0W', c) in lr]
    res[c] = wilcoxon(d, alt[c])
hc = holm([res[c][0] for c in CAND]); h2 = holm([res[c][1] for c in CAND])
for c, a, b in zip(CAND, hc, h2):
    p, p2, n, w = res[c]
    print(f'{c:16s} n={n} W+={w} p({alt[c]})={p:.4f} p_two={p2:.4f} · Holm 주장방향 {a:.4f} · Holm 양측 {b:.4f}')
