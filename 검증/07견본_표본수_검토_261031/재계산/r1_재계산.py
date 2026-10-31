# -*- coding: utf-8 -*-
"""방법 검토 r1 — 표본수 틀 · 07 견본의 예시 숫자 재계산 + 원자료에서 S_w · ICC · 분포 모양 확인.
저장소 코드(code/24)를 import 하지 않고 따로 짠다. 저장소 data/ 는 읽기만 한다.
출력: r1_재계산.json · 화면
"""
import csv, json, math, os, sys, time, tracemalloc
import numpy as np
from scipy import stats
from statsmodels.stats.power import TTestIndPower

sys.stdout.reconfigure(encoding='utf-8')
T0 = time.perf_counter(); tracemalloc.start()
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = r'C:\Users\109so\OneDrive\Desktop\학교\micro3\새주제_v11'
RAW = os.path.join(REPO, 'data')
CAND = ['Cutibacterium', 'Staphylococcus', 'Corynebacterium']
EPS = 1e-5
z = stats.norm.ppf
out = {}

# ── 1. 상수 ─────────────────────────────────────────────
zb = z(0.80)
C = {m: 2 * (z(1 - 0.05 / (2 * m)) + zb) ** 2 for m in (1, 3)}
C_one = {m: 2 * (z(1 - 0.05 / m) + zb) ** 2 for m in (1, 3)}   # 한쪽(주장 방향) α 0.05 — code/24 E1 관례
tost = 2 * (z(0.95) + z(0.90)) ** 2
out['constants'] = {'2(z+z)^2 m1': C[1], 'm3': C[3], 'one-sided m1': C_one[1], 'one-sided m3': C_one[3],
                    'TOST': tost, 'ln1.5': math.log(1.5), 'ln2': math.log(2), 'ln3': math.log(3),
                    'check d0.5/0.8/1.0': [math.ceil(C[1] / d ** 2) for d in (0.5, 0.8, 1.0)]}

# ── 2. 문서의 S_w (JSON) ───────────────────────────────
J = json.load(open(os.path.join(REPO, 'results', '효능판정_시범.json'), encoding='utf-8'))
K = 1.96 * math.sqrt(2)
SW = {}
for c in CAND:
    n15 = J['noise_15season']['genera'][c]
    i07 = J['cohort07']['CR_internal_2W4W']['raw'][c]
    SW[c] = {'15_1': n15['S_w_oneway'], '15_2': n15['S_w_twoway'],
             '07_1': i07['oneway'] / K, '07_2': i07['twoway'] / K}
out['S_w_from_json'] = SW

# ── 3. 5장 표 · 6장 표 재계산 ─────────────────────────
def n0(sw, delta, m=3, r=2, const=None):
    k = const if const is not None else C[m]
    return k * sw ** 2 * (1 + 1 / r) / delta ** 2

tab5 = {}
for src in ('15_2', '07_2'):
    for c in CAND:
        row = {}
        for lab, dl in (('x1.5', math.log(1.5)), ('x2', math.log(2)), ('x3', math.log(3))):
            row[lab] = {'exact_delta': math.ceil(n0(SW[c][src], dl)),
                        'rounded_delta': math.ceil(n0(SW[c][src], round(dl, 2)))}
        tab5[f'{src}|{c}'] = row
out['table5'] = tab5
doc5 = {'15_2|Cutibacterium': (189, 65, 26), '15_2|Staphylococcus': (184, 63, 25),
        '15_2|Corynebacterium': (123, 42, 17), '07_2|Cutibacterium': (27, 10, 4),
        '07_2|Staphylococcus': (83, 29, 12), '07_2|Corynebacterium': (44, 15, 6)}
cmp5 = []
for k, v in doc5.items():
    mine = tuple(tab5[k][l]['exact_delta'] for l in ('x1.5', 'x2', 'x3'))
    mine_r = tuple(tab5[k][l]['rounded_delta'] for l in ('x1.5', 'x2', 'x3'))
    cmp5.append({'row': k, 'doc': v, 'mine_exact_delta': mine, 'mine_rounded_delta': mine_r,
                 'match_exact': mine == v, 'match_rounded': mine_r == v})
out['table5_compare'] = cmp5

def chain(sw, m=3, r=2, s=0.95, T=3, q=0.10, are=0.864, delta=math.log(2)):
    a = n0(sw, delta, m, r)
    b = math.ceil(a / are)
    rec = math.ceil(b / (s ** T * (1 - q)))
    return {'n0': round(a, 3), 'n1': b, 'recruit_per_arm': rec, 'total': 2 * rec, 'samples': 2 * rec * T}

sec6 = {}
for src in ('07_2', '15_2'):
    sw = SW['Cutibacterium'][src]
    sec6[src] = {'m3_s95': chain(sw), 'm3_s90': chain(sw, s=0.90), 'm1_s95': chain(sw, m=1),
                 'm3_s775': chain(sw, s=0.775), 'm3_s95_rounded_delta': chain(sw, delta=0.69)}
# 셋 다 주장하면 가장 많이 필요한 균(5장 '가장 많이 필요한 균에 맞춘다')
sec6['07_2_all3_max'] = {c: chain(SW[c]['07_2']) for c in CAND}
sec6['15_2_all3_max'] = {c: chain(SW[c]['15_2']) for c in CAND}
out['section6'] = sec6

# 비교: t 분포(작은 n)로 푼 군당 인원 — 같은 α/m · 검정력
pw = TTestIndPower()
tsol = {}
for src in ('07_2', '15_2'):
    for c in CAND:
        sd = SW[c][src] * math.sqrt(1.5)
        d = math.log(2) / sd
        tsol[f'{src}|{c}'] = {'d': round(d, 3), 'normal_n0': round(C[3] / d ** 2, 2),
                              't_n': round(pw.solve_power(effect_size=d, alpha=0.05 / 3, power=0.8,
                                                           ratio=1, alternative='two-sided'), 2)}
out['t_based_n'] = tsol

# 비율 · 배수 문장
r2 = {c: SW[c]['15_2'] / SW[c]['07_2'] for c in CAND}
out['ratios'] = {'Sw_ratio_15over07_2way': r2, 'n_ratio': {c: v ** 2 for c, v in r2.items()},
                 'x1.5_vs_x3_n_ratio': (math.log(3) / math.log(1.5)) ** 2,
                 'samples_ratio_15_vs_07': [588 / 90, 690 / 102],
                 'waste_usd': [75 * 63, 75 * 278]}

# ── 4. 원자료에서 직접: 로그비 · S_w · ICC · 분포 ─────────
def load(tag, site, tpkey):
    si = {r['run']: r for r in csv.DictReader(open(os.path.join(RAW, tag + '_sample_info.tsv'), encoding='utf-8-sig'), delimiter='\t')}
    rd = csv.reader(open(os.path.join(RAW, tag + '_genus_pct.tsv'), encoding='utf-8'), delimiter='\t')
    hdr = next(rd)
    rows = [r for r in rd if r[0] in si and si[r[0]].get('body_site') == site and si[r[0]].get('group') != '음성대조']
    G = hdr[1:]
    X = np.array([[float(x) if x not in ('', 'NA') else 0.0 for x in r[1:]] for r in rows])
    meta = [(si[r[0]]['subject'], si[r[0]][tpkey], r[0]) for r in rows]
    return G, X, meta, si

def lr(G, X, excl=()):
    ci = [G.index(c) for c in CAND]
    rest_mask = np.array([g not in CAND and g not in excl for g in G])
    rest = X[:, rest_mask].sum(1)
    return {c: np.log((X[:, G.index(c)] + EPS) / (rest + EPS)) for c in CAND}

def wide(meta, v, tps):
    d = {}
    for (s, t, _), x in zip(meta, v):
        d.setdefault(s, {})[t] = x
    subj = sorted(k for k in d if all(t in d[k] for t in tps))
    return subj, np.array([[d[k][t] for t in tps] for k in subj])

def sw1(M):
    n, k = M.shape
    return math.sqrt(((M - M.mean(1, keepdims=True)) ** 2).sum() / (n * (k - 1)))

def sw2(M):
    n, k = M.shape
    R = M - M.mean(1, keepdims=True) - M.mean(0, keepdims=True) + M.mean()
    return math.sqrt((R ** 2).sum() / ((n - 1) * (k - 1)))

def icc_2way(M):
    """사람 분산 / (사람 + 잔차) — 시점 평균 이동은 뺀 일치도(consistency) ICC(3,1)"""
    n, k = M.shape
    msr = k * ((M.mean(1) - M.mean()) ** 2).sum() / (n - 1)
    R = M - M.mean(1, keepdims=True) - M.mean(0, keepdims=True) + M.mean()
    mse = (R ** 2).sum() / ((n - 1) * (k - 1))
    return (msr - mse) / (msr + (k - 1) * mse), max((msr - mse) / k, 0), mse

G7, X7, M7, SI7 = load('07_skincare_16S', 'Cheek', 'timepoint')
G15, X15, M15, SI15 = load('15_season_16S', 'Forehead', 'collection_date')
L7, L15 = lr(G7, X7), lr(G15, X15)
tp15 = sorted({t for _, t, _ in M15})
own = {}
dist = {}
for c in CAND:
    s7, W7 = wide(M7, L7[c], ['2W', '4W'])
    s15, W15 = wide(M15, L15[c], tp15)
    s7all, W7all = wide(M7, L7[c], ['0W', '2W', '4W'])
    i15 = icc_2way(W15); i7 = icc_2way(W7)
    d7 = W7[:, 1] - W7[:, 0]
    R15 = (W15 - W15.mean(1, keepdims=True) - W15.mean(0, keepdims=True) + W15.mean()).ravel()
    own[c] = {'07_n': len(s7), '07_sw1': sw1(W7), '07_sw2': sw2(W7), '15_n': len(s15), '15_k': len(tp15),
              '15_sw1': sw1(W15), '15_sw2': sw2(W15),
              'ICC_15_2way': i15[0], 'ICC_07_2W4W_2way': i7[0],
              'sigma_b2_15': i15[1], 'sw2_sq_15': i15[2],
              # S_w 추정의 95% CI(카이제곱 · 자유도) — 계획값의 불확실성
              'Sw07_2_CI95': [sw2(W7) * math.sqrt((len(s7) - 1) / stats.chi2.ppf(q, len(s7) - 1)) for q in (0.975, 0.025)],
              'Sw15_2_CI95': [sw2(W15) * math.sqrt(((len(s15) - 1) * (len(tp15) - 1)) / stats.chi2.ppf(q, (len(s15) - 1) * (len(tp15) - 1))) for q in (0.975, 0.025)],
              'Sw07_2_UCL80': sw2(W7) * math.sqrt((len(s7) - 1) / stats.chi2.ppf(0.20, len(s7) - 1)),
              }
    dist[c] = {'07_d_2W4W_kurtosis_excess': float(stats.kurtosis(d7)), '07_d_skew': float(stats.skew(d7)),
               '07_d_shapiro_p': float(stats.shapiro(d7).pvalue),
               '15_resid_kurtosis_excess': float(stats.kurtosis(R15)), '15_resid_shapiro_p': float(stats.shapiro(R15).pvalue),
               '07_zero_samples_(pct==0)': int((X7[:, G7.index(c)] == 0).sum()),
               '15_zero_samples_(pct==0)': int((X15[:, G15.index(c)] == 0).sum()),
               '07_min_logratio': float(L7[c].min()), '15_min_logratio': float(L15[c].min())}
out['own_from_raw'] = own
out['distribution'] = dist

# ANCOVA 대 변화량(기준 r=2) 분산 비 — 2원 모형 성분으로
anc = {}
for c in CAND:
    sb2, se2 = own[c]['sigma_b2_15'], own[c]['sw2_sq_15']
    for r in (1, 2):
        change = se2 * (1 + 1 / r)
        post = sb2 + se2
        ancova = post - sb2 ** 2 / (sb2 + se2 / r)
        anc[f'{c}|r{r}'] = {'ICC': sb2 / (sb2 + se2), 'var_change': change, 'var_post_only': post,
                            'var_ancova': ancova, 'ancova_over_change': ancova / change}
out['ancova_vs_change_15season'] = anc

# ── 5. 07 견본 A-5 숫자: 6속 · 50% 초과 수 · 검출률 ──────────
six = J['E0']['07']['confined_grid']['1%']['genera']
six_names = [g['genus'] for g in six]
out['six_genera_shares'] = {g['genus']: g['share'] for g in six}
idx = [G7.index(g) for g in six_names]
tps7 = [t for _, t, _ in M7]
sums0 = np.array([X7[i, idx].sum() for i in range(len(M7)) if tps7[i] == '0W'])
out['six_sum_0W'] = {'n': len(sums0), 'median': float(np.median(sums0)), 'min': float(sums0.min()),
                     'max': float(sums0.max()), 'n_over_50': int((sums0 > 50).sum())}
for t in ('2W', '4W'):
    s = np.array([X7[i, idx].sum() for i in range(len(M7)) if tps7[i] == t])
    out['six_sum_0W'][f'median_{t}'] = float(np.median(s))
# 관찰 OTU 중앙값(시료표의 observed_otus_paper · Observed)
for col in ('observed_otus_paper', 'Observed'):
    try:
        out[f'median_{col}'] = {t: float(np.median([float(SI7[r][col]) for _, tt, r in M7 if tt == t])) for t in ('0W', '2W', '4W')}
    except Exception as e:
        out[f'median_{col}'] = str(e)
# 06 참고 수치(A-5) — JSON 에서 그대로 읽어 대조
pc = J['cohort07']['comparisons']['0W→4W']['versions']
ref = {}
for v in ('raw', 'excl_1pct', 'excl_0.1pct', 'excl_0pct', 'excl_union'):
    ref[v] = {c: {'hl': pc[v][c]['hl'], 'ci': pc[v][c]['hl_ci95'], 'p_holm': pc[v][c]['p_holm'],
                  'p_two_holm': pc[v][c]['p_two_sided_holm']} for c in CAND}
out['A5_reference_json'] = ref

out['mem_time'] = {'sec': time.perf_counter() - T0, 'peak_MB': tracemalloc.get_traced_memory()[1] / 1e6}
json.dump(out, open(os.path.join(HERE, 'r1_재계산.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1, default=float)
print(json.dumps(out, ensure_ascii=False, indent=1, default=float))
