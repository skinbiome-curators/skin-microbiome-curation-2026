#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
피지 때문 아닌가 — Malassezia 착시 점검 ⑤ (2026-09-30 추가, 1회차 멘토링 제안)

⚠️ 멘토링 정리 범위를 넘은 추가 분석이다(Claude 판단으로 수행). 경위·결과는 온라인멘토링_정리.md §6.

1회차 멘토: "같은 사람에게서 Malassezia 와 피지 양을 기기로 잰 자료가 있을 것이다.
8장(착시 점검)에 보강 자료로 쓸 수 있다."

Malassezia 는 지방산 합성효소가 없어 피지에 의존한다. 그러니 민감군의 Malassezia 감소가
"민감군 피지가 적어서" 생긴 것이라면, 민감피부 신호가 아니라 피지 신호다. 세 가지를 본다.

  (1) 두 군의 피지가 다른가           03 민감 42명 · sebumeter_sebum_rigth (오른쪽 볼, 기기 측정)
  (2) 피지와 Malassezia 가 같이 움직이나   같은 42명 · Spearman
  (3) 피지를 넣어도 Malassezia 가 더하나   로지스틱 — 9번(나이 보정)과 같은 방식
        (가) 피지만 → (나) 피지 + Malassezia      증가분
        (다) 나이+피지 → (라) 나이+피지+Malassezia  증가분
  (4) 교차 확인                        01 노화 58명 · sebum · 전체 / 볼 30명만
  (5) scipy 교차 확인 (선택)            (1)(2)를 scipy 로 다시 계산 · scipy 없으면 건너뜀 · JSON 미저장

(5)는 대화 중 즉석 계산(scipy)을 옮겨 두었던 피지_즉석계산_20260930.py 의 역할을 흡수한 것이다
(2026-09-30 · 그 파일은 삭제).
'외부 라이브러리 0' 규칙은 (1)~(4)와 결과 JSON 에 적용되고, (5)는 독립 구현과 대조하는 출력일 뿐이다.

(3)의 부트스트랩은 사람 단위, seed 20260918 (다른 스크립트와 동일).
재대입 AUC 는 변수를 더하면 거의 항상 오르므로 한 명씩 빼는 교차검증(LOOCV)을 같이 낸다.

⚠️ 한계: 피지는 한 번, 볼 한 곳에서 잰 값이다. 01 노화는 다른 연구이므로 '재현'이 아니라
   '조건이 다른 자료에서의 일관성'이다.

결과: results/피지_점검.json
"""
import os, sys, csv, math, random, statistics as st
# 이 파일은 온라인멘토링_정리.md 옆(새주제_v11/온라인멘토링/)에 둔다 — 공유 통계 함수는 ../code/stats_util.py 에서 불러온다
# (2026-09-30 폴더 최상위 → 온라인멘토링/ 으로 옮기며 한 단계 위를 보도록 수정 · 2026-10-07 폴더 이름 코드 → code)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'code'))
from stats_util import save_json, RAW, TAB, f, spearman, auc_p, auc_ci

SEED, N_BOOT = 20260918, 4000
MAL = 'g__Malassezia'


def read_tsv(path):
    with open(path, encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def load(tag, sebum_col):
    info = read_tsv(os.path.join(RAW, f'{tag}_sample_info.tsv'))
    pct = {r['run']: r for r in read_tsv(os.path.join(RAW, f'{tag}_genus_pct.tsv'))}
    rows = []
    for r in info:
        seb, mal = f(r.get(sebum_col)), f(pct[r['run']].get(MAL))
        if seb is None or mal is None:
            continue
        grp = str(r.get('group', ''))
        y = 0 if grp.startswith('Non') else (1 if grp == 'Sensitive' else None)
        rows.append(dict(run=r['run'], y=y, age=f(r.get('age')), sebum=seb,
                         mal=mal, site=r.get('body_site')))
    return rows


# ── 9번과 같은 로지스틱·AUC·LOOCV (9번은 실행 스크립트라 import 하지 않고 같은 구현을 둔다) ──
def auc(y, x):
    pos = [x[i] for i in range(len(y)) if y[i] == 1]
    neg = [x[i] for i in range(len(y)) if y[i] == 0]
    if not pos or not neg:
        return float('nan')
    s = 0.0
    for p in pos:
        for n in neg:
            s += 1.0 if p > n else (0.5 if p == n else 0.0)
    return s / (len(pos) * len(neg))


def zscore(col):
    m = sum(col) / len(col)
    v = sum((c - m) ** 2 for c in col) / max(1, len(col) - 1)
    return [(c - m) / (math.sqrt(v) if v > 0 else 1.0) for c in col]


def logistic(X, y, iters=200, lr=0.1):
    n, k = len(X), len(X[0])
    w, b = [0.0] * k, 0.0
    for _ in range(iters):
        gw, gb = [0.0] * k, 0.0
        for i in range(n):
            z = max(-30.0, min(30.0, b + sum(w[j] * X[i][j] for j in range(k))))
            e = 1.0 / (1.0 + math.exp(-z)) - y[i]
            gb += e
            for j in range(k):
                gw[j] += e * X[i][j]
        b -= lr * gb / n
        for j in range(k):
            w[j] -= lr * gw[j] / n
    return w, b


def fit(rows, ix, use):
    y = [rows[i]['y'] for i in ix]
    if len(set(y)) < 2:
        return float('nan')
    z = [zscore([rows[i][c] for i in ix]) for c in use]
    X = [[z[j][t] for j in range(len(use))] for t in range(len(ix))]
    w, b = logistic(X, y)
    return auc(y, [b + sum(w[j] * r[j] for j in range(len(w))) for r in X])


def loocv(rows, use):
    """표준화도 학습 폴드 안에서만 계산한다 (9번과 동일)."""
    n = len(rows)
    y = [r['y'] for r in rows]
    oof = [0.0] * n
    for h in range(n):
        tr = [i for i in range(n) if i != h]
        mu, sd = [], []
        for c in use:
            v = [rows[i][c] for i in tr]
            m = sum(v) / len(v)
            mu.append(m)
            sd.append(math.sqrt(sum((x - m) ** 2 for x in v) / max(1, len(v) - 1)) or 1.0)
        X = [[(rows[i][c] - mu[j]) / sd[j] for j, c in enumerate(use)] for i in tr]
        w, b = logistic(X, [y[i] for i in tr])
        xh = [(rows[h][c] - mu[j]) / sd[j] for j, c in enumerate(use)]
        oof[h] = b + sum(w[j] * xh[j] for j in range(len(w)))
    return auc(y, oof)


def increment(rows, base_use, label):
    n = len(rows)
    idx = list(range(n))
    full_use = base_use + ['mal']
    base, full = fit(rows, idx, base_use), fit(rows, idx, full_use)
    rng = random.Random(SEED)
    boots = []
    for _ in range(N_BOOT):
        ix = [rng.randrange(n) for _ in range(n)]
        if len(set(rows[i]['y'] for i in ix)) < 2:
            continue
        boots.append(fit(rows, ix, full_use) - fit(rows, ix, base_use))
    s = sorted(v for v in boots if v == v)
    lo, hi = s[int(0.025 * len(s))], s[int(0.975 * len(s))]
    base_cv, full_cv = loocv(rows, base_use), loocv(rows, full_use)
    d, d_cv = full - base, full_cv - base_cv
    print(f'  {label:<28} 기준 {base:.3f} → +Malassezia {full:.3f}  증가분 {d:+.3f} '
          f'[{lo:+.3f}, {hi:+.3f}]  · 교차검증 {base_cv:.3f} → {full_cv:.3f} ({d_cv:+.3f})')
    return {'base_vars': base_use, 'auc_base': round(base, 3), 'auc_plus_malassezia': round(full, 3),
            'delta': round(d, 3), 'ci': [round(lo, 3), round(hi, 3)], 'crosses_zero': lo <= 0 <= hi,
            'auc_base_loocv': round(base_cv, 3), 'auc_plus_malassezia_loocv': round(full_cv, 3),
            'delta_loocv': round(d_cv, 3), 'loocv_agrees': (d_cv > 0.0) == (not (lo <= 0 <= hi))}


def rho_block(rows):
    r, p = spearman([x['sebum'] for x in rows], [x['mal'] for x in rows])
    return {'n': len(rows), 'rho': round(r, 3), 'p': round(p, 4)}


print('=' * 88)
print('피지 때문 아닌가 — Malassezia 착시 점검 ⑤')
print('=' * 88)
out = {'seed': SEED, 'n_boot': N_BOOT,
       'question': '민감군의 Malassezia 감소가 피지 차이로 설명되는가',
       'limitation': '피지는 한 번·볼 한 곳에서 잰 기기값. 01 노화는 다른 연구라 재현이 아니라 일관성 확인.'}

# (1)(2)(3) 03 민감 42명
r3 = load('03_sensitive_ITS', 'sebumeter_sebum_rigth')
S = [r['sebum'] for r in r3 if r['y'] == 1]
N = [r['sebum'] for r in r3 if r['y'] == 0]
a, p = auc_p(S, N)
lo, hi = auc_ci(S, N)
out['cohort03'] = {'n': len(r3), 'sensitive': len(S), 'non_sensitive': len(N),
                   'site': '오른쪽 볼', 'sebum_col': 'sebumeter_sebum_rigth'}
out['sebum_by_group'] = {'median_sensitive': st.median(S), 'median_non_sensitive': st.median(N),
                         'auc_sensitive_higher': round(a, 3), 'ci': [round(lo, 3), round(hi, 3)],
                         'p_mann_whitney_normal_approx': round(p, 3)}
print(f'\n(1) 피지 — 민감 중앙 {st.median(S)} vs 비민감 {st.median(N)} · '
      f'민감군이 클 확률 AUC {a:.3f} [{lo:.3f}, {hi:.3f}] · p={p:.3f}')

out['sebum_vs_malassezia_03'] = {'all': rho_block(r3),
                                 'sensitive': rho_block([r for r in r3 if r['y'] == 1]),
                                 'non_sensitive': rho_block([r for r in r3 if r['y'] == 0])}
print('(2) 피지 ~ Malassezia (Spearman)')
for k, v in out['sebum_vs_malassezia_03'].items():
    print(f"  {k:<14} n={v['n']:<3} rho={v['rho']:+.3f}  p={v['p']:.4f}")

print('(3) 피지를 넣어도 Malassezia 가 더하나')
out['adjusted'] = {'sebum': increment(r3, ['sebum'], '(가) 피지만'),
                   'age_sebum': increment(r3, ['age', 'sebum'], '(다) 나이+피지')}

# (4) 01 노화 58명 교차 확인
r1 = load('01_aging_ITS', 'sebum')
out['cross_check_01_aging'] = {'all': rho_block(r1),
                               'cheek_only': rho_block([r for r in r1 if r['site'] == 'Cheek']),
                               'forehead_only': rho_block([r for r in r1 if r['site'] == 'Forehead']),
                               'note': '민감 라벨 없음. 피지와 Malassezia 가 같이 움직이는지만 본다. 볼만 따로 본 것은 03(오른쪽 볼)과 부위를 맞추기 위해서.'}
print('(4) 01 노화 — 피지 ~ Malassezia (Spearman)')
for k in ('all', 'cheek_only', 'forehead_only'):
    v = out['cross_check_01_aging'][k]
    print(f"  {k:<14} n={v['n']:<3} rho={v['rho']:+.3f}  p={v['p']:.4f}")

# (5) scipy 교차 확인 — 피지_즉석계산_20260930.py(삭제)의 역할을 흡수 (2026-09-30)
#     JSON 에는 넣지 않는다 — scipy 유무에 따라 결과 파일이 달라지면 안 되므로.
#     Malassezia 는 즉석 계산처럼 다른 출처(sample_info 의 g_g__Malassezia 열)에서 읽어 입력 경로도 대조한다.
try:
    from scipy import stats as sps
except ImportError:
    sps = None
print('\n(5) scipy 교차 확인 — JSON 에는 저장하지 않음')
if sps is None:
    print('  scipy 없음 → 건너뜀 (위 (1)~(4)와 JSON 에는 영향 없음)')
else:
    mal_si = {r['run']: f(r['g_g__Malassezia'])
              for r in read_tsv(os.path.join(RAW, '03_sensitive_ITS_sample_info.tsv'))}
    print(f"  Malassezia 출처 대조  genus_pct vs sample_info(g_g__Malassezia)  "
          f"최대 차이 {max(abs(r['mal'] - mal_si[r['run']]) for r in r3):g}")
    u = sps.mannwhitneyu(S, N, alternative='two-sided')
    print(f"  (1) AUC stats_util {a:.3f} · scipy {u.statistic / (len(S) * len(N)):.3f}   "
          f"p stats_util {p:.3f} · scipy {u.pvalue:.3f}  ← p 차이는 scipy 기본의 연속성 보정")
    for k, sub in (('all', r3), ('sensitive', [r for r in r3 if r['y'] == 1]),
                   ('non_sensitive', [r for r in r3 if r['y'] == 0])):
        v = out['sebum_vs_malassezia_03'][k]
        rho_s, p_s = sps.spearmanr([r['sebum'] for r in sub], [mal_si[r['run']] for r in sub])
        print(f"  (2) {k:<14} rho stats_util {v['rho']:+.3f} · scipy {rho_s:+.3f}   "
              f"p stats_util {v['p']:.4f} · scipy {p_s:.4f}")

print('\n읽는 법')
print('  (1) 두 군 피지가 같으면 → Malassezia 차이를 "피지가 적어서"로 설명할 수 없다')
print('  (3) 증가분 CI 가 0 을 제외하고 교차검증도 같은 방향 → 피지를 넣어도 Malassezia 가 더한다')
save_json(out, os.path.join(TAB, '피지_점검.json'))
print('\n저장: results/피지_점검.json')
