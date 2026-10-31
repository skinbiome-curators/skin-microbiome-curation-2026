# -*- coding: utf-8 -*-
"""방법 검토 r2 — 표본수 틀이 주장하는 검정력(80%)이 정해진 n 에서 실제로 나오는지 모의실험.

설계(틀 6장 · 견본 B-1): 사용군 · 무처리군 1:1 · 기준 두 번(b1, b2) · 4주(post) · 변화량 Δ = post − (b1+b2)/2
E3 = 두 군의 Δ 를 Mann–Whitney(Wilcoxon 순위합)로 비교 · 양측 α/m (m = 3 이면 Holm 첫 단계 = Bonferroni)
효과: 사용군 post 에 δ = ln 2 더함(대조보다 ×2 — 후보/나머지 비)

잡음 모형
  normal   : 측정마다 N(0, S_w²) 독립
  emp07    : 07 2W→4W 변화량 d_i(평균 뺌)를 ±d/√2 로 측정 오차 풀을 만들어 S_w 로 맞춤(꼬리 모양 유지)
  emp15    : 15계절 2원 잔차 풀을 S_w 로 맞춤
  drift    : S_w² 의 절반은 '상태'(1주 간격 기준 두 번이 공유 · 4주 뒤는 새로 뽑음) · 절반은 측정 오차
  het      : 사람마다 효과가 다름 — N(δ, (δ/2)²) · 또는 반응자 50%(2δ) / 비반응자 50%(0)
  Sw_unc   : 참 S_w 를 07 추정의 불확실성(자유도 24)에서 뽑음 → 기대 검정력(assurance)

추가
  - 모집 경로: 군당 15명 모집 · 시료 성공 s(독립) · 사람 탈락 q → 분석 n 이 우연히 정해짐
      strict = 세 시료 모두 성공해야 분석 / lenient = post + 기준 하나 이상
  - 결정표 5번(근거로 쓸 수 있음)까지: E3 유의 AND E2(같은 시험 잡음) 통과 확률 — E2 숫자 기준은 빈칸(5-6 #1)이라
      후보 두 개를 시험: (a) 사용군 Δ_i > CR̂ 인 사람 k 의 정확 95% CI 하한 > 2.5% · (b) 사용군 평균 Δ > CR̂/√n
      CR̂ = 2.77·Ŝ_w, Ŝ_w = 두 기준값(b1, b2)의 2원 S_w(모든 사람)
  - 세 균 함께(07 경험 분포에서 사람 단위로 세 균을 같이 뽑음 · Holm) — 각 균 유의 · 셋 모두 유의 확률
  - 1종 오류(δ = 0)
출력: r2_모의실험.json · 화면 진행 표시
"""
import csv, json, math, os, sys, time, tracemalloc
import numpy as np
from scipy import stats
from joblib import Parallel, delayed

sys.stdout.reconfigure(encoding='utf-8')
T0 = time.perf_counter(); tracemalloc.start()
HERE = os.path.dirname(os.path.abspath(__file__))
RAW = r'C:\Users\109so\OneDrive\Desktop\학교\micro3\새주제_v11\data'
CAND = ['Cutibacterium', 'Staphylococcus', 'Corynebacterium']
EPS = 1e-5
DELTA = math.log(2)
SEED = 20261031
NREP = int(os.environ.get('NREP', 10000))
CR_K = 1.96 * math.sqrt(2)
SW = {'07': {'Cutibacterium': 0.37346782647158877, 'Staphylococcus': 0.6558415395505228, 'Corynebacterium': 0.4777299998200134},
      '15': {'Cutibacterium': 0.9924, 'Staphylococcus': 0.9797, 'Corynebacterium': 0.8},
      '07_1': {'Cutibacterium': 0.4424468145138683, 'Staphylococcus': 0.713925310862275, 'Corynebacterium': 0.4972475900558258}}


# ── 경험 잡음 풀 ────────────────────────────────────────
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

def wide(meta, v, tps):
    d = {}
    for (s, t), x in zip(meta, v):
        d.setdefault(s, {})[t] = x
    subj = sorted(s for s in d if all(t in d[s] for t in tps))
    return np.array([[d[s][t] for t in tps] for s in subj])

m7, L7 = load('07_skincare_16S', 'Cheek', 'timepoint')
m15, L15 = load('15_season_16S', 'Forehead', 'collection_date')
tp15 = sorted({t for _, t in m15})
POOL = {'07': {}, '15': {}}
D07 = {}
for c in CAND:
    W = wide(m7, L7[c], ['2W', '4W']); d = W[:, 1] - W[:, 0]; d = d - d.mean()
    D07[c] = d
    e = np.concatenate([d, -d]) / math.sqrt(2)
    POOL['07'][c] = e / e.std() * SW['07'][c]
    W = wide(m15, L15[c], tp15)
    R = (W - W.mean(1, keepdims=True) - W.mean(0, keepdims=True) + W.mean()).ravel()
    POOL['15'][c] = R / R.std() * SW['15'][c]
D07_joint = np.column_stack([D07[c] for c in CAND])            # 사람 × 3균(같은 사람에서 함께)


# ── 한 시나리오 ─────────────────────────────────────────
def draw_err(rng, model, src, genus, size, sw):
    if model in ('normal', 'het_norm', 'het_resp', 'Sw_unc', 'recruit_strict', 'recruit_lenient', 'drift'):
        return rng.normal(0, sw, size)
    if model == 'emp07':
        return rng.choice(POOL['07'][genus], size) * (sw / SW['07'][genus])
    if model == 'emp15':
        return rng.choice(POOL['15'][genus], size) * (sw / SW['15'][genus])
    raise ValueError(model)

def mw_p(x, y, exact):
    return stats.mannwhitneyu(x, y, alternative='two-sided', method='exact' if exact else 'asymptotic', axis=1).pvalue

def e2_pass(b1, b2, delta_use, n_all):
    """같은 시험 잡음 E2 후보 두 가지. b1, b2: (rep, N) 모든 사람 기준값. delta_use: (rep, n) 사용군 Δ."""
    dd = b2 - b1
    dd = dd - dd.mean(1, keepdims=True)
    sw_hat = np.sqrt((dd ** 2).sum(1) / (2 * (dd.shape[1] - 1)))      # 2원(k=2) S_w
    cr = CR_K * sw_hat
    n = delta_use.shape[1]
    k = (delta_use > cr[:, None]).sum(1)
    lo = np.where(k > 0, stats.beta.ppf(0.025, np.maximum(k, 1), n - k + 1), 0.0)
    a = lo > 0.025
    b = delta_use.mean(1) > cr / math.sqrt(n)
    return a, b, k

def scenario(name, model, src, genus, n, m, delta=DELTA, nrep=NREP, seed=0, exact=None, s=0.95, q=0.10, T=3, n_recruit=None):
    rng = np.random.default_rng([SEED, seed])
    sw = SW[src][genus]
    alpha = 0.05 / m
    if exact is None:
        exact = n <= 30
    t0 = time.perf_counter()
    if model in ('recruit_strict', 'recruit_lenient'):
        # 모집 n_recruit · 시료 성공 · 탈락 → 분석 n 이 반복마다 다름 → 반복마다 따로 검정
        sig = 0; ns = []
        for r in range(nrep):
            arms = []
            for arm in (0, 1):
                N = n_recruit
                stay = rng.random(N) > q
                ok = rng.random((N, 3)) < s                  # b1, b2, post
                if model == 'recruit_strict':
                    use = stay & ok.all(1)
                    base = None
                else:
                    use = stay & ok[:, 2] & (ok[:, 0] | ok[:, 1])
                b = rng.normal(0, sw, (N, 2)); p = rng.normal(0, sw, N) + (delta if arm == 1 else 0)
                bm = np.where(ok[:, 0] & ok[:, 1], b.mean(1), np.where(ok[:, 0], b[:, 0], b[:, 1]))
                arms.append((p - bm)[use])
            ns.append((len(arms[0]), len(arms[1])))
            if min(len(arms[0]), len(arms[1])) < 2:
                continue
            pv = stats.mannwhitneyu(arms[1], arms[0], alternative='two-sided',
                                    method='exact' if max(len(arms[0]), len(arms[1])) <= 30 else 'asymptotic').pvalue
            sig += pv < alpha
        ns = np.array(ns)
        return {'name': name, 'power_MW': sig / nrep, 'n_analysed_mean': ns.mean(0).tolist(),
                'P(n_analysed>=target_each_arm)': float(((ns >= n).all(1)).mean()), 'sec': time.perf_counter() - t0}
    # 기본 경로 — 행렬로 한 번에
    sw_true = sw
    if model == 'Sw_unc':
        df = 24
        sw_true = sw * np.sqrt(df / rng.chisquare(df, nrep))[:, None]
    shape = (nrep, n)
    def meas(size):
        if model == 'Sw_unc':
            return rng.normal(0, 1, size) * sw_true
        return draw_err(rng, model, src, genus, size, sw)
    res = {}
    if model == 'drift':
        se = sw / math.sqrt(2)                                # 측정 오차 절반 · 상태 절반
        out = []
        for arm in (0, 1):
            ub = rng.normal(0, se, shape); up = rng.normal(0, se, shape)
            b1 = ub + rng.normal(0, se, shape); b2 = ub + rng.normal(0, se, shape)
            p = up + rng.normal(0, se, shape) + (delta if arm == 1 else 0)
            out.append((b1, b2, p))
    else:
        out = []
        for arm in (0, 1):
            b1, b2, p = meas(shape), meas(shape), meas(shape)
            if arm == 1:
                if model == 'het_norm':
                    p = p + rng.normal(delta, delta / 2, shape)
                elif model == 'het_resp':
                    p = p + np.where(rng.random(shape) < 0.5, 2 * delta, 0.0)
                else:
                    p = p + delta
            out.append((b1, b2, p))
    dC = out[0][2] - (out[0][0] + out[0][1]) / 2
    dU = out[1][2] - (out[1][0] + out[1][1]) / 2
    pmw = mw_p(dU, dC, exact)
    pt = stats.ttest_ind(dU, dC, axis=1, equal_var=False).pvalue
    sigE3 = pmw < alpha
    b1all = np.concatenate([out[0][0], out[1][0]], 1); b2all = np.concatenate([out[0][1], out[1][1]], 1)
    a, b, k = e2_pass(b1all, b2all, dU, 2 * n)
    res.update(name=name, model=model, src=src, genus=genus, n_per_arm=n, m=m, alpha=alpha, delta=delta,
               power_MW=float(sigE3.mean()), power_Welch_t=float((pt < alpha).mean()),
               P_E3_and_E2a=float((sigE3 & a).mean()), P_E3_and_E2b=float((sigE3 & b).mean()),
               P_E2a=float(a.mean()), P_E2b=float(b.mean()), k_mean=float(k.mean()),
               sec=time.perf_counter() - t0)
    return res


def joint3(name, n, nrep=NREP, seed=0, delta=DELTA):
    """세 균 함께 · 07 경험 분포(사람 단위로 세 균 같이) · Holm(양측) — 각 균에 δ = ln2"""
    rng = np.random.default_rng([SEED, seed])
    t0 = time.perf_counter()
    def draw(shape):
        idx = rng.integers(0, D07_joint.shape[0], shape)
        sg = rng.choice([-1, 1], shape)
        return D07_joint[idx] * sg[..., None] / math.sqrt(2)     # 측정 하나의 오차(세 균 상관 유지)
    sh = (nrep, n)
    dd = []
    for arm in (0, 1):
        b1, b2, p = draw(sh), draw(sh), draw(sh)
        if arm == 1:
            p = p + delta
        dd.append(p - (b1 + b2) / 2)
    P = np.stack([stats.mannwhitneyu(dd[1][..., g], dd[0][..., g], alternative='two-sided', method='exact', axis=1).pvalue
                  for g in range(3)], 1)
    # Holm
    order = np.argsort(P, 1); Ps = np.take_along_axis(P, order, 1)
    adj = np.maximum.accumulate(np.minimum(1, Ps * np.array([3, 2, 1])), 1)
    rej_sorted = adj < 0.05
    rej = np.zeros_like(rej_sorted); np.put_along_axis(rej, order, rej_sorted, 1)
    return {'name': name, 'n_per_arm': n, 'P_sig_each': dict(zip(CAND, rej.mean(0).tolist())),
            'P_all3': float(rej.all(1).mean()), 'sec': time.perf_counter() - t0}


SC = [
    # name, model, src, genus, n, m, kwargs
    ('A 07정규 Cuti n11 m3', 'normal', '07', 'Cutibacterium', 11, 3, {}),
    ('A0 07정규 Cuti n11 m3 δ=0(1종오류)', 'normal', '07', 'Cutibacterium', 11, 3, {'delta': 0.0}),
    ('B 07정규 Cuti n8 m1', 'normal', '07', 'Cutibacterium', 8, 1, {}),
    ('C 15정규 Cuti n75 m3', 'normal', '15', 'Cutibacterium', 75, 3, {}),
    ('D 15정규 Cuti n56 m1', 'normal', '15', 'Cutibacterium', 56, 1, {}),
    ('E 07경험 Cuti n11 m3', 'emp07', '07', 'Cutibacterium', 11, 3, {}),
    ('F 15경험 Cuti n75 m3', 'emp15', '15', 'Cutibacterium', 75, 3, {}),
    ('G1 07 사람마다 효과 N(δ,δ/2) n11', 'het_norm', '07', 'Cutibacterium', 11, 3, {}),
    ('G2 07 반응자 50% n11', 'het_resp', '07', 'Cutibacterium', 11, 3, {}),
    ('H 07 S_w 불확실(df24) n11', 'Sw_unc', '07', 'Cutibacterium', 11, 3, {}),
    ('I 07 표류(기준 두 번이 상태 공유) n11', 'drift', '07', 'Cutibacterium', 11, 3, {}),
    ('I15 15 표류 n75', 'drift', '15', 'Cutibacterium', 75, 3, {}),
    ('J1 모집15 strict s.95 q.1', 'recruit_strict', '07', 'Cutibacterium', 11, 3, {'n_recruit': 15, 'nrep': 4000}),
    ('J2 모집15 lenient s.95 q.1', 'recruit_lenient', '07', 'Cutibacterium', 11, 3, {'n_recruit': 15, 'nrep': 4000}),
    ('K Staph 07정규 n11(Cuti 기준 인원)', 'normal', '07', 'Staphylococcus', 11, 3, {}),
    ('K2 Coryne 07정규 n11', 'normal', '07', 'Corynebacterium', 11, 3, {}),
    ('I0 07 표류 δ=0(E2 거짓 통과)', 'drift', '07', 'Cutibacterium', 11, 3, {'delta': 0.0}),
    ('A1 07 1원 S_w 가 참이면 n11', 'normal', '07_1', 'Cutibacterium', 11, 3, {}),
    ('L 07정규 Cuti n11 m3 · 15 잡음이 참이면', 'normal', '15', 'Cutibacterium', 11, 3, {}),
]

if __name__ == '__main__':
    print(f'모의실험 {len(SC) + 1}개 시나리오 · 반복 {NREP} · 시드 {SEED} · joblib 22코어')
    jobs = [delayed(scenario)(nm, mo, sr, ge, n, m, seed=i, **kw) for i, (nm, mo, sr, ge, n, m, kw) in enumerate(SC)]
    jobs.append(delayed(joint3)('M 세 균 함께 07경험 n11 Holm', 11, seed=99))
    jobs.append(delayed(joint3)('M2 세 균 함께 07경험 n33 Holm', 33, seed=98))
    res = Parallel(n_jobs=22, verbose=10)(jobs)
    out = {'seed': SEED, 'nrep': NREP, 'delta': DELTA, 'results': res,
           'mem_time': {'sec': time.perf_counter() - T0, 'peak_MB_main': tracemalloc.get_traced_memory()[1] / 1e6}}
    json.dump(out, open(os.path.join(HERE, 'r2_모의실험.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1, default=float)
    for r in res:
        print(json.dumps(r, ensure_ascii=False, default=float))
    print(out['mem_time'])
