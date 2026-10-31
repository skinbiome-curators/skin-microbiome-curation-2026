# -*- coding: utf-8 -*-
"""파일럿 검토 p1 — 순차 모집 블라인드 내부 파일럿(무리 = 런 블록)의 1종 오류 · 검정력 · 인원 분포 모의실험.

설계(고친 표본수 틀 10/31 판 기준)
  사용 · 무처리 1:1 · 기준 두 번(b1 = -1주, b2 = 0주) · 4주(post) · 변화량 Δ = post - (b1+b2)/2
  주 균 하나(m = 1) · α 0.05 양측 · 검정력 0.80 · δ = ln 2
  인원 규칙(틀 3장): 비중심 t 로 푼 군당 n_t → ⌈n_t / 0.864⌉ (분석 인원 · 모집 환산은 따로)
  확인: 07 2원 S_w 0.374 → 10 · 80% 상한 0.431 → 13 · 15계절 2원 0.992 → 58 (틀 6장과 같아야 함)

비교하는 설계
  F10 · F13 · F58         : 처음부터 고정(07 2원 · 07 80% 상한 · 15계절 2원)  ← 대안 (가)
  P{n1}-base              : 첫 무리 n1/군 → 기준 두 번 차이로 S_w 추정(군 무관) → σ̂Δ² = 1.5 Ŝ_w²
                            (대안 (나) '기준 시료만 먼저 읽기'도 쓸 수 있는 정보가 이것뿐 → 같은 통계 성질)
  P{n1}-lump              : 첫 무리 모든 사람의 Δ 를 군 무시하고 분산(Kieser–Friede 의 단순 블라인드 추정)
  P15-adj                 : lump - δ_plan²·n/(4(n-1)) (효과 몫을 빼는 보정 · 음수면 바닥)
  P15-lump-cap30          : 상한 30/군(예산 상한)
  공통: 하한 = n1(첫 무리는 되돌릴 수 없음) · 상한 58/군(15계절 고정 설계와 같은 크기) · 상한에 닿으면 상한까지 진행
최종 분석(무리 = 층)
  T   : OLS Δ ~ 군 + 무리 (t 검정)
  VE  : van Elteren 층별 Wilcoxon(정규 근사)  ← 순위 검정을 고를 때
  MWp : 층을 무시한 Mann–Whitney(정규 근사)   ← 층을 빼먹으면
  AN  : OLS post ~ 군 + 기준 평균 + 무리 (ANCOVA)
참 시나리오(잡음 · 무리 효과 · 표류 · 반응 차이) — 아래 TRUTH
출력: p1_결과.json · 로그(진행 표시) · 시간 · peak 메모리
"""
import csv, json, math, os, sys, time, tracemalloc
import numpy as np
from scipy import stats
from scipy.optimize import brentq
from joblib import Parallel, delayed

sys.stdout.reconfigure(encoding='utf-8')
T0 = time.perf_counter(); tracemalloc.start()
HERE = os.path.dirname(os.path.abspath(__file__))
RAW = r'C:\Users\109so\OneDrive\Desktop\학교\micro3\새주제_v11\data'
CAND = ['Cutibacterium', 'Staphylococcus', 'Corynebacterium']
EPS = 1e-5
DELTA = math.log(2)
ALPHA, POWER, ARE = 0.05, 0.80, 0.864
SEED = 20261031
NREP0 = int(os.environ.get('NREP0', 20000))   # δ = 0(1종 오류)
NREP1 = int(os.environ.get('NREP1', 5000))    # δ = ln 2(검정력)
SW = {'07': 0.37346782647158877, '07u': 0.431, 'mid': 0.65, '15': 0.9924}
NMAX = 58


# ── 인원 규칙(틀 3장 · t 분포 · ÷0.864) ─────────────────────────────────
def build_rule(nhi=400):
    # n 4 부터(설계의 첫 무리 ≥ 10 이라 더 작은 n 은 쓰이지 않음 · 아주 작은 자유도에서 scipy nct 가 NaN 을 냄)
    ns = np.arange(4, nhi + 1)
    vmax = np.empty(len(ns))
    for k, n in enumerate(ns):
        df = 2 * n - 2
        tc = stats.t.ppf(1 - ALPHA / 2, df)
        g = lambda ncp: float(np.nan_to_num(stats.nct.sf(tc, df, ncp), nan=1.0)) + float(np.nan_to_num(stats.nct.cdf(-tc, df, ncp), nan=0.0)) - POWER
        ncp = brentq(g, 0.01, 12)
        vmax[k] = n * DELTA ** 2 / (2 * ncp ** 2)        # 이 분산 이하면 n 으로 충분
    return ns, vmax

NS, VMAX = build_rule()
ZC = (stats.norm.ppf(1 - ALPHA / 2) + stats.norm.ppf(POWER)) ** 2

def n_rule(var):
    k = np.searchsorted(VMAX, var)                       # VMAX 는 n 에 따라 증가
    if k < len(NS):
        nt = NS[k]
    else:
        nt = math.ceil(2 * ZC * var / DELTA ** 2) + 1
    return math.ceil(nt / ARE)

CHECK = {k: n_rule(1.5 * SW[k] ** 2) for k in SW}


# ── 경험 잡음 풀(15계절 2원 잔차 · r2 와 같은 방법) ─────────────────────
def load(tag, site, tpkey):
    si = {r['run']: r for r in csv.DictReader(open(os.path.join(RAW, tag + '_sample_info.tsv'), encoding='utf-8-sig'), delimiter='\t')}
    rd = csv.reader(open(os.path.join(RAW, tag + '_genus_pct.tsv'), encoding='utf-8'), delimiter='\t')
    hdr = next(rd); G = hdr[1:]
    rows = [r for r in rd if r[0] in si and si[r[0]].get('body_site') == site and si[r[0]].get('group') != '음성대조']
    X = np.array([[float(x) if x not in ('', 'NA') else 0.0 for x in r[1:]] for r in rows])
    meta = [(si[r[0]]['subject'], si[r[0]][tpkey]) for r in rows]
    rest = X[:, [g not in CAND for g in G]].sum(1)
    return meta, np.log((X[:, G.index('Cutibacterium')] + EPS) / (rest + EPS))

m15, L15 = load('15_season_16S', 'Forehead', 'collection_date')
tp15 = sorted({t for _, t in m15})
d = {}
for (s, t), x in zip(m15, L15):
    d.setdefault(s, {})[t] = x
W = np.array([[d[s][t] for t in tp15] for s in sorted(d) if all(t in d[s] for t in tp15)])
R = (W - W.mean(1, keepdims=True) - W.mean(0, keepdims=True) + W.mean()).ravel()
POOL15 = R / R.std()                                     # 단위 SD 로 맞춘 풀(모양만)
POOL_INFO = {'n_subject': int(W.shape[0]), 'n_tp': int(W.shape[1]), 'excess_kurtosis': float(stats.kurtosis(POOL15)),
             'skew': float(stats.skew(POOL15))}


# ── 참 시나리오 ────────────────────────────────────────────────────────
# sw: 무리1 S_w · sw2_mult: 무리2 S_w 배수(계절) · f: 표류(기준 두 번이 공유하는 '상태' 몫) · emp: 경험 분포
# g_sd: 무리마다 공통 변화(계절 · 두 군 같음) · run_sd: 무리(런)마다 모든 측정에 더해지는 값 · d2_mult: 무리2 효과 배수
# tau: 'none' | 'norm'(사람마다 N(δ, δ/2)) | 'resp'(반응자 50% 2δ)
TRUTH = {
    'S07':        dict(sw=SW['07'], f=0.0),
    'S07u':       dict(sw=SW['07u'], f=0.0),
    'Smid':       dict(sw=SW['mid'], f=0.0),
    'S15':        dict(sw=SW['15'], f=0.0),
    'S15emp':     dict(sw=SW['15'], f=0.0, emp=True),
    'S07drift':   dict(sw=SW['07'], f=0.5),
    'Smiddrift':  dict(sw=SW['mid'], f=0.5),
    'S15drift':   dict(sw=SW['15'], f=0.5),
    'S07season':  dict(sw=SW['07'], f=0.0, g_sd=0.5, run_sd=1.0, sw2_mult=1.3),
    'Smidseason': dict(sw=SW['mid'], f=0.0, g_sd=0.5, run_sd=1.0, sw2_mult=1.3),
    'S15season':  dict(sw=SW['15'], f=0.0, g_sd=0.5, run_sd=1.0, sw2_mult=1.3),
    'S07tau':     dict(sw=SW['07'], f=0.0, tau='norm'),
    'S07resp':    dict(sw=SW['07'], f=0.0, tau='resp'),
    'Smidmod':    dict(sw=SW['mid'], f=0.0, g_sd=0.5, run_sd=1.0, d2_mult=0.5),
}
SIGMA_B = 1.0   # 사람 사이 SD(ICC 를 정할 뿐 · 변화량 분석에는 영향 없음)

DESIGNS = {
    'F10': dict(kind='fixed', n=10), 'F13': dict(kind='fixed', n=13), 'F58': dict(kind='fixed', n=58),
    'P10-base': dict(kind='pilot', n1=10, est='base'), 'P15-base': dict(kind='pilot', n1=15, est='base'), 'P20-base': dict(kind='pilot', n1=20, est='base'),
    'P10-lump': dict(kind='pilot', n1=10, est='lump'), 'P15-lump': dict(kind='pilot', n1=15, est='lump'), 'P20-lump': dict(kind='pilot', n1=20, est='lump'),
    'P15-adj': dict(kind='pilot', n1=15, est='adj'),
    'P15-lump-cap30': dict(kind='pilot', n1=15, est='lump', nmax=30),
}


def gen(rng, n, tr, cohort, delta):
    """한 무리: 군당 n 명 · 반환 arm, b1, b2, post (길이 2n)."""
    sw = tr['sw'] * (tr.get('sw2_mult', 1.0) if cohort == 1 else 1.0)
    f = tr.get('f', 0.0)
    se, ss = sw * math.sqrt(1 - f), sw * math.sqrt(f)
    N = 2 * n
    arm = np.repeat([0, 1], n)
    mu = rng.normal(0, SIGMA_B, N)
    run = rng.normal(0, tr.get('run_sd', 0.0))          # 이 무리(런)의 모든 측정에 같은 값
    g = rng.normal(0, tr.get('g_sd', 0.0))              # 이 무리의 공통 변화(계절 · post 에만)
    if tr.get('emp'):
        e = rng.choice(POOL15, (3, N)) * se
    else:
        e = rng.normal(0, 1, (3, N)) * se
    sb, sp = rng.normal(0, 1, N) * ss, rng.normal(0, 1, N) * ss
    b1 = mu + run + sb + e[0]
    b2 = mu + run + sb + e[1]
    dc = delta * (tr.get('d2_mult', 1.0) if cohort == 1 else 1.0)
    eff = np.zeros(N)
    if dc != 0:
        tau = tr.get('tau', 'none')
        if tau == 'norm':
            eff[arm == 1] = rng.normal(dc, dc / 2, n)
        elif tau == 'resp':
            eff[arm == 1] = np.where(rng.random(n) < 0.5, 2 * dc, 0.0)
        else:
            eff[arm == 1] = dc
    post = mu + run + g + sp + e[2] + eff
    return arm, b1, b2, post


def p_norm(z):
    return math.erfc(abs(z) / math.sqrt(2))


def tests(arm, coh, b1, b2, post):
    D = post - (b1 + b2) / 2
    bm = (b1 + b2) / 2
    cs = np.unique(coh)
    out = {}
    # T: Δ ~ 군 + 무리
    X = [np.ones(len(D)), arm.astype(float)] + [(coh == c).astype(float) for c in cs[1:]]
    X = np.column_stack(X)
    for key, y, XX in (('T', D, X), ('AN', post, np.column_stack([X, bm]))):
        beta, *_ = np.linalg.lstsq(XX, y, rcond=None)
        r = y - XX @ beta
        df = len(y) - XX.shape[1]
        s2 = r @ r / df
        v = s2 * np.linalg.inv(XX.T @ XX)[1, 1]
        out[key] = 2 * stats.t.sf(abs(beta[1] / math.sqrt(v)), df)
    # VE: 층별 Wilcoxon
    T = V = 0.0
    for c in cs:
        m = coh == c
        y, a = D[m], arm[m]
        Nn = len(y); n1 = a.sum(); n0 = Nn - n1
        rk = np.empty(Nn); rk[np.argsort(y)] = np.arange(1, Nn + 1)
        T += (rk[a == 1].sum() - n1 * (Nn + 1) / 2) / (Nn + 1)
        V += n1 * n0 / (12 * (Nn + 1))
    out['VE'] = p_norm(T / math.sqrt(V))
    # MWp: 층 무시
    Nn = len(D); n1 = arm.sum(); n0 = Nn - n1
    rk = np.empty(Nn); rk[np.argsort(D)] = np.arange(1, Nn + 1)
    out['MWp'] = p_norm((rk[arm == 1].sum() - n1 * (Nn + 1) / 2) / math.sqrt(n1 * n0 * (Nn + 1) / 12))
    return out


def run_cell(dname, tname, delta, nrep, seed):
    tm0 = time.perf_counter(); tracemalloc.start()
    rng = np.random.default_rng([SEED, seed])
    ds, tr = DESIGNS[dname], TRUTH[tname]
    sig = {k: 0 for k in ('T', 'AN', 'VE', 'MWp')}
    Ns, ratio, single, cap = [], [], 0, 0
    f = tr.get('f', 0.0)
    var_true = tr['sw'] ** 2 * (1.5 * (1 - f) + 2 * f)   # 무리1 대조군 Δ 의 참 분산(효과 · τ 제외)
    for _ in range(nrep):
        if ds['kind'] == 'fixed':
            arm, b1, b2, post = gen(rng, ds['n'], tr, 0, delta)
            coh = np.zeros(len(arm), int); Ntot = ds['n']
        else:
            n1, nmax = ds['n1'], ds.get('nmax', NMAX)
            a1, x1, y1, p1 = gen(rng, n1, tr, 0, delta)
            if ds['est'] == 'base':
                dd = y1 - x1
                V = 1.5 * dd.var(ddof=1) / 2
            else:
                D1 = p1 - (x1 + y1) / 2
                V = D1.var(ddof=1)
                if ds['est'] == 'adj':
                    nn = len(D1)
                    V = max(V - DELTA ** 2 * nn / (4 * (nn - 1)), 1e-6)
            ratio.append(V / var_true)
            nt = n_rule(V)
            if nt > nmax:
                cap += 1
            Ntot = min(max(nt, n1), nmax)
            n2 = Ntot - n1
            if n2 > 0:
                a2, x2, y2, p2 = gen(rng, n2, tr, 1, delta)
                arm = np.concatenate([a1, a2]); b1 = np.concatenate([x1, x2]); b2 = np.concatenate([y1, y2]); post = np.concatenate([p1, p2])
                coh = np.concatenate([np.zeros(2 * n1, int), np.ones(2 * n2, int)])
            else:
                arm, b1, b2, post, coh = a1, x1, y1, p1, np.zeros(2 * n1, int)
                single += 1
        Ns.append(Ntot)
        ps = tests(arm, coh, b1, b2, post)
        for k in sig:
            sig[k] += ps[k] < ALPHA
    Ns = np.array(Ns)
    res = {'design': dname, 'truth': tname, 'delta': delta, 'nrep': nrep,
           'rej': {k: v / nrep for k, v in sig.items()},
           'N_per_arm_mean': float(Ns.mean()), 'N_q10_50_90': [int(x) for x in np.percentile(Ns, [10, 50, 90])],
           'P_single_cohort': single / nrep, 'P_cap_hit': cap / nrep,
           'Vhat_over_true_mean': float(np.mean(ratio)) if ratio else None,
           'sec': time.perf_counter() - tm0, 'peak_MB_worker': tracemalloc.get_traced_memory()[1] / 1e6}
    tracemalloc.stop()
    return res


if __name__ == '__main__':
    print('인원 규칙 확인(틀 6장: 07 10 · 07u 13 · 15 58):', CHECK, flush=True)
    assert CHECK['07'] == 10 and CHECK['07u'] == 13 and CHECK['15'] == 58, CHECK
    print('15계절 경험 풀:', POOL_INFO, flush=True)
    cells = []
    i = 0
    for dn in DESIGNS:
        for tn in TRUTH:
            for dl, nr in ((0.0, NREP0), (DELTA, NREP1)):
                cells.append((dn, tn, dl, nr, i)); i += 1
    print(f'칸 {len(cells)}개 · 반복 δ=0 {NREP0} · δ=ln2 {NREP1} · 시드 {SEED} · joblib 22코어', flush=True)
    res = Parallel(n_jobs=22, verbose=10)(delayed(run_cell)(*c) for c in cells)
    out = {'seed': SEED, 'nrep0': NREP0, 'nrep1': NREP1, 'delta': DELTA, 'alpha': ALPHA, 'power': POWER, 'ARE': ARE,
           'rule_check': CHECK, 'pool15': POOL_INFO, 'truth': TRUTH, 'designs': DESIGNS, 'results': res,
           'mem_time': {'sec_total': time.perf_counter() - T0, 'peak_MB_main': tracemalloc.get_traced_memory()[1] / 1e6,
                        'peak_MB_worker_max': max(r['peak_MB_worker'] for r in res)}}
    json.dump(out, open(os.path.join(HERE, 'p1_결과.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1, default=float)
    print(out['mem_time'], flush=True)
