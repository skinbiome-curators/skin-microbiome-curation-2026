"""3~6번 스크립트가 공유하는 통계 함수. 외부 라이브러리를 쓰지 않는다."""
import csv, math, os, random, statistics as st

# ── OS 이식성 (2026-09-18 추가) ──────────────────────────────────────────────
# 맥에서 만든 결과를 Windows 에서 그대로 재현하기 위한 두 가지.
#   1) stdout 을 UTF-8 로 고정 — 한글 Windows 기본 cp949 는 em-dash(—) 를 못 쓴다
#   2) save_json — encoding·newline 을 못 박아 OS 와 무관하게 같은 바이트를 낸다
import sys as _sys, json as _json
try:
    _sys.stdout.reconfigure(encoding='utf-8')
    _sys.stderr.reconfigure(encoding='utf-8')
except Exception:
    pass

def save_json(obj, path):
    with open(path, 'w', encoding='utf-8', newline=chr(10)) as _fh:
        _json.dump(obj, _fh, ensure_ascii=False, indent=2)
# ────────────────────────────────────────────────────────────────────────────


BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB  = os.path.join(BASE, '분석표')
RAW  = os.path.join(BASE, '원자료')

CAND = ['Cutibacterium', 'Staphylococcus', 'Corynebacterium']

# 구강·타액 상재로 알려진 속. 후보 선정에서 제외 여부를 판정하기 위한 사전 목록이며,
# 이 자료에서 구강 유래임을 증명한 목록이 아니다.
ORAL = ['Streptococcus', 'Rothia', 'Actinomyces', 'Haemophilus', 'Neisseria',
        'Veillonella', 'Gemella', 'Granulicatella', 'Lautropia', 'Prevotella']


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def rank(v):
    idx = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and v[idx[j + 1]] == v[idx[i]]:
            j += 1
        rr = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[idx[k]] = rr
        i = j + 1
    return r


# ── t 분포 양측 p (2026-09-18) ──────────────────────────────────────────────
# 이전에는 정규근사(erf)를 썼는데 n 이 작을 때 p 를 과소평가한다.
#   n=42, rho=0.385 → 정규근사 0.0083 · t분포 0.0118 (29% 과소)
# 표준 라이브러리만으로 불완전 베타 함수를 구현해 t 분포를 쓴다. scipy 와 1e-13 이내 일치.
import math

def _betacf(a, b, x, itmax=300, eps=3e-14, fpmin=1e-300):
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < fpmin: d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, itmax + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin: d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin: c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin: d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin: c = fpmin
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < eps:
            break
    return h

def _betai(a, b, x):
    if x <= 0.0: return 0.0
    if x >= 1.0: return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                  + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b

def t_p2(t, df):
    """t 통계량의 양측 p값."""
    if df <= 0: return float('nan')
    return _betai(0.5 * df, 0.5, df / (df + t * t))


def spearman(x, y):
    """스피어만 상관계수와 양측 p값 (t 분포, df = n-2)."""
    rx, ry = rank(x), rank(y)
    n = len(x)
    mx, my = st.mean(rx), st.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    rho = num / den if den else 0.0
    t = rho * math.sqrt((n - 2) / max(1e-12, 1 - rho ** 2))
    return rho, t_p2(t, n - 2)


def auc_p(a, b):
    """Mann-Whitney U 기반 AUC와 정규근사 양측 p값. a가 양성군."""
    allv = sorted([(v, 1) for v in a] + [(v, 0) for v in b])
    ranked = []
    i = 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        rr = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranked.append((rr, allv[k][1]))
        i = j + 1
    n1, n2 = len(a), len(b)
    R = sum(r for r, g in ranked if g == 1)
    U = R - n1 * (n1 + 1) / 2
    mu = n1 * n2 / 2
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    z = (U - mu) / sd
    return U / (n1 * n2), 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))


def auc_ci(a, b, N=4000, seed=20260918):
    """사람 단위 부트스트랩 95% 신뢰구간."""
    rng = random.Random(seed)
    s = [auc_p([rng.choice(a) for _ in a], [rng.choice(b) for _ in b])[0] for _ in range(N)]
    s.sort()
    return s[int(.025 * N)], s[int(.975 * N)]


def bh(pvals):
    """Benjamini-Hochberg q값. 입력 순서대로 반환."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for rk, i in enumerate(reversed(order), start=1):
        val = pvals[i] * m / (m - rk + 1)
        prev = min(prev, val)
        q[i] = min(1.0, prev)
    return q


def load_595():
    """595명 속 상대풍부도(비율)와 메타데이터를 함께 반환."""
    gr = list(csv.DictReader(open(os.path.join(TAB, '595_genus_relabund.tsv'), encoding='utf-8'), delimiter='\t'))
    md = {r['metadata_id']: r for r in csv.DictReader(open(os.path.join(TAB, '595_metadata.tsv'), encoding='utf-8'), delimiter='\t')}
    qc = {r['metadata_id']: r for r in csv.DictReader(open(os.path.join(TAB, '595_qc_flags.tsv'), encoding='utf-8'), delimiter='\t')}
    ids = [r['metadata_id'] for r in gr]
    genera = [k for k in gr[0] if k not in ('metadata_id', 'mb_sample_id')]
    val = {g: [f(r.get(g)) or 0.0 for r in gr] for g in genera}
    return ids, genera, val, md, qc


def load_42():
    """42명 속 상대비율(%)과 시료 정보를 반환."""
    si = {r['run']: r for r in csv.DictReader(open(os.path.join(RAW, '03_sensitive_16S_sample_info.tsv'), encoding='utf-8'), delimiter='\t')}
    rd = csv.reader(open(os.path.join(RAW, '03_sensitive_16S_genus_pct.tsv'), encoding='utf-8'), delimiter='\t')
    hdr = next(rd)
    rows = [r for r in rd]
    col = {g: i for i, g in enumerate(hdr)}
    info = [si[r[0]] for r in rows]
    val = {g: [f(r[col[g]]) or 0.0 for r in rows] for g in hdr[1:]}
    return info, val
