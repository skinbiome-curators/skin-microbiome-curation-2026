#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
24 효능 판정 시범 — 07 화장품(변화) + 15계절(잡음 · 규칙 시험)  (2026-10-17 작성 · 팀 합의)

규칙: ../경진대회/효능판정기준_초안_261017.md v2 — ★ 결정 1–8 · 2장 판정 항목 · 3장 결정표 · 5-2 사양.

목적(5-2):
  (1) 품질 점검(E0)이 07 에서 '시점과 배치가 겹침'을 자동으로 잡아 결정표 1번('아직 안 됨')을 내는지
  (2) 결정표를 무처리 자료(15계절 계절 쌍 6개)에 적용하면 어떤 등급이 나오는지
      — 대조가 없으므로 '근거로 쓸 수 있음'은 구조상 나올 수 없다. 볼 것은 '조건부'가 얼마나 나오나다.
  E1 · E2 · E0-③ 수치는 **참고**다 — 07 의 등급은 규칙(E0)으로 정해진다(결정 5 · 7).

입력
  data/07_skincare_16S_{sample_info,genus_pct}.tsv   볼 · 25명 × 0W · 2W · 4W · 음성대조 없음
  data/15_season_16S_{sample_info,genus_pct}.tsv     이마 · 24명 × 4계절(2014-01 · 04 · 08 · 11) · 음성대조 4(모두 2014-01)
  FASTQ 머리줄(ENA) → results/효능판정_시범_런번호.tsv 에 저장하고, 파일이 있으면 다시 받지 않는다.
     07: Illumina 머리줄 `기기:런:플로셀:레인:` 에서 시퀀싱 런을 읽는다(파일마다 첫 3읽기).
     15: 머리줄이 `시료이름|번호|…` 로 바뀌어 있어 런 정보가 없다 → '정보 없음'.

방법(문서 5-2 번호)
  #2 E0-② (a) 두 시점의 시퀀싱 런 — crossed(두 시점이 쓴 모든 런에 두 시점 시료가 함께 있음) /
              partial(일부 런만 공유) / disjoint(공유 없음) / unknown. crossed 가 아니면 '겹침'.
          (b) 한 시점에만 몰린 속 — 검출 문턱 1 · 0.1 · 0 % 격자, 한 시점 ≥ 80 % · 다른 모든 시점 ≤ 20 %.
              E0 상태 판정에는 1 % 문턱만 쓰고 나머지는 보고만 한다. 문턱별 목록은 포함 관계가 아닐 수 있다.
     E0-④ 음성대조 — 비교하는 두 시점 모두에 있어야 한다(시점별로 센다).
  #3 E0-③ 문턱마다(+ 세 목록의 합집합) 표시된 속을 로그비 분모(rest)에서 뺀 판 — 참고
  #4 로그비 = log((후보 + EPS) / (rest + EPS)), rest = 후보 셋이 아닌 속의 합(Unassigned 포함),
     EPS = 1e-5 (% 척도) — code/12 와 같다. 보조: 후보끼리의 로그비(분모 없음)
  #5 E1 Wilcoxon 부호순위(scipy · method 명시) · Hodges–Lehmann 추정치 · 정확 95 % CI(부호순위 귀무분포 직접 계산)
     · Holm(3속). 주장 방향 검정과 양측 검정을 함께 낸다 — 양측으로 유의하고 방향이 반대면 '주장과 반대 방향'.
     주장 방향: Cutibacterium · Staphylococcus = 증가(원 논문) · Corynebacterium = 양쪽(원 논문 주장 없음)
  #6 E2 15계절 S_w(1원 = 판정용 / 2원 = 참고) → CR = 1.96·√2·S_w (≈ 2.77·S_w, Vaz 2013).
     07 안의 대리 잡음(2W ↔ 4W, 1원 · 평균 이동을 뺀 2원)도 함께. 'CR 을 넘은 사람 k / n (95 % 정확 CI)'.
  #7 15계절 계절 쌍 6개에 결정표 전체를 적용

⚠️ 같이 말할 것
  1. 07 의 0W · 2W · 4W 는 서로 다른 MiSeq 런이다(머리줄). 추출 · PCR 배치는 모른다.
  2. 주 비교 0W → 4W 는 시점별 중앙값을 본 뒤 정했다(부분 사후). 방법 검토가 일부 수치를 미리 봤다(문서 ★ v1 → v2).
  3. 15계절 잡음은 다른 부위 · 영역 · 간격이다 → E2 는 참고다(결정 7). 15계절 규칙 시험의 E2 는 같은 자료로 잰 잡음이다.
  4. 대조가 있는 경우의 결정표(2–6번)는 실제 자료가 없어 합성 사례 자가검사로만 확인한다.
  5. E0-①(심도 · 오염 의심 표시 · 이름표 신뢰도)은 이 시범에서 계산하지 않는다(5-2 범위 밖).
  6. 같은 시험 안의 잡음으로 E2 를 판정할 숫자 기준은 아직 없다(문서 «빈칸») → 이 시범에서는 e2_gate = None.
  7. 부호순위의 정확 CI 는 0 · 동점이 없다고 가정한다. 0 · 동점이 있으면 JSON 에 표시한다(이번 자료는 없음).

변경 기록
  2026-10-17 첫 작성.
  2026-10-17 독립 검수 반영 —
    · 런 겹침 판정: 공유 런이 하나라도 있으면 통과 → crossed 일 때만 통과(partial 도 '겹침')
    · E1: 양측 p 도 계산해 결정표 8번 '주장과 반대 방향'에 닿게 함 · wilcoxon method 명시 · 0 · 동점 표시
    · 음성대조: 코호트 전체 → 시점별(두 시점 모두에 있어야 함)
    · 15계절 규칙 시험: 'passed' → 'structurally_guaranteed' + 조건부 비율 보고
    · E0-③: 문턱별 목록의 포함 관계 · 교집합 기록 + 합집합을 뺀 판 추가
    · E2 양쪽 주장: max(up, down) 대 2.5 % → up + down 대 5 %
    · 머리줄: 첫 1읽기 → 첫 3읽기 · 런 정보가 없는 행은 status 'ok_no_runfield'
    · 07 대리 잡음: 평균 이동을 뺀 2원 S_w 추가 · E4(기준 시점 반복측정) 기록
    · JSON 을 그림보다 먼저 저장 · 시점 정렬을 숫자 기준으로 · 쓰지 않던 변수 삭제 · 그림 범례 · 이름 고침

결과: results/효능판정_시범.json · 효능판정_시범_사람별변화.tsv(행 = 사람, 열 = 판|비교|균 — matrix) ·
      효능판정_시범_런번호.tsv · 효능판정_시범_그림.{png,pdf} · 효능판정_시범_그림데이터.tsv (TSV 는 기존 results 처럼 BOM 없음)
의존성: scipy · matplotlib (conda bio_proj) + ./stats_util.py · ./common42_261007.py(mem_time_line)
"""
import csv, datetime, hashlib, math, os, re, statistics as st, sys, time, tracemalloc, zlib
import urllib.request
from concurrent.futures import ThreadPoolExecutor

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

T0 = time.perf_counter()
tracemalloc.start()

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from stats_util import save_json, TAB, RAW, CAND           # noqa: E402
from common42_261007 import mem_time_line                  # noqa: E402

from scipy import stats                                    # noqa: E402

EPS = 1e-5                     # 로그비의 0 처리 — code/12 · code/4 와 같은 값(% 척도)
THRESH = (1.0, 0.1, 0.0)       # E0-② 검출 문턱 격자(%) — 0.0 은 '0 초과'
PRIMARY_THR = 1.0              # E0 상태 판정에 쓰는 문턱
HIGH, LOW = 0.8, 0.2           # 한 시점 ≥ 80 % · 다른 모든 시점 ≤ 20 %
CR_K = 1.96 * math.sqrt(2)     # 반복성 계수 CR = 1.96·√2·S_w (Vaz 2013 의 2.77)
ALPHA = 0.05
N_HEAD = 3                     # 파일마다 확인하는 머리줄 수
CLAIM = {'Cutibacterium': 'greater', 'Staphylococcus': 'greater', 'Corynebacterium': 'two-sided'}
PAIRS = (('Cutibacterium', 'Staphylococcus'), ('Cutibacterium', 'Corynebacterium'),
         ('Staphylococcus', 'Corynebacterium'))
PRIMARY = ('0W', '4W')         # 07 주 비교(부분 사후 — 문서 5-2 #5)

RUNS_TSV = os.path.join(TAB, '효능판정_시범_런번호.tsv')
OUT_JSON = os.path.join(TAB, '효능판정_시범.json')
OUT_DELTA = os.path.join(TAB, '효능판정_시범_사람별변화.tsv')
OUT_FIG = os.path.join(TAB, '효능판정_시범_그림')
OUT_FIGDATA = os.path.join(TAB, '효능판정_시범_그림데이터.tsv')

COHORTS = {
    '07': dict(tag='07_skincare_16S', site='Cheek', tp_key='timepoint', project='PRJEB44885'),
    '15': dict(tag='15_season_16S', site='Forehead', tp_key='collection_date', project='PRJNA390040'),
}


def tp_sorted(tps):
    """시점 정렬 — '10W' 가 '2W' 보다 앞에 오지 않게 숫자 W 는 숫자로, 나머지(날짜)는 문자열로."""
    def key(t):
        m = re.match(r'^(\d+)W$', t)
        return (0, int(m.group(1)), '') if m else (1, 0, t)
    return sorted(tps, key=key)


def cname(t1, t2):
    return f'{t1}→{t2}'


# ── 입력 ────────────────────────────────────────────────────────────────────
def read_tsv(p):
    with open(p, encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def input_fingerprints():
    """입력 지문 — data/SHA256SUMS 와 같은지 본다(재현 패키지 G)."""
    listed = {}
    with open(os.path.join(RAW, 'SHA256SUMS'), encoding='utf-8') as fh:
        for line in fh:
            parts = line.split()
            if len(parts) == 2:
                listed[parts[1]] = parts[0]
    out = {}
    for c in COHORTS.values():
        for kind in ('sample_info', 'genus_pct'):
            name = f"{c['tag']}_{kind}.tsv"
            got = sha256(os.path.join(RAW, name))
            out[name] = {'sha256': got, 'matches_SHA256SUMS': listed.get(name) == got}
    return out


def load(key):
    """시료(사람 · 시점) 목록과 속별 % 행렬. 음성대조는 빼고 시점별 개수만 센다."""
    c = COHORTS[key]
    si_rows = read_tsv(os.path.join(RAW, c['tag'] + '_sample_info.tsv'))
    neg_by_tp = {}
    for r in si_rows:
        if r.get('group') == '음성대조':
            neg_by_tp[r[c['tp_key']]] = neg_by_tp.get(r[c['tp_key']], 0) + 1
    si = {r['run']: r for r in si_rows}
    with open(os.path.join(RAW, c['tag'] + '_genus_pct.tsv'), encoding='utf-8', newline='') as fh:
        rd = csv.reader(fh, delimiter='\t')
        hdr = next(rd)
        rows = [r for r in rd if r[0] in si and si[r[0]].get('body_site') == c['site']
                and si[r[0]].get('group') != '음성대조']
    samples = [{'run': r[0], 'subject': si[r[0]]['subject'], 'tp': si[r[0]][c['tp_key']]} for r in rows]
    genera = hdr[1:]
    X = {g: [float(r[j]) if r[j] not in ('', 'NA') else 0.0 for r in rows] for j, g in enumerate(genera, 1)}
    tps = tp_sorted({s['tp'] for s in samples})
    return samples, genera, X, tps, neg_by_tp


# ── FASTQ 머리줄 → 시퀀싱 런 ────────────────────────────────────────────────
ILLUMINA = re.compile(r'^@\S+\s+([A-Za-z0-9_-]+):(\d+):([A-Za-z0-9-]+):(\d+):')
RUN_FIELDS = ['cohort', 'run', 'subject', 'timepoint', 'status', 'n_headers_checked', 'first_header',
              'instrument', 'run_number', 'flowcell', 'run_id', 'fetched_on']


def fetch_headers(run, n=N_HEAD):
    url = ('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=' + run +
           '&result=read_run&fields=fastq_ftp&format=tsv')
    with urllib.request.urlopen(url, timeout=60) as r:
        lines = r.read().decode('utf-8').strip().split('\n')
    ftp = lines[1].split('\t')[-1].split(';')[0]
    req = urllib.request.Request('https://' + ftp, headers={'Range': 'bytes=0-16383'})
    with urllib.request.urlopen(req, timeout=60) as r:
        raw = r.read()
    text = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(raw).decode('utf-8', 'replace').split('\n')
    return [text[4 * i].strip() for i in range(n) if 4 * i < len(text) and text[4 * i].startswith('@')]


def parse_header(h):
    m = ILLUMINA.match(h or '')
    if not m:
        return None
    inst, runno, fc, _lane = m.groups()
    return inst, runno, fc, f'{inst}:{runno}:{fc}'


def run_table(all_samples):
    """런번호 TSV 를 읽고, 없거나 실패했거나 확인한 머리줄이 모자란 행만 다시 받는다."""
    cached = {}
    if os.path.exists(RUNS_TSV):
        for r in read_tsv(RUNS_TSV):
            if r['status'] in ('ok', 'ok_no_runfield') and r.get('n_headers_checked', '') == str(N_HEAD):
                cached[r['run']] = r
    todo = [(k, s) for k, ss in all_samples.items() for s in ss if s['run'] not in cached]

    def job(item):
        k, s = item
        row = {'cohort': k, 'run': s['run'], 'subject': s['subject'], 'timepoint': s['tp'],
               'fetched_on': datetime.date.today().isoformat()}
        try:
            hs = fetch_headers(s['run'])
            parsed = [parse_header(h) for h in hs]
            ids = {p[3] for p in parsed if p}
            row.update(n_headers_checked=str(len(hs)), first_header=hs[0][:120] if hs else '')
            if len(ids) == 1 and all(parsed):
                inst, runno, fc, rid = parsed[0]
                row.update(status='ok', instrument=inst, run_number=runno, flowcell=fc, run_id=rid)
            elif len(ids) > 1:                                   # 한 파일에 여러 런 — 따로 표시
                row.update(status='ok', instrument='', run_number='', flowcell='',
                           run_id='mixed:' + '+'.join(sorted(ids)))
            else:
                row.update(status='ok_no_runfield', instrument='', run_number='', flowcell='', run_id='')
        except Exception as e:                                   # 네트워크 실패 → 정보 없음
            row.update(status='fetch_failed: ' + type(e).__name__, n_headers_checked='0', first_header='',
                       instrument='', run_number='', flowcell='', run_id='')
        return row

    if todo:
        print(f'  FASTQ 머리줄 받기 — {len(todo)}개(ENA · 각 파일 앞 16 KB · 첫 {N_HEAD}읽기)')
        with ThreadPoolExecutor(max_workers=8) as ex:
            for row in ex.map(job, todo):
                cached[row['run']] = row
    order = [s['run'] for k in all_samples for s in all_samples[k]]
    with open(RUNS_TSV, 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=RUN_FIELDS, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for run in order:
            w.writerow({k: cached[run].get(k, '') for k in RUN_FIELDS})
    return {run: cached[run].get('run_id', '') for run in order}


# ── E0-② · ④ ────────────────────────────────────────────────────────────────
def run_crosstab(samples, runid, tps):
    tab = {}
    for s in samples:
        rid = runid.get(s['run']) or '(정보 없음)'
        tab.setdefault(rid, {t: 0 for t in tps})[s['tp']] += 1
    return tab


def run_overlap(samples, runid, t1, t2):
    """crossed = 두 시점이 쓴 모든 런에 두 시점 시료가 함께 있음. 공유 런에 속한 시료 비율도 낸다."""
    sel = [s for s in samples if s['tp'] in (t1, t2)]
    if any(not runid.get(s['run']) or runid[s['run']].startswith('mixed:') for s in sel):
        return {'status': 'unknown', 'frac_samples_in_shared_runs': None}
    r1 = {runid[s['run']] for s in sel if s['tp'] == t1}
    r2 = {runid[s['run']] for s in sel if s['tp'] == t2}
    shared = r1 & r2
    frac = sum(1 for s in sel if runid[s['run']] in shared) / len(sel)
    status = 'disjoint' if not shared else 'crossed' if r1 == r2 else 'partial'
    return {'status': status, 'frac_samples_in_shared_runs': round(frac, 3)}


def confined(samples, genera, X, tps, thr):
    out = []
    for g in genera:
        share = {}
        for t in tps:
            v = [x for x, s in zip(X[g], samples) if s['tp'] == t]
            share[t] = sum((x >= thr) if thr > 0 else (x > 0) for x in v) / len(v)
        for t in tps:
            if share[t] >= HIGH and all(share[o] <= LOW for o in tps if o != t):
                out.append({'genus': g, 'timepoint': t, 'share': {k: round(v, 3) for k, v in share.items()}})
    return out


def e0_state(overlap, flagged_primary, neg_by_tp, t1, t2):
    reasons = []
    hit = [f['genus'] for f in flagged_primary if f['timepoint'] in (t1, t2)]
    if overlap['status'] == 'disjoint':
        reasons.append('두 시점의 시퀀싱 런이 하나도 겹치지 않음')
    elif overlap['status'] == 'partial':
        reasons.append(f"두 시점이 시퀀싱 런을 일부만 공유함(공유 런 시료 {overlap['frac_samples_in_shared_runs']:.0%})")
    if hit:
        reasons.append(f'한 시점에만 몰린 속 {len(hit)}개({PRIMARY_THR:g} % 문턱)')
    if reasons:
        return '혼동 의심', reasons
    missing = []
    if overlap['status'] == 'unknown':
        missing.append('시퀀싱 런 정보 없음')
    no_neg = [t for t in (t1, t2) if neg_by_tp.get(t, 0) == 0]
    if no_neg:
        missing.append('음성대조 없음(' + ' · '.join(no_neg) + ')')
    if missing:
        return '정보 없음', missing
    return '통과', []


# ── 표현 ────────────────────────────────────────────────────────────────────
def logratio(X, genera, n, excl=frozenset()):
    excl = set(excl) - set(CAND)
    rest = [sum(X[g][i] for g in genera if g not in CAND and g not in excl) for i in range(n)]
    return {c: [math.log((X[c][i] + EPS) / (rest[i] + EPS)) for i in range(n)] for c in CAND}


def pairratio(X, n):
    return {f'{a}/{b}': [math.log((X[a][i] + EPS) / (X[b][i] + EPS)) for i in range(n)] for a, b in PAIRS}


def by_subject(samples, vals):
    d = {}
    for s, v in zip(samples, vals):
        d.setdefault(s['subject'], {})[s['tp']] = v
    return d


def paired(samples, vals, t1, t2):
    d = by_subject(samples, vals)
    subj = sorted(k for k, v in d.items() if t1 in v and t2 in v)
    return subj, [d[k][t2] - d[k][t1] for k in subj]


# ── E1: 부호순위 · Hodges–Lehmann · 정확 CI ─────────────────────────────────
_CDF = {}


def signrank_cdf(n):
    """귀무가설에서 W+ 의 누적분포(동점 · 0 없음 가정) — 부분합 개수를 동적계획으로 센다."""
    if n not in _CDF:
        M = n * (n + 1) // 2
        c = [0] * (M + 1)
        c[0] = 1
        for k in range(1, n + 1):
            for w in range(M, k - 1, -1):
                c[w] += c[w - k]
        tot, acc, cdf = 2 ** n, 0, []
        for x in c:
            acc += x
            cdf.append(acc / tot)
        _CDF[n] = cdf
    return _CDF[n]


def hl_ci(d, alpha=ALPHA):
    n = len(d)
    walsh = sorted((d[i] + d[j]) / 2 for i in range(n) for j in range(i, n))
    M = len(walsh)
    cdf = signrank_cdf(n)
    ks = [w for w in range(M + 1) if cdf[w] <= alpha / 2]
    if not ks:
        return st.median(walsh), None, None
    k = ks[-1]                        # P(W+ ≤ k) ≤ α/2 인 가장 큰 k
    return st.median(walsh), walsh[k], walsh[M - 1 - k]


def e1(d, alt):
    n = len(d)
    n_zero = sum(1 for x in d if x == 0)
    ties = len({abs(x) for x in d if x != 0}) < n - n_zero
    method = 'exact' if (n <= 50 and not ties and n_zero == 0) else 'approx'
    p_claim = float(stats.wilcoxon(d, alternative=alt, method=method).pvalue)
    p_two = p_claim if alt == 'two-sided' else float(stats.wilcoxon(d, alternative='two-sided',
                                                                    method=method).pvalue)
    hl, lo, hi = hl_ci(d)
    return {'n': n, 'median_delta': round(st.median(d), 4), 'hl': round(hl, 4),
            'hl_ci95': [round(lo, 4), round(hi, 4)] if lo is not None else None,
            'p': p_claim, 'alternative': alt, 'p_two_sided': p_two, 'p_method': method,
            'n_zero': n_zero, 'ties': ties, 'ci_assumes_no_zero_no_ties': True}


def holm(ps):
    m = len(ps)
    order = sorted(range(m), key=lambda i: ps[i])
    adj, run = [0.0] * m, 0.0
    for rank_, i in enumerate(order):
        run = max(run, min(1.0, (m - rank_) * ps[i]))
        adj[i] = run
    return adj


def e1_decision(rows, genera):
    """Holm 을 주장 방향 p 와 양측 p 에 각각 건다. 반환: {균: (유의?, 방향 맞음?)}"""
    hc = holm([rows[c]['p'] for c in genera])
    h2 = holm([rows[c]['p_two_sided'] for c in genera])
    out = {}
    for c, pc, p2 in zip(genera, hc, h2):
        rows[c]['p_holm'], rows[c]['p_two_sided_holm'] = pc, p2
        alt = rows[c]['alternative']
        if alt == 'two-sided':
            out[c] = (pc < ALPHA, True)
            continue
        same = (rows[c]['hl'] > 0) == (alt == 'greater')
        if pc < ALPHA and same:
            out[c] = (True, True)
        elif p2 < ALPHA and not same:                            # 결정표 8번 — 주장과 반대 방향
            out[c] = (True, False)
        else:
            out[c] = (False, True)
    return out


# ── E2: 개인 내 SD · CR ─────────────────────────────────────────────────────
def sw_oneway(groups):
    n, k = len(groups), len(groups[0])
    msw = sum((x - sum(g) / k) ** 2 for g in groups for x in g) / (n * (k - 1))
    return math.sqrt(msw)


def sw_twoway(groups):
    n, k = len(groups), len(groups[0])
    gm = sum(x for g in groups for x in g) / (n * k)
    rm = [sum(g) / k for g in groups]
    cm = [sum(g[j] for g in groups) / n for j in range(k)]
    sse = sum((groups[i][j] - rm[i] - cm[j] + gm) ** 2 for i in range(n) for j in range(k))
    return math.sqrt(sse / ((n - 1) * (k - 1)))


def beyond(d, cr, alt):
    n = len(d)
    up, down = sum(1 for x in d if x > cr), sum(1 for x in d if x < -cr)

    def ci(k):
        c = stats.binomtest(k, n).proportion_ci(confidence_level=0.95, method='exact')
        return [round(c.low, 3), round(c.high, 3)]

    out = {'CR': round(cr, 4), 'n': n, 'k_up': up, 'k_up_ci95': ci(up), 'k_down': down, 'k_down_ci95': ci(down)}
    # 참고: 잡음만이라면 한쪽으로 CR 을 넘을 확률 ≈ 2.5 %, 양쪽 합 ≈ 5 %(정규 · 같은 S_w 가정)
    if alt == 'two-sided':
        k, p0 = up + down, 0.05
    else:
        k, p0 = (up, 0.025) if alt == 'greater' else (down, 0.025)
    out.update(k_claim=k, k_claim_ci95=ci(k), binom_null=p0,
               binom_p=float(stats.binomtest(k, n, p0, alternative='greater').pvalue))
    return out


# ── 3장 결정표 ──────────────────────────────────────────────────────────────
def grade(has_control, e0, e1_sig=False, e1_dir_ok=True, e2_gate=None, e4=False,
          group_batch_overlap=False, e3_sig=False, e3_dir_ok=True):
    """문서 3장 결정표 1–8번. e2_gate: None = 같은 시험 잡음 없음(참고만) · True/False = 관문 결과."""
    if not has_control:
        if e0 == '혼동 의심':
            return '아직 안 됨', 1, '시점과 배치가 겹쳐 제품 효과와 처리 차이를 가를 수 없음'
        if e1_sig and e1_dir_ok and e2_gate is not False:
            return '조건부', 7, '전후 변화는 관찰됐지만 원인을 가를 대조가 없음'
        why = ('전후 변화가 유의하지 않음' if not e1_sig else
               '주장과 반대 방향' if not e1_dir_ok else '같은 시험 잡음 기준을 넘지 못함')
        return '아직 안 됨', 8, why
    if group_batch_overlap:
        return '아직 안 됨', 2, '군 차이와 처리 차이를 가를 수 없음'
    if e3_sig and not e3_dir_ok:
        return '아직 안 됨', 3, '주장과 반대 방향의 변화'
    if not e3_sig:
        return '아직 안 됨', 4, '대조와 비교해 유의한 차이를 보이지 않음'
    missing = []
    if e0 != '통과':
        missing.append(f'E0 {e0}')
    if not e4:
        missing.append('기준 시점 반복측정 없음')
    if e2_gate is not True:
        missing.append('같은 시험 잡음 기준 없음' if e2_gate is None else '같은 시험 잡음 기준을 넘지 못함')
    if not missing:
        return '근거로 쓸 수 있음', 5, ''
    return '조건부', 6, ' · '.join(missing)


def _selftest():
    """결정표 · 정확 CI · 런 겹침 · 반대 방향 판정 검사. 대조가 있는 줄(2–6)은 여기서만 확인한다."""
    cases = [
        (dict(has_control=False, e0='혼동 의심', e1_sig=True), 1),
        (dict(has_control=False, e0='정보 없음', e1_sig=True), 7),
        (dict(has_control=False, e0='통과', e1_sig=True, e2_gate=False), 8),
        (dict(has_control=False, e0='통과', e1_sig=False), 8),
        (dict(has_control=False, e0='통과', e1_sig=True, e1_dir_ok=False), 8),
        (dict(has_control=True, e0='통과', group_batch_overlap=True, e3_sig=True), 2),
        (dict(has_control=True, e0='통과', e3_sig=True, e3_dir_ok=False), 3),
        (dict(has_control=True, e0='통과', e3_sig=False), 4),
        (dict(has_control=True, e0='통과', e3_sig=True, e4=True, e2_gate=True), 5),
        (dict(has_control=True, e0='정보 없음', e3_sig=True, e4=True, e2_gate=True), 6),
        (dict(has_control=True, e0='통과', e3_sig=True, e4=False, e2_gate=None), 6),
    ]
    for kw, row in cases:
        got = grade(**kw)[1]
        assert got == row, (kw, got, row)
    # n = 25 · 양측 0.05 의 부호순위 임계값은 89(표준 표값) — P(W+ ≤ 89) ≤ 0.025 < P(W+ ≤ 90)
    cdf = signrank_cdf(25)
    assert cdf[89] <= 0.025 < cdf[90], (cdf[89], cdf[90])
    # 정확 분포의 p 가 scipy(정확)와 같은지 — 동점 없는 고정 자료
    d = [0.31, -0.12, 0.55, 0.08, -0.41, 0.77, 0.19, 0.63, -0.05, 0.92, 0.27, 0.48]
    wplus = sum(r for r, x in zip(stats.rankdata([abs(x) for x in d]), d) if x > 0)
    p_mine = 1 - signrank_cdf(len(d))[int(round(wplus)) - 1]
    p_scipy = stats.wilcoxon(d, alternative='greater', method='exact').pvalue
    assert abs(p_mine - p_scipy) < 1e-12, (p_mine, p_scipy)
    # 런 겹침: 0W 24개는 A, 1개는 4W 와 같은 B → partial(통과 아님)
    smp = [{'run': f'a{i}', 'tp': '0W'} for i in range(25)] + [{'run': f'b{i}', 'tp': '4W'} for i in range(25)]
    rid = {f'a{i}': ('B' if i == 0 else 'A') for i in range(25)} | {f'b{i}': 'B' for i in range(25)}
    assert run_overlap(smp, rid, '0W', '4W')['status'] == 'partial'
    rid2 = {f'a{i}': ('A' if i % 2 else 'B') for i in range(25)} | {f'b{i}': ('A' if i % 2 else 'B') for i in range(25)}
    assert run_overlap(smp, rid2, '0W', '4W')['status'] == 'crossed'
    # 반대 방향: 증가를 주장했는데 모두 줄었으면 '주장과 반대 방향'(8번)
    dd = [-0.5 - 0.01 * i for i in range(25)]
    rows = {c: e1(dd, CLAIM[c]) for c in CAND}
    dec = e1_decision(rows, CAND)
    assert dec['Cutibacterium'] == (True, False), dec
    assert grade(False, '통과', *dec['Cutibacterium'])[2] == '주장과 반대 방향'
    return {'decision_table_cases': len(cases), 'signrank_crit_n25': 89, 'exact_p_matches_scipy': True,
            'run_overlap_partial_and_crossed': True, 'opposite_direction_reaches_row8': True, 'passed': True}


# ── 그림 ────────────────────────────────────────────────────────────────────
def figure(samples, X, flagged_genera, runid, tps):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    runs = sorted({runid[s['run']] or 'unknown' for s in samples})
    pal = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']        # dataviz 기본 범주 1–4(검사 통과)
    mk = ['o', 's', '^', 'D']                                   # 색 + 모양 이중 표시(대비 경고 보완)
    assert len(runs) <= len(pal), f'런이 {len(runs)}개 — 색을 돌려 쓰지 않는다(dataviz 규칙)'
    style = {r: (pal[i], mk[i]) for i, r in enumerate(runs)}
    fsum = [sum(X[g][i] for g in flagged_genera) for i in range(len(samples))]
    cuti = X['Cutibacterium']
    renorm = [c / (100 - f) * 100 if f < 100 else float('nan') for c, f in zip(cuti, fsum)]
    xpos = {t: i for i, t in enumerate(tps)}
    # 점 흩뿌림 — hash() 는 실행마다 달라 재현이 안 되므로 md5 로 고정
    jit = [((int(hashlib.md5(s['subject'].encode()).hexdigest(), 16) % 97) / 97 - 0.5) * 0.36 for s in samples]

    plt.rcParams.update({'font.size': 13, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.9))
    leg = dict(frameon=False, handletextpad=0.4, handlelength=1.0, borderaxespad=0.3)
    a = ax[0]
    for r in runs:
        idx = [i for i, s in enumerate(samples) if (runid[s['run']] or 'unknown') == r]
        col, m = style[r]
        lab = (f'run {r.split(":")[1]}' if ':' in r else 'run unknown') + f' (n = {len(idx)})'
        a.scatter([xpos[samples[i]['tp']] + jit[i] for i in idx], [fsum[i] for i in idx], s=70, c=col,
                  marker=m, edgecolors='black', linewidths=0.9, label=lab, zorder=3)
    for t in tps:
        v = [fsum[i] for i, s in enumerate(samples) if s['tp'] == t]
        a.hlines(st.median(v), xpos[t] - 0.28, xpos[t] + 0.28, colors='black', linewidth=2.2, zorder=4)
    a.set_xticks(range(len(tps)), tps)
    a.set_ylabel(f'Sum of the {len(flagged_genera)} genera (% of reads)')
    a.set_title(f'A  {len(flagged_genera)} genera detected almost only at baseline', loc='left',
                fontweight='bold')
    a.legend(title='MiSeq run', loc='upper right', **leg)
    a.grid(axis='y', color='#e4e3df', linewidth=0.8, zorder=0)

    b = ax[1]
    ymax = 0.0
    for vals, col, m, lab, off in ((cuti, '#8a8983', 'o', 'all genera (original %)', -0.15),
                                   (renorm, '#6250d6', 's', f'{len(flagged_genera)} baseline-only genera removed',
                                    0.15)):
        med = [st.median([v for v, s in zip(vals, samples) if s['tp'] == t]) for t in tps]
        lab = f'{lab}: median ' + ' → '.join(f'{x:.1f}' for x in med)     # 끝 숫자는 범례에 — 어느 선인지 분명하게
        b.scatter([xpos[s['tp']] + off + jit[i] * 0.4 for i, s in enumerate(samples)], vals, s=55, c=col,
                  marker=m, edgecolors='black', linewidths=0.8, label=lab, zorder=3)
        b.plot([xpos[t] + off for t in tps], med, color=col, linewidth=2.2, zorder=4)
        ymax = max(ymax, max(v for v in vals if v == v))
    b.set_ylim(0, ymax * 1.32)                                   # 범례 자리 — 점과 겹치지 않게
    b.set_xticks(range(len(tps)), tps)
    b.set_ylabel(r'$\it{Cutibacterium}$ (% of reads)')
    b.set_title(r'B  $\bf{\it{Cutibacterium}}$, line = median', loc='left', fontweight='bold')
    b.legend(loc='upper left', fontsize=11.5, **leg)
    b.grid(axis='y', color='#e4e3df', linewidth=0.8, zorder=0)
    fig.text(0.01, 0.005, 'Public cosmetic trial data (cheek, n = 25, no untreated control). '
             'Sequencing run read from FASTQ headers. Panel B is on the % scale; '
             'the analysis uses log-ratios.', fontsize=10, color='#6b6a64')
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    for ext in ('png', 'pdf'):
        fig.savefig(f'{OUT_FIG}.{ext}', dpi=300)
    plt.close(fig)
    with open(OUT_FIGDATA, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['run', 'subject', 'timepoint', 'sequencing_run', 'baseline_only_genera_sum_pct',
                    'cutibacterium_pct_original', 'cutibacterium_pct_baseline_only_removed'])
        for i, s in enumerate(samples):
            w.writerow([s['run'], s['subject'], s['tp'], runid[s['run']], round(fsum[i], 4),
                        round(cuti[i], 4), round(renorm[i], 4)])


# ── 실행 ────────────────────────────────────────────────────────────────────
def main():
    out = {'question': '화장품 시험의 균 변화를 근거로 말할 수 있나 — 효능 판정 규칙 v2 시범',
           'rules': '경진대회/효능판정기준_초안_261017.md v2 (결정 1–8 · 3장 결정표 · 5-2)',
           'note': 'E1 · E2 · E0-③ 은 참고. 07 의 등급은 E0 로 정해진다(결정 5 · 7).',
           'params': {'EPS': EPS, 'thresholds_pct': THRESH, 'primary_threshold_pct': PRIMARY_THR,
                      'confined_rule': f'한 시점 ≥ {HIGH:.0%} · 다른 모든 시점 ≤ {LOW:.0%}',
                      'run_overlap_rule': 'crossed(모든 런에 두 시점이 함께) 일 때만 통과 — partial · disjoint 는 겹침',
                      'CR': '1.96·√2·S_w', 'alpha': ALPHA, 'claim_direction': CLAIM,
                      'E2_gate': 'None — 같은 시험 잡음의 숫자 기준은 아직 정하지 않음(문서 빈칸)'},
           'inputs': input_fingerprints()}
    out['selftest'] = _selftest()
    print('자가검사 통과 — 결정표 11개 · n=25 임계값 89 · 정확 p = scipy · 런 partial/crossed · 반대 방향 → 8번')

    D = {k: load(k) for k in COHORTS}
    runid = run_table({k: D[k][0] for k in D})

    # ── 15계절: 잡음 기준(S_w · CR)
    s15, g15, X15, t15, neg15 = D['15']
    lr15 = logratio(X15, g15, len(s15))
    noise = {}
    for c in CAND:
        d = by_subject(s15, lr15[c])
        groups = [[d[k][t] for t in t15] for k in sorted(d) if all(t in d[k] for t in t15)]
        sw1, sw2 = sw_oneway(groups), sw_twoway(groups)
        noise[c] = {'n_subject': len(groups), 'n_repeat': len(t15), 'S_w_oneway': round(sw1, 4),
                    'S_w_twoway': round(sw2, 4), 'CR_oneway': round(CR_K * sw1, 4),
                    'CR_twoway': round(CR_K * sw2, 4)}
    out['noise_15season'] = {'note': '홍콩 중국인 11가구 24명 · 이마 · V4 · 2014-01·04·08·11 · 로그비 원자료',
                             'genera': noise}

    # ── 코호트별 E0
    e0info = {}
    for k, (ss, gg, XX, tt, nn) in D.items():
        grid, sets = {}, {}
        for thr in THRESH:
            fl = confined(ss, gg, XX, tt, thr)
            sets[thr] = {f['genus'] for f in fl}
            sums = {}
            for t in tt:
                idx = [i for i, s in enumerate(ss) if s['tp'] == t]
                sums[t] = round(st.median([sum(XX[f['genus']][i] for f in fl) for i in idx]), 3) if fl else 0.0
            grid[f'{thr:g}%'] = {'n_genera': len(fl), 'genera': fl, 'sum_median_pct_by_timepoint': sums}
        hi, mid, lo = (sets[t] for t in THRESH)
        e0info[k] = {'timepoints': tt, 'n_samples': len(ss), 'negative_controls_by_timepoint': nn,
                     'run_crosstab': run_crosstab(ss, runid, tt), 'confined_grid': grid,
                     'confined_nesting': {'nested': hi <= mid <= lo,
                                          'intersection_1_and_0.1': sorted(hi & mid),
                                          'intersection_1_and_0': sorted(hi & lo),
                                          'union_n': len(hi | mid | lo)},
                     'E4_baseline_repeat': False,
                     'E4_note': '바르기 전(기준 시점) 반복측정이 설계에 없다' if k == '07'
                                else '관찰 연구 — 기준 시점 반복측정 개념이 없다'}
    out['E0'] = e0info

    # ── 07: 비교 · 판 · 균
    s07, g07, X07, t07, neg07 = D['07']
    n07 = len(s07)
    fl_primary = e0info['07']['confined_grid'][f'{PRIMARY_THR:g}%']['genera']
    versions = {'raw': logratio(X07, g07, n07)}
    union = set()
    for thr in THRESH:
        excl = {f['genus'] for f in e0info['07']['confined_grid'][f'{thr:g}%']['genera']}
        union |= excl
        versions[f'excl_{thr:g}pct'] = logratio(X07, g07, n07, excl)
    versions['excl_union'] = logratio(X07, g07, n07, union)
    pair = pairratio(X07, n07)
    comps = [('0W', '2W'), PRIMARY, ('2W', '4W')]

    cr_int = {}
    for v, lr in versions.items():
        cr_int[v] = {}
        for c in CAND:
            d = by_subject(s07, lr[c])
            groups = [[d[s]['2W'], d[s]['4W']] for s in sorted(d) if '2W' in d[s] and '4W' in d[s]]
            cr_int[v][c] = {'oneway': CR_K * sw_oneway(groups), 'twoway': CR_K * sw_twoway(groups)}

    res07, delta_cols, subjects = {}, {}, set()
    for t1, t2 in comps:
        cn = cname(t1, t2)
        ov = run_overlap(s07, runid, t1, t2)
        state, why = e0_state(ov, fl_primary, neg07, t1, t2)
        block = {'run_overlap': ov, 'E0_state': state, 'E0_reasons': why, 'versions': {}}
        for v, lr in versions.items():
            rows = {}
            for c in CAND:
                subj, d = paired(s07, lr[c], t1, t2)
                subjects |= set(subj)
                delta_cols[f'{v}|{cn}|{c}'] = dict(zip(subj, d))
                r = e1(d, CLAIM[c])
                r['E2_borrowed_oneway'] = beyond(d, noise[c]['CR_oneway'], CLAIM[c])
                r['E2_borrowed_twoway'] = beyond(d, noise[c]['CR_twoway'], CLAIM[c])
                if (t1, t2) != ('2W', '4W'):                 # 대리 잡음을 만든 비교 자체에는 쓰지 않는다
                    r['E2_internal_2W4W_oneway'] = beyond(d, cr_int[v][c]['oneway'], CLAIM[c])
                    r['E2_internal_2W4W_twoway'] = beyond(d, cr_int[v][c]['twoway'], CLAIM[c])
                rows[c] = r
            for c, (sig, dir_ok) in e1_decision(rows, CAND).items():
                rows[c]['grade'] = dict(zip(('grade', 'row', 'reason'),
                                            grade(False, state, e1_sig=sig, e1_dir_ok=dir_ok, e2_gate=None)))
            block['versions'][v] = rows
        prow = {}
        for name, vals in pair.items():
            subj, d = paired(s07, vals, t1, t2)
            delta_cols[f'pair|{cn}|{name}'] = dict(zip(subj, d))
            prow[name] = e1(d, 'two-sided')
        for name, pa in zip(pair, holm([prow[x]['p'] for x in pair])):
            prow[name]['p_holm'] = pa
        block['candidate_pair_logratio'] = prow
        res07[cn] = block
    pc = cname(*PRIMARY)
    out['cohort07'] = {'design': '볼 · 25명 × 0W · 2W · 4W · 대조 없음 · 음성대조 없음',
                       'primary_comparison': f'{pc} (부분 사후 — 시점별 중앙값을 본 뒤 정함)',
                       'excl_union_genera': sorted(union),
                       'CR_internal_2W4W': {v: {c: {m: round(x, 4) for m, x in dd.items()} for c, dd in d.items()}
                                            for v, d in cr_int.items()},
                       'comparisons': res07,
                       'grade_primary': {c: res07[pc]['versions']['raw'][c]['grade'] for c in CAND}}

    # ── 15계절: 결정표 규칙 시험(계절 쌍 6개)
    rule = {'pairs': {}, 'grade_counts': {}}
    fl15 = e0info['15']['confined_grid'][f'{PRIMARY_THR:g}%']['genera']
    for i in range(len(t15)):
        for j in range(i + 1, len(t15)):
            t1, t2 = t15[i], t15[j]
            ov = run_overlap(s15, runid, t1, t2)
            state, why = e0_state(ov, fl15, neg15, t1, t2)
            rows = {}
            for c in CAND:
                subj, d = paired(s15, lr15[c], t1, t2)
                r = e1(d, 'two-sided')
                r['E2_in_sample_oneway'] = beyond(d, noise[c]['CR_oneway'], 'two-sided')
                rows[c] = r
            for c, (sig, dir_ok) in e1_decision(rows, CAND).items():
                g = grade(False, state, e1_sig=sig, e1_dir_ok=dir_ok, e2_gate=None)
                rows[c]['grade'] = dict(zip(('grade', 'row', 'reason'), g))
                rule['grade_counts'][g[0]] = rule['grade_counts'].get(g[0], 0) + 1
            rule['pairs'][cname(t1, t2)] = {'run_overlap': ov, 'E0_state': state, 'E0_reasons': why,
                                            'genera': rows}
    n_all = sum(rule['grade_counts'].values())
    rule['structurally_guaranteed'] = ("대조가 없으므로 결정표 1 · 7 · 8번만 나올 수 있다 — "
                                       "'근거로 쓸 수 있음' 0건은 데이터와 무관하게 보장된다")
    rule['conditional_rate'] = f"{rule['grade_counts'].get('조건부', 0)} / {n_all}"
    rule['reading'] = "무처리 자료에서 나온 '조건부'는 '변화가 관찰됨 · 원인은 모름'이다(효능의 약한 근거가 아님)."
    out['rule_test_15season'] = rule

    # ── 저장: JSON 을 먼저(그림에서 실패해도 결과는 남게) → 사람별 변화 → 그림
    out['mem_time'] = mem_time_line(T0)
    save_json(out, OUT_JSON)
    cols = list(delta_cols)
    with open(OUT_DELTA, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['subject'] + cols)
        for s in sorted(subjects):
            w.writerow([s] + [round(delta_cols[c][s], 6) if s in delta_cols[c] else '' for c in cols])
    figure(s07, X07, [f['genus'] for f in fl_primary], runid, t07)

    # ── 화면 요약
    print('\n■ 07 화장품 — E0')
    for r, c in e0info['07']['run_crosstab'].items():
        print(f'  런 {r}: {c}')
    for thr in THRESH:
        gq = e0info['07']['confined_grid'][f'{thr:g}%']
        print(f'  한 시점에만 몰린 속({thr:g} % 문턱): {gq["n_genera"]}개 · 합계 중앙값 {gq["sum_median_pct_by_timepoint"]}')
    print(f"  문턱별 목록 포함 관계: {e0info['07']['confined_nesting']}")
    b = res07[pc]
    print(f'  {pc}: 런 {b["run_overlap"]} · E0 {b["E0_state"]} — {" · ".join(b["E0_reasons"])}')
    print(f'\n■ 07 등급(주 비교 {pc})')
    for c in CAND:
        g = out['cohort07']['grade_primary'][c]
        print(f'  {c:16s} {g["grade"]}(결정표 {g["row"]}번) — {g["reason"]}')
    print(f'\n■ 07 참고 수치({pc}) — 판마다 HL 추정치 [95% CI] · Holm p(주장 방향) · 빌린 CR(1원) 넘은 사람')
    for v in versions:
        for c in CAND:
            r = b['versions'][v][c]
            e2 = r['E2_borrowed_oneway']
            print(f'  {v:12s} {c:16s} HL {r["hl"]:+.3f} {r["hl_ci95"]} · p_holm {r["p_holm"]:.2g} · '
                  f'k {e2["k_claim"]}/{e2["n"]}')
    print('\n■ 15계절 규칙 시험 — 계절 쌍 6개 × 3속')
    print(f'  등급 집계: {rule["grade_counts"]} · 조건부 {rule["conditional_rate"]}')
    print(f'  E0: {sorted({p["E0_state"] for p in rule["pairs"].values()})}')
    print('\n' + out['mem_time'])
    print(f'저장: {OUT_JSON}\n      {OUT_DELTA}\n      {RUNS_TSV}\n      {OUT_FIG}.png/.pdf · {OUT_FIGDATA}')


if __name__ == '__main__':
    main()
