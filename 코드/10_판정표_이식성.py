"""판정표가 우리 네 편 밖에서도 작동하는지 확인한다 — 중간발표 피드백 ⑤ 정면 답변.

피드백 원문: *"4편에만 해당하는. 마이크로바이옴이 아니라면 활용불가하겠네? …
4개의 논문 이외에 아예 쓸 수 없는건가 그럼 방향이 맞을까?"*

⚠️ **이것은 재현성 검정이 아니다.** 어떤 균의 질환 연관이 다른 자료에서도 반복되는지는
   묻지 않는다. 그건 다른 질문이고 우리 자료로 답할 수 없다. 여기서 묻는 것은 하나다 —

       **5관문 중 무엇이 새 자료에서 라벨 없이 판정되는가?**

   ①검출 안정성·②교란 취약성은 민감/비민감 라벨이 필요 없다. 시료와 심도만 있으면 된다.
   ③분류 신뢰성은 종 수준 분류표가, ④판별력은 군 라벨이, ⑤방향성은 종 수준 계수가 있어야 한다.
   그 경계를 자료로 보이는 것이 이 스크립트의 산출물이고, 그대로 **「자료 요건표」**가 된다.

쓰는 자료 — 질환 코호트는 쓰지 않는다
   16S: 595명(얼굴 세척액) · 42명 03민감(볼) · 02건강(이마/볼) · 07화장품(볼) ·
        15계절(이마) · 16비듬(두피)
   ITS: 42명 03민감(볼) · 58명 01노화(볼·이마) · 15계절(이마) · 16비듬(두피)

결과: 분석표/판정표_이식성.json
"""
import csv, os, statistics as st
from stats_util import save_json, TAB, RAW, CAND, f, spearman, load_595, load_42

# ── 관문 문턱 ───────────────────────────────────────────────────────────────
# 15차 판정표에 쓴 판정을 문턱으로 명문화한 것이다. 사후 규칙이며 사전등록이 아니다.
#   ① 검출 안정성 — 595 Cutibacterium 99.16% 통과 / Corynebacterium 90.42% 주의.
#      심도 상관도 ① 쪽이다. 15차 표가 Corynebacterium 을 「주의 90.42% · 저심도 81.4% ·
#      42명 심도 rho=-0.385」로 묶어 적은 그대로다. **음의 상관만** 문제 삼는다 —
#      얕게 읽을수록 비율이 커지는 것이 시약 배경의 특징이기 때문이다(6_오염점검.py 와 같은 방향).
#      양의 상관은 "깊게 읽을수록 더 잡힌다"는 정상 검출이므로 강등하지 않는다.
#   ② 교란 취약성 — 나이다. 595 Cutibacterium 이 세 표현 모두 -0.28~-0.38 일 때 「취약」이었다.
#      교란은 방향과 무관하므로 절댓값으로 본다.
DET_PASS, DET_WARN = 95.0, 90.0
RHO_WEAK, RHO_WARN = 0.30, 0.15
NOT_GENUS = {'unidentified', 'Unassigned', 'metagenome'}

MAL = 'g__Malassezia'          # 수업 ITS 표의 Malassezia 열 이름
FUNGI_CAND = [MAL]


def gate1(det, rho_depth):
    """① 검출 안정성. 검출률이 1차 기준이고, 심도와의 음의 상관이 있으면 강등한다."""
    if det < DET_WARN:
        return '불가'
    if rho_depth <= -RHO_WEAK:
        return '주의 — 심도 의존'
    if det < DET_PASS:
        return '주의'
    return '통과'


def gate2(rho):
    """② 교란 취약성. 나이 상관의 절댓값."""
    a = abs(rho)
    return '취약' if a >= RHO_WEAK else ('주의' if a >= RHO_WARN else '통과')


# ── 코호트 적재 ─────────────────────────────────────────────────────────────
def load_class(tag, keep=None):
    """수업 공개 코호트 하나를 읽는다. keep(시료정보 dict) -> bool 로 부분집합을 고른다."""
    si = {r['run']: r for r in csv.DictReader(
        open(os.path.join(RAW, tag + '_sample_info.tsv'), encoding='utf-8'), delimiter='\t')}
    rd = csv.reader(open(os.path.join(RAW, tag + '_genus_pct.tsv'), encoding='utf-8'), delimiter='\t')
    hdr = next(rd)
    col = {g: i for i, g in enumerate(hdr)}
    rows = [r for r in rd if r[0] in si and (keep is None or keep(si[r[0]]))]
    info = [si[r[0]] for r in rows]
    val = {g: [f(r[col[g]]) or 0.0 for r in rows] for g in hdr[1:]}
    depth = [f(r.get('read_depth_final')) or 0.0 for r in info]
    age = [f(r.get('age')) for r in info]
    return {'info': info, 'val': val, 'depth': depth,
            'age': age if any(a is not None for a in age) else None, 'n': len(rows)}


def load_595_as_cohort():
    ids, genera, val, md, qc = load_595()
    return {'info': None, 'n': len(ids),
            'val': {g: [100 * x for x in val[g]] for g in genera},   # 비율 -> %
            'depth': [f(qc[i]['non_chimeric_reads']) or 0.0 for i in ids],
            'age': [f(md[i]['Age']) for i in ids]}


def load_42_as_cohort():
    info, val = load_42()
    return {'info': info, 'n': len(info), 'val': val,
            'depth': [f(r['read_depth_final']) or 0.0 for r in info],
            'age': [f(r['age']) for r in info]}


site = lambda s: (lambda r: r.get('body_site') == s and r.get('group') != '음성대조')

# 나이 상관을 그대로 읽으면 안 되는 코호트. 설계가 양 극단만 뽑았으면 |rho| 가 부풀려진다.
AGE_DESIGN = {
    '02건강 · 이마': '설계상 Old(49-67) 52 · Young(21-36) 50 두 군. 중간 연령이 없어 |rho| 가 커지는 쪽으로 치우친다',
    '02건강 · 볼': '설계상 Old(49-67) 52 · Young(21-36) 50 두 군. 중간 연령이 없어 |rho| 가 커지는 쪽으로 치우친다',
    '58명 01노화 · 볼·이마': '설계상 19-29세 26명 · 60-63세 32명. 30-50대가 0명이다. '
                            '연속 상관으로 읽으면 안 되고 두 군 비교로 봐야 한다. '
                            '15차 §6.2 가 이미 고연령 65.96% 대 저연령 84.71% 로 보고한 것과 같은 현상이다',
}

COHORTS_16S = [
    ('595명 · 얼굴 세척액', 'V3-V4', '여 595', '우리 주 자료', load_595_as_cohort),
    ('42명 · 오른쪽 볼', 'V4-V5', '여 42', '우리 판별 자료', load_42_as_cohort),
    ('02건강 · 이마', 'V3-V4', '한국 남녀', '새 자료', lambda: load_class('02_healthy_16S', site('Forehead'))),
    ('02건강 · 볼', 'V3-V4', '한국 남녀', '새 자료', lambda: load_class('02_healthy_16S', site('Cheek'))),
    ('07화장품 · 볼', 'V3-V4', '여', '새 자료', lambda: load_class('07_skincare_16S', site('Cheek'))),
    ('15계절 · 이마', 'V4', '미기재', '새 자료', lambda: load_class('15_season_16S', site('Forehead'))),
    ('16비듬 · 두피', 'V3', '여', '새 자료 · 다른 부위', lambda: load_class('16_dandruff_16S', site('Scalp'))),
]

COHORTS_ITS = [
    ('42명 · 오른쪽 볼', 'ITS1', '여 42', '우리 판별 자료',
     lambda: load_class('03_sensitive_ITS')),
    ('58명 01노화 · 볼·이마', 'ITS1', '여 58', '우리 교차 자료',
     lambda: load_class('01_aging_ITS')),
    ('15계절 · 이마', 'ITS1', '미기재', '새 자료',
     lambda: load_class('15_season_ITS', site('Forehead'))),
    ('16비듬 · 두피', 'ITS1', '여', '새 자료 · 다른 부위',
     lambda: load_class('16_dandruff_ITS', site('Scalp'))),
]


# ── 관문 ①② 판정 ───────────────────────────────────────────────────────────
def judge(coh, genus):
    v = coh['val'].get(genus)
    if v is None:
        return None
    n = len(v)
    det = 100 * sum(1 for x in v if x > 0) / n
    d = coh['depth']
    k = max(1, n // 10)
    lo_idx = sorted(range(n), key=lambda i: d[i])[:k]          # 심도 하위 10%
    hi_idx = sorted(range(n), key=lambda i: d[i])[k:]
    r_dep, p_dep = spearman(d, v)
    row = {'genus': genus, 'n': n,
           'detect_pct': round(det, 2), 'median_pct': round(st.median(v), 3),
           'detect_lowdepth_pct': round(100 * sum(1 for i in lo_idx if v[i] > 0) / len(lo_idx), 1),
           'detect_rest_pct': round(100 * sum(1 for i in hi_idx if v[i] > 0) / max(1, len(hi_idx)), 1),
           'rho_depth': round(r_dep, 3), 'p_depth': round(p_dep, 4),
           'gate1': gate1(det, r_dep)}
    if coh['age'] and sum(1 for a in coh['age'] if a is not None) >= 10:
        idx = [i for i in range(n) if coh['age'][i] is not None]
        r_age, p_age = spearman([coh['age'][i] for i in idx], [v[i] for i in idx])
        row.update({'rho_age': round(r_age, 3), 'p_age': round(p_age, 4),
                    'n_age': len(idx), 'gate2_age': gate2(r_age)})
    else:
        row.update({'rho_age': None, 'p_age': None, 'n_age': 0, 'gate2_age': '판정 불가 — 나이 없음'})
    return row


out = {'question': '5관문 중 무엇이 새 자료에서 라벨 없이 판정되는가',
       'not_a_replication_test':
           '질환 연관이 독립 자료에서 반복되는지를 묻는 검정이 아니다. 도구의 적용 가능성만 본다. '
           '질환 코호트(04·08·09·10·11·12·13·14·17)는 쓰지 않았다.',
       'thresholds': {'gate1_pass': DET_PASS, 'gate1_warn': DET_WARN,
                      'gate2_weak': RHO_WEAK, 'gate2_warn': RHO_WARN,
                      'note': '15차 판정표의 판정을 문턱으로 명문화한 사후 규칙'},
       'bacteria': {}, 'fungi': {}}

print('=' * 100)
print('판정표 이식성 — 5관문 중 무엇이 새 자료에서 라벨 없이 판정되는가')
print('=' * 100)

for kind, cohorts, cands, store in (('세균(16S)', COHORTS_16S, CAND, out['bacteria']),
                                    ('진균(ITS)', COHORTS_ITS, FUNGI_CAND, out['fungi'])):
    print(f'\n■ {kind}')
    print(f"  {'코호트':22s} {'영역':7s} {'n':>5s}  {'후보':16s} {'검출률':>8s} {'중앙값':>9s} "
          f"{'저심도검출':>9s} {'rho(심도)':>10s} {'① 검출 안정성':>14s} {'rho(나이)':>10s} {'② 교란':>8s}")
    for name, region, sex, role, loader in cohorts:
        try:
            coh = loader()
        except Exception as e:
            print(f"  {name:22s} — 읽지 못함 ({e})")
            continue
        entry = {'region': region, 'subjects': sex, 'role': role, 'n': coh['n'],
                 'depth_median': round(st.median(coh['depth']), 0) if coh['depth'] else None,
                 'has_age': bool(coh['age']), 'genera': {}}
        if name in AGE_DESIGN:
            entry['age_design_caveat'] = AGE_DESIGN[name]
        # 속 수준 미분류 비율 — ③의 대리지표일 뿐 ③ 판정이 아니다
        unk = [g for g in coh['val'] if g in NOT_GENUS or g.startswith('g__unidentified')]
        if unk:
            tot = [sum(coh['val'][g][i] for g in unk) for i in range(coh['n'])]
            entry['unassigned_median_pct'] = round(st.median(tot), 2)
        for g in cands:
            row = judge(coh, g)
            if row is None:
                print(f"  {name:22s} {region:7s} {coh['n']:5d}  {g[:16]:16s}  — 이 자료에 없는 속")
                continue
            entry['genera'][g] = row
            ra = f"{row['rho_age']:+10.3f}" if row['rho_age'] is not None else f"{'-':>10s}"
            g2 = row['gate2_age'] + ('*' if name in AGE_DESIGN else '')
            print(f"  {name:22s} {region:7s} {coh['n']:5d}  {g[:16]:16s} {row['detect_pct']:7.2f}% "
                  f"{row['median_pct']:8.3f}% {row['detect_lowdepth_pct']:8.1f}% "
                  f"{row['rho_depth']:+10.3f} {row['gate1']:>14s} {ra} {g2:>8s}")
        store[name] = entry
    if any(n in AGE_DESIGN for n, *_ in cohorts):
        print('  * 나이 상관을 그대로 읽으면 안 되는 코호트 — 설계가 양 극단만 뽑았다. JSON 의 age_design_caveat 참조')


# ── 관문별 적용 가능성 요약 = 자료 요건표 ──────────────────────────────────
REQ = [
    ('① 검출 안정성', '시료별 read 심도',             '모든 코호트',            True),
    ('② 교란 취약성', '심도 + 교란변수(나이·피지)',   '심도는 모두 · 나이는 일부', True),
    ('③ 분류 신뢰성', '종 수준 분류표 또는 ASV',      '속 수준 표만 있으면 불가', False),
    ('④ 판별력',      '군 라벨(민감/비민감)',          '라벨 있는 코호트만',      False),
    ('⑤ 방향성',      '종 수준 상대량 + 표현형 연계',  '속 수준 표만 있으면 불가', False),
]

n16 = len(out['bacteria'])
nits = len(out['fungi'])
n_age = sum(1 for e in out['bacteria'].values() if e['has_age'])

out['requirement_table'] = [
    {'gate': g, 'needs': need, 'availability': avail, 'label_free': lf} for g, need, avail, lf in REQ]
out['summary'] = {
    'n_cohort_16S': n16, 'n_cohort_ITS': nits,
    'gate1_judged': f'{n16 + nits}/{n16 + nits} 코호트',
    'gate2_depth_judged': f'{n16 + nits}/{n16 + nits} 코호트',
    'gate2_age_judged': f'{n_age}/{n16} 세균 코호트 (나이가 기재된 자료만)',
    # 2026-09-19 — 이 두 문자열에 «속 요약표 기준» 이라는 단서를 넣었다.
    # 같은 수업 자료의 asv/ 표를 열어 보니 ⑤ 는 진균에서 «부분 가능» 이었다(11_진균_종수준.py).
    # 분석표/판정표_이식성.json 의 같은 필드도 똑같이 고쳐 두었다.
    'gate3_judged': '0 (속 요약표 11개 기준) — 수업 공개 자료의 속 수준 표로는 종 수준 분류 신뢰성을 볼 수 없다. 단 같은 자료의 asv/ 표를 쓰면 달라질 수 있다 — 2026-09-19 현재 ③은 ASV로 다시 재지 않았다',
    'gate4_judged': '2 — 03민감(민감/비민감) 뿐. 나머지는 군 라벨이 없다',
    'gate5_judged': '0 (속 요약표 11개 기준) — 속 수준 계수가 없다. 단 같은 자료의 asv/ 표로 내려가니 진균 ⑤는 부분 가능이었다(11_진균_종수준.py · M. restricta AUC 0.165 q=0.0024). 세균은 그대로 불가',
}

print('\n' + '=' * 100)
print('관문별 적용 가능성 = 자료 요건표')
print('=' * 100)
print(f"  {'관문':16s} {'필요한 것':30s} {'라벨 불필요':>11s}  가용성")
for g, need, avail, lf in REQ:
    print(f"  {g:16s} {need:30s} {'예' if lf else '아니오':>11s}  {avail}")
print()
for k, v in out['summary'].items():
    if k.startswith('gate'):
        print(f"  {k:20s} {v}")

print('\n  → ①과 ②(심도)는 새 자료 %d개 전부에서 판정됐다. 판정표는 우리 네 편에 묶여 있지 않다.' % (n16 + nits))
print('  → 속 요약표를 쓰는 한 ③⑤는 어느 자료에서도 판정되지 않는다. 우리 방법이 아니라 자료 형식의 한계다.')
print('     ⚠️ 2026-09-19 — 같은 자료의 asv/ 표로 내려가니 진균 ⑤ 는 «부분 가능» 이었다(11_진균_종수준.py).')
print('        «③⑤ 는 어느 공개 자료에서도 판정되지 않는다» 고 말하면 틀린 말이다. ③ 은 아직 ASV 로 다시 재지 않았다.')
print('  → ④만 군 라벨을 요구한다. 그래서 「판별력을 모른다」가 자료의 문제라는 우리 결론과 맞는다.')


# ── 2차 규칙(코호트 간 배수)을 새 짝에 적용 ────────────────────────────────
# 595 는 얼굴 세척액·V3-V4 다. 같은 V3-V4 한국인 얼굴(02건강 볼)과 비교하면 배수가 작아야 하고,
# 채취법·영역이 다른 자료(42명 볼 swab V4-V5, 16비듬 두피 V3)와는 커야 한다.
print('\n' + '=' * 100)
print('2차 규칙(두 자료 중앙값 배수)을 새 짝에 적용 — 배수는 균의 양이 아니라 방법 차이를 잰다')
print('=' * 100)
base = load_595_as_cohort()
folds = {}
print(f"  {'상대 코호트':22s} {'영역':7s} " + ' '.join(f'{g[:14]:>15s}' for g in CAND))
for name, region, sex, role, loader in COHORTS_16S[1:]:
    try:
        coh = loader()
    except Exception:
        continue
    row = {}
    cells = []
    for g in CAND:
        m5 = st.median(base['val'][g]) if g in base['val'] else None
        mx = st.median(coh['val'][g]) if g in coh['val'] else None
        if m5 and mx:
            fold = max(m5 / mx, mx / m5)
            row[g] = round(fold, 1)
            cells.append(f'{fold:14.1f}배')
        else:
            row[g] = None
            cells.append(f"{'-':>15s}")
    folds[name] = {'region': region, 'fold_vs_595': row}
    print(f"  {name:22s} {region:7s} " + ' '.join(cells))
out['fold_vs_595'] = {
    'note': '595(얼굴 세척액·V3-V4) 중앙값 대비 배수. 2차 규칙의 문턱은 5배다.',
    'threshold': 5.0, 'pairs': folds}


# ── 속별 코호트 간 폭 — 절대 문턱을 쓸 수 없다는 근거 ──────────────────────
# 예상과 달랐던 것: 같은 V3-V4 한국인 얼굴인 02건강조차 595 기준 5배 문턱을 넘는다
# (Cutibacterium 이마 7.1배·볼 5.4배). 595 는 얼굴 전체 **세척액**이고 나머지는 **면봉**이다.
# 그런데 면봉끼리도 벌어진다. 그래서 결론은 "세척액이 특이하다"가 아니라 아래 한 줄이다.
print('\n' + '=' * 100)
print('속별 코호트 간 중앙값 폭 — 패널에 절대 문턱을 쓸 수 없는 이유')
print('=' * 100)
spread = {}
print(f"  {'후보':18s} {'최저 코호트':28s} {'최고 코호트':28s} {'폭':>9s}")
for g in CAND + FUNGI_CAND:
    vals = [(name, e['genera'][g]['median_pct'])
            for store in (out['bacteria'], out['fungi'])
            for name, e in store.items() if g in e['genera']]
    if len(vals) < 2:
        continue
    lo = min(vals, key=lambda t: t[1]); hi = max(vals, key=lambda t: t[1])
    ratio = hi[1] / lo[1] if lo[1] > 0 else None
    spread[g] = {'n_cohort': len(vals), 'min': {'cohort': lo[0], 'median_pct': lo[1]},
                 'max': {'cohort': hi[0], 'median_pct': hi[1]},
                 'max_over_min': round(ratio, 1) if ratio else None}
    print(f"  {g[:18]:18s} {lo[0][:22] + f' {lo[1]:.2f}%':28s} "
          f"{hi[0][:22] + f' {hi[1]:.2f}%':28s} {ratio:8.1f}배")
out['cross_cohort_spread'] = {
    'note': '같은 속의 중앙값이 연구마다 몇 배까지 벌어지는가. 채취 방식·영역·부위·파이프라인이 섞인 폭이다. '
            '어느 한쪽이 틀린 것이 아니라 셋 다 그 자료 안에서는 맞는 값이다.',
    'implication': '패널 기준을 "Cutibacterium 20% 이상" 같은 절대 문턱으로 쓸 수 없다. '
                   '기준은 반드시 같은 채취 프로토콜 안에서 정의돼야 하고, 프로토콜이 바뀌면 다시 정해야 한다.',
    'genera': spread}
print('  → 어느 한쪽이 틀린 값이 아니다. 각 자료 안에서는 셋 다 맞다.')
print('  → 그래서 패널 기준을 "Cutibacterium 20% 이상" 같은 **절대 문턱**으로 쓸 수 없다.')
print('     기준은 같은 채취 프로토콜 안에서 정의돼야 하고, 프로토콜이 바뀌면 다시 정해야 한다.')

save_json(out, os.path.join(TAB, '판정표_이식성.json'))
print(f"\n저장: {os.path.join(TAB, '판정표_이식성.json')}")
