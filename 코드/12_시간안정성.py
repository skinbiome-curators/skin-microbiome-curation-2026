"""후보균이 **같은 사람 안에서 시간이 지나도 유지되는가** — 2026-09-19 신규.

15차 §10 의 2급 한계 하나를 정면으로 잰다.

    "한 시점 한 번 측정이다. 같은 사람 재측정이 없어 지표의 시간 안정성을 모른다."

**모르는 게 아니라 안 쟀던 것이다.** 수업 공개 자료에 반복측정 설계가 둘 있다.

  · `15_season` 이마 — **24명 × 4계절**(2014-01·04·08·11). 관찰 연구, 1년 추적
  · `07_skincare` 볼 — **25명 × 3시점**(0W·2W·4W). **화장품 중재** 연구

재는 것: **ICC(1)** — 전체 변동 중 «사람 사이» 차이가 차지하는 몫.

    ICC = (MSB - MSW) / (MSB + (k-1) * MSW)

  ICC 가 높다  = 그 사람의 값이 **그 사람의 특징**이다 → 패널 마커로 쓸 수 있다
  ICC 가 낮다  = 같은 사람도 잴 때마다 흔들린다 → 한 번 재서 판정하면 안 된다

조성 착시를 피하려고 **원래 %** 와 **로그비**(후보를 뺀 나머지 대비) 두 표현으로 각각 낸다.
§5 에서 ②를 세 표현으로 다시 물었던 것과 같은 이유다.

⚠️ 반드시 같이 말할 것
  1. **두 코호트 다 우리 대상이 아니다.** 15계절은 미국·중국 가구 이마(V4), 07 은 화장품 중재다.
     한국인 여성 얼굴에서 잰 값이 아니므로 **참고값**이다.
  2. **ICC 가 낮다고 «측정이 나쁘다»가 아니다.** 계절에 따라 진짜로 변했을 수도 있다.
     어느 쪽이든 **패널에는 같은 문제** — 한 번 재서 판정할 수 없다는 뜻이다.
  3. 07 은 **중재 연구**라 사람 안 변동에 처치 효과가 섞인다. 그래서 07 의 ICC 는 **하한**이다.

결과: 분석표/시간안정성.json
"""
import csv, math, os, statistics as st
from stats_util import save_json, TAB, RAW, CAND, f

MAL = 'g__Malassezia'
EPS = 1e-5          # 로그비의 0 처리 — 4_조성효과_점검.py 와 같은 값


def read_tsv(p):
    with open(p, encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def load(tag, site):
    si = {r['run']: r for r in read_tsv(os.path.join(RAW, tag + '_sample_info.tsv'))}
    rd = csv.reader(open(os.path.join(RAW, tag + '_genus_pct.tsv'), encoding='utf-8'), delimiter='\t')
    hdr = next(rd)
    col = {g: i for i, g in enumerate(hdr)}
    rows = [r for r in rd if r[0] in si
            and si[r[0]].get('body_site') == site and si[r[0]].get('group') != '음성대조']
    info = [si[r[0]] for r in rows]
    val = {g: [f(r[col[g]]) or 0.0 for r in rows] for g in hdr[1:]}
    return info, val


def icc1(groups):
    """일원 변량 ICC. groups = [[같은 사람의 반복값...], ...] (반복수 동일 가정)"""
    k = len(groups[0])
    n = len(groups)
    allv = [x for g in groups for x in g]
    gm = sum(allv) / len(allv)
    msb = sum(k * (sum(g) / k - gm) ** 2 for g in groups) / (n - 1)
    msw = sum((x - sum(g) / k) ** 2 for g in groups for x in g) / (n * (k - 1))
    den = msb + (k - 1) * msw
    if den <= 0:
        return None
    return (msb - msw) / den


def icc_ci(groups, B=2000, seed=20260918):
    """사람 단위 부트스트랩 95% CI. n=24 라 ICC 추정이 부정확하므로 반드시 같이 낸다."""
    import random
    rng = random.Random(seed)
    n = len(groups)
    vals = []
    for _ in range(B):
        samp = [groups[rng.randrange(n)] for _ in range(n)]
        v = icc1(samp)
        if v is not None:
            vals.append(v)
    if len(vals) < 100:
        return None, None
    vals.sort()
    return vals[int(0.025 * len(vals))], vals[int(0.975 * len(vals))]


def _selftest():
    """구현이 맞는지 극단 두 경우로 확인한다."""
    same = [[5.0, 5.0, 5.0], [9.0, 9.0, 9.0], [1.0, 1.0, 1.0]]      # 사람 안 변동 0 -> 1.0
    noise = [[1.0, 9.0, 5.0], [5.0, 1.0, 9.0], [9.0, 5.0, 1.0]]     # 사람 사이 차이 0 -> <=0
    a, b = icc1(same), icc1(noise)
    assert a is not None and abs(a - 1.0) < 1e-9, a
    assert b is not None and b <= 0.01, b
    return round(a, 3), round(b, 3)


def logratio(v, others_sum):
    """후보를 뺀 나머지 대비 로그비."""
    return [math.log((a + EPS) / (b + EPS)) for a, b in zip(v, others_sum)]


def analyse(label, tag, site, cands, note):
    info, val = load(tag, site)
    subj = [r['subject'] for r in info]
    order = {}
    for i, s in enumerate(subj):
        order.setdefault(s, []).append(i)
    ks = sorted({len(v) for v in order.values()})
    if len(ks) != 1:
        return {'error': f'반복수가 균일하지 않다: {ks}'}
    k = ks[0]
    n = len(order)
    # 후보 외 나머지 합 (로그비의 분모)
    rest = [sum(val[g][i] for g in val if g not in cands) for i in range(len(subj))]

    entry = {'cohort': tag, 'site': site, 'n_subject': n, 'n_repeat': k,
             'n_sample': len(subj), 'note': note, 'genera': {}}
    print(f'\n■ {label} — {n}명 × {k}회 = {len(subj)}시료 · {note}')
    print(f"  {'후보':18s} {'중앙값':>9s} {'ICC(원래%)':>7s} {'95% CI':>13s} "
          f"{'ICC(로그비)':>8s} {'95% CI':>13s} {'변동폭':>11s}  판정")
    for g in cands:
        if g not in val:
            print(f'  {g[:18]:18s} — 이 자료에 없는 속')
            continue
        v = val[g]
        lr = logratio(v, rest)
        gp_raw = [[v[i] for i in ix] for ix in order.values()]
        gp_lr = [[lr[i] for i in ix] for ix in order.values()]
        a, b = icc1(gp_raw), icc1(gp_lr)
        a_lo, a_hi = icc_ci(gp_raw)
        b_lo, b_hi = icc_ci(gp_lr)
        # 같은 사람 안에서 최대/최소 비 (0 은 EPS 로 보정)
        folds = [max(x + EPS for x in gpr) / min(x + EPS for x in gpr) for gpr in gp_raw]
        fold_med = st.median(folds)
        worst = max(a or 0.0, b or 0.0)
        verdict = ('안정' if worst >= 0.75 else '보통' if worst >= 0.5
                   else '취약' if worst >= 0.4 else '불안정')
        entry['genera'][g] = {
            'median_pct': round(st.median(v), 3),
            'icc_raw': round(a, 3) if a is not None else None,
            'icc_raw_ci': [round(a_lo, 3), round(a_hi, 3)] if a_lo is not None else None,
            'icc_logratio': round(b, 3) if b is not None else None,
            'icc_logratio_ci': [round(b_lo, 3), round(b_hi, 3)] if b_lo is not None else None,
            'ci_upper_below_0_5': (a_hi is not None and b_hi is not None
                                   and max(a_hi, b_hi) < 0.5),
            'within_subject_fold_median': round(fold_med, 1),
            'within_subject_fold_max': round(max(folds), 1),
            'verdict': verdict}
        flag = '  CI 상한도 0.5 미만' if (a_hi is not None and b_hi is not None
                                          and max(a_hi, b_hi) < 0.5) else ''
        print(f"  {g[:18]:18s} {st.median(v):8.2f}% {a:7.3f} [{a_lo:.2f}-{a_hi:.2f}] "
              f"{b:7.3f} [{b_lo:.2f}-{b_hi:.2f}] {fold_med:9.1f}배  {verdict}{flag}")
    return entry


out = {'question': '후보균 값이 같은 사람 안에서 시간이 지나도 유지되는가',
       'metric': 'ICC(1) — 전체 변동 중 사람 사이 차이의 몫. 높을수록 그 사람의 특징이다.',
       'thresholds': {'안정': '>=0.75', '보통': '0.50-0.75', '취약': '0.40-0.50', '불안정': '<0.40'},
       'threshold_note': '문턱은 재현성 연구의 관례를 따랐고 사전등록한 것이 아니다.',
       'caveats': [
           '두 코호트 다 우리 대상(한국인 여성 얼굴)이 아니다. 참고값이다.',
           'ICC 가 낮다고 측정이 나쁘다는 뜻이 아니다 — 계절에 따라 진짜로 변했을 수 있다. '
           '어느 쪽이든 패널에는 같은 문제다: 한 번 재서 판정할 수 없다.',
           '07 은 화장품 중재 연구라 사람 안 변동에 처치 효과가 섞인다. 07 의 ICC 는 하한이다.'],
       'cohorts': {}}

_st = _selftest()
out['implementation_selftest'] = {
    'note': '사람 안 변동이 0 이면 ICC=1, 사람 사이 차이가 0 이면 ICC<=0 이어야 한다.',
    'no_within_variation': _st[0], 'no_between_variation': _st[1], 'passed': True}
print(f'구현 자체 검사 — 사람안 변동 0 → ICC {_st[0]} · 사람사이 차이 0 → ICC {_st[1]}  (통과)')
print('=' * 104)
print('시간 안정성 — 같은 사람을 여러 번 재면 후보균 값이 유지되는가')
print('=' * 104)

out['cohorts']['15계절 · 이마 (세균)'] = analyse(
    '15계절 · 이마 · 세균', '15_season_16S', 'Forehead', CAND,
    '2014년 1·4·8·11월 · 관찰 연구')
out['cohorts']['15계절 · 이마 (진균)'] = analyse(
    '15계절 · 이마 · 진균', '15_season_ITS', 'Forehead', [MAL],
    '2014년 1·4·8·11월 · 관찰 연구')
out['cohorts']['07화장품 · 볼 (세균)'] = analyse(
    '07화장품 · 볼 · 세균', '07_skincare_16S', 'Cheek', CAND,
    '0W·2W·4W · ⚠️ 중재 연구이므로 ICC 는 하한')

print()
print('읽는 법')
print('  ICC 가 높다 → 그 사람의 값이 그 사람의 특징이다. 한 번 재도 그 사람을 대표한다')
print('  ICC 가 낮다 → 같은 사람도 잴 때마다 흔들린다. 한 번 재서 패널 판정을 내리면 안 된다')
print('  원래 % 와 로그비 중 하나만 높으면 "표현 의존" 이다 — ② 관문과 같은 성질')

save_json(out, os.path.join(TAB, '시간안정성.json'))
print(f"\n저장: {os.path.join(TAB, '시간안정성.json')}")
