"""진균(ITS) 후보의 5관문 판정.

세균 3속은 ④ 판별력에서 전부 막혔다. 같은 42명의 진균 자료를 같은 문턱 «숫자»로 판정한다.

⚠️ "같은 규칙"이 아니다. 문턱 상수(검출 90% · 중앙 0.5% · 배수 5배)는 같지만 적용 범위가 다르다.
   세균: 1차 = 양쪽 코호트(595+42) 검출 AND 595 중앙값  →  2차 = 두 코호트 배수 < 5
   진균: 1차 = 42명 검출 AND 42명 중앙값                 →  2차 없음 (595에 진균 자료가 없어 배수를 못 잼)
   따라서 진균 후보는 세균과 같은 심사를 통과한 것이 아니다.
사람은 늘지 않으므로 신뢰구간이 좁아지지는 않는다. 측정 축이 하나 더 생기는 것이다.

⚠️ 채취 조건이 자료마다 다르다. 어디에 영향이 있고 없는지 구분해서 쓴다.

  자료              부위              영역      성별
  03 민감피부 ITS   오른쪽 볼 swab    ITS1      여 42
  01 노화 ITS       볼 30 / 이마 28   ITS1      여 58
  595명 (AAML)      얼굴 전체 세척액  16S V3-V4 여 595   ← 진균 자료 없음

  ④ 판별력: 민감군·비민감군이 같은 연구·같은 부위·같은 영역이므로 채취 차이의 영향을 받지 않는다.
  ② 나이   : 03 내부(같은 조건)로 먼저 보고, 01 노화로 교차 확인한다. 01은 다른 연구이므로
             '재현'이 아니라 '조건이 다른 자료에서의 일관성'이며, 부위를 볼로 맞춰 비교한다.
  ①③⑤    : 595명에 진균이 없어 세균만큼 넓게 볼 수 없다. 이 한계를 판정에 적는다.

결과: results/진균_판별력.json
"""
import csv, json, math, os, sys, statistics as st
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_util import save_json, TAB, RAW, f, spearman, auc_p, auc_ci, bh

DET_MIN, MED_MIN, FOLD_MAX = 90.0, 0.5, 5.0
EPS = 1e-3          # ITS 상대비율(%) 자료의 0 처리

def load(tag, table):
    si = {r['run']: r for r in csv.DictReader(open(os.path.join(RAW, f'{tag}_sample_info.tsv')), delimiter='\t')}
    rd = csv.reader(open(os.path.join(RAW, table)), delimiter='\t')
    hdr = next(rd); rows = [r for r in rd]
    col = {g: i for i, g in enumerate(hdr)}
    return [si[r[0]] for r in rows], {g: [f(r[col[g]]) or 0.0 for r in rows] for g in hdr[1:]}, hdr

info3, v3, hdr3 = load('03_sensitive_ITS', '03_sensitive_ITS_genus_pct.tsv')
info1, v1, hdr1 = load('01_aging_ITS',     '01_aging_ITS_genus_pct.tsv')
n3 = len(info3)
grp  = [r['group'] for r in info3]
age3 = [f(r['age']) for r in info3]
dep3 = [f(r['read_depth_final']) for r in info3]

out = {'cohort': {'n': n3, 'sensitive': grp.count('Sensitive'), 'non_sensitive': grp.count('Non-sensitive'),
                  'site': '오른쪽 볼 swab', 'region': 'ITS1', 'sex': '여성'},
       'sampling_note': '④는 같은 연구·부위·영역 안의 군 비교이므로 채취 차이의 영향을 받지 않는다. '
                        '01 노화와의 비교는 다른 연구이므로 일관성 확인이며, 부위를 볼로 맞췄다.',
       'limitation': '595명에는 진균 자료가 없다. 세균 후보와 달리 ①②③⑤를 595명으로 넓게 볼 수 없다.',
       'candidates': {}}

# ---------- 1차 규칙 통과 진균 ----------
stage1 = []
for g, v in v3.items():
    if g in ('Unassigned', 'unidentified'):
        continue
    det = 100 * sum(1 for x in v if x > 0) / n3
    if det >= DET_MIN and st.median(v) >= MED_MIN:
        stage1.append(g)

pv = []
for g in stage1:
    v = v3[g]
    S = [v[k] for k in range(n3) if grp[k] == 'Sensitive']
    N = [v[k] for k in range(n3) if grp[k] == 'Non-sensitive']
    a, p = auc_p(S, N); lo, hi = auc_ci(S, N)

    # ① 검출·심도
    det = 100 * sum(1 for x in v if x > 0) / n3
    r_dep, p_dep = spearman(dep3, v)

    # ② 나이 — 03 내부
    r_age, p_age = spearman(age3, v)
    bands = {}
    for blo, bhi in [(20, 30), (31, 55)]:
        idx = [k for k in range(n3) if blo <= (age3[k] or 0) <= bhi]
        bS = [v[k] for k in idx if grp[k] == 'Sensitive']
        bN = [v[k] for k in idx if grp[k] == 'Non-sensitive']
        if len(bS) >= 4 and len(bN) >= 4:
            ba, bp = auc_p(bS, bN)
            bands[f'{blo}-{bhi}세'] = {'auc': round(ba, 3), 'p': round(bp, 4), 'n_S': len(bS), 'n_N': len(bN)}

    # ② 조성 — 로그비
    ref = [sum(v3[o][k] for o in v3 if o not in (g, 'Unassigned', 'unidentified')) for k in range(n3)]
    lr = [math.log10((v[k] + EPS) / (ref[k] + EPS)) for k in range(n3)]
    lS = [lr[k] for k in range(n3) if grp[k] == 'Sensitive']
    lN = [lr[k] for k in range(n3) if grp[k] == 'Non-sensitive']
    la, lp = auc_p(lS, lN); llo, lhi = auc_ci(lS, lN)

    # 방법 의존성 + 나이 교차 확인 — 01 노화 ITS, 부위를 볼로 맞춤
    cross = None
    if g in v1:
        cheek = [k for k, r in enumerate(info1) if r['body_site'] == 'Cheek']
        vc = [v1[g][k] for k in cheek]
        gc = [info1[k]['group'] for k in cheek]
        O = [vc[k] for k in range(len(cheek)) if gc[k] == 'O']
        Y = [vc[k] for k in range(len(cheek)) if gc[k] == 'Y']
        m1 = st.median(vc); m3 = st.median(v)
        fold = max(m1 / m3, m3 / m1) if m1 > 0 and m3 > 0 else float('inf')
        ca, cp = auc_p(O, Y) if O and Y else (None, None)
        cross = {'site': 'Cheek만', 'n': len(cheek), 'median_01': round(m1, 3), 'median_03': round(m3, 3),
                 # ⚠️ 관문이 아니다. 진균은 2차 규칙을 «치르지 않았고» 이 배수는 선정이 끝난 뒤 붙인
                 #    참고값이다. 옛 이름 method_stable 은 판정처럼 읽혀 2026-09-20 에 바꿨다.
                 'fold': round(fold, 2), 'fold_under_5_reference_only': fold < FOLD_MAX,
                 'auc_old_vs_young': round(ca, 3) if ca is not None else None,
                 'p': round(cp, 4) if cp is not None else None,
                 'median_O': round(st.median(O), 2) if O else None,
                 'median_Y': round(st.median(Y), 2) if Y else None}

    out['candidates'][g] = {
        'gate1': {'detect_pct': round(det, 1), 'median_pct': round(st.median(v), 3),
                  'rho_depth': round(r_dep, 3), 'p_depth': round(p_dep, 3),
                  'contamination_flag': r_dep <= -0.30},
        'gate2_age': {'rho_in_03': round(r_age, 3), 'p': round(p_age, 3), 'age_bands': bands},
        'gate2_composition': {'auc_logratio': round(la, 3), 'ci': [round(llo, 3), round(lhi, 3)],
                              'p': round(lp, 4), 'crosses_half': llo <= 0.5 <= lhi},
        'gate4': {'auc': round(a, 3), 'ci': [round(lo, 3), round(hi, 3)], 'p': round(p, 4),
                  'crosses_half': lo <= 0.5 <= hi,
                  'median_sensitive': round(st.median(S), 2), 'median_non_sensitive': round(st.median(N), 2)},
        'cross_cohort_01_aging': cross,
        # 2026-09-19 — 이 문자열이 바뀌었다. 옛 값은
        #   '불가 — 속 수준 표만 제공됨. 원저는 M. restricta 감소를 보고했으나
        #    restricta·globosa·sympodialis가 한 칸에 합쳐져 방향을 말할 수 없음'
        # 이었는데, 같은 자료의 ASV 표에 종이 나뉘어 있어 전제가 틀렸다. 11_진균_종수준.py 참조.
        # 분석표/진균_판별력.json 의 같은 필드도 이 문자열과 똑같이 고쳐 두었다.
        # ⚠️ 2026-09-19 — 후보별로 분기한다. 이 문자열을 모든 후보에 공유하면
        #    종 수준 분석을 하지 않은 속까지 「부분 가능」이 된다.
        #    11_진균_종수준.py 는 GENUS='g__Malassezia' 하나만 본다.
        #    부정 판정은 공유해도 되지만 긍정 판정은 공유하면 안 된다.
        'gate5': (
            '부분 가능(진균 한정) — 속 요약표에서만 한 칸이었다. 같은 자료의 ASV 표에는 종이 나뉘어 있고'
            '(Malassezia ASV 116개 중 74개에 종 라벨) 원저가 보고한 M. restricta 감소가 재현된다'
            '(AUC 0.165 [0.050-0.307] q=0.0024). 통과는 아니다 — 종 라벨의 36%가 여전히 미분류이고 '
            'restricta 가 속(0.153)보다 나은 지표라는 뜻도 아니다. 생성기 11_진균_종수준.py'
        ) if g == 'g__Malassezia' else (
            '불가 — 이 속은 종 수준 분석을 하지 않았다. 속 요약표만 있다. '
            'ASV 표로 내려간 것은 Malassezia 뿐이다(11_진균_종수준.py).'
        ),
    }
    pv.append(p)

for g, q in zip(stage1, bh(pv)):
    out['candidates'][g]['gate4']['q_bh'] = round(q, 4)

# ---------- 출력 ----------
print(f"03 민감피부 ITS — 민감 {grp.count('Sensitive')} / 비민감 {grp.count('Non-sensitive')}  ·  오른쪽 볼 swab · ITS1")
print(f"1차 규칙(검출 ≥{DET_MIN:.0f}% AND 중앙값 ≥{MED_MIN}%) 통과 진균 {len(stage1)}속\n")
for g in stage1:
    c = out['candidates'][g]
    g4, g1, g2 = c['gate4'], c['gate1'], c['gate2_composition']
    print(f"■ {g}")
    print(f"   ① 검출 {g1['detect_pct']}%  중앙값 {g1['median_pct']}%  심도 ρ={g1['rho_depth']:+.3f} "
          f"→ {'오염 의심' if g1['contamination_flag'] else '해당 없음'}")
    print(f"   ② 나이(03 내부) ρ={c['gate2_age']['rho_in_03']:+.3f} (p={c['gate2_age']['p']})")
    for k, b in c['gate2_age']['age_bands'].items():
        print(f"        {k}: AUC {b['auc']} (p={b['p']}, 민감 {b['n_S']}/비민감 {b['n_N']})")
    print(f"   ② 조성(로그비) AUC {g2['auc_logratio']} [{g2['ci'][0]}-{g2['ci'][1]}] "
          f"→ 0.5 {'걸침' if g2['crosses_half'] else '안 걸침'}")
    print(f"   ④ 판별력 AUC {g4['auc']} [{g4['ci'][0]}-{g4['ci'][1]}] p={g4['p']} q={g4['q_bh']} "
          f"→ {'근거 부족' if g4['crosses_half'] else '근거 있음'}")
    print(f"        민감군 {g4['median_sensitive']}% vs 비민감군 {g4['median_non_sensitive']}%")
    cc = c['cross_cohort_01_aging']
    if cc:
        print(f"   ↔ 01 노화 ITS (볼 {cc['n']}개): 중앙값 {cc['median_01']}% vs 03의 {cc['median_03']}% "
              f"= {cc['fold']}배 (참고값 — 관문 아님, 5배 미만 {cc['fold_under_5_reference_only']})")
        print(f"        고연령 {cc['median_O']}% vs 저연령 {cc['median_Y']}%  AUC(고>저) {cc['auc_old_vs_young']} (p={cc['p']})")
    print(f"   ⑤ {c['gate5'][:52]}…\n")

save_json(out, os.path.join(TAB, '진균_판별력.json'))
print(f"저장: {os.path.join(TAB, '진균_판별력.json')}")
