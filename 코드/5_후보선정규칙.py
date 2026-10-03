"""후보 3속을 어떤 규칙으로 골랐는지 코드로 고정하고, 규칙을 통과한 속을 전부 공개한다.

⚠️ 이 규칙은 결과를 본 뒤에 명문화한 사후 규칙이다. 사전등록된 규칙이 아니다.
   그렇게 표시하지 않고 "미리 정했다"고 말하면 안 된다.

규칙
  1차 (존재)   양쪽 코호트 검출률 ≥ 90%  AND  595명 중앙값 ≥ 0.5%
               'unidentified' / 'Unassigned'는 속이 아니므로 제외
  2차 (방법)   두 코호트의 중앙값 배수가 5배 미만
               595는 얼굴 전체 세척액·16S V3-V4, 42명은 오른쪽 볼 swab·16S V4-V5다.
               채취 부위·시퀀싱 영역·파이프라인이 모두 다르다. 이 차이를 건너 값이
               유지되지 않는 속은 '균의 양'이 아니라 '측정 방법'을 재고 있는 것이므로
               패널 후보에서 뺀다. 구강 유래라고 단정하지 않는다 — 세 요인이 섞여 있다.

1차를 통과했으나 2차에서 탈락한 속도 판별력(AUC)까지 계산해 전부 기록한다.
유의해 보이는 것을 빼고 보고하지 않는 것이 선택적 보고이기 때문이다.

결과: 분석표/후보선정규칙.json
"""
import json, os, statistics as st
from stats_util import save_json, TAB, CAND, f, auc_p, auc_ci, bh, load_595, load_42

DET_MIN    = 0.90   # 양쪽 검출률 하한
MED_MIN    = 0.50   # 595 중앙값 하한 (%)
FOLD_MAX   = 5.0    # 두 코호트 중앙값 배수 상한
NOT_GENUS  = {'unidentified', 'Unassigned', 'metagenome'}

ids, genera, val595, md, qc = load_595()
info42, val42 = load_42()
n5, n4 = len(ids), len(info42)
grp = [r['group'] for r in info42]

def pct(v):        return [100 * x for x in v]
def det(v):        return 100 * sum(1 for x in v if x > 0) / len(v)

# ---------- 1차: 존재 규칙 ----------
stage1 = []
for g in genera:
    if g in NOT_GENUS or g not in val42:
        continue
    v5, v4 = pct(val595[g]), val42[g]
    d5, d4 = det(v5), det(v4)
    m5, m4 = st.median(v5), st.median(v4)
    if d5 >= DET_MIN * 100 and d4 >= DET_MIN * 100 and m5 >= MED_MIN:
        fold = max(m5 / m4, m4 / m5) if m4 > 0 and m5 > 0 else float('inf')
        stage1.append({'genus': g, 'detect_595': round(d5, 1), 'detect_42': round(d4, 1),
                       'median_595': round(m5, 3), 'median_42': round(m4, 3),
                       'fold': round(fold, 1)})
stage1.sort(key=lambda r: -r['median_595'])

# ---------- 2차: 방법 의존성 ----------
for r in stage1:
    r['stage2_pass'] = r['fold'] < FOLD_MAX
    r['verdict'] = '후보 채택' if r['stage2_pass'] else f"제외 — 두 자료 중앙값 {r['fold']}배 차이"

# ---------- 1차 통과 속 전부의 판별력 (선택적 보고 금지) ----------
pv = []
for r in stage1:
    v = val42[r['genus']]
    S = [v[k] for k in range(n4) if grp[k] == 'Sensitive']
    N = [v[k] for k in range(n4) if grp[k] == 'Non-sensitive']
    a, p = auc_p(S, N)
    lo, hi = auc_ci(S, N)
    r['auc_42'] = round(a, 3); r['ci_42'] = [round(lo, 3), round(hi, 3)]; r['p_42'] = round(p, 4)
    pv.append(p)
for r, q in zip(stage1, bh(pv)):
    r['q_bh_42'] = round(q, 4)

passed  = [r['genus'] for r in stage1 if r['stage2_pass']]
dropped = [r for r in stage1 if not r['stage2_pass']]

# ---------- 1차 탈락했지만 595에서 흔한 속도 공개한다 ----------
# 595 중앙값은 기준을 넘지만 42명 검출률이 모자라 1차에서 떨어진 속.
# 판별력이 유의해 보이는 것이 섞여 있을 수 있으므로 숨기지 않고 기록한다.
seen = {r['genus'] for r in stage1}
disclosed = []
for g in genera:
    if g in NOT_GENUS or g in seen or g not in val42:
        continue
    v5, v4 = pct(val595[g]), val42[g]
    if st.median(v5) < MED_MIN:
        continue
    S = [v4[k] for k in range(n4) if grp[k] == 'Sensitive']
    N = [v4[k] for k in range(n4) if grp[k] == 'Non-sensitive']
    a, p = auc_p(S, N)
    lo, hi = auc_ci(S, N)
    m5d, m4d = st.median(v5), st.median(v4)
    fold_d = max(m5d / m4d, m4d / m5d) if m4d > 0 and m5d > 0 else float('inf')
    disclosed.append({'genus': g, 'detect_595': round(det(v5), 1), 'detect_42': round(det(v4), 1),
                      'median_595': round(m5d, 3), 'median_42': round(m4d, 3),
                      'fold': round(fold_d, 1),
                      'auc_42': round(a, 3), 'ci_42': [round(lo, 3), round(hi, 3)], 'p_42': round(p, 4),
                      'why_dropped': ' · '.join(
                          ([f"595 검출률 {det(v5):.1f}% < {DET_MIN*100:.0f}%"] if det(v5) < DET_MIN*100 else [])
                          + ([f"42명 검출률 {det(v4):.1f}% < {DET_MIN*100:.0f}%"] if det(v4) < DET_MIN*100 else []))})
# 1차 통과 4속 + 공개 대상을 합쳐 다중검정 보정을 다시 건다
allp = [r['p_42'] for r in stage1] + [r['p_42'] for r in disclosed]
allq = bh(allp)
for r, q in zip(stage1 + disclosed, allq):
    r['q_bh_all'] = round(q, 4)

out = {
    'rule': {'note': '결과를 본 뒤 명문화한 사후 규칙. 사전등록 규칙이 아님.',
             'stage1': f'양쪽 검출률 >= {DET_MIN*100:.0f}% AND 595 중앙값 >= {MED_MIN}%',
             'stage2': f'두 코호트 중앙값 배수 < {FOLD_MAX}배',
             'excluded_labels': sorted(NOT_GENUS),
             'cohort_difference': '595 = 얼굴 전체 세척액 · 16S V3-V4 / 42명 = 오른쪽 볼 swab · 16S V4-V5'},
    'n_genus_595': len(genera), 'n_genus_42': len(val42),
    'stage1_passed': len(stage1), 'stage2_passed': len(passed),
    'selected': passed,
    'matches_hardcoded_candidates': sorted(passed) == sorted(CAND),
    'table': stage1,
    'disclosed_stage1_fail': disclosed,
    'multiplicity': f'1차 통과 {len(stage1)}속 + 공개 {len(disclosed)}속 = {len(stage1)+len(disclosed)}가설에 BH 보정',
}
save_json(out, os.path.join(TAB, '후보선정규칙.json'))

print(f'1차 통과 {len(stage1)}속 → 2차 통과 {len(passed)}속')
print(f'기존 하드코딩 3속과 일치: {"예" if out["matches_hardcoded_candidates"] else "아니오"}\n')
print(f"{'속':18s} {'595검출':>7s} {'595중앙':>8s} {'42중앙':>8s} {'배수':>7s} {'AUC(42)':>8s} {'q(BH)':>7s}  판정")
for r in stage1:
    print(f"{r['genus']:18s} {r['detect_595']:6.1f}% {r['median_595']:8.3f} {r['median_42']:8.3f} "
          f"{r['fold']:7.1f} {r['auc_42']:8.3f} {r['q_bh_42']:7.3f}  {r['verdict']}")
print('\n제외된 속:', ', '.join(f"{r['genus']}({r['fold']}배)" for r in dropped) or '없음')
if disclosed:
    print('\n[공개] 1차 탈락했으나 595에서 흔한 속 — 숨기지 않고 기록')
    for r in disclosed:
        print(f"  {r['genus']:16s} 595 {r['median_595']:7.3f}% / 42 {r['median_42']:6.3f}%  "
              f"AUC {r['auc_42']:.3f} [{r['ci_42'][0]:.3f}-{r['ci_42'][1]:.3f}] p={r['p_42']:.4f} "
              f"q={r['q_bh_all']:.3f}  ← {r['why_dropped']}")
