"""② 교란 취약성 판정이 조성 효과(합이 1로 고정된 데이터의 착시)인지 확인한다.

상대풍부도는 사람마다 합이 1이다. 그래서 어떤 속이 늘면 나머지는 생물학적 이유 없이
기계적으로 줄어든다. ② 관문의 판정이 이 착시에서 온 것인지 가리려면, 같은 질문을
여러 표현으로 다시 물어봐야 한다. 세 가지로 물었다.

  표현 1  원래 상대비율               (현재 판정표가 쓰는 값)
  표현 2  로그비                      분모 = 후보·구강의심·미분류를 뺀 나머지 전체
  표현 3  구강 비중이 낮은 사람만      구강 10속 합계 < 20%인 부분집합

세 표현에서 방향과 크기가 유지되면 견고, 부호가 바뀌거나 사라지면 표현 의존이다.

결과: 분석표/조성효과_점검.json
"""
import json, math, os, statistics as st
from stats_util import save_json, TAB, CAND, ORAL, f, spearman, auc_p, auc_ci, load_595, load_42

EPS = 1e-5          # 로그비의 0 처리. 상대비율(합=1) 자료이므로 count용 0.5를 쓰지 않는다.
EPS_LADDER = [1e-5, 1e-4, 1e-3]   # 0 처리 민감도. 값을 바꿔도 ② 판정이 유지되는지 본다.
ORAL_CUT = 0.20     # 표현 3의 부분집합 기준

ids, genera, val, md, qc = load_595()
n = len(ids)
age = [f(md[i]['Age']) for i in ids]
oil = [f(md[i]['Oil_avg']) for i in ids]

oral_present = [g for g in ORAL if g in val]
osum = [sum(val[g][k] for g in oral_present) for k in range(n)]

# 피부 기준집합 = 전체 - 후보 3속 - 구강의심 - 미분류
EXC = set(CAND) | set(oral_present) | {'unidentified'}
ref = [sum(val[g][k] for g in genera if g not in EXC) for k in range(n)]

sub = [k for k in range(n) if osum[k] < ORAL_CUT]

def logratio(g, eps=EPS):
    return [math.log10((val[g][k] + eps) / (ref[k] + eps)) for k in range(n)]

out = {
    'n': n,
    'oral_genera_used': oral_present,
    'oral_sum_pct': {'median': round(100 * st.median(osum), 2),
                     'mean': round(100 * st.mean(osum), 2),
                     'max': round(100 * max(osum), 2),
                     'n_over_50pct': sum(1 for x in osum if x > 0.5)},
    'skin_reference_pct': {'median': round(100 * st.median(ref), 3),
                           'n_zero': sum(1 for x in ref if x == 0)},
    'subset_oral_under_20pct_n': len(sub),
    'eps': EPS,
    'gate2': {},
}

r, p = spearman(age, osum);          out['oral_sum_vs_age'] = {'rho': round(r, 3), 'p': p}
r, p = spearman(osum, val['Cutibacterium']); out['oral_sum_vs_cutibacterium'] = {'rho': round(r, 3), 'p': p}

print(f"구강의심 {len(oral_present)}속 합계: 중앙값 {out['oral_sum_pct']['median']}%  최대 {out['oral_sum_pct']['max']}%  "
      f"50% 초과 {out['oral_sum_pct']['n_over_50pct']}명")
print(f"피부 기준집합(분모) 중앙값 {out['skin_reference_pct']['median']}%\n")
print(f"{'후보':18s} {'표현1 원래%':>16s} {'표현2 로그비':>16s} {'표현3 구강<20%':>17s}  판정")

for g in CAND:
    lr = logratio(g)
    r1, p1 = spearman(age, val[g])
    r2, p2 = spearman(age, lr)
    r3, p3 = spearman([age[k] for k in sub], [val[g][k] for k in sub])
    o1, _ = spearman(oil, val[g])
    o2, op2 = spearman(oil, lr)
    o3, op3 = spearman([oil[k] for k in sub], [val[g][k] for k in sub])

    rs = [r1, r2, r3]
    same_sign = all(x > 0 for x in rs) or all(x < 0 for x in rs)
    strong = all(abs(x) >= 0.15 for x in rs)
    if same_sign and strong:
        verdict = '견고 — 세 표현 모두 같은 방향·비슷한 크기'
    elif not same_sign:
        verdict = '표현 의존 — 부호가 뒤집힘'
    else:
        verdict = '표현 의존 — 표현에 따라 사라짐'

    out['gate2'][g] = {
        'age': {'raw': {'rho': round(r1, 3), 'p': p1},
                'logratio': {'rho': round(r2, 3), 'p': p2},
                'subset': {'rho': round(r3, 3), 'p': p3, 'n': len(sub)}},
        'sebum': {'raw': {'rho': round(o1, 3)},
                  'logratio': {'rho': round(o2, 3), 'p': op2},
                  'subset': {'rho': round(o3, 3), 'p': op3}},
        'verdict': verdict,
    }
    print(f'{g:18s} {r1:+8.3f}(p={p1:5.1g}) {r2:+8.3f}(p={p2:5.1g}) {r3:+8.3f}(p={p3:5.1g})  {verdict}')

# ---- 0 처리(eps) 민감도 — 값을 바꿔도 ② 판정이 유지되는가 ----
print()
print('② 0 처리 민감도 — 로그비 eps 사다리')
print(f"{'후보':18s}" + ''.join(f'{("eps=" + str(e)):>16s}' for e in EPS_LADDER) + '   판정')
out['gate2_eps_ladder'] = {'eps_values': EPS_LADDER, 'genera': {}}
for g in CAND:
    rr = []
    for e in EPS_LADDER:
        r, pv = spearman(age, logratio(g, e))
        rr.append({'eps': e, 'rho': round(r, 3), 'p': pv})
    rhos = [x['rho'] for x in rr]
    stable = (all(x > 0 for x in rhos) or all(x < 0 for x in rhos)) and              (max(rhos) - min(rhos)) < 0.05
    out['gate2_eps_ladder']['genera'][g] = {'by_eps': rr,
                                            'rho_range': round(max(rhos) - min(rhos), 4),
                                            'stable': stable}
    print(f'{g:18s}' + ''.join(f'{x["rho"]:+16.3f}' for x in rr)
          + ('   불변' if stable else '   ** 달라짐 **'))
print('  → eps 를 100배 바꿔도 ② 판정이 유지되면 0 처리 선택이 결론을 만들지 않았다는 뜻이다.')

# ---- ④ 판별력도 로그비로 다시 본다 (42명) ----
info42, v42 = load_42()
grp = [r['group'] for r in info42]
m = len(info42)
exc42 = set(CAND) | set(g for g in ORAL if g in v42) | {'Unassigned'}
ref42 = [sum(v42[g][k] for g in v42 if g not in exc42) for k in range(m)]

out['gate4_logratio'] = {}
print('\n④ 판별력 — 로그비로 다시 계산 (42명)')
for g in CAND:
    lr = [math.log10((v42[g][k] + EPS) / (ref42[k] + EPS)) for k in range(m)]
    S = [lr[k] for k in range(m) if grp[k] == 'Sensitive']
    N = [lr[k] for k in range(m) if grp[k] == 'Non-sensitive']
    a, p = auc_p(S, N)
    lo, hi = auc_ci(S, N)
    out['gate4_logratio'][g] = {'auc': round(a, 3), 'ci': [round(lo, 3), round(hi, 3)], 'p': round(p, 4),
                               'crosses_half': lo <= 0.5 <= hi}
    print(f'  {g:18s} AUC {a:.3f} [{lo:.3f}-{hi:.3f}] p={p:.3f}  0.5 걸침 {"예" if lo <= 0.5 <= hi else "아니오"}')

save_json(out, os.path.join(TAB, '조성효과_점검.json'))
print(f"\n저장: {os.path.join(TAB, '조성효과_점검.json')}")
