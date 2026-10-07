"""설문 라벨(Baumann 민감/저항)이 실제 피부 측정값과 연결되는지 검증한다.

HTML 그림 1의 값을 만드는 스크립트. 결과는 results/설문라벨_검증.json.

핵심 질문: 595명의 설문 라벨을 민감피부 진단의 정답표로 쓸 수 있는가?
→ 쓸 수 없다. 붉은기(a*)를 전혀 가르지 못하고 나이만 강하게 가른다.
"""
import json, os, statistics as st
from stats_util import save_json, TAB, f, auc_p, bh, load_595

# (열 이름, 표시 이름) — 붉은기는 홍반계(Mexameter) 측정이 아니라 색측정 a*축이다.
MEASURES = [
    ('Age',                  '나이'),
    ('Oil_avg',              '피지'),
    ('TEWL_Cheek_avg',      '수분손실(TEWL, 볼)'),
    ('Hydration_Cheek_avg',  '수분(볼)'),
    ('Color_Cheek_ITA',      '피부 밝기(ITA)'),
    ('Color_Cheek_a_value',  '붉은기(CIELAB a*)'),
]

ids, genera, val, md, qc = load_595()
sr = [md[i]['baumann_sr'] for i in ids]

out = {'n': len(ids),
       'n_S': sr.count('S'), 'n_R': sr.count('R'),
       'note': '붉은기는 홍반계 측정이 아니라 색측정 기기의 CIELAB a* 축이다.',
       'measures': {}}

pv = []
for col, label in MEASURES:
    S = [f(md[i][col]) for i in ids if md[i]['baumann_sr'] == 'S' and f(md[i][col]) is not None]
    R = [f(md[i][col]) for i in ids if md[i]['baumann_sr'] == 'R' and f(md[i][col]) is not None]
    a, p = auc_p(S, R)
    pv.append(p)
    out['measures'][label] = {
        'column': col, 'auc': round(a, 4), 'discrimination': round(abs(a - 0.5), 4),
        'p': p, 'median_S': round(st.median(S), 3), 'median_R': round(st.median(R), 3),
        'n_S': len(S), 'n_R': len(R),
    }

for (col, label), q in zip(MEASURES, bh([m['p'] for m in out['measures'].values()])):
    out['measures'][label]['q_bh'] = round(q, 4)

# 붉은기가 나이 층화 후에도 구별되지 않는지 확인
strat = {}
for lo, hi in [(20, 39), (40, 54), (55, 79)]:
    sub = [i for i in ids if lo <= (f(md[i]['Age']) or 0) <= hi]
    S = [f(md[i]['Color_Cheek_a_value']) for i in sub if md[i]['baumann_sr'] == 'S']
    R = [f(md[i]['Color_Cheek_a_value']) for i in sub if md[i]['baumann_sr'] == 'R']
    S = [x for x in S if x is not None]; R = [x for x in R if x is not None]
    a, p = auc_p(S, R)
    strat[f'{lo}-{hi}세'] = {'auc': round(a, 4), 'p': round(p, 4),
                            'median_S': round(st.median(S), 3), 'median_R': round(st.median(R), 3),
                            'n_S': len(S), 'n_R': len(R)}
out['redness_by_age_band'] = strat

save_json(out, os.path.join(TAB, '설문라벨_검증.json'))

print(f"설문 라벨: 민감 S {out['n_S']}명 / 저항 R {out['n_R']}명\n")
print(f"{'측정값':22s} {'AUC':>7s} {'구별력':>8s} {'p':>10s} {'q(BH)':>8s}")
for _, label in MEASURES:
    m = out['measures'][label]
    print(f"{label:22s} {m['auc']:7.3f} {m['discrimination']:8.3f} {m['p']:10.2g} {m['q_bh']:8.3g}")
print('\n붉은기 — 나이 층화 후')
for k, v in strat.items():
    print(f"  {k}  S {v['median_S']:.2f} (n={v['n_S']})  R {v['median_R']:.2f} (n={v['n_R']})  AUC {v['auc']:.3f}  p={v['p']:.3f}")
