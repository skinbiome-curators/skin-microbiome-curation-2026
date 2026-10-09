#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
④ 정규성 점검 — ④ 판별력을 Mann–Whitney(순위 기반)로 잰 근거 (2026-10-07 작성 · 사용자 합의)

왜: 지금까지 코드에 정규성 검정이 없어 「정규성을 확인했다」고 말할 수 없었다(발표 금칙).
    「균 비율은 0이 많고 한쪽으로 쏠려 평균 비교(t 검정)보다 순위 비교가 맞다」를 수치로 남긴다.

방법 — 9번과 같은 42명(민감 23 · 비민감 19) · 후보 4속의 비율(%) · 군별로
  - Shapiro–Wilk 검정(scipy.stats.shapiro) : p < 0.05 면 「정규분포라는 가정을 기각」
  - 0 의 비율 · 왜도(scipy.stats.skew, 표본 보정) · 중앙값 · 평균
  - 참고로 나이도 같은 방식으로 본다(로지스틱의 입력)
  - 여러 검정이므로 BH 보정 q 값도 함께 적는다
⚠️ 정규성을 기각하지 못했다고 정규분포라는 뜻은 아니다(n=19–23 이라 검정력이 낮다).

결과: ../results/정규성_점검_261007.json · ../results/정규성_점검_261007.tsv (행 = 변수 × 군, matrix 형태)
의존성: scipy (conda bio_proj) + ./common42_261007.py · ./stats_util.py · 2026-10-09 최종발표_준비/추가분석/ → code/ 로 옮김
"""
import os, sys, time, statistics, tracemalloc

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

T0 = time.perf_counter()
tracemalloc.start()

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common42_261007 import load, save_json, TAB, CAND, FUNGI, mem_time_line   # noqa: E402
from stats_util import bh                                              # noqa: E402

from scipy import stats                                                # noqa: E402

GROUPS = ((1, '민감'), (0, '비민감'))


def check(name, x):
    w, p = stats.shapiro(x)
    return {'variable': name, 'n': len(x), 'shapiro_W': round(float(w), 4), 'shapiro_p': float(p),
            'zero_n': sum(1 for v in x if v == 0), 'zero_pct': round(100 * sum(1 for v in x if v == 0) / len(x), 1),
            'skewness': round(float(stats.skew(x, bias=False)), 3),
            'median': round(statistics.median(x), 3), 'mean': round(statistics.mean(x), 3)}


rows_all = []
for prefix, cols in (('16S', CAND), ('ITS', FUNGI)):
    rows = load(prefix, cols)
    for c in cols:
        for g, gname in GROUPS:
            r = check(c, [r[c] for r in rows if r['y'] == g])
            r['group'] = gname
            r['kind'] = 'candidate'
            rows_all.append(r)
    if prefix == '16S':   # 나이는 한 번만(16S·ITS 같은 42명)
        for g, gname in GROUPS:
            r = check('age', [r['age'] for r in rows if r['y'] == g])
            r['group'] = gname
            r['kind'] = 'reference'
            rows_all.append(r)

cand_rows = [r for r in rows_all if r['kind'] == 'candidate']
for r, q in zip(cand_rows, bh([r['shapiro_p'] for r in cand_rows])):
    r['shapiro_q_bh'] = round(q, 4)

n_reject = sum(1 for r in cand_rows if r['shapiro_p'] < 0.05)
summary = (f'후보 4속 × 2군 = {len(cand_rows)}개 분포 중 {n_reject}개에서 정규분포 가정을 기각(Shapiro–Wilk p<0.05). '
           f'→ 평균 비교(t 검정) 대신 순위 비교(Mann–Whitney)를 쓴 근거.')

print('=' * 96)
print('④ 정규성 점검 — 후보 비율(%) 군별 Shapiro–Wilk · 03 민감피부 42명')
print('=' * 96)
print(f"  {'변수':<18}{'군':<6}{'n':>4}{'W':>8}{'p':>11}{'q(BH)':>9}{'0 비율':>8}{'왜도':>8}{'중앙값':>9}{'평균':>9}")
for r in rows_all:
    q = r.get('shapiro_q_bh', '')
    print(f"  {r['variable']:<18}{r['group']:<6}{r['n']:>4}{r['shapiro_W']:>8.3f}{r['shapiro_p']:>11.2e}"
          f"{(f'{q:.4f}' if q != '' else '—'):>9}{r['zero_pct']:>7.1f}%{r['skewness']:>8.2f}{r['median']:>9.3f}{r['mean']:>9.3f}")
print(f'\n  {summary}')
print('  ⚠️ 기각하지 못한 분포가 정규분포라는 뜻은 아니다(n=19–23).')

out = {'method': 'scipy.stats.shapiro(군별) · 0 비율 · 표본 보정 왜도. q 는 후보 8개 분포에 대한 BH 보정.',
       'summary': summary,
       'caution': '정규성을 기각하지 못했다고 정규분포인 것은 아니다(n=19–23 이라 검정력이 낮다). 나이는 참고(로지스틱 입력).',
       'rows': rows_all}
save_json(out, os.path.join(TAB, '정규성_점검_261007.json'))

tsv = os.path.join(TAB, '정규성_점검_261007.tsv')
cols_t = ['n', 'shapiro_W', 'shapiro_p', 'shapiro_q_bh', 'zero_n', 'zero_pct', 'skewness', 'median', 'mean']
with open(tsv, 'w', encoding='utf-8', newline='\n') as fh:
    fh.write('variable_group\t' + '\t'.join(cols_t) + '\n')
    for r in rows_all:
        fh.write(f"{r['variable']}:{r['group']}\t" + '\t'.join(str(r.get(c, '')) for c in cols_t) + '\n')

print('\n저장: results/정규성_점검_261007.json · results/정규성_점검_261007.tsv')
print(mem_time_line(T0))
tracemalloc.stop()
