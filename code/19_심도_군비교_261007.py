#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
민감군 vs 비민감군 — 시퀀싱 심도 · 오염 의심 진균 비교 (2026-10-03 작성 · 사용자 합의)

왜: Rhodotorula 는 우리 자료에서 오염 의심(심도와 음의 상관)인데, Lu 2024(2bRAD-M)는 민감 안면에서
    증가를 보고했다. 민감군 시료의 심도가 낮다면, 오염처럼 행동하는 균이 민감군에서 늘어 보이는
    착시일 수 있다. 질문 대비용 점검이며 **판정 등급은 바꾸지 않는다.**

  ① 42명 16S · ITS 의 read_count(원 read) · read_depth_final(처리 후) — 두 군 비교
  ② 42명 ITS 의 g__Rhodotorula · g__Naganishia 비율 — 두 군 비교 (g__Malassezia 는 기준으로 함께 표시, BH 제외)
  검정: stats_util.auc_p (Mann–Whitney, 정규근사) · AUC = 민감군이 클 확률 · ①②의 6개 검정에 BH

⚠️ read 수는 생물량(DNA 양)이 아니다 — 원자료에 DNA 농도가 없어 **대리 지표**로만 쓴다.
⚠️ 결과를 본 뒤 고른 점검이 아니라 사전에 합의한 점검이지만, 판정 이후의 **사후 점검**이다.

군 판정 규칙: 9_나이보정_판별력.py 와 같다('Non' 이 있으면 비민감, 'ensitive' 이면 민감).
입력: ../data/03_sensitive_{16S,ITS}_sample_info.tsv · 03_sensitive_ITS_genus_pct.tsv
결과: ../results/심도_군비교_261003.json
의존성: Python 3 표준 라이브러리만 (+ ./stats_util.py) · 2026-10-09 최종발표_준비/추가분석/ → code/ 로 옮김
"""
import os, sys, csv, time, statistics, tracemalloc

try:   # 한글 Windows 콘솔(cp949) 대응
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

T0 = time.perf_counter()
tracemalloc.start()

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)   # 2026-10-09 code/ 로 옮김 — stats_util.py 가 같은 폴더
from stats_util import save_json, TAB, RAW, auc_p, bh, f   # noqa: E402

DEPTH_COLS = ('read_count', 'read_depth_final')
TAXA = ('g__Rhodotorula', 'g__Naganishia')
REF = 'g__Malassezia'


def read_tsv(path):
    with open(path, encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def label(grp):
    g = str(grp)
    if 'Non' in g or 'non' in g:
        return 0
    return 1 if 'ensitive' in g else None


def compare(s, n):
    a, p = auc_p(s, n)
    return {'n_sensitive': len(s), 'n_non_sensitive': len(n),
            'median_sensitive': round(statistics.median(s), 4),
            'median_non_sensitive': round(statistics.median(n), 4),
            'auc_sensitive_higher': round(a, 3), 'p_mann_whitney_normal_approx': round(p, 4)}


tests, out = [], {'depth': {}, 'taxa_its': {}, 'reference_its': {}}

for assay in ('16S', 'ITS'):
    info = read_tsv(os.path.join(RAW, f'03_sensitive_{assay}_sample_info.tsv'))
    for col in DEPTH_COLS:
        s = [f(r[col]) for r in info if label(r['group']) == 1 and r.get(col)]
        n = [f(r[col]) for r in info if label(r['group']) == 0 and r.get(col)]
        res = compare(s, n)
        out['depth'][f'{assay}_{col}'] = res
        tests.append(('depth', f'{assay}_{col}'))

info_its = {r['run']: label(r['group']) for r in read_tsv(os.path.join(RAW, '03_sensitive_ITS_sample_info.tsv'))}
pct = read_tsv(os.path.join(RAW, '03_sensitive_ITS_genus_pct.tsv'))
for g in TAXA + (REF,):
    s = [f(r[g]) for r in pct if info_its.get(r['run']) == 1]
    n = [f(r[g]) for r in pct if info_its.get(r['run']) == 0]
    res = compare(s, n)
    res['detect_pct_sensitive'] = round(100 * sum(v > 0 for v in s) / len(s), 1)
    res['detect_pct_non_sensitive'] = round(100 * sum(v > 0 for v in n) / len(n), 1)
    if g == REF:
        out['reference_its'][g] = res
    else:
        out['taxa_its'][g] = res
        tests.append(('taxa_its', g))

q = bh([out[k][name]['p_mann_whitney_normal_approx'] for k, name in tests])
for (k, name), qq in zip(tests, q):
    out[k][name]['q_bh_6tests'] = round(qq, 4)

out['note'] = ('read 수는 생물량이 아니다(DNA 농도 없음 · 대리 지표). 판정 등급은 바꾸지 않는다. '
               'AUC = 민감군이 클 확률. BH 는 심도 4 + 진균 2 = 6개 검정. Malassezia 는 기준으로만 표시.')
save_json(out, os.path.join(TAB, '심도_군비교_261003.json'))

print(f"{'항목':<28}{'민감 중앙':>12}{'비민감 중앙':>12}{'AUC(민감>)':>11}{'p':>8}{'q':>8}")
for k in ('depth', 'taxa_its', 'reference_its'):
    for name, r in out[k].items():
        qv = r.get('q_bh_6tests', float('nan'))
        print(f"{name:<28}{r['median_sensitive']:>12}{r['median_non_sensitive']:>12}"
              f"{r['auc_sensitive_higher']:>11.3f}{r['p_mann_whitney_normal_approx']:>8.4f}{qv:>8.4f}")
print('\n저장: results/심도_군비교_261003.json')

cur, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()
print(f'메모리(Python 할당 기준) 현재 {cur / 1024 ** 3:.4f} GB · peak {peak / 1024 ** 3:.4f} GB · '
      f'소요 시간 {(time.perf_counter() - T0) / 3600:.6f} Hour')
