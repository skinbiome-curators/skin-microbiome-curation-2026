#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2단계 검증의 표본 크기 — AUC 기반 검정력 · 정밀도 (2026-10-03 작성)

팀 합의(2026-10-03): 기대 AUC 0.65 / 0.70 / 0.75 / 0.80 · 검정력 80% · 90% · 피부 증폭 실패 반영.
지금 「군당 65명」의 근거는 손 계산(AUC 신뢰구간 반폭 ±0.10, Hanley–McNeil)뿐이라 코드로 남긴다.

  ① 검정력 기준  H0: AUC = 0.5 를 양측 α=0.05 로 기각할 군당 인원(군 1:1)
                 z(1-α/2)·SE0(n) + z(power)·SE1(n) ≤ |AUC1 − 0.5|  를 만족하는 가장 작은 n
                 SE 는 Hanley & McNeil (1982) 분산식. H0 쪽 분산은 Mann–Whitney 귀무 분산과 같다.
  ② 정밀도 기준  95% 신뢰구간 반폭 = 1.96·SE(AUC, n) — 옛 손 계산(AUC=0.5 가정)을 재현
  ③ 증폭 실패 반영  필요 인원 ÷ 증폭 성공 비율(93/120 — Ahle 2022 *Commun Biol*, 피부 tuf2 앰플리콘)
  ④ 독립 검수  이항정규 모의실험(양성 N(δ,1) · 음성 N(0,1), δ = √2·Φ⁻¹(AUC1)) 으로
                 ①의 n 에서 실제 검정력을 잰다. 검정은 stats_util.auc_p(Mann–Whitney 정규근사).
                 ⚠️ Hanley–McNeil 은 지수형 분포를 가정한 근사라 이항정규 모의실험과 조금 다를 수 있다.

방향: AUC 0.65 는 반대 방향 0.35 와 같은 크기다(*Malassezia* 는 민감군에서 적어 AUC < 0.5).
⚠️ 기대 AUC 는 가정값이다. 42명 *Malassezia* AUC 0.153(=0.847)은 작은 표본이라 부풀려졌을 수 있어 쓰지 않았다.

결과: ../results/검정력_표본크기_261003.json · ../results/검정력_표본크기_261003.tsv (행 = 기대 AUC, matrix 형태)
의존성: Python 3 표준 라이브러리만 (+ ./stats_util.py) · 2026-10-09 최종발표_준비/추가분석/ → code/ 로 옮김
"""
import os, sys, math, time, random, statistics, tracemalloc

try:   # 한글 Windows 콘솔(cp949) 대응 — stats_util.py 와 같은 처리
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

T0 = time.perf_counter()
tracemalloc.start()

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)   # 2026-10-09 code/ 로 옮김 — stats_util.py 가 같은 폴더
from stats_util import save_json, TAB, auc_p   # noqa: E402

ALPHA = 0.05
POWERS = (0.80, 0.90)
AUCS = (0.65, 0.70, 0.75, 0.80)
SUCCESS = 93 / 120          # 피부 시료 증폭 성공 비율 (Ahle 2022)
HALF_TARGET = 0.10          # 옛 손 계산의 목표 반폭
N_NOW = (23, 19)            # 지금 42명 (민감, 비민감)
SEED, N_SIM = 20260918, 2000
Z = statistics.NormalDist()


def hm_var(A, n1, n2):
    """Hanley & McNeil (1982) AUC 분산."""
    q1 = A / (2 - A)
    q2 = 2 * A * A / (1 + A)
    return (A * (1 - A) + (n1 - 1) * (q1 - A * A) + (n2 - 1) * (q2 - A * A)) / (n1 * n2)


def n_power(a1, power):
    """검정력 기준 군당 인원(1:1)."""
    za, zb = Z.inv_cdf(1 - ALPHA / 2), Z.inv_cdf(power)
    d = abs(a1 - 0.5)
    for n in range(2, 5001):
        if za * math.sqrt(hm_var(0.5, n, n)) + zb * math.sqrt(hm_var(a1, n, n)) <= d:
            return n
    return None


def half_width(A, n1, n2):
    return Z.inv_cdf(0.975) * math.sqrt(hm_var(A, n1, n2))


def n_precision(target, A=0.5):
    for n in range(2, 5001):
        if half_width(A, n, n) <= target:
            return n
    return None


def min_detectable_auc(n, power):
    """군당 n 명으로 검정력 power 를 내는 가장 작은 AUC(0.5 위쪽)."""
    for k in range(501, 1000):
        a = k / 1000
        m = n_power(a, power)
        if m is not None and m <= n:
            return a
    return None


def simulate_power(a1, n, rng):
    delta = math.sqrt(2) * Z.inv_cdf(a1)
    hits = 0
    for _ in range(N_SIM):
        pos = [rng.gauss(delta, 1) for _ in range(n)]
        neg = [rng.gauss(0, 1) for _ in range(n)]
        if auc_p(pos, neg)[1] < ALPHA:
            hits += 1
    return hits / N_SIM


rng = random.Random(SEED)
rows = []
print('기대 AUC × 검정력 → 군당 인원 (양측 α=0.05, 1:1)')
print(f"{'AUC':>5} {'검정력':>6} {'군당 n':>7} {'실패 반영':>9} {'총 인원(반영)':>13} {'모의 검정력':>11}")
for a in AUCS:
    for pw in POWERS:
        n = n_power(a, pw)
        n_adj = math.ceil(n / SUCCESS)
        sim = simulate_power(a, n, rng)
        rows.append({'auc': a, 'power': pw, 'n_per_group': n, 'n_per_group_with_failure': n_adj,
                     'total_with_failure': 2 * n_adj, 'simulated_power_at_n': round(sim, 3)})
        print(f'{a:>5.2f} {pw:>6.0%} {n:>7d} {n_adj:>9d} {2 * n_adj:>13d} {sim:>11.3f}')

prec = {'assumed_auc': 0.5,
        'half_width_now_23v19': round(half_width(0.5, *N_NOW), 4),
        'n_per_group_for_half_width_0.10': n_precision(HALF_TARGET),
        'half_width_at_65': round(half_width(0.5, 65, 65), 4)}
mda = {f'{int(pw * 100)}%': min_detectable_auc(65, pw) for pw in POWERS}
mda_adj = {f'{int(pw * 100)}%': min_detectable_auc(math.floor(65 * SUCCESS), pw) for pw in POWERS}

print('\n정밀도 기준(AUC=0.5 가정 — 옛 손 계산 재현)')
print(f"  지금 23 vs 19 반폭 ±{prec['half_width_now_23v19']:.3f} · ±0.10 이 되는 군당 n = "
      f"{prec['n_per_group_for_half_width_0.10']} · 65명 반폭 ±{prec['half_width_at_65']:.3f}")
print(f'군당 65명으로 검출 가능한 최소 AUC: {mda}  (증폭 실패 반영 시 군당 {math.floor(65 * SUCCESS)}명: {mda_adj})')

out = {'settings': {'alpha_two_sided': ALPHA, 'powers': list(POWERS), 'aucs': list(AUCS),
                    'amplification_success': round(SUCCESS, 4),
                    'amplification_source': 'Ahle 2022 Commun Biol PMC9452508 — tuf2 93/120 증폭',
                    'variance': 'Hanley & McNeil 1982',
                    'simulation': f'binormal equal-variance, N_SIM={N_SIM}, seed={SEED}, test=stats_util.auc_p',
                    'note': '기대 AUC 는 가정값. 42명 Malassezia 0.847 은 부풀려졌을 수 있어 쓰지 않음. '
                            'AUC 0.65 는 반대 방향 0.35 와 같은 크기.'},
       'power_table': rows,
       'precision': prec,
       'min_detectable_auc_at_65_per_group': mda,
       'min_detectable_auc_at_65_after_failure': {'n_per_group': math.floor(65 * SUCCESS), **mda_adj}}
save_json(out, os.path.join(TAB, '검정력_표본크기_261003.json'))

# TSV — matrix 형태(행 = 기대 AUC, 열 = 검정력별 값)
tsv = os.path.join(TAB, '검정력_표본크기_261003.tsv')
cols = ['n_per_group', 'n_per_group_with_failure', 'simulated_power_at_n']
with open(tsv, 'w', encoding='utf-8', newline='\n') as fh:
    fh.write('auc\t' + '\t'.join(f'{c}_power{int(p * 100)}' for p in POWERS for c in cols) + '\n')
    for a in AUCS:
        vals = []
        for p in POWERS:
            r = next(x for x in rows if x['auc'] == a and x['power'] == p)
            vals += [str(r[c]) for c in cols]
        fh.write(f'{a}\t' + '\t'.join(vals) + '\n')
print(f'\n저장: results/검정력_표본크기_261003.json · results/검정력_표본크기_261003.tsv')

cur, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()
print(f'메모리(Python 할당 기준) 현재 {cur / 1024 ** 3:.4f} GB · peak {peak / 1024 ** 3:.4f} GB · '
      f'소요 시간 {(time.perf_counter() - T0) / 3600:.6f} Hour')
