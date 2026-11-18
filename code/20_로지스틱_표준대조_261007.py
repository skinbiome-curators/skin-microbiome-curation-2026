#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
① 9번 나이 보정 로지스틱을 표준 패키지와 대조 + 우도비 검정 (2026-10-07 작성 · 팀 합의)

왜: 9번은 고정 200회 경사하강(직접 구현)이라 수렴을 확인하지 않았다. 발표에서
    「나이만 0.653 → 나이+Malassezia 0.883 (+0.23 [0.07–0.41])」를 말하려면 표준 도구로도 같은지 봐야 한다.

무엇을 하나 — 9번과 같은 42명 · 같은 변수(나이 + 후보 비율) · 같은 z-표준화
  (가) 표준 MLE   statsmodels Logit(Newton) 으로 다시 적합 → 계수 · 재대입 AUC · 수렴 여부
  (나) 또 다른 표준 도구  sklearn LogisticRegression(벌점 없음) 계수 — (가)와 같아야 한다
  (다) 우도비 검정  「나이만」 대 「나이 + 후보」 — LR = 2·(llf₁ − llf₀), 자유도 1, χ² 분포
  (라) LOOCV    한 명씩 빼고 MLE 로 맞춰 뺀 사람을 예측(표준화는 학습 41명 안에서) — 9번 LOOCV 와 대조
  (마) 증가분 부트스트랩  9번과 같은 seed · 같은 재표집 순서 4,000회를 MLE 로 — 9번 CI 와 대조
       Newton 이 실패한 재표집(특이행렬 등)은 sklearn(벌점 없음)으로 다시 맞추고, 그래도 실패한 것만 버린다.
       완전분리 · 미수렴은 버리지 않고 세어서 밝힌다(AUC 는 순위만 쓰므로 값이 정해진다).
       (2026-10-07 독립 검수 반영 — 처음 판은 Newton 실패 17건을 버려 증가분이 큰 쪽이 빠졌다)

결과: ../results/로지스틱_표준대조_261007.json · ../results/로지스틱_표준대조_261007.tsv (행 = 모형, matrix 형태)
의존성: numpy · scipy · statsmodels · scikit-learn (conda bio_proj) + ./common42_261007.py · 2026-10-09 최종발표_준비/추가분석/ → code/ 로 옮김
"""
import os, sys, json, math, time, random, warnings, tracemalloc

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

T0 = time.perf_counter()
tracemalloc.start()

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common42_261007 import (load, auc, zscore, logistic_gd, score, save_json, TAB, CAND, FUNGI,   # noqa: E402
                      SEED, mem_time_line)

import numpy as np                                    # noqa: E402
import statsmodels.api as sm                          # noqa: E402
from scipy import stats                               # noqa: E402
from sklearn.linear_model import LogisticRegression   # noqa: E402

N_BOOT = 4000          # 9번과 같음
NINE = json.load(open(os.path.join(TAB, '나이보정_판별력.json'), encoding='utf-8'))


def design(rows, ix, use):
    """9번 fit() 과 같은 표준화 — 선택된 사람(ix) 안에서 열마다 z-점수."""
    z = [zscore([rows[i][c] for i in ix]) for c in use]
    return [[z[j][t] for j in range(len(use))] for t in range(len(ix))]


def mle(X, y):
    """statsmodels Logit(Newton) 으로 맞춘다. Newton 이 실패하면(특이행렬 등) sklearn(벌점 없음)으로 다시 맞춘다.
    (독립 검수 2026-10-07: 진균 부트스트랩에서 Newton 실패 17건을 버리면 증가분이 큰 쪽이 빠져 CI 위쪽이 낮아진다.)
    완전분리 경고가 나거나 수렴하지 않으면 converged=False 로 표시한다 — AUC 는 순위만 쓰므로 값은 정해진다.
    반환: dict(params, llf, converged, method, res) — 둘 다 실패하면 None."""
    Xa = np.asarray(X, dtype=float)
    Xc = sm.add_constant(Xa, has_constant='add')
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter('always')
        try:
            res = sm.Logit(np.asarray(y, dtype=float), Xc).fit(disp=0, method='newton', maxiter=200)
            if np.all(np.isfinite(res.params)):
                sep = any('separation' in str(w.message).lower() for w in wl)
                conv = bool(res.mle_retvals.get('converged', False)) and not sep
                return {'params': res.params, 'llf': float(res.llf), 'converged': conv, 'method': 'newton', 'res': res}
        except Exception:
            pass
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            skl = LogisticRegression(penalty=None, solver='lbfgs', max_iter=100000).fit(Xa, np.asarray(y, dtype=int))
        except Exception:
            return None
    params = np.concatenate([skl.intercept_, skl.coef_[0]])
    if not np.all(np.isfinite(params)):
        return None
    return {'params': params, 'llf': float('nan'), 'converged': False, 'method': 'sklearn_fallback', 'res': None}


def lin(X, params):
    Xc = sm.add_constant(np.asarray(X, dtype=float), has_constant='add')
    return list(Xc @ params)


def loocv_mle(rows, use):
    """반환: (LOOCV AUC, Newton 이 아닌 적합 수)."""
    n = len(rows)
    y = [r['y'] for r in rows]
    oof = [0.0] * n
    n_alt = 0
    for h in range(n):
        tr = [i for i in range(n) if i != h]
        mu, sd = [], []
        for c in use:
            v = [rows[i][c] for i in tr]
            m = sum(v) / len(v)
            s = math.sqrt(sum((x - m) ** 2 for x in v) / max(1, len(v) - 1)) or 1.0
            mu.append(m); sd.append(s)
        X = [[(rows[i][c] - mu[j]) / sd[j] for j, c in enumerate(use)] for i in tr]
        fit = mle(X, [y[i] for i in tr])
        if fit is None:
            return float('nan'), n_alt
        n_alt += fit['method'] != 'newton' or not fit['converged']
        xh = [[(rows[h][c] - mu[j]) / sd[j] for j, c in enumerate(use)]]
        oof[h] = lin(xh, fit['params'])[0]
    return auc(y, oof), n_alt


def run(label, rows, cols, nine):
    n = len(rows)
    idx = list(range(n))
    y = [r['y'] for r in rows]
    out = {'n': n, 'n_sensitive': sum(y), 'models': {}}

    def one(use):
        X = design(rows, idx, use)
        w_gd, b_gd = logistic_gd(X, y)
        a_gd = auc(y, score(X, w_gd, b_gd))
        fit = mle(X, y)
        if fit is None or fit['method'] != 'newton':
            raise RuntimeError(f'전체 {len(y)}명 MLE(Newton) 적합 실패: {use}')
        params, llf, conv, res = fit['params'], fit['llf'], fit['converged'], fit['res']
        a_mle = auc(y, lin(X, params))
        loo_auc, loo_alt = loocv_mle(rows, list(use))
        skl = LogisticRegression(penalty=None, solver='lbfgs', max_iter=100000, tol=1e-12).fit(np.asarray(X), y)
        coef_skl = [float(skl.intercept_[0])] + [float(c) for c in skl.coef_[0]]
        return {'use': list(use),
                'auc_gd_9번방식': round(a_gd, 4), 'auc_mle': round(a_mle, 4),
                'coef_gd_9번방식': [round(b_gd, 4)] + [round(w, 4) for w in w_gd],
                'coef_mle_statsmodels': [round(float(p), 4) for p in params],
                'coef_sklearn': [round(c, 4) for c in coef_skl],
                'max_abs_diff_statsmodels_vs_sklearn': round(max(abs(float(p) - c) for p, c in zip(params, coef_skl)), 6),
                'converged': conv, 'llf': round(llf, 4),
                'wald_p_last': float(res.pvalues[-1]),
                'loocv_auc_mle': round(loo_auc, 4),
                'loocv_fits_not_newton_or_not_converged': int(loo_alt)}, llf

    base, llf0 = one(['age'])
    base['auc_9번_json'] = nine['auc_age_only']
    base['loocv_9번_json'] = nine['auc_age_only_loocv']
    out['models']['age'] = base
    print(f'\n== {label} — n={n} ==')
    print(f"  {'모형':<26}{'9번(GD)':>9}{'MLE':>8}{'수렴':>6}{'GD 계수(후보)':>14}{'MLE 계수(후보)':>15}"
          f"{'LOOCV 9번':>11}{'LOOCV MLE':>11}")
    print(f"  {'나이만':<26}{base['auc_gd_9번방식']:>9.3f}{base['auc_mle']:>8.3f}{str(base['converged']):>6}"
          f"{'':>14}{'':>15}{nine['auc_age_only_loocv']:>11.3f}{base['loocv_auc_mle']:>11.3f}")

    for c in cols:
        m, llf1 = one(['age', c])
        lr = 2.0 * (llf1 - llf0)
        m['lr_test_vs_age'] = {'LR': round(lr, 4), 'df': 1, 'p': float(stats.chi2.sf(lr, 1))}
        nc = nine['candidates'][c]
        m['auc_9번_json'] = nc['auc_age_plus']
        m['delta_9번_json'] = nc['delta']
        m['ci_9번_json'] = nc['ci']
        m['loocv_9번_json'] = nc['auc_age_plus_loocv']
        m['delta_loocv_9번_json'] = nc['delta_loocv']
        m['delta_mle'] = round(m['auc_mle'] - base['auc_mle'], 4)
        m['delta_loocv_mle'] = round(m['loocv_auc_mle'] - base['loocv_auc_mle'], 4)
        out['models'][c] = m
        print(f"  {'나이 + ' + c:<26}{m['auc_gd_9번방식']:>9.3f}{m['auc_mle']:>8.3f}{str(m['converged']):>6}"
              f"{m['coef_gd_9번방식'][-1]:>+14.3f}{m['coef_mle_statsmodels'][-1]:>+15.3f}"
              f"{nc['auc_age_plus_loocv']:>11.3f}{m['loocv_auc_mle']:>11.3f}"
              f"   LR={lr:.2f} p={m['lr_test_vs_age']['p']:.4f}")

    # (마) 부트스트랩 — 9번과 같은 seed · 같은 재표집 순서
    rng = random.Random(SEED)
    deltas = {c: [] for c in cols}
    fails = {c: 0 for c in cols}       # Newton · sklearn 모두 실패 → 버림
    alt = {c: 0 for c in cols}         # Newton 실패 → sklearn 으로 다시 맞춤
    notconv = {c: 0 for c in cols}     # 완전분리 경고 · 미수렴(다시 맞춘 것 포함) — 버리지 않고 센다
    used = 0
    for _ in range(N_BOOT):
        ix = [rng.randrange(n) for _ in range(n)]
        yy = [rows[i]['y'] for i in ix]
        if len(set(yy)) < 2:
            continue
        used += 1
        X0 = design(rows, ix, ['age'])
        f0 = mle(X0, yy)
        for c in cols:
            X1 = design(rows, ix, ['age', c])
            f1 = mle(X1, yy)
            if f0 is None or f1 is None:
                fails[c] += 1
                continue
            alt[c] += (f0['method'] != 'newton') or (f1['method'] != 'newton')
            notconv[c] += (not f0['converged']) or (not f1['converged'])
            deltas[c].append(auc(yy, lin(X1, f1['params'])) - auc(yy, lin(X0, f0['params'])))
    for c in cols:
        s = sorted(deltas[c])
        lo, hi = s[int(0.025 * len(s))], s[int(0.975 * len(s))]
        out['models'][c]['boot_mle'] = {'n_resamples_used': used, 'n_mle_failed': fails[c],
                                        'n_refit_with_sklearn': alt[c], 'n_not_converged_or_separation': notconv[c],
                                        'ci_delta': [round(lo, 4), round(hi, 4)],
                                        'crosses_zero': lo <= 0 <= hi}
        print(f"  부트스트랩 {c:<16} MLE 증가분 CI [{lo:+.3f}, {hi:+.3f}]  "
              f"(9번 GD CI {out['models'][c]['ci_9번_json']}) · 버림 {fails[c]} · sklearn 재적합 {alt[c]} · "
              f"미수렴/분리 {notconv[c]} / {used}")
    return out


print('=' * 96)
print('① 9번 나이 보정 로지스틱 — 표준 MLE(statsmodels · sklearn) 대조 + 우도비 검정 · 03 민감피부 42명')
print('=' * 96)

out = {'seed': SEED, 'n_boot': N_BOOT,
       'note': '9번(경사하강 200회)과 같은 입력 · 표준화로 표준 MLE 를 다시 적합했다. 재대입 AUC 는 낙관 편향이 있어 '
               'LOOCV 를 함께 본다. 우도비 검정은 「나이만」 대비 후보를 더했을 때의 적합도 개선이다(n=42).',
       'bacteria': run('세균 (16S V4-V5)', load('16S', CAND), CAND, NINE['bacteria']),
       'fungi': run('진균 (ITS1)', load('ITS', FUNGI), FUNGI, NINE['fungi'])}

save_json(out, os.path.join(TAB, '로지스틱_표준대조_261007.json'))

tsv = os.path.join(TAB, '로지스틱_표준대조_261007.tsv')
cols_t = ['auc_9번_json', 'auc_gd_9번방식', 'auc_mle', 'converged', 'delta_9번_json', 'delta_mle',
          'loocv_9번_json', 'loocv_auc_mle', 'delta_loocv_9번_json', 'delta_loocv_mle',
          'lr', 'lr_p', 'wald_p_last', 'boot_ci_lo_mle', 'boot_ci_hi_mle', 'boot_mle_failed',
          'boot_refit_sklearn', 'boot_not_converged_or_separation']
with open(tsv, 'w', encoding='utf-8', newline='\n') as fh:
    fh.write('model\t' + '\t'.join(cols_t) + '\n')
    for grp in ('bacteria', 'fungi'):
        for name, m in out[grp]['models'].items():
            row = dict(m)
            row['lr'] = m.get('lr_test_vs_age', {}).get('LR', '')
            row['lr_p'] = m.get('lr_test_vs_age', {}).get('p', '')
            row['boot_ci_lo_mle'] = m.get('boot_mle', {}).get('ci_delta', ['', ''])[0]
            row['boot_ci_hi_mle'] = m.get('boot_mle', {}).get('ci_delta', ['', ''])[1]
            row['boot_mle_failed'] = m.get('boot_mle', {}).get('n_mle_failed', '')
            row['boot_refit_sklearn'] = m.get('boot_mle', {}).get('n_refit_with_sklearn', '')
            row['boot_not_converged_or_separation'] = m.get('boot_mle', {}).get('n_not_converged_or_separation', '')
            label = f'{grp}:age' if name == 'age' else f'{grp}:age+{name}'
            fh.write(label + '\t' + '\t'.join(str(row.get(c, '')) for c in cols_t) + '\n')

print('\n저장: results/로지스틱_표준대조_261007.json · results/로지스틱_표준대조_261007.tsv')
print(mem_time_line(T0))
tracemalloc.stop()
