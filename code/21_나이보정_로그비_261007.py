#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
② 나이 보정 모형에 균 「비율」 대신 「로그비」를 넣으면 증가분이 달라지는가 (2026-10-07 작성 · 사용자 합의)

왜: 비율은 사람마다 합이 100% 로 묶여 있어 「비율을 그대로 회귀에 넣어도 되나」 질문이 나온다(발표문구 §1-4).
    9번은 비율을 넣었다. 같은 계산을 로그비로 다시 해서 결론(증가분 · CI · LOOCV)이 유지되는지 본다.

방법 — 9번 run() 과 같다. 바꾼 것은 후보 열 하나뿐이다(비율 → 로그비).
  (가) 나이만  (나) 나이 + 후보 로그비  (다) 증가분 = (나) − (가)
  로지스틱은 9번과 같은 경사하강 200회(common42_261007.logistic_gd) · z-표준화 · 재대입 AUC
  부트스트랩 4,000회(9번과 같은 seed · 같은 재표집 순서) · LOOCV(표준화는 학습 폴드 안에서)
  로그비 정의: 세균 = 4번 ④(분모 = 후보·구강 의심·미분류를 뺀 나머지, EPS 1e-5) · 진균 = 8번 ②(분모 = 그 속·미분류를 뺀 나머지, EPS 1e-3)
  자체 점검: 로그비 「단독」 AUC 가 4번 · 8번 저장값과 같아야 같은 정의다.

결과: ../results/나이보정_로그비_261007.json · ../results/나이보정_로그비_261007.tsv (행 = 후보, matrix 형태)
의존성: Python 3 표준 라이브러리만 (+ ./common42_261007.py · ./stats_util.py) · 2026-10-09 최종발표_준비/추가분석/ → code/ 로 옮김
"""
import os, sys, json, math, time, random, tracemalloc

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

T0 = time.perf_counter()
tracemalloc.start()

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common42_261007 import (load, add_logratio, auc, zscore, logistic_gd, score, save_json, TAB, CAND, FUNGI,   # noqa: E402
                      SEED, mem_time_line)

N_BOOT = 4000
NINE = json.load(open(os.path.join(TAB, '나이보정_판별력.json'), encoding='utf-8'))
FOUR = json.load(open(os.path.join(TAB, '조성효과_점검.json'), encoding='utf-8'))['gate4_logratio']
EIGHT = json.load(open(os.path.join(TAB, '진균_판별력.json'), encoding='utf-8'))['candidates']


def loocv(rows, use):
    """9번 loocv() 와 같다."""
    n = len(rows)
    y = [r['y'] for r in rows]
    oof = [0.0] * n
    for h in range(n):
        tr = [i for i in range(n) if i != h]
        mu, sd = [], []
        for c in use:
            v = [rows[i][c] for i in tr]
            m = sum(v) / len(v)
            sdev = math.sqrt(sum((x - m) ** 2 for x in v) / max(1, len(v) - 1)) or 1.0
            mu.append(m); sd.append(sdev)
        X = [[(rows[i][c] - mu[j]) / sd[j] for j, c in enumerate(use)] for i in tr]
        w, b = logistic_gd(X, [y[i] for i in tr])
        xh = [(rows[h][c] - mu[j]) / sd[j] for j, c in enumerate(use)]
        oof[h] = b + sum(w[j] * xh[j] for j in range(len(w)))
    return auc(y, oof)


TOL = 0.01


def sign_tol(d):
    """LOOCV 증가분의 부호 — |d| < 0.01 이면 0(사실상 변화 없음)."""
    return 0 if abs(d) < TOL else (1 if d > 0 else -1)


def run(label, prefix, cols, nine):
    rows = add_logratio(prefix, load(prefix, cols), cols)
    n = len(rows)
    idx = list(range(n))
    y = [r['y'] for r in rows]

    def fit(ix, use):
        yy = [rows[i]['y'] for i in ix]
        if len(set(yy)) < 2:
            return float('nan')
        z = [zscore([rows[i][c] for i in ix]) for c in use]
        X = [[z[j][t] for j in range(len(use))] for t in range(len(ix))]
        w, b = logistic_gd(X, yy)
        return auc(yy, score(X, w, b))

    base = fit(idx, ['age'])
    base_cv = loocv(rows, ['age'])
    rng = random.Random(SEED)
    boots = {c: [] for c in cols}
    for _ in range(N_BOOT):
        ix = [rng.randrange(n) for _ in range(n)]
        if len(set(rows[i]['y'] for i in ix)) < 2:
            continue
        bb = fit(ix, ['age'])
        for c in cols:
            boots[c].append(fit(ix, ['age', 'lr_' + c]) - bb)

    out = {'n': n, 'n_sensitive': sum(y), 'auc_age_only': round(base, 3), 'auc_age_only_loocv': round(base_cv, 3),
           'auc_age_only_9번_json': nine['auc_age_only'], 'candidates': {}}
    print(f'\n== {label} — n={n} · 나이만 AUC {base:.3f} (9번 {nine["auc_age_only"]}) · 교차검증 {base_cv:.3f} ==')
    print(f"  {'후보':<18}{'로그비 단독':>11}{'저장값':>8}{'나이+로그비':>12}{'증가분':>9}  증가분 95% CI{'':>8}"
          f"{'LOOCV 증가분':>13}   │ 9번(비율) 증가분 · CI · LOOCV")
    for c in cols:
        lr_alone = auc(y, [r['lr_' + c] for r in rows])
        stored = FOUR[c]['auc'] if prefix == '16S' else EIGHT[c]['gate2_composition']['auc_logratio']
        full = fit(idx, ['age', 'lr_' + c])
        d = full - base
        s = sorted(v for v in boots[c] if v == v)
        lo, hi = s[int(0.025 * len(s))], s[int(0.975 * len(s))]
        full_cv = loocv(rows, ['age', 'lr_' + c])
        d_cv = full_cv - base_cv
        nc = nine['candidates'][c]
        sg_lr, sg_9 = sign_tol(d_cv), sign_tol(nc['delta_loocv'])
        out['candidates'][c] = {
            'auc_logratio_alone': round(lr_alone, 3), 'auc_logratio_alone_stored_4or8': stored,
            'definition_matches_4or8': round(lr_alone, 3) == stored,
            'auc_age_plus_logratio': round(full, 3), 'delta_logratio': round(d, 3),
            'ci_logratio': [round(lo, 3), round(hi, 3)], 'crosses_zero_logratio': lo <= 0 <= hi,
            'auc_age_plus_logratio_loocv': round(full_cv, 3), 'delta_logratio_loocv': round(d_cv, 3),
            'ratio_9번': {'auc_age_plus': nc['auc_age_plus'], 'delta': nc['delta'], 'ci': nc['ci'],
                         'crosses_zero': nc['crosses_zero'], 'delta_loocv': nc['delta_loocv']},
            # 2026-10-07 독립 검수 반영 — CI 결론과 LOOCV 부호를 따로 본다. 부호는 |증가분| < 0.01 이면 0 으로(경계값 흔들림 방지)
            'ci_conclusion_same_as_9번': (lo <= 0 <= hi) == nc['crosses_zero'],
            'loocv_sign_logratio': sg_lr, 'loocv_sign_9번': sg_9,
            'loocv_sign_same_as_9번': sg_lr == sg_9}
        m = out['candidates'][c]
        verdict = ('CI·LOOCV 모두 같음' if m['ci_conclusion_same_as_9번'] and m['loocv_sign_same_as_9번']
                   else ('⚠️ CI 결론 다름' if not m['ci_conclusion_same_as_9번'] else '⚠️ LOOCV 부호 다름'))
        print(f"  {c:<18}{lr_alone:>11.3f}{stored:>8.3f}{full:>12.3f}{d:>+9.3f}  [{lo:+.3f}, {hi:+.3f}]"
              f"{d_cv:>+13.3f}   │ {nc['delta']:+.3f} {nc['ci']} {nc['delta_loocv']:+.3f}  {verdict}")
    return out


print('=' * 100)
print('② 나이 보정 — 균 비율 대신 로그비 · 03 민감피부 42명 (9번과 같은 방법)')
print('=' * 100)

out = {'seed': SEED, 'n_boot': N_BOOT,
       'note': '9번과 같은 계산에서 후보 열만 비율 → 로그비로 바꿨다. 증가분 CI 가 0 을 포함하는지와 LOOCV 증가분의 '
               '부호(|증가분|<0.01 이면 0)가 9번(비율)과 같으면 「비율을 그대로 넣었기 때문」이라는 반론은 결론을 바꾸지 못한다. '
               '둘을 따로 저장한다.',
       'logratio_definition': {'bacteria': '4번 ④ 와 같음 — log10((후보+1e-5)/(후보·구강 의심·미분류를 뺀 나머지 합+1e-5))',
                               'fungi': '8번 ② 와 같음 — log10((속+1e-3)/(그 속·미분류를 뺀 나머지 합+1e-3))'},
       'bacteria': run('세균 (16S V4-V5)', '16S', CAND, NINE['bacteria']),
       'fungi': run('진균 (ITS1)', 'ITS', FUNGI, NINE['fungi'])}

save_json(out, os.path.join(TAB, '나이보정_로그비_261007.json'))

tsv = os.path.join(TAB, '나이보정_로그비_261007.tsv')
cols_t = ['auc_logratio_alone', 'auc_logratio_alone_stored_4or8', 'auc_age_plus_logratio', 'delta_logratio',
          'ci_lo_logratio', 'ci_hi_logratio', 'delta_logratio_loocv',
          'ratio_9번_delta', 'ratio_9번_ci_lo', 'ratio_9번_ci_hi', 'ratio_9번_delta_loocv',
          'ci_conclusion_same_as_9번', 'loocv_sign_logratio', 'loocv_sign_9번', 'loocv_sign_same_as_9번']
with open(tsv, 'w', encoding='utf-8', newline='\n') as fh:
    fh.write('candidate\t' + '\t'.join(cols_t) + '\n')
    for grp in ('bacteria', 'fungi'):
        for c, m in out[grp]['candidates'].items():
            row = dict(m)
            row['ci_lo_logratio'], row['ci_hi_logratio'] = m['ci_logratio']
            row['ratio_9번_delta'] = m['ratio_9번']['delta']
            row['ratio_9번_ci_lo'], row['ratio_9번_ci_hi'] = m['ratio_9번']['ci']
            row['ratio_9번_delta_loocv'] = m['ratio_9번']['delta_loocv']
            fh.write(c + '\t' + '\t'.join(str(row.get(k, '')) for k in cols_t) + '\n')

print('\n저장: results/나이보정_로그비_261007.json · results/나이보정_로그비_261007.tsv')
print(mem_time_line(T0))
tracemalloc.stop()
