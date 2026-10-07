#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
나이를 보정하면 균이 무엇을 더하는가 — 발표 질문 대비

④ 관문은 균 하나만 써서 AUC 를 냈다. 그런데 42명에서 민감군은 29세,
비민감군은 40세이고 **나이 자체의 AUC 0.347 이 세균 세 속을 전부 이긴다.**
그래서 발표에서 이 질문이 반드시 나온다.

  "그럼 나이만 물어보면 되지 왜 균을 재나?"

답하려면 균 단독 AUC 가 아니라 **나이 위에 얼마나 더하는지**를 재야 한다.
  (가) 나이만          로지스틱
  (나) 나이 + 후보     로지스틱
  (다) 증가분          (나) - (가)

n=42 라 검정력이 거의 없다. 이 계산의 목적은 "얼마나 모르는가" 를 재는 것이다.
부트스트랩은 사람 단위, seed 20260918 (다른 스크립트와 동일).

결과: results/나이보정_판별력.json
"""
import os, sys, csv, math, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_util import save_json, BASE, RAW, TAB, CAND, f

SEED, N_BOOT = 20260918, 4000
FUNGI = ['g__Malassezia']


def read_tsv(path):
    with open(path, encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def auc(y, x):
    pos = [x[i] for i in range(len(y)) if y[i] == 1]
    neg = [x[i] for i in range(len(y)) if y[i] == 0]
    if not pos or not neg:
        return float('nan')
    s = 0.0
    for p in pos:
        for n in neg:
            s += 1.0 if p > n else (0.5 if p == n else 0.0)
    return s / (len(pos) * len(neg))


def zscore(col):
    m = sum(col) / len(col)
    v = sum((c - m) ** 2 for c in col) / max(1, len(col) - 1)
    return [(c - m) / (math.sqrt(v) if v > 0 else 1.0) for c in col]


def logistic(X, y, iters=200, lr=0.1):
    n, k = len(X), len(X[0])
    w, b = [0.0] * k, 0.0
    for _ in range(iters):
        gw, gb = [0.0] * k, 0.0
        for i in range(n):
            z = max(-30.0, min(30.0, b + sum(w[j] * X[i][j] for j in range(k))))
            e = 1.0 / (1.0 + math.exp(-z)) - y[i]
            gb += e
            for j in range(k):
                gw[j] += e * X[i][j]
        b -= lr * gb / n
        for j in range(k):
            w[j] -= lr * gw[j] / n
    return w, b


def score(X, w, b):
    return [max(-30.0, min(30.0, b + sum(w[j] * r[j] for j in range(len(w))))) for r in X]


def load(prefix, cols):
    """03_sensitive 의 16S 또는 ITS 를 읽어 (행: 사람) 로 정리."""
    info = read_tsv(os.path.join(RAW, f'03_sensitive_{prefix}_sample_info.tsv'))
    pct = read_tsv(os.path.join(RAW, f'03_sensitive_{prefix}_genus_pct.tsv'))
    by = {r['run']: r for r in pct}
    rows = []
    for r in info:
        run, grp, age = r.get('run'), str(r.get('group', '')), r.get('age')
        if run not in by or not age:
            continue
        y = 0 if 'Non' in grp or 'non' in grp else (1 if 'ensitive' in grp else None)
        if y is None:
            continue
        vals = {}
        ok = True
        for c in cols:
            if c not in by[run]:
                ok = False; break
            vals[c] = f(by[run][c])
        if not ok:
            continue
        rows.append(dict(run=run, y=y, age=f(age), **vals))
    return rows


def loocv(rows, use):
    """한 사람을 빼고 41명으로 맞춘 뒤 뺀 사람을 예측한다. 42명 전원에게 한 번씩.

    fit() 은 같은 사람으로 맞추고 같은 사람으로 재는 **재대입(resubstitution)** 이라
    변수를 하나 더하면 AUC 가 거의 항상 오른다 — 낙관 편향이다. 그 편향을 걷어낸 값을 같이 낸다.
    표준화도 학습 폴드 안에서만 계산한다. 전체로 계산하면 빼놓은 사람의 정보가 새어 든다.
    """
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
            mu.append(m)
            sd.append(sdev)
        X = [[(rows[i][c] - mu[j]) / sd[j] for j, c in enumerate(use)] for i in tr]
        w, b = logistic(X, [y[i] for i in tr])
        xh = [(rows[h][c] - mu[j]) / sd[j] for j, c in enumerate(use)]
        oof[h] = b + sum(w[j] * xh[j] for j in range(len(w)))
    return auc(y, oof)


def run(label, rows, cols):
    n = len(rows)
    idx = list(range(n))

    def fit(ix, use):
        y = [rows[i]['y'] for i in ix]
        if len(set(y)) < 2:
            return float('nan')
        z = [zscore([rows[i][c] for i in ix]) for c in use]
        X = [[z[j][t] for j in range(len(use))] for t in range(len(ix))]
        w, b = logistic(X, y)
        return auc(y, score(X, w, b))

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
            boots[c].append(fit(ix, ['age', c]) - bb)

    out = {'n': n, 'n_sensitive': sum(r['y'] for r in rows),
           'auc_age_only': round(base, 3), 'auc_age_only_loocv': round(base_cv, 3),
           'loocv_note': '재대입(fit)은 같은 42명으로 맞추고 같은 42명으로 잰다. 변수를 더하면 거의 항상 오르므로 증가분이 후하게 나온다. LOOCV 는 그 편향을 걷어낸 값이다. 둘의 방향이 같아야 결론을 믿을 수 있다.',
           'candidates': {}}
    print(f'\n== {label} — n={n} · 나이만 AUC {base:.3f} (교차검증 {base_cv:.3f}) ==')
    print(f"  {'후보':<22}{'AUC':>8}{'증가분':>10}   증가분 95% CI{'':>10}{'교차검증':>9}{'증가분':>9}")
    for c in cols:
        full = fit(idx, ['age', c])
        d = full - base
        s = sorted(v for v in boots[c] if v == v)
        lo, hi = s[int(0.025 * len(s))], s[int(0.975 * len(s))]
        crosses = lo <= 0 <= hi
        full_cv = loocv(rows, ['age', c])
        d_cv = full_cv - base_cv
        out['candidates'][c] = {'auc_age_plus': round(full, 3), 'delta': round(d, 3),
                                'ci': [round(lo, 3), round(hi, 3)],
                                'crosses_zero': crosses,
                                'auc_age_plus_loocv': round(full_cv, 3),
                                'delta_loocv': round(d_cv, 3),
                                'loocv_agrees': (d_cv > 0.0) == (not crosses)}
        mark = '0 포함' if crosses else '** 0 제외 **'
        print(f'  {c:<22}{full:>8.3f}{d:>+10.3f}   [{lo:+.3f}, {hi:+.3f}]  {mark:<12}'
              f'{full_cv:>9.3f}{d_cv:>+9.3f}')
    return out


print('=' * 88)
print('나이를 보정하면 후보가 무엇을 더하는가 — 03 민감피부 42명')
print('=' * 88)

out = {'seed': SEED, 'n_boot': N_BOOT,
       'note': '증가분 CI 가 0 을 포함하면 그 후보가 나이 위에 무엇을 더한다고 말할 수 없다. '
               'n=42 라 검정력이 거의 없으므로 "모른다" 의 크기를 재는 계산이다.'}

out['bacteria'] = run('세균 (16S V4-V5)', load('16S', CAND), CAND)
out['fungi'] = run('진균 (ITS1)', load('ITS', FUNGI), FUNGI)

print('\n읽는 법')
print('  증가분 CI 가 0 을 포함 → 나이 위에 더한다고 말할 수 없다')
print('  동시에 상한이 크면 → 안 더한다고도 말할 수 없다')
print('  교차검증 증가분 → 낙관 편향을 걷어낸 값. 재대입과 방향이 같아야 결론을 믿을 수 있다')

save_json(out, os.path.join(TAB, '나이보정_판별력.json'))
print('\n저장: results/나이보정_판별력.json')
