#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
③ 랜덤포레스트로 「나이만」 대 「나이 + 균」 — 반복 교차검증 비교 (2026-10-07 작성 · 사용자 합의)

왜: 데이터 멘토가 랜덤포레스트를 제안했다. 지금까지는 「42명에서 과적합 위험이 커 쓰지 않았다(해 보지 않았다)」고만
    답할 수 있었다. 실제로 돌려 보고, 같은 분할에서 로지스틱과 나란히 놓는다.

방법
  - 입력: 9번과 같은 42명 · 나이 + 후보 비율(%) — 나무 모형이라 표준화하지 않는다
  - 교차검증: 층화 5겹 × 50회 반복(seed 고정). 반복마다 5겹의 예측을 모아(out-of-fold) AUC 하나를 낸다
    → 반복 50개의 AUC 분포. 모든 모형이 **같은 분할**을 쓴다(증가분은 반복마다 짝지어 뺀다)
  - 랜덤포레스트: 나무 300그루 · 나머지는 scikit-learn 기본값 · random_state 고정
  - 비교 기준: 로지스틱(벌점 없음, 학습 폴드 안에서 표준화) — 같은 분할
  ⚠️ 50회 반복의 2.5–97.5 백분위는 「분할을 바꿨을 때의 흔들림」이지 모집단 신뢰구간이 아니다.

  병렬(2026-10-07 수정): 처음 판은 나무 300그루를 나눠 키우게(n_jobs=-1) 했으나 표본이 33명이라 사실상 코어 1개만 썼고
    25분을 넘겨도 끝나지 않았다. 반복 50회를 코어에 나눠 동시에 돌린다(joblib). 분할 · seed · 모형 설정은 같아
    결과는 순서대로 돌린 것과 같다. 메모리 기록은 주 프로세스 기준이다(작업 프로세스 제외).

결과: ../../results/랜덤포레스트_비교_261007.json · ../../results/랜덤포레스트_비교_261007.tsv (행 = 후보 × 모형, matrix 형태)
의존성: numpy · scikit-learn · joblib (conda bio_proj) + ./common42_261007.py
"""
import os, sys, time, warnings, tracemalloc

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

T0 = time.perf_counter()
tracemalloc.start()

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common42_261007 import load, save_json, TAB, CAND, FUNGI, SEED, mem_time_line   # noqa: E402

import numpy as np                                                        # noqa: E402
from sklearn.ensemble import RandomForestClassifier                       # noqa: E402
from sklearn.linear_model import LogisticRegression                       # noqa: E402
from sklearn.pipeline import make_pipeline                                # noqa: E402
from sklearn.preprocessing import StandardScaler                          # noqa: E402
from sklearn.model_selection import RepeatedStratifiedKFold               # noqa: E402
from sklearn.metrics import roc_auc_score                                 # noqa: E402
from joblib import Parallel, delayed                                      # noqa: E402

N_SPLITS, N_REPEATS, N_TREES = 5, 50, 300
N_JOBS = max(1, (os.cpu_count() or 2) - 2)     # 코어 2개는 남겨 둔다


def make(kind):
    if kind == 'rf':   # 작업 프로세스 안에서는 한 코어로 — 병렬은 반복 단위에서 한다
        return RandomForestClassifier(n_estimators=N_TREES, random_state=SEED, n_jobs=1)
    return make_pipeline(StandardScaler(), LogisticRegression(penalty=None, solver='lbfgs', max_iter=10000))


def one_repeat(X, y, rsplits, kind):
    """반복 하나: 5겹의 out-of-fold 확률을 모아 AUC 하나."""
    oof = np.zeros(len(y))
    for tr, te in rsplits:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            m = make(kind).fit(X[tr], y[tr])
        oof[te] = m.predict_proba(X[te])[:, 1]
    return roc_auc_score(y, oof)


def cv_aucs(X, y, splits, kind):
    """반복 수만큼의 AUC 목록 — 반복을 코어에 나눠 동시에 계산(결과 순서는 반복 순서 그대로)."""
    jobs = [splits[r * N_SPLITS:(r + 1) * N_SPLITS] for r in range(N_REPEATS)]
    return np.array(Parallel(n_jobs=N_JOBS)(delayed(one_repeat)(X, y, s, kind) for s in jobs))


def summ(a):
    return {'mean': round(float(a.mean()), 3), 'p2_5': round(float(np.percentile(a, 2.5)), 3),
            'p97_5': round(float(np.percentile(a, 97.5)), 3)}


def run(label, rows, cols):
    # 16S 와 ITS 는 행 순서가 달라 그대로 두면 같은 42명이 다른 분할로 나뉜다(2026-10-07 독립 검수).
    # 사람(subject) 순서로 정렬해 두 자료가 같은 분할을 쓰게 한다 → 「나이만」 값도 두 블록에서 같아진다.
    rows = sorted(rows, key=lambda r: str(r['subject']))
    y = np.array([r['y'] for r in rows])
    age = np.array([[r['age']] for r in rows], dtype=float)
    splits = list(RepeatedStratifiedKFold(n_splits=N_SPLITS, n_repeats=N_REPEATS, random_state=SEED).split(age, y))
    out = {'n': int(len(y)), 'n_sensitive': int(y.sum()), 'candidates': {}}
    base = {k: cv_aucs(age, y, splits, k) for k in ('rf', 'logistic')}
    out['age_only'] = {k: summ(v) for k, v in base.items()}
    print(f"\n== {label} — n={len(y)} · 나이만  RF {out['age_only']['rf']['mean']:.3f} "
          f"· 로지스틱 {out['age_only']['logistic']['mean']:.3f} (반복 {N_REPEATS}회 평균) ==")
    print(f"  {'후보':<18}{'모형':<10}{'나이+균':>9}{'증가분 평균':>12}  증가분 2.5–97.5%{'':>4}{'증가분>0 비율':>14}")
    for c in cols:
        X = np.array([[r['age'], r[c]] for r in rows], dtype=float)
        res = {}
        for k in ('rf', 'logistic'):
            full = cv_aucs(X, y, splits, k)
            d = full - base[k]
            res[k] = {'auc_age_plus': summ(full), 'delta': summ(d),
                      'share_delta_positive': round(float((d > 0).mean()), 3)}
            print(f"  {c:<18}{k:<10}{res[k]['auc_age_plus']['mean']:>9.3f}{res[k]['delta']['mean']:>+12.3f}"
                  f"  [{res[k]['delta']['p2_5']:+.3f}, {res[k]['delta']['p97_5']:+.3f}]"
                  f"{res[k]['share_delta_positive']:>14.2f}")
        out['candidates'][c] = res
    return out


def main():
    print('=' * 92)
    print(f'③ 랜덤포레스트 대 로지스틱 — 「나이만」 대 「나이 + 균」 · 층화 {N_SPLITS}겹 × {N_REPEATS}회 · 03 민감피부 42명'
          f' · 병렬 {N_JOBS}코어')
    print('=' * 92)

    out = {'seed': SEED, 'n_splits': N_SPLITS, 'n_repeats': N_REPEATS, 'n_trees': N_TREES, 'n_jobs': N_JOBS,
           'note': '반복마다 out-of-fold AUC 하나. 범위(2.5–97.5%)는 분할을 바꿨을 때의 흔들림이지 모집단 신뢰구간이 아니다. '
                   '모든 모형이 같은 분할을 쓰고 증가분은 반복마다 짝지어 뺐다. 두 자료를 사람(subject) 순서로 정렬해 같은 분할을 쓴다. '
                   'n=42 라 랜덤포레스트는 흔들림이 크다.',
           'bacteria': run('세균 (16S V4-V5)', load('16S', CAND), CAND),
           'fungi': run('진균 (ITS1)', load('ITS', FUNGI), FUNGI)}

    save_json(out, os.path.join(TAB, '랜덤포레스트_비교_261007.json'))

    tsv = os.path.join(TAB, '랜덤포레스트_비교_261007.tsv')
    with open(tsv, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('candidate\tmodel\tauc_age_only_mean\tauc_age_plus_mean\tdelta_mean\tdelta_p2_5\tdelta_p97_5\tshare_delta_positive\n')
        for grp in ('bacteria', 'fungi'):
            for c, res in out[grp]['candidates'].items():
                for k, m in res.items():
                    fh.write(f"{c}\t{k}\t{out[grp]['age_only'][k]['mean']}\t{m['auc_age_plus']['mean']}\t{m['delta']['mean']}\t"
                             f"{m['delta']['p2_5']}\t{m['delta']['p97_5']}\t{m['share_delta_positive']}\n")

    print('\n저장: results/랜덤포레스트_비교_261007.json · results/랜덤포레스트_비교_261007.tsv')
    print(mem_time_line(T0) + ' (주 프로세스 기준 — 병렬 작업 프로세스 제외)')
    tracemalloc.stop()


if __name__ == '__main__':   # Windows 병렬 작업 프로세스가 본문을 다시 실행하지 않게
    main()
