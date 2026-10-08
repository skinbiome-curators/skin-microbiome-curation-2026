#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
추가분석 20–23번이 함께 쓰는 42명 읽기 · 로그비 · 9번 로지스틱 (2026-10-07 작성)

같은 42명을 같은 규칙으로 읽어야 9번 결과와 나란히 비교할 수 있다. 그래서 새로 짜지 않고 옮겼다.
  load()        ← code/9_나이보정_판별력.py 의 load() 그대로 (같은 사람 · 같은 군 판정 · 같은 값)
  zscore() · logistic_gd() · score() · auc()
                ← code/9_나이보정_판별력.py 의 zscore · logistic · score · auc 그대로 (경사하강 200회 · lr 0.1)
  add_logratio()  세균 ← code/4_조성효과_점검.py 의 ④ 로그비(분모 = 후보 · 구강 의심 · 미분류를 뺀 나머지, EPS 1e-5)
                  진균 ← code/8_진균_판별력.py 의 ② 로그비(분모 = 그 속 · 미분류를 뺀 나머지, EPS 1e-3)
군 판정: 'Non' 이 있으면 비민감(y=0), 'ensitive' 이면 민감(y=1) — 9번과 같다.
"""
import csv, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'code'))
from stats_util import RAW, TAB, CAND, ORAL, f, save_json   # noqa: E402,F401

SEED = 20260918
FUNGI = ['g__Malassezia']
EPS_16S = 1e-5      # 4번 ④ 로그비와 같음
EPS_ITS = 1e-3      # 8번 ② 로그비와 같음


def read_tsv(path):
    with open(path, encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def load(prefix, cols):
    """03_sensitive 의 16S 또는 ITS 를 읽어 (행: 사람) 로 정리 — 9번 load() 와 같다."""
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
        # 9번과 다른 점은 subject 키 하나뿐(2026-10-07) — 16S 와 ITS 는 run 번호와 행 순서가 달라서,
        # 22번이 두 자료를 같은 사람 순서로 정렬해 같은 교차검증 분할을 쓰게 하려고 더했다. 계산에는 쓰지 않는다.
        rows.append(dict(run=run, subject=r.get('subject'), y=y, age=f(age), **vals))
    return rows


def add_logratio(prefix, rows, cols):
    """rows 의 각 사람에 'lr_<속>' 을 더한다. 분모 정의는 4번(세균) · 8번(진균)과 같다."""
    pct = read_tsv(os.path.join(RAW, f'03_sensitive_{prefix}_genus_pct.tsv'))
    allv = {r['run']: {k: (f(v) or 0.0) for k, v in r.items() if k != 'run'} for r in pct}
    for r in rows:
        v = allv[r['run']]
        if prefix == '16S':
            exc = set(CAND) | {g for g in ORAL if g in v} | {'Unassigned'}
            ref = sum(x for g, x in v.items() if g not in exc)
            for c in cols:
                r['lr_' + c] = math.log10((v[c] + EPS_16S) / (ref + EPS_16S))
        else:
            for c in cols:
                ref = sum(x for g, x in v.items() if g not in (c, 'Unassigned', 'unidentified'))
                r['lr_' + c] = math.log10((v[c] + EPS_ITS) / (ref + EPS_ITS))
    return rows


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


def logistic_gd(X, y, iters=200, lr=0.1):
    """9번 logistic() 그대로 — 고정 200회 경사하강. 수렴을 확인하지 않는다(20번이 표준 MLE 와 대조한다)."""
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


def mem_time_line(t0):
    """Grok 규칙 — 메모리(GB) · 소요 시간(Hour) 기록."""
    import time, tracemalloc
    cur, peak = tracemalloc.get_traced_memory()
    return (f'메모리(Python 할당 기준) 현재 {cur / 1024 ** 3:.4f} GB · peak {peak / 1024 ** 3:.4f} GB · '
            f'소요 시간 {(time.perf_counter() - t0) / 3600:.6f} Hour')
