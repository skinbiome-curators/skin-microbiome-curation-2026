#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
42명 16S·ITS ASV 표에서 속별로 «종 라벨이 붙은 read 비율»을 센다 (2026-10-02 작성)

발표준비_종수준한계_후보생물학근거_20260930.md §7 이 이 스크립트를 생성기로 가리키는데 파일이 없었다.
그 표를 같은 정의로 다시 계산해 재현할 수 있게 한다.

  종 라벨 = asv_taxonomy 의 Species 칸이 비어 있지 않은 ASV (ITS 는 's__' 접두어를 떼고 표시)
  read 비율 = 42명 count 합 기준 — (그 속에서 종 라벨 ASV 의 read 합) / (그 속 전체 read 합) × 100
  대상 = §7 표의 7개 속 + 자료 전체

⚠️ 종 라벨은 16S·ITS 짧은 구간 기반 추정이다. 라벨이 붙었다고 종이 확정된 것이 아니다
   (예: Staphylococcus 최다 ASV 는 BLAST 에서 여러 종과 동률 — 같은 폴더 점검 문서 §5).

입력: ../../data/03_sensitive_{16S,ITS}_asv_{taxonomy,count}.tsv
결과: 같은 폴더 속별_종라벨비율.tsv (행 = 자료·속, matrix 형태) + 화면 출력
의존성: Python 3 표준 라이브러리만
"""
import os, sys, csv, time, tracemalloc

try:   # 한글 Windows 콘솔(cp949)에서 출력이 깨지지 않게 — stats_util.py 와 같은 처리
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

T0 = time.perf_counter()
tracemalloc.start()

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'data')   # 2026-10-07 폴더 이름 변경: 원자료 → data
OUT = os.path.join(HERE, '속별_종라벨비율.tsv')
TOP = 4   # 표에 남길 상위 종 수

TARGETS = {'16S': ['Cutibacterium', 'Staphylococcus', 'Corynebacterium', 'Streptococcus', 'Lactobacillus'],
           'ITS': ['Malassezia', 'Cladosporium']}


def read_tsv(path):
    with open(path, encoding='utf-8-sig', newline='') as fh:
        return list(csv.reader(fh, delimiter='\t'))


def load(tag):
    tax = read_tsv(os.path.join(RAW, f'03_sensitive_{tag}_asv_taxonomy.tsv'))
    cnt = read_tsv(os.path.join(RAW, f'03_sensitive_{tag}_asv_count.tsv'))
    h = tax[0]
    gi, si = h.index('Genus'), h.index('Species')
    reads = {r[0]: sum(float(x) for x in r[1:] if x != '') for r in cnt[1:]}
    asvs = []
    for r in tax[1:]:
        genus = r[gi][3:] if r[gi].startswith('g__') else r[gi]
        sp = r[si][3:] if r[si].startswith('s__') else r[si]
        asvs.append((r[0], genus, sp.strip(), reads.get(r[0], 0.0)))
    missing = sum(1 for r in tax[1:] if r[0] not in reads)
    return asvs, missing


def summarize(asvs, label):
    total = sum(a[3] for a in asvs)
    named = [a for a in asvs if a[2]]
    by_sp = {}
    for a in named:
        by_sp[a[2]] = by_sp.get(a[2], 0.0) + a[3]
    top = sorted(by_sp.items(), key=lambda kv: (-kv[1], kv[0]))[:TOP]
    row = {'genus': label, 'asv_labeled': len(named), 'asv_total': len(asvs),
           'reads_labeled': int(round(sum(a[3] for a in named))), 'reads_total': int(round(total)),
           'labeled_read_pct': round(100 * sum(a[3] for a in named) / total, 2) if total else 0.0}
    for k in range(TOP):
        sp, rd = top[k] if k < len(top) else ('', 0.0)
        row[f'top{k + 1}_species'] = sp
        row[f'top{k + 1}_pct'] = round(100 * rd / total, 2) if (sp and total) else ''
    return row


rows = []
for tag, genera in TARGETS.items():
    asvs, missing = load(tag)
    if missing:
        print(f'⚠️ {tag}: taxonomy 에는 있으나 count 표에 없는 ASV {missing}개 → read 0 으로 처리')
    for g in genera:
        sub = [a for a in asvs if a[1] == g]
        rows.append({'data': tag, **summarize(sub, g)})
    rows.append({'data': tag, **summarize(asvs, '(전체)')})

cols = ['data', 'genus', 'asv_labeled', 'asv_total', 'reads_labeled', 'reads_total', 'labeled_read_pct'] + \
       [f'top{k + 1}_{x}' for k in range(TOP) for x in ('species', 'pct')]
with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(cols)
    for r in rows:
        w.writerow([r[c] for c in cols])

print(f"{'자료':4s} {'속':16s} {'ASV(라벨/전체)':>14s} {'종 라벨 read %':>14s}   상위 종")
for r in rows:
    tops = ' · '.join(f"{r[f'top{k + 1}_species']} {r[f'top{k + 1}_pct']}" for k in range(TOP) if r[f'top{k + 1}_species'])
    print(f"{r['data']:4s} {r['genus']:16s} {str(r['asv_labeled']) + '/' + str(r['asv_total']):>14s} {r['labeled_read_pct']:>14.2f}   {tops}")
print(f'\n저장: {os.path.relpath(OUT, HERE)}')

cur, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()
print(f'메모리(Python 할당 기준) 현재 {cur / 1024 ** 3:.4f} GB · peak {peak / 1024 ** 3:.4f} GB · '
      f'소요 시간 {(time.perf_counter() - T0) / 3600:.6f} Hour')
