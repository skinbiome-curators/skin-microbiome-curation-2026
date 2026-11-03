# 덱 v12 — 결과 3(공개 화장품 시험) p값 독립 재계산 (2026-11-03)

- **왜:** 덱 9장 «독립 코드로 재계산한 결과가 발표 수치와 일치». 결과 3의 p 0.18 · 0.75(8장 꼬리말 · 부록 A2-2)는 독립 재계산 기록이 없었음(`경진대회/효능판정기준_초안_261017.md` 202행 · 덱 v12 검증에서 지적)
- **어떻게:** `code/24` 를 import 하지 않고, 명세만 보고 표준 라이브러리로 새로 짬. 명세는 덱 8장 꼬리말과 효능판정기준 초안 #4 로그비
  - 6속 제외: *Wolbachia* · *Sphingomonas* · *Akkermansia* · *Acetobacter* · *Bacteroides* · *Lachnospiraceae NK4A136 group*
  - 로그비 = log((후보 + 1e-5) / (나머지 + 1e-5)) — 나머지 = 후보 3속 · 6속을 뺀 모든 속
  - 검정: 0W → 4W 사람별 차이로 정확 Wilcoxon 부호순위(동점 · 0 없음 확인). 주장 방향은 *Cutibacterium* · *Staphylococcus* 증가 = 한쪽, *Corynebacterium* = 양측
  - Holm 보정: 3속
- **실행:** `python -I 검증/덱v12_07재계산_261103/recalc_07.py <저장소 루트>` → [`결과.txt`](결과.txt)
- **결과:** `results/효능판정_시범.json` (cohort07 · 0W→4W · excl_1pct)과 일치

| 속 | 주장 방향 Holm p | 양측 Holm p | 저장소 값 |
|---|---|---|---|
| *Cutibacterium* | **0.1817** | 0.3635 | p_holm 0.1817 · p_two_sided_holm 0.3635 |
| *Staphylococcus* | **0.7461** | 0.5249 | 덱 0.75 |

- **참고(덱에 없음):** 같은 계산에서 *Corynebacterium* 의 0W → 4W 변화는 증가 쪽으로 양측 Holm p 0.0007(W+ 290 / 325). 원 보고가 주장한 균이 아니고, 회차 겹침 때문에 등급(근거 불충분)은 그대로임
