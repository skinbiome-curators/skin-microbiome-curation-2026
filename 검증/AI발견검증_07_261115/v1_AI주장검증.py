"""AI 발견 검증 — 07 (2026-11-15)

AI 비교 시험(검증/AI비교시험_261113)에서 AI 가 07 에 대해 낸 주장 14개를
원본 파일로 독립 계산해 판정한다. AI 가 쓴 코드는 보지 않는다.

판정(코드 전 사용자 합의 · 11/15)
- 맞음: 보고된 자릿수까지 일치(인원 · 개수는 정확히)
- 근사: 방향 · 크기는 맞고 값이 조금 다름(계산법 차이) — 차이를 적음
- 틀림
- 확인 불가: 자료로 가를 수 없는 해석이거나 주장의 정의가 없음 — 사실 부분만 적음

입력: data/07_skincare_16S_sample_info.tsv · data/07_skincare_16S_genus_pct.tsv · results/효능판정_시범_런번호.tsv
출력: 이 폴더의 결과.tsv · 결과.md
실행: python -I 검증/AI발견검증_07_261115/v1_AI주장검증.py
"""
import sys
import time
import tracemalloc
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import braycurtis
from scipy.stats import kruskal, pearsonr, wilcoxon

sys.stdout.reconfigure(encoding="utf-8")
T0 = time.time()
tracemalloc.start()

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
info = pd.read_csv(REPO / "data/07_skincare_16S_sample_info.tsv", sep="\t", dtype={"subject": str})
gen = pd.read_csv(REPO / "data/07_skincare_16S_genus_pct.tsv", sep="\t").set_index("run")
runs = pd.read_csv(REPO / "results/효능판정_시범_런번호.tsv", sep="\t", dtype=str)
runs = runs[runs["cohort"] == "07"].set_index("run")

assert len(info) == 75 and set(info["run"]) == set(gen.index) == set(runs.index), "07 시료 75개가 세 파일에서 같아야 함"
info = info.set_index("run")
info["run_number"] = runs.loc[info.index, "run_number"].astype(int)
gen = gen.loc[info.index]
rowsum = gen.sum(axis=1)
assert np.allclose(rowsum, 100, atol=0.5), f"속 비율 합이 100 이 아님: {rowsum.min():.2f}–{rowsum.max():.2f}"

TP = ["0W", "2W", "4W"]
SIX = ["Wolbachia", "Sphingomonas", "Akkermansia", "Acetobacter", "Lachnospiraceae NK4A136 group", "Bacteroides"]
R = []


def add(vid, claim, src, ours, verdict, note=""):
    R.append(dict(번호=vid, AI주장=claim, 답=src, 우리_계산=ours, 판정=verdict, 비고=note))


def num(ours, claimed, nd):
    """보고 자릿수 nd 에서 같으면 맞음 · 10% 안이면 근사 · 아니면 틀림"""
    if round(ours, nd) == round(claimed, nd):
        return "맞음"
    return "근사" if abs(ours - claimed) <= max(10 ** -nd, 0.1 * abs(claimed)) else "틀림"


def cnt(ours, claimed):
    return "맞음" if ours == claimed else ("근사" if abs(ours - claimed) == 1 else "틀림")


def worst(*vs):
    order = {"맞음": 0, "근사": 1, "틀림": 2}
    return max(vs, key=lambda v: order[v])


def wide(col_values):
    """시료 단위 값 → 사람 × 시점 표"""
    df = pd.DataFrame({"subject": info["subject"], "tp": info["timepoint"], "v": col_values})
    return df.pivot(index="subject", columns="tp", values="v")[TP]


# V1 대상자 · 나이 열
subs = sorted(info["subject"].unique())
expect = [f"S{i:02d}" for i in range(1, 27) if i != 10]
age_empty = info["age"].isna().all() and info["age_group"].isna().all()
add("V1", "25명 = S01–S26 중 S10 없음 · age · age_group 모두 비어 있음", "07_A_1 · 07_D_1",
    f"대상자 {len(subs)}명 · {'S01–S26 중 S10 없음' if subs == expect else subs} · age 빈칸 {age_empty}",
    "맞음" if subs == expect and age_empty else "틀림")

# V2 런 × 시점
ct = pd.crosstab(info["timepoint"], info["run_number"])
s244 = sorted(info[(info["timepoint"] == "2W") & (info["run_number"] == 244)]["subject"])
ok2 = (ct.loc["0W", 242] == 25 and ct.loc["2W", 241] == 22 and ct.loc["2W", 244] == 3 and ct.loc["4W", 243] == 25
       and s244 == ["S08", "S13", "S15"])
add("V2", "0주 242(25) · 2주 241(22) + 244(3: S08 · S13 · S15) · 4주 243(25)", "07_C_2 · 07_D_1",
    "; ".join(f"{t}: " + ", ".join(f"{r}({int(n)})" for r, n in ct.loc[t].items() if n) for t in TP) + f" · 244 = {s244}",
    "맞음" if ok2 else "틀림")

# V3 0주 여섯 속
z = info["timepoint"] == "0W"
g0 = gen[z]
claimed = {"Wolbachia": 27.6, "Sphingomonas": 17.7, "Akkermansia": 4.7, "Acetobacter": 4.6,
           "Lachnospiraceae NK4A136 group": 3.2, "Bacteroides": 2.3}
vs = [num(g0[k].mean(), v, 1) for k, v in claimed.items()]
add("V3a", "0주 평균: Wolbachia 27.6 · Sphingomonas 17.7 · Akkermansia 4.7 · Acetobacter 4.6 · Lachnospiraceae NK4A136 3.2 · Bacteroides 2.3 %",
    "07_A_2 · 07_C_2 · 07_D_2", " · ".join(f"{k.split()[0]} {g0[k].mean():.2f}" for k in claimed), worst(*vs))
six = g0[SIX].sum(axis=1)
add("V3b", "여섯 속 합: 0주 평균 60% · 중앙값 58%", "07_A_3",
    f"평균 {six.mean():.1f} · 중앙값 {six.median():.1f}", worst(num(six.mean(), 60, 0), num(six.median(), 58, 0)))
w0 = g0["Wolbachia"]
top_wo_un = (g0.drop(columns="Unassigned").idxmax(axis=1) == "Wolbachia").sum()
top_with_un = (g0.idxmax(axis=1) == "Wolbachia").sum()
add("V3c", "Wolbachia 0주 25명 전원 17.6–36.6% · 20명에게서 가장 많은 속", "07_A_3 · 07_C_2 · 07_D_1",
    f"범위 {w0.min():.1f}–{w0.max():.1f} · 검출 {int((w0 > 0).sum())}/25 · 최다 속인 사람 {top_wo_un}명(Unassigned 제외) / {top_with_un}명(포함)",
    worst(num(w0.min(), 17.6, 1), num(w0.max(), 36.6, 1),
          "맞음" if top_with_un == 20 and top_wo_un == 20 else ("근사" if 20 in (top_with_un, top_wo_un) else cnt(int(top_wo_un), 20))),
    "🔧 첫 실행 뒤 고침: ‘최다 속’에 Unassigned 를 넣으면 20명(AI 와 같음) · 빼면 23명 — 정의 차이라 근사(합의한 기준)")

# V4 Wolbachia : Sphingomonas
ratio = g0["Wolbachia"] / g0["Sphingomonas"]
cv = ratio.std(ddof=1) / ratio.mean() * 100
r = pearsonr(g0["Wolbachia"], g0["Sphingomonas"])[0]
add("V4", "Wolbachia/Sphingomonas 비 25명 모두 1.42–1.76 · 변동계수 약 5% · r = 0.985", "07_A_2 · 07_C_1 · 07_C_2 · 07_A_3",
    f"비 {ratio.min():.2f}–{ratio.max():.2f} · 중앙값 {ratio.median():.2f} · 변동계수 {cv:.1f}% · r = {r:.3f}",
    worst(num(ratio.min(), 1.42, 2), num(ratio.max(), 1.76, 2), num(cv, 5, 0), num(r, 0.985, 3)),
    "‘한 출처의 혼합(오염)’은 해석 — 비가 일정하다는 사실만 판정")

# V5 깊이 · Observed
med = info.groupby("timepoint")[["read_count", "read_depth_final", "Observed"]].median().loc[TP]
kw = kruskal(*[info.loc[info["timepoint"] == t, "read_count"] for t in TP]).pvalue
v5 = [cnt(int(med.loc[t, "read_count"]), c) for t, c in zip(TP, [148418, 131874, 214290])]
v5 += [cnt(int(med.loc[t, "read_depth_final"]), c) for t, c in zip(TP, [29878, 24442, 47764])]
v5 += [cnt(int(med.loc[t, "Observed"]), c) for t, c in zip(TP, [455, 270, 223])]
add("V5", "중앙값 read_count 148,418 → 131,874 → 214,290(KW p ≈ 3×10⁻¹⁰) · read_depth_final 29,878 → 24,442 → 47,764 · Observed 455 → 270 → 223",
    "07_D_2 · 07_A_3 · 07_D_1 · 07_D_3",
    " / ".join(f"{t}: {int(med.loc[t, 'read_count']):,} · {int(med.loc[t, 'read_depth_final']):,} · {int(med.loc[t, 'Observed'])}" for t in TP) + f" · KW p = {kw:.1e}",
    worst(*v5, "맞음" if kw < 1e-9 and kw > 1e-11 else "근사"))

# V6 · V7 2주 → 4주
st = wide(gen["Staphylococcus"])
d = st["4W"] - st["2W"]
p6 = wilcoxon(st["4W"], st["2W"]).pvalue
add("V6", "Staphylococcus 2→4주: 25명 중 19명 감소 · 6명 증가(p ≈ 0.017–0.02)", "07_A_1 · 07_A_2 · 07_C_1 · 07_D_3",
    f"감소 {int((d < 0).sum())} · 증가 {int((d > 0).sum())} · 같음 {int((d == 0).sum())} · Wilcoxon p = {p6:.3f}",
    worst(cnt(int((d < 0).sum()), 19), cnt(int((d > 0).sum()), 6), "맞음" if 0.01 <= p6 <= 0.025 else "근사"),
    "원 비율(%) 기준")
cu = wide(gen["Cutibacterium"])
add("V7", "Cutibacterium 2→4주: 25명 중 21명 증가", "07_A_2",
    f"증가 {int(((cu['4W'] - cu['2W']) > 0).sum())}", cnt(int(((cu["4W"] - cu["2W"]) > 0).sum()), 21))

# V8 Corynebacterium 대비 로그비 0 → 4주
co = gen["Corynebacterium"]
zeros = {k: int((gen[k] == 0).sum()) for k in ["Cutibacterium", "Staphylococcus", "Corynebacterium"]}
lcc = wide(np.log(gen["Cutibacterium"] / co))
lsc = wide(np.log(gen["Staphylococcus"] / co))
up_c, up_s = int((lcc["4W"] > lcc["0W"]).sum()), int((lsc["4W"] > lsc["0W"]).sum())
p_c, p_s = wilcoxon(lcc["4W"], lcc["0W"]).pvalue, wilcoxon(lsc["4W"], lsc["0W"]).pvalue
add("V8", "Corynebacterium 75개 모두 검출 · log(Cuti/Coryne) 0→4주 6/25 증가 · log(Staph/Coryne) 4/25 증가(p ≈ 0.004)", "07_C_2 · 07_C_3",
    f"0 인 시료 {zeros} · Cuti/Coryne 증가 {up_c}/25(p = {p_c:.3f}) · Staph/Coryne 증가 {up_s}/25(p = {p_s:.3f})",
    worst("맞음" if zeros["Corynebacterium"] == 0 else "틀림", cnt(up_c, 6), cnt(up_s, 4),
          num(p_c, 0.004, 3), num(p_s, 0.004, 3)))

# V9 Unassigned · Pseudomonas
un = info.groupby("timepoint").apply(lambda x: gen.loc[x.index, "Unassigned"].agg(["mean", "median"]))
ps = info.groupby("timepoint").apply(lambda x: gen.loc[x.index, "Pseudomonas"].agg(["mean", "median"]))
psw = wide(gen["Pseudomonas"])
up_ps = int((psw["4W"] > psw["0W"]).sum())
# 🔧 첫 실행 뒤 고침: 07_C_3 의 값은 「덩어리를 빼도 런 차이가 남습니다」 다음 문장 — 0주 덩어리를 빼고 다시 100% 로 맞춘 값이다.
#    AI 는 덩어리를 정확히 정의하지 않았다 → 우리 여섯 속을 빼고 다시 맞춰 근사 계산한다(맞음은 줄 수 없음).
rest = gen.drop(columns=SIX)
ren = rest.div(rest.sum(axis=1), axis=0) * 100
un_r = info.groupby("timepoint").apply(lambda x: ren.loc[x.index, "Unassigned"].agg(["mean", "median"]))
ps_r = info.groupby("timepoint").apply(lambda x: ren.loc[x.index, "Pseudomonas"].agg(["mean", "median"]))
best = lambda *v: min(v, key=lambda u: {"맞음": 0, "근사": 1, "틀림": 2}[u])
v9 = worst(best(num(un_r.loc["0W", "mean"], 31.0, 1), num(un_r.loc["0W", "median"], 31.0, 1)),
           best(num(un_r.loc["4W", "mean"], 17.2, 1), num(un_r.loc["4W", "median"], 17.2, 1)),
           best(num(ps_r.loc["0W", "mean"], 3.1, 1), num(ps_r.loc["0W", "median"], 3.1, 1)),
           best(num(ps_r.loc["4W", "mean"], 5.0, 1), num(ps_r.loc["4W", "median"], 5.0, 1)),
           cnt(up_ps, 25))
v9 = "확인 불가"   # 🔧 두 번째 실행 뒤: 여섯 속을 빼도 0주 값이 맞지 않음(27.6 대 31.0) — AI 의 덩어리 정의가 없어 재현 불가 · 4주 값과 증가 인원은 일치
add("V9", "(덩어리를 빼고 다시 맞추면) Unassigned 0주 31.0% → 4주 17.2% · Pseudomonas 3.1% → 5.0% · Pseudomonas 0→4주 25/25명 증가(원 비율)", "07_C_3 · 07_A_3",
    f"여섯 속을 빼고 다시 맞춘 값 — Unassigned 평균 {un_r.loc['0W', 'mean']:.1f} → {un_r.loc['4W', 'mean']:.1f} (중앙값 {un_r.loc['0W', 'median']:.1f} → {un_r.loc['4W', 'median']:.1f}) · "
    f"Pseudomonas 평균 {ps_r.loc['0W', 'mean']:.1f} → {ps_r.loc['4W', 'mean']:.1f} (중앙값 {ps_r.loc['0W', 'median']:.1f} → {ps_r.loc['4W', 'median']:.1f}) · "
    f"원 비율 Pseudomonas 증가 {up_ps}/25 · (참고: 원 비율 Unassigned 평균 {un.loc['0W', 'mean']:.1f} → {un.loc['4W', 'mean']:.1f})",
    v9, "🔧 첫 실행 뒤 고침: 원 비율로 비교해 ‘틀림’ 이 났으나 답의 문맥상 덩어리를 뺀 값 — AI 의 덩어리 정의를 몰라 최대 ‘근사’")

# V10 Bray–Curtis
P = gen.div(gen.sum(axis=1), axis=0).to_numpy()
ids = list(info.index)
same_run_diff_sub, same_sub_diff_run = [], []
for i in range(len(ids)):
    for j in range(i + 1, len(ids)):
        a, b = info.iloc[i], info.iloc[j]
        dij = braycurtis(P[i], P[j])
        if a["run_number"] == b["run_number"] and a["subject"] != b["subject"]:
            same_run_diff_sub.append(dij)
        if a["subject"] == b["subject"] and a["run_number"] != b["run_number"]:
            same_sub_diff_run.append(dij)
m1, m2 = float(np.median(same_run_diff_sub)), float(np.median(same_sub_diff_run))
add("V10", "Bray–Curtis: 같은 런의 다른 사람 중앙값 0.434 < 같은 사람의 다른 시점 0.578", "07_C_3 · 07_C_2",
    f"같은 런 · 다른 사람 {m1:.3f}(쌍 {len(same_run_diff_sub)}) · 같은 사람 · 다른 런 {m2:.3f}(쌍 {len(same_sub_diff_run)})",
    worst(num(m1, 0.434, 3), num(m2, 0.578, 3)) if m1 < m2 else "틀림",
    "모든 속(Unassigned 포함) · 비율을 합 1 로 맞춤 · ‘같은 사람의 다른 시점’은 런도 다른 쌍만")

# V11 런 번호 순서
add("V11", "0주 런(242)이 2주 런(241)보다 번호가 뒤", "07_C_2 · 07_D_1",
    "242 > 241 · 기기 M04177 같음", "맞음",
    "‘0주 시료를 나중에 읽었다’는 해석 — 런 번호가 시간 순서인지는 자료로 확인 불가")

# V12 Lactobacillus
lb = gen["Lactobacillus"]
lb_mean = info.groupby("timepoint").apply(lambda x: lb.loc[x.index].mean())
lb_hi = lb[(info["timepoint"] != "0W") & (lb >= 10)]
add("V12", "Lactobacillus 0주 평균 0.19% → 4주 2.35% · 몇 명은 11–18%", "07_A_1 · 07_C_1",
    f"평균 0W {lb_mean['0W']:.2f} · 2W {lb_mean['2W']:.2f} · 4W {lb_mean['4W']:.2f} · 2 · 4주 10% 이상 {len(lb_hi)}개({lb_hi.min():.1f}–{lb_hi.max():.1f}%)" if len(lb_hi) else
    f"평균 0W {lb_mean['0W']:.2f} · 4W {lb_mean['4W']:.2f} · 2 · 4주 10% 이상 0개",
    worst(num(lb_mean["0W"], 0.19, 2), num(lb_mean["4W"], 2.35, 2),
          "맞음" if len(lb_hi) and lb_hi.min() >= 10.5 and lb_hi.max() < 18.5 else ("근사" if len(lb_hi) else "틀림")),
    "‘에센스 제품의 DNA’는 해석 — 확인 불가(제품 시퀀싱 필요)")

# V13 여러 속이 함께 증가 0 → 4주
named = ["Pseudomonas", "Streptococcus", "Corynebacterium", "Haemophilus", "Enhydrobacter"]
ups = {}
for g in named:
    w = wide(gen[g])
    ups[g] = int((w["4W"] > w["0W"]).sum())
add("V13", "Pseudomonas(25/25) · Streptococcus · Corynebacterium · Haemophilus · Enhydrobacter 가 0→4주에 함께 증가", "07_A_3",
    " · ".join(f"{g} {n}/25" for g, n in ups.items()),
    "맞음" if all(n >= 13 for n in ups.values()) and ups["Pseudomonas"] == 25 else ("근사" if sum(n >= 13 for n in ups.values()) >= 4 else "틀림"),
    "‘함께 증가’ = 25명 중 13명 이상(과반)에서 증가로 정함(11/15 · 계산 전)")

# V14 0주에만 있다가 사라지는 속
det0 = (g0 > 0).sum()
g4 = gen[info["timepoint"] == "4W"]
det4 = (g4 > 0).sum()
gone = det0[(det0 == 25) & (det4 <= 5)].index.tolist()
add("V14", "0주 25개에 모두 있다가 이후 거의 사라지는 속 20개 · 합 63%", "07_C_2",
    f"정의(0주 25/25 검출 · 4주 검출 5명 이하)로 {len(gone)}개 · 0주 합 평균 {g0[gone].sum(axis=1).mean():.1f}%",
    "확인 불가", "AI 가 ‘거의 사라짐’을 정의하지 않음 — 우리 정의로 낸 값은 참고")

out = pd.DataFrame(R)
out.to_csv(HERE / "결과.tsv", sep="\t", index=False)
cnts = out["판정"].value_counts().to_dict()
md = ["| 번호 | AI 주장 | 답 | 우리 계산 | 판정 | 비고 |", "|---|---|---|---|:-:|---|"]
md += [f"| {x.번호} | {x.AI주장} | {x.답} | {x.우리_계산} | {x.판정} | {x.비고} |" for x in out.itertuples()]
(HERE / "결과.md").write_text("\n".join(md) + "\n", encoding="utf-8")
print(out[["번호", "판정", "우리_계산"]].to_string(index=False))
print("판정 합계:", cnts)
cur, peak = tracemalloc.get_traced_memory()
print(f"메모리(Python 할당 기준) 현재 {cur / 1e9:.4f} GB · peak {peak / 1e9:.4f} GB · 소요 시간 {(time.time() - T0) / 3600:.6f} Hour")
