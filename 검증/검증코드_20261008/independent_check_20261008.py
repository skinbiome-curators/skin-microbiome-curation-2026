"""핵심 수치 독립 재계산 (클로드 검증 v4, 2026-10-08)

팀 코드(code/, 추가분석/)를 import 하지 않고 data/ 원자료에서 발표 핵심 수치를 다시 계산한다.
정의(AUC 방향 · 대상 시료)만 팀 문서에서 가져오고 계산은 scipy · statsmodels · scikit-learn 으로 새로 짰다.
- AUC = P(민감 > 비민감) + 0.5·P(같음). 0.5 미만이면 민감군에서 적다.
- 부트스트랩은 군 안에서 사람 단위로 다시 뽑는다(층화), 4,000회, 시드 20261008, 백분위 구간.
- 595명 검출 인원만은 원 xlsx 대신 팀이 만든 results/595_genus_relabund.tsv 를 읽는다(반독립 — 표에 표시).
출력: independent_values_20261008.tsv (matrix: 행 = 지표, 열 = 보고값 · 재계산값 · 차이 · 허용오차 · 일치)

실행: conda activate HH_python3 후  python independent_check_20261008.py --repo <저장소> --out <출력폴더>
"""
import argparse, math, os, sys, time
import numpy as np
import pandas as pd
import psutil
from scipy import stats
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score

T0 = time.time()
SEED, NBOOT = 20261008, 4000
S, N = "Sensitive", "Non-sensitive"


def auc(x_s, x_n):
    """P(민감 > 비민감) + 0.5·동점 = U(민감) / (n_s·n_n)."""
    u = stats.mannwhitneyu(x_s, x_n, alternative="two-sided", method="asymptotic").statistic
    return u / (len(x_s) * len(x_n))


def boot_ci(x_s, x_n, rng):
    x_s, x_n = np.asarray(x_s, float), np.asarray(x_n, float)
    v = [auc(rng.choice(x_s, len(x_s)), rng.choice(x_n, len(x_n))) for _ in range(NBOOT)]
    return np.percentile(v, [2.5, 97.5])


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def logit_auc(y, X):
    """z-표준화한 X 로 로지스틱 MLE → 재대입 AUC, 로그우도."""
    Xz = (X - X.mean(0)) / X.std(0, ddof=0)
    m = sm.Logit(y, sm.add_constant(Xz)).fit(disp=0)
    return roc_auc_score(y, m.predict(sm.add_constant(Xz))), m.llf


def loocv_auc(y, X):
    """한 명씩 빼고 학습(학습 자료 안에서 표준화) → 빠진 사람 예측 → 예측을 모아 AUC 하나."""
    pred = np.empty(len(y))
    for i in range(len(y)):
        tr = np.arange(len(y)) != i
        mu, sd = X[tr].mean(0), X[tr].std(0, ddof=0)
        m = sm.Logit(y[tr], sm.add_constant((X[tr] - mu) / sd)).fit(disp=0)
        xi = np.concatenate([[1.0], (X[i] - mu) / sd])
        pred[i] = 1 / (1 + math.exp(-float(xi @ m.params)))
    return roc_auc_score(y, pred)


def hm_se(a, n1, n2):
    """Hanley–McNeil 1982 AUC 표준오차."""
    q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
    return math.sqrt((a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n2 - 1) * (q2 - a * a)) / (n1 * n2))


def icc1(df, value):
    """일원 랜덤효과 ICC(1): 사람(subject) × 시점, 모든 시점이 있는 사람만."""
    # 검수 반영: 사람×계절 중복이 있으면 조용히 평균내지 않고 멈춘다
    assert not df.duplicated(["subject", "season"]).any(), "사람×계절 중복 시료"
    w = df.pivot_table(index="subject", columns="season", values=value, aggfunc="mean").dropna()
    k, n = w.shape[1], w.shape[0]
    arr = w.to_numpy(float)
    grand = arr.mean()
    msb = k * ((arr.mean(1) - grand) ** 2).sum() / (n - 1)
    msw = ((arr - arr.mean(1)[:, None]) ** 2).sum() / (n * (k - 1))   # numpy 배열에서만 [:, None] 사용
    return (msb - msw) / (msb + (k - 1) * msw), n, k


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    D = lambda f: os.path.join(a.repo, "data", f)
    rng = np.random.default_rng(SEED)
    rows = []   # (지표, 보고값, 재계산값, 허용오차, 비고)

    # ── 03 민감 42명 · ITS ──────────────────────────────────────────────
    its = pd.read_csv(D("03_sensitive_ITS_genus_pct.tsv"), sep="\t", index_col="run")
    info = pd.read_csv(D("03_sensitive_ITS_sample_info.tsv"), sep="\t", index_col="run").loc[its.index]
    g = info["group"].to_numpy()
    rows.append(("42명 ITS 시료 수 (민감/비민감)", "23/19", f"{(g == S).sum()}/{(g == N).sum()}", None, ""))
    mal = its["g__Malassezia"].to_numpy(float)
    a_mal = auc(mal[g == S], mal[g == N])
    lo, hi = boot_ci(mal[g == S], mal[g == N], rng)
    rows += [("Malassezia AUC", 0.153, a_mal, 0.001, ""),
             ("Malassezia AUC 95% CI 하한", 0.048, lo, 0.01, "몬테카를로 오차"),
             ("Malassezia AUC 95% CI 상한", 0.286, hi, 0.01, "몬테카를로 오차"),
             ("Malassezia 중앙값 비민감 %", 84.2, np.median(mal[g == N]), 0.05, ""),
             ("Malassezia 중앙값 민감 %", 36.4, np.median(mal[g == S]), 0.05, "")]
    genera = [c for c in its.columns if c != "Unassigned" and (its[c] > 0).any()]
    p_all = [stats.mannwhitneyu(its.loc[g == S, c], its.loc[g == N, c], method="asymptotic",
                                use_continuity=False).pvalue for c in genera]
    q_all = dict(zip(genera, bh(p_all)))
    rows.append((f"Malassezia BH q — 검출 진균 {len(genera)}속 전체", 0.030, q_all["g__Malassezia"], 0.002,
                 "팀 보고 q=0.0003은 1차 통과 2속 범위"))

    # 종 수준 — ASV read ÷ 시료 전체 read
    cnt = pd.read_csv(D("03_sensitive_ITS_asv_count.tsv"), sep="\t", index_col="ASV")[its.index]
    tax = pd.read_csv(D("03_sensitive_ITS_asv_taxonomy.tsv"), sep="\t", index_col="ASV").loc[cnt.index]
    is_mal = (tax["Genus"] == "g__Malassezia").to_numpy()
    is_res = is_mal & (tax["Species"] == "s__restricta").to_numpy()
    tot = cnt.sum(0).to_numpy(float)
    res = cnt[is_res].sum(0).to_numpy(float) / tot * 100
    malr = cnt[is_mal].sum(0).to_numpy(float)
    frac = np.divide(cnt[is_res].sum(0).to_numpy(float), malr, out=np.full_like(malr, np.nan), where=malr > 0)
    ok = ~np.isnan(frac)
    rows += [("M. restricta AUC", 0.165, auc(res[g == S], res[g == N]), 0.002, ""),
             ("M. restricta 중앙값 민감 / 비민감 %", "30.2/66.9",
              f"{np.median(res[g == S]):.1f}/{np.median(res[g == N]):.1f}", None, ""),
             ("Malassezia read 중 종 이름 있는 비율 %", 91.2,
              cnt[is_mal & tax["Species"].notna().to_numpy()].to_numpy().sum() / malr.sum() * 100, 0.1, ""),
             ("속 안 restricta 비율 AUC (군간)", 0.444, auc(frac[ok & (g == S)], frac[ok & (g == N)]), 0.002,
              "재검증 문서 값")]

    # 진균 다양성 — 속 235열(Unassigned 포함) Shannon(ln)
    p = its.to_numpy(float)
    p = p / p.sum(1, keepdims=True)
    shannon = -(np.where(p > 0, p * np.log(np.where(p > 0, p, 1)), 0)).sum(1)
    rows.append(("ρ(Malassezia, 진균 Shannon)", -0.979, stats.spearmanr(mal, shannon).statistic, 0.001, ""))

    # 나이 보정 로지스틱
    y = (g == S).astype(float)
    age = info["age"].to_numpy(float)
    X1, X2 = age[:, None], np.column_stack([age, mal])
    a1, ll1 = logit_auc(y, X1)
    a2, ll2 = logit_auc(y, X2)
    lr = 2 * (ll2 - ll1)
    rows += [("나이만 AUC (재대입)", 0.653, a1, 0.001, ""),
             ("나이 + Malassezia AUC (재대입)", 0.883, a2, 0.001, ""),
             ("우도비 검정 LR", 18.50, lr, 0.05, ""),
             ("우도비 검정 p", 1.7e-5, stats.chi2.sf(lr, 1), 2e-6, ""),
             ("LOOCV 나이만 (MLE, 예측 모아 AUC)", 0.581, loocv_auc(y, X1), 0.005, ""),
             ("LOOCV 나이 + Malassezia (MLE)", 0.847, loocv_auc(y, X2), 0.005, "팀 경사하강 구현은 0.854")]

    # 정규성 — 후보 4속 × 2군
    b16 = pd.read_csv(D("03_sensitive_16S_genus_pct.tsv"), sep="\t", index_col="run")
    i16 = pd.read_csv(D("03_sensitive_16S_sample_info.tsv"), sep="\t", index_col="run").loc[b16.index]
    g16 = i16["group"].to_numpy()
    dists = [b16[c].to_numpy(float)[g16 == grp] for c in ["Cutibacterium", "Staphylococcus", "Corynebacterium"]
             for grp in (S, N)] + [mal[g == S], mal[g == N]]
    rows.append(("정규성 기각 수 (Shapiro p<0.05, 8개 분포)", 7, sum(stats.shapiro(d).pvalue < 0.05 for d in dists), 0, ""))
    rows.append(("후보 4속에서 0인 사람 수", 0, int(sum((d == 0).sum() for d in dists)), 0, ""))

    # 세균 3속 AUC (42명 16S)
    for c, rep in [("Cutibacterium", 0.487), ("Staphylococcus", 0.449), ("Corynebacterium", 0.389)]:
        x = b16[c].to_numpy(float)
        rows.append((f"{c} AUC", rep, auc(x[g16 == S], x[g16 == N]), 0.001, ""))

    # 595명 검출 인원 (반독립: 팀이 만든 분석표)
    t595 = pd.read_csv(os.path.join(a.repo, "results", "595_genus_relabund.tsv"), sep="\t", index_col=0)
    if t595.shape[0] != 595:            # 방향이 속 × 사람이면 뒤집는다
        t595 = t595.T
    for c, rep in [("Cutibacterium", 590), ("Staphylococcus", 591), ("Corynebacterium", 538)]:
        col = [k for k in t595.columns if str(k).split(";")[-1].strip() in (c, f"g__{c}")]
        rows.append((f"595명 {c} 검출 인원 [반독립]", rep, int((t595[col[0]] > 0).sum()) if col else "열 없음", 0, ""))

    # 반복측정 ICC — 15계절 이마
    s16 = pd.read_csv(D("15_season_16S_genus_pct.tsv"), sep="\t", index_col="run")
    m16 = pd.read_csv(D("15_season_16S_sample_info.tsv"), sep="\t", index_col="run").loc[s16.index]
    f16 = m16["body_site"].eq("Forehead") & m16["group"].ne("음성대조")   # 검수 반영: 팀 12번과 같은 제외 조건
    d16 = pd.DataFrame({"subject": m16.loc[f16, "subject"],
                        "season": m16.loc[f16, "collection_date"].str[5:7]})
    sIT = pd.read_csv(D("15_season_ITS_genus_pct.tsv"), sep="\t", index_col="run")
    mIT = pd.read_csv(D("15_season_ITS_sample_info.tsv"), sep="\t", index_col="run").loc[sIT.index]
    fIT = mIT["body_site"].eq("Forehead") & mIT["group"].ne("음성대조")
    dIT = pd.DataFrame({"subject": mIT.loc[fIT, "subject"].str.upper(),
                        "season": mIT.loc[fIT, "isolation_source"].str.split().str[-1]})
    for c, rep in [("Cutibacterium", 0.312), ("Staphylococcus", 0.102), ("Corynebacterium", 0.381)]:
        v, n, k = icc1(d16.assign(v=s16.loc[f16, c].to_numpy(float)), "v")
        rows.append((f"ICC 15계절 이마 {c} (n={n}, k={k})", rep, v, 0.01, ""))
    v, n, k = icc1(dIT.assign(v=sIT.loc[fIT, "g__Malassezia"].to_numpy(float)), "v")
    rows.append((f"ICC 15계절 이마 Malassezia (n={n}, k={k})", 0.195, v, 0.01, "ITS 계절은 isolation_source 끝 단어"))

    # 검정력 — Hanley–McNeil, 양측 α 0.05, 80%
    z_a, z_b = stats.norm.ppf(0.975), stats.norm.ppf(0.80)
    n57 = next(n for n in range(5, 1000) if 0.15 >= z_a * hm_se(0.5, n, n) + z_b * hm_se(0.65, n, n))
    n65 = next(n for n in range(5, 1000) if z_a * hm_se(0.5, n, n) <= 0.10)
    rows += [("CI 반폭 지금 (23/19, AUC 0.5 가정)", 0.1775, z_a * hm_se(0.5, 23, 19), 0.0005, ""),
             ("군당 n — AUC 0.65 · 검정력 80%", 57, n57, 0, ""),
             ("군당 모집 n — 57 ÷ (93/120)", 74, math.ceil(n57 / (93 / 120)), 0, ""),
             ("군당 n — CI 반폭 ≤ 0.10", 65, n65, 0, ""),
             ("군당 모집 n — 65 ÷ (93/120)", 84, math.ceil(n65 / (93 / 120)), 0, "README 「65 → 74」 서술 점검")]

    # 저장 — matrix 형태
    out = []
    for name, rep, val, tol, note in rows:
        if tol is None or isinstance(rep, str) or isinstance(val, str):
            match = "확인" if str(rep) == str(val) else "다름"
            out.append((name, rep, val, "", "", match, note))
        else:
            d = abs(float(val) - float(rep))
            # 검수 반영: 작은 p값이 반올림으로 뭉개지지 않게 유효숫자 6자리, 정수 지표는 정수로 저장
            shown = int(round(float(val))) if tol == 0 else float(f"{float(val):.6g}")
            out.append((name, rep, shown, float(f"{d:.6g}"), tol, "일치" if d <= tol + 1e-12 else "불일치", note))
    os.makedirs(a.out, exist_ok=True)
    pd.DataFrame(out, columns=["metric", "reported", "recomputed", "abs_diff", "tolerance", "match", "note"]) \
        .to_csv(os.path.join(a.out, "independent_values_20261008.tsv"), sep="\t", index=False)
    for r in out:
        print(f"{r[5]:4s} | {r[0]} | 보고 {r[1]} | 재계산 {r[2]}")
    print(f"total {(time.time() - T0) / 3600:.4f} Hour · peak {psutil.Process().memory_info().peak_wset / 1024**3:.3f} GB "
          f"· python {sys.version.split()[0]} · numpy {np.__version__} · scipy {sys.modules['scipy'].__version__}")


if __name__ == "__main__":
    main()
