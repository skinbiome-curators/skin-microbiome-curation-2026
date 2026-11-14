"""AI 비교 시험 — 집계 (설계서 v2 · 7절 · 8절)

채점자 2명의 JSON 을 합쳐, 비교 대상(A · C · D)별로 k/n 과 95% Clopper–Pearson 구간을 낸다.
채점자끼리 다르면 둘 다 남기고, AI 쪽에 유리한 판정을 채택한다(우리 주장에 보수적으로).

실행: python -I 검증/AI비교시험_261113/스크립트/s4_집계.py 07|설계
"""
import csv
import json
import re
import sys
import time
import tracemalloc
from collections import defaultdict
from pathlib import Path

from scipy.stats import beta

sys.stdout.reconfigure(encoding="utf-8")
T0 = time.time()
tracemalloc.start()
HERE = Path(__file__).resolve().parent.parent
SC = HERE / "채점"
P_RANK = {"P0": 0, "P1": 1, "P2": 2}
ALLOWED = {  # 설계서 6절 허용 판정
    "가": {"받칠 수 없음"}, "나": {"받칠 수 없음"}, "다": {"받칠 수 없음", "조건부"},
    "라": {"받칠 수 있음", "조건부"}, "마": {"받칠 수 없음", "조건부"},
}
FLAWS = {"가": ["가1"], "나": ["나1", "나2", "나3"], "다": ["다1", "다2"], "라": [], "마": ["마1", "마2"]}


def cp(k, n):
    if n == 0:
        return "—"
    lo = 0.0 if k == 0 else beta.ppf(0.025, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(0.975, k + 1, n - k)
    return f"{k}/{n} ({lo * 100:.0f}–{hi * 100:.0f}%)"


def scorer_json(test, k):
    p = SC / "원답" / f"채점_{test}_{k}.jsonl"
    txt = ""
    for line in p.read_text(encoding="utf-8").splitlines():
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            continue
        if d.get("type") == "result":
            txt = d.get("result") or ""
    m = re.search(r"```json\s*(\[.*?\])\s*```", txt, flags=re.S)
    if not m:
        raise SystemExit(f"채점자 {k}: JSON 블록을 찾지 못함")
    return {o["답"]: o for o in json.loads(m.group(1))}


def read_key(test):
    with (SC / f"{test}_키.tsv").open(encoding="utf-8", newline="") as f:
        return {r["답"]: (r["job"], r["비교대상"]) for r in csv.DictReader(f, delimiter="\t")}


def agg07():
    key, s1, s2 = read_key("07"), scorer_json("07", 1), scorer_json("07", 2)
    with (SC / "07_유효성.tsv").open(encoding="utf-8", newline="") as f:
        pscript = {r["job"]: r["P스크립트"] for r in csv.DictReader(f, delimiter="\t")}
    items = ["문장1", "문장2", "P", "S1", "S2", "S3", "S4", "X", "기억"]
    agree = defaultdict(lambda: [0, 0])
    per = []
    for aid, (job, arm) in key.items():
        a, b = s1.get(aid), s2.get(aid)
        if a is None or b is None:
            raise SystemExit(f"{aid}: 채점자 한쪽에 없음")
        for it in items:
            agree[it][1] += 1
            agree[it][0] += int(a.get(it) == b.get(it))
        # AI 쪽에 유리하게: P 는 높은 쪽 · S 는 참 쪽 · X 는 거짓 쪽 · 문장 판정은 엔진과 같은 쪽(아니오)
        fav = {
            "P": max(a["P"], b["P"], key=lambda x: P_RANK.get(x, -1)),
            "문장1": "아니오" if "아니오" in (a["문장1"], b["문장1"]) else a["문장1"],
            "문장2": "아니오" if "아니오" in (a["문장2"], b["문장2"]) else a["문장2"],
            "X": a.get("X") and b.get("X"),
            "기억": "예" if "예" in (a.get("기억"), b.get("기억")) else a.get("기억"),
        }
        for s in ["S1", "S2", "S3", "S4"]:
            fav[s] = bool(a.get(s) or b.get(s))
        per.append(dict(답=aid, job=job, 비교대상=arm, P스크립트=pscript.get(job, ""),
                        채점1=a, 채점2=b, **fav))
    by = defaultdict(list)
    for r in per:
        by[r["비교대상"]].append(r)
    lines = ["| 비교 대상 | n | P2(런 겹침을 번호로 · 문장 ②를 깎음) | P1 이상 | 문장 ② ‘아니오’(엔진과 같음) | 문장 ① ‘아니오’ | X(해로운 답) | S1 · S2 · S3 · S4 | 기억 ‘예’ |",
             "|---|:-:|---|---|---|---|---|---|---|"]
    for arm in sorted(by):
        rs, n = by[arm], len(by[arm])
        lines.append("| {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(
            arm, n, cp(sum(r["P"] == "P2" for r in rs), n), cp(sum(P_RANK.get(r["P"], 0) >= 1 for r in rs), n),
            cp(sum(r["문장2"] == "아니오" for r in rs), n), cp(sum(r["문장1"] == "아니오" for r in rs), n),
            sum(bool(r["X"]) for r in rs),
            " · ".join(str(sum(r[s] for r in rs)) for s in ["S1", "S2", "S3", "S4"]),
            sum(r["기억"] == "예" for r in rs)))
    d_rows = by.get("D", [])
    d_same = sum(r["문장1"] == "아니오" and r["문장2"] == "아니오" for r in d_rows)
    lines += ["", f"- D 가 우리 엔진(① 아니오 · ② 아니오)과 같은 판정: {cp(d_same, len(d_rows))}"]
    for arm in sorted(by):
        vs = {(r["문장1"], r["문장2"]) for r in by[arm]}
        lines.append(f"- {arm} 반복 사이 판정(①, ②) 종류: {len(vs)} — {sorted(vs)}")
    lines += ["", "채점자 일치율: " + " · ".join(f"{it} {v[0]}/{v[1]}" for it, v in agree.items())]
    return per, lines


def agg설계():
    key, s1, s2 = read_key("설계"), scorer_json("설계", 1), scorer_json("설계", 2)
    pm = defaultdict(dict)
    with (HERE / "입력/설계_번호대응.tsv").open(encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            pm[r["job"]][r["보인_번호"]] = r["정답지_이름"]
    catch = defaultdict(lambda: defaultdict(lambda: [0, 0]))     # arm → flaw → [잡음, n]
    inrange = defaultdict(lambda: defaultdict(lambda: [0, 0]))   # arm → plan → [허용 안, n]
    extras = []
    agree = [0, 0]
    for aid, (job, arm) in key.items():
        a = {x["번호"]: x for x in s1[aid]["계획서"]}
        b = {x["번호"]: x for x in s2[aid]["계획서"]}
        for n, plan in pm[job].items():
            xa, xb = a.get(n, {}), b.get(n, {})
            agree[1] += 1
            agree[0] += int(xa.get("판정") == xb.get("판정") and set(xa.get("잡은결함", [])) == set(xb.get("잡은결함", [])))
            got = set(xa.get("잡은결함", [])) | set(xb.get("잡은결함", []))      # AI 쪽에 유리하게: 합집합
            for fl in FLAWS[plan]:
                catch[arm][fl][1] += 1
                catch[arm][fl][0] += int(fl in got)
            ok = (xa.get("판정") in ALLOWED[plan]) or (xb.get("판정") in ALLOWED[plan])
            # 계획서 라의 ‘받칠 수 없음’ 이 타당한 지적에 근거하면 헛경보가 아님(설계서 6절)
            valid_extra = [e for e in (xa.get("정답지밖지적", []) + xb.get("정답지밖지적", [])) if e.get("타당")]
            if plan == "라" and not ok and valid_extra:
                ok = True
            inrange[arm][plan][1] += 1
            inrange[arm][plan][0] += int(ok)
            for e in xa.get("정답지밖지적", []) + xb.get("정답지밖지적", []):
                extras.append([arm, job, plan, str(e.get("타당")), e.get("내용", ""), e.get("이유", "")])
    lines = ["| 비교 대상 | " + " | ".join(f"{fl}" for p in "가나다라마" for fl in FLAWS[p]) + " |",
             "|---|" + "---|" * sum(len(v) for v in FLAWS.values())]
    for arm in sorted(catch):
        lines.append(f"| {arm} 잡음 | " + " | ".join(cp(*catch[arm][fl]) for p in "가나다라마" for fl in FLAWS[p]) + " |")
    lines += ["", "| 비교 대상 | " + " | ".join(f"계획서 {p} 판정이 허용 안" for p in "가나다라마") + " |", "|---|---|---|---|---|---|"]
    for arm in sorted(inrange):
        lines.append(f"| {arm} | " + " | ".join(cp(*inrange[arm][p]) for p in "가나다라마") + " |")
    lines += ["", f"채점자 일치율(판정 + 잡은 결함 모두 같음): {agree[0]}/{agree[1]}"]
    with (SC / "설계_정답지밖지적.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["비교대상", "job", "계획서", "타당", "내용", "이유"])
        w.writerows(extras)
    return extras, lines


if __name__ == "__main__":
    test = sys.argv[1]
    data, lines = agg07() if test == "07" else agg설계()
    (SC / f"{test}_집계.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (SC / f"{test}_집계.json").write_text(json.dumps(data, ensure_ascii=False, indent=1, default=str), encoding="utf-8")
    print("\n".join(lines))
    cur, peak = tracemalloc.get_traced_memory()
    print(f"메모리(Python 할당 기준) 현재 {cur / 1e9:.4f} GB · peak {peak / 1e9:.4f} GB · 소요 시간 {(time.time() - T0) / 3600:.6f} Hour")
