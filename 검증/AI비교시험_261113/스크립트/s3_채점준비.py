"""AI 비교 시험 — 채점 준비 (설계서 v2 · 5절 무효 · 7절)

하는 일
1. 원답/<job>.jsonl 에서 마지막 답(result)을 꺼낸다
2. 무효를 정한다 — 실행 검사(_검사.json)의 invalid · 사람 확인표(원답/검사_사람확인.tsv) · 우리 용어
3. 07: 런 번호 · 플로셀 ID 를 찾아 P2 후보를 1차 표시(스크립트)
4. 비교 대상 이름을 지우고 시드 고정으로 섞은 답모음 + 숨긴 키를 만든다
5. 채점 실행 목록(채점자 2명)을 만든다

실행: python -I 검증/AI비교시험_261113/스크립트/s3_채점준비.py 07|설계
"""
import csv
import json
import random
import re
import sys
import time
import tracemalloc
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
T0 = time.time()
tracemalloc.start()

HERE = Path(__file__).resolve().parent.parent
RAW = HERE / "원답"
SC = HERE / "채점"
SEED = "20261113-채점"
TERMS_ALL = ["관문", "E0", "결정표", "AI비교시험", "새주제_v11", "판정 엔진"]
TERMS_D_OK = {"E0", "결정표"}          # D 는 규칙 문서에 이 말이 있어 허용
P_PAT = re.compile(r"(?<!\d)24[1-4](?!\d)|CW4R9|CW4RL|J2766|J276W")

# 설계 정답지(설계서 6절) — 채점자에게 답마다 붙여 준다
KEY = {
    "가": ("받칠 수 없음", [("가1", "군과 시퀀싱 런 · 추출 플레이트가 겹침(제품군 = 런 1 · 무처리군 = 런 2)")]),
    "나": ("받칠 수 없음", [("나1", "무처리 대조가 없음"), ("나2", "사용 전 측정이 1회뿐"),
                          ("나3", "주 균 · 방향을 결과를 본 뒤 정함(사전 등록 없음)")]),
    "다": ("받칠 수 없음 · 조건부(문장을 ‘포도상구균 속 비율’로 바꾸는 조건)도 허용",
           [("다1", "비율 자료로 ‘양’을 주장(절대량 측정 없음)"),
            ("다2", "16S V3–V4 로는 포도상구균 종(S. epidermidis)을 가르기 어려움")]),
    "라": ("받칠 수 있음 · 조건부 모두 허용(‘받칠 수 없음’은 헛경보 — 정답지밖지적이 타당하면 예외)", []),
    "마": ("받칠 수 없음 · 조건부도 허용",
           [("마1", "채취 부위가 시점마다 다름(사용 전 볼 → 4주 이마)"),
            ("마2", "함께 쓰는 세안제 · 화장품을 통제 · 기록하지 않음")]),
}


def final_text(p: Path) -> str:
    out = ""
    for line in p.read_text(encoding="utf-8").splitlines():
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            continue
        if d.get("type") == "result":
            out = d.get("result") or ""
    return out


def human_checks():
    p = RAW / "검사_사람확인.tsv"
    if not p.exists():
        return {}
    with p.open(encoding="utf-8", newline="") as f:
        return {r["job"]: (r["판정"], r["이유"]) for r in csv.DictReader(f, delimiter="\t")}


def main(test: str):
    SC.mkdir(exist_ok=True)
    hc = human_checks()
    jobs = sorted(p.stem for p in RAW.glob(f"{test}_*.jsonl"))
    rows, answers, problems = [], [], []
    for job in jobs:
        arm = job.split("_")[1]
        chk = json.loads((RAW / f"{job}_검사.json").read_text(encoding="utf-8"))
        text = final_text(RAW / f"{job}.jsonl")
        why = []
        human_ok = job in hc and hc[job][0] == "유효"   # 11/13 · 사람 감사(s2b)가 자동 검사를 덮어씀 — 설계서 «바꾼 것»
        if chk.get("invalid") and not human_ok:
            why.append("폴더 밖 접근: " + chk["invalid"][:200])
        if chk.get("suspect") and job not in hc:
            problems.append(f"{job}: 의심 표시가 있으나 사람 확인표에 없음 — 확인 뒤 다시 실행")
        if job in hc and hc[job][0] == "무효":
            why.append("사람 확인: " + hc[job][1])
        if chk.get("is_error") or not text.strip():
            why.append("오류 · 빈 답")
        used = [t for t in TERMS_ALL if t in text and not (arm == "D" and t in TERMS_D_OK)]
        if used:
            why.append("우리 용어: " + ",".join(used))
        p_script = bool(P_PAT.search(text)) if test == "07" else None
        rows.append([job, arm, "무효" if why else "유효", " / ".join(why), "" if p_script is None else str(p_script)])
        if not why:
            answers.append((job, arm, text))
    if problems:
        for x in problems:
            print("✗", x)
        sys.exit(1)

    rng = random.Random(f"{SEED}-{test}")
    rng.shuffle(answers)
    ids = [f"R{rng.randrange(1000, 9999)}" for _ in answers]
    assert len(set(ids)) == len(ids)
    key_rows, parts = [], [f"# 답모음 ({len(answers)}개)\n"]
    plan_map = {}
    if test == "설계":
        with (HERE / "입력/설계_번호대응.tsv").open(encoding="utf-8", newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                plan_map.setdefault(r["job"], {})[r["보인_번호"]] = r["정답지_이름"]
    for aid, (job, arm, text) in zip(ids, answers):
        clean = text.replace("판정규칙.md", "[규칙 문서]")
        parts.append(f"\n---\n\n## 답 {aid}\n")
        if test == "설계":
            plan_text = (HERE / f"입력/설계/{job}_계획서.md").read_text(encoding="utf-8")
            parts.append("### 이 답이 본 계획서\n\n" + plan_text.replace("# 계획서 다섯", "").strip() + "\n")
            parts.append("### 정답지(이 답의 번호 기준)\n")
            for n in sorted(plan_map[job]):
                exp, flaws = KEY[plan_map[job][n]]
                fl = "; ".join(f"{i} {d}" for i, d in flaws) or "심은 결함 없음"
                parts.append(f"- 계획서 {n}: 기대 판정 {exp} · 결함 {fl}")
            parts.append("\n### 답\n")
        parts.append(clean + "\n")
        key_rows.append([aid, job, arm])
    (SC / f"{test}_답모음.md").write_text("\n".join(parts), encoding="utf-8")
    with (SC / f"{test}_키.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["답", "job", "비교대상"])
        w.writerows(key_rows)
    with (SC / f"{test}_유효성.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["job", "비교대상", "판정", "이유", "P스크립트"])
        w.writerows(rows)
    # 채점 실행 목록(채점자 2명 · 같은 답모음)
    lst = SC / f"{test}_채점목록.tsv"
    with lst.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["job", "prompt", "files"])
        for k in (1, 2):
            w.writerow([f"채점_{test}_{k}", f"채점/프롬프트_{test}.txt", f"채점/{test}_답모음.md>답모음.md"])
    cur, peak = tracemalloc.get_traced_memory()
    print(f"{test}: 답 {len(jobs)}개 · 유효 {len(answers)} · 무효 {len(jobs) - len(answers)} → 채점/{test}_답모음.md · 키 · 채점목록")
    print(f"메모리(Python 할당 기준) 현재 {cur / 1e9:.4f} GB · peak {peak / 1e9:.4f} GB · 소요 시간 {(time.time() - T0) / 3600:.6f} Hour")


if __name__ == "__main__":
    main(sys.argv[1])
