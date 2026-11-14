"""AI 비교 시험 — 격리 감사 (s2 의 자동 검사가 '무효 · 의심'을 낸 회차를 사람이 보기 쉽게 풀어 씀)

회차마다 도구 호출의 경로를 셋으로 나눈다
- 작업 폴더(입력이 있는 곳)
- 그 세션 전용 임시 폴더(Claude Code 가 세션마다 새로 만드는 scratchpad — 우리 자료 없음)
- 그 밖 → 하나라도 있으면 무효
그리고 우리 세션 · 저장소 · 메모리 흔적, 거부된 호출을 적는다.

실행: python -I 검증/AI비교시험_261113/스크립트/s2b_감사.py <job 접두어> [--출력 원답 폴더]
"""
import argparse
import glob
import json
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent.parent
OURS = ["928799ca", "fe8fc9b252b6b757", "onedrive", "새주제", "micro3", ".claude/projects", "/memory"]
PAD = re.compile(r"C:[\\/]Users[\\/]109so[\\/]AppData[\\/]Local[\\/]Temp[\\/]claude[\\/][0-9a-f]+[\\/][0-9a-f-]+[\\/]scratchpad", re.I)
ABS = re.compile(r"[A-Za-z]:[\\/][^\s\"'`)\]]+")


def n(p: str) -> str:
    return p.replace("\\", "/").lower()


def audit(jsonl: Path, outdir: Path):
    job = jsonl.stem
    chk = json.loads((outdir / f"{job}_검사.json").read_text(encoding="utf-8"))
    wd = n(chk["workdir"])
    paths, texts, denied = [], [], []
    for line in jsonl.read_text(encoding="utf-8").splitlines():
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            continue
        m = d.get("message") or {}
        c = m.get("content") if isinstance(m, dict) else None
        if not isinstance(c, list):
            continue
        for x in c:
            if x.get("type") == "tool_use":
                inp = x.get("input") or {}
                for k in ("file_path", "path"):
                    if isinstance(inp.get(k), str):
                        paths.append((x["name"], inp[k]))
                for k in ("command", "content"):
                    if isinstance(inp.get(k), str):
                        texts.append((x["name"], inp[k]))
            if x.get("type") == "tool_result" and x.get("is_error"):
                t = x.get("content")
                t = t if isinstance(t, str) else json.dumps(t, ensure_ascii=False)
                denied.append(t[:100].replace("\n", " "))
    blob = " ".join(v for _, v in texts + paths)
    pads = {n(p) for p in PAD.findall(blob)}

    def ok(p):
        q = n(p)
        return q.startswith(wd) or any(q.startswith(x) for x in pads)

    other = [f"{a}:{v}" for a, v in paths if re.match(r"[A-Za-z]:[\\/]|/", v) and not ok(v)]
    other += [f"(명령·스크립트 안) {a}" for _, v in texts for a in ABS.findall(v) if not ok(a)]
    other += [f"(상위 폴더 '..') {a}: …{v[max(0, v.find('..') - 30):v.find('..') + 30]}…" for a, v in texts + paths if ".." in v and not re.search(r"\.\.\.", v)]
    hits = sorted({o for o in OURS if o in n(blob)})
    return job, len(pads), sorted(set(other)), hits, denied


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix")
    ap.add_argument("--출력", default=str(HERE / "원답"))
    a = ap.parse_args()
    outdir = Path(a.출력)
    for p in sorted(outdir.glob(f"{a.prefix}*.jsonl")):
        job, npad, other, hits, denied = audit(p, outdir)
        print(f"== {job}: 전용 임시 폴더 {npad}개 · 그 밖 경로 {other or '없음'} · 우리 흔적 {hits or '없음'}")
        for dn in denied:
            print("   거부 · 오류:", dn)


if __name__ == "__main__":
    main()
