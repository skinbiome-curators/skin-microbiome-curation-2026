"""AI 비교 시험 — 실행 (설계서 v2 · 5절)

실행목록.tsv 의 job 을 저장소 밖 격리 폴더에서 Claude Code 로 따로 돌리고,
stream-json 기록 전체를 원답/ 에 남긴 뒤, 도구 기록으로 격리를 검사한다.

격리 (11/13 권한 시험으로 확인)
- 회차마다 무작위 이름 부모 폴더 아래 새 작업 폴더 · 입력만 복사
- --setting-sources user(허용 목록 · 훅 없음) · --strict-mcp-config(MCP 없음) · dontAsk
- 폴더 밖 Read · cat · head · ls 는 막힘. python 으로 여는 것은 막을 수 없음 → 기록 검사로 거름
- 끝나면 작업 폴더의 파일을 원답/<job>_작업파일/ 로 옮기고 폴더를 지운다

실행: python -I 검증/AI비교시험_261113/스크립트/s2_실행.py <job 접두어…> [--목록 파일] [--동시 3]
예:   ... s2_실행.py 07_A 07_C        → 07_A_1–3 · 07_C_1–3
"""
import argparse
import csv
import json
import re
import secrets
import shutil
import subprocess
import sys
import time
import tracemalloc
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
T0 = time.time()
tracemalloc.start()

HERE = Path(__file__).resolve().parent.parent
OUT = HERE / "원답"
CLAUDE = Path(r"C:\Users\109so\.vscode\extensions\anthropic.claude-code-2.1.322-win32-x64\resources\native-binary\claude.exe")
MODEL = "claude-opus-5-5"
ROOT = Path(r"C:\Users\109so\AppData\Local\Temp\aitest_261113")
ALLOWED = "Write,Bash(python *),Bash(python3 *),Bash(py *)"
DISALLOWED = "WebSearch,WebFetch,Agent,Task,PowerShell"
TIMEOUT_S = 45 * 60

# 폴더 밖을 가리킬 수 있는 흔적(파이썬 코드 · 명령 안)
SUSPECT = [r"\.\.", r"[A-Za-z]:[\\/]", r"(^|[\s'\"(])/(c|Users|home|tmp)/", r"~[\\/]",
           r"os\.walk", r"listdir", r"scandir", r"\.parent", r"dirname\(", r"expanduser", r"getcwd\(\)\s*\.\s*\w", r"chdir"]


def norm(p: str) -> str:
    return p.replace("\\", "/").lower().rstrip("/")


# 11/13 · 1차 뒤: Claude Code 가 세션마다 새로 만드는 전용 임시 폴더(scratchpad)는 작업 폴더와 같게 본다
# (우리 자료 없음 · s2b 감사로 확인 · 설계서 «바꾼 것»)
PAD = re.compile(r"C:[\\/]Users[\\/]109so[\\/]AppData[\\/]Local[\\/]Temp[\\/]claude[\\/][0-9a-f]+[\\/][0-9a-f-]+[\\/]scratchpad[^\s\"'`)]*", re.I)


def scan(jsonl: Path, workdir: Path):
    """도구 호출을 훑어 폴더 밖 접근 흔적을 찾는다 → (무효 사유, 의심 사유, 도구 수, 거부 수)"""
    wd = norm(str(workdir))
    invalid, suspect, n_tool, n_denied = [], [], 0, 0
    for line in jsonl.read_text(encoding="utf-8").splitlines():
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            continue
        msg = d.get("message") or {}
        content = msg.get("content") if isinstance(msg, dict) else None
        if not isinstance(content, list):
            continue
        for c in content:
            if c.get("type") == "tool_result" and c.get("is_error"):
                t = c.get("content")
                t = t if isinstance(t, str) else json.dumps(t, ensure_ascii=False)
                if "denied" in t:
                    n_denied += 1
            if c.get("type") != "tool_use":
                continue
            n_tool += 1
            name, inp = c.get("name"), c.get("input") or {}
            for key in ("file_path", "path", "notebook_path"):
                v = inp.get(key)
                if isinstance(v, str) and v and not PAD.match(v):
                    if (re.match(r"[A-Za-z]:[\\/]|/", v) and not norm(v).startswith(wd)) or ".." in v:
                        invalid.append(f"{name}.{key}={v}")
            texts = [inp.get("command"), inp.get("content"), inp.get("pattern")]
            for t in texts:
                if not isinstance(t, str):
                    continue
                t2 = t.replace(str(workdir), "<WD>").replace(str(workdir).replace("\\", "/"), "<WD>")
                t2 = PAD.sub("<PAD>", t2)
                for pat in SUSPECT:
                    m = re.search(pat, t2, flags=re.M)
                    if m:
                        suspect.append(f"{name}: '{pat}' → …{t2[max(0, m.start()-40):m.end()+40]!r}…")
    return invalid, suspect, n_tool, n_denied


def result_of(jsonl: Path):
    res = {}
    for line in jsonl.read_text(encoding="utf-8").splitlines():
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            continue
        if d.get("type") == "system" and d.get("subtype") == "init":
            res["model_init"] = d.get("model")
            res["memory_paths"] = json.dumps(d.get("memory_paths"), ensure_ascii=False)
            res["mcp"] = json.dumps(d.get("mcp_servers"), ensure_ascii=False)
        if d.get("type") == "result":
            res.update(result=d.get("result", ""), is_error=d.get("is_error"), num_turns=d.get("num_turns"),
                       duration_ms=d.get("duration_ms"), cost_usd=d.get("total_cost_usd"),
                       models=",".join((d.get("modelUsage") or {}).keys()))
    return res


def run_job(job, prompt_rel, files):
    OUT.mkdir(exist_ok=True)
    if (OUT / f"{job}.jsonl").exists():
        return job, "건너뜀(이미 있음)", {}
    parent = ROOT / secrets.token_hex(8)
    wd = parent / "w"
    wd.mkdir(parents=True)
    for spec in files.split(";"):
        src, dst = spec.split(">")
        shutil.copyfile(HERE / src, wd / dst)
    prompt = (HERE / prompt_rel).read_text(encoding="utf-8")
    (OUT / f"{job}_보낸프롬프트.txt").write_text(prompt, encoding="utf-8")
    cmd = [str(CLAUDE), "-p", prompt, "--model", MODEL, "--output-format", "stream-json", "--verbose",
           "--setting-sources", "user", "--strict-mcp-config", "--permission-mode", "dontAsk",
           "--allowedTools", ALLOWED, "--disallowedTools", DISALLOWED]
    t = time.time()
    with (OUT / f"{job}.jsonl").open("wb") as fo, (OUT / f"{job}.err").open("wb") as fe:
        try:
            subprocess.run(cmd, cwd=wd, stdout=fo, stderr=fe, timeout=TIMEOUT_S)
            status = "끝"
        except subprocess.TimeoutExpired:
            status = "시간 초과"
    inv, sus, n_tool, n_den = scan(OUT / f"{job}.jsonl", wd)
    # 작업 파일 보관(입력 복사본 제외) · 폴더 지우기
    keep = OUT / f"{job}_작업파일"
    inputs = {s.split(">")[1] for s in files.split(";")}
    for p in wd.rglob("*"):
        if p.is_file() and p.name not in inputs:
            dst = keep / p.relative_to(wd)
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(p, dst)
    shutil.rmtree(parent, ignore_errors=True)
    info = result_of(OUT / f"{job}.jsonl")
    info.update(status=status, minutes=round((time.time() - t) / 60, 1), n_tool=n_tool, n_denied=n_den,
                invalid=" | ".join(inv), suspect=" | ".join(sus[:12]), workdir=str(wd))
    (OUT / f"{job}_검사.json").write_text(json.dumps(info, ensure_ascii=False, indent=1), encoding="utf-8")
    verdict = "무효" if inv else ("의심(사람 확인)" if sus else "유효")
    return job, f"{status} · {verdict}", info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefixes", nargs="+")
    ap.add_argument("--목록", default=str(HERE / "실행목록.tsv"))
    ap.add_argument("--동시", type=int, default=3)
    ap.add_argument("--출력", default=str(HERE / "원답"))  # 실행기 검수 때만 바꿈
    a = ap.parse_args()
    global OUT
    OUT = Path(a.출력)
    with open(a.목록, encoding="utf-8", newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))[1:]
    jobs = [r for r in rows if any(r[0].startswith(p) for p in a.prefixes)]
    print(f"실행 {len(jobs)}개 · 동시 {a.동시} · 모델 {MODEL}", flush=True)
    done = 0
    with ThreadPoolExecutor(max_workers=a.동시) as ex:
        futs = {ex.submit(run_job, *r): r[0] for r in jobs}
        for fu in as_completed(futs):
            done += 1
            try:
                job, st, info = fu.result()
                extra = f" · {info.get('minutes')}분 · 도구 {info.get('n_tool')} · 거부 {info.get('n_denied')} · 모델 {info.get('models')}" if info else ""
                print(f"[{done}/{len(jobs)}] {job}: {st}{extra}", flush=True)
            except Exception as e:  # 실행기 자체 오류는 기록하고 계속
                print(f"[{done}/{len(jobs)}] {futs[fu]}: 실행기 오류 {e!r}", flush=True)
    cur, peak = tracemalloc.get_traced_memory()
    print(f"메모리(Python 할당 기준) 현재 {cur / 1e9:.4f} GB · peak {peak / 1e9:.4f} GB · 소요 시간 {(time.time() - T0) / 3600:.4f} Hour", flush=True)


if __name__ == "__main__":
    main()
