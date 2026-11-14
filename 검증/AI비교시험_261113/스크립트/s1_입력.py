"""AI 비교 시험 — 입력 만들기 (설계서 v2 · 3절 · 6절)

하는 일
1. 07 자료에서 출처 식별자를 지운 가공본을 만든다
   - sample_info: run(ENA) · project · *_paper 세 열을 지우고, FASTQ 첫 줄 열(fastq_header_line1)을 더한다
   - genus_pct: 첫 열 ENA 번호 → 시료 ID(sample_id)
2. 설계 시험의 계획서를 회차마다 시드 고정으로 섞어 '계획서 1–5' 한 파일로 만든다
3. 실행 목록(실행목록.tsv)을 만든다
4. 스스로 검사한다(식별자 남음 · 행 수 · 열 · 머리줄)

실행: python -I 스크립트/s1_입력.py   (저장소 새주제_v11 안 어디서든)
"""
import csv
import hashlib
import random
import re
import sys
import time
import tracemalloc
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
T0 = time.time()
tracemalloc.start()

HERE = Path(__file__).resolve().parent.parent          # 검증/AI비교시험_261113
REPO = HERE.parent.parent                              # 새주제_v11
IN = HERE / "입력"
SRC_INFO = REPO / "data/07_skincare_16S_sample_info.tsv"
SRC_GENUS = REPO / "data/07_skincare_16S_genus_pct.tsv"
SRC_HEAD = REPO / "results/효능판정_시범_런번호.tsv"
PLANS = IN / "계획서_원본.md"
DROP = ["run", "project", "reads_paper", "observed_otus_paper", "shannon_paper"]
SEED = 20261113
REPS = 3


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read_tsv(p: Path):
    with p.open(encoding="utf-8", newline="") as f:
        r = csv.reader(f, delimiter="\t")
        rows = list(r)
    return rows[0], rows[1:]


def write_tsv(p: Path, header, rows):
    with p.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


# ── 1. 07 가공본 ──────────────────────────────────────────────
h_info, info = read_tsv(SRC_INFO)
h_gen, gen = read_tsv(SRC_GENUS)
h_head, head = read_tsv(SRC_HEAD)

i_run, i_sid = h_info.index("run"), h_info.index("sample_id")
ena2sid = {r[i_run]: r[i_sid] for r in info}
assert len(ena2sid) == 75 and len(set(ena2sid.values())) == 75, "07 시료 75개 · 시료 ID 중복 없음이어야 함"

# FASTQ 첫 줄: 07 행만 · 앞의 '@ERR….1' 을 '@<시료 ID>.1' 로
hc, hr, hf = h_head.index("cohort"), h_head.index("run"), h_head.index("first_header")
first = {}
for r in head:
    if r[hc] != "07":
        continue
    ena = r[hr]
    line = r[hf]
    tok, rest = line.split(" ", 1)
    m = re.fullmatch(r"@(ERR\d+)\.(\d+)", tok)
    assert m and m.group(1) == ena, f"머리줄 앞부분이 예상과 다름: {tok}"
    first[ena] = f"@{ena2sid[ena]}.{m.group(2)} {rest}"
assert set(first) == set(ena2sid), "07 시료마다 머리줄이 하나씩 있어야 함"

keep = [i for i, c in enumerate(h_info) if c not in DROP]
new_h_info = [h_info[i] for i in keep] + ["fastq_header_line1"]
new_info = [[r[i] for i in keep] + [first[r[i_run]]] for r in info]

assert h_gen[0] == "run"
new_h_gen = ["sample_id"] + h_gen[1:]
new_gen = [[ena2sid[r[0]]] + r[1:] for r in gen]

out07 = IN / "07"
out07.mkdir(parents=True, exist_ok=True)
write_tsv(out07 / "sample_info.tsv", new_h_info, new_info)
write_tsv(out07 / "genus_pct.tsv", new_h_gen, new_gen)
write_tsv(IN / "대응표_ENA_시료ID.tsv", ["ena_run", "sample_id"], sorted(ena2sid.items()))

# ── 2. 설계 계획서 섞기 ──────────────────────────────────────
text = PLANS.read_text(encoding="utf-8")
blocks = re.split(r"^## (가|나|다|라|마)\s*$", text, flags=re.M)
plans = {blocks[i]: blocks[i + 1].strip() for i in range(1, len(blocks), 2)}
assert list(plans) == ["가", "나", "다", "라", "마"], "계획서 원본은 가–마 다섯이어야 함"
lens = {k: len(v.splitlines()) for k, v in plans.items()}
assert len(set(lens.values())) == 1, f"계획서마다 항목 수가 같아야 함: {lens}"

outd = IN / "설계"
outd.mkdir(parents=True, exist_ok=True)
key_rows = []
jobs_design = []
for arm in ["A", "C"]:
    for rep in range(1, REPS + 1):
        job = f"설계_{arm}_{rep}"
        order = list(plans)
        random.Random(f"{SEED}-{job}").shuffle(order)
        parts = ["# 계획서 다섯\n"]
        for n, k in enumerate(order, 1):
            parts.append(f"## 계획서 {n}\n{plans[k]}\n")
            key_rows.append([job, str(n), k])
        (outd / f"{job}_계획서.md").write_text("\n".join(parts), encoding="utf-8")
        jobs_design.append([job, f"프롬프트/설계_{arm}.txt", f"입력/설계/{job}_계획서.md>계획서.md"])
write_tsv(IN / "설계_번호대응.tsv", ["job", "보인_번호", "정답지_이름"], key_rows)

# ── 3. 실행 목록 (07 9개 → 설계 6개) ──────────────────────────
jobs07 = []
for arm in ["A", "C", "D"]:
    for rep in range(1, REPS + 1):
        files = ["입력/07/genus_pct.tsv>genus_pct.tsv", "입력/07/sample_info.tsv>sample_info.tsv"]
        if arm == "D":
            files.append("입력/판정규칙.md>판정규칙.md")
        jobs07.append([f"07_{arm}_{rep}", f"프롬프트/07_{arm}.txt", ";".join(files)])
jobs = jobs07 + jobs_design
write_tsv(HERE / "실행목록.tsv", ["job", "prompt", "files(원본>폴더 안 이름)"], jobs)

# ── 4. 스스로 검사 ────────────────────────────────────────────
problems = []
for p in [out07 / "sample_info.tsv", out07 / "genus_pct.tsv"] + sorted(outd.glob("*.md")):
    t = p.read_text(encoding="utf-8")
    for bad in ["ERR59", "PRJEB", "_paper"]:
        if bad in t:
            problems.append(f"{p.name}: '{bad}' 남음")
hi, ri = read_tsv(out07 / "sample_info.tsv")
hg, rg = read_tsv(out07 / "genus_pct.tsv")
if len(ri) != 75 or len(rg) != 75:
    problems.append("행 수가 75가 아님")
if len(hi) != len(h_info) - len(DROP) + 1:
    problems.append("sample_info 열 수가 맞지 않음")
if any(not r[-1].startswith("@P") or ":24" not in r[-1] for r in ri):
    problems.append("머리줄 열 형식이 이상함")
# 머리줄의 기기:런:플로셀 부분이 원본과 그대로인지
orig_mid = {r[hr]: r[hf].split(" ", 1)[1] for r in head if r[hc] == "07"}
sid2ena = {v: k for k, v in ena2sid.items()}
for r in ri:
    sid = r[hi.index("sample_id")]
    if r[-1].split(" ", 1)[1] != orig_mid[sid2ena[sid]]:
        problems.append(f"{sid}: 머리줄 뒷부분이 원본과 다름")
# genus 값이 원본과 같은지(시료 ID 만 바뀌었는지)
orig_gen = {r[0]: r[1:] for r in gen}
for r in rg:
    if r[1:] != orig_gen[sid2ena[r[0]]]:
        problems.append(f"{r[0]}: 속 비율 값이 원본과 다름")

# SHA
sha_rows = [[str(p.relative_to(HERE)).replace("\\", "/"), sha(p)] for p in
            [out07 / "sample_info.tsv", out07 / "genus_pct.tsv", IN / "판정규칙.md", PLANS] + sorted(outd.glob("*.md"))]
sha_rows += [[f"(원본) {p.relative_to(REPO)}".replace("\\", "/"), sha(p)] for p in [SRC_INFO, SRC_GENUS, SRC_HEAD]]
write_tsv(IN / "SHA256.tsv", ["파일", "sha256"], sha_rows)

cur, peak = tracemalloc.get_traced_memory()
print(f"07 가공본: 시료 {len(ri)} · sample_info 열 {len(h_info)} → {len(hi)} · genus 열 {len(hg)}")
print(f"설계 회차 {len(list(outd.glob('*.md')))}개 · 실행 목록 {len(jobs)}개")
print("검사:", "통과" if not problems else f"문제 {len(problems)}개")
for x in problems:
    print("  ✗", x)
print(f"메모리(Python 할당 기준) 현재 {cur / 1e9:.4f} GB · peak {peak / 1e9:.4f} GB · 소요 시간 {(time.time() - T0) / 3600:.6f} Hour")
sys.exit(1 if problems else 0)
