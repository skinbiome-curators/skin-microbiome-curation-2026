"""D1 — 기능성화장품 보고품목 집계(10/17 규칙과 같게) + 다른 방법으로 다시 세기.

규칙(10/17 `references/수요추정_원문확인_20261017.md` 2절과 같음):
- 품목명(ITEM_NAME)에 검색어가 그대로 들어간 것만 센다 — 마이크로바이옴 · 프로바이오틱 · 유산균
- 제품 1개 = 보고일련번호(COSMETIC_REPORT_SEQ) 1개 · 회사 = 업체일련번호(ENTP_SEQ) 고유 수 · 연도 = 보고일(REPORT_DATE) 앞 4자리
- 취소(CANCEL_APPROVAL_YN = 'Y')는 따로 센다
더 보는 것:
- 연도별 전체 기능성 보고 건수 · 보고 회사 수(분모) → 비율
- 상위 5개 회사 이름 · 해마다 다시 보고하는 회사 · 제품 1개뿐인 회사
- 보조 검색어(바이옴 · biome · probiotic · 포스트바이오틱 · 락토 · 비피다)로 빠진 것이 얼마인지 — 주 규칙에는 넣지 않음
다시 세기(독립 방법):
- ① API 자체 품목명 검색(item_name=검색어)의 totalCount 와 내려받은 전체에서 센 수를 맞춘다
- ② 10/17 숫자(372 제품 · 137 회사 · 2026-10-15 까지)와 같은 기간으로 맞춘다
결과: 결과_d1/ (저장소 안) · 원자료는 micro3/data/업계데이터_261116/d1/

사용: python -I s6_d1_집계.py
"""
import sys, json, gzip, glob, re, csv, time, tracemalloc, urllib.request, urllib.parse
from collections import Counter, defaultdict
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
tracemalloc.start(); T0 = time.time()
HERE = Path(__file__).resolve().parent
MICRO3 = HERE.parents[2]
RAW = MICRO3 / "data" / "업계데이터_261116" / "d1"
OUT = HERE / "결과_d1"; OUT.mkdir(exist_ok=True)
KW = ["마이크로바이옴", "프로바이오틱", "유산균"]
KW_EXTRA = re.compile(r"바이옴|biome|probiotic|포스트바이오틱|락토|비피다", re.I)

# --- 읽기
rows = []
for f in sorted(glob.glob(str(RAW / "page_*.json.gz"))):
    j = json.loads(gzip.decompress(Path(f).read_bytes()))
    it = j["body"].get("items") or []
    if isinstance(it, dict): it = it.get("item", [])
    rows += it
meta = json.loads((RAW / "받기_기록.json").read_text(encoding="utf-8"))
seqs = Counter(r["COSMETIC_REPORT_SEQ"] for r in rows)
dup = sum(c - 1 for c in seqs.values() if c > 1)
uniq = {}
for r in rows:
    uniq.setdefault(r["COSMETIC_REPORT_SEQ"], r)
rows = list(uniq.values())
print(f"받은 행 {sum(seqs.values()):,} · 중복 {dup} · 고유 {len(rows):,} · API totalCount {meta['totalCount']:,}")

yr = lambda r: (r.get("REPORT_DATE") or "")[:4]
canc = lambda r: (r.get("CANCEL_APPROVAL_YN") or "") == "Y"
name = lambda r: r.get("ITEM_NAME") or ""
hit = lambda r: any(k in name(r) for k in KW)

# --- 전체(분모)
all_y = defaultdict(lambda: [0, set(), 0])  # 제품, 회사, 취소
for r in rows:
    a = all_y[yr(r)]; a[0] += 1; a[1].add(r["ENTP_SEQ"]); a[2] += canc(r)

# --- 검색어 집합
sel = [r for r in rows if hit(r)]
first_year = {}
for r in sorted(sel, key=lambda r: r.get("REPORT_DATE") or ""):
    first_year.setdefault(r["ENTP_SEQ"], yr(r))
sel_y = defaultdict(lambda: [0, set(), 0])
for r in sel:
    a = sel_y[yr(r)]; a[0] += 1; a[1].add(r["ENTP_SEQ"]); a[2] += canc(r)
new_y = Counter(first_year.values())

years = [str(y) for y in range(2017, 2027)]
with open(OUT / "연도별.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["연도", "검색어 제품", "검색어 보고 회사", "처음 보고한 회사", "검색어 취소", "전체 기능성 제품", "전체 보고 회사", "제품 비율(%)", "회사 비율(%)"])
    for y in years:
        s, a = sel_y[y], all_y[y]
        w.writerow([y, s[0], len(s[1]), new_y.get(y, 0), s[2], a[0], len(a[1]),
                    round(100 * s[0] / a[0], 3) if a[0] else "", round(100 * len(s[1]) / len(a[1]), 2) if a[1] else ""])
print("\n연도 | 검색어 제품 · 회사 · 첫 보고 | 전체 제품 · 회사 | 제품 % · 회사 %")
for row in csv.reader(open(OUT / "연도별.tsv", encoding="utf-8"), delimiter="\t"):
    print(" | ".join(row))

# --- 키워드별 · 10/17 대조
def window(rs, end="20261015"):
    return [r for r in rs if (r.get("REPORT_DATE") or "") <= end]
print("\n검색어별(전체 기간) · (2026-10-15 까지)")
kw_rows = []
for k in KW + ["셋 중 하나 이상"]:
    rs = sel if k.startswith("셋") else [r for r in rows if k in name(r)]
    wv = window(rs)
    kw_rows.append([k, len(rs), len({r['ENTP_SEQ'] for r in rs}), len(wv), len({r['ENTP_SEQ'] for r in wv}), sum(canc(r) for r in rs)])
    print(f"  {k}: 제품 {len(rs)} · 회사 {kw_rows[-1][2]} | 10/15까지 제품 {len(wv)} · 회사 {kw_rows[-1][4]} | 취소 {kw_rows[-1][5]}")
since24 = [r for r in sel if yr(r) >= "2024"]
print(f"  2024년 이후: 제품 {len(since24)} · 회사 {len({r['ENTP_SEQ'] for r in since24})}")
print("  10/17 기록: 마이크로바이옴 143 · 55 / 프로바이오틱 214 · 75 / 유산균 15 · 10 / 셋 372 · 137 / 2024 이후 139 · 67")

# --- 회사 분포 · 상위 5 · 재보고
by_co = defaultdict(list)
for r in sel: by_co[r["ENTP_SEQ"]].append(r)
co_rows = []
for c, rs in by_co.items():
    names = Counter(r["ENTP_NAME"] for r in rs).most_common(1)[0][0]
    ys = sorted({yr(r) for r in rs})
    co_rows.append([names, c, len(rs), len(ys), ys[0], ys[-1], " ".join(ys)])
co_rows.sort(key=lambda x: (-x[2], x[0]))
with open(OUT / "회사별.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["회사(가장 많이 쓴 이름)", "업체일련번호", "검색어 제품 수", "보고한 해 수", "처음", "마지막", "보고한 해"])
    w.writerows(co_rows)
top5 = co_rows[:5]
one = sum(1 for x in co_rows if x[2] == 1)
multi_year = sum(1 for x in co_rows if x[3] >= 2)
print(f"\n회사 {len(co_rows)}곳 · 제품 1개뿐 {one}곳 · 2개 해 이상 보고 {multi_year}곳")
print(f"상위 5개 회사가 {sum(x[2] for x in top5)}/{len(sel)} 제품:")
for x in top5: print(f"  {x[0]} · {x[2]}개 · {x[4]}–{x[5]} · {x[3]}개 해")

# --- 보조 검색어로 빠진 것
extra = [r for r in rows if not hit(r) and KW_EXTRA.search(name(r))]
ec = Counter(m.group(0).lower() for r in extra for m in [KW_EXTRA.search(name(r))])
print(f"\n보조 검색어에만 걸리는 제품 {len(extra)}개 · {dict(ec)}")
with open(OUT / "보조검색어_제품.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["보고일", "회사", "품목명"])
    for r in sorted(extra, key=lambda r: r.get("REPORT_DATE") or ""):
        w.writerow([r.get("REPORT_DATE"), r.get("ENTP_NAME"), name(r)])

# --- 넓힌 규칙(11/16 추가): 품목명에 '바이옴'(마이크로바이옴 · 스킨바이옴 · ○○바이옴 등) 또는 주 검색어
#     '락토'는 대부분 갈락토미세스(효모 발효물)라 넣지 않는다 · '비피다'(발효 용해물)도 넣지 않는다
wide = [r for r in rows if hit(r) or "바이옴" in name(r)]
wide_first = {}
for r in sorted(wide, key=lambda r: r.get("REPORT_DATE") or ""):
    wide_first.setdefault(r["ENTP_SEQ"], yr(r))
wy = defaultdict(lambda: [0, set()])
for r in wide:
    a = wy[yr(r)]; a[0] += 1; a[1].add(r["ENTP_SEQ"])
wnew = Counter(wide_first.values())
wco = Counter()
wname = {}
for r in wide:
    wco[r["ENTP_SEQ"]] += 1; wname.setdefault(r["ENTP_SEQ"], r["ENTP_NAME"])
print(f"\n넓힌 규칙(주 검색어 + 바이옴): 제품 {len(wide)} · 회사 {len(wco)} · 취소 {sum(canc(r) for r in wide)}")
with open(OUT / "연도별_넓힌규칙.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["연도", "제품", "보고 회사", "처음 보고한 회사", "전체 기능성 제품", "전체 보고 회사", "제품 비율(%)", "회사 비율(%)"])
    for y in years:
        a, t = wy[y], all_y[y]
        w.writerow([y, a[0], len(a[1]), wnew.get(y, 0), t[0], len(t[1]),
                    round(100 * a[0] / t[0], 3) if t[0] else "", round(100 * len(a[1]) / len(t[1]), 2) if t[1] else ""])
        print(f"  {y}: 제품 {a[0]} · 회사 {len(a[1])} · 처음 {wnew.get(y, 0)} · 회사 비율 {100*len(a[1])/len(t[1]):.2f}%")
top = wco.most_common(5)
print(f"  상위 5개 회사 {sum(n for _, n in top)}/{len(wide)}: " + " · ".join(f"{wname[c]} {n}" for c, n in top))
print(f"  제품 1개뿐 {sum(1 for n in wco.values() if n == 1)}곳 · 2024년 이후 회사 {len({r['ENTP_SEQ'] for r in wide if yr(r) >= '2024'})}곳")

# --- 다시 세기 ① API 자체 검색
key = (MICRO3 / ".secrets" / "data_go_kr_key.txt").read_text(encoding="utf-8-sig").strip()
EP = "https://apis.data.go.kr/1471000/FtnltCosmRptPrdlstInfoService01/getRptPrdlstInq01"
print("\n① API 자체 품목명 검색 vs 내려받은 전체")
chk = []
for k in KW + ["바이옴"]:
    url = f"{EP}?serviceKey={key}&pageNo=1&numOfRows=1&type=json&item_name={urllib.parse.quote(k)}"
    try:
        j = json.loads(urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"}), timeout=120).read())
        api_n = int(j["body"]["totalCount"])
    except Exception as e:
        api_n = None; print("  API 실패", str(e).replace(key, "***")[:80])
    local_n = sum(1 for r in rows if k in name(r))
    chk.append([k, api_n, local_n, "OK" if api_n == local_n else "다름"])
    print(f"  {k}: API {api_n} · 내려받은 전체 {local_n} → {chk[-1][3]}")
with open(OUT / "대조.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["검색어", "API totalCount", "내려받은 전체에서 셈", "결과"]); w.writerows(chk)
    w.writerow([]); w.writerow(["검색어", "전체 기간 제품", "전체 기간 회사", "10/15까지 제품", "10/15까지 회사", "취소"]); w.writerows(kw_rows)

cur, peak = tracemalloc.get_traced_memory()
print(f"\n→ {OUT}")
print(f"peak {peak/1e9:.3f} GB · {(time.time()-T0)/3600:.4f} Hour")
