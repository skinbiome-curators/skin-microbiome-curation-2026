"""D3 — 받은 사업보고서에서 수치 · 문장을 꺼낸다.

1) 주요 계정 API → 매출액 · 영업이익 · 당기순이익(당기)
2) 직원 현황 API → 행 그대로 + 합계
3) 보고서 본문(XML) → 텍스트로 바꿔 장(TITLE)별로 나눔
   - 장 '사업의 개요' · '주요 제품 및 서비스' · '매출 및 수주상황' 전문을 텍스트 파일로(원자료 폴더)
   - 본문 전체에서 검색어가 든 줄을 뽑음(마이크로바이옴 · 미생물 · 인체적용시험 · 고객 · 수출 · 단가 · 가격)
결과: 저장소 안 결과_d3/ (표 · 검색어 줄) · 장 전문은 원자료 폴더 텍스트/

사용: python -I s2_dart_추출.py
"""
import sys, json, re, html, time, tracemalloc, zipfile, csv
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
tracemalloc.start(); T0 = time.time()
HERE = Path(__file__).resolve().parent
MICRO3 = HERE.parents[2]
RAW = MICRO3 / "data" / "업계데이터_261116" / "dart"
TXT = RAW.parent / "텍스트"; TXT.mkdir(exist_ok=True)
OUT = HERE / "결과_d3"; OUT.mkdir(exist_ok=True)
LOG = json.loads((RAW / "받기_기록.json").read_text(encoding="utf-8"))

KW = ["마이크로바이옴", "미생물", "인체적용시험", "인체 적용시험", "고객", "수출", "단가", "가격"]
SECTIONS = ["사업의 개요", "주요 제품 및 서비스", "매출 및 수주상황"]


def xml_to_text(s: str) -> str:
    s = re.sub(r"<TITLE[^>]*>(.*?)</TITLE>", lambda m: "\n\n## " + re.sub(r"<[^>]+>", "", m.group(1)).strip() + "\n", s, flags=re.S)
    s = re.sub(r"</(TD|TH|TE|TU)>", " | ", s)
    s = re.sub(r"</TR>|<BR\s*/?>|</P>|</TABLE>", "\n", s)
    s = re.sub(r"<[^>]+>", "", s)
    s = html.unescape(s)
    s = re.sub(r"[ \t ]+", " ", s)
    s = re.sub(r"\n\s*\n+", "\n", s)
    return s


def split_sections(t: str) -> dict:
    parts = re.split(r"\n## ", t)
    out = {}
    for p in parts[1:]:
        title, _, body = p.partition("\n")
        out.setdefault(title.strip(), []).append(body)
    return {k: "\n".join(v) for k, v in out.items()}


def num(x):
    x = (x or "").replace(",", "").strip()
    return int(x) if re.fullmatch(r"-?\d+", x) else None


fin_rows, emp_rows, kw_rows, sec_index = [], [], [], []
for r in LOG:
    cc, yr, rn, co = r["corp_code"], r["연도"], r["rcept_no"], r["회사"]
    base = RAW / cc / f"{yr}_{rn}"
    # 1) 주요 계정
    fj = json.loads(Path(f"{base}_fin.json").read_text(encoding="utf-8"))
    if fj.get("status") == "000":
        for fs in ("CFS", "OFS"):
            rows = [x for x in fj["list"] if x.get("fs_div") == fs]
            if rows: break
        for x in rows:
            if x["account_nm"] in ("매출액", "영업수익", "영업이익", "영업이익(손실)", "당기순이익", "당기순이익(손실)"):
                fin_rows.append([co, yr, fs, x["account_nm"], num(x.get("thstrm_amount")), num(x.get("frmtrm_amount")), x.get("currency", "")])
    else:
        fin_rows.append([co, yr, "-", f"API 상태 {fj.get('status')} {fj.get('message','')}", None, None, ""])
    # 2) 직원 현황
    ej = json.loads(Path(f"{base}_emp.json").read_text(encoding="utf-8"))
    tot = 0
    for x in ej.get("list", []):
        emp_rows.append([co, yr, x.get("fo_bbm"), x.get("sexdstn"), x.get("rgllbr_co"), x.get("cnttk_co"), x.get("sm")])
        if (x.get("fo_bbm") or "") not in ("합계", "총계") and (x.get("sexdstn") or "") not in ("합계", "총계"):
            tot += num(x.get("sm")) or 0
    emp_rows.append([co, yr, "(행 합)", "", "", "", tot])
    # 3) 본문
    z = zipfile.ZipFile(f"{base}_doc.zip")
    main = f"{rn}.xml"
    t = xml_to_text(z.read(main).decode("utf-8", "replace"))
    (TXT / f"{co}_{yr}_{rn}_본문.txt").write_text(t, encoding="utf-8")
    secs = split_sections(t)
    for want in SECTIONS:
        hits = [k for k in secs if want in k.replace(" ", "").replace("　", "") or want.replace(" ", "") in k.replace(" ", "")]
        body = "\n".join(f"## {k}\n{secs[k]}" for k in hits)
        f = TXT / f"{co}_{yr}_{want.replace(' ', '')}.txt"
        f.write_text(body, encoding="utf-8")
        sec_index.append([co, yr, want, "; ".join(hits), len(body)])
    for ln in t.split("\n"):
        for k in KW:
            if k in ln:
                kw_rows.append([co, yr, k, ln.strip()[:400]])
                break

def w(name, head, rows):
    with open(OUT / name, "w", encoding="utf-8", newline="") as f:
        cw = csv.writer(f, delimiter="\t"); cw.writerow(head); cw.writerows(rows)

w("재무_주요계정.tsv", ["회사", "연도", "재무제표", "계정", "당기(원)", "전기(원)", "통화"], fin_rows)
w("직원현황.tsv", ["회사", "연도", "사업부문", "성별", "정규직", "계약직", "합계"], emp_rows)
w("검색어_줄.tsv", ["회사", "연도", "검색어(첫 일치)", "줄"], kw_rows)
w("장_찾기.tsv", ["회사", "연도", "찾는 장", "찾은 장 제목", "글자 수"], sec_index)

# 검색어 줄 수 요약
from collections import Counter
c = Counter((r[0], r[1], r[2]) for r in kw_rows)
for co in sorted({r[0] for r in kw_rows}):
    for yr in sorted({r[1] for r in kw_rows if r[0] == co}):
        print(co, yr, {k: c[(co, yr, k)] for k in KW if c[(co, yr, k)]})
cur, peak = tracemalloc.get_traced_memory()
print(f"\n→ {OUT} · 장 전문 {TXT}")
print(f"peak {peak/1e9:.3f} GB · {(time.time()-T0)/3600:.4f} Hour")
