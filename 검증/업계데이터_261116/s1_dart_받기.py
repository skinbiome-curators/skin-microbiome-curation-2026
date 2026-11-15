"""D3 — 인체적용시험기관 2곳의 사업보고서를 DART OpenAPI 로 받는다(받기만, 분석 없음).

대상: 지에프씨생명과학(01221752 · 388610) · 피엔케이피부임상연구센타(01290381 · 347740)
받는 것(사업보고서마다): 공시 목록 · 단일회사 주요 계정(매출 · 영업이익 · 순이익) · 직원 현황 · 보고서 원문(zip)
키: micro3/.secrets/dart_api_key.txt (화면 · 파일에 남기지 않음)
원자료: micro3/data/업계데이터_261116/dart/ (저장소 밖)

사용: python -I s1_dart_받기.py
"""
import sys, json, time, tracemalloc, urllib.request, urllib.parse, re
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
tracemalloc.start(); T0 = time.time()

MICRO3 = Path(__file__).resolve().parents[3]
KEY = (MICRO3 / ".secrets" / "dart_api_key.txt").read_text(encoding="utf-8-sig").strip()
RAW = MICRO3 / "data" / "업계데이터_261116" / "dart"
RAW.mkdir(parents=True, exist_ok=True)

CORPS = {"지에프씨생명과학": "01221752", "피엔케이피부임상연구센타": "01290381"}
API = "https://opendart.fss.or.kr/api/"


def mask(s: str) -> str:
    return s.replace(KEY, "***")


def call(ep: str, **p) -> bytes:
    url = API + ep + "?" + urllib.parse.urlencode({"crtfc_key": KEY, **p})
    for i in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"}), timeout=120) as r:
                b = r.read()
            time.sleep(0.4)  # 호출 간격
            return b
        except Exception as e:
            print("  재시도", i + 1, mask(str(e)))
            time.sleep(2)
    raise RuntimeError(f"{ep} 실패")


log = []
for name, cc in CORPS.items():
    d = RAW / cc; d.mkdir(exist_ok=True)
    # 1) 정기공시 목록 — 최종 보고서만(정정 반영본)
    rows, page = [], 1
    while True:
        j = json.loads(call("list.json", corp_code=cc, bgn_de="20150101", end_de="20261116",
                            pblntf_ty="A", last_reprt_at="Y", page_no=page, page_count=100))
        if j.get("status") != "000":
            print(name, "목록 상태", j.get("status"), j.get("message")); break
        rows += j["list"]
        if page >= int(j.get("total_page", 1)): break
        page += 1
    (d / "list.json").write_text(json.dumps(rows, ensure_ascii=False, indent=1), encoding="utf-8")
    annual = [r for r in rows if r["report_nm"].strip().replace(" ", "").find("사업보고서") >= 0]
    # [첨부정정] 은 감사보고서 등 첨부만 고친 것 — 본문 원문이 없다(document.xml 014).
    # 그 경우 정정 전 목록(last_reprt_at=N)에서 같은 기간의 본문 보고서 중 가장 늦은 것을 쓴다.
    if any(r["report_nm"].startswith("[첨부정정]") for r in annual):
        allrows = json.loads(call("list.json", corp_code=cc, bgn_de="20150101", end_de="20261116",
                                  pblntf_ty="A", last_reprt_at="N", page_no=1, page_count=100)).get("list", [])
        (d / "list_정정전포함.json").write_text(json.dumps(allrows, ensure_ascii=False, indent=1), encoding="utf-8")
        fixed = []
        for r in annual:
            if r["report_nm"].startswith("[첨부정정]"):
                per = re.search(r"\(\d{4}\.\d{2}\)", r["report_nm"]).group(0)
                body = [x for x in allrows if per in x["report_nm"] and "사업보고서" in x["report_nm"]
                        and not x["report_nm"].startswith("[첨부정정]")]
                body = max(body, key=lambda x: x["rcept_no"])
                print(f"  {r['report_nm']} {r['rcept_no']} → 본문 {body['report_nm']} {body['rcept_no']}")
                r = body
            fixed.append(r)
        annual = fixed
    print(f"\n{name}: 정기공시 {len(rows)}건 · 사업보고서 {len(annual)}건")
    for r in sorted(annual, key=lambda r: r["rcept_dt"]):
        m = re.search(r"\((\d{4})\.\d{2}\)", r["report_nm"])
        yr = m.group(1) if m else "?"
        rn = r["rcept_no"]
        # 2) 주요 계정 · 3) 직원 현황 · 4) 원문
        fin = call("fnlttSinglAcnt.json", corp_code=cc, bsns_year=yr, reprt_code="11011")
        emp = call("empSttus.json", corp_code=cc, bsns_year=yr, reprt_code="11011")
        doc = call("document.xml", rcept_no=rn)
        (d / f"{yr}_{rn}_fin.json").write_bytes(fin)
        (d / f"{yr}_{rn}_emp.json").write_bytes(emp)
        (d / f"{yr}_{rn}_doc.zip").write_bytes(doc)
        fs = json.loads(fin).get("status"); es = json.loads(emp).get("status")
        ok_doc = doc[:2] == b"PK"
        print(f"  {yr} · {r['report_nm'].strip()} · 접수 {r['rcept_dt']} · {rn} · 계정 {fs} · 직원 {es} · 원문 {'zip' if ok_doc else '없음'} {len(doc):,}B")
        log.append({"회사": name, "corp_code": cc, "연도": yr, "보고서": r["report_nm"].strip(),
                    "접수일": r["rcept_dt"], "rcept_no": rn, "계정상태": fs, "직원상태": es, "원문": ok_doc})

(RAW / "받기_기록.json").write_text(json.dumps(log, ensure_ascii=False, indent=1), encoding="utf-8")
cur, peak = tracemalloc.get_traced_memory()
print(f"\n받은 사업보고서 {len(log)}건 → {RAW}")
print(f"peak {peak/1e9:.3f} GB · {(time.time()-T0)/3600:.4f} Hour")
