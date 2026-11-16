"""D1 — 식약처 「기능성화장품 보고품목정보」 전체를 받는다(받기만).

API: https://apis.data.go.kr/1471000/FtnltCosmRptPrdlstInfoService01/getRptPrdlstInq01
     (공공데이터포털 15095680 페이지 원문 주소 · 11/16 확인. 버전 없는 옛 주소는 '등록되지 않은 서비스키'를 냄)
키: micro3/.secrets/data_go_kr_key.txt (Encoding 키 — 그대로 붙인다 · 화면 · 파일에 남기지 않음)
원자료: micro3/data/업계데이터_261116/d1/page_XXXX.json.gz (저장소 밖) — 이미 받은 쪽은 건너뛴다(다시 실행 가능)
한 쪽 500행 · 동시 6개 · 쪽마다 3번까지 다시 시도

사용: python -I s5_d1_받기.py
"""
import sys, json, gzip, math, time, tracemalloc, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
tracemalloc.start(); T0 = time.time()
HERE = Path(__file__).resolve().parent
MICRO3 = HERE.parents[2]
KEY = (MICRO3 / ".secrets" / "data_go_kr_key.txt").read_text(encoding="utf-8-sig").strip()
RAW = MICRO3 / "data" / "업계데이터_261116" / "d1"; RAW.mkdir(parents=True, exist_ok=True)
EP = "https://apis.data.go.kr/1471000/FtnltCosmRptPrdlstInfoService01/getRptPrdlstInq01"
N = 500


def get(page):
    url = f"{EP}?serviceKey={KEY}&pageNo={page}&numOfRows={N}&type=json"
    last = ""
    for i in range(3):
        try:
            b = urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"}), timeout=180).read()
            j = json.loads(b)
            if j.get("header", {}).get("resultCode") == "00":
                return j
            last = str(j.get("header"))[:120]
        except Exception as e:
            last = str(e).replace(KEY, "***")[:120]
        time.sleep(3 * (i + 1))
    raise RuntimeError(f"쪽 {page} 실패: {last}")


first = get(1)
total = int(first["body"]["totalCount"])
pages = math.ceil(total / N)
(RAW / "page_0001.json.gz").write_bytes(gzip.compress(json.dumps(first, ensure_ascii=False).encode()))
meta = {"받은_시각_시작": time.strftime("%Y-%m-%d %H:%M:%S"), "totalCount": total, "쪽": pages, "쪽당": N, "주소": EP}
print(f"전체 {total:,}행 · {pages}쪽", flush=True)

todo = [p for p in range(2, pages + 1) if not (RAW / f"page_{p:04d}.json.gz").exists()]
done, fail = pages - len(todo), []


def work(p):
    j = get(p)
    (RAW / f"page_{p:04d}.json.gz").write_bytes(gzip.compress(json.dumps(j, ensure_ascii=False).encode()))
    return p, len(j["body"].get("items") or [])


with ThreadPoolExecutor(max_workers=6) as ex:
    futs = {ex.submit(work, p): p for p in todo}
    for f in as_completed(futs):
        try:
            f.result(); done += 1
        except Exception as e:
            fail.append(futs[f]); print("  ", e, flush=True)
        if done % 20 == 0 or done == pages:
            el = time.time() - T0
            print(f"  {done}/{pages}쪽 · {el/60:.1f}분 · 남은 예상 {el/max(done,1)*(pages-done)/60:.1f}분", flush=True)

meta.update({"받은_시각_끝": time.strftime("%Y-%m-%d %H:%M:%S"), "실패한_쪽": fail})
(RAW / "받기_기록.json").write_text(json.dumps(meta, ensure_ascii=False, indent=1), encoding="utf-8")
cur, peak = tracemalloc.get_traced_memory()
print(f"\n끝 · 실패 {len(fail)}쪽 {fail[:10]} → {RAW}")
print(f"peak {peak/1e9:.3f} GB · {(time.time()-T0)/3600:.4f} Hour")
