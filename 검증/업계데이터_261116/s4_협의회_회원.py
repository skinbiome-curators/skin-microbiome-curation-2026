"""I-06 확인 — 한국인체적용시험기관협의회 회원기관이 홈페이지에 마이크로바이옴 시험을 내세우는가.

1) 협의회 '회원기관 명단'(www.kahsrc.or.kr/member/member.php) → 기관 이름 · 소개문 · 홈페이지
2) 기관마다 홈페이지 첫 화면 + 같은 도메인 링크(최대 MAXP 쪽)를 받아 검색어를 찾는다
   검색어: 마이크로바이옴 · microbiome · 상재균 · 16S · 미생물 군집 · 균총
3) 결과: 결과_i06/회원기관.tsv (기관 · 홈페이지 · 받은 쪽 · 검색어가 나온 쪽 · 첫 문장 · 그 주소)
원자료(HTML): micro3/data/업계데이터_261116/kahsrc/ (저장소 밖)

한계: 홈페이지에 없다고 그 시험을 안 한다는 뜻은 아니다(‘못 찾음’). 첫 화면에서 한 단계 링크까지만 본다.
사용: python -I s4_협의회_회원.py
"""
import sys, re, html, json, time, tracemalloc, hashlib, csv
import urllib.request, urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
tracemalloc.start(); T0 = time.time()
HERE = Path(__file__).resolve().parent
RAW = HERE.parents[2] / "data" / "업계데이터_261116" / "kahsrc"; RAW.mkdir(parents=True, exist_ok=True)
OUT = HERE / "결과_i06"; OUT.mkdir(exist_ok=True)
MAXP = 40
KW = re.compile(r"마이크로바이옴|microbiome|상재균|16S|미생물\s*군집|균총", re.I)
UA = {"User-Agent": "Mozilla/5.0"}


def fetch(url, timeout=20):
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
            raw = r.read(3_000_000); ct = r.headers.get("Content-Type", ""); final = r.geturl()
    except Exception as e:
        return None, None, str(e)[:80]
    m = re.search(rb'charset=["\']?([\w-]+)', raw[:4000]) or re.search(r"charset=([\w-]+)", ct)
    enc = (m.group(1).decode() if isinstance(m.group(1), bytes) else m.group(1)) if m else "utf-8"
    try:
        s = raw.decode(enc, "replace")
    except LookupError:
        s = raw.decode("utf-8", "replace")
    return s, final, None


def text_of(s):
    s = re.sub(r"(?is)<(script|style|noscript)[^>]*>.*?</\1>", " ", s)
    s = html.unescape(re.sub(r"<[^>]+>", "\n", s))
    return [l.strip() for l in re.sub(r"[ \t ]+", " ", s).split("\n") if l.strip()]


# 1) 회원 명단
src, _, err = fetch("http://www.kahsrc.or.kr/member/member.php")
if err: sys.exit("명단 받기 실패 " + err)
(RAW / "member.php.html").write_text(src, encoding="utf-8")
lines = text_of(src)
i = next(k for k, l in enumerate(lines) if "가다나 순으로" in l)
members, cur = [], None
for k in range(i + 1, len(lines)):
    l = lines[k]
    if l.startswith("대표자") and len(l) < 40:
        cur = {"기관": lines[k - 1], "대표자": l.replace("대표자", "").strip(), "소개": []}
        members.append(cur)
        if len(members) > 1:
            members[-2]["소개"] = members[-2]["소개"][:-1]  # 다음 기관 이름 줄 빼기
    elif cur is not None:
        cur["소개"].append(l)
for m in members:
    t = " ".join(m["소개"])
    u = re.findall(r"https?://[^\s]+", t)
    m["홈페이지"] = u[0].rstrip("/") if u else ""
    m["소개문_검색어"] = "; ".join(sorted(set(KW.findall(t))))
print(f"회원기관 {len(members)}곳")


# 2) 홈페이지 + 한 단계 링크
def crawl(m):
    home = m["홈페이지"]
    if not home: return m, [], 0, "홈페이지 없음"
    s, final, err = fetch(home)
    if err: return m, [], 0, "첫 화면 실패: " + err
    dom = urllib.parse.urlparse(final).netloc
    links = []
    for h in re.findall(r'href=["\']([^"\'#]+)', s):
        u = urllib.parse.urljoin(final, html.unescape(h))
        p = urllib.parse.urlparse(u)
        if p.netloc == dom and p.scheme in ("http", "https") and not re.search(r"\.(jpg|jpeg|png|gif|pdf|zip|hwp|hwpx|docx?|xlsx?|pptx?|mp4)$", p.path, re.I):
            if u not in links and u != final: links.append(u)
    pages = [(final, s)]
    for u in links[:MAXP]:
        s2, f2, e2 = fetch(u)
        if s2: pages.append((f2, s2))
    hits = []
    for u, s in pages:
        name = hashlib.md5(u.encode()).hexdigest()[:10]
        (RAW / f"{dom}_{name}.html").write_text(s, encoding="utf-8")
        for l in text_of(s):
            if KW.search(l) and len(l) < 600:
                hits.append((u, l)); break
    return m, hits, len(pages), ""


rows = []
with ThreadPoolExecutor(max_workers=8) as ex:
    futs = [ex.submit(crawl, m) for m in members]
    for n, f in enumerate(as_completed(futs), 1):
        m, hits, np_, note = f.result()
        rows.append([m["기관"], m["홈페이지"], m["소개문_검색어"], np_, len(hits),
                     hits[0][1][:300] if hits else "", hits[0][0] if hits else "", note])
        print(f"  [{n}/{len(members)}] {m['기관'][:20]:20s} 쪽 {np_:2d} · 검색어 쪽 {len(hits)} {note}", flush=True)

rows.sort(key=lambda r: r[0])
with open(OUT / "회원기관.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["기관", "홈페이지", "협의회 소개문 검색어", "받은 쪽", "검색어가 나온 쪽", "첫 문장", "그 주소", "비고"])
    w.writerows(rows)
hit = [r for r in rows if r[4] or r[2]]
print(f"\n회원 {len(rows)}곳 · 홈페이지나 소개문에 검색어가 있는 곳 {len(hit)} · 첫 화면 실패 {sum(1 for r in rows if r[7])}")
cur_, peak = tracemalloc.get_traced_memory()
print(f"peak {peak/1e9:.3f} GB · {(time.time()-T0)/3600:.4f} Hour")
