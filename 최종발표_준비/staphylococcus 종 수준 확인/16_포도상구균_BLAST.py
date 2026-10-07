"""42명 16S Staphylococcus ASV 20개를 NCBI BLAST(웹 URL API)에 보낸다 (2026-09-30).

DB 두 개:
  rRNA_typestrains/16S_ribosomal_RNA — NCBI 16S 기준균주(type strain) DB. 종 동정의 표준 참조
  core_nt           — 미분류 진균 BLAST(results/미분류진균_BLAST_*)와 같은 DB
결과: 포도상구균_BLAST_입력.fa · 포도상구균_BLAST_원본_<db>.json · 포도상구균_BLAST_요약_<db>.tsv
NCBI DB 가 갱신되면 hit 가 달라질 수 있다.

실행
  python 16_포도상구균_BLAST.py                   NCBI 에 제출(또는 RID_<db> 환경변수로 이어받기) → 원본 JSON → 요약
  python 16_포도상구균_BLAST.py --summarize-only  제출하지 않고, 저장된 원본 JSON 에서 요약 TSV 만 다시 만든다
                                                  (FASTA 도 다시 쓰지 않는다 · 원본 JSON 이 없는 DB 는 건너뛴다)

변경 기록 (2026-10-02)
  - 요약 단계를 폴더에 실제로 있는 요약 파일(포도상구균_BLAST_요약_core_nt.tsv)의 형식으로 바꿨다.
    이전 코드는 다른 이름(포도상구균_BLAST_요약.tsv)·다른 열(db, top_ident …)로 써서 그 파일을 다시 만들 수 없었다.
    새 요약 단계로 원본 JSON 에서 다시 만든 결과가 기존 파일과 바이트 단위로 같음을 확인했다.
  - --summarize-only 추가 (RID 는 하루 남짓 뒤 만료되므로, 제출 없이 요약을 재현할 수 있어야 한다)
  - species() 삭제 — 새 요약 단계는 label() 을 쓴다
  - 점검 문서 §5-3 의 RID BTY8VWS5014(ASV0006 단독) · BTYZGGFV016(Staphylococcaceae 한정, hit 500)은
    이 스크립트로 제출한 검색이 아니다(이 스크립트는 항상 ASV 20개 · hit 50 · 조건 없음으로 제출한다)
"""
import os, sys, csv, json, time, re, io, zipfile, tracemalloc, urllib.request, urllib.parse
from collections import Counter

try:   # 한글 Windows 콘솔(cp949)에서 출력이 깨지지 않게
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
T0 = time.perf_counter()
tracemalloc.start()

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'data')   # 2026-10-07 폴더 이름 변경: 원자료 → data
URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
DBS = ['core_nt', 'rRNA_typestrains/16S_ribosomal_RNA']
SUMMARIZE_ONLY = '--summarize-only' in sys.argv[1:]

tax = list(csv.reader(open(os.path.join(RAW, '03_sensitive_16S_asv_taxonomy.tsv'), encoding='utf-8'), delimiter='\t'))
h = tax[0]
asvs = [(r[0], r[h.index('Species')], r[h.index('sequence')]) for r in tax[1:] if r[h.index('Genus')] == 'Staphylococcus']
fa = ''.join(f'>{a}\n{s}\n' for a, _, s in asvs)
if not SUMMARIZE_ONLY:
    open(os.path.join(HERE, '포도상구균_BLAST_입력.fa'), 'w', encoding='utf-8', newline='\n').write(fa)

# 요약의 read 비중·검출 인원용 — 42명 ASV count (15_포도상구균_종수준.py 와 같은 계산)
cnt = list(csv.reader(open(os.path.join(RAW, '03_sensitive_16S_asv_count.tsv'), encoding='utf-8'), delimiter='\t'))
cmap = {r[0]: [int(float(x)) for x in r[1:]] for r in cnt[1:]}
SILVA = {a: sp for a, sp, _ in asvs}
STAPH_READS = sum(sum(cmap.get(a, [])) for a in SILVA)

def post(params):
    data = urllib.parse.urlencode(params).encode()
    return urllib.request.urlopen(urllib.request.Request(URL, data=data), timeout=120).read()

def run(db):
    rid = os.environ.get('RID_' + db.split('/')[-1])   # 이미 제출한 검색을 이어받을 때
    if not rid:
        txt = post({'CMD': 'Put', 'PROGRAM': 'blastn', 'MEGABLAST': 'on', 'DATABASE': db,
                    'QUERY': fa, 'HITLIST_SIZE': '50', 'TOOL': 'student_project'}).decode()
        rid = re.search(r'RID = (\S+)', txt).group(1)
        time.sleep(30)
    print(db, 'RID', rid, flush=True)
    while True:
        st = urllib.request.urlopen(f'{URL}?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}', timeout=120).read().decode()
        m = re.search(r'Status=(\w+)', st).group(1)
        print(' ', m, flush=True)
        if m == 'READY': break
        if m in ('FAILED', 'UNKNOWN'): raise SystemExit(f'{db} {m}')
        time.sleep(20)
    raw = urllib.request.urlopen(f'{URL}?CMD=Get&FORMAT_TYPE=JSON2&RID={rid}', timeout=300).read()
    z = zipfile.ZipFile(io.BytesIO(raw))
    reports = []
    for n in z.namelist():
        d = json.loads(z.read(n))
        if 'BlastOutput2' in d:
            reports += [x['report'] for x in (d['BlastOutput2'] if isinstance(d['BlastOutput2'], list) else [d['BlastOutput2']])]
    json.dump({'rid': rid, 'db': db, 'reports': reports},
              open(os.path.join(HERE, f'포도상구균_BLAST_원본_{db.split("/")[-1]}.json'), 'w', encoding='utf-8'), ensure_ascii=False)
    return reports

def label(title):
    """hit 제목 → 종 이름. 'Uncultured <소문자 일반명>'(bacterium·organism 등) = None(uncultured),
    'Uncultured <대문자 속명> …' = 그 속명부터 두 단어(예: Staphylococcus sp.), 그 외 = 제목 첫 두 단어."""
    w_ = title.split()
    if w_ and w_[0].lower() == 'uncultured':
        return None if len(w_) < 2 or w_[1][:1].islower() else ' '.join(w_[1:3])
    return ' '.join(w_[:2])


def summarize(db, reports):
    """원본 JSON 보고서 → 요약 TSV. 동률 = 질의 길이 95% 이상 정렬된 hit 중 최고 일치도와 같은 hit."""
    rows = []
    for rep in reports:
        s = rep['results']['search']
        q = s['query_title'].split()[0]; qlen = s['query_len']
        hits = [(100 * x['hsps'][0]['identity'] / x['hsps'][0]['align_len'],
                 100 * x['hsps'][0]['align_len'] / qlen, x['description'][0]['title']) for x in s.get('hits', [])]
        c = cmap.get(q, [])
        share, det = round(100 * sum(c) / STAPH_READS, 3), sum(1 for x in c if x > 0)
        ok = [x for x in hits if x[1] >= 95]
        if ok:
            best = max(x[0] for x in ok)
            labs = [label(x[2]) for x in ok if x[0] >= best - 1e-9]
            named = Counter(l for l in labs if l)
            rows.append((share, [q, qlen, share, SILVA.get(q, ''), det, f'{best:.2f}', len(labs), labs.count(None),
                                 '; '.join(f'{k}({v})' for k, v in named.most_common()), 'Y' if len(labs) == 50 else '']))
        else:
            rows.append((share, [q, qlen, share, SILVA.get(q, ''), det, '', 0, 0, '', '']))
    path = os.path.join(HERE, f'포도상구균_BLAST_요약_{db.split("/")[-1]}.tsv')
    with open(path, 'w', encoding='utf-8', newline='\n') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['asv', 'len', 'share_of_staph_reads_pct', 'silva_species', 'n_detected_42', 'best_ident_cov95',
                    'n_hits_at_best', 'n_uncultured_at_best', 'named_species_at_best(count)', 'all_50_hits_tied'])
        for _, r in sorted(rows, key=lambda x: -x[0]):
            w.writerow(r)
    print(f'요약 저장: {os.path.basename(path)} ({len(rows)} ASV)')


for db in DBS:
    raw_json = os.path.join(HERE, f'포도상구균_BLAST_원본_{db.split("/")[-1]}.json')
    if SUMMARIZE_ONLY:
        if not os.path.exists(raw_json):
            print(f'{db}: 원본 JSON 없음 → 건너뜀')
            continue
        reps = json.load(open(raw_json, encoding='utf-8'))['reports']
    else:
        reps = run(db)
    summarize(db, reps)

cur, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()
print(f'메모리(Python 할당 기준) 현재 {cur / 1024 ** 3:.4f} GB · peak {peak / 1024 ** 3:.4f} GB · '
      f'소요 시간 {(time.perf_counter() - T0) / 3600:.6f} Hour')
