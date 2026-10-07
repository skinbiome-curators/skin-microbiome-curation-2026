import csv, json, os, sys, statistics, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xlsx import rows as xrows


# ── OS 이식성 (2026-09-18 추가) ──────────────────────────────────────────────
# 맥에서 만든 결과를 Windows 에서 그대로 재현하기 위한 두 가지.
#   1) stdout 을 UTF-8 로 고정 — 한글 Windows 기본 cp949 는 em-dash(—) 를 못 쓴다
#   2) save_json — encoding·newline 을 못 박아 OS 와 무관하게 같은 바이트를 낸다
import sys as _sys, json as _json
try:
    _sys.stdout.reconfigure(encoding='utf-8')
    _sys.stderr.reconfigure(encoding='utf-8')
except Exception:
    pass

def save_json(obj, path):
    with open(path, 'w', encoding='utf-8', newline=chr(10)) as _fh:
        _json.dump(obj, _fh, ensure_ascii=False, indent=2)
# ────────────────────────────────────────────────────────────────────────────

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EV   = os.path.join(BASE, 'data')      # 2026-10-07 폴더 이름 변경: 원자료 → data
OUT  = os.path.join(BASE, 'results')   # 2026-10-07 폴더 이름 변경: 분석표 → results
os.makedirs(OUT, exist_ok=True)

def f(x):
    try: return float(x)
    except: return None

# ---------- 1. 595명 선정 (20-79세) ----------
meta_rows = list(csv.DictReader(open(f'{EV}/AAML_merged_data.csv', newline='', encoding='utf-8-sig')))
sel = [r for r in meta_rows if f(r['Age']) is not None and 20 <= f(r['Age']) <= 79]
assert len(sel) == 595, len(sel)
sid2meta = {r['mb.selected_sample-id']: r for r in sel}
assert len(sid2meta) == 595

# ---------- 2. 속 수준 상대풍부도 ----------
it = xrows(f'{EV}/AAML_마이크로바이옴 분석 데이터_final.xlsx')
hdr = next(it)
samp_cols = {hdr[i]: i for i in range(8, len(hdr)) if hdr[i]}
use = {s: i for s, i in samp_cols.items() if s in sid2meta}
missing = set(sid2meta) - set(use)

genus_sum = {}          # genus -> {sample: value}
mismatch = []
domain_cnt = {}
n_taxa = 0
SKIP = {'uncultured', 'unidentified', 'metagenome', 'bacterium', 'sp.', 'unculturedbacterium'}
for r in it:
    n_taxa += 1
    dom, gen, spe = r[1], r[6], r[7]
    domain_cnt[dom] = domain_cnt.get(dom, 0) + 1
    g = genus_sum.setdefault(gen, {})
    for s, i in use.items():
        v = f(r[i]) if i < len(r) else None
        if v: g[s] = g.get(s, 0.0) + v
    # 속-종 라벨 불일치 점검
    tok = (spe or '').split('_')[0]
    if gen not in ('unidentified', '') and tok and tok.lower() not in SKIP and tok[:1].isupper() \
       and tok.lower() != gen.split('_')[0].lower():
        tot = sum(f(r[i]) or 0 for i in use.values())
        det = sum(1 for i in use.values() if (f(r[i]) or 0) > 0)
        mismatch.append({'row_index': r[0], 'Domain': dom, 'Phylum': r[2], 'Genus': gen,
                         'Species': spe, 'n_detected_595': det, 'sum_relabund_595': round(tot, 6)})

genera = sorted(genus_sum)
order = [r['mb.selected_sample-id'] for r in sel]

with open(f'{OUT}/595_genus_relabund.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['metadata_id', 'mb_sample_id'] + genera)
    for s in order:
        m = sid2meta[s]
        w.writerow([m['metadata-id'], s] + [f'{genus_sum[g].get(s, 0.0):.8g}' for g in genera])

col_sum = {s: sum(genus_sum[g].get(s, 0.0) for g in genera) for s in order}
n_det = {s: sum(1 for g in genera if genus_sum[g].get(s, 0.0) > 0) for s in order}

# ---------- 3. 메타데이터 ----------
MCOLS = ['metadata-id', 'mb.selected_sample-id', 'Age', 'Age_segment', 'Sex',
         'Hydration.Cheek.avg', 'Hydration.Forehead.avg', 'TEWL.Cheek.avg', 'TEWL.Forehead.avg',
         'Oil.Cheek', 'Oil.Forehead', 'Oil.Nose', 'Oil.avg',
         'Elasticity.Cheek.R2', 'Elasticity.Cheek.R7', 'Pore.Cheek.avg',
         'Color.Cheek.L-value', 'Color.Cheek.a-value', 'Color.Cheek.b-value', 'Color.Cheek.ITA',
         'baumann.sr', 'baumann.sr_score', 'baumann.od', 'baumann.od_score',
         'baumann.pn', 'baumann.wt', 'baumann',
         'Moist.ss', 'TE.ss', 'Oil.ss', 'ITA.ss', 'R7.ss', 'OM_group', 'TE_group', 'Aging_group',
         'mb.cluster', 'mb.subcluster', 'mb_useData']
with open(f'{OUT}/595_metadata.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow([c.replace('.', '_').replace('-', '_') for c in MCOLS] + ['hydration_face_avg', 'tewl_face_avg'])
    for s in order:
        m = sid2meta[s]
        h = [f(m['Hydration.Cheek.avg']), f(m['Hydration.Forehead.avg'])]
        t = [f(m['TEWL.Cheek.avg']), f(m['TEWL.Forehead.avg'])]
        w.writerow([m[c] for c in MCOLS]
                   + [round(sum(h) / 2, 4) if all(x is not None for x in h) else '',
                      round(sum(t) / 2, 4) if all(x is not None for x in t) else ''])

# ---------- 4. 품질·제외 기록 ----------
q = {}
qit = xrows(f'{EV}/AAML_피부 임상 측정.xlsx', sheet='xl/worksheets/sheet5.xml')
qh = next(qit)
for r in qit:
    if r and r[0]:
        q[r[0]] = dict(zip(qh, r))
qmiss = [s for s in order if s not in q]
nc = sorted(f(q[s]['non-chimeric']) for s in order if s in q and f(q[s]['non-chimeric']) is not None)
p10 = nc[int(len(nc) * 0.10)]

with open(f'{OUT}/595_qc_flags.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['metadata_id', 'mb_sample_id', 'input_reads', 'filtered_reads', 'merged_reads',
                'non_chimeric_reads', 'pct_input_non_chimeric', 'analysis_status',
                'genus_relabund_sum', 'n_genus_detected', 'flag_low_depth_p10', 'flag_qc_row_missing'])
    for s in order:
        m = sid2meta[s]; qq = q.get(s, {})
        v = f(qq.get('non-chimeric', ''))
        w.writerow([m['metadata-id'], s, qq.get('input', ''), qq.get('filtered', ''), qq.get('merged', ''),
                    qq.get('non-chimeric', ''), qq.get('percentage of input non-chimeric', ''),
                    qq.get('analysis', ''), f'{col_sum[s]:.6f}', n_det[s],
                    'Y' if (v is not None and v < p10) else 'N', 'Y' if s not in q else 'N'])

with open(f'{OUT}/595_taxonomy_mismatch.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, delimiter='\t',
                       fieldnames=['row_index', 'Domain', 'Phylum', 'Genus', 'Species',
                                   'n_detected_595', 'sum_relabund_595'])
    w.writeheader()
    for d in sorted(mismatch, key=lambda x: -x['sum_relabund_595']): w.writerow(d)

# ---------- 5. 요약 ----------
cs = sorted(col_sum.values())
summ = {
    'n_people': len(order), 'n_taxa_rows': n_taxa, 'n_genera': len(genera),
    'samples_missing_in_microbiome_table': sorted(missing),
    'genus_relabund_colsum': {'min': round(cs[0], 6), 'median': round(statistics.median(cs), 6), 'max': round(cs[-1], 6)},
    'n_genus_detected': {'min': min(n_det.values()), 'median': statistics.median(n_det.values()), 'max': max(n_det.values())},
    'non_chimeric_reads': {'n': len(nc), 'min': nc[0], 'p10_cutoff': p10, 'median': statistics.median(nc), 'max': nc[-1]},
    'qc_rows_missing': qmiss,
    'taxonomy_mismatch_rows': len(mismatch),
    'taxonomy_mismatch_total_relabund_share': round(sum(d['sum_relabund_595'] for d in mismatch) / sum(cs), 6),
    'domains': domain_cnt,
    'baumann_sr': {g: sum(1 for s in order if sid2meta[s]['baumann.sr'] == g) for g in ('S', 'R')},
}
for g in ['Cutibacterium', 'Staphylococcus', 'Corynebacterium']:
    vals = [genus_sum.get(g, {}).get(s, 0.0) for s in order]
    summ[g] = {'detected': sum(1 for v in vals if v > 0),
               'median_relabund': round(statistics.median(vals), 6)}
save_json(summ, f'{OUT}/build_summary.json')
print(json.dumps(summ, ensure_ascii=False, indent=2))
