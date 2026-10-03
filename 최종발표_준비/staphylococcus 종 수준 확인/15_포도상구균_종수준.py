"""Staphylococcus 종 수준 점검 (2026-09-30) — 595명 SILVA 종 라벨 + 42명 16S ASV 종 라벨을 센다. 결과: 같은 폴더 staph_out.json"""
import os, sys, csv, json
CODE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), '코드')
sys.path.insert(0, CODE)
from stats_util import RAW, TAB, f, load_595
from xlsx import rows as xrows

out = {}
# ---------- 595 AAML (SILVA138v_Species) ----------
ids, genera, val595, md, qc = load_595()
sel = set(md[i]['mb_selected_sample_id'] for i in ids)
it = xrows(os.path.join(RAW, 'AAML_마이크로바이옴 분석 데이터_final.xlsx'))
hdr = next(it)
print('HDR', hdr[:8])
cols = [i for i in range(8, len(hdr)) if hdr[i] in sel]
staph = []; other_with_staph_species = []
for r in it:
    gen, spe = r[6], r[7]
    vals = [f(r[i]) or 0.0 for i in cols]
    if gen == 'Staphylococcus':
        staph.append((spe, vals, r[:8]))
    elif spe and 'Staphylococcus' in str(spe):
        other_with_staph_species.append((gen, spe, sum(vals)))
n = len(cols)
gtot = [sum(v[1][j] for v in staph) for j in range(n)]
T = sum(gtot)
print('595 people', n, 'staph rows', len(staph), 'detected people', sum(1 for x in gtot if x > 0))
res = []
for spe, vals, meta in sorted(staph, key=lambda x: -sum(x[1])):
    s = sum(vals); det = sum(1 for v in vals if v > 0)
    # per-person share of genus
    shares = [vals[j]/gtot[j] for j in range(n) if gtot[j] > 0]
    res.append({'species': spe, 'share_pct': round(100*s/T, 3), 'n_detected': det,
                'max_share_in_person_pct': round(100*max(shares), 1) if shares else 0,
                'max_relabund_pct': round(100*max(vals), 3)})
    print(f"  {str(spe):40s} {100*s/T:7.3f}%  det {det:3d}  maxshare {100*max(shares):5.1f}%  maxrel {100*max(vals):.3f}%  meta={meta[:6]}")
# per person: species-labelled share of Staph
named = [v for v in staph if v[0] not in ('unidentified', 'uncultured_bacterium', '', None)]
pp = []
for j in range(n):
    if gtot[j] > 0:
        pp.append(sum(v[1][j] for v in named)/gtot[j])
pp.sort()
import statistics as st
print('per-person named share: median', st.median(pp), 'mean', st.mean(pp), 'n with any named', sum(1 for x in pp if x > 0), 'n >50%', sum(1 for x in pp if x > .5))
print('other genus rows w/ Staph species', other_with_staph_species)
out['595'] = {'n_people': n, 'n_detected': sum(1 for x in gtot if x > 0), 'rows': res,
              'per_person_named_median': st.median(pp), 'per_person_named_mean': st.mean(pp),
              'n_any_named': sum(1 for x in pp if x > 0), 'n_named_gt50': sum(1 for x in pp if x > .5),
              'n_with_staph': len(pp)}

# taxonomy mismatch rows mentioning Staph
with open(os.path.join(TAB, '595_taxonomy_mismatch.tsv'), encoding='utf-8') as fh:
    lines = [l.rstrip('\n') for l in fh]
print('MISMATCH HDR', lines[0])
for l in lines[1:]:
    if 'Staphylo' in l: print('  MM', l[:200])

# ---------- 03 sensitive 16S ASV ----------
def tsv(p):
    with open(p, encoding='utf-8') as fh:
        return list(csv.reader(fh, delimiter='\t'))
tax = tsv(os.path.join(RAW, '03_sensitive_16S_asv_taxonomy.tsv'))
cnt = tsv(os.path.join(RAW, '03_sensitive_16S_asv_count.tsv'))
info = tsv(os.path.join(RAW, '03_sensitive_16S_sample_info.tsv'))
ih = info[0]; grp = {r[0]: r[ih.index('group')] for r in info[1:]}
th = tax[0]
taxd = {r[0]: r for r in tax[1:]}
ch = cnt[0]; samples = ch[1:]
print('03 samples in count', len(samples), 'in info', len(grp), 'overlap', len(set(samples) & set(grp)))
cmap = {r[0]: [int(float(x)) for x in r[1:]] for r in cnt[1:]}
st_asv = [a for a in taxd if taxd[a][th.index('Genus')] == 'Staphylococcus']
tot_reads = [sum(cmap[a][k] for a in cmap) for k in range(len(samples))]
asv_rows = []
G = [sum(cmap[a][k] for a in st_asv if a in cmap) for k in range(len(samples))]
GT = sum(G)
for a in st_asv:
    c = cmap.get(a, [0]*len(samples))
    sp = taxd[a][th.index('Species')]
    det = sum(1 for x in c if x > 0)
    detS = sum(1 for k, x in enumerate(c) if x > 0 and grp.get(samples[k]) == 'Sensitive')
    asv_rows.append({'asv': a, 'species': sp, 'reads': sum(c), 'share_pct': round(100*sum(c)/GT, 3),
                     'n_detected': det, 'n_det_sensitive': detS, 'seq': taxd[a][th.index('sequence')]})
asv_rows.sort(key=lambda x: -x['reads'])
for x in asv_rows:
    print(f"  {x['asv']} {x['species'] or '-':15s} reads {x['reads']:7d} {x['share_pct']:7.3f}% det {x['n_detected']:2d} (S {x['n_det_sensitive']})")
named03 = [x for x in asv_rows if x['species']]
pp3 = [sum(cmap[x['asv']][k] for x in named03)/G[k] for k in range(len(samples)) if G[k] > 0]
print('03 people with staph', sum(1 for g in G if g > 0), 'of', len(samples))
print('03 per-person named share median', st.median(pp3), 'mean', st.mean(pp3), 'any', sum(1 for x in pp3 if x > 0))
print('groups', {g: sum(1 for s in samples if grp.get(s) == g) for g in set(grp.values())})
out['03'] = {'asvs': asv_rows, 'n_samples': len(samples), 'n_with_staph': sum(1 for g in G if g > 0),
             'pp_named_median': st.median(pp3), 'pp_named_mean': st.mean(pp3), 'pp_any': sum(1 for x in pp3 if x > 0),
             'staph_reads': GT, 'total_reads': sum(tot_reads)}
json.dump(out, open(os.path.join(os.path.dirname(__file__), 'staph_out.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
