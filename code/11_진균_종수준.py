"""⑤ 방향성을 진균에서 **종 수준으로 내려** 다시 판정한다 — 2026-09-19 신규.

15차 §10-2 는 *Malassezia* 의 ⑤를 이렇게 적었다.

    "restricta·globosa·sympodialis 가 한 칸에 합쳐져 방향을 말할 수 없음"

**그 전제가 틀렸다.** 수업 공개 자료의 `asv/03_sensitive_ITS_asv_taxonomy.tsv` 에는
*Malassezia* ASV 116개 중 **74개에 종 라벨이 붙어 있다**(restricta 46 · globosa 15 · sympodialis 2 …).
속 수준 요약표(`*_genus_pct.tsv`)만 쓰다 보니 못 본 것이지 자료가 없던 것이 아니다.

여기서 하는 것
  1. ASV 계수 → 사람별 **종 수준 상대비율** (시료 합 100%)
  2. 종마다 민감/비민감 AUC + 사람 단위 부트스트랩 CI + BH 보정
  3. 원저(Keum 2020)가 보고한 ***M. restricta* 감소**가 우리 자료에서 재현되는가
  4. 속 합계(66.4%)와 비교해 **종을 나누면 신호가 어디서 오는가**

⚠️ 이것이 ⑤를 통과시키는 것은 아니다. 종 라벨은 **ITS 기반 추정**이고 36%는 여전히 미분류다.
   판정은 「불가」에서 **「부분 가능 — 진균에 한정」** 으로 바뀔 뿐이다. 세균은 그대로 불가다
   (16S ASV 에서 *Cutibacterium* 종 라벨 **9%** · *Staphylococcus* **10%**).

결과: results/진균_종수준.json
"""
import csv, json, os, statistics as st
from stats_util import save_json, TAB, RAW, f, auc_p, auc_ci, bh

# 2026-09-19 — ASV 표 4개를 원자료/ 로 옮겼다. 이제 패키지 밖을 보지 않는다.
# (옛 경로: <패키지의 부모>/data/microbiome-2026-main/asv/ — 폴더를 떼어 옮기면 안 돌았다)
# 없으면 옛 자리를 찾아본다 — 이전 판올림과 호환을 위해서다.
ASV = RAW
if not os.path.exists(os.path.join(ASV, '03_sensitive_ITS_asv_taxonomy.tsv')):
    ASV = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                       'data', 'microbiome-2026-main', 'asv')
GENUS = 'g__Malassezia'
MIN_DETECT = 50.0     # 종 수준 판정에 올릴 최소 검출률(%)


def read_tsv(p):
    with open(p, encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


# ── 1. ASV → 사람별 종 수준 상대비율 ───────────────────────────────────────
tax = {r['ASV']: r for r in read_tsv(os.path.join(ASV, '03_sensitive_ITS_asv_taxonomy.tsv'))}
rd = csv.reader(open(os.path.join(ASV, '03_sensitive_ITS_asv_count.tsv'), encoding='utf-8'), delimiter='\t')
hdr = next(rd)
runs = hdr[1:]
counts = {r[0]: [f(x) or 0.0 for x in r[1:]] for r in rd}

info = {r['run']: r for r in read_tsv(os.path.join(RAW, '03_sensitive_ITS_sample_info.tsv'))}
keep = [i for i, run in enumerate(runs) if run in info]
runs = [runs[i] for i in keep]
grp = [info[r]['group'] for r in runs]
S = [k for k in range(len(runs)) if grp[k] == 'Sensitive']
N = [k for k in range(len(runs)) if grp[k] == 'Non-sensitive']

total = [0.0] * len(runs)
for asv, v in counts.items():
    for k, i in enumerate(keep):
        total[k] += v[i]

def label(asv):
    t = tax.get(asv, {})
    if GENUS not in (t.get('Genus') or ''):
        return None
    sp = (t.get('Species') or '').strip()
    return sp.replace('s__', '') if sp else '(종 미분류)'

species = {}
for asv, v in counts.items():
    lab = label(asv)
    if lab is None:
        continue
    row = species.setdefault(lab, [0.0] * len(runs))
    for k, i in enumerate(keep):
        row[k] += v[i]

pct = {sp: [100 * row[k] / total[k] if total[k] > 0 else 0.0 for k in range(len(runs))]
       for sp, row in species.items()}
genus_total = [sum(pct[sp][k] for sp in pct) for k in range(len(runs))]

out = {'source': 'data/microbiome-2026-main/asv/03_sensitive_ITS_asv_*',
       'n': len(runs), 'n_sensitive': len(S), 'n_non_sensitive': len(N),
       'n_asv_total': len(counts), 'n_asv_malassezia': sum(1 for a in counts if label(a)),
       'why': '15차 §10-2 는 restricta·globosa·sympodialis 가 한 칸이라 ⑤ 불가라고 적었다. '
              'ASV 표에는 종이 나뉘어 있다. 속 요약표만 쓴 탓이지 자료가 없던 것이 아니다.',
       'caveat': '종 라벨은 ITS 기반 추정이고 일부는 미분류다. ⑤를 통과시키는 것이 아니라 '
                 '「불가」를 「부분 가능 — 진균에 한정」으로 바꾼다. 세균은 16S ASV 에서도 '
                 'Cutibacterium 종 라벨 9% · Staphylococcus 10% 라 그대로 불가다.',
       'genus_total_median': {'sensitive': round(st.median([genus_total[k] for k in S]), 2),
                              'non_sensitive': round(st.median([genus_total[k] for k in N]), 2)},
       'species': {}}

rows, pv = [], []
for sp, v in sorted(pct.items(), key=lambda kv: -st.median(kv[1])):
    det = 100 * sum(1 for x in v if x > 0) / len(v)
    a, p = auc_p([v[k] for k in S], [v[k] for k in N])
    lo, hi = auc_ci([v[k] for k in S], [v[k] for k in N])
    rows.append({'species': sp, 'detect_pct': round(det, 1),
                 'median_sensitive': round(st.median([v[k] for k in S]), 3),
                 'median_non_sensitive': round(st.median([v[k] for k in N]), 3),
                 'auc': round(a, 3), 'ci': [round(lo, 3), round(hi, 3)], 'p': round(p, 4),
                 'crosses_half': lo <= 0.5 <= hi, 'screened': det >= MIN_DETECT})
    pv.append(p)
for r, q in zip(rows, bh(pv)):
    r['q_bh'] = round(q, 4)
out['species'] = rows
out['multiplicity'] = f'검출률 무관 전체 {len(rows)}개 종/그룹에 BH 보정'

print('=' * 96)
print('⑤ 방향성 — 진균을 종 수준으로 내려 다시 본다 (03 민감피부 42명 · ITS ASV)')
print('=' * 96)
print(f"  ASV {out['n_asv_total']}개 중 Malassezia {out['n_asv_malassezia']}개 · "
      f"속 합계 중앙 민감 {out['genus_total_median']['sensitive']}% 대 "
      f"비민감 {out['genus_total_median']['non_sensitive']}%")
print()
print(f"  {'종':22s} {'검출':>7s} {'민감중앙':>9s} {'비민감중앙':>10s} {'AUC':>7s} {'95% CI':>16s} {'q':>8s}")
for r in rows:
    if not r['screened']:
        continue
    print(f"  {r['species'][:22]:22s} {r['detect_pct']:6.1f}% {r['median_sensitive']:9.2f} "
          f"{r['median_non_sensitive']:10.2f} {r['auc']:7.3f} "
          f"[{r['ci'][0]:.3f}-{r['ci'][1]:.3f}]{'':>3s} {r['q_bh']:8.4f}")
low = [r for r in rows if r['screened'] and not r['crosses_half']]
print(f"\n  검출률 {MIN_DETECT:.0f}% 이상 {sum(1 for r in rows if r['screened'])}개 중 "
      f"CI 가 0.5 를 넘지 않는 것 {len(low)}개: " + (', '.join(r['species'] for r in low) or '없음'))

res = next((r for r in rows if r['species'] == 'restricta'), None)
if res:
    out['restricta_vs_original_paper'] = {
        'original': 'Keum 2020 이 민감군에서 M. restricta 감소를 보고했다',
        'ours': {k: res[k] for k in ('median_sensitive', 'median_non_sensitive', 'auc', 'ci', 'q_bh')},
        'reproduced': res['auc'] < 0.5 and not res['crosses_half']}
    print(f"\n  [원저 대조] Keum 2020 은 민감군에서 M. restricta 감소를 보고했다.")
    print(f"    우리 자료: 민감 {res['median_sensitive']}% 대 비민감 {res['median_non_sensitive']}% · "
          f"AUC {res['auc']} [{res['ci'][0]}-{res['ci'][1]}] · q={res['q_bh']}")
    print(f"    재현 여부: {'재현됨 (감소 방향 · CI 가 0.5 를 넘지 않음)' if out['restricta_vs_original_paper']['reproduced'] else '방향은 같으나 CI 가 0.5 를 걸침'}")

print('\n  ⚠️ 이것으로 ⑤가 통과가 되는 것이 아니다. 「불가」 → 「부분 가능 — 진균에 한정」 이다.')
print('     세균은 16S ASV 에서도 종 라벨이 Cutibacterium 9% · Staphylococcus 10% 라 그대로 불가다.')

save_json(out, os.path.join(TAB, '진균_종수준.json'))
print(f"\n저장: {os.path.join(TAB, '진균_종수준.json')}")
