"""후보 속 안에서 '좋은 균 / 나쁜 균' 방향을 말할 수 있는지 확인한다.

중간발표 피드백: "병원균이 있고 상쇄균이 있는데 그렇게라도 나눴으면 좋겠다"

같은 속 안에 방향이 반대인 종이 섞여 있으면 속 수준 관찰로는 방향을 말할 수 없다.
대표적으로 Staphylococcus 안에는 S. aureus(민감·아토피 악화)와 S. epidermidis(장벽 보조)가
함께 있다. 원자료의 종 라벨을 실제로 세어 이 한계를 수치로 보인다.

주의: 여기 나오는 종 라벨은 16S 기반 추정이며, 이 자료의 종 라벨에는 이미 오류가 확인됐다
(분석표/595_taxonomy_mismatch.tsv). 종 수준 결론의 근거로 쓰지 않고, '종을 가릴 수 없다'는
판정의 근거로만 쓴다.

결과: 분석표/종수준_한계.json
"""
import json, os, sys, statistics as st
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_util import save_json, BASE, TAB, RAW, CAND, f, load_595
from xlsx import rows as xrows

# 방향성이 알려진 대표 종. 문헌 근거의 방향이며 이 자료에서 검증한 것이 아니다.
KNOWN = {
    'Staphylococcus_aureus':      ('악화 방향', '아토피·민감피부 악화와 연관 보고'),
    'Staphylococcus_epidermidis': ('보호 방향', '장벽 지질·세라마이드 관련 기전 보고'),
    'Cutibacterium_acnes':        ('양방향',    '균주(phylotype)에 따라 역할이 갈림'),
}

ids, genera, val595, md, qc = load_595()
sel = set(r for r in (md[i]['mb_selected_sample_id'] for i in ids))

it = xrows(os.path.join(RAW, 'AAML_마이크로바이옴 분석 데이터_final.xlsx'))
hdr = next(it)
cols = {hdr[i]: i for i in range(8, len(hdr)) if hdr[i] in sel}

acc = {g: [] for g in CAND}
for r in it:
    gen, spe = r[6], r[7]
    if gen not in acc:
        continue
    tot = sum(f(r[i]) or 0.0 for i in cols.values())
    det = sum(1 for i in cols.values() if (f(r[i]) or 0.0) > 0)
    acc[gen].append({'species': spe, 'sum_relabund': tot, 'n_detected': det})

out = {'n_people': len(cols), 'note': '종 라벨은 16S 기반 추정이며 이 자료에서 오류가 확인됐다.',
       'genera': {}}

print(f'595명 종 라벨 구성 (사람 {len(cols)}명)\n')
for g in CAND:
    rows = sorted(acc[g], key=lambda x: -x['sum_relabund'])
    tot = sum(x['sum_relabund'] for x in rows) or 1e-12
    named = [x for x in rows if x['species'] not in ('unidentified', 'uncultured_bacterium', '')]
    top = [{'species': x['species'], 'share_of_genus': round(100 * x['sum_relabund'] / tot, 2),
            'n_detected': x['n_detected'],
            'known_direction': KNOWN.get(f"{g}_{x['species'].split('_')[-1]}", KNOWN.get(x['species'], ('미상', '')))[0]}
           for x in rows[:6]]
    unid = round(100 * sum(x['sum_relabund'] for x in rows
                           if x['species'] in ('unidentified', 'uncultured_bacterium', '')) / tot, 2)

    # 방향을 말할 수 있는가 — 문헌이 아니라 이 자료의 종 라벨 분포로 판정한다.
    known_here = {x['species']: KNOWN[x['species']] for x in rows if x['species'] in KNOWN
                  and x['n_detected'] > 0}
    if unid >= 50:
        verdict = f'불가 — 속 내 {unid}%가 종 수준 미분류'
    elif len({v[0] for v in known_here.values()}) >= 2:
        verdict = '불가 — 방향이 반대인 종이 한 속에 섞여 있음'
    elif any(v[0] == '양방향' for v in known_here.values()):
        verdict = '조건부 — 우점종이 균주에 따라 역할이 갈림'
    else:
        verdict = '조건부 — 방향이 알려진 종이 우점하지 않음'
    out_known = {k: {'direction': v[0], 'basis': v[1],
                     'n_detected': next(x['n_detected'] for x in rows if x['species'] == k)}
                 for k, v in known_here.items()}
    missing = [k for k in KNOWN if k.startswith(g + '_') and k not in known_here]

    out['genera'][g] = {'n_species_rows': len(rows), 'n_named_species': len(named),
                        'unidentified_share_pct': unid, 'top_species': top,
                        'known_species_found': out_known,
                        'known_species_absent': missing,
                        'direction_verdict': verdict}
    print(f'{g} — 종 라벨 {len(rows)}행 (이름 붙은 것 {len(named)}개), 미분류 비중 {unid}%')
    for t in top:
        print(f"   {t['species'][:38]:40s} 속 내 비중 {t['share_of_genus']:6.2f}%  검출 {t['n_detected']:3d}명  {t['known_direction']}")
    if missing:
        print(f"   방향이 알려졌으나 이 자료에 라벨이 없는 종: {', '.join(missing)}")
    print(f'   → ⑤ 방향성 판정: {verdict}\n')

save_json(out, os.path.join(TAB, '종수준_한계.json'))
print(f"저장: {os.path.join(TAB, '종수준_한계.json')}")
