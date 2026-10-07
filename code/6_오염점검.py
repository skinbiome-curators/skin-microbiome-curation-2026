"""시약 유래 오염(kitome) 흔적을 찾고, 후보 3속이 그 패턴을 보이는지 확인한다.

DNA 추출 키트·시약에는 미량의 세균 DNA가 들어 있다. 시료의 DNA가 적을수록 이 배경이
차지하는 비율이 커지므로, **시퀀싱 심도와 강한 음의 상관**을 보이는 속은 오염을 의심한다.
우리 코호트에는 음성대조 시료가 없으므로 이 상관을 대체 지표로 쓴다. 오염을 확정하는
검사가 아니라 의심 신호를 찾는 검사다.

**음성대조 참조 (2026-09-19 신규).** 수업 공개 자료 `15_season` 에 블랭크가 16S 4개·ITS 1개 들어
있다. 이것은 **다른 연구(PRJNA390040)·다른 키트·다른 실험실**의 블랭크이므로 우리 런의 음성대조가
아니다. 우리 시료의 오염을 확정하지 못한다. 쓰는 용도는 둘뿐이다 —
  (1) 시약·환경 유래로 알려진 속의 **경험적 참조 목록**을 우리 자료 안에서 만든다
  (2) 심도-상관이 지목한 속이 실제 블랭크에도 나오는지 봄으로써 **대체 지표가 작동하는지 점검**한다
n=4 이고 ITS 는 n=1 이라 비율은 참고값이다.

세균(16S)뿐 아니라 진균(ITS)도 훑는다. Malassezia 가 ④를 통과한 유일한 후보이므로
**그것이 줄어든 자리에 무엇이 들어왔는지**를 같이 본다. 한 속이 조성의 66%를 차지하는
자료에서는 그 속이 줄면 나머지가 자동으로 오르기 때문이다.

결과: results/오염점검.json
"""
import csv, json, os, statistics as st
from stats_util import save_json, TAB, RAW, CAND, f, spearman, bh, auc_p, load_595, load_42

RHO_CUT  = -0.30   # 심도와 이보다 강한 음의 상관이면 의심
DET_MIN  = 50.0    # 의심 목록에 올릴 최소 검출률(%)

out = {'method': '상대비율과 시퀀싱 심도의 스피어만 상관. 강한 음의 상관 = 시약 오염 의심.',
       'rho_threshold': RHO_CUT, 'cohorts': {}}

def scan(name, val, depth, n, top=8):
    rows = []
    for g, v in val.items():
        if g in ('unidentified', 'Unassigned', 'metagenome'):
            continue
        det = 100 * sum(1 for x in v if x > 0) / n
        if det < DET_MIN:
            continue
        r, p = spearman(depth, v)
        rows.append({'genus': g, 'detect_pct': round(det, 1),
                     'median_pct': round(100 * st.median(v) if max(v) <= 1.0 else st.median(v), 4),
                     'rho_depth': round(r, 3), 'p': p})
    for r, q in zip(rows, bh([x['p'] for x in rows])):
        r['q_bh'] = round(q, 4)
    suspect = sorted([r for r in rows if r['rho_depth'] <= RHO_CUT and r['q_bh'] < 0.05],
                     key=lambda r: r['rho_depth'])
    cand = {r['genus']: r for r in rows if r['genus'] in CAND}
    out['cohorts'][name] = {'n': n, 'n_genus_screened': len(rows),
                            'suspect': suspect[:top], 'candidates': cand}
    print(f'\n[{name}] 검출률 {DET_MIN:.0f}% 이상 {len(rows)}속 중 의심 {len(suspect)}속')
    if suspect:
        print(f"  {'속':18s} {'검출':>7s} {'중앙값':>9s} {'rho(심도)':>10s} {'q':>8s}")
        for r in suspect[:top]:
            print(f"  {r['genus']:18s} {r['detect_pct']:6.1f}% {r['median_pct']:9.3f} {r['rho_depth']:+10.3f} {r['q_bh']:8.4f}")
    print('  후보 3속 —')
    for g in CAND:
        r = cand.get(g)
        if r:
            flag = '의심' if r['rho_depth'] <= RHO_CUT else '해당 없음'
            print(f"    {g:18s} rho={r['rho_depth']:+.3f} (q={r['q_bh']:.3f})  {flag}")

# ---- 42명 (오른쪽 볼 swab) ----
info42, val42 = load_42()
depth42 = [f(r['read_depth_final']) for r in info42]
scan('42명', val42, depth42, len(info42))

# ---- 595명 (얼굴 전체 세척액) ----
ids, genera, val595, md, qc = load_595()
depth595 = [f(qc[i]['non_chimeric_reads']) for i in ids]
scan('595명', val595, depth595, len(ids))


# ---- 진균(ITS) ----------------------------------------------------------
# Malassezia 가 ④를 통과한 유일한 후보다. 그 반대편을 보지 않으면 판정이 반쪽이다.

def load_its(tag, table):
    si = {r['run']: r for r in csv.DictReader(
        open(os.path.join(RAW, tag + '_sample_info.tsv'), encoding='utf-8'), delimiter='	')}
    rd = csv.reader(open(os.path.join(RAW, table), encoding='utf-8'), delimiter='	')
    hdr = next(rd); rows = [r for r in rd]
    col = {g: i for i, g in enumerate(hdr)}
    return [si[r[0]] for r in rows], {g: [f(r[col[g]]) or 0.0 for r in rows] for g in hdr[1:]}

i3, v3 = load_its('03_sensitive_ITS', '03_sensitive_ITS_genus_pct.tsv')
i1, v1 = load_its('01_aging_ITS', '01_aging_ITS_genus_pct.tsv')
scan('42명 진균(ITS)', v3, [f(r['read_depth_final']) for r in i3], len(i3))
scan('58명 진균(01 노화 ITS)', v1, [f(r['read_depth_final']) for r in i1], len(i1))

# ---- Malassezia 가 줄어든 자리에 무엇이 오르는가 -------------------------
grp = [r['group'] for r in i3]
S = [k for k in range(len(i3)) if grp[k] == 'Sensitive']
N = [k for k in range(len(i3)) if grp[k] == 'Non-sensitive']
rise = []
for g, v in v3.items():
    det = 100 * sum(1 for x in v if x > 0) / len(v)
    if det < DET_MIN:
        continue
    a_, p_ = auc_p([v[k] for k in S], [v[k] for k in N])
    rise.append({'genus': g, 'detect_pct': round(det, 1),
                 'median_S': round(st.median([v[k] for k in S]), 3),
                 'median_N': round(st.median([v[k] for k in N]), 3),
                 'auc': round(a_, 3), 'p': p_,
                 'rho_depth': round(spearman([f(r['read_depth_final']) for r in i3], v)[0], 3)})
for r, q in zip(rise, bh([x['p'] for x in rise])):
    r['q_bh'] = round(q, 4)
rise.sort(key=lambda r: -r['auc'])
out['fungal_displacement_42'] = {
    'note': 'Malassezia 가 ITS 조성의 66%를 차지한다. 그것이 줄면 나머지가 자동으로 오른다. '
            '민감군에서 오르는 속이 피부 상재 진균인지 환경·미분류인지 본다.',
    'n_screened': len(rise), 'top_rising': rise[:8],
    'malassezia': next((r for r in rise if 'Malassezia' in r['genus']), None)}

print()
print('[42명 진균] Malassezia 가 줄어든 자리에 무엇이 오르는가')
print('  %-26s %8s %10s %10s %8s %9s' % ('속', '검출', '민감중앙', '비민감중앙', 'AUC', 'rho(심도)'))
for r in rise[:8]:
    print('  %-26s %7.1f%% %10.2f %10.2f %8.3f %+9.3f'
          % (r['genus'][:26], r['detect_pct'], r['median_S'], r['median_N'], r['auc'], r['rho_depth']))
m = out['fungal_displacement_42']['malassezia']
if m:
    print('  %-26s %7.1f%% %10.2f %10.2f %8.3f %+9.3f   ← 유일한 ④ 통과 후보'
          % (m['genus'][:26], m['detect_pct'], m['median_S'], m['median_N'], m['auc'], m['rho_depth']))
print('  → 오르는 쪽이 피부 상재가 아니면, "Malassezia 감소" 와 "시료 품질 저하" 가 구별되지 않는다.')

# ---- 음성대조 참조 (15_season 블랭크) — 2026-09-19 신규 ---------------------
# ⚠️ 다른 연구(PRJNA390040)·다른 키트·다른 실험실의 블랭크다. 우리 런의 음성대조가 아니므로
#    우리 시료의 오염을 확정하지 못한다. 경험적 참조 목록이자 대체 지표의 작동 점검이다.

BLANK = '음성대조'

def load_cohort(tag):
    """수업 공개 코호트의 시료정보와 속 상대비율(%)을 같이 읽는다."""
    si = {r['run']: r for r in csv.DictReader(
        open(os.path.join(RAW, tag + '_sample_info.tsv'), encoding='utf-8'), delimiter='	')}
    rd = csv.reader(open(os.path.join(RAW, tag + '_genus_pct.tsv'), encoding='utf-8'), delimiter='	')
    hdr = next(rd); rows = [r for r in rd]
    col = {g: i for i, g in enumerate(hdr)}
    return [si[r[0]] for r in rows], {g: [f(r[col[g]]) or 0.0 for r in rows] for g in hdr[1:]}

def blank_profile(tag, top=12):
    info, val = load_cohort(tag)
    B = [k for k in range(len(info)) if info[k]['group'] == BLANK]
    R = [k for k in range(len(info)) if info[k]['group'] != BLANK]
    prof = {}
    for g, v in val.items():
        mb = st.median([v[k] for k in B])
        mr = st.median([v[k] for k in R])
        prof[g] = {'genus': g,
                   'median_blank': round(mb, 4), 'median_sample': round(mr, 4),
                   'detect_blank_pct': round(100 * sum(1 for k in B if v[k] > 0) / len(B), 1),
                   'blank_over_sample': round(mb / mr, 2) if mr > 0 else None}
    ranked = sorted(prof.values(), key=lambda r: -r['median_blank'])
    return {'n_blank': len(B), 'n_sample': len(R), 'n_genus': len(prof),
            'top_in_blank': ranked[:top]}, prof

neg = {'source': '15_season (PRJNA390040) 의 음성대조 시료',
       'caveat': '다른 연구·다른 키트·다른 실험실·소표본(16S n=4, ITS n=1). '
                 '우리 코호트(595·42·58)의 음성대조가 아니며 우리 시료의 오염을 확정하지 못한다. '
                 '경험적 참조 목록으로만 쓴다.',
       'cohorts': {}}

prof16 = profITS = None
try:
    neg['cohorts']['15_season_16S'], prof16 = blank_profile('15_season_16S')
except Exception as e:
    neg['cohorts']['15_season_16S'] = {'error': str(e)}
try:
    neg['cohorts']['15_season_ITS'], profITS = blank_profile('15_season_ITS')
except Exception as e:
    neg['cohorts']['15_season_ITS'] = {'error': str(e)}

print()
print('[음성대조 참조] 15_season 블랭크 — 다른 연구·다른 키트. 우리 런의 대조가 아니다.')
for tag, prof in (('15_season_16S', prof16), ('15_season_ITS', profITS)):
    c = neg['cohorts'].get(tag, {})
    if prof is None or 'error' in c:
        print(f"  {tag}: 읽지 못함 ({c.get('error')})")
        continue
    print()
    print(f"  {tag} — 블랭크 {c['n_blank']}개 대 시료 {c['n_sample']}개, {c['n_genus']}속")
    print(f"    {'속':26s} {'블랭크중앙':>10s} {'시료중앙':>10s} {'배수':>8s}")
    for r in c['top_in_blank'][:8]:
        fold = f"{r['blank_over_sample']:8.1f}" if r['blank_over_sample'] is not None else '     inf'
        print(f"    {r['genus'][:26]:26s} {r['median_blank']:9.2f}% {r['median_sample']:9.2f}% {fold}")

# (2) 심도-상관이 지목한 속이 실제 블랭크에도 나오는가 — 대체 지표의 작동 점검
check = []
if prof16 or profITS:
    # 세균(16S)과 진균(ITS)을 각각 해당 블랭크와 맞춘다. 서로 다른 라이브러리이므로 섞지 않는다.
    SHORT = {'42명': '42명', '595명': '595명',
             '42명 진균(ITS)': '42명·ITS', '58명 진균(01 노화 ITS)': '58명·ITS'}
    for assay, prof, cohorts in (('16S', prof16, ('42명', '595명')),
                                 ('ITS', profITS, ('42명 진균(ITS)', '58명 진균(01 노화 ITS)'))):
        if prof is None:
            continue
        seen = []
        for cname in cohorts:
            for r in out['cohorts'].get(cname, {}).get('suspect', []):
                if r['genus'] not in seen:
                    seen.append(r['genus'])
        for g in seen:
            pr = prof.get(g)
            check.append({'assay': assay, 'genus': g,
                          'flagged_in': [c for c in cohorts
                                         if any(x['genus'] == g for x in out['cohorts'].get(c, {}).get('suspect', []))],
                          'in_blank': pr is not None,
                          'median_blank': pr['median_blank'] if pr else None,
                          'median_sample': pr['median_sample'] if pr else None,
                          'detect_blank_pct': pr['detect_blank_pct'] if pr else None,
                          'blank_over_sample': pr['blank_over_sample'] if pr else None})
    hit = [c for c in check if c['in_blank'] and c['detect_blank_pct'] == 100.0]
    neg['proxy_check'] = {
        'note': '심도-상관으로 의심한 속이 블랭크에서도 검출되면 대체 지표가 같은 것을 잡았다는 뜻이다. '
                '블랭크에 없다고 해서 오염이 아니라는 뜻은 아니다 — 키트가 다르기 때문이다. '
                'ITS 블랭크는 1개뿐이라 검출 여부만 보고 비율은 참고값이다.',
        'n_flagged': len(check), 'n_in_all_blanks': len(hit), 'detail': check}
    print()
    print('  [대체 지표 작동 점검] 심도-상관 의심 속이 블랭크에도 나오는가')
    print(f"    {'':4s}{'속':24s} {'지목한 코호트':>22s} {'블랭크검출':>10s} {'블랭크중앙':>10s} {'블랭크/시료':>11s}")
    for c in check:
        d = f"{c['detect_blank_pct']:9.0f}%" if c['detect_blank_pct'] is not None else '         -'
        m = f"{c['median_blank']:9.3f}%" if c['median_blank'] is not None else '         -'
        if c['blank_over_sample'] is not None:
            fo = f"{c['blank_over_sample']:10.1f}배"          # 블랭크 중앙 / 시료 중앙
        elif c['in_blank'] and c['median_blank'] > 0:
            fo = '       inf'                                 # 시료 중앙이 0인데 블랭크에는 있다
        else:
            fo = '         -'                                 # 양쪽 다 중앙 0 또는 목록에 없음
        tags = '/'.join(SHORT.get(x, x) for x in c['flagged_in'])
        print(f"    {c['assay']:4s}{c['genus'][:24]:24s} {tags[:22]:>22s} {d} {m} {fo}")
    print(f"    → 의심 {len(check)}속 중 {len(hit)}속이 해당 블랭크 전부에서 검출된다.")

    # Cladosporium — 우리가 심도-상관으로 떨어뜨린 진균 후보가 ITS 블랭크 1위였다
    if profITS:
        cl = next((profITS[g] for g in profITS if 'Cladosporium' in g), None)
        if cl:
            neg['cladosporium_note'] = {
                'why': '§6.4 에서 Cladosporium 은 진균 1차를 통과했으나 심도 상관 rho=-0.378 로 '
                       '오염이 의심돼 후보에서 뺐다. 그 판정을 블랭크가 독립적으로 뒷받침한다.',
                **cl,
                'caveat': 'ITS 블랭크는 1개다. 방향은 일치하지만 한 시료의 값이므로 근거 강도는 약하다.'}
            print()
            print(f"    [Cladosporium] ITS 블랭크 1위 — 블랭크 {cl['median_blank']:.2f}% / 시료 "
                  f"{cl['median_sample']:.2f}% ({cl['blank_over_sample']}배)")
            print('      심도-상관(rho=-0.378)으로 이미 뺀 후보다. 블랭크가 같은 방향을 가리킨다. '
                  '단 블랭크 1개이므로 보강 근거일 뿐이다.')

# (3) 후보 3속 + Malassezia 의 블랭크 수준
cand_blank = {}
if prof16:
    for g in CAND:
        if g in prof16:
            cand_blank[g] = prof16[g]
if profITS:
    for g in profITS:
        if 'Malassezia' in g:
            cand_blank[g] = profITS[g]
neg['candidates_in_blank'] = cand_blank
if cand_blank:
    print()
    print('  [후보 속의 블랭크 수준] — 낮다고 오염이 아니고, 있다고 오염도 아니다')
    for g, r in cand_blank.items():
        fold = f"{r['blank_over_sample']:.3f}배" if r['blank_over_sample'] is not None else 'inf'
        print(f"    {g[:24]:24s} 블랭크 {r['median_blank']:7.3f}% / 시료 {r['median_sample']:7.3f}%  ({fold})")
    print('    → Corynebacterium 은 Salter 2014 의 키트 오염균 목록에도 있다. 판정 「주의」는 유지된다.')

out['negative_control_reference'] = neg

save_json(out, os.path.join(TAB, '오염점검.json'))
print(f"\n저장: {os.path.join(TAB, '오염점검.json')}")
