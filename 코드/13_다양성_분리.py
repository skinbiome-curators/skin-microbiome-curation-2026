"""1급 한계 — «*Malassezia* 축이 진균 다양성과 분리되지 않는다» 를 검증하고 종 수준으로 확장한다.

## 왜 이 스크립트가 필요했나

15차 판정표·README·HTML·덱 6장·문헌 2편이 **ρ = −0.979** 와 **구별력 0.347 대 0.324** 를 인용하는데,
**`분석표/` 어디에도 없고 어느 코드도 계산하지 않았다.** 재현할 수 없는 숫자였다.
(팀원이 2026-09-19 에 지운 «근거 없는 ICC 0.714 / 0.176» 과 같은 성격이다.)

그래서 이 스크립트는 **검증 + 확장** 둘 다 한다.

  A. 검증 — ρ(*Malassezia* 속, Shannon) 을 처음부터 다시 계산한다. −0.979 가 맞는가
  B. 확장 — 종 수준(*restricta*)에서도 같은 벽이 있는가
  C. 구별력 — 어느 축이 더 잘 가르는가
  D. **핵심** — *restricta* 가 다양성 **위에** 정보를 더하는가 (LOOCV)
  E. 조성 보정 — 로그비로 바꿔도 같은가
  F. 음성 대조 — 신호 없는 *globosa* 에 같은 검정. 패턴이 안 나와야 정상

## 해석 규칙 — **돌리기 전에 정했다**

우리 방법론의 약점이 «사후 규칙» 이었으므로 이번엔 미리 못 박는다.

  · *restricta* 가 Shannon 위에 더하고 역은 아니다  → **분리된다.** 한계 #6 완화
  · 둘 다 서로 위에 못 더한다                      → **분리 안 된다.** 한계 유지·강화
  · Shannon 이 *restricta* 위에 더한다             → **다양성 쪽이 축.** Malassezia 근거 약화
  · 증가분 CI 가 전부 0 을 걸친다                   → **n=42 로는 판단 불가**

## 자료 출처 — 섞여 있다는 것을 밝힌다

  · Shannon · *Malassezia* 속  → `원자료/03_sensitive_ITS_genus_pct.tsv` (패키지 안, 해시 등록됨)
  · *restricta* · *globosa*    → `원자료/03_sensitive_ITS_asv_*.tsv` (2026-09-19 부터 **패키지 안**)

Shannon 정의: **속 수준 상대비율의 자연로그 기반 H = −Σ p·ln p.** 미분류 라벨도 한 칸으로 센다.
(원자료 `sample_info` 의 `Shannon` 열은 산출 방식이 달라 값이 다르다 — 패키지 README 경고 참조)

결과: 분석표/다양성_분리.json
"""
import csv, math, os, random, statistics as st
from stats_util import save_json, TAB, RAW, f, spearman, auc_p, auc_ci, bh

SEED, N_BOOT = 20260918, 2000
# 2026-09-19 — ASV 표를 원자료/ 로 옮겼다 (11번과 같다).
ASV = RAW
if not os.path.exists(os.path.join(ASV, '03_sensitive_ITS_asv_taxonomy.tsv')):
    ASV = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                       'data', 'microbiome-2026-main', 'asv')
GEN = 'g__Malassezia'
NOT_GENUS = {'unidentified', 'Unassigned', 'metagenome'}


def read_tsv(p):
    with open(p, encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def shannon(props):
    """H = -Σ p ln p. 입력은 합이 1 이 되도록 정규화된 비율."""
    return -sum(p * math.log(p) for p in props if p > 0)


# ── 속 표에서 Shannon · Malassezia ─────────────────────────────────────────
info = {r['run']: r for r in read_tsv(os.path.join(RAW, '03_sensitive_ITS_sample_info.tsv'))}
rd = csv.reader(open(os.path.join(RAW, '03_sensitive_ITS_genus_pct.tsv'), encoding='utf-8'), delimiter='\t')
hdr = next(rd)
rows = [r for r in rd if r[0] in info]
col = {g: i for i, g in enumerate(hdr)}
runs = [r[0] for r in rows]
grp = [info[r]['group'] for r in runs]
y = [1 if g == 'Sensitive' else 0 for g in grp]
S = [i for i in range(len(runs)) if y[i] == 1]
N = [i for i in range(len(runs)) if y[i] == 0]

genera = [g for g in hdr[1:]]
val = {g: [f(r[col[g]]) or 0.0 for r in rows] for g in genera}

mal = val[GEN]
sh_all, sh_excl, rich = [], [], []
for k in range(len(rows)):
    v = [val[g][k] for g in genera]
    tot = sum(v)
    sh_all.append(shannon([x / tot for x in v]) if tot > 0 else 0.0)
    ve = [val[g][k] for g in genera if g != GEN]
    te = sum(ve)
    sh_excl.append(shannon([x / te for x in ve]) if te > 0 else 0.0)
    rich.append(sum(1 for g in genera if g not in NOT_GENUS and val[g][k] > 0))

# ── ASV 표에서 종 수준 ─────────────────────────────────────────────────────
species, asv_ok = {}, True
try:
    tax = {r['ASV']: r for r in read_tsv(os.path.join(ASV, '03_sensitive_ITS_asv_taxonomy.tsv'))}
    ard = csv.reader(open(os.path.join(ASV, '03_sensitive_ITS_asv_count.tsv'), encoding='utf-8'), delimiter='\t')
    ah = next(ard)
    cnt = {r[0]: [f(x) or 0.0 for x in r[1:]] for r in ard}
    idx = {run: i for i, run in enumerate(ah[1:])}
    keep = [idx[r] for r in runs if r in idx]
    tot = [sum(v[i] for v in cnt.values()) for i in keep]
    for asv, v in cnt.items():
        t = tax.get(asv, {})
        if GEN not in (t.get('Genus') or ''):
            continue
        sp = (t.get('Species') or '').strip().replace('s__', '') or '(미분류)'
        row = species.setdefault(sp, [0.0] * len(keep))
        for k, i in enumerate(keep):
            row[k] += v[i]
    species = {sp: [100 * r[k] / tot[k] if tot[k] > 0 else 0.0 for k in range(len(keep))]
               for sp, r in species.items()}
except Exception as e:
    asv_ok = False
    print(f'⚠️ ASV 표를 못 읽었다 ({e}) — A·C 의 속 수준 부분만 계산한다')

res = species.get('restricta')
glo = species.get('globosa')


# ── 로지스틱 · LOOCV ───────────────────────────────────────────────────────
def logistic(X, yy, iters=200, lr=0.1):
    n, k = len(X), len(X[0])
    w, b = [0.0] * k, 0.0
    for _ in range(iters):
        gw, gb = [0.0] * k, 0.0
        for i in range(n):
            z = max(-30.0, min(30.0, b + sum(w[j] * X[i][j] for j in range(k))))
            e = 1.0 / (1.0 + math.exp(-z)) - yy[i]
            gb += e
            for j in range(k):
                gw[j] += e * X[i][j]
        b -= lr * gb / n
        for j in range(k):
            w[j] -= lr * gw[j] / n
    return w, b


def loocv(cols, yy, ix=None):
    """한 사람을 빼고 맞춘 뒤 뺀 사람을 예측. 표준화도 학습 폴드 안에서만."""
    ix = list(range(len(yy))) if ix is None else ix
    n = len(ix)
    oof = [0.0] * n
    for h in range(n):
        tr = [ix[i] for i in range(n) if i != h]
        mu, sd = [], []
        for c in cols:
            v = [c[i] for i in tr]
            m = sum(v) / len(v)
            s = math.sqrt(sum((x - m) ** 2 for x in v) / max(1, len(v) - 1)) or 1.0
            mu.append(m); sd.append(s)
        X = [[(c[i] - mu[j]) / sd[j] for j, c in enumerate(cols)] for i in tr]
        w, b = logistic(X, [yy[i] for i in tr])
        xh = [(c[ix[h]] - mu[j]) / sd[j] for j, c in enumerate(cols)]
        oof[h] = b + sum(w[j] * xh[j] for j in range(len(w)))
    return auc_p([oof[i] for i in range(n) if yy[ix[i]] == 1],
                 [oof[i] for i in range(n) if yy[ix[i]] == 0])[0]


out = {'why': 'ρ=−0.979 와 구별력 0.324 가 문서에만 있고 생성기가 없었다. 검증하고 종 수준으로 확장한다.',
       'shannon_def': 'H = -Σ p·ln p · 속 수준 상대비율 · 미분류 라벨도 한 칸. '
                      'sample_info 의 Shannon 열과는 산출 방식이 달라 값이 다르다.',
       'decision_rule_prespecified':
           'restricta 가 Shannon 위에 더하고 역은 아니다 → 분리된다 / '
           '둘 다 못 더한다 → 분리 안 된다 / Shannon 이 더한다 → 다양성이 축 / '
           'CI 가 전부 0 을 걸친다 → n=42 로 판단 불가',
       'n': len(runs), 'n_sensitive': len(S), 'n_non_sensitive': len(N),
       'asv_available': asv_ok, 'seed': SEED, 'n_boot': N_BOOT}

print('=' * 92)
print('1급 한계 검증 — Malassezia 축은 진균 다양성과 분리되는가')
print('=' * 92)
print(f'  n={len(runs)} (민감 {len(S)} · 비민감 {len(N)}) · ASV 표 {"있음" if asv_ok else "없음"}')

# ── A·B. 상관 구조 ─────────────────────────────────────────────────────────
print('\n[A·B] 상관 구조 — 문서의 ρ=−0.979 가 맞는가, 종으로 내려가면 달라지는가')
pairs = [('Malassezia 속', mal, 'Shannon(전체)', sh_all),
         ('Malassezia 속', mal, 'Shannon(속 제외)', sh_excl)]
if res:
    pairs += [('restricta', res, 'Shannon(전체)', sh_all),
              ('restricta', res, 'Shannon(속 제외)', sh_excl),
              ('restricta', res, 'Malassezia 속', mal)]
if glo:
    pairs += [('globosa(음성대조)', glo, 'Shannon(전체)', sh_all)]
corr = []
print(f"  {'축 A':20s} {'축 B':20s} {'rho':>8s} {'p':>10s}")
for na, a, nb, b in pairs:
    r, p = spearman(a, b)
    corr.append({'a': na, 'b': nb, 'rho': round(r, 3), 'p': p})
    print(f'  {na:20s} {nb:20s} {r:+8.3f} {p:10.2e}')
out['correlations'] = corr
doc = next((c for c in corr if c['a'] == 'Malassezia 속' and c['b'] == 'Shannon(전체)'), None)
out['documented_rho_check'] = {
    'documented': -0.979, 'recomputed': doc['rho'] if doc else None,
    'matches': (doc is not None and abs(doc['rho'] - (-0.979)) < 0.005)}
print(f"  → 문서의 −0.979 {'재현됨' if out['documented_rho_check']['matches'] else '재현 안 됨 ←'}"
      f" (다시 계산 {doc['rho'] if doc else '?'})")

# ── C. 구별력 ──────────────────────────────────────────────────────────────
print('\n[C] 구별력 — 어느 축이 더 잘 가르는가 (AUC · |AUC−0.5| · 사람 단위 부트스트랩 CI)')
axes = [('Malassezia 속', mal), ('Shannon(전체)', sh_all),
        ('Shannon(속 제외)', sh_excl), ('관찰 속 수', [float(x) for x in rich])]
if res: axes.insert(1, ('restricta', res))
if glo: axes.append(('globosa(음성대조)', glo))
disc, pv = [], []
print(f"  {'축':20s} {'AUC':>7s} {'95% CI':>17s} {'구별력':>8s} {'p':>9s}")
for name, v in axes:
    a, p = auc_p([v[i] for i in S], [v[i] for i in N])
    lo, hi = auc_ci([v[i] for i in S], [v[i] for i in N], N=N_BOOT, seed=SEED)
    disc.append({'axis': name, 'auc': round(a, 3), 'ci': [round(lo, 3), round(hi, 3)],
                 'discrimination': round(abs(a - 0.5), 3), 'p': round(p, 4),
                 'crosses_half': lo <= 0.5 <= hi})
    pv.append(p)
    print(f'  {name:20s} {a:7.3f} [{lo:6.3f},{hi:6.3f}] {abs(a-0.5):8.3f} {p:9.4f}')
for d, q in zip(disc, bh(pv)):
    d['q_bh'] = round(q, 4)
out['discrimination'] = disc
out['multiplicity'] = f'이 표의 {len(disc)}개 축에 BH 보정. 42명에 가설을 더 얹은 것이다.'

# ── D. 핵심 — 서로 위에 정보를 더하는가 ────────────────────────────────────
if res:
    print('\n[D] 핵심 — 한쪽이 다른 쪽 «위에» 정보를 더하는가 (LOOCV · 재대입 둘 다)')

    def resub(cols):
        zc = []
        for c in cols:
            m = sum(c) / len(c)
            s = math.sqrt(sum((x - m) ** 2 for x in c) / max(1, len(c) - 1)) or 1.0
            zc.append([(x - m) / s for x in c])
        X = [[zc[j][i] for j in range(len(cols))] for i in range(len(y))]
        w, b = logistic(X, y)
        sc = [b + sum(w[j] * X[i][j] for j in range(len(w))) for i in range(len(y))]
        return auc_p([sc[i] for i in S], [sc[i] for i in N])[0]

    models = {'Shannon 만': [sh_all], 'Shannon + restricta': [sh_all, res],
              'restricta 만': [res], 'restricta + Shannon': [res, sh_all]}
    got = {k: {'resub': round(resub(v), 3), 'loocv': round(loocv(v, y), 3)} for k, v in models.items()}
    for k, v in got.items():
        print(f'  {k:22s} 재대입 {v["resub"]:.3f} · 교차검증 {v["loocv"]:.3f}')

    # 증가분 + 부트스트랩 CI (재대입 기준 — 9번과 같은 방식)
    rng = random.Random(SEED)
    boots = {'shannon_plus_restricta': [], 'restricta_plus_shannon': []}
    n = len(y)
    for _ in range(N_BOOT):
        ix = [rng.randrange(n) for _ in range(n)]
        if len(set(y[i] for i in ix)) < 2:
            continue
        def rs(cols):
            sub = [[c[i] for i in ix] for c in cols]
            zc = []
            for c in sub:
                m = sum(c) / len(c)
                s = math.sqrt(sum((x - m) ** 2 for x in c) / max(1, len(c) - 1)) or 1.0
                zc.append([(x - m) / s for x in c])
            X = [[zc[j][i] for j in range(len(sub))] for i in range(n)]
            yy = [y[i] for i in ix]
            w, b = logistic(X, yy)
            sc = [b + sum(w[j] * X[i][j] for j in range(len(w))) for i in range(n)]
            return auc_p([sc[i] for i in range(n) if yy[i] == 1],
                         [sc[i] for i in range(n) if yy[i] == 0])[0]
        boots['shannon_plus_restricta'].append(rs([sh_all, res]) - rs([sh_all]))
        boots['restricta_plus_shannon'].append(rs([res, sh_all]) - rs([res]))
    deltas = {}
    print(f"\n  {'무엇 위에 무엇을 더하나':32s} {'증가분':>8s} {'95% CI':>18s} {'교차검증 증가분':>14s}")
    for key, base, full in (('shannon_plus_restricta', 'Shannon 만', 'Shannon + restricta'),
                            ('restricta_plus_shannon', 'restricta 만', 'restricta + Shannon')):
        s_ = sorted(boots[key])
        lo, hi = s_[int(.025 * len(s_))], s_[int(.975 * len(s_))]
        d = got[full]['resub'] - got[base]['resub']
        dcv = got[full]['loocv'] - got[base]['loocv']
        deltas[key] = {'base': base, 'full': full, 'delta_resub': round(d, 3),
                       'ci': [round(lo, 3), round(hi, 3)], 'crosses_zero': lo <= 0 <= hi,
                       'delta_loocv': round(dcv, 3)}
        lab = f'{base} → +{full.split("+")[1].strip()}'
        print(f'  {lab:32s} {d:+8.3f} [{lo:+.3f},{hi:+.3f}] {dcv:+14.3f}')
    out['added_information'] = {'models': got, 'deltas': deltas}

    a, b = deltas['shannon_plus_restricta'], deltas['restricta_plus_shannon']
    if a['crosses_zero'] and b['crosses_zero']:
        verdict = ('판단 불가 — 양쪽 증가분 CI 가 모두 0 을 걸친다. n=42 로는 어느 쪽이 축인지 말할 수 없다. '
                   '한계 #6 은 그대로 유지된다.')
    elif not a['crosses_zero'] and b['crosses_zero']:
        verdict = 'restricta 가 다양성 위에 정보를 더하고 역은 아니다 → 분리된다. 한계 #6 완화.'
    elif a['crosses_zero'] and not b['crosses_zero']:
        verdict = 'Shannon 이 restricta 위에 정보를 더한다 → 다양성 쪽이 축이다. Malassezia 근거가 약해진다.'
    else:
        verdict = '양쪽 다 서로 위에 더한다 — 두 축이 다른 정보를 담고 있다. 해석에 주의.'
    out['verdict'] = verdict
    print(f'\n  판정(사전 규칙 적용): {verdict}')

# ── E. 조성 보정 ───────────────────────────────────────────────────────────
if res:
    EPS = 1e-5
    rest_lr = [math.log((res[k] + EPS) / (100 - mal[k] + EPS)) for k in range(len(runs))]
    a, p = auc_p([rest_lr[i] for i in S], [rest_lr[i] for i in N])
    lo, hi = auc_ci([rest_lr[i] for i in S], [rest_lr[i] for i in N], N=N_BOOT, seed=SEED)
    out['logratio'] = {'note': 'restricta 를 «Malassezia 속을 뺀 나머지» 대비 로그비로 바꿔 다시 잰다. '
                               '«restricta 가 낮다» 와 «나머지가 많다» 를 가르기 위한 것.',
                       'auc': round(a, 3), 'ci': [round(lo, 3), round(hi, 3)],
                       'crosses_half': lo <= 0.5 <= hi, 'eps': EPS}
    print(f'\n[E] 조성 보정 — restricta 로그비 AUC {a:.3f} [{lo:.3f},{hi:.3f}]'
          f'{"  (0.5 를 걸침)" if lo <= 0.5 <= hi else ""}')

print('\n[한계] n=42(23/19) 라 증가분 신뢰구간이 넓다. 이 계산은 «모른다의 크기» 를 재는 것이다.')
print('       restricta·globosa 는 원자료/ 의 ASV 표에서 왔다 (15.13 해소 · 해시 등록됨).')

save_json(out, os.path.join(TAB, '다양성_분리.json'))
print(f"\n저장: {os.path.join(TAB, '다양성_분리.json')}")
