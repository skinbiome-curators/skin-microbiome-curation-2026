# -*- coding: utf-8 -*-
"""p1_결과.json 을 표로 요약(읽기 전용)."""
import json, os, sys
sys.stdout.reconfigure(encoding='utf-8')
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 'p1_결과.json'), encoding='utf-8'))
print('mem_time', d['mem_time'], '| rule', d['rule_check'], '| pool', d['pool15'])
s = sum(r['sec'] for r in d['results']); print('cell-sec total %.1f · per rep ms %.3f' % (s, 1000 * s / sum(r['nrep'] for r in d['results'])))
R = {(r['design'], r['truth'], round(r['delta'], 3)): r for r in d['results']}
designs = list(d['designs']); truths = list(d['truth'])
for dl in (0.0, round(d['delta'], 3)):
    print('\n== delta', dl, '(T / VE / MWp / AN) · N/arm q10-50-90 · P(one cohort) · P(cap) · Vhat/true')
    for t in truths:
        print('--', t)
        for dn in designs:
            r = R.get((dn, t, dl))
            if not r: continue
            j = r['rej']
            print('  %-15s %.4f %.4f %.4f %.4f | N %5.1f %s | one %.2f cap %.2f | V %.2f' % (
                dn, j['T'], j['VE'], j['MWp'], j['AN'], r['N_per_arm_mean'], r['N_q10_50_90'], r['P_single_cohort'], r['P_cap_hit'],
                r['Vhat_over_true_mean'] if r['Vhat_over_true_mean'] is not None else float('nan')))
