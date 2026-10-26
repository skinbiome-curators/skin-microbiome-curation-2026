"""인디발판 3절 '사람 수 산수' 재계산 + 대안 설계 비교.
가정: 결과 변수 = 개인 변화량(사후 - 기준), 분산 s2 는 군마다 같다.
비교 하나(브랜드 군 vs 대조)의 분산을 '혼자 하는 2군 시험(군마다 n명)'과 같게 맞춘다.
"""
import math, sys
sys.stdout.reconfigure(encoding='utf-8')

print('A. 문서 표 그대로(군마다 n명 · 동등 배정)')
for k in [2, 3, 4, 6, 10, 100]:
    per = (k + 1) / k
    print(f'  k={k:>3}: 브랜드당 {per:.4f}n · 절감 {1 - per / 2:.1%}')

print('\nB. 대조군을 √k 배로 키우되, 비교 하나의 분산은 혼자 2군 시험과 같게 유지(Dunnett 형 최적 배정 · 계산)')
# 1/nt + 1/nc = 2/n, 총합 k*nt + nc 최소 -> nc = sqrt(k)*nt
for k in [2, 3, 4, 6, 10, 100]:
    nt = (1 + 1 / math.sqrt(k)) / 2      # n 단위
    nc = math.sqrt(k) * nt
    total = k * nt + nc
    per = total / k
    print(f'  k={k:>3}: 군당 {nt:.3f}n · 대조 {nc:.3f}n · 브랜드당 {per:.4f}n · 절감 {1 - per / 2:.1%}  (동등 배정 {1 - (k + 1) / k / 2:.1%})')

print('\nC. 반쪽 얼굴(split-face) · 브랜드 혼자 · 두 쪽 변화량 상관 rho (계산)')
# 짝 차이 분산 = 2 s2 (1 - rho). 같은 분산(2 s2 / n) -> m = n (1 - rho) 명
for rho in [0.0, 0.3, 0.5, 0.7]:
    m = 1 - rho
    print(f'  rho={rho}: 사람 {m:.2f}n (절감 {1 - m / 2:.0%}) · 시료(두 쪽) {2 * m:.2f}n x 시점 (혼자 2군 2n x 시점 대비 {1 - m:.0%} 절감)')

print('\nD. 시료 수 · 시퀀싱만 어림 (팀 문서 간이견적 O2 시료당 $63-$278 · 16S 만이면 더 낮음 · 원 출처 미확인)')
for n in [20, 40]:
    for k in [2, 4]:
        persons = (k + 1) * n / k
        T = 3  # 기준 2회 + 4주
        samples = persons * T
        print(f'  n={n} k={k}: 브랜드당 {persons:.0f}명 x {T}시점 = {samples:.0f}시료 -> ${samples * 63:,.0f}-${samples * 278:,.0f} (16S 낮은 끝 $42 이면 ${samples * 42:,.0f})')

print('\nE. 라운드 채우기 (확인된 숫자 x 팀 시나리오 가정)')
for firms in [23, 29]:
    for reach in [0.1, 0.2, 0.3]:
        print(f'  해마다 {firms}곳 x 도달률 {reach:.0%} = {firms * reach:.1f}곳/년  -> 분기 라운드(4회) x k=4 에 필요한 16곳의 {firms * reach / 16:.0%}')

print('\nF. 공유 대조 비교 사이 상관 (동등 배정, 변화량 분산 같음)')
for k in [2, 4]:
    # cov(T1-C, T2-C) = var(C) ; corr = (s2/n)/(2 s2/n) = 0.5
    print(f'  k={k}: 두 브랜드의 대조 대비 추정치 상관 = 0.5 (대조군 우연 편차가 모든 브랜드에 같은 방향으로 실림)')
