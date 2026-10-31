import math
from scipy.stats import norm, t
d = math.log(2)
def solve(Sw):
    var = 1.5 * Sw**2
    n0 = 2 * (norm.ppf(0.975) + norm.ppf(0.8))**2 * var / d**2
    # t: smallest n with n >= 2(t_{.975,df}+t_{.8,df})^2 var/d^2, df=2n-2
    for n in range(2, 500):
        df = 2 * n - 2
        if n >= 2 * (t.ppf(0.975, df) + t.ppf(0.8, df))**2 * var / d**2:
            nt = n; break
    n1 = math.ceil(nt / 0.864)
    out = [round(n0, 2), nt, n1]
    for s in (0.95, 0.90):
        per = math.ceil(n1 / (s**3 * 0.9))
        out += [per, 2 * per, 2 * per * 3]
    return out
for name, Sw in [("07 2way", 1.0352 / (1.96 * math.sqrt(2))), ("07 80%up", 0.431), ("15 2way", 0.9924)]:
    print(name, round(Sw, 4), solve(Sw))
p = 0.236
print("ratio x2 -> proportion fold (other unchanged, 3 cands excluded approx):", ((2*p/(1-p))/(1+2*p/(1-p)))/p)
