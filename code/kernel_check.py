"""Numerical check of the negative-jump kernel on the explicit critical sequence of Section 6.

q_1 = 0, q_k = C k^-4 Z^(1-k) / binom(2k-1, k-1) for k >= 2, with C = 1/(zeta(3)-1), Z = (zeta(3)-1)/(zeta(3)-zeta(4)).
Checks:
  1. the critical equations f(Z) = 1 - 1/Z and Z^2 f'(Z) = 1 (40 digits);
  2. nu(-k) computed two ways agree: from the partition function, nu(-k) = (2/c) d(k-1) I_{k-1} with the closed form of
     I_j, and from the kernel nu(-k) = d(k)/(2k-1) * T(k) (40 digits for k <= 50; a direct double-precision sum with
     an explicit tail for k = 200, 1000, where extrapolated summation of T loses digits);
  3. nu has total mass one, the negative part through (2/c) int_0^1 J'(r) (1-r)^(-1/2) dr;
  4. h(j) = 2 j d(j) is nu-harmonic at p = 1, ..., 6 (h = 0 on j <= 0).
"""
import mpmath as mp
import numpy as np
from scipy.special import gammaln

mp.mp.dps = 40
z2, z3, z4 = mp.zeta(2), mp.zeta(3), mp.zeta(4)
C = 1 / (z3 - 1)
Z = (z3 - 1) / (z3 - z4)
c = 4 * Z


def d(j):
    return mp.binomial(2 * j, j) / mp.mpf(4) ** j


def a(k):                       # a_k = q_k binom(2k-1,k-1) Z^(k-1)
    return C / mp.mpf(k) ** 4 if k >= 2 else mp.mpf(0)


def nu_pos(m):                  # nu(m) = q_{m+1} c^m, m >= 0
    k = m + 1
    return a(k) * mp.mpf(4) ** m / mp.binomial(2 * k - 1, k - 1)


print("1. critical equations")
f_Z = C * (z4 - 1)
fp_Z = C * (z3 - z4) / Z
print("   f(Z) - (1 - 1/Z)  = %s" % mp.nstr(f_Z - (1 - 1 / Z), 5))
print("   Z^2 f'(Z) - 1     = %s" % mp.nstr(Z ** 2 * fp_Z - 1, 5))
print("   same, by summing the series: %s, %s" % (
    mp.nstr(mp.nsum(a, [2, mp.inf]) - (1 - 1 / Z), 5),
    mp.nstr(Z * mp.nsum(lambda k: (k - 1) * a(k), [2, mp.inf]) - 1, 5)))


def I(j):                       # closed form of int_0^1 r^j J'(r) dr; I_0 = 1
    if j == 0:
        return mp.mpf(1)
    return Z * C * (-z3 / (j * (j + 1)) + z2 / mp.mpf(j) ** 2 - mp.harmonic(j) / mp.mpf(j) ** 3)


def nu_neg(k):                  # nu(-k) from the partition function
    return 2 / c * d(k - 1) * I(k - 1)


def w(m):
    return m * (2 * m + 1) * d(m) * nu_pos(m)


def T(x):
    return mp.nsum(lambda m: w(m) / (m + x), [1, mp.inf])


print("2. nu(-k): partition function against the kernel d(k) T(k)/(2k-1)")
for k in [1, 2, 3, 5, 10, 50]:
    A, B = nu_neg(k), d(k) / (2 * k - 1) * T(k)
    print("   k=%5d  partition %s  kernel %s  rel.diff %s" % (k, mp.nstr(A, 15), mp.nstr(B, 15), mp.nstr(abs(A / B - 1), 3)))
M = 10_000_000
m = np.arange(1, M + 1, dtype=np.float64)
logd = gammaln(2 * m + 1) - 2 * gammaln(m + 1) - m * np.log(4.0)
nu_m = float(C) / (m + 1) ** 4 * np.exp(m * np.log(4.0) - (gammaln(2 * m + 2) - gammaln(m + 1) - gammaln(m + 2)))
w_m = m * (2 * m + 1) * np.exp(logd) * nu_m
Aw = w_m[-1] * M ** 2                                   # w_m ~ Aw m^-2
for k in [200, 1000]:
    Tk = np.sum(w_m / (m + k)) + Aw * (1.0 / (k * M) - np.log1p(k / M) / k ** 2)   # int_M^inf Aw/(x^2(x+k)) dx
    B = float(d(k)) / (2 * k - 1) * Tk
    A = float(nu_neg(k))
    print("   k=%5d  partition %.15g  kernel %.15g  rel.diff %.2e  (direct sum)" % (k, A, B, abs(A / B - 1)))

print("3. total mass")
Jp = lambda r: Z * C * (z3 - mp.polylog(3, r) / r) if r > 0 else Z * C * (z3 - 1)
neg_mass = 2 / c * mp.quad(lambda r: Jp(r) / mp.sqrt(1 - r), [0, 0.5, 0.9, 1])
An = nu_m[-1] * M ** 3.5                                # nu(m) ~ An m^-7/2
pos_mass = np.sum(nu_m) + An * (M ** -2.5 / 2.5 - M ** -3.5 / 2)   # Euler-Maclaurin tail
print("   negative part %s, positive part %.15f, total - 1 = %.2e" % (mp.nstr(neg_mass, 15), pos_mass, float(neg_mass) + pos_mass - 1))


def h(j):
    return 2 * j * d(j) if j >= 1 else mp.mpf(0)


print("4. harmonicity of h")
for p in range(1, 7):
    s = mp.nsum(lambda mm: nu_pos(mm) * h(p + mm), [0, mp.inf])
    s += mp.fsum(nu_neg(k) * h(p - k) for k in range(1, p))
    print("   p=%d  sum_j nu(j) h(p+j) / h(p) - 1 = %s" % (p, mp.nstr(s / h(p) - 1, 5)))
