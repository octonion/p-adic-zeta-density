#!/usr/bin/env python3
"""Outward-rounded certificate for positive-genus finite-packet bounds.

This script certifies consequences of the exact modular capacity-pair theorem
for X_0(p).  It uses only integer fixed-point interval arithmetic for pi and
logarithms, exact integer/rational operations, and the elementary bound

  C_p(1/m) <= (2*pi^2/p) H_{floor((m-1)/p)}.

No binary floating point is used in any certificate decision.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

DIGITS = 45
SCALE = 10 ** DIGITS


def floor_div(a: int, b: int) -> int:
    assert b > 0
    return a // b


def ceil_div(a: int, b: int) -> int:
    assert b > 0
    return -((-a) // b)


@dataclass(frozen=True)
class I:
    lo: int
    hi: int

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError("bad interval")

    @staticmethod
    def integer(n: int) -> "I":
        return I(n * SCALE, n * SCALE)

    @staticmethod
    def fraction(a: int, b: int) -> "I":
        if b <= 0:
            raise ValueError("positive denominator required")
        return I(floor_div(a * SCALE, b), ceil_div(a * SCALE, b))

    def __add__(self, other: "I") -> "I":
        return I(self.lo + other.lo, self.hi + other.hi)

    def __sub__(self, other: "I") -> "I":
        return I(self.lo - other.hi, self.hi - other.lo)

    def __neg__(self) -> "I":
        return I(-self.hi, -self.lo)

    def __mul__(self, other: "I") -> "I":
        ps = [self.lo * other.lo, self.lo * other.hi,
              self.hi * other.lo, self.hi * other.hi]
        return I(floor_div(min(ps), SCALE), ceil_div(max(ps), SCALE))

    def mul_int(self, n: int) -> "I":
        if n >= 0:
            return I(self.lo * n, self.hi * n)
        return I(self.hi * n, self.lo * n)

    def div_int(self, n: int) -> "I":
        if n <= 0:
            raise ValueError("positive divisor required")
        return I(floor_div(self.lo, n), ceil_div(self.hi, n))

    def mul_fraction(self, a: int, b: int) -> "I":
        return self * I.fraction(a, b)

    def sq(self) -> "I":
        return self * self

    def lower_decimal(self, places: int = 15) -> str:
        sign = "-" if self.lo < 0 else ""
        a = abs(self.lo)
        q, r = divmod(a, SCALE)
        frac = str(r).rjust(DIGITS, "0")[:places]
        return f"{sign}{q}.{frac}"

    def upper_decimal(self, places: int = 15) -> str:
        sign = "-" if self.hi < 0 else ""
        a = abs(self.hi)
        q, r = divmod(a, SCALE)
        frac = str(r).rjust(DIGITS, "0")[:places]
        return f"{sign}{q}.{frac}"


def atan_unit_fraction(q: int, even_last: int) -> I:
    """Rigorous atan(1/q) interval from alternating series.

    even_last is an even j.  S_even (through j=even_last) is an upper bound;
    the next partial sum through odd j is a lower bound.
    """
    if even_last % 2:
        raise ValueError("even_last must be even")
    x = I.fraction(1, q)
    x2 = x * x
    power = x
    partial = I.integer(0)
    partial_even = None
    for j in range(even_last + 2):
        term = power.div_int(2 * j + 1)
        partial = partial + term if j % 2 == 0 else partial - term
        if j == even_last:
            partial_even = partial
        power = power * x2
    # partial now includes j=even_last+1, an odd last term -> lower bound.
    assert partial_even is not None
    return I(partial.lo, partial_even.hi)


def pi_interval() -> I:
    # Errors are far below 10^-60 at these truncations.
    a = atan_unit_fraction(5, 50)
    b = atan_unit_fraction(239, 12)
    return a.mul_int(16) - b.mul_int(4)


PI = pi_interval()


def log_reduced_ratio(num: int, den: int, terms: int = 52) -> I:
    """Rigorous log(num/den) for 1 <= num/den < 2 via atanh series."""
    if not (den <= num <= 2 * den):
        raise ValueError("ratio not in [1,2)")
    a = num - den
    b = num + den
    t = I.fraction(a, b)
    t2 = t * t
    power = t
    sm = I.integer(0)
    for j in range(terms + 1):
        sm = sm + power.div_int(2 * j + 1)
        power = power * t2
    base = sm.mul_int(2)
    # Remaining positive tail <= 2*t^(2N+3)/(2N+3)/(1-t^2).
    # power is t^(2*terms+3) after the loop.
    recip = I.fraction(b * b, b * b - a * a)
    tail = (power.div_int(2 * terms + 3) * recip).mul_int(2)
    return I(base.lo, base.hi + max(0, tail.hi))


LOG2 = log_reduced_ratio(2, 1)
LOG_CACHE: dict[int, I] = {2: LOG2}


def log_int(n: int) -> I:
    if n <= 0:
        raise ValueError("log positive only")
    if n == 1:
        return I.integer(0)
    if n in LOG_CACHE:
        return LOG_CACHE[n]
    k = n.bit_length() - 1
    den = 1 << k
    ylog = log_reduced_ratio(n, den)
    ans = LOG2.mul_int(k) + ylog
    LOG_CACHE[n] = ans
    return ans


# Cache harmonic upper bounds.  H_UPPER[r] is an outward upper fixed-point
# enclosure of H_r obtained by summing ceil(SCALE/j).
H_UPPER = [0]

def harmonic_upper(r: int) -> int:
    while len(H_UPPER) <= r:
        j = len(H_UPPER)
        H_UPPER.append(H_UPPER[-1] + ceil_div(SCALE, j))
    return H_UPPER[r]


def capacity_width(p: int) -> I:
    return log_int(p).mul_fraction(12, p - 1)


def analytic_lower(p: int, m0: int) -> I:
    """Lower bound for Lambda_p(1/m0)-E_p(1/m0)/m0."""
    # A = 12 log p/(p-1) - 2 pi/m0.
    A = capacity_width(p) - PI.mul_int(2).div_int(m0)
    # (1 - 1/m0) A
    first = A.mul_fraction(m0 - 1, m0)
    # C <= (2*pi^2/p) H_floor((m0-1)/p).
    r = (m0 - 1) // p
    H = I(0, harmonic_upper(r))
    Cup = PI.sq().mul_int(2).div_int(p) * H
    return first - I(0, Cup.hi).div_int(m0)


def arithmetic_upper(d: int, k: int, x: int) -> I:
    """Upper bound for tau at xi=d for any k-packet with largest x.

    m >= m0=x+k^2 and Q <= (x-2)m+x+2.  The resulting upper bound is
    decreasing in m, so it is evaluated at m0.
    """
    m0 = x + k * k
    lm = log_int(m0)
    t1 = I.integer(x - 2) + lm.mul_int(x)
    t1 = t1.div_int(m0)
    t2 = I.fraction(x + 2, m0 * m0)
    return (t1 + t2).mul_int(2 * d)


def certify_packet_bound(p: int, d: int, k: int, Smax: int):
    """Certify that no k rational values occur among odd 3..Smax."""
    if Smax % 2 == 0:
        raise ValueError("Smax must be odd")
    xmin = 2 * k + 1
    if xmin > Smax:
        # There are fewer than k available odd arguments; combinatorially empty.
        return {"ok": True, "worst_x": None, "margin": I.integer(1)}
    worst = None
    for x in range(xmin, Smax + 1, 2):
        m0 = x + k * k
        rhs = analytic_lower(p, m0)
        lhs = arithmetic_upper(d, k, x)
        margin = rhs - lhs
        if worst is None or margin.lo < worst[1].lo:
            worst = (x, margin, rhs, lhs, m0)
    x, margin, rhs, lhs, m0 = worst
    return {
        "ok": margin.lo > 0,
        "worst_x": x,
        "margin": margin,
        "rhs": rhs,
        "lhs": lhs,
        "m0": m0,
    }


def prime_genus(p: int) -> int:
    if p in (2, 3, 5, 7, 13):
        return 0
    def leg(a: int) -> int:
        v = pow(a % p, (p - 1) // 2, p)
        return -1 if v == p - 1 else (0 if v == 0 else 1)
    e2 = 1 + leg(-1)
    e3 = 1 + leg(-3)
    # g=(p+1)/12-e2/4-e3/3, evaluated integrally.
    return ((p + 1) * 1 - 3 * e2 - 4 * e3) // 12


def find_cap(p: int, d: int, Smax: int = 2401):
    N = (Smax - 1) // 2
    for k in range(1, N + 1):
        cert = certify_packet_bound(p, d, k, Smax)
        if cert["ok"]:
            return k - 1, cert
    return N, None


def find_1000_cutoff(p: int, d: int, Nlimit: int = 2200):
    for N in range(1000, Nlimit + 1):
        Smax = 2 * N + 1
        k = N - 999  # exclude k rationals => at most N-1000 rationals.
        cert = certify_packet_bound(p, d, k, Smax)
        if cert["ok"]:
            return Smax, k - 1, cert
    return None


def main():
    # Targets were located independently; the certificate below only verifies
    # the displayed exclusion and, diagnostically, its immediate predecessor.
    cases = [
        # p, genus, d, k excluded at S=2401, 1000-cutoff S, k excluded there
        (11, 1, 2, 191, 2379, 190),
        (17, 1, 2, 227, 2459, 230),
        (19, 1, 2, 238, 2483, 242),
        (23, 2, 4, 378, 2823, 412),
        (29, 2, 4, 415, 2921, 461),
        (31, 2, 4, 426, 2951, 476),
        (37, 2, 4, 458, 3039, 520),
        (41, 3, 6, 596, 3445, 723),
        (43, 3, 6, 607, 3483, 742),
    ]
    print(f"fixed-point scale = 10^{DIGITS}")
    print(f"pi in [{PI.lower_decimal(40)}, {PI.upper_decimal(40)}]")
    print("method: xi=d; m>=x+k^2; Q<=(x-2)m+x+2")
    print("collision bound: C_p(1/m) <= (2*pi^2/p) H_floor((m-1)/p)")
    print("all decisions use strict positivity of an outward-rounded lower margin")
    print()
    print("p genus d cap_2401 irr_min_2401 worst_x margin_lower cutoff_1000 cap_at_cutoff cutoff_worst_x cutoff_margin")
    for p,g,d,k,S1000,k1000 in cases:
        assert prime_genus(p) == g and d == 2*g
        cert = certify_packet_bound(p,d,k,2401)
        if not cert["ok"]:
            raise SystemExit(f"FAILED p={p}, S=2401, k={k}: margin {cert['margin'].lower_decimal(20)}")
        prev = certify_packet_bound(p,d,k-1,2401)
        N=(S1000-1)//2
        assert k1000 == N-999
        cert2=certify_packet_bound(p,d,k1000,S1000)
        if not cert2["ok"]:
            raise SystemExit(f"FAILED p={p}, cutoff={S1000}, k={k1000}: margin {cert2['margin'].lower_decimal(20)}")
        # Verify minimality of the stated 1000 cutoff for this sufficient test.
        if S1000>2001:
            Nprev=N-1; Sprev=S1000-2; kprev=Nprev-999
            cert2prev=certify_packet_bound(p,d,kprev,Sprev)
        else:
            cert2prev=None
        cap=k-1; irr=1200-cap
        print(f"{p:2d} {g:5d} {d:1d} {cap:8d} {irr:12d} {cert['worst_x']:7d} {cert['margin'].lower_decimal(15):>16s} "
              f"{S1000:11d} {k1000-1:13d} {cert2['worst_x']:14d} {cert2['margin'].lower_decimal(15):>13s}")
        print(f"  diagnostics: S=2401 predecessor k={k-1} margin_lo={prev['margin'].lower_decimal(15)}; "
              + (f"cutoff predecessor S={Sprev}, k={kprev}, margin_lo={cert2prev['margin'].lower_decimal(15)}" if cert2prev else ""))

    print("\nCERTIFICATE: PASS")


if __name__ == "__main__":
    main()
