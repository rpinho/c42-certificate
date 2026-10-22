#!/usr/bin/env python3
"""Certify the analytic inequalities for the step-profile Turan construction.

All input decimals are exact rationals. All transcendental operations use
outward-rounded intervals, and every verdict compares exact dyadic endpoints.
The conversion from these inequalities to C_42 is proved in notes/turan.md.
"""
import argparse
import json
from fractions import Fraction as Q
from pathlib import Path

from mpmath import iv
from mpmath.libmp import to_rational

if not __debug__:
    raise RuntimeError("run without -O: assertions are part of certificate verification")


def rational(x):
    if not isinstance(x, (str, int)):
        raise ValueError("certificate numbers must be exact strings or integers")
    return Q(x)


def interval(x):
    x = Q(x)
    return iv.mpf(x.numerator) / iv.mpf(x.denominator)


def endpoints(x):
    return tuple(Q(*to_rational(t)) for t in x._mpi_)


def show(x, digits=30):
    lo, hi = endpoints(x)
    scale = 10**digits
    a = lo.numerator * scale // lo.denominator
    b = -((-hi.numerator * scale) // hi.denominator)
    def dec(k):
        sign = "-" if k < 0 else ""
        k = abs(k)
        return f"{sign}{k//scale}.{k%scale:0{digits}d}"
    return [dec(a), dec(b)]


def error_box(bound):
    _, hi = endpoints(bound)
    e = iv.mpf([-interval(hi).b, interval(hi).b])
    return iv.mpc(e, e)


def power(x, alpha):
    return iv.exp(alpha * iv.log(x))


def integral(alpha, xq, terms):
    x = interval(xq)
    a = alpha.real
    s = iv.mpc(0)
    for k in range(terms):
        s += x**k / (alpha + k)
    tail = x ** (a + terms) / ((a + terms) * (1-x))
    return power(x, alpha) * s + error_box(tail)


def pair_integral(alpha, tq, rq, terms):
    lq = 1-tq-rq
    if lq <= 0:
        return iv.mpc(0)
    t, r, l = map(interval, (tq, rq, lq))
    bt, br = 1-t, 1-r
    a = alpha.real
    c0 = iv.log(bt*br/(t*r))
    c = c0
    s = c/alpha
    for k in range(1, terms):
        c -= (bt**(-k) + br**(-k))/k
        s += c * l**k / (alpha+k)
    tail = l**a/(a+terms) * (
        c0*l**terms/(1-l)
        + (l/bt)**terms/(t*(1-l/bt))
        + (l/br)**terms/(r*(1-l/br))
    )
    return power(l, alpha)*s + error_box(tail)


def certify(data, terms=240):
    c = rational(data["C"])
    alphaq = tuple(map(rational, data["alpha"]))
    cuts = list(map(rational, data["breaks"]))
    eta = [tuple(map(rational, z)) for z in data["eta"]]
    assert len(cuts) == len(eta)+1 and len(eta) > 0
    assert Q(1, 3) <= cuts[0] < Q(1, 2)
    assert cuts[-1] == 1-cuts[0]
    assert all(a < b for a, b in zip(cuts, cuts[1:]))
    assert alphaq[0] > 0 and 0 < c < Q("0.6906538")
    sq = (1-alphaq[0], -alphaq[1])
    margins = [c*c-z[0]*z[0]-z[1]*z[1] for z in [sq]+eta]
    assert min(margins) > 0, "profile exceeds C"
    alpha = iv.mpc(*map(interval, alphaq))
    s = iv.mpc(*map(interval, sq))
    k = integral(alpha, cuts[0], terms)
    dcomplex = integral(iv.mpc(interval(alphaq[0]), 0), cuts[0], terms)
    d = dcomplex.real
    assert endpoints(d)[0] > 0
    y = 1+s*k
    prev = s
    steps = []
    for t, z in zip(cuts, eta):
        new = iv.mpc(*map(interval, z))
        step = new-prev
        steps.append(step)
        y -= step*(integral(alpha, 1-t, terms)-k)
        prev = new
    for i in range(len(steps)):
        for j in range(i+1):
            coefficient = interval(Q(1, 2)) if i == j else interval(1)
            y += coefficient*steps[i]*steps[j]*pair_integral(alpha, cuts[i], cuts[j], terms)
    normsq = y.real*y.real+y.imag*y.imag
    slack = interval(c)**2*d*d-normsq
    assert endpoints(slack)[0] > 0, "limiting inequality fails"
    ratio = iv.sqrt(normsq)/d
    assert endpoints(ratio)[1] < c
    return {
        "certified_upper_bound": str(c),
        "decimal_upper_bound": data["C"],
        "terms": terms,
        "interval_dps": iv.dps,
        "profile_squared_radius_margins": [str(x) for x in margins],
        "D": show(d),
        "Y_real": show(y.real),
        "Y_imag": show(y.imag),
        "abs_Y_over_D": show(ratio),
        "C_squared_D_squared_minus_abs_Y_squared": show(slack),
        "improvement_over_listed_bound": str(Q("0.6906538")-c),
        "status": "PASS: analytic inequalities certified; asymptotic implication in notes/turan.md",
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("certificate", nargs="?", default="certificates/turan.json")
    ap.add_argument("--dps", type=int, default=100)
    ap.add_argument("--terms", type=int, default=240)
    ap.add_argument("--output")
    args = ap.parse_args()
    assert args.dps >= 40 and args.terms >= 40
    iv.dps = args.dps
    result = certify(json.loads(Path(args.certificate).read_text()), args.terms)
    rendered = json.dumps(result, indent=2)+"\n"
    print(rendered, end="")
    if args.output:
        Path(args.output).write_text(rendered)
