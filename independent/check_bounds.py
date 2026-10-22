#!/usr/bin/env python3
"""Independent FLINT/Arb checkers; written from mathematical notes only.

No floating-point search output is trusted. Inputs are exact rational strings.
C42 uses hypergeometric incomplete integrals and an explicit J-series tail.
"""
import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys

from flint import arb, acb, fmpq, ctx


def rational(s):
    if not isinstance(s, (str, int, Fraction)):
        raise ValueError("exact rational strings or integers required")
    return Fraction(s)


def ball(s):
    x = rational(s)
    return arb(fmpq(x.numerator, x.denominator))


def ibeta(alpha, x):
    """Integral_0^x v^(alpha-1)/(1-v) dv, by hypergeometric identity."""
    x = acb(x)
    return x**alpha / alpha * x.hypgeom_2f1(1, alpha, alpha + 1)


def j_integral(alpha, t, r, terms):
    """Original triangular integral, via independently bounded log series.

    c_k = log((1-t)(1-r)/(tr)) - sum_{q=1}^k ((1-t)^-q+(1-r)^-q)/q.
    Its absolute value is <= c0 + k*((1-t)^-k+(1-r)^-k).
    We sum the resulting k-weighted geometric tails explicitly.
    """
    if t + r >= 1:
        return acb(0)
    t, r = ball(t), ball(r)
    L, b, d = 1-t-r, 1-t, 1-r
    c = (b*d/(t*r)).log()
    c0 = c
    value = acb(0)
    p, ib, id_ = arb(1), arb(1), arb(1)
    for k in range(terms):
        if k:
            ib /= b
            id_ /= d
            c -= (ib + id_) / k
        value += c*p/(alpha + k)
        p *= L
    value *= acb(L)**alpha
    tail = c0 * L**terms/(1-L)
    for q in (L/b, L/d):
        tail += q**terms * (terms/(1-q) + q/(1-q)**2)
    tail *= L**alpha.real / (alpha.real + terms)
    if not tail > 0:
        raise ArithmeticError("J tail positivity failed")
    error = arb(0, tail.upper())
    return value + acb(error, error)


def check_c42(data, terms):
    Cq = rational(data["C"])
    aq, bq = map(rational, data["alpha"])
    ts = list(map(rational, data["breaks"]))
    es = [(rational(x), rational(y)) for x, y in data["eta"]]
    if not (0 < Cq < 1 and aq > 0 and len(es) == len(ts)-1):
        raise ValueError("invalid C42 dimensions or parameters")
    if not (Fraction(1,3) <= ts[0] < Fraction(1,2)
            and ts[-1] == 1-ts[0]
            and all(x < y for x,y in zip(ts, ts[1:]))):
        raise ValueError("C42 requires 1/3 <= tau < 1/2 and B=1-tau")
    exact_margins = [Cq*Cq - ((1-aq)**2 + bq*bq)]
    exact_margins += [Cq*Cq - x*x-y*y for x,y in es]
    if not all(m > 0 for m in exact_margins):
        raise ArithmeticError("exact profile radius inequality fails")
    alpha = acb(ball(aq), ball(bq))
    s = 1-alpha
    eta = [acb(ball(x), ball(y)) for x,y in es]
    delta = [eta[0]-s] + [eta[j]-eta[j-1] for j in range(1,len(eta))]
    tau = ball(ts[0])
    K = ibeta(alpha,tau)
    Dcomplex = ibeta(acb(ball(aq)),tau)
    if not Dcomplex.imag.contains(0):
        raise ArithmeticError("real normalization developed imaginary part")
    D = Dcomplex.real
    if not D > 0:
        raise ArithmeticError("normalization is not proven positive")
    Y = 1+s*K
    for t, e in zip(ts[:-1],delta):
        Y -= e*(ibeta(alpha,ball(1-t))-K)
    for i,ti in enumerate(ts[:-1]):
        for j,tj in enumerate(ts[:i+1]):
            mult = 1 if i != j else Fraction(1,2)
            Y += ball(mult)*delta[i]*delta[j]*j_integral(alpha,ti,tj,terms)
    margin = ball(Cq)**2 * D**2 - abs(Y)**2
    if not margin > 0:
        raise ArithmeticError("C42 analytic inequality not certified: " + str(margin))
    return {"bound":str(Cq), "D":str(D), "Y":str(Y),
            "ratio":str(abs(Y)/D), "analytic_margin":str(margin),
            "smallest_exact_profile_margin":str(min(exact_margins)),
            "terms":terms}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("certificate",type=Path)
    p.add_argument("--bits",type=int,default=512)
    p.add_argument("--terms",type=int,default=320)
    p.add_argument("--out",type=Path)
    a=p.parse_args()
    ctx.prec=a.bits
    input_bytes=a.certificate.read_bytes()
    data=json.loads(input_bytes)
    try:
        if data["constant"]=="C_42": result=check_c42(data,a.terms)
        else: raise ValueError("unknown constant")
        result.update(status="PASS",bits=a.bits,pid=os.getpid(),
                      certificate=str(a.certificate),
                      certificate_sha256=hashlib.sha256(input_bytes).hexdigest(),
                      implementation="independent python-flint/Arb")
        status=0
    except Exception as e:
        result={"status":"FAIL", "reason":str(e), "pid":os.getpid()}
        status=1
    output=json.dumps(result,indent=2)+"\n"
    if a.out: a.out.write_text(output)
    print(output,end="")
    return status


if __name__=="__main__": sys.exit(main())
