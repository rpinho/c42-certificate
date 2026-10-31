# Adapted from Bakeas/turan-c42-bound/independent_check.py (MIT).
# Original license: third_party/LICENSE.bakeas-turan-c42-bound. See BAKEAS_ADAPTATION.md.
# Additional local adaptation: integrated scalar Volterra residual majorant,
# derived in VOLTERRA_VERIFIER_MATH.md. No import from our verifier/search.
# This variant additionally uses the local contraction before rational cutoffs.
# Knot geometry, the untrusted trial solve, threading, receipt output and both
# propagation bounds change. The refusal under python -O and the finite checks
# follow version 1.1.0 of the original.
"""Independent rigorous check of F < T for a C42 profile (python-flint / Arb ball arithmetic).

Written from NOTES.md (Theorem A, Lemma 3, bound (8)) only; it does not reuse rigorous_limit.py.

Usage:  python independent_check.py c42_profile_0.6865.json M      (grid N = 48*M)

Method
------
* rho~ is piecewise linear on the grid t_i = i/N (jump allowed at 1/2). Its nodal values come
  from a floating-point collocation R(t_i) = 0. rho~ may be any function: rigour comes only
  from the a-posteriori residual bound (Lemma 3).
* rho~ is written in "ramp" form
        rho~(v) = r0 + s0*v + sum_k Del_k (v - t_k)_+ + J*1[v > 1/2],
  so that, with D1' = d, D2' = D1 (D1(0)=D2(0)=0), G1(s) = s^(1-g)/(1-g), G2(s) = s^(2-g)/((1-g)(2-g)):
        int_0^t d(t-v) rho~(v) dv     = r0 D1(t) + s0 D2(t) + sum_k Del_k D2(t-t_k) + J D1(t-1/2)
        int_0^t (t-v)^(-g) rho~(v) dv = r0 G1(t) + s0 G2(t) + sum_k Del_k G2(t-t_k) + J G1(t-1/2).
  The sums over k are discrete convolutions (acb_poly products).
* Residual R(t) = d(t) + int_0^t d(u) rho~(t-u) du - t rho~(t).  On each cell:
        sup |R/t| <= (max |R at the two ends| + h^2/8 * sup_cell |R''|) / t_i      (i >= 1),
        sup |R/t| <= |R'(0+)| + h/2 * sup_cell |R''|                               (i = 0),
  where sup_cell |R''| is an interval enclosure (d, d', d'' enclosed on each cell by a
  centre value plus a Lipschitz radius), not a global third-derivative bound.
* kappa = sup|d| = T is verified;  eps = sup|R/t| / (1 - kappa)  (Lemma 3).
* F <= (|p~(1)| + eps/(1-a)) / (Den - eps * I),  I >= int_0^{1/2} t^(1-a)/((1-a)(1-t)) dt,
  bounded with 1/(1-t) <= 1 + 2t on [0, 1/2].  Den is a rigorous lower bound for
  int_0^{1/2} |p~(t)| dt/(1-t):  near 0 by |p~| >= t^-a - P* t^(1-a)/(1-a); on the other cells by
  |int L/(1-t)| - (h^2/8) int|p~''| / (1 - t_{i+1}),  L = linear interpolant of p~.
"""
import json
import hashlib
from pathlib import Path
from fractions import Fraction
from scipy.signal import fftconvolve
import math
import sys
import time

import numpy as np
from flint import acb, acb_poly, arb, arb_poly, ctx

if not __debug__:   # the checks below are assert statements, which python -O would remove
    sys.exit("run this program without python -O")

PREC = 192
NPIECES = 24          # knots j/48, j = 0..24, cover [0, 1/2]


def series_f(z, c, terms=30):
    """sum_{n>=0} z^n/(n+c)!  as a rigorous ball, for |z| <= 1/2."""
    s, zn, f = acb(0), acb(1), arb.fac_ui(c)
    for n in range(terms):
        s += zn / f
        zn *= z
        f *= n + 1 + c
    az = z.abs_upper()
    assert az < arb(1) / 2
    tail = (zn.abs_upper() / (f * (1 - az))).upper()      # sum_{n>=terms} |z|^n/(n+c)!
    return s + acb(arb(0, tail), arb(0, tail))


def box(r):
    """Complex ball centred at 0 that contains the disk of radius r."""
    r = r.upper()
    return acb(arb(0, r), arb(0, r))


def to_c(x):
    return complex(float(x.real.mid()), float(x.imag.mid()))


def main():
    path = sys.argv[1]
    M = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    ctx.prec = PREC
    ctx.threads = 2
    t0 = time.time()
    code_sha=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    certificate_raw=Path(path).read_bytes()
    prof = json.loads(certificate_raw)
    T = arb(prof["T"])
    psi = [arb(s) for s in prof["phases"]]
    NPIECES=len(psi)-1
    assert len(psi) == NPIECES + 1
    knots = list(map(Fraction, prof['breaks'])) if isinstance(prof.get('breaks'), list) else [Fraction(j,48) for j in range(25)]
    assert len(knots)==len(psi) and knots[0]==0 and knots[-1]==Fraction(1,2)
    claim_q=Fraction(prof.get('C',prof['T']))
    assert Fraction(prof['T'])<=claim_q<1
    claim=arb(str(claim_q))
    assert all(x<y for x,y in zip(knots,knots[1:]))
    beta = [(psi[j+1]-psi[j])/arb(str(knots[j+1]-knots[j])) for j in range(NPIECES)]
    I = acb(0, 1)
    gam = T * (I * psi[0]).exp()
    a = gam.real

    # ---- hypotheses of Theorem A (checked rigorously) -------------------------------------
    # On (0,1/2]: |d| = T|1 - e^{i(phi-psi0)}| <= T * max_j |psi_j - psi_0|  (phi piecewise linear).
    dev = max(abs(p - psi[0]).abs_upper() for p in psi)
    assert dev < arb(1), "need max|psi_j - psi_0| < 1 so that kappa = T"
    assert T>0 and T < arb(1), "need 0<T < 1 (kappa = T < 1)"
    assert a > arb(0) and a < arb(1), "need 0 < a < 1"
    kappa = T                       # |d| = T on (1/2, 1], and |d| < T on (0, 1/2]

    N = 48 * M
    H = N // 2
    assert all((x*N).denominator==1 for x in knots)
    starts = [int(x*N) for x in knots]
    h = arb(1) / N
    t = [arb(i) / N for i in range(N + 1)]
    bmax = max(abs(b).abs_upper() for b in beta)

    # ---- per-piece and per-cell data for d ---------------------------------------------------
    f1 = [series_f(I * beta[j] * h, 1) for j in range(NPIECES)]
    f2 = [series_f(I * beta[j] * h, 2) for j in range(NPIECES)]
    m0 = [None] * N             # int_0^h d(t_m + x)(1 - x/h) dx   (only for the float solve)
    m1 = [None] * N             # int_0^h d(t_m + x)(x/h) dx
    D1 = [acb(0)] * (N + 1)     # D1(t_m)
    D2 = [acb(0)] * (N + 1)     # D2(t_m)
    dc0, dc1, dc2 = [None] * N, [None] * N, [None] * N      # d, d', d'' at cell centres
    r0, r1, r2 = [arb(0)] * N, [arb(0)] * N, [arb(0)] * N   # Lipschitz radii on the cell
    dplus = [None] * (N + 1)    # d(t_i^+)
    dminus = [None] * (N + 1)   # d(t_i^-)
    for m in range(N):
        if m < H:
            j = int(np.searchsorted(starts,m,side="right")-1)
            th = psi[j] + beta[j] * ((m - starts[j]) * h)
            E = (I * th).exp()
            Ec = (I * (th + beta[j] * h / 2)).exp()
            intd = gam * h - T * E * h * f1[j]
            intd2 = gam * h * h / 2 - T * E * h * h * f2[j]      # int_0^h (h-x) d(t_m+x) dx
            m0[m] = intd2 / h
            m1[m] = intd - m0[m]
            b = beta[j]
            dc0[m] = gam - T * Ec
            dc1[m] = -I * T * b * Ec
            dc2[m] = T * b * b * Ec
            ab = abs(b)
            r0[m] = T * ab * h / 2
            r1[m] = T * ab * ab * h / 2
            r2[m] = T * ab * ab * ab * h / 2
            dplus[m] = gam - T * E
            thr = psi[j] + beta[j] * ((m + 1 - starts[j]) * h)     # phase at t_{m+1}, same piece
            dminus[m + 1] = gam - T * (I * thr).exp()
        else:
            intd = gam * h
            intd2 = gam * h * h / 2
            m0[m] = intd2 / h
            m1[m] = intd - m0[m]
            dc0[m], dc1[m], dc2[m] = gam, acb(0), acb(0)
            dplus[m] = gam
            dminus[m + 1] = gam
        D1[m + 1] = D1[m] + intd
        D2[m + 1] = D2[m] + h * D1[m] + intd2
    t1 = time.time()

    # ---- floating-point collocation for rho~ (any rho~ is allowed) ------------------------
    m0f = np.array([to_c(x) for x in m0])
    m1f = np.array([to_c(x) for x in m1])
    dmf = np.array([0j] + [to_c(dminus[i]) for i in range(1, N + 1)])
    gf, Tf = to_c(gam), float(T.mid())
    # ADAPTATION: fast trial generation only. All rigorous bounds below are unchanged.
    jump = 2*float(T.mid())*np.exp(1j*float(psi[-1].mid()))
    origin = -1j*float(beta[0].mid())*gf
    times = np.arange(N+1)/N
    rhs = dmf.copy(); rhs[1:] += m1f*origin
    d1f = np.array([to_c(z) for z in D1])
    rhs[H+1:] += jump*(d1f[1:H+1]-times[H+1:])
    weights = np.zeros(N+1,complex); weights[1:N] = m0f[1:]+m1f[:-1]
    diag = times-m0f[0]
    vv = np.zeros(N+1,complex); vv[0] = origin
    def solve(lo,hi):
        if hi-lo<=48:
            for i in range(lo,hi):
                if i>lo: rhs[i] += vv[lo:i]@weights[i-lo:0:-1]
                vv[i] = rhs[i]/diag[i]
            return
        mid=(lo+hi)//2;solve(lo,mid)
        cv=fftconvolve(vv[lo:mid],weights[:hi-lo])
        rhs[mid:hi]+=cv[mid-lo:hi-lo];solve(mid,hi)
    solve(1,N+1)
    rm=vv.copy();rm[H+1:]+=jump
    rp=vv.copy();rp[H:]+=jump
    t2 = time.time()

    # ---- exact binary values of rho~ and its ramp coefficients -----------------------------
    RP = [acb(complex(z)) for z in rp[:N]]        # exact (doubles)
    RM = [acb(0)] + [acb(complex(z)) for z in rm[1:]]
    for i in range(1, N):
        if i != H:
            assert rm[i] == rp[i]
    sig = [(RM[i + 1] - RP[i]) * N for i in range(N)]          # slope on cell i
    Del = [acb(0)] + [sig[k] - sig[k - 1] for k in range(1, N)]
    J = RP[H] - RM[H]
    r0t = RP[0]
    Apoly = acb_poly(Del)

    # ---- (1) nodal residuals ------------------------------------------------------------
    C = (Apoly * acb_poly(D2)).coeffs()
    C += [acb(0)] * (N + 1 - len(C))
    Rm = [None] * (N + 1)
    Rp = [None] * (N + 1)
    for i in range(1, N + 1):
        base = r0t * D1[i] + sig[0] * D2[i] + C[i]
        if i > H:
            base += J * D1[i - H]
        Rm[i] = dminus[i] + base - t[i] * RM[i]
        if i < N:
            Rp[i] = dplus[i] + base - t[i] * RP[i]

    # ---- (2) enclosure of R'' on each cell -------------------------------------------------
    #   R'' = d'' + r0 d' + s0 d + sum_{k<=i} Del_k d(t - t_k) + [i>=H] J d'(t - 1/2) - 2 sig_i
    Dcell = [dc0[m] + box(r0[m]) for m in range(N)]
    B2 = (Apoly * acb_poly(Dcell)).coeffs()
    B2 += [acb(0)] * (N - len(B2))
    ar0, as0, aJ = abs(r0t).abs_upper(), abs(sig[0]).abs_upper(), abs(J).abs_upper()
    U = [None] * N
    for i in range(N):
        val = dc2[i] + r0t * dc1[i] + sig[0] * dc0[i] + B2[i] - 2 * sig[i]
        rad = r2[i] + ar0 * r1[i] + as0 * r0[i]
        if i >= H:
            val += J * dc1[i - H]
            rad += aJ * r1[i - H]
        U[i] = (val.abs_upper() + rad).upper()

    # ---- (3) sup |R/t| cell by cell ------------------------------------------------------
    Rd0 = -I * beta[0] * gam - r0t                # R'(0+) = d'(0+) - rho~(0)
    res = [(Rd0.abs_upper() + h / 2 * U[0]).upper()]
    for i in range(1, N):
        e = max(Rp[i].abs_upper(), Rm[i + 1].abs_upper())
        res.append(((e + h * h / 8 * U[i]) / t[i]).upper())
    assert all(r.is_finite() for r in res), "a cell bound is not finite"
    imax = max(range(N), key=lambda i: res[i])
    supR = res[imax]
    eps = (supR / (1 - kappa)).upper()
    # Independent implementation of the integrated-residual majorant.
    # A(x)=int_0^x res(v)*v^(-T)dv bounds the integrated rho error.
    Awhole=arb(0);Ahalf=arb(0);Rhalf=arb(0);Rnum=arb(0)
    oldA=arb(0);oldN=arb(1)
    for cell in range(N):
        newA=t[cell+1]**(1-T)
        newN=(1-t[cell+1])**(1-a) if cell+1<N else arb(0)
        area=res[cell]*(newA-oldA)/(1-T)
        Awhole+=area
        Rnum+=res[cell]*(oldN-newN)/(1-a)
        if cell<H:
            Ahalf+=area;Rhalf+=h*res[cell]
        oldA=newA;oldN=newN
    betaN=T.gamma()*(1-a).gamma()/(T+1-a).gamma()
    betaD=T.gamma()*(2-a).gamma()/(T+2-a).gamma()
    new_num_error=(Rnum+T*Awhole*betaN).upper()
    new_den_error=(2/(1-a)*(arb('0.5')**(1-a)*Rhalf+T*Ahalf*arb('0.5')**(T+1-a)*betaD)).upper()
    local_errors=[]
    for cut_q in (Fraction(1,256),Fraction(1,64),Fraction(1,16)):
        if (cut_q*N).denominator!=1:continue
        exact_psi=list(map(Fraction,prof['phases']))
        excursions=[abs(p-exact_psi[0]) for p,u in zip(exact_psi,knots) if u<=cut_q]
        for j in range(NPIECES):
            if knots[j]<=cut_q<=knots[j+1]:
                interpolated=exact_psi[j]+(exact_psi[j+1]-exact_psi[j])*(cut_q-knots[j])/(knots[j+1]-knots[j])
                excursions.append(abs(interpolated-exact_psi[0]));break
        k_q=Fraction(prof['T'])*max(excursions)
        assert 0<=k_q<Fraction(prof['T'])
        cut=arb(str(cut_q));small=arb(str(k_q));count=int(cut_q*N)
        local_weighted=arb(0);early_weighted=arb(0);early_plain=arb(0);early_num=arb(0)
        old_small=arb(0);old_global=arb(0);old_num=arb(1)
        for j in range(count):
            u=t[j+1]**(1-small);v=t[j+1]**(1-T);w=(1-t[j+1])**(1-a)
            local_weighted+=res[j]*(u-old_small)/(1-small)
            early_weighted+=res[j]*(v-old_global)/(1-T)
            early_num+=res[j]*(old_num-w)/(1-a);early_plain+=h*res[j]
            old_small=u;old_global=v;old_num=w
        Ecut=cut**small*local_weighted
        effective_whole=cut**(-T)*Ecut+Awhole-early_weighted
        effective_half=cut**(-T)*Ecut+Ahalf-early_weighted
        num_cut=((1-cut)**(-a)*Ecut+Rnum-early_num+T*effective_whole*betaN).upper()
        den_cut=(2/(1-a)*(arb('0.5')**(1-a)*(Ecut+Rhalf-early_plain)+T*effective_half*arb('0.5')**(T+1-a)*betaD)).upper()
        local_errors.append(dict(delta=str(cut_q),kappa=str(k_q),E_delta=Ecut.str(20),numerator_error=num_cut.str(20),denominator_error=den_cut.str(20)))
        new_num_error=min(new_num_error,num_cut);new_den_error=min(new_den_error,den_cut)
    t3 = time.time()

    # ---- (4) p~ : numerator and nodal values ------------------------------------------------
    one_g, two_g = 1 - gam, 2 - gam
    tpow = [acb(0)] + [(-gam * t[m].log()).exp() for m in range(1, N + 1)]    # t_m^(-g)
    G2 = [acb(0)] + [t[m] * t[m] * tpow[m] / (one_g * two_g) for m in range(1, N + 1)]
    C3 = (Apoly * acb_poly(G2)).coeffs()
    C3 += [acb(0)] * (N + 1 - len(C3))
    half = arb(1) / 2
    p1 = 1 + r0t / one_g + sig[0] / (one_g * two_g) + C3[N] \
        + J * (-gam * half.log()).exp() * half / one_g
    num = (abs(p1).abs_upper() + min(new_num_error,(eps/(1-a)).upper())).upper()

    def ptil(i):     # p~(t_i), 1 <= i <= H
        return tpow[i] + r0t * t[i] * tpow[i] / one_g + sig[0] * G2[i] + C3[i]

    # ---- (5) denominator lower bound ----------------------------------------------------
    istar = max(2, math.ceil(0.25 * N ** (1 / 3)))
    ts = t[istar]
    Pst = max([abs(RP[0]).abs_upper()] + [abs(RM[i]).abs_upper() for i in range(1, istar + 1)])
    den0 = ts ** (1 - a) / (1 - a) - Pst * ts ** (2 - a) / ((1 - a) * (2 - a))
    absg2 = abs(gam * (gam + 1)).abs_upper()
    absgr = abs(gam * r0t).abs_upper()
    hh = h ** (1 - a) / (1 - a)
    cm = [hh * (arb(m + 1) ** (1 - a) - arb(m) ** (1 - a)) for m in range(H + 1)]
    C4 = (arb_poly([abs(x) for x in Del[:H + 1]]) * arb_poly(cm)).coeffs()
    C4 += [arb(0)] * (H + 1 - len(C4))
    den = den0
    pcur = ptil(istar)
    for i in range(istar, H):
        pnext = ptil(i + 1)
        ell = (h / (1 - t[i + 1])).log1p()                  # int_cell dt/(1-t)
        intLw = pcur * ell + (pnext - pcur) * ((1 - t[i]) * ell - h) / h
        ipp = absg2 * (t[i] ** (-a - 1) - t[i + 1] ** (-a - 1)) / (a + 1) \
            + absgr * (t[i] ** (-a) - t[i + 1] ** (-a)) / a \
            + as0 * (t[i + 1] ** (1 - a) - t[i] ** (1 - a)) / (1 - a) \
            + C4[i]
        den += abs(intLw) - h * h / 8 * ipp / (1 - t[i + 1])
        pcur = pnext
    den_p = den.lower()
    Iw = (half ** (2 - a) / (2 - a) + 2 * half ** (3 - a) / (3 - a)) / (1 - a)
    den_lo = (den - min(new_den_error,(eps*Iw).upper())).lower()
    assert num.is_finite() and den_lo.is_finite()
    assert den_lo > 0
    Fup = (num / den_lo).upper()
    t4 = time.time()

    ftil = float(abs(p1).mid()) / float(den.mid())
    print(f"profile        : {path}")
    print(f"M = {M}, N = {N}, prec = {PREC} bits, i* = {istar}")
    print(f"a = Re(gamma)  : {a.str(10)}   kappa = T = {T.str(6)}")
    print(f"sup|R/t|       : {float(supR):.4e}   (worst cell t = {imax / N:.6f})")
    print(f"eps            : {float(eps):.4e}")
    print(f"p~(1)          : {p1.str(20, radius=False)}")
    print(f"|p(1)| <=      : {float(num):.7f}")
    print(f"int|p|/(1-t) >=: {float(den_lo):.7f}   (for p~ alone: {float(den_p):.7f})")
    print(f"F~ (approx)    : {ftil:.8f}")
    print(f"F <=           : {Fup.str(10, radius=False)}")
    print(f"time           : setup {t1 - t0:.1f}s, float solve {t2 - t1:.1f}s, "
          f"residual {t3 - t2:.1f}s, F bound {t4 - t3:.1f}s, total {t4 - t0:.1f}s")
    if len(sys.argv)>3:
        receipt=dict(verdict='PASS' if Fup<claim else 'UNPROVED',C=str(claim_q),T=prof['T'],N=N,bits=PREC,threads=2,
                     F_upper=Fup.str(30),F_upper_rational=str(Fup.fmpq()),p_tilde_1=p1.str(40),
                     residual_upper=supR.str(16),rho_error_upper=eps.str(16),
                     integrated_residual=dict(A1=Awhole.str(20),Ah=Ahalf.str(20),RN=Rnum.str(20),Rh=Rhalf.str(20),
                                              numerator_error=new_num_error.str(20),denominator_error=new_den_error.str(20),local_contractions=local_errors),
                     numerator_upper=num.str(30),denominator_lower=den_lo.str(30),seconds=t4-t0,
                     certificate_sha256=hashlib.sha256(certificate_raw).hexdigest(),
                     verifier_sha256=code_sha,
                     adaptation='Bakeas independent_check: rational knots, fast trial solve, original residual/denominator enclosures, integrated residual plus local contraction')
        Path(sys.argv[3]).write_text(json.dumps(receipt,indent=2)+'\n')
    if Fup < claim:
        print("PASS  (max(T,F) <= C certified)")
    else:
        print("FAIL  (bound not below T)")
        sys.exit(1)


if __name__ == "__main__":
    main()
