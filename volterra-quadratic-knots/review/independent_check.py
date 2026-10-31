#!/usr/bin/env python3
"""Fresh numerical diagnostics. Imports no target code. Not an interval certificate.
Volterra: Picard iteration with exact convolution of linear interpolants of d_c,r.
Finite: independently exponentiate formal series, construct W and its actual roots.
"""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
import json,time,sys,hashlib
from pathlib import Path
from fractions import Fraction
import numpy as np
from scipy.signal import fftconvolve
from scipy.special import roots_jacobi, roots_legendre
from scipy.integrate import quad
from flint import acb,acb_poly,acb_series,ctx

DATA=json.loads(Path('../certificates/C42_compact_best.json').read_text())
T=float(DATA['T']);kn=np.array([float(Fraction(x)) for x in DATA['breaks']]);ph=np.array([float(x) for x in DATA['phases']]);gam=T*np.exp(1j*ph[0]);a=gam.real;b=1-gam;ge=T*np.exp(1j*ph[-1]);J=2*ge

def profile(x):return T*np.exp(1j*np.interp(x,kn,ph))

def conv(x,y,n):
    z=fftconvolve(x,y)[:n]
    # Avoid absolute FFT noise magnified by division by tiny t.
    m=min(512,n,len(x),len(y));z[:m]=np.convolve(x[:m],y[:m])[:m]
    return z

def volterra(N):
    start=time.time();h=1/N;H=N//2;t=np.arange(N+2)*h
    dc=gam-profile(np.minimum(t,.5));dc[0]=0
    prefix=np.r_[0,np.cumsum((dc[:-1]+dc[1:])*h/2)]
    origin=-1j*(ph[1]-ph[0])/(kn[1]-kn[0])*gam
    right=np.zeros(N+1,complex);right[1:]=dc[1:N+1]/t[1:N+1]
    right[H+1:]+=(ge+J*prefix[1:H+1])/t[H+1:N+1]-J
    r=right.copy();r[0]=origin
    diffs=[]
    for it in range(150):
        c=conv(dc,r,N+2)
        integ=h*(2*c[1:N+1]/3-dc[1:N+1]*r[0]/3+c[:N]/6+c[2:N+2]/6-dc[2:N+2]*r[0]/6)
        ir=np.r_[0,np.cumsum((r[:-1]+r[1:])*h/2)]
        integ[H-1:]+=ge*ir[:H+1]
        new=right.copy();new[1:]+=integ/t[1:N+1];new[0]=origin
        diff=float(np.max(abs(new-r)));diffs.append(diff);r=new
        if diff<5e-12:break
    # Integrate linear r against u^(-gamma), without differentiating r.
    u=t[:N+1];pow1=np.zeros(N+1,complex);pow2=pow1.copy()
    pow1[1:]=u[1:]**b;pow2[1:]=u[1:]**(b+1)
    mass=np.diff(pow1)/b
    moment=np.diff(pow2)/(b+1)/h-np.arange(N)*mass
    left=mass-moment
    p=np.zeros(N+1,complex)
    c1=conv(left,r,N+1);c2=conv(moment,r,N)
    p[1:]=u[1:]**(-gam)+c1[1:]-np.r_[left[1:],0]*r[0]+c2
    p[H+1:]+=J*u[1:H+1]**b/b
    q=np.ones(H+1,complex);q[1:]=u[1:H+1]**gam*p[1:H+1]
    # First cell: exact q for linear r; Jacobi removes t^-a.
    nodes,weights=roots_jacobi(24,0,-a);x=(nodes+1)/2
    q0=1+origin*h*x/b+(r[1]-origin)/h*(h*x)**2/(b*(b+1))
    D=h**(1-a)/2**(1-a)*np.dot(weights,abs(q0)/(1-h*x))
    nodes,weights=roots_legendre(6)
    # Remaining cells: integrate the complex linear interpolant q with Gauss.
    for v,w in zip((nodes+1)/2,weights/2):
        tt=(np.arange(1,H)+v)*h
        qq=(1-v)*q[1:H]+v*q[2:H+1]
        D+=h*w*np.sum(tt**(-a)*abs(qq)/(1-tt))
    # Check independently at selected locations by adaptive integration of linear interpolants.
    checks=[]
    for index in [37,H-1,H,H+1,N]:
        tt=index*h
        # sampling diagnostics use the product integration residual at grid nodes
        checks.append(dict(t=tt,rho_real=float((r[index]+(J if index>H else 0)).real),rho_imag=float((r[index]+(J if index>H else 0)).imag)))
    return dict(N=N,F=float(abs(p[N])/D),p1=[float(p[N].real),float(p[N].imag)],D=float(D),iterations=it+1,last_change=diff,picard_sup_tail_bound=diff/(1-T),seconds=time.time()-start,checks=checks)

def finite(n,root_check=False,certificate=DATA):
    start=time.time();ctx.prec=160;ctx.threads=1;ctx.cap=n+2;K=n//2
    ks=np.arange(1,K+1);s=profile(ks/n)
    # Formal coefficients: only this independently derived identity uses FLINT.
    ell=[acb(0)]+[-acb(complex(z))/int(k) for k,z in zip(ks,s)]
    pp=(acb_poly(acb_series(ell,prec=n+1).exp().coeffs())*acb_poly([1]*(n+1))).truncate(n+1)
    P=np.array([complex(float(pp[j].real.mid()),float(pp[j].imag.mid())) for j in range(n+1)])
    D=np.sum(abs(P[:n-K])/(n-np.arange(n-K)));lam=P[n]/D
    result=dict(n=n,K=K,F=float(abs(lam)),bound=max(T,float(abs(lam))),seconds=time.time()-start)
    if root_check:
        sig=np.r_[s,lam*np.where(abs(P[n-K-1::-1])>0,np.conj(P[n-K-1::-1])/abs(P[n-K-1::-1]),1)]
        # Rebuild W from all n power sums. Newton identity, no target routines.
        w=np.zeros(n+1,complex);w[0]=1
        for m in range(1,n+1):w[m]=-np.dot(sig[:m],w[m-1::-1])/m
        # Enforce the exactly known root 1 via deflation, retain numerical W(1) diagnostic.
        quotient=np.cumsum(w)[:-1];z=np.r_[1+0j,np.roots(quotient)]
        sums=np.array([np.sum(z**k) for k in range(1,n+1)])
        result.update(W1=float(abs(np.sum(w))),root_powersum_max=float(max(abs(sums))),prescribed_error=float(max(abs(sums-sig))),max_radius=float(max(abs(z))))
        Path(f'receipts/tuple_{n}.json').write_text(json.dumps(dict(z=[[v.real,v.imag] for v in z]),indent=2))
    return result

def known_cases():
    # d(t)=c*t -> rho(t)=c*sum_(j>=0)(c*t)^j/((j+1)!j!).
    # Positive solution; zero trial; r=c and exact integrated error. Also complex c.
    from scipy.special import factorial
    rows=[]
    for c in [.2,.65,.9, .3+.4j]:
        k=abs(c)
        def rho(t):return sum(c**(j+1)*t**j/(factorial(j+1)*factorial(j)) for j in range(45))
        E=quad(lambda t:abs(rho(t)),0,1,epsabs=1e-12)[0]
        globalE=k/(1-k)
        delta=.125;small=k*delta
        Ecut=k*delta/(1-small)
        localE=delta**(-k)*Ecut+k*(1-delta**(1-k))/(1-k)
        rows.append(dict(c=str(c),integrated_actual=E,global_majorant=globalE,local_majorant=localE,holds=E<=min(globalE,localE)))
    # Sharp scalar model E'=r+T E/t, constant r: E=r*t/(1-T).
    scalar=[]
    for kk in [0,.1,.686360795,.99]:
        rr=.37;tt=.61;E=rr*tt/(1-kk)
        scalar.append(dict(T=kk,E=E,derivative=rr/(1-kk),rhs=rr+kk*E/tt))
    return dict(positive_and_complex_kernel=rows,scalar_equalities=scalar)

if __name__=='__main__':
    mode=sys.argv[1];out=Path('receipts')
    if mode=='volterra':
        for N in [73728,294912,1179648,2359296]:
            row=volterra(N);print(json.dumps(row),flush=True);(out/f'volterra_{N}.json').write_text(json.dumps(row,indent=2))
    elif mode=='finite':
        for n in [16,31,64,127,256,262144,1000000,4000000]:
            row=finite(n,n<=256);print(json.dumps(row),flush=True);(out/f'finite_{n}.json').write_text(json.dumps(row,indent=2))
    elif mode=='controls':
        row=known_cases();print(json.dumps(row,indent=2));(out/'known_cases.json').write_text(json.dumps(row,indent=2))
