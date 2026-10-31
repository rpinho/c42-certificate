#!/usr/bin/env python3
"""Independent Arb continuous-phase certificate checker; see VOLTERRA_VERIFIER_MATH.md.

No optimizer/search or public-checker imports. Floating collocation produces
only an arbitrary dyadic trial function; a posteriori residual bounds certify it.
"""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
import argparse,hashlib,json,math,sys,time
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.signal import fftconvolve
from flint import acb,acb_poly,arb,arb_poly,ctx,fmpq

if not __debug__:sys.exit('run this program without python -O')

def R(q):
 q=Fraction(q);return arb(fmpq(q.numerator,q.denominator))
def exact_upper(x):return arb(x.upper())
def exact_lower(x):return arb(x.lower())
def midpoint(z):return complex(float(z.real.mid()),float(z.imag.mid()))
def profile(data,N):
 if data.get('constant','C_42')!='C_42':raise ValueError('wrong constant')
 T=Fraction(data['T']);ph=list(map(Fraction,data['phases']));K=len(ph)-1
 if not 0<T<1 or K<1:raise ValueError('profile dimensions or T')
 if isinstance(data.get('breaks'),list):ts=list(map(Fraction,data['breaks']))
 else:
  if K!=24:raise ValueError('explicit rational breaks required for a new profile')
  ts=[Fraction(j,48) for j in range(25)]
 if len(ts)!=len(ph) or ts[0]!=0 or ts[-1]!=Fraction(1,2) or any(t>=u for t,u in zip(ts,ts[1:])):raise ValueError('phase grid')
 if N<16 or N%2 or N>16000000 or any((t*N).denominator!=1 for t in ts):raise ValueError('N must align all profile knots and 1/2')
 if max(abs(x-ph[0]) for x in ph)>=1:raise ValueError('this checker uses the sufficient kappa=T condition max phase excursion<1')
 slopes=[(y-x)/(u-t) for x,y,t,u in zip(ph,ph[1:],ts,ts[1:])]
 return T,ph,ts,slopes

def float_witness(A,B,dnodes,prefix,gamma,g_end,N):
 """Lower-triangular linear collocation; this is not a proof step."""
 af=np.array([midpoint(z) for z in A]);bf=np.array([midpoint(z) for z in B])
 df=np.array([midpoint(z) for z in dnodes]);pf=np.array([midpoint(z) for z in prefix])
 gf=midpoint(gamma);J=2*midpoint(g_end);H=N//2
 # d'(0) is supplied separately by the caller through an attribute.
 v0=float_witness.origin
 t=np.arange(N+1)/N
 rhs=df.copy();rhs[1:]+=af*v0
 rhs[H+1:]+=J*(pf[1:H+1]-t[H+1:])
 kernel=np.zeros(N+1,complex);kernel[1:N]=bf[1:]+af[:-1]
 denominator=t-bf[0]
 v=np.zeros(N+1,complex);v[0]=v0
 def solve(left,right):
  if right-left<=48:
   for i in range(left,right):
    if i>left:rhs[i]+=np.dot(v[left:i],kernel[i-left:0:-1])
    v[i]=rhs[i]/denominator[i]
   return
  middle=(left+right)//2;solve(left,middle)
  product=fftconvolve(v[left:middle],kernel[:right-left])
  rhs[middle:right]+=product[middle-left:right-left]
  solve(middle,right)
 solve(1,N+1)
 if not np.isfinite(v).all():raise ArithmeticError('nonfinite trial')
 return v,J

def verify(data,N=30720,bits=192,threads=1,verbose=False,weighted_residual=False):
 start=time.time();ctx.prec=bits;ctx.threads=threads
 Tq,ph,ts,slopes=profile(data,N);T=R(Tq);H=N//2;h=R(Fraction(1,N))
 Cq=Fraction(data.get('C',data['T']))
 if not Tq<=Cq<1:raise ValueError('claim must satisfy T<=C<1')
 def progress(stage):
  if verbose:print(stage,round(time.time()-start,3),flush=True,file=sys.stderr)
 gamma=T*acb(0,R(ph[0])).exp();a=gamma.real;br=1-a;b=1-gamma;den=b*(b+1)
 if not a>0 or not br>0:raise ArithmeticError('need certified 0<Re gamma<1')
 beta_cap=R(max(abs(s) for s in slopes));g_end=T*acb(0,R(ph[-1])).exp()
 dprime0=-acb(0,1)*R(slopes[0])*gamma
 A=[acb(0)]*N;B=[acb(0)]*N;dn=[acb(0)]*(N+1)
 dc=[acb(0)]*N;d1=[acb(0)]*N;d2=[acb(0)]*N
 for seg,betaq in enumerate(slopes):
  left=int(ts[seg]*N);right=int(ts[seg+1]*N);beta=R(betaq)
  z=acb(0,beta*h);step=z.exp();half=(z/2).exp();power=acb(1);e1=acb(0);e2=acb(0)
  for k in range(16):
   if k:power*=z/k
   e1+=power/(k+1);e2+=power/(k+2)
  err=T*h*abs(z)**16*abs(z).exp()/math.factorial(16)
  er=arb(0,err.upper());unc=acb(er,er)
  g=T*acb(0,R(ph[seg])).exp()
  zero=(betaq==0 and ph[seg]==ph[0])
  for i in range(left,right):
   if zero:
    A[i]=B[i]=dc[i]=d1[i]=d2[i]=dn[i]=acb(0)
   else:
    A[i]=gamma*h/2-g*h*e2+unc;B[i]=gamma*h/2-g*h*(e1-e2)+unc
    dn[i]=gamma-g;gc=g*half;dc[i]=gamma-gc;d1[i]=-acb(0,1)*beta*gc;d2[i]=beta**2*gc
   g*=step
 dn[0]=acb(0);dn[H]=gamma-g_end
 for i in range(H,N):
  A[i]=B[i]=gamma*h/2;dc[i]=gamma
  if i>H:dn[i]=gamma
 dn[N]=gamma
 prefix=[acb(0)]
 for x,y in zip(A,B):prefix.append(prefix[-1]+x+y)
 progress('moments')
 float_witness.origin=midpoint(dprime0)
 vf,Jf=float_witness(A,B,dn,prefix,gamma,g_end,N)
 trial_hash=hashlib.sha256(np.asarray(vf,dtype='<c16').tobytes()+np.asarray([Jf],dtype='<c16').tobytes()).hexdigest()
 V=[acb(float(z.real),float(z.imag)).mid() for z in vf];J=acb(float(Jf.real),float(Jf.imag)).mid()
 sigma=[(V[i+1]-V[i])*N for i in range(N)]
 delta=[acb(0)]+[sigma[i]-sigma[i-1] for i in range(1,N)]
 progress('dyadic_trial')
 vp=acb_poly(V)
 av=(acb_poly(A)*vp).truncate(N+1);bv=(acb_poly(B)*vp).truncate(N+1)
 cm=(acb_poly(delta)*acb_poly(dc)).truncate(N)
 t=[R(Fraction(i,N)) for i in range(N+1)]
 rn=[arb(0)]*(N+1)
 for i in range(1,N+1):
  conv=bv[i]-(B[i]*V[0] if i<N else 0)+av[i-1]
  jump=J*prefix[i-H] if i>H else acb(0)
  value=dn[i]+conv+jump-t[i]*(V[i]+(J if i>H else 0))
  if i==1:r_first=value
  rn[i]=exact_upper(abs(value))
  if i==H:rh_plus=exact_upper(abs(value+g_end-t[i]*J))
 total_delta=sum((abs(z) for z in delta[1:]),arb(0))
 B3=T*beta_cap**3+T*beta_cap**2*abs(V[0])+T*beta_cap*(abs(sigma[0])+total_delta)+abs(J)*T*beta_cap**2
 residual=arb(0);worst=0
 A1=arb(0);Ah=arb(0);Rh=arb(0);RN=arb(0)
 left_A=arb(0);left_N=arb(1)
 cuts=[]
 if weighted_residual:
  for dq in (Fraction(1,256),Fraction(1,64),Fraction(1,16)):
   if (dq*N).denominator!=1:continue
   phase_values=[ph[j] for j,z in enumerate(ts) if z<=dq]
   for j in range(len(slopes)):
    if ts[j]<=dq<=ts[j+1]:phase_values.append(ph[j]+slopes[j]*(dq-ts[j]));break
   kq=Tq*max(abs(x-ph[0]) for x in phase_values)
   assert 0<=kq<Tq
   cuts.append(dict(delta=dq,index=int(dq*N),kappa=kq,local=arb(0),previous=arb(0),
                    earlyA=arb(0),earlyRN=arb(0),earlyR=arb(0)))
 for i in range(N):
  rcc=d2[i]+d1[i]*V[0]+sigma[0]*dc[i]+cm[i]-2*sigma[i]
  if i>=H:rcc+=J*d1[i-H]
  second=abs(rcc)+h*B3/2
  if i==0:
   bound=abs(dprime0-V[0])+h*second/2
   beta0=R(slopes[0])
   first_third=T*abs(acb(0,1)*beta0**3+beta0**2*V[0]-acb(0,1)*beta0*sigma[0])
   endpoint=exact_upper(abs(dprime0-V[0]));right=exact_upper(abs(r_first)/h)
   endpoint=endpoint if endpoint>right else right
   hermite=endpoint+h*h*first_third/24
   if hermite.upper()<bound.upper():bound=hermite
  else:
   edge=rh_plus if i==H else rn[i]
   edge=exact_upper(edge if edge>rn[i+1] else rn[i+1])
   bound=(edge+h*h*second/8)/t[i]
  bound=exact_upper(bound)
  if bound>residual:residual=bound;worst=i
  if weighted_residual:
   right_A=t[i+1]**(1-T)
   right_N=(1-t[i+1])**br if i+1<N else arb(0)
   moment_A=(right_A-left_A)/(1-T)
   moment_N=(left_N-right_N)/br
   A1+=bound*moment_A
   RN+=bound*moment_N
   if i<H:
    Ah+=bound*moment_A;Rh+=bound*h
   left_A=right_A;left_N=right_N
   for cut in cuts:
    if i>=cut['index']:continue
    kk=R(cut['kappa']);new=t[i+1]**(1-kk)
    cut['local']+=bound*(new-cut['previous'])/(1-kk);cut['previous']=new
    cut['earlyA']+=bound*moment_A;cut['earlyRN']+=bound*moment_N;cut['earlyR']+=bound*h
 eps=exact_upper(residual/(1-T))
 progress('residual_certified')
 # Closed antiderivatives; p(1), and q=t^gamma*p on the denominator interval.
 Z=[acb(0)]+[acb(t[i])**b for i in range(1,N+1)]
 G2=[acb(0)]+[t[i]*Z[i]/den for i in range(1,N+1)]
 ramp=(acb_poly(delta)*acb_poly(G2)).truncate(N+1)
 p1=1+V[0]/b+sigma[0]/den+ramp[N]+J*Z[H]/b
 q=[acb(1)]+[1+V[0]*t[i]/b+sigma[0]*t[i]**2/den+(t[i]/Z[i])*ramp[i] for i in range(1,H+1)]
 progress('trial_integrated')
 P=[arb(0)]+[abs(Z[i]) for i in range(1,H+1)] # exactly |t^b|=t^(1-Re gamma), enclosed.
 ipoly=arb_poly([1/(br+k) for k in range(96)])
 I=[arb(0)]
 for i in range(1,H+1):
  tail=P[i]*t[i]**96/((br+96)*(1-t[i]))
  I.append(P[i]*ipoly(t[i])+arb(0,tail.upper()))
 lower=arb(0)
 for i in range(H):
  W0=I[i+1]-I[i];W1=W0-(P[i+1]-P[i])/br
  hat_right=(W1-t[i]*W0)/h;hat_left=W0-hat_right
  lower+=exact_lower(abs(hat_left*q[i]+hat_right*q[i+1]))
 # Integral bounds on q'' via three positive real convolutions.
 ad=arb_poly([arb(0)]+[exact_upper(abs(delta[i])) for i in range(1,H)])
 K0=[];K1=[];K2=[]
 for i in range(H):
  K0.append((P[i+1]-P[i])/br)
  K1.append((t[i+1]*P[i+1]-t[i]*P[i])/(br+1))
  K2.append((t[i+1]**2*P[i+1]-t[i]**2*P[i])/(br+2))
 c0=(ad*arb_poly(K0)).truncate(H);c1=(ad*arb_poly(K1)).truncate(H);c2=(ad*arb_poly(K2)).truncate(H)
 coeff0=abs(sigma[0]/den);coeff1=abs(gamma/(b+1));coeff2=abs(2*gamma/b)
 interp=coeff0*h*h*P[1]/((1-h)*(br+1)*(br+2))
 for i in range(1,H):
  second_integral=2*coeff0*h+coeff1*c2[i]/(t[i]*P[i])+coeff2*c1[i]/P[i]+t[i+1]*c0[i]/P[i+1]
  weight=P[i]/(t[i]*(1-t[i+1]))
  interp+=h*h*weight*second_integral/8
 lower_trial=exact_lower(lower-interp)
 correction=exact_upper(eps/br*(I[H]-P[H]/br))
 numerator_error=exact_upper(eps/br)
 integrated={}
 if weighted_residual:
  beta1=T.gamma()*br.gamma()/(T+br).gamma()
  beta2=T.gamma()*(br+1).gamma()/(T+br+1).gamma()
  EN=exact_upper(RN+T*A1*beta1)
  ED=exact_upper(2/br*(t[H]**br*Rh+T*Ah*t[H]**(T+br)*beta2))
  integrated=dict(A1=A1.str(20),Ah=Ah.str(20),Rh=Rh.str(20),RN=RN.str(20),
                  numerator_error=EN.str(20),denominator_error=ED.str(20))
  local_rows=[]
  for cut in cuts:
   dd=R(cut['delta']);kk=R(cut['kappa']);Edelta=dd**kk*cut['local']
   B1=dd**(-T)*Edelta+A1-cut['earlyA'];Bh=dd**(-T)*Edelta+Ah-cut['earlyA']
   local_EN=exact_upper((1-dd)**(-a)*Edelta+RN-cut['earlyRN']+T*B1*beta1)
   local_ED=exact_upper(2/br*(t[H]**br*(Edelta+Rh-cut['earlyR'])+T*Bh*t[H]**(T+br)*beta2))
   local_rows.append(dict(delta=str(cut['delta']),kappa=str(cut['kappa']),
                          E_delta=Edelta.str(20),numerator_error=local_EN.str(20),denominator_error=local_ED.str(20)))
   if local_EN<EN:EN=local_EN
   if local_ED<ED:ED=local_ED
  integrated['local_contractions']=local_rows
  integrated['best_numerator_error']=EN.str(20);integrated['best_denominator_error']=ED.str(20)
  if EN<numerator_error:numerator_error=EN
  if ED<correction:correction=ED
 D=exact_lower(lower_trial-correction)
 if not D>0:raise ArithmeticError('denominator lower bound not positive')
 numerator=exact_upper(abs(p1)+numerator_error);upper=exact_upper(numerator/D)
 Fq=Fraction(str(upper.fmpq()));verdict='PASS' if Fq<Cq else 'UNPROVED'
 progress('denominator_certified')
 return dict(verdict=verdict,C=str(Cq),T=str(Tq),N=N,bits=bits,threads=threads,a=a.str(20),kappa=str(Tq),
             residual_upper=residual.str(16),rho_error_upper=eps.str(16),worst_cell=worst,B3=B3.str(12),
             p_tilde_1=p1.str(40),numerator_upper=numerator.str(30),denominator_lower=D.str(30),
             denominator_trial_lower=lower_trial.str(30),denominator_interpolation_error=interp.str(14),
             weighted_residual=weighted_residual,integrated_residual=integrated,
             F_upper=upper.str(30),F_upper_rational=str(Fq),exact_margin=str(Cq-Fq),
             trial_sha256=trial_hash,seconds=time.time()-start)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('certificate');p.add_argument('--N',type=int,default=30720);p.add_argument('--bits',type=int,default=192)
 p.add_argument('--threads',type=int,default=1);p.add_argument('--verbose',action='store_true')
 p.add_argument('--weighted-residual',action='store_true');args=p.parse_args()
 raw=Path(args.certificate).read_bytes();codehash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 try:
  if args.bits<128 or not 1<=args.threads<=8:raise ValueError('precision or threads')
  result=verify(json.loads(raw),args.N,args.bits,args.threads,args.verbose,args.weighted_residual)
 except Exception as e:result=dict(verdict='REFUSED',reason=str(e))
 result.update(certificate_sha256=hashlib.sha256(raw).hexdigest(),verifier_sha256=codehash)
 print(json.dumps(result,indent=2));sys.exit(0 if result['verdict']=='PASS' else 1)
