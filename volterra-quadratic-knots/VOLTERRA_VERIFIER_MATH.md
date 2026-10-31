# Independent continuous-phase verifier: derivation before new search

Written from Bakeas, turan-c42-bound/NOTES.md (Theorem A, Lemma 3), after
reproducing its new public T=.6865 certificate with both public checkers.
No code from either checker is imported. The new phase optimizer has not
yet been written. This is an unrefereed extension of the computation, not
a new claim about the validity of Bakeas's analytic theorem.

Let gamma=T exp(i phi(0)), a=Re gamma, b=1-gamma, br=1-a. Assume 0<a<1,
T<1, and max|phi-phi(0)|<1. Then kappa=T by |1-exp(i theta)|<=|theta|.
The profile is continuous piecewise linear in phase on [0,1/2], constant
modulus T. Theorem A and its proof work for any finite rational knot grid;
all knots and 1/2 must align with the verifier's uniform grid.

## Trial rho and rigorous residual

Use any exact dyadic continuous piecewise-linear function v plus a dyadic
jump J at 1/2: rho_tilde=v+J 1_(t>1/2). A floating collocation constructs
only a witness; every subsequent bound is outward-rounded Arb. Let
sigma_i=(v_(i+1)-v_i)/h and Delta_i=sigma_i-sigma_(i-1).

For each cell, A_j=int d(u)*(u-t_j)/h du and
B_j=int d(u)*(t_(j+1)-u)/h du. They are explicit exponential moments;
a finite exponential series has an explicit exponential tail bound.
Convolutions give nodal residuals
R_i=d(t_i)+(B*v)_i-B_i v_0+(A*v)_(i-1)
    +J int_0^(t_i-1/2) d -t_i(v_i+J 1_(i>H)).
At H both one-sided values are used. Inside cell i,
R''=d''+d'v_0+sigma_0 d+sum Delta_k d(t-t_k)
      +J d'(t-1/2)-2sigma_i.
The derivative bound B3 from NOTES.md bounds R''' there. Linear
interpolation gives sup|R/t| on all cells; the first cell uses R(0)=0 and
R'(0)=d'(0)-v_0. Lemma 3 gives eps=sup|R/t|/(1-T).

## Exact integration of the trial

Write G1(t)=t^b/b and G2(t)=t^(b+1)/(b(b+1)). Then
p_tilde(t)=t^-gamma+v_0 G1(t)+sigma_0 G2(t)
    +sum Delta_k G2(t-t_k)+J G1(t-1/2).
All nodal values follow from one ordinary polynomial convolution.
The numerator is bounded by |p_tilde(1)|+eps/br.

## Denominator with its singular factor removed

Put q(t)=t^gamma p_tilde(t). On [0,1/2],
q=1+v_0 t/b+sigma_0 t^2/(b(b+1))
       +sum Delta_k t^gamma (t-t_k)^(2-gamma)/(b(b+1)).
The denominator is int t^-a |q(t)|/(1-t) dt. On each cell let L be the
linear interpolant of q at its endpoints. A rigorous lower bound is
|int L(t) W(t) dt|-int|q-L|W, W=t^-a/(1-t).
The two weighted hat moments use I_br(t) and
I_(br+1)(t)=I_br(t)-t^br/br. I_br is evaluated by a positive geometric
series with explicit tail t^br*t^K/((br+K)(1-t)).

For t>t_k, the nonpolynomial term has second derivative
Delta_k[-gamma/(b+1)t^(gamma-2)(t-t_k)^(2-gamma)
 +2gamma/b*t^(gamma-1)(t-t_k)^(1-gamma)
 +t^gamma*(t-t_k)^(-gamma)].
On cells i>=1, integrate its modulus by bounding the t factors at their
appropriate endpoints. The remaining integrals are convolutions of
|Delta_k| with exact power moments for exponents 2-a,1-a,-a. All are
integrable since a<1. Add 2|sigma_0/(b(b+1))|h. The Peano-kernel bound
int_cell|q-L| <= h^2/8 int_cell|q''| then applies. Multiply by sup_cell W.
On the first cell q is exactly quadratic; its weighted error is at most
|sigma_0/(b(b+1))| h^(3-a)/((1-h)(2-a)(3-a)).

Finally subtract eps/br * I_(br+1)(1/2) to pass from p_tilde to p.
If the resulting denominator lower bound is positive, divide the rigorous
numerator upper bound by it. PASS requires F_upper<T, compared exactly as
rationals. The certified C42 bound is T, never the smaller F value.

Factoring t^-gamma avoids the near-zero denominator interpolation loss in
the public checkers. All interpolation and residual inequalities are the
complex-valued Peano-kernel inequalities, not a real-valued Rolle argument.

## Sharpened first-cell residual (added before adaptive-profile search)

For f=R on [0,h], f(0)=0. Hermite interpolation through f'(0) and f(h)
has remainder bounded by sup|f'''|*t^2(h-t)/6. This is a complex-valued
Peano-kernel bound: the kernel is everywhere nonpositive, since
(t-s)^2 <= (t/h)^2(h-s)^2 for 0<=s<=t. Dividing by t yields

 sup|R/t| <= max(|R'(0)|,|R(h)|/h) + h^2 sup|R'''|/24.

On the first cell, before any shifted ramp or jump can enter,
R'''(t)=g(t)*(i beta_0^3+beta_0^2 v_0-i beta_0 sigma_0).
Since |g(t)|=T, this gives an exact local third-derivative bound. The
verifier takes the smaller certified upper endpoint of this bound and
the original first-cell bound. This sharpening is independent of search.

## Integrated residual bound (derived before implementation, 20:36 UTC)

Let r(t) be any nonnegative piecewise-constant certified upper bound for
|R(t)/t|, and put kappa=T<1. For e=|rho-rho_tilde|, the residual identity
gives e(t)<=r(t)+T/t int_0^t e(v)dv. Since e is bounded, an integrating
factor for E(t)=int_0^t e gives

 E(t)<=t^T int_0^t r(v)v^(-T)dv,
 e(t)<=r(t)+T*t^(T-1) int_0^t r(v)v^(-T)dv.

The zero-endpoint boundary term vanishes because E(t)=O(t) and T<1.
Write A_x=int_0^x r(v)v^(-T)dv and R_x=int_0^x r(v)dv. Then the error
in p(1) is at most

 EN=int_0^1 r(v)(1-v)^(-a)dv + T*A_1*Beta(T,1-a).

For the denominator, bound 1/(1-t)<=2 on [0,1/2] and use Tonelli on
the nonnegative majorant. With br=1-a,

 ED <= (2/br)*(1/2)^br*R_(1/2)
       +(2*T/br)*A_(1/2)*(1/2)^(T+br)*Beta(T,br+1).

All r-cell moments are explicit positive power differences; their endpoint
singularities are integrable because T<1 and a<1. Evaluate EN and ED
outward in Arb, and take the smaller certified upper endpoint of EN and
eps/br (respectively ED and the original denominator correction). No
floating residual samples enter this argument. This is an optional
sharpening; the original sup-norm proof remains available. Its intended
purpose is to account for residuals concentrated near small phase knots.

The certificate may optionally state a distinct exact rational claim C,
with T<=C<1. Theorem A concludes C42<=max(T,F), so F_upper<C suffices.
T remains the modulus of the prescribed profile; it is never replaced by
the smaller F value. If C is absent the checker uses C=T as before.

## Local contraction before a rational cutoff (derived before code)

For an aligned rational 0<delta<=1/2, let k_delta be an exact upper bound
for |d(u)| on [0,delta], with 0<=k_delta<=T. For our linear phase this can
be taken as T times the largest phase excursion up to delta, including
the interpolated phase at delta. The same scalar inequality on this first
interval gives

 E_delta=delta^k_delta int_0^delta r(v)v^(-k_delta)dv >= int_0^delta e.

For t>=delta, integrating the global inequality from delta gives

 int_0^t e <= t^T [delta^(-T) E_delta + int_delta^t r(v)v^(-T)dv].

Define B_x=delta^(-T) E_delta+int_delta^x r(v)v^(-T)dv. Split the
numerator error integral at delta and use (1-v)^(-a)<=(1-delta)^(-a)
on the early part. Valid error bounds are

 EN_delta=(1-delta)^(-a)E_delta
          +int_delta^1 r(v)(1-v)^(-a)dv+T B_1 Beta(T,br),
 ED_delta=(2/br)[(1/2)^br(E_delta+int_delta^(1/2) r(v)dv)
          +T B_(1/2)(1/2)^(T+br)Beta(T,br+1)].

These remain valid when k_delta=0; no Gamma(k_delta) is used. The program
may take the minimum certified upper endpoint over several rational delta
values and the original bounds. This is useful because d(0)=0 and the
global estimate |d|<=T can be very wasteful for errors near the origin.
