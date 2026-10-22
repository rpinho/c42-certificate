# An eight-step profile for Turán's pure power-sum problem

This note gives a proposed computer-assisted upper bound

$$
C_{42}=\limsup_{n\to\infty}\min_{\max_j|z_j|=1}
             \max_{1\le k\le n}\left|\sum_{j=1}^n z_j^k\right|
\le 0.688983.
$$

Author: André Röhrig. September 9, 2026.

Exact rational parameters and outward-rounded interval calculations certify
the strict inequalities used below. The analytic argument is included in
full. This is an asymptotic construction, and no explicit finite threshold
$N$ is claimed. The author reviewed the proof and supporting material
before publication; the result has not been peer reviewed.

## Attribution and change from the previous construction

The coefficient-cancellation construction is adapted from Sebastian Griego's
[two-block certificate](https://github.com/sebastian-griego/turan-c42-certificate/tree/v1.0.0),
whose stated bound is $0.6906538$. The earlier literature includes A. Biró,
*An upper estimate in Turán's pure power sum problem*, Indagationes
Mathematicae 11 (2000), 499–508. We use Griego's construction as the starting
point and replace its single constant middle profile by eight separate
complex values. The additional degrees of freedom give the improvement.
No correctness claim about an unpublished source is needed as an extra
hypothesis: the argument for the generalized construction follows below.

The difference from the [listed bound](https://teorth.github.io/optimizationproblems/constants/42a.html)
is exactly $0.0016708=4177/2500000$.

## Exact parameters and the limiting inequality

All displayed terminating decimals are exact rationals. Set

$$
 C=0.688983,\qquad \tau=1/3,\qquad B=2/3,
 \qquad \alpha=0.6460641528+0.5911224477i,
 \qquad s=1-\alpha.
$$

Let $t_j=(8+j)/24$, for $0\le j\le8$. On each middle interval
$t_j<u<t_{j+1}$, use the following value $\eta_j$:

| $j$ | $\operatorname{Re}\eta_j$ | $\operatorname{Im}\eta_j$ |
| --- | --- | --- |
| 0 | 0.4826503951 | -0.4916756331 |
| 1 | 0.5360145885 | -0.4328795363 |
| 2 | 0.5637805037 | -0.3960402450 |
| 3 | 0.5881806899 | -0.3588033000 |
| 4 | 0.6126301773 | -0.3152468841 |
| 5 | 0.6317490513 | -0.2749353892 |
| 6 | 0.6531184739 | -0.2193914107 |
| 7 | 0.6719343778 | -0.1523171824 |

Define the jumps $\delta_0=\eta_0-s$ and
$\delta_j=\eta_j-\eta_{j-1}$ for $j\ge1$. For $a=\operatorname{Re}\alpha>0$, put

$$
 I_\alpha(x)=\int_0^x\frac{v^{\alpha-1}}{1-v}\,dv,
 \qquad K=I_\alpha(\tau),\qquad D=I_a(\tau)>0.
$$

Positive real bases use their real logarithms in complex powers. Define

$$
 J_\alpha(t,r)=
 \int_{\substack{u\ge t,\ v\ge r\\u+v\le1}}
 \frac{(1-u-v)^{\alpha-1}}{uv}\,du\,dv.
$$

Set $J_\alpha(t,r)=0$ when $t+r\ge1$. If $t+r<1$, the substitution
$z=1-u-v$, followed by partial fractions in $u$, gives

$$
 J_\alpha(t,r)=\int_0^{1-t-r}\frac{z^{\alpha-1}}{1-z}
 \log\frac{(1-t-z)(1-r-z)}{tr}\,dz.
$$

In particular $J_\alpha(\tau,\tau)$ is Griego's double-integral term.
The endpoint $B=1-\tau$ never contributes to these double integrals,
since every other correction index is at least $\tau$.

The complex number to be certified is

$$
Y=1+sK-\sum_{j=0}^7\delta_j\bigl(I_\alpha(1-t_j)-K\bigr)
  +\frac12\sum_{i,j=0}^7\delta_i\delta_jJ_\alpha(t_i,t_j).
$$

The checker establishes, by exact rational comparisons of interval endpoints,

$$
 |s|<C,\quad |\eta_j|<C\ (0\le j\le7),\quad |Y|<CD.
$$

For readability, short outward enclosures from the calculation are

$$
\begin{aligned}
 D&\in[0.888190650452511954178223193231,
        0.888190650452511954178223193232],\\
 |Y|/D&\in[0.688982098504060748505120997583,
             0.688982098504060748505120997584],\\
 C^2D^2-|Y|^2&\in[0.000000979973624636182168032566,
                      0.000000979973624636182168032567].
\end{aligned}
$$

Here $0.688982098504\ldots$ is an evaluated term in the limiting
criterion. It is **not** the claimed upper bound: the profile magnitudes
must be bounded as well. The claimed rational bound is $0.688983$.

## A finite construction for all sufficiently large $n$

Put $A=\lfloor n/3\rfloor$, and partition $1,\dots,n$ into

$$
 E=\{1,\dots,A\},\quad M=\{A+1,\dots,n-A-1\},
 \quad F=\{n-A,\dots,n\}.
$$

For $m\in E$, set $S_m=s$. For $m\in M$, set

$$
 w_m=\sum_{j=0}^7\delta_j\mathbf1_{m>\lfloor nt_j\rfloor},
 \qquad S_m=s+w_m.
$$

The sum telescopes, so each middle value is one of the $\eta_j$.
Let $\beta_l=(\alpha)_l/l!$, the coefficients of $(1-z)^{-\alpha}$, and

$$
 P_m=\frac{n}{m}\beta_{n-m}\quad(m\in F).
$$

All $P_m$ are nonzero since $\operatorname{Re}\alpha>0$. Define

$$
 Y_n=n\beta_n
 -n\sum_{m\in M}\frac{w_m\beta_{n-m}}m
 +\frac n2\sum_{\substack{u,v\in M\\u+v\le n}}
             \frac{w_uw_v\beta_{n-u-v}}{uv}
 +s\sum_{m\in F}P_m.
$$

We will prove $|Y_n|/\sum_{m\in F}|P_m|\to|Y|/D<C$.
Consequently, for all sufficiently large $n$, the values

$$
 S_m=\frac{Y_n}{\sum_{r\in F}|P_r|}
             \frac{\overline{P_m}}{|P_m|}\qquad(m\in F)
$$

satisfy $|S_m|\le C$ and $\sum_{m\in F}S_mP_m=Y_n$.

Consider coefficients $b_l$ defined through degree $n$ by

$$
 \mathcal B(z)=(1-z)^{-1}
       \exp\left(-\sum_{m=1}^n\frac{S_mz^m}{m}\right)
       =\sum_{l\ge0}b_lz^l.
$$

Through degree $n$, factor this as $(1-z)^{-\alpha}$ times the
exponential of the corrections $S_m-s$, all at indices $m>A$.
Since $3(A+1)>n$ **for every integer $n$**, at most two corrections
can enter $b_n$. A final correction has index at least $n-A$, so
its sum with any other correction exceeds $n$. It follows exactly that

$$
 n b_n=Y_n-\sum_{m\in F}S_mP_m=0.
$$

This also explains why $\tau=1/3$ is allowed, despite the strict
inequality $\tau>1/3$ sufficient for the same correction expansion.

Let $y_2,\dots,y_n$ be the roots, with multiplicity, of

$$
 p_n(Z)=Z^{n-1}+b_1Z^{n-2}+\cdots+b_{n-1},
 \qquad y_1=1.
$$

Logarithmic differentiation of $\mathcal B$, or Newton's identities,
shows that for $1\le k\le n$, its coefficients and the power sums
$Q_k=\sum_{j=1}^ny_j^k$ satisfy the same triangular recurrence:

$$
 Q_k+\sum_{j=1}^{k-1}b_jQ_{k-j}
   =1+\sum_{j=1}^{k-1}b_j-kb_k
   =S_k+\sum_{j=1}^{k-1}b_jS_{k-j}.
$$

For $k=n$, use $b_n=0$, because $p_n$ has degree $n-1$.
Induction gives $Q_k=S_k$. Finally set
$\Lambda=\max_j|y_j|\ge1$ and $z_j=y_j/\Lambda$. Then

$$
 \max_j|z_j|=1,\qquad
 \left|\sum_jz_j^k\right|=\Lambda^{-k}|S_k|\le C
 \quad(1\le k\le n).
$$

## Justification of the limit, including the singular endpoint

The gamma-ratio asymptotic gives

$$
 \beta_l=\frac{l^{\alpha-1}}{\Gamma(\alpha)}(1+O(l^{-1})).
$$

For the fixed parameter $a=\operatorname{Re}\alpha>0$, it also gives
$|\beta_l|\le K_0(l+1)^{a-1}$ for all $l\ge0$, with a finite constant
$K_0$. Restricting first to $l\ge\varepsilon n$, the asymptotic is
uniform and the relevant sums converge as Riemann sums. There are finitely
many step discontinuities, which cause no difficulty on that restricted
domain. Integer rounding of the breakpoints does not change these limits.

To justify subsequently sending $\varepsilon$ to zero, observe that

$$
 \sum_{0\le l<\varepsilon n}|\beta_l|
       =O((\varepsilon n+1)^a).
$$

In every single sum, the remaining denominator is $m\asymp n$. In
the double sum, group by $l=n-u-v$; for each fixed $l$, there are
at most $n$ terms and $1/(uv)=O(n^{-2})$. Thus, after normalization
by $n^{1-a}$, the absolute contribution of either type of endpoint
region is $O(\varepsilon^a+n^{-a})$. The middle weights are uniformly
bounded, so they preserve this estimate. The same argument applies to
the sum of absolute values over $F$. This proves convergence including
$l=0$, rather than treating a singular endpoint as an ordinary Riemann sum.

With $w(u)=\sum_j\delta_j\mathbf1_{u>t_j}$ on $(\tau,B)$, the limits are

$$
\begin{aligned}
\Gamma(\alpha)n^{1-\alpha}\beta_n&\longrightarrow1,\\
\Gamma(\alpha)n^{1-\alpha}
 \sum_{m\in M}\frac{w_m\beta_{n-m}}m
 &\longrightarrow\sum_j\delta_j(I_\alpha(1-t_j)-K),\\
\Gamma(\alpha)n^{1-\alpha}
 \sum_{\substack{u,v\in M\\u+v\le n}}
 \frac{w_uw_v\beta_{n-u-v}}{uv}
 &\longrightarrow\sum_{i,j}\delta_i\delta_jJ_\alpha(t_i,t_j),\\
\Gamma(\alpha)n^{1-\alpha}
 \sum_{m\in F}\frac{\beta_{n-m}}m&\longrightarrow K,\\
|\Gamma(\alpha)|n^{1-a}
 \sum_{m\in F}\frac{|\beta_{n-m}|}m&\longrightarrow D.
\end{aligned}
$$

Therefore $\Gamma(\alpha)n^{-\alpha}Y_n\to Y$, and
$|\Gamma(\alpha)|n^{-a}\sum_{m\in F}|P_m|\to D>0$. Taking absolute
values and their quotient proves the required limit. The strict certified
inequality implies $R_n\le C$ eventually and hence $C_{42}\le C$.

## The rigorous calculation

The script `scripts/certify_turan.py` computes the elementary integrals by
convergent series with explicit tails. For $0<x<1$,

$$
 I_\alpha(x)=\sum_{k=0}^{N-1}\frac{x^{\alpha+k}}{\alpha+k}+E,
 \qquad |E|\le\frac{x^{a+N}}{(a+N)(1-x)}.
$$

For $J_\alpha(t,r)$, set $L=1-t-r>0$, $b=1-t$, $d=1-r$,

$$
 c_0=\log\frac{bd}{tr},\qquad
 c_k=c_0-\sum_{q=1}^k\frac{b^{-q}+d^{-q}}q.
$$

Expanding the one-dimensional integral gives

$$
 J_\alpha(t,r)=\sum_{k=0}^{N-1}c_k\frac{L^{\alpha+k}}{\alpha+k}+E_J,
$$

with the valid absolute tail bound

$$
 |E_J|\le\frac{L^a}{a+N}
 \left(\frac{c_0L^N}{1-L}
 +\frac{(L/b)^N}{t(1-L/b)}
 +\frac{(L/d)^N}{r(1-L/d)}\right).
$$

Indeed $c_0>0$ and
$|c_k|\le c_0+b^{-k}/t+d^{-k}/r$; summing three geometric series proves
the bound. Every real and imaginary component of each tail is enclosed
by the corresponding symmetric error interval.

Reproduce the proof's inequalities with:

```sh
python3 scripts/certify_turan.py
python3 scripts/certify_turan.py --dps 200 --terms 480
```

The exact certificate is `certificates/turan.json`; the default verification
transcript is `certificates/turan_verification.json`. Both precisions pass.

This result, its verification code, and the exposition were developed with
AI assistance using Codex. André Röhrig reviewed the mathematical arguments,
references, and supporting verification material before publication.
The verification commands were executed and checked during preparation.
