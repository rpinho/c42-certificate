# Adapted public cross-checker

`check_bakeas_adapted.py` is Bakeas's MIT-licensed `independent_check.py`,
with these explicit changes (full diff in bakeas_checker_adaptation.diff):

1. Read explicit rational phase knots instead of hardcoding j/48; require
   every knot to align with the uniform integration grid.
2. Replace its quadratic-time FLOATING trial solve by a divide-and-conquer
   convolution solve. This only constructs a trial function. Every original
   a posteriori residual, interpolation, and integral error bound is retained.
   No proof bound or optimizer is imported from our verifier/search code.
3. Use two FLINT threads and emit a machine-readable receipt.

The exact residual and denominator algorithms are unchanged. This checker
therefore provides a separate arithmetic implementation of the proof bounds,
using the public author's original denominator treatment. It is not described
as an unmodified public checker. Both original public checkers were separately
run unchanged on the published .6865 profile.

`check_bakeas_integrated.py` is a separately retained further variant,
whose complete original-to-variant diff is bakeas_integrated_adaptation.diff.
It keeps the public author's residual and denominator enclosures, and adds
our integrated scalar Volterra majorant from VOLTERRA_VERIFIER_MATH.md.
This variant therefore DOES change the propagation of the residual bound.
It also permits any finite rational knot count and a separate exact claim
C>=T. It imports neither our verifier nor search code. Both source and
certificate hashes are captured at the beginning of each run.

`check_bakeas_local.py` additionally implements the local-contraction
cutoffs from that same derivation, again using the public residual and
denominator enclosures. Full diff: bakeas_local_adaptation.diff. Both new
variants are our mathematical extensions of public code; neither is
described as an unmodified public checker or a new independent author.

The final `check_bakeas_local` also follows version 1.1.0 of the original on
three points: it refuses to run under `python -O`, it checks that the cell
bounds, numerator and denominator are finite, and it exits non-zero on
failure. `verifier_versions/GUARD_CHANGE.diff` shows this against the earlier
version.
