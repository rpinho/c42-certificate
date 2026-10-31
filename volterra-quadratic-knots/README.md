# C42 <= 0.6863607950 on Bakeas's Theorem A

Unrefereed. This certifies an upper bound on Turán's pure power-sum constant, 42a in [Tao's optimization-constants table](https://teorth.github.io/optimizationproblems/constants/42a.html).

This is Vasileios Bakeas's method, not a new construction. Theorem A, the residual lemma and the phase construction are theirs ([turan-c42-bound v1.1.0](https://doi.org/10.5281/zenodo.23836302)), building on Harcos's and Biró's finite linearization. Their notes already report about 0.686361 as the floating optimum of this family, so this mostly certifies that value.

New here are the non-uniform grid, 49 rational phases on knots j²/4608, and two error-propagation bounds, an integrated residual estimate and a smaller kernel cap near 0.

Theorem A is stated for knots j/48. Its proof uses only a finite piecewise-linear phase, with |d(u)| ≤ TΛu for Λ the largest slope, so it applies verbatim to these knots. Bakeas confirmed this. It also assumes Re γ ≥ 0. Here a = Re γ = 0.1118, and both programs check a > 0.

The claim is max(T,F) ≤ T = 0.6863607950. Two Arb verifiers give F ≤ 0.6863607742 and F ≤ 0.6863607773.

## What's here

- `certificates/C42_compact_best.json`: the profile, with exact modulus T. SHA-256 `ac8d052d567069c6f2b32583791515628fde0fdfe3a1c15e8b2967d45e898765`.
- `verifier_versions/verify_volterra_7134ee6184e563f2.py`: a verifier written from Theorem A and Lemma 3.
- `verifier_versions/check_bakeas_local_3d07e55ffcc2e49a.py`: Bakeas's `independent_check.py`, adapted. `BAKEAS_ADAPTATION.md` lists the changes, and `third_party/check_bakeas_local.diff` is the full diff from version 1.0.0 of the original.
- `verifier_versions/GUARD_CHANGE.diff`: the only change since the versions that Bakeas and the review read. Both programs now refuse to run under `python -O`, which would remove their `assert` checks. The adapted checker also gets the finite checks and the non-zero exit on failure from Bakeas's version 1.1.0, and a corrected header.
- `VOLTERRA_VERIFIER_MATH.md`: the derivation, including both new bounds.
- `logs/`: receipts of both final runs.
- `review/`: a separate AI review's direct Volterra solve.

File names carry the first 16 hex digits of each file's SHA-256. `SHA256SUMS` covers the folder.

## Reproduce

Needs python-flint 0.9.0, numpy, scipy and mpmath.

    ./reproduce_best.sh

It took about 40 minutes on 32 vCPUs and needs up to about 95 GB of RAM. `reproduce_best.sh` is what produced the receipts in `logs/`. Set `PYTHON` to change the interpreter.

## What this doesn't show

The two verifiers are separate implementations, but they share the new propagation bounds, so a mistake there wouldn't be caught by running both.

Bakeas checked the bundle and found no error. They re-derived every new bound in `VOLTERRA_VERIFIER_MATH.md`, read `verify_volterra` line by line against it, and read the diff of the adapted checker against their own. They didn't have the memory for the full N = 14,155,776 run. At N = 294,912 and 589,824 the verifier gave F ≤ 0.6864311 and F ≤ 0.6863783, unproved as expected, and h² scaling to the full N predicts F ≤ 0.68636077, which matches the receipts. A separate floating-point solver of theirs converges to F ≈ 0.686360744.

The review re-derived Theorem A and the new bounds, found no mathematical gap, and got F ≈ 0.68636074 from its own solve. That number is a diagnostic, not a certificate, and the review is AI, so it isn't a second party.

This says nothing about optimality. A wider search for anything below this family's optimum has found nothing so far.

## AI use

Analysis and proof drafting by AI agents under my direction; I reviewed the results at summary level. AI agents also prepared the certificate, verifiers and this text. This does not claim a line-by-line human verification.

Code MIT. Bakeas's MIT notice stays with the adapted checker, in `third_party/LICENSE.bakeas-turan-c42-bound`.
