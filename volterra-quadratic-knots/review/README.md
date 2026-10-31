# Separate AI review: direct Volterra solve

`independent_check.py` imports neither verifier nor the search. It solves the original Volterra equation by Picard iteration, with exact convolution of piecewise-linear interpolants. It's a numerical diagnostic, not an interval certificate.

The only change from the reviewed file is the certificate path on line 17. Run it from this folder with `python independent_check.py volterra`. Needs python-flint, numpy and scipy.

`receipts/` holds its output at N = 73,728 to 2,359,296. The last two give F ≈ 0.6863607440163906 and 0.6863607440081801.
