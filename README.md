# A 256-step certificate for Turán's pure power-sum constant

Ricardo Pinho · [ORCID 0009-0006-2744-2599](https://orcid.org/0009-0006-2744-2599)

This unrefereed note gives an analytic argument and an exact-input certificate for
**C₄₂ ≤ 215276991/312500000 = 0.6888863712**, where

\[
R_n=\min_{\max_j|z_j|=1}\max_{1\le k\le n}\left|\sum_{j=1}^n z_j^k\right|,
\qquad C_{42}=\limsup_{n\to\infty}R_n.
\]

The construction uses 256 middle blocks in the coefficient-cancellation method
of Sebastian Griego, André Röhrig and Kerrigan Plata. The strongest public
proposal located before publication was Plata's **0.688970**; the improvement
is **0.0000836288**. See [FINAL_LIVE_CHECK.md](FINAL_LIVE_CHECK.md) for the
dated search, primary sources and limits of coverage.

Read the [paper (PDF)](paper/c42.pdf) and [TeX source](paper/c42.tex). The proof
gives the bound for all sufficiently large n without an explicit threshold.
The programs verify the strict finite inequalities; the asymptotic implication
is an analytic proof requiring mathematical review. No optimality, human expert
peer review or proof-assistant formalization is claimed.

## Reproduce

Python 3.12.13, mpmath 1.3.0 and python-flint 0.9.0 were used.
From a fresh checkout:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
shasum -a 256 -c SHA256SUMS
PYTHON=python sh reproduce.sh
```

The script runs the interval checker at 100 decimal places with 240 terms,
the separately written Arb checker at 768 bits with 400 terms, and four
corruptions through each checker (eight expected rejections). These cover a
too-small radius, a zero middle profile that fails the analytic inequality,
repeated thresholds, and the false near-threshold bound 0.6888863709. Runtime
is several minutes, mostly in the interval implementation. Run without Python's
`-O` flag: the interval engine uses assertions and explicitly rejects that mode.

Individual commands, as in the paper:

```sh
python verify_primary.py certificates/c42_256.json
python independent/check_bounds.py certificates/c42_256.json --bits 768 --terms 400
python rejection_controls.py
```

The certificate contains exact rational strings. Its SHA-256 is
`705a729155eb508878d17ebea0322627a0ae02451c044a346a99d2e81c7afcab`.
Floating-point search diagnostics in that file are unused by both checkers.
The outward enclosure for |Y|/D lies in
`[0.688886370983465645822767803129, 0.688886370983465645822767803130]`.
The exact minimum squared profile margin is
`1491623179869746287/5000000000000000000000000000 > 0`.

Archived clean-copy receipts are under [verification/](verification/).
Reproduction writes fresh results to the ignored `reproduction-output/` folder.
`SHA256SUMS` covers all distributed files except itself.

## Credit and provenance

Griego's generating function, final-block cancellation, Newton identities and
asymptotic reduction are the underlying method. Röhrig generalized the middle
profile, and Plata supplied a six-level refinement. This note retains that
mechanism and supplies finer rational profile data. The paper credits all three
and the classical work of Biró and Harcos.

- [Griego, v1.0.0](https://github.com/sebastian-griego/turan-c42-certificate/tree/v1.0.0): 0.6906538.
- [Röhrig, pinned source](https://github.com/AndreRoehrig/optimization-constants-certificates/tree/ed6387399cf6de02eaa4cf10c43787633496619a): 0.688983.
- [Plata, pinned review draft](https://github.com/krplatz/turan-c42-six-level-review/tree/d3650caf944e8602e867e2e81a439e037a7b91a9): 0.688970.

`sources/roehrig/scripts/certify_turan.py` is Röhrig's unchanged public interval
engine; `verify_primary.py` is a command wrapper. `notes/turan.md` preserves his
analytic note. [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) records their
provenance and licensing status. The Arb checker was written by a separate AI
agent from the mathematical notes without reading the search or interval code.
Its C42 mathematical functions are unchanged from the approved work package;
unrelated routines were omitted from this distribution. Neither checker imports
the other or uses an optimizer.

## AI use and author review

Analysis and proof drafting by AI agents under my direction; I reviewed the
results at summary level. AI agents also prepared the certificate, checkers and
publication materials. This does not claim a line-by-line human verification.
The separate AI-agent check is not independent human expert review.

The approved note's historical references to a local, unsubmitted draft describe
its preparation state. This repository makes that note public. The sole paper
edit for publication adds the summary-level AI review disclosure; its mathematical
content is unchanged. The note refers to the original local audit files; the
publication audit is [FINAL_LIVE_CHECK.md](FINAL_LIVE_CHECK.md).

Original repository software and certificate: [MIT](LICENSE).
The approved paper retains its stated CC BY 4.0 license. Third-party material
retains its upstream rights, as described in the notices.
