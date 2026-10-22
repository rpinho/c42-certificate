# Final live check — C42

Final refresh before publication: requests initiated from 2026-10-22T17:25:12.822156+00:00 to 2026-10-22T17:25:13.783806+00:00 (UTC).
An initial pass at 17:17 UTC gave the same relevant bounds and repository heads.

**Outcome: no bound below 0.6888863712 was located in the sources checked.**
This is a bounded public-source audit, not exhaustive priority certification.

| Source | Live result |
| --- | --- |
| [Tao's C42 page](https://teorth.github.io/optimizationproblems/constants/42a.html) and [raw source](https://github.com/teorth/optimizationproblems/blob/46dbf03669219cb2bcf33fd5bbf65209d137f994/constants/42a.md) | Griego's 0.6906538* remains the best listed upper bound. Head `46dbf03669219cb2bcf33fd5bbf65209d137f994`. |
| [Tao's open PRs](https://github.com/teorth/optimizationproblems/pulls) | All 29 open PRs fetched (page limit 100). The relevant pending proposal is [PR184](https://github.com/teorth/optimizationproblems/pull/184), Röhrig's C42 bound 0.688983. No stronger C42 proposal found. A separate issue/PR search returned 8 matches, inspected for relevance. |
| [Plata's README](https://github.com/krplatz/turan-c42-six-level-review/tree/d3650caf944e8602e867e2e81a439e037a7b91a9) and [proof](https://github.com/krplatz/turan-c42-six-level-review/blob/d3650caf944e8602e867e2e81a439e037a7b91a9/proof.md) | Still proposes 68897/100000 = 0.688970; explicitly seeks independent mathematical review. Head `d3650caf944e8602e867e2e81a439e037a7b91a9`. This is the strongest located public comparison, and is absent from Tao's current bound table. |
| [Röhrig's repository](https://github.com/AndreRoehrig/optimization-constants-certificates/tree/ed6387399cf6de02eaa4cf10c43787633496619a) | C42 proposal remains 0.688983. Head `ed6387399cf6de02eaa4cf10c43787633496619a`. |
| [GitHub repository search](https://github.com/search?q=turan-c42&type=repositories&s=updated&o=desc) | `turan-c42`, sorted by update, returns only the Plata and Griego repositories before this publication. |
| [arXiv API](https://export.arxiv.org/api/query?search_query=all:%22pure%20power%20sum%22%20OR%20(all:Turan%20AND%20all:%22power%20sum%22)&sortBy=submittedDate&sortOrder=descending&max_results=100) | All 4 results inspected by title and abstract; none gives a competing C42 bound. Returned IDs: 2112.11228v2, 1008.4535v3, 0706.4131v3, math/0607238v2. These cover other power-sum settings or an unrelated ξ-function representation. |
| [arXiv web search](https://arxiv.org/search/?query=Turan+pure+power+sum&searchtype=all&abstracts=show&order=-announced_date_first&size=50) | One result, Andersson's 2006/2007 inf-max construction, concerning different modulus and exponent constraints. No newer competing bound found. |

Supplemental web queries included `"Turan" "0.688886"`,
`"power sum" "2026" "upper bound"`, `site:arxiv.org "Turán" "power sum"`,
`site:arxiv.org "pure power" "2026"`, `"0.688970"`,
`"Turan" "0.6888"`, and `"Turán" "Griego"`.
They yielded no additional relevant smaller bound. Search indexing is incomplete;
this is why live repository and arXiv API reads were also used.

All 13 direct requests succeeded. Response SHA-256 receipts are in
[verification/live_source_receipts.json](verification/live_source_receipts.json).
Full source snapshots were retained in the local publication workspace.
No private sources, author contact, or deposit were used for this audit.

The proposed exact value is 215276991/312500000 = 0.6888863712.
Compared with Plata's unrefereed 0.688970, the improvement is exactly
0.0000836288. The current table convention calls for retaining historical
rows, appending an asterisked claim, and updating the README's bound and recent
progress; the absent 0.688970 comparison is explicitly cited in the PR text.

The arithmetic publication gate is separate from this literature check:
clean-copy interval and Arb runs and their rejection controls are recorded in
[verification/](verification/). The analytic proof remains unrefereed;
there is no optimality, formalization or exhaustive priority claim.
