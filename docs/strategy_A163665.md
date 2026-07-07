# Computational Strategy — OEIS A163665

## 0. Inputs And Confirmation Handoff

Inputs:

- OEIS snapshot: `A163665 - OEIS.md`;
- external triage: `input_A163665.md`;
- b-file: `data/b163665.txt`;
- local methodological PDFs: `Webster-Purdum.pdf`, `rpb272-final.pdf`;
- confirmation dossier: `docs/confirmation_A163665.md`;
- NotebookLM meta-compendium notebook: `OEIS_META_COMPENDIUM`
  (`4f1bfbce-3e9f-4a35-9548-aadbfa2d7d04`).

Confirmation verdict: `CONFIRMED-GO`.

Handoff constraints:

- target predicate: find `n` such that `phi(sigma(n)) = n`, then report
  `k = prime(n)`;
- old bound: no further `k <= 5*10^9`, corresponding to
  `n <= pi(5*10^9) ~= 234954223`;
- confirmed pruning: skip odd `n > 1`;
- unconfirmed pruning: divisibility by 4, 6, 12, or stronger congruence classes;
- expected bottleneck: exact factorization of `sigma(n)` after computing
  `sigma(n)`;
- architecture ceiling: one segmented C17 source, one Makefile, fixed range
  contract, one manifest and one hit file per segment, no external verifier or
  textual verdict parser in the hot loop;
- current authorization: strategy and bounded smoke tests only; production
  extension above the old bound requires benchmark evidence and Carlo approval.

## 1. Meta-Compendium Consultation

NotebookLM consultation status:

- source list succeeded; `OEIS_META_COMPENDIUM` has 20 ready sources;
- first `notebooklm ask` attempt stalled without output and was interrupted;
- compact retry succeeded outside the sandbox;
- retry conversation ID: `7518c18f-25f7-4c12-abb1-5f9c9f90a804`;
- answer used as methodological memory, not as mathematical proof.

NotebookLM precedents judged reusable:

- use a filter/reject cascade before high-cost arithmetic, analogous to
  A138044 and A065992;
- use segmented scanning with atomic per-segment manifests;
- perform mandatory b-file/known-term replay before extension claims;
- require thread-count invariance if parallelism is enabled;
- use `uint64_t` storage with `__uint128_t` intermediates and static numeric
  assertions;
- make pilot throughput a veto gate before campaign authorization;
- avoid undocumented discovery-probability claims.

NotebookLM precedents judged unsafe or weakened for A163665:

- bottom-up or inductive generation analogies from A096242/A224540 do not
  transfer cleanly because `n -> phi(sigma(n))` is not reversible enough;
- digit/base filters from A029483/A359813 are irrelevant;
- a suggested "skip unresolved factorization" policy is not acceptable for
  certified negative ranges. If `sigma(n)` cannot be factored exactly, the
  segment is unresolved, not certified.

Influence on this strategy:

- choose a top-down segmented scan, not inverse generation;
- keep all hot-loop arithmetic local and exact;
- design manifests to certify empty intervals;
- make the first implementation production-track `v1`, with benchmark gates
  before any long run;
- keep the convenience index blocked until `P_hit` and wall time are grounded.

Local PDF review:

- `Webster-Purdum.pdf` and `rpb272-final.pdf` concern algorithms for the
  multiplication-table problem `M(n)`, not the map `n -> phi(sigma(n))`.
  They do not provide a direct pruning theorem for A163665.
- Reusable methodological lessons: segmentation is a first-class design choice,
  cache/space behavior can dominate asymptotic simplicity, wheel/modulo methods
  may buy speed only when their residue rules are mathematically valid, and
  independent implementations are useful for validating exact computations.
- Non-transferable lessons: multiplication-table rectangle geometry,
  `delta(n)`, `tau+(n)` bounds, and Monte Carlo estimation do not certify
  negative ranges for A163665.
- Strategy impact: add smoothness/factorization telemetry and keep a future
  wheel/modulo benchmark as an optional post-`v1` refinement, not a production
  requirement or current pruning rule.

## 2. Computational Problem Synthesis

The search variable is the prime index `n`, not the prime value `k`.

For each candidate `n`:

1. compute the exact factorization of `n`;
2. compute `S = sigma(n)` from the factorization;
3. compute `phi(S)` from the exact factorization of `S`;
4. if `phi(S) == n`, record an index hit and later compute `k = prime(n)`.

The useful scientific outputs are:

- new terms `k = prime(n)`, if any;
- a certified negative extension of the known bound;
- segment manifests supporting replay and publication audit.

The known hits have highly smooth-looking indices, but this is not yet a proof
or a safe filter. The scan must remain complete over all even `n` in each
claimed segment.

## 3. Engineering Surface Budget

| Component | Scientifically necessary? | Trust boundary or state count | Main failure mode | Simpler alternative |
|---|---:|---:|---|---|
| C17 predicate kernel | yes | 1 trusted codebase | formula bug, overflow | none |
| Segmented factorization of `n` | yes | internal state arrays | missed prime power, segment boundary error | slower direct factoring |
| Factorization of `sigma(n)` | yes | internal deterministic 64-bit routines | incomplete factorization, rare Pollard-Rho failure | trial division only for small probes |
| Deterministic Miller-Rabin for 64-bit | recommended | internal function | wrong base set, overflow | full trial division, too slow |
| Pollard-Rho for 64-bit cofactors | recommended | internal function with fixed seeds | nontermination or repeated failure | mark segment failed, never skip |
| Prime table | yes | generated local data | insufficient limit | static assertion and manifest limit |
| `prime(n)` conversion | cold path only | optional external cross-check | wrong reported OEIS term | report `n` plus two independent `prime(n)` checks |
| CLI modes | yes, minimal | 3 modes: selftest, segment, summarize | mode drift | compile-time constants, no general framework |
| Threads | optional | 1 or fixed N states | race in hit/log output | segment-level parallelism only |
| Manifest | yes | one durable artifact per segment | incomplete negative certificate | no publication-grade negative claim |
| Checksums | yes | manifest fields | artifact mismatch | weaker audit trail |
| Retry/timeout machinery | no | 0 automated states | hidden partial reruns | explicit manual rerun of failed segment |
| Textual verdict parser | no | 0 | parser ambiguity | prohibit |
| External services | no | 0 in hot loop | transport/auth failure | prohibit |

Engineering surface class: `MODERATE`, bounded and acceptable.

## 4. Minimal Architecture And Complexity Veto

### Minimal Fixed Production Workflow

Expected files:

- planned versioned production source, normalized in the public package as
  `src/a163665.c`;
- `src/Makefile`.

Supported modes:

- `selftest`: known index hits, b-file mapping, overflow-guard tests, small
  independent range;
- `segment --lo N --hi M --threads T --out PREFIX`: inclusive segmented scan;
- `summarize --manifest-glob ...`: check completed segment continuity and hit
  list integrity.

Durable artifacts per segment:

- `PREFIX.manifest.tsv`;
- `PREFIX.hits.tsv`;
- optional `PREFIX.stderr.log` only for diagnostics, not as a parsed verdict.

Manifest minimum fields:

```text
seq_id, version, git_hash_or_source_sha256, lo, hi, parity_policy,
threads, segment_size, candidates_raw, candidates_tested, hits,
unresolved_count, started_utc, ended_utc, wall_seconds,
binary_sha256, hit_file_sha256, status
```

Completion statuses:

- `OK_EMPTY`;
- `OK_HITS`;
- `FAILED`;
- `UNRESOLVED`.

Only `OK_EMPTY` and `OK_HITS` certify coverage. No automatic retry tree is
included in `v1`.

### Generalized Monolithic Implementation

A broader framework with selectable algorithms, arbitrary numeric ranges,
external nth-prime tools, retry policies, and multiple parser states is rejected.
It would add engineering burden without improving the immediate mathematical
claim.

Architecture ceiling:

- one source, one Makefile;
- at most three CLI modes;
- no external parser/verifier in the hot loop;
- no more than four manifest statuses;
- fixed segmented inclusive ranges;
- exact 64-bit arithmetic only, with `__uint128_t` intermediates;
- failed or unresolved segments are not certified.

Complexity gate: `PASS`.

## 5. Authorization Envelope And Pre-Proposal Gate

Extracted authorization envelope:

| Limit | Value |
|---|---:|
| Unconditional production scan above old bound | none |
| Strategy smoke-test range | at most `1e7` candidate indices per probe |
| Strategy smoke-test wall time | at most 30 minutes total |
| Old bound replay | not authorized until projected by benchmark |
| Extension to `n <= 1e10` | outside current authorization |
| Extension to `n <= 1e11` | outside current authorization |
| Extension to `n <= 1e12` | outside current authorization |
| Expansion authority | Carlo only |

Pre-proposal gate:

| Proposed action | Estimated range | Estimated wall time | Estimated success probability | Triage limit | Authorized now? | Unlock condition |
|---|---:|---:|---:|---:|---|---|
| Implement production-track `v1` | no production scan | development only | not applicable | architecture ceiling | yes | keep complexity gate passing |
| `selftest` and known-term replay | small fixed data | minutes expected | deterministic | 30 min smoke total | yes | none |
| one benchmark segment near old bound | `<= 1e7` indices | `<= 30 min` | not a hit campaign | 30 min smoke total | yes, if bounded | exact command and manifest |
| full old-bound replay to `pi(5e9)` | `~2.35e8` indices | unknown | deterministic | not authorized | no | benchmark projection plus Carlo approval |
| extension to `n <= 1e10` | about 40x old index bound | unknown | not estimable | outside envelope | no | benchmark, validation, runtime budget, Carlo approval |
| extension to `n <= 1e11` | about 400x old index bound | unknown | not estimable | outside envelope | no | separate explicit approval |

Authorization gate: `PASS` for implementation and bounded smoke tests only.

Production campaign gate: `BLOCKED BY TRIAGE` until benchmark evidence and Carlo
approval exist.

## 6. Recommended Main Algorithm

Use a segmented exact-factorization scan over even `n`.

For a segment `[lo, hi]`:

1. normalize `lo` to include only `n = 1, 2` special cases or even `n > 1`;
2. allocate arrays for the represented even candidates:
   - `rem[i] = n`;
   - `sigma[i] = 1`;
   - optional compact status/diagnostic counters;
3. generate primes up to a verified limit at least `sqrt(hi)` for factoring
   `n`, and up to the selected `sqrt(max_sigma_bound)` if trial division is
   used before Miller-Rabin/Pollard-Rho;
4. for each prime `p`, visit multiples of `p` in the segment, extract the
   exponent `a`, update `rem`, and multiply `sigma` by
   `(p^(a+1)-1)/(p-1)` using `__uint128_t` guarded arithmetic;
5. after the prime pass, if `rem[i] > 1`, multiply `sigma[i]` by `rem[i] + 1`;
6. factor `S = sigma[i]` exactly:
   - trial divide by a small-prime prefix;
   - as soon as a prime factor `q` of `S` is found, reject if `(q - 1)` does
     not divide `n`;
   - finish remaining cofactors with deterministic Miller-Rabin and fixed-seed
     64-bit Pollard-Rho;
7. compute `phi(S)` exactly and compare to `n`;
8. write hits as index witnesses:

```text
n, sigma_n, phi_sigma_n, factor_n, factor_sigma_n
```

9. convert `n` to `prime(n)` only after a hit, preferably with two independent
   cold-path checks before public reporting.

The early `(q - 1) | n` reject is safe because every prime factor `q` of an
integer `S` with `phi(S) = n` contributes a factor `q - 1` to `phi(S)`.

## 7. Alternative Strategies

Segmented divisor-sum sieve:

- advantage: conceptually simple;
- problem: many divisor additions over very large ranges;
- verdict: useful as a small independent validator, not the production kernel.

Direct on-the-fly factoring of each `n`:

- advantage: smallest code surface;
- problem: repeats work and loses segment locality;
- verdict: acceptable only as a selftest validator over small ranges.

Inverse-totient-assisted generation:

- advantage: might reduce candidate volume if a strong theory emerges;
- problem: `m = sigma(n)` still has to be matched to `n`, and inverse totient
  branching is not obviously smaller;
- verdict: not production `v1`.

Smooth-family search:

- advantage: known hits have smooth-looking index factorizations;
- problem: incomplete and publication-risky;
- verdict: may inform heuristics, never certify negative intervals.

Wheel/modulo refinement:

- advantage: the multiplication-table PDFs show that residue-class generation
  can produce major speedups in sieve-like computations;
- problem: A163665 has no proven residue-class exclusion beyond even `n > 1`,
  and importing mod-6, mod-12, mod-60, or similar wheel rules would be a false
  analogy;
- verdict: not `v1`; benchmark only after a rigorous congruence proof or as a
  purely diagnostic implementation with no certification role.

## 8. Data Structures And Numeric Types

Primary types:

- `uint64_t` for `n`, `sigma(n)`, `phi`, prime factors, counters where bounded;
- `__uint128_t` for multiplication, prime powers, sigma products, and
  Pollard-Rho modular multiplication;
- `uint32_t` acceptable only for prime table entries below `2^32`, never for
  segment endpoints or global counters.

Numeric-domain ledger:

| Quantity | Planned type | Guard |
|---|---|---|
| `n` | `uint64_t` | `hi <= configured_n_max` |
| `sigma(n)` | `uint64_t` | analytic/configured bound plus runtime overflow check |
| intermediate `p^(a+1)` | `__uint128_t` | checked before cast |
| `phi(sigma(n))` | `uint64_t` | computed by divide-then-multiply where possible |
| segment counters | `uint64_t` | manifest consistency checks |
| modular products | `__uint128_t` | required for 64-bit Miller-Rabin/Pollard-Rho |

Supported maxima should be compile-time or startup-validated. Numeric type
capacity does not authorize a campaign range.

## 9. Hardware Feasibility

Expected platform compatibility: 64 GB installed RAM and a 16-thread CPU class
are adequate for bounded segment scans, provided segment length is tuned.

The likely pressure points are:

- CPU time in exact factorization of `sigma(n)`;
- cache behavior of `rem`/`sigma` arrays;
- load imbalance if Pollard-Rho is needed frequently;
- manifest overhead if segments are made too small.

Initial segment size should be chosen empirically, not hard-coded from triage.
Start with a small range such as `2^20` to `2^22` represented even candidates,
then measure.

## 10. Throughput, Chunking, And Telemetry

Segment contract:

- ranges are inclusive `[lo, hi]`;
- `lo` and `hi` are recorded exactly in every manifest;
- odd `n > 1` are excluded by policy and counted separately from tested
  candidates;
- segment artifacts are immutable after `OK_*` completion;
- failed or unresolved segments must be rerun from scratch.

Required telemetry:

- raw range width;
- even candidates represented;
- candidates rejected by small factor `(q - 1) | n` test;
- candidates requiring Miller-Rabin;
- candidates requiring Pollard-Rho;
- largest prime factor profile for `n` and sampled/final `sigma(n)` values;
- counts by factorization path, including trial-division-only, Miller-Rabin, and
  Pollard-Rho-assisted completions;
- unresolved/failure count;
- hit count;
- wall seconds and candidates per second;
- thread count and segment size;
- artifact checksums.

Chunking recommendation:

- use fixed-size segments for initial validation;
- allow only a compile-time or documented run-time segment length parameter;
- do not split a campaign to evade cumulative authorization limits.

## 11. Transfer-Risk And Performance Audit

Risks imported from prior projects:

- 64-bit fit is not equivalent to feasible runtime;
- b-file replay can still share a common-mode formula bug unless an independent
  small validator is present;
- hit density is not stable enough to estimate `P_hit`;
- manifest language must distinguish proven empty segments from failed or
  unresolved segments;
- thread determinism is mandatory if hits can be written concurrently.
- multiplication-table algorithms support the general engineering lesson that
  segmentation, residue classes, and factored inputs can matter, but their
  specific geometry and Monte Carlo methods do not transfer to this exact
  fixed-point search.

Audit actions before campaign approval:

- independent Python/SymPy or PARI/GP small-range cross-check up to at least
  `n = 100000`;
- C selftest against all 20 b-file terms and their index witnesses;
- C vs independent validator on random small intervals;
- one high-region timing segment near `n ~= 2.35e8`;
- overflow tests at the configured maximum range;
- 1-thread vs N-thread equality on the same segment.
- inspect smoothness telemetry from the timing segment before choosing larger
  campaign chunks; this may tune segment size but must not become a skip rule.

## 12. Minimum Validation Plan

Pre-build or build-time:

- compile with warnings enabled;
- run sanitizer selftests on small ranges;
- verify `phi(sigma(n)) = n` for all known index hits;
- verify b-file terms map to `prime(n)`;
- verify no extra hits in a small fully enumerated range.

Pre-campaign:

- reproduce a bounded interval ending at the old OEIS bound, if benchmark
  authorizes it;
- run segmented and monolithic-small modes over the same small ranges and
  compare identical hit sets;
- record source and binary checksums in manifests.

Post-hit:

- recompute `phi(sigma(n)) = n` independently;
- compute `prime(n)` using two independent tools or one internal plus one
  external method;
- store factorization witnesses for `n` and `sigma(n)`.

Negative-result certification:

- only contiguous `OK_EMPTY`/`OK_HITS` manifests certify coverage;
- failed or unresolved segments create holes and must be excluded from any
  claimed bound.

## 13. Convenience Index

A numeric Type A convenience index is not responsibly computable yet.

Reasons:

- `P_hit` is not estimable from the current sparse and irregular data;
- `T_wall` is unknown until the exact `sigma(n)`/`phi(sigma(n))` kernel is
  benchmarked;
- the triage's `bound_factor` suggestions are candidate ambitions, not
  authorized campaigns;
- the meta-compendium explicitly warns against undocumented discovery
  probabilities.

Required inputs for the first valid index computation:

| Scenario | Bound factor candidate | Required measured input |
|---|---:|---|
| conservative | about `2x` to `40x` old index bound | exact projected wall time and Carlo-approved `P_hit` prior |
| recommended | about `40x` to `400x` | same, plus validation rerun cost |
| ambitious | `400x+` | same, plus explicit high-runtime authorization |

Until those exist, the convenience-index verdict is:

```text
DEFERRED — do not use a numerical score for campaign authorization.
```

This is a deliberate veto, not an omission.

## 14. Final Strategy Verdict And Next Action

Strategy verdict: `STRATEGY-GO` for production-track `v1` implementation and
bounded smoke tests.

Campaign verdict: production extension is `BLOCKED BY TRIAGE` until the
authorization gate is rerun with benchmark data.

Recommended next action:

```text
Implement the minimal segmented C17 production artifact described above
(planned at this stage under a versioned source name and normalized in the
public package as `src/a163665.c`) together with `src/Makefile`. Run only
selftest, small independent validation, and one bounded timing probe before
requesting approval for any old-bound replay or extension campaign.
```

Non-negotiable constraints for `v1`:

- segmented/checkpointable execution from the start;
- exact factorization or unresolved segment failure;
- no external parser/verifier in the hot loop;
- no generalized framework;
- no campaign recommendation from unmeasured throughput;
- no numerical convenience index until `P_hit` and `T_wall` are defensible.
