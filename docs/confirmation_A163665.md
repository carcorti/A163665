# Triage Confirmation — OEIS A163665

## 0. Input Summary

### 0.1 OEIS Snapshot

[OEIS] A163665 is defined as the numbers `k` such that

```text
k = prime(phi(sigma(pi(k)))).
```

The listed terms are:

```text
2, 3, 19, 37, 719, 1511, 5443, 69709, 386093, 1907819,
10777931, 17819101, 18653749, 56125547, 60163267, 98911811,
272887613, 567611663, 989060309, 2444540149
```

[OEIS] The comments state that all terms are primes and that there is no further
term up to `5*10^9`. The Mathematica line searches prime indices up to
`235000000`, consistent with `pi(5*10^9) ~= 234954223`.

### 0.2 External Triage

The Gemini triage gives a favorable `PROMETTENTE` verdict. Its main proposal is
to invert the search variable: write `n = pi(k)`, test
`n = phi(sigma(n))`, and output `k = prime(n)` only for hits. It identifies the
problem as a rare-hit Type A search, expects native `uint64_t` arithmetic to be
sufficient for the proposed exploratory ranges, and proposes a segmented
computation with careful throughput measurement before any extended campaign.

### 0.3 Consistency Check

Confirmed. The OEIS definition is equivalent to the triage's working predicate:

```text
k = prime(phi(sigma(pi(k))))
n = pi(k)
k = prime(n)
therefore n = phi(sigma(n)).
```

Small diagnostic checks recovered the first index hits

```text
1, 2, 8, 12, 128, 240, 720, 6912, 32768
```

which map to the first OEIS terms via `prime(n)`. SymPy checks also verified
the full local b-file list against `phi(sigma(pi(k))) = pi(k)`.

No conflict was found between the OEIS snapshot, b-file, and triage.

### 0.4 Bibliographic Search

No external bibliographic search has been integrated at this stage. A focused
search is recommended before publication work, not before strategy design:

- OEIS history and author notes for A163665 and A163666;
- prior computations of fixed points of `phi(sigma(n))`;
- implementation precedents for segmented multiplicative-function evaluation.

Candidate references should not be imported into the project record until Carlo
approves them.

## 1. Mathematical Confirmation

### 1.1 Formal Definition

[OEIS] A term is a prime `k = p_n` whose index `n` satisfies:

```text
n = phi(sigma(n)).
```

The computational task is therefore to find index hits `n`, then convert each
hit to `p_n`.

### 1.2 Mathematical Nature

[PLAUSIBLE] This is a sparse fixed-point search for the composite arithmetic
map `n -> phi(sigma(n))`. It is not a primality search in the hot loop. The
prime `k` is a reporting consequence of the index equality, not the primary
object to test.

[OEIS] The known finite data and the `no further term up to 5*10^9` comment
make the immediate scientific objective a certification extension beyond the
2009 bound, with or without a new hit.

### 1.3 Key Arithmetic Properties

[OEIS] Every reported `k` is prime because `k = prime(n)`.

[PLAUSIBLE] For `n > 1`, `sigma(n) > 2`, and Euler's totient is even for every
integer greater than 2. Therefore any solution with `n > 1` must have even
`n`. This is a rigorous half-domain reduction.

[TO-VERIFY] Stronger congruence restrictions are not established. A small check
up to `n = 100000` found hits with residues incompatible with simple mandatory
divisibility by 4, 6, or 12, so the triage's suggested follow-up modular
filters should be treated as open, not as available pruning.

[PLAUSIBLE] If `q` divides `sigma(n)` and `phi(sigma(n)) = n`, then `q-1`
divides `n`. This gives a useful validation and possible pruning invariant
after factorization of `sigma(n)`, but it does not by itself remove the need to
compute or factor `sigma(n)`.

### 1.4 Computational Reduction

[OEIS] The Mathematica line already searches over the prime index variable.
The correct production problem is:

```text
scan n, preferably even n only after n = 2;
compute sigma(n);
compute phi(sigma(n));
if phi(sigma(n)) == n, report prime(n).
```

[PLAUSIBLE] The best initial production kernel is not a naive divisor-sum loop.
It should factor all `n` in a segment, compute `sigma(n)` from the factorization,
then factor `sigma(n)` with deterministic 64-bit routines to compute `phi`.

### 1.5 Qualitative Complexity

Mathematical difficulty: `MODERATE`.

The problem is easy to define and validate locally, but no density law or
strong constructive structure is visible from the snapshot. New-hit probability
is not responsibly estimable from the current data.

Computational difficulty: `MODERATE`.

The loop is local, parallelizable, and native-width for the proposed ranges, but
the decisive throughput depends on the cost of repeated exact factorization of
`sigma(n)`.

### 1.6 Pruning Levers

Confirmed or immediately usable:

- scan index `n`, not candidate prime `k`;
- exclude odd `n > 1`;
- report `prime(n)` only after a confirmed hit;
- use known-term replay as a deterministic regression suite.

Open or secondary:

- modular restrictions beyond parity;
- inverse-totient generation as a replacement for scanning;
- precomputed or cached factorization of `sigma(n)` components;
- special handling of smooth `n` families observed among known hits.

### 1.7 Preliminary Feasibility Judgment

The external triage is confirmed in direction but weakened in scale. Extending
the 2009 bound is plausible and worth strategy design. Claims that ranges such
as `n <= 10^11` or `n <= 10^12` are easily affordable must wait for measured
throughput on the exact kernel.

### 1.8 Software And Operational Complexity

Dependency inventory:

| Item | Unavoidable? | Trust boundary | Notes |
|---|---:|---:|---|
| C17 production executable | yes | local compiled code | main trusted artifact |
| Prime table up to `sqrt(max_sigma)` | yes | internal generated data | deterministic and checksumable |
| `prime(n)` reporting tool/library | yes, but cold path only | external or internal | may be isolated from search kernel |
| Textual verdict parser | no | none | prohibit in v1 |
| External long-running service | no | none | prohibit |
| Checkpoint/resume manifests | yes for production campaign | durable local state | must be small and structured |
| Retry/conflict machinery | no for v1 | none | failed segment should be rerun explicitly |

Software/operational complexity class: `MODERATE`.

The intrinsic sequence is computationally nontrivial, but the engineering
surface can remain small if the production artifact is a single segmented C17
program with fixed modes and no external verifier/parser loop.

### 1.9 Minimal Viable Computation And Architecture Ceiling

Smallest credible computation:

1. replay all known index hits and the OEIS b-file correspondence;
2. verify no hits in a small interval above the last known hit;
3. benchmark one fixed segment near the old bound using the intended exact
   kernel;
4. only then decide whether a production extension beyond `pi(5*10^9)` is
   authorized.

Prohibited in the first production design:

- multiple algorithm families exposed as CLI modes;
- external textual parsers or verdict labels;
- arbitrary user-selected arithmetic domains;
- automatic retry trees;
- unbounded checkpoint formats;
- generalized campaign framework.

Architecture ceiling:

- one production C17 source, originally expected during planning under a
  versioned source name and normalized in the public package as `src/a163665.c`;
- one Makefile, expected as `src/Makefile`;
- segmented inclusive `n` ranges;
- one manifest per segment;
- one hit file per segment, usually empty;
- deterministic known-term and segment-vs-monolithic regressions;
- optional cold-path prime-index conversion isolated from the hot predicate.

Production-track `v1` can plausibly start after smoke tests and CLI timing
measurements. No family of exploratory prototypes is required.

### 1.10 Pre-Development Stop Criteria

Stop before production implementation if any of the following occurs:

- the exact segmented factorization kernel cannot replay the b-file terms;
- `sigma(n)` or `phi(sigma(n))` requires a numeric domain wider than `uint64_t`
  for the selected range;
- computing `phi(sigma(n))` needs an external factorization service;
- a safe design requires textual parsing, retry/conflict machinery, and complex
  checkpoint state;
- a benchmark near the old bound implies that even a modest extension, including
  validation reruns, exceeds Carlo's approved runtime budget;
- no compact segmented artifact contract can be specified.

### 1.11 Authorization Envelope

No production search campaign has been quantitatively authorized by the input
files. The following envelope is therefore a hard constraint for the strategy
stage:

| Class | Range / runtime | Status |
|---|---:|---|
| Confirmation diagnostics | small local checks only | completed |
| Strategy smoke tests | known-term replay and one or two bounded timing probes, each no more than `1e7` candidate indices and no more than 30 minutes total wall time | conditionally acceptable for strategy design |
| Production extension beyond old bound | any cumulative scan above `n = pi(5*10^9)` | not yet authorized |
| Ranges `n <= 1e10`, `1e11`, or `1e12` | triage suggestions only | outside current authorization until benchmark evidence and Carlo approval |

Required unlock evidence:

- exact kernel benchmark in the target numeric region;
- known-term replay and negative-range replay matching OEIS;
- projected cumulative wall time including validation reruns;
- Carlo's explicit approval for the selected cumulative range and runtime.

Only Carlo may authorize expansion beyond this envelope.

## 2. Triage Verdict

### 2.1 Confirmed Points

- The index-variable reduction is correct.
- The OEIS snapshot, b-file, and external triage are mutually consistent.
- All nontrivial index hits must be even.
- The problem has a credible native-integer, local, segmented route.
- The engineering surface can be kept below the prohibitive threshold.

### 2.2 Rejected Or Weakened Points

- The claim that very large ranges are cheaply reachable is not confirmed
  before benchmark evidence.
- Stronger modular pruning than parity is not available yet.
- A naive segmented divisor-sum sieve is not automatically the most efficient
  route; segmented factorization of `n` plus 64-bit factorization of
  `sigma(n)` is the stronger baseline candidate.

### 2.3 Open Verification Points

- exact throughput near and above the old bound;
- best hot-loop choice among segmented factorization, divisor-sum sieving, and
  inverse-totient-assisted pruning;
- whether smoothness patterns among known hits can produce a safe filter;
- exact supported maximum range under `uint64_t` after bounding `sigma(n)`.

### 2.4 Recommendation For Strategy Stage

Verdict: `CONFIRMED-GO`.

The sequence is suitable for Computational Strategy Design. The strategy must
query `OEIS_META_COMPENDIUM`, then choose the smallest segmented production
architecture. A production campaign must remain blocked until benchmark-based
authorization passes.

Strategy handoff:

```text
exact computational target: find n > pi(2444540149) and especially n > pi(5*10^9)
such that phi(sigma(n)) = n; report k = prime(n)
unresolved mathematical risks: no density model; no confirmed pruning beyond parity
expected bottleneck: exact factorization of sigma(n) after computing sigma(n)
candidate scale: old bound n ~= 2.35e8; larger ranges require benchmark unlock
required validation: b-file replay, known negative interval replay, independent
  small-range cross-check, segmented = monolithic regression
software/operational complexity class: MODERATE
external dependency/trust-boundary count: 0 in hot loop; 1 optional cold-path prime(n) tool
allowed architecture ceiling: one C17 source, one Makefile, fixed segmented range
production-track v1: can start after smoke tests and CLI measurements
forbidden generalizations: multi-mode framework, external verdict parsers,
  arbitrary domains, automated retry trees
parser/subprocess complexity: not acceptable in the hot loop
pre-development kill criteria: see Section 1.10
unconditional campaign limits: no production search above old bound authorized yet
conditional limits: bounded smoke tests up to 1e7 candidate indices and 30 minutes total
unlock evidence: exact kernel throughput, replay validation, wall-time projection,
  explicit Carlo approval
query OEIS_META_COMPENDIUM: yes, before final strategy
```
