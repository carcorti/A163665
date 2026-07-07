# A163665 Campaign Final Validation

Date: 2026-07-06

Package congruence re-audit: 2026-07-07

## Verdict

VALIDATED.

The official campaign runs `A163665_run_001` through `A163665_run_015` form a contiguous certified extension from

```text
n = 234954224
```

through

```text
n = 405234954223
```

with no unresolved or failed segments.

## Scope

Excluded from the official synthesis:

- `A163665_run_000`, pilot/overlap run.
- Any interrupted old-binary artifacts.

Included in the official local campaign synthesis:

- `runs/campaign/A163665_run_001.*` through `runs/campaign/A163665_run_015.*`.

These local run artifacts are not included as files in the raw public package;
their aggregate certified content is represented by `validation/run_manifest.tsv`
and this validation summary.

Previous OEIS-style bound:

```text
pi(5*10^9) = 234954223
```

Final campaign bound:

```text
prime(405234954223) = 11776117514717
```

Therefore the campaign certifies the range after the old `k <= 5*10^9` bound up to

```text
k <= 11776117514717
```

with exactly the new hits listed below.

## Aggregate Counts

```text
official runs:          15
certified subsegments:  40500
raw n-width:            405000000000
candidates tested:      202500000000
unresolved segments:    0
failed segments:        0
new hits:               8
script elapsed total:   22:12:58
```

The script elapsed total is computed from the `started_unix` and `ended_unix` fields in the official run logs.

## New Terms

```tsv
n	k	sigma_n	phi_sigma_n	factor_n	factor_sigma_n
2147483648	50685770167	4294967295	2147483648	2^31	3*5*17*257*65537
3889036800	94206075193	16331433888	3889036800	2^9*3^4*5^2*11^2*31	2^5*3*7*11^3*19*31^2
4389396480	106882008719	18377794080	4389396480	2^13*3^7*5*7^2	2^5*3^3*5*19*41*43*127
21946982400	571264143487	94951936080	21946982400	2^13*3^7*5^2*7^2	2^4*3^2*5*19*31*41*43*127
47416320000	1272242879459	204721968000	47416320000	2^13*3^3*5^4*7^3	2^7*3*5^3*11*43*71*127
92177326080	2536961181761	386940247200	92177326080	2^13*3^8*5*7^3	2^5*3^2*5^2*13*43*127*757
133145026560	3715374607207	601662398400	133145026560	2^13*3^6*5*7^3*13	2^6*3^2*5^2*7*43*127*1093
331914240000	9576704442589	1433565580920	331914240000	2^13*3^3*5^4*7^4	2^3*3*5*11*43*71*127*2801
```

Candidate b-file extension:

```text
21 50685770167
22 94206075193
23 106882008719
24 571264143487
25 1272242879459
26 2536961181761
27 3715374607207
28 9576704442589
```

The last discovered term is

```text
a(28) = 9576704442589
```

The campaign also certifies that there is no further term in

```text
9576704442589 < k <= 11776117514717
```

## Validation Checks

Global run-summary validation:

```text
bad=0
segments_ok=40500
segments_unresolved=0
segments_failed=0
total_hits=8
global_lo=234954224
global_hi=405234954223
width=405000000000
```

Submanifest count:

```text
40500 official submanifest files
```

Final executable selftest:

```text
selftest OK: 20 b-file terms, known indices, low-range cross-check, segmented=direct
```

Final high-hit one-point reproduction:

```text
OK_HITS [331914240000,331914240000] tested=1 hits=1 unresolved=0 failures=0
segments_ok=1	segments_unresolved=0	segments_failed=0	total_hits=1	coverage=CONTIGUOUS	lo=331914240000	hi=331914240000
```

The reproduced row was:

```tsv
331914240000	9576704442589	1433565580920	331914240000	2^13*3^3*5^4*7^4	2^3*3*5*11*43*71*127*2801
```

## Public Consolidated Files

```text
validation/run_manifest.tsv
results/A163665_campaign_hits.tsv
data/b163665.txt
data/certified_terms.tsv
```

SHA-256:

```text
7948c44af7e6c692e180db04e0c5a9edb47d880f056f1489a80a1808d62a95ce  validation/run_manifest.tsv
352f0e414ec085030a5268e96ddc1fd6f6728e02df6d1c7d567c3a885ddbab81  results/A163665_campaign_hits.tsv
afd499cd98241e57e76166160cb763a5bce52eb76327c7f2c52828fcb4fb889e  data/b163665.txt
c07764a2724bc6493399abc13c0d21a9bd0e963de54729151a1dd01dc006d368  data/certified_terms.tsv
```

## GitHub Package Validation

The raw GitHub package was validated after copying and public-package normalization.

Passed commands:

```bash
make -C src clean all
make -C src test THREADS=4
make -C src sanitize
gp -q validation/verify_hits.gp
validation/verify_prime_mapping.sh
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=GitHub/paper GitHub/paper/A163665.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=GitHub/paper GitHub/paper/A163665.tex
make -C src clean
```

Data checks:

```text
data/b163665.txt has a real final blank line.
The last eight k-values in data/b163665.txt match data/certified_terms.tsv.
results/A163665_campaign_hits.tsv matches data/certified_terms.tsv, ignoring comment/header lines.
paper/A163665.tex compiles cleanly after two LaTeX passes and produces paper/A163665.pdf.
No scratch directories named tmp, temp, build, dist, .cache, __pycache__, or runs remain under the package root.
```

Independent positive-hit verification:

```text
PARI_GP_ALL_28_HIT_VERIFICATION=OK
PRIME_MAPPING_ALL_28_VERIFICATION=OK
```

These independent scripts verify all 28 consolidated positive terms, including the eight new discoveries. They do not re-certify the negative interval.

GitHub package SHA-256, excluding validation/validation_summary.md itself to avoid a self-referential checksum:

```text
1a1dbe176bc233b499d35a57db7513f2941c99ab9759f177830c9149be99005b  .gitattributes
8c32e9785d7f0919943d6fb4a053dc619bd2521993b4d8ed3885de0fc9b0b6b2  .gitignore
3b4d6b8cc861620cfb2ffd0ad61b048ab13700eb773061c2432e4731c065a7ca  README.md
b5241b675e5d0915f5b3140f794f38e88fd9f75c504967dd7fb6106c9e854013  CITATION.cff
1db7ea53590c61186548979643d1892af163e23c8254a198be6b1e50ae189911  LICENSE
afd499cd98241e57e76166160cb763a5bce52eb76327c7f2c52828fcb4fb889e  data/b163665.txt
c07764a2724bc6493399abc13c0d21a9bd0e963de54729151a1dd01dc006d368  data/certified_terms.tsv
2f8cdca9f510622e6f858e0fd49accf3a70195a12812377f69ce28dd0ac581f0  docs/confirmation_A163665.md
50df4ab13938d7fee53201bf814b7d00c639c4296eac98901ac3d9c674f95305  docs/strategy_A163665.md
f9d09fe8e9d30da75c7d88c63f73949d26d6423087324df86bf2e4e5a4c82603  paper/A163665.pdf
4fd44ff26223a43e53b4142ddb3cad431012c8e62d9a9646fe8b2916d3099169  paper/A163665.tex
352f0e414ec085030a5268e96ddc1fd6f6728e02df6d1c7d567c3a885ddbab81  results/A163665_campaign_hits.tsv
af8cf84c1be09fa6beab60f9ec223cbd32e7a7b20797ab8bb048e324290791ca  src/Makefile
b72b98eb016f59ac6f44262e2cef7df86bdaf40a8274d1c63c457e0e45a3a21f  src/a163665.c
7948c44af7e6c692e180db04e0c5a9edb47d880f056f1489a80a1808d62a95ce  validation/run_manifest.tsv
c51b002eb9ce57ef47d804da9afb5338e4d6d2ea8002a382c5e53b98b751d5c9  validation/validation_notes.md
7c351924614a5300e2a520bbcd8afa6701f6fc1756f30557fd4a30d74cdd3d15  validation/verify_hits.gp
f883479a7ade6e139527a6417c1e5c88fd73925f3869ef2f62b8d7c1646ba7ea  validation/verify_prime_mapping.sh
```
