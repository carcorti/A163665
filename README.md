# Computational Extension of OEIS Sequence A163665

Author: Carlo Corti.

This repository contains the C17 segmented search code, data tables, selected validation artifacts, and the pre-release manuscript supporting a computational extension of OEIS A163665.

OEIS A163665 is the sequence of integers \(k\) such that

\[
k = \operatorname{prime}(\varphi(\sigma(\pi(k)))).
\]

The computation uses the exact reduction to prime indices: a prime-index \(n\) is a hit precisely when

\[
\varphi(\sigma(n)) = n,
\]

and the corresponding OEIS term is \(k=\operatorname{prime}(n)\).

## Main result

The production campaign searched the contiguous prime-index interval

\[
234954224 \le n \le 405234954223.
\]

It found eight new hit indices:

| New OEIS index | Hit index \(n\) | Term \(k=\operatorname{prime}(n)\) |
|---:|---:|---:|
| 21 | 2147483648 | 50685770167 |
| 22 | 3889036800 | 94206075193 |
| 23 | 4389396480 | 106882008719 |
| 24 | 21946982400 | 571264143487 |
| 25 | 47416320000 | 1272242879459 |
| 26 | 92177326080 | 2536961181761 |
| 27 | 133145026560 | 3715374607207 |
| 28 | 331914240000 | 9576704442589 |

Subject to the correctness of the validated C17 production kernel and to the integrity of the recorded official run manifests, the computation certifies that there are no further terms through

\[
\operatorname{prime}(405234954223)=11776117514717.
\]

The consolidated candidate b-file is:

- `data/b163665.txt`

The certified-term table for the newly reported terms is:

- `data/certified_terms.tsv`

## Repository structure

```text
.
├── CITATION.cff
├── data
│   ├── b163665.txt
│   └── certified_terms.tsv
├── docs
│   ├── confirmation_A163665.md
│   └── strategy_A163665.md
├── .gitignore
├── LICENSE
├── README.md
├── paper
│   └── A163665_v6.tex
├── results
│   └── A163665_campaign_hits.tsv
├── src
│   ├── a163665.c
│   └── Makefile
└── validation
    ├── run_manifest.tsv
    ├── validation_notes.md
    ├── validation_summary.md
    ├── verify_hits.gp
    └── verify_prime_mapping.sh
```

The manuscript filename under `paper/` is intentionally pre-release at this stage. It should be finalized only after the first GitHub upload, Zenodo activation, GitHub release, and Zenodo archive creation.

## Method

The search is performed over prime indices \(n\), not directly over candidate terms \(k\). The production code applies:

- the exact prime-index reduction \(k=\operatorname{prime}(n)\);
- the parity filter for \(n>1\);
- multiplicative computation of \(\sigma(n)\);
- factorization of \(\sigma(n)\);
- the exact rejection condition \((q-1)\mid n\) for every prime divisor \(q\) of \(\sigma(n)\);
- exact computation of \(\varphi(\sigma(n))\);
- conversion from hit index \(n\) to sequence term \(k=\operatorname{prime}(n)\).

The C source is in:

- `src/a163665.c`

The build file is:

- `src/Makefile`

## Validation

The public package records the final validation state in:

- `validation/validation_summary.md`

The official run-level manifest is:

- `validation/run_manifest.tsv`

The reported aggregate validation state is:

| Quantity | Value |
|---|---:|
| Official runs | 15 |
| Certified subsegments | 40500 |
| Raw \(n\)-width | 405000000000 |
| Candidates tested | 202500000000 |
| Unresolved segments | 0 |
| Failed segments | 0 |
| New hits | 8 |
| Reconstructed elapsed total | 22h:12m:58s |

The positive rows are independently checked by:

- `validation/verify_hits.gp`, which recomputes \(\varphi(\sigma(n))=n\) and checks primality of the corresponding \(k\)-values;
- `validation/verify_prime_mapping.sh`, which verifies \(k=\operatorname{prime}(n)\) with `primecount`.

These scripts are independent positive-hit checks. They do not constitute a second exhaustive scan of the non-hit interval. The negative certification rests on the validated segmented C17 campaign and its recorded manifests.

## Reproducibility

The package validation environment recorded in the manuscript was:

- Linux Mint 22.3;
- AMD Ryzen 9 7940HS, 8 cores, 16 hardware threads;
- 64 GB DDR5 installed RAM;
- GCC 13.3.0;
- ISO/IEC 9899:2018 C / C17;
- OpenMP;
- PARI/GP 2.15.4;
- primecount 7.10.

The manuscript records the following validation commands as passed:

```sh
make -C src clean all
make -C src test THREADS=4
make -C src sanitize
gp -q validation/verify_hits.gp
validation/verify_prime_mapping.sh
make -C src clean
```

The large local subsegment directories generated during the campaign are not included in this public package. Their aggregate certified content is represented by the run manifest and validation summary.

## Paper

The accompanying pre-release manuscript is:

- `paper/A163665_v6.tex`

Current title:

> Computational Extension of OEIS Sequence A163665: Eight New Terms via a Segmented Prime-Index Search

The public paper filename is intentionally not finalized yet. After the first GitHub upload, Zenodo activation, GitHub release, and Zenodo archive creation, the repository metadata and manuscript references should be updated with the final release tag and Zenodo DOI.

## Data and result files

- `data/b163665.txt` — candidate extended OEIS b-file.
- `data/certified_terms.tsv` — certified new terms with witnesses.
- `results/A163665_campaign_hits.tsv` — consolidated campaign hit table.
- `validation/run_manifest.tsv` — official run-level manifest.
- `validation/validation_summary.md` — final validation summary.
- `validation/validation_notes.md` — validation notes.
- `docs/strategy_A163665.md` — computational strategy notes.
- `docs/confirmation_A163665.md` — confirmation notes.

## GitHub and Zenodo status

Repository URL intended for publication:

- `https://github.com/carcorti/A163665`

Zenodo DOI placeholder:

- `10.5281/zenodo.xxxxxxxx`

The Zenodo DOI is intentionally still a placeholder in this pre-release package. After Zenodo archives the first GitHub release, replace the placeholder DOI in the paper, `README.md`, and `CITATION.cff`.

## Citation

Please cite the repository using the metadata in:

- `CITATION.cff`

Until the Zenodo archive has been created, the DOI entry remains a placeholder and must not be treated as an active DOI.

## License

This project is released under the MIT License. See:

- `LICENSE`
