# Validation Notes

This directory contains the archival validation summary for the official A163665 campaign.

The public package includes run-level telemetry in `run_manifest.tsv` and the final validation report in `validation_summary.md`. The many local subsegment files under `runs/campaign/A163665_run_*.tmp/` are intentionally not included in this raw GitHub package.

Certification semantics:

- every included official run has `coverage=CONTIGUOUS`;
- every included official run has `segments_unresolved=0`;
- every included official run has `segments_failed=0`;
- only `A163665_run_001` through `A163665_run_015` are part of the official synthesis;
- `A163665_run_000` was a pilot overlap run and is excluded.

The final package validation repeated:

```text
selftest OK
one-point reproduction of the largest discovered hit
global run-summary contiguity check
official submanifest count check
PARI/GP positive-hit verification
primecount prime-index mapping verification
```

See `validation_summary.md` for the numerical details and SHA-256 hashes.
