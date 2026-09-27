# Arabidopsis UGT cross-screen profiles

Analysis and figure code for **Cross-screen substrate profiles of Arabidopsis
UDP-glycosyltransferases diverge by chemical series**.

This repository reproduces the processed, matched comparison of plant UGT
screens reported by [Yang et al. (2018)](https://doi.org/10.1038/s41589-018-0154-9)
and [Sirirungruang et al. (2025)](https://doi.org/10.1038/s41467-025-61530-6).
It contains 40 enzyme constructs, 15 acceptors and 600 reaction units. The 22
Y18 `NOT_TESTED` units remain explicit; 578 tested pairs enter binary comparisons.

## Reproduce the numerical results

Use Python 3.12. Install the numerical requirements, then run the entry point
from the repository root:

```sh
python -m venv .venv
# Activate .venv using the command appropriate for your shell.
python -m pip install -r 05_REPRODUCIBILITY/requirements-numerics.txt
python reproduce.py
```

The entry point sets all numerical thread limits to one and runs calculation
followed by a separate acceptance check. It replaces files in
`05_REPRODUCIBILITY/generated/`; use a fresh clone to retain the committed snapshot.
No GPU, large model, raw-MS download or account credentials are required.

The calculation reads only `S1_EnzymeIdentity.csv`, `S2_AcceptorIdentity.csv`,
`S3_ReactionUnits600.csv` and `recompute_config.json`. Frozen expected values
are read exclusively by the separate validator. Validation uses tolerance
`1e-10`, checks six effect estimates and intervals, and retains all 1,000
bootstrap rows (500 per analysis set; seed 20260919).

## Rebuild figures and supplementary tables

```sh
python -m pip install -r 05_REPRODUCIBILITY/requirements-analysis.txt
python reproduce.py --figures --supplement
```

Main figures are written to `02_FIGURES/`; supplementary figures and readable
tables go to `03_SUPPLEMENT/`. Main figures use a 180-mm canvas and are intended
for full-width placement of at least 169 mm. The build exports vector PDF/SVG,
300-dpi PNG and TIFF. TIFFs are generated locally and are not tracked by Git.
Figure 3 uses explicit scaffold/substitution drawings checked against recorded
SMILES/connectivity keys; RDKit is not required.

`--supplement` also regenerates the 24-sheet workbook data as
`05_REPRODUCIBILITY/generated/workbook_content.json`. A complete XLSX is included.
The optional original formatting exporter `build_workbook.mjs` requires the
environment-specific `@oai/artifact-tool` 2.8.74 package; it is not required for
numerical reproduction and is not represented as a public npm dependency.
Where that package is available, run:

```sh
node 05_REPRODUCIBILITY/build_workbook.mjs . internal_checks/workbook_previews
```

## Files and evidence scope

| Location | Contents |
| --- | --- |
| `05_REPRODUCIBILITY/input/` | Processed reaction units, identity records, sensitivity inputs and source provenance |
| `05_REPRODUCIBILITY/generated/` | Numerical summaries, full resamples, and workbook data |
| `05_REPRODUCIBILITY/*.py` | Statistics, verification, plotting and supplementary-table code |
| `02_FIGURES/` | Current Figures 1–3 |
| `03_SUPPLEMENT/` | Readable tables, source-data workbook, captions and supplementary figures |
| `06_QC/` | Current bibliography and reference verification |
| `provenance/` | Source-file hashes and repository adaptation record |

`TESTED_ACTIVE`, `TESTED_INACTIVE`, `NOT_TESTED` and `AMBIGUOUS` remain distinct.
Missing or ambiguous units never become negatives. Multiple product rows do not
constitute independent enzyme–acceptor reactions. Computation starts from
processed source calls and recorded identity matches; it does not reconstruct
raw mass spectra or sequence alignments. Historical inputs remain provenance
records; current numerical results are in `generated/` and current references
are in `06_QC/REFERENCES_CURRENT.md`.

The HCA series contracts from 38/117 to 1/117 positive pairs. Coumarins retain
87/116 positives in each screen, with 65 shared positives and 22 unique to each
screen. These are cross-screen observations, not isolated causal effects of
pH, pooling or other jointly changing experimental conditions.

Third-party article PDFs, raw mass spectra, manuscript drafts, author forms,
credentials, and historical algorithm benchmarks are outside this repository.
The input source manifest records original public records and available hashes.

## Provenance, citation and rights

The initial repository snapshot derives from the locally verified manuscript
package dated 26 September 2026. Numerical and input files are preserved;
repository-only portability edits are listed in `provenance/ADAPTATIONS.md`.
`SHA256SUMS.txt` records the committed payload bytes. Local validation results
are recorded in `provenance/VALIDATION.json`.

Cite the two original studies for their experimental data. For this software,
record the repository URL and exact Git commit. The public archive does not imply manuscript acceptance. Author names, order and affiliations
were supplied from Manuscript(1).docx; see `AUTHORS.md` and `CITATION.cff`. See `RIGHTS.md`
for the current licensing status.

## Public release and licensing

Version 1.0.0 preserves the scientific code, reaction inputs and numerical
results of the verified snapshot. Original code and software documentation
are licensed under MIT; original derived data, annotations and compilation
are licensed under CC BY 4.0. Third-party source material retains its original
rights. See LICENSE, LICENSE-DATA.md and RIGHTS.md.

Version 1.0.0 is archived at [10.5281/zenodo.22994959](https://doi.org/10.5281/zenodo.22994959).
The record contains the five manuscript authors and three affiliations.
The archive preserves release commit `adff639bea9b2ff64660b80c65d4b4379b2ca472`.
Cite this version DOI when referring to the results reproduced by this release.
