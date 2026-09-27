# Repository preparation, 27 September 2026

Source: the 26 September 2026 FINAL_SUBMISSION_PACKAGE. Original inputs,
statistical calculation, acceptance check, figure code and committed numerical
outputs were copied unchanged. Source and destination hashes are recorded in
SOURCE_FILES.json; the original package was not edited.

Two portability-only edits were made:

1. `build_supplement.py` reads the same reference records from
   `06_QC/REFERENCES_CURRENT.md`, removing its dependency on the full manuscript.
2. `build_workbook.mjs` defaults to the repository root (`.`) instead of a local
   submission-directory name; an explicit root argument remains supported.

Repository documentation, a single-thread entry point, numerical-only dependency
list, Git exclusions and checksum manifests were added. Main-text/document
builders and submission-package validators were excluded because their required
manuscript and author forms are not part of this code release.

No author identity, open license, new experiment, statistical result or DOI was
created during repository preparation.

## Supplied authorship metadata

On 27 September 2026 the user designated Manuscript(1).docx as the source for
author names, order and affiliations. AUTHORS.md and CITATION.cff transcribe
those records; the contribution statement is also present in that document.
No scientific code, reaction inputs, estimates or plots changed in this update.

## Public release metadata

On 27 September 2026 the project owner authorized public repository access,
DOI archiving, MIT licensing for the original code, and CC BY 4.0 for the
project's own derived data. License files, versioned citation metadata and
.zenodo.json were added. Third-party rights were retained. Scientific inputs,
analysis scripts, numerical outputs and figures were unchanged.
