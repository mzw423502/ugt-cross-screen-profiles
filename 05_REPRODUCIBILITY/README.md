# Scientific reproduction scripts

See the repository-root README for installation, commands, scope and licensing.
Keep the numbered directories together. `python reproduce.py` runs the numerical
analysis and its separate validator with one numerical thread. Add `--figures`
and/or `--supplement` to rebuild the corresponding outputs.

The matched reaction inputs and original numerical code are byte-identical to
the frozen submission package. Only the supplement reference-file path and the
optional workbook exporter's default root were adapted for this repository.
See `provenance/ADAPTATIONS.md` and `provenance/SOURCE_FILES.json`.
