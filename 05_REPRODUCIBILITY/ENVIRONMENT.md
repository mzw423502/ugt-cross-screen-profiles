# Recorded environment

Numerical/plotting baseline: Windows 11, Python 3.12.10, NumPy 2.5.1,
pandas 3.0.5, Matplotlib 3.11.1, Pillow 12.3.0 and PyMuPDF 1.28.2.
The pinned requirements describe the environment used for the frozen output.

All entry-point subprocesses use OMP_NUM_THREADS, MKL_NUM_THREADS,
OPENBLAS_NUM_THREADS and NUMEXPR_NUM_THREADS set to 1. CPU only, one process
at a time; no model training, docking, MD or raw-MS reanalysis.

Arial is preferred for plots with DejaVu Sans as fallback. PDF fonts are
embedded or converted to vector paths. Exact rendering may vary by font/runtime;
numerical reproduction is checked separately at tolerance 1e-10.

Optional workbook formatting used Node.js with @oai/artifact-tool 2.8.74 from
the Codex runtime. That environment-specific exporter is not needed for the
analysis; CSV and JSON versions expose every data cell independently.
