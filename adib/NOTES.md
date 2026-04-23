# Adib's Pipeline

**What it does:** End-to-end pharmacogenomics reporting pipeline. Takes a VCF
file, runs PharmCAT 3.2.0 (with a custom two-pass preprocessing step that fills
reference calls at missing PGx positions so GIAB-style benchmark VCFs produce
usable results), and generates two reports: a plain-language patient report and
a technical clinician report.

**Input:** A VCF file (plain `.vcf` or gzipped `.vcf.gz`) using GRCh38
coordinates. Tested on GIAB benchmark samples HG001 (NA12878), HG002 (NA24385),
and HG005 — the pipeline handles the FORMAT/sample field mismatches common in
those files.

**Output:** Both reports are rendered as HTML and (when WeasyPrint is
installed) PDF, written to `output/reports/`:
- `patient_report.html` / `.pdf` — plain-language, colorblind-safe design
  (symbols `▲ ◆ ✓ —` + labels + color), genes grouped by functional category.
- `clinician_report.html` / `.pdf` — same theme, full technical detail
  (diplotypes, activity scores, star alleles, drug recommendations).

Intermediate artifacts are also kept: cleaned VCF, preprocessed VCF, and the
raw PharmCAT JSON output.

**Main script to run:**
```bash
python pgx_report.py <path-to-vcf> -o ./output
```
Example: `python pgx_report.py sample.vcf.gz -o ./output --no-pdf`

**Key dependencies:**
- Python 3.10+
- Java 17+ (required by PharmCAT) — a Temurin 17 JDK was bundled in `lib/jdk/`
  for local testing, but it's excluded from the repo to keep size down
- PharmCAT JAR 3.2.0+ — downloaded on demand via `download_pharmcat.py`
- `jinja2` (HTML templating)
- `weasyprint` (optional, for PDF output)

**Pipeline stages (6 in total):**
1. Validate the VCF.
2. Initial PharmCAT pass — discovers all PGx positions PharmCAT checks
   (or reuses a cached `phenotype.json` via `--reference-phenotype`, useful
   when batch-processing multiple samples).
3. Screen coverage + preprocess the VCF (fills `0/0` reference calls at
   missing PGx positions, skipping HLA-A/HLA-B/MT-RNR1).
4. Final PharmCAT pass on the preprocessed VCF (with `-research cyp2d6`
   enabled by default so CYP2D6 also gets called).
5. Parse PharmCAT JSON into a structured model.
6. Render both reports via Jinja2 templates.

**Notes for merging:**
- **Heavy artifacts excluded from the repo:** `lib/` (bundled JDK + PharmCAT
  JAR) and `output/` are not committed — they're too big for GitHub. The JAR
  can be fetched on demand with `python download_pharmcat.py`. A Java 17+
  runtime must be installed separately on the machine running the pipeline.
- **Two-pass preprocessing is the key trick.** GIAB benchmark VCFs only
  contain variant calls — reference-matching positions are absent, which
  causes PharmCAT to report "No Result" for 15/23 PGx genes. Preprocessing
  lifts this to ~19–20/23 genes. If Marco's pipeline is also based on
  GIAB-style data, we should keep this step.
- **Colorblind-safe design is load-bearing.** Both reports use shape + label
  + color (not color alone): `▲ Action · ◆ Review · ✓ Normal · — No Data`.
  Gene results are grouped into 5 functional categories (Phase I Metabolism /
  Phase II Metabolism / Drug Transporters / Immune Markers (HLA) / Other).
  This was a deliberate accessibility choice — any merged report should keep
  this system.
- **Known uncallable genes:** HLA-A and HLA-B need specialized HLA typing
  (not solvable from normal VCF); MT-RNR1 is mitochondrial and absent from
  most chr1-22 benchmark VCFs. These are documented as known limitations in
  both reports rather than silently missing.
- **CYP2D6 via PharmCAT research mode.** Passing `-research cyp2d6` lets
  PharmCAT attempt a call, but results carry a reliability caveat. Consider
  how Marco's pipeline handles CYP2D6 — we may want to align on one approach.
- **Repo layout (in `adib/`):**
  - `pgx_report.py` — CLI entry point
  - `download_pharmcat.py` — fetches the JAR from GitHub releases
  - `src/pharmcat/` — VCF validator, PharmCAT subprocess runner, output
    parser (with gene categories + action symbols)
  - `src/screening/` — coverage screening
  - `src/preprocessing/` — `clean_vcf` (fixes FORMAT mismatches) + reference-fill
  - `src/reports/` — patient + clinician report generators
  - `src/templates/` — Jinja2 HTML templates (shared design system)
  - `requirements.txt`, `README.md`
