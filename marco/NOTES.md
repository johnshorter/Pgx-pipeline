# NOTES

## What it does
End-to-end pharmacogenomics reporting pipeline: takes a patient VCF,
runs PharmCAT to call pharmacogene diplotypes/phenotypes, and produces
two audience-specific reports (one for the patient, one for the
clinician/researcher) in HTML and PDF.

## Input
- A single VCF file (`.vcf` or `.vcf.gz`), aligned to **GRCh38**.
- Chromosome names must use the `chr1, chr2, ..., chrX` convention
  (not `1, 2, ...`).
- Can be a small targeted VCF or a full whole-genome VCF (tested with
  GIAB HG005, ~2.3 GB).
- Optional CLI flags: `--sample-id`, `--output`, `--preprocess`,
  `--java-memory`, `--timeout`, `--keep-intermediate`, `--skip-validation`.

## Output
Per sample, under `output/<sample-id>/`:
- `patient_report.html` / `patient_report.pdf`
  - Lay-friendly: traffic-light gene overview cards (green / orange /
    red / grey) with color-blind-safe symbols, phenotype explanation,
    and grouped medication lists.
- `clinician_report.html` / `clinician_report.pdf`
  - Technical: diplotypes, star alleles, phenotypes, CPIC/DPWG
    recommendations, activity scores, and per-gene variant tables.
- `pharmcat/` subfolder with raw PharmCAT JSON (only kept with
  `--keep-intermediate`).

PDFs are generated via WeasyPrint; if WeasyPrint is not installed the
pipeline still produces the HTML versions.

## Main script to run
```
python generate_reports.py path/to/sample.vcf
```
Whole-genome VCF (recommended flags):
```
python generate_reports.py sample.vcf.gz --preprocess --java-memory 4g
```
Run `python generate_reports.py --help` for all options.

## Key dependencies
- **Python 3.10+**
- **Java 17+** (for PharmCAT; Eclipse Adoptium tested)
- **PharmCAT JAR** — must be placed at `lib/pharmcat.jar`
  (download from https://github.com/PharmGKB/PharmCAT/releases;
  not tracked in git)
- Python packages (see `requirements.txt`):
  - `jinja2` — HTML templating
  - `weasyprint` — HTML to PDF conversion (optional but recommended)
  - `pysam` / native `gzip` — VCF reading

## Notes for merging
- **No bcftools dependency.** The official PharmCAT preprocessor relies
  on bcftools, which is painful to install on Windows. Instead,
  `pharmcat_wrapper/preprocessor.py` is a pure-Python stream filter: it
  reads PharmCAT's built-in PGx position list straight out of the JAR
  (the `org/pharmgkb/pharmcat/definition/alleles/*_translation.json`
  files) and keeps only matching VCF lines. Use `--preprocess` on
  whole-genome VCFs — it cuts PharmCAT runtime from ~30–60 min down to
  seconds.
- **Off-by-one tolerance.** The preprocessor keeps positions
  `pos-1 / pos / pos+1` to be forgiving about indel-anchor differences
  between VCF conventions.
- **Windows-specific Java detection.** `config/settings.py` auto-finds
  Eclipse Adoptium JDKs under `C:\Program Files\Eclipse Adoptium`
  before falling back to `java` on PATH. On Linux/macOS it just uses
  `java` from PATH — no changes needed.
- **PharmCAT JAR is NOT in the repo.** Licensing is fine (MPL-2.0) but
  the JAR is large and versioned externally. Users download it.
- **Large VCFs and memory.** Default PharmCAT timeout is 3600s. For
  whole-genome VCFs without `--preprocess`, pass `--java-memory 4g` or
  `8g`. With `--preprocess`, defaults are usually fine.
- **Sample ID = output subfolder name.** If not passed, derived from
  the VCF basename (suffixes `.vcf` / `.vcf.gz` stripped).
- **Cleanup behavior.** By default the `pharmcat/` raw JSON folder is
  deleted after reports are generated. Pass `--keep-intermediate` to
  keep it (useful for debugging or re-generating reports without
  re-running PharmCAT).
- **Test fixtures.** `tests/fixtures/pharmcat_output/` contains real
  PharmCAT JSON output for **HG005 (GIAB public reference)** — safe to
  distribute. Do not replace these with files from actual patient VCFs
  unless you re-review for identifiability.
- **Known limitations.**
  - GRCh37 VCFs are not supported (PharmCAT requires GRCh38).
  - CYP2D6 copy-number / structural variants are not called from
    short-read VCFs; results for CYP2D6 are best-effort and flagged
    accordingly in the report.
  - Chromosome naming must be `chrN` style.
