# PGx Pipeline — Genetic Variant Interpretation for Personalized Medicine

A command-line pipeline that turns a single-sample genome VCF into **plain-language
pharmacogenomic reports**. It runs [PharmCAT](https://pharmcat.org/) to call pharmacogene
diplotypes and metabolizer phenotypes, maps them to **CPIC, DPWG, and FDA** drug-label
guidance, and renders two self-contained HTML reports:

- a **patient report** — plain-language summary of how your genes may affect medications;
- a **clinician/researcher report** — detailed diplotypes, phenotypes, activity scores, and
  drug-by-drug recommendations.

Every gene call is grouped into a colourblind-safe, four-level risk model so the most
important findings are easy to spot:

| Symbol | Level | Meaning |
|:------:|-------|---------|
| ▲ | **Action** | Atypical result that may warrant a change in drug or dose (e.g. poor/ultrarapid metabolizer). |
| ◆ | **Review** | Worth reviewing (e.g. intermediate/rapid metabolizer, decreased function). |
| ✓ | **Normal** | Typical function — standard prescribing usually applies. |
| — | **No Data** | No genotype could be called for this gene. |

> **Disclaimer.** This software is for **research and educational use only**. It is **not a
> medical device** and must **not** be used for diagnosis or treatment decisions. Reports are
> produced by computational analysis and must be confirmed by a CLIA/CLIA-equivalent certified
> laboratory before any clinical use. Always consult a qualified healthcare professional.

---

## Project context

Developed as a **bachelor project at Roskilde University (RUC), Spring 2026** —
*"Building a Genetic Variant Interpretation Pipeline for Personalized Medicine"* —
by **Adib Vishtal Ahmad** and **Marco Antonio Román Faraldo**, supervised by **John Shorter**.

The pipeline began as two independent reporters (`adib/` and `marco/`) which were later merged
into the unified codebase under [`src/`](src/). The `adib/` and `marco/` folders are kept as
**legacy reference implementations only** — use `src/pgx_report.py` for all real runs.

---

## What it reports

The pipeline reports on the pharmacogenes PharmCAT supports, organised into functional
categories:

- **Phase I metabolism (CYP enzymes):** CYP2B6, CYP2C9, CYP2C19, CYP2D6, CYP3A4, CYP3A5, CYP4F2
- **Phase II metabolism:** DPYD, NAT2, TPMT, NUDT15, UGT1A1
- **Drug transporters:** ABCG2, SLCO1B1
- **Immune markers (HLA):** HLA-A, HLA-B
- **Other pharmacogenes:** CACNA1S, CFTR, F2, F5, G6PD, IFNL3, MT-RNR1, RYR1, VKORC1

Recommendations are sourced from **CPIC**, **DPWG**, and **FDA** drug labels via PharmCAT's
reporter.

---

## See an example report (no install needed)

The repository ships pre-generated example output you can open in a browser. Because GitHub
does not render raw HTML files inline, use one of these options:

**Option A — view online via htmlpreview:**

- Patient report (HG005): <https://htmlpreview.github.io/?https://github.com/johnshorter/Pgx-pipeline/blob/main/output/HG005/patient_report_v2.html>
- Clinician report (HG005): <https://htmlpreview.github.io/?https://github.com/johnshorter/Pgx-pipeline/blob/main/output/HG005/clinician_report_v2.html>

**Option B — clone and open locally:**

- [`output/HG005/patient_report_v2.html`](output/HG005/patient_report_v2.html) — a real public
  benchmark genome (Genome in a Bottle **HG005**), showing a mix of Action / Review / Normal calls.
- [`output/HG005/clinician_report_v2.html`](output/HG005/clinician_report_v2.html) — the matching detailed report.
- [`output/reference_test/`](output/reference_test/) — an all-reference baseline (everything Normal),
  equivalent to what the [Quick start](#quick-start) below produces.

> HG005 is a **public reference genome**, not personal data.

---

## Prerequisites

| Requirement | Version | Check |
|-------------|---------|-------|
| **Python**  | 3.10 or newer | `python --version` |
| **Java (JDK)** | 17 or newer (PharmCAT 3.2.0 needs JDK 17+) | `java -version` |
| **PharmCAT JAR** | 3.2.0+ | downloaded in the install step below |

Works on **Windows, macOS, and Linux**. Java does **not** need to be on your `PATH` — the
pipeline auto-detects a JDK from a bundled `lib/jdk/`, from Eclipse Adoptium under
`C:\Program Files\Eclipse Adoptium\` on Windows, or from `java` on `PATH`. You can also point
at a specific binary with `--java-executable`.

---

## Installation

```bash
# 1. Clone
git clone https://github.com/johnshorter/Pgx-pipeline.git
cd Pgx-pipeline

# 2. Create and activate a virtual environment
python -m venv .venv
#   macOS / Linux:
source .venv/bin/activate
#   Windows (PowerShell):
.venv\Scripts\Activate.ps1

# 3. Install Python dependencies
pip install -r requirements.txt

# 4. Download the PharmCAT JAR into lib/ (~32 MB)
python src/download_pharmcat.py --version 3.2.0
```

**Note on PDFs.** Reports are always written as **HTML**. PDF export is *best-effort*: it relies
on WeasyPrint, which needs GTK system libraries that are often missing on Windows. If PDF
generation fails you will see a WeasyPrint warning and a `*.pdf_not_generated.txt` placeholder —
the HTML reports are unaffected. Pass `--no-pdf` to skip PDF entirely and silence the warning.

---

## Quick start

A tiny example VCF is included so you can produce a full report immediately:

```bash
python src/pgx_report.py examples/reference_baseline.vcf --ref-fill -o output
```

This writes reports to `output/reference_baseline/`:

```
output/reference_baseline/
├── patient_report.html        ├── patient_report_v2.html      (risk-grouped)
└── clinician_report.html      └── clinician_report_v2.html    (risk-grouped)
```

Open any `*.html` file in a web browser. The `_v2` reports are the current, risk-grouped layout
and are the recommended ones to read.

> [`examples/reference_baseline.vcf`](examples/reference_baseline.vcf) is an empty GRCh38 VCF;
> `--ref-fill` fills every pharmacogene position with the reference genotype, producing an
> "everything Normal" baseline. To see a report with real variation, run the pipeline on a
> genome VCF of your own (see below) or view the bundled HG005 example above.

---

## Usage

```
python src/pgx_report.py <VCF> [options]
```

The input is a **single-sample VCF aligned to GRCh38** (`.vcf` or `.vcf.gz`). The two
preprocessing flags handle the common input types:

| Your input | Recommended command |
|------------|---------------------|
| **Whole-genome / whole-exome VCF** (millions of variants) | `--filter` (subset to PharmCAT positions first — much faster) |
| **GIAB benchmark VCF** (e.g. HG005) or any variant-only VCF where absence means reference | `--ref-fill` |
| **Large WGS, and you also want missing positions treated as reference** | `--filter --ref-fill` |
| **Small targeted PGx VCF already covering the positions** | *(no preprocessing flags)* |

### Examples

```bash
# Whole-genome VCF: filter to PharmCAT positions, give Java more heap
python src/pgx_report.py wgs.vcf.gz --filter --java-memory 4g -o output

# GIAB benchmark VCF: reference-fill missing positions
python src/pgx_report.py HG005.vcf.gz --ref-fill -o output

# Both (filter first, then fill the small filtered VCF)
python src/pgx_report.py wgs.vcf.gz --filter --ref-fill -o output

# Batch: reuse one phenotype.json to skip per-sample position discovery
python src/pgx_report.py HG002.vcf.gz --ref-fill \
    --reference-phenotype output/HG001/pharmcat_raw/HG001.phenotype.json -o output
```

### All options

| Option | Description |
|--------|-------------|
| `vcf` (positional) | Input VCF file (GRCh38), `.vcf` or `.vcf.gz`. |
| `-o, --output-dir DIR` | Output root directory (default: `./output`). |
| `--sample-id ID` | Name of the per-sample output subfolder (default: VCF filename without extension). |
| `--jar PATH` | Path to the PharmCAT JAR (default: auto-detect in `lib/`). |
| `--filter` | Pre-filter the VCF to PharmCAT-relevant positions (WGS fast-path). |
| `--ref-fill` | Fill missing PGx positions with reference `0/0` calls (GIAB rescue). |
| `--skip-validation` | Skip the VCF validation step. |
| `--reference-phenotype PATH` | Existing `phenotype.json` for position discovery; with `--ref-fill`, skips the initial PharmCAT pass. |
| `--research FLAGS` | Comma-separated PharmCAT research-mode flags (default: none). See the CYP2D6 caveat below. |
| `--java-memory SIZE` | Java `-Xmx` heap size, e.g. `4g` (recommended for large VCFs without `--filter`). |
| `--java-executable PATH` | Java binary to use (default: auto-detected). |
| `--timeout SECONDS` | PharmCAT subprocess timeout (default: 3600). |
| `--no-pdf` | Skip PDF generation (HTML only). |
| `--keep-intermediate` | Keep PharmCAT raw output and intermediate VCFs (default: deleted after reports). |
| `-v, --verbose` | Enable DEBUG logging. |

---

## Input requirements

The VCF is validated before processing (skip with `--skip-validation`):

- **Format:** must be a valid VCF with `##fileformat=` and a `#CHROM` header line.
- **Single sample:** one sample column is expected.
- **Genome build: GRCh38 / hg38.** GRCh37/hg19 inputs are **rejected** — lift them over to
  GRCh38 first (e.g. with Picard `LiftoverVcf` or `bcftools +liftover`).
- **Chromosome naming:** chr-prefixed contigs (`chr1`, `chr2`, …) are expected. Numeric contigs
  (`1`, `2`, …) can cause PharmCAT to return all no-calls; rename them to `chr`-prefixed first.
- **Size:** up to 3000 MB (sized for whole-genome inputs).

---

## Output

Each run creates `output/<sample-id>/` containing:

| File | Audience |
|------|----------|
| `patient_report_v2.html` | **Patient** — risk-grouped, plain language *(recommended)* |
| `clinician_report_v2.html` | **Clinician/researcher** — risk-grouped, detailed *(recommended)* |
| `patient_report.html`, `clinician_report.html` | Original (v1) layouts, kept for comparison |
| `*.pdf` *(best-effort)* | PDF copies if WeasyPrint is available; otherwise a `*.pdf_not_generated.txt` note |

With `--keep-intermediate`, the run also keeps `pharmcat/` (final PharmCAT JSON output),
`pharmcat_raw/` (the position-discovery pass used by `--ref-fill`), and `intermediates/`
(filtered / reference-filled VCFs).

---

## How it works

```
VCF ─► validate ─► [filter] ─► [position discovery] ─► [reference-fill] ─► PharmCAT ─► parse ─► render
                   (--filter)   (--ref-fill, pass 1)     (--ref-fill, pass 2)                    HTML reports
```

1. **Validate** the VCF (format, single sample, GRCh38 build).
2. **Filter** *(optional, `--filter`)* — subset a whole-genome VCF to PharmCAT positions.
3. **Position discovery** *(when `--ref-fill` is set)* — an initial PharmCAT pass enumerates
   every position PharmCAT checks.
4. **Reference-fill** *(optional, `--ref-fill`)* — positions absent from the VCF are written as
   reference `0/0` so PharmCAT can assign diplotypes. HLA-A, HLA-B, and MT-RNR1 are never filled.
5. **Run PharmCAT** to call diplotypes, phenotypes, and drug recommendations.
6. **Parse** PharmCAT's JSON into the four-level risk model.
7. **Render** the patient and clinician reports.

---

## Known limitations & troubleshooting

| Symptom | Cause / fix |
|---------|-------------|
| **HLA-A, HLA-B, MT-RNR1 always show "No Data"** | These need specialised typing (HLA imputation, mitochondrial calling) that PharmCAT does not perform from a standard VCF. Expected behaviour. |
| **CYP2D6 missing or flagged** | CYP2D6 is unreliable from short-read whole-genome sequencing due to complex structural variation. Use dedicated clinical CYP2D6 testing for actionable decisions. `--research cyp2d6` enables best-effort calling **but PharmCAT 3.2.0+ then disables the full reporter, dropping all drug recommendations** — avoid unless you only need the raw call. |
| **`VCF appears to use genome build GRCh37`** | The input is GRCh37/hg19. Lift it over to GRCh38 before running. |
| **PharmCAT returns all no-calls** | Often numeric contigs (`1` instead of `chr1`). Rename contigs to `chr`-prefixed GRCh38 names. |
| **A whole-genome run "looks Normal" everywhere with `--ref-fill`** | `--ref-fill` assumes any absent position equals reference. This is valid for whole-genome / GIAB inputs but **not** for variant-only panels with limited coverage, where it can produce false "Normal" calls. |
| **WeasyPrint / `libgobject` errors, no PDF** | PDF export is best-effort and often fails on Windows. HTML reports are still produced; add `--no-pdf` to silence the warning. |
| **`No JAR asset found` / no PharmCAT JAR** | Run `python src/download_pharmcat.py --version 3.2.0`, or pass `--jar path/to/pharmcat.jar`. |
| **Run is very slow on a whole-genome VCF** | Add `--filter` (and optionally `--java-memory 4g`). |

---

## Repository layout

```
Pgx-pipeline/
├── src/                       # The unified pipeline (use this)
│   ├── pgx_report.py          #   CLI entry point
│   ├── download_pharmcat.py   #   PharmCAT JAR downloader
│   ├── config/settings.py     #   Paths, gene categories, risk model, disclaimers
│   ├── pharmcat/              #   Runner, VCF validator, output parser
│   ├── preprocessing/        #   clean / filter / reference-fill VCF steps
│   ├── reports/              #   Report generators + Jinja2 HTML templates
│   └── screening/            #   PGx position handling & coverage screening
├── examples/                  # Small example input VCF for the quick start
├── output/                    # Pre-generated example reports (HG005, reference_test)
├── docs/                      # Report figures
├── adib/, marco/              # Legacy standalone reporters (reference only)
├── lib/                       # PharmCAT JAR (downloaded; not in git)
└── requirements.txt
```

---

## Acknowledgements

This pipeline is a front-end to **[PharmCAT](https://pharmcat.org/)** (Pharmacogenomics Clinical
Annotation Tool), developed by PharmGKB / Stanford University. Pharmacogenomic guidance comes from
**[CPIC](https://cpicpgx.org/)**, **[DPWG](https://www.knmp.nl/)**, and FDA drug labels.
