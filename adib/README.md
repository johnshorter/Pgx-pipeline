# PGx Reporter — Standalone

Standalone pharmacogenomics reporting pipeline.

**Input:** a VCF file (GRCh38 coordinates).
**Output:** two reports — one for the patient (plain language, colorblind-safe
design) and one for the clinician (full technical detail, PharmCAT-backed). Both
are produced as HTML and, when WeasyPrint is installed, as PDF.

Internally the tool wraps the [PharmCAT](https://pharmcat.org) JAR (developed
by PharmGKB/Stanford) and extends it with a VCF preprocessing step that makes
it work cleanly on GIAB-style benchmark VCFs (which otherwise leave most PGx
genes uncallable).

---

## Features

- Single-command CLI: `python pgx_report.py sample.vcf`
- Colorblind-safe design system (symbols + shapes, not just colors):
  - `▲` Action · `◆` Review · `✓` Normal · `—` No Data
- Gene results grouped by functional category:
  - Phase I Metabolism (CYP Enzymes)
  - Phase II Metabolism
  - Drug Transporters
  - Immune Markers (HLA)
  - Other Pharmacogenes
- Two-pass VCF preprocessing — fills reference calls at missing PGx positions
  so that **all genotype-matchable genes** get called, not just the variants
  that happened to be in the VCF.
- CYP2D6 supported via PharmCAT research mode.
- Two report tiers sharing one theme (patient + clinician).

---

## Requirements

| Dependency | Minimum |
| --- | --- |
| Python | 3.10 |
| Java | 17 (for PharmCAT) |
| PharmCAT JAR | 3.2.0+ (download with the helper script below) |
| Python packages | `jinja2`, `weasyprint` (optional, for PDFs) |

On macOS / Linux:
```bash
# Java
brew install openjdk@17    # or: sudo apt install openjdk-17-jdk

# Python deps
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# PharmCAT JAR
python download_pharmcat.py
```

If you don’t install WeasyPrint, the pipeline still runs — it just emits HTML
only. Pass `--no-pdf` to skip PDF generation explicitly.

---

## Quick start

```bash
python pgx_report.py path/to/sample.vcf --output-dir ./output
```

The output directory will contain:
```
output/
├── cleaned/              # VCF with FORMAT/sample field mismatches fixed
├── preprocessed/         # VCF with reference calls filled for all PGx positions
├── pharmcat_raw/         # Output of the initial PharmCAT pass (position discovery)
├── pharmcat_final/       # Output of the final PharmCAT pass
└── reports/
    ├── patient_report.html       (+ .pdf)
    └── clinician_report.html     (+ .pdf)
```

---

## CLI reference

```
python pgx_report.py VCF [options]

Positional:
  VCF                          Input VCF file (GRCh38).

Options:
  -o, --output-dir DIR         Directory for all outputs (default: ./output).
      --jar PATH               Path to PharmCAT JAR (default: auto-detect in lib/).
      --no-pdf                 Skip PDF generation (HTML only).
      --skip-preprocess        Skip reference-fill preprocessing.
      --reference-phenotype P  Use an existing phenotype.json to discover PGx
                               positions (skips initial PharmCAT pass — useful
                               when batch-processing multiple samples).
      --research FLAGS         Comma-separated PharmCAT research flags
                               (default: "cyp2d6"; pass "" to disable).
      --timeout SEC            PharmCAT timeout (default: 300).
  -v, --verbose                Enable DEBUG logging.
```

---

## Pipeline stages

1. **Validate VCF** — format, chromosomes, reference build.
2. **Position discovery** — run PharmCAT once to enumerate all PGx positions
   it will check, or reuse a cached `phenotype.json` via
   `--reference-phenotype`.
3. **Screen + preprocess** — report coverage per gene, then write a new VCF
   with synthetic `0/0` reference calls filling the missing PGx positions
   (skipping genes that need specialized typing: HLA-A, HLA-B, MT-RNR1).
4. **Final PharmCAT run** — against the preprocessed VCF, with research mode
   enabled for CYP2D6 by default.
5. **Parse** — convert PharmCAT JSON into a structured model of gene results
   and drug recommendations.
6. **Render reports** — Jinja2 templates produce both report tiers.

Typical coverage on GIAB benchmark VCFs: **≈19–20 of 23 genes called** per
sample. Remaining uncallable genes (HLA-A, HLA-B, MT-RNR1) are documented as
known limitations in both reports.

---

## Batch processing

PGx positions are sample-independent — you can run the initial PharmCAT pass
once and reuse its `phenotype.json` for every subsequent sample:

```bash
# First sample: full two-pass run
python pgx_report.py HG001.vcf -o out/HG001

# Subsequent samples: reuse the cached position map
python pgx_report.py HG002.vcf -o out/HG002 \
    --reference-phenotype out/HG001/pharmcat_raw/HG001.phenotype.json
python pgx_report.py HG005.vcf -o out/HG005 \
    --reference-phenotype out/HG001/pharmcat_raw/HG001.phenotype.json
```

---

## Project layout

```
pgx-reporter-standalone/
├── pgx_report.py             # CLI entry point (this is what you run)
├── download_pharmcat.py      # Helper: download PharmCAT JAR from GitHub
├── requirements.txt
├── README.md
├── lib/                      # PharmCAT JAR goes here
└── src/
    ├── pharmcat/
    │   ├── vcf_validator.py       # VCF format & build checks
    │   ├── runner.py              # PharmCAT subprocess invocation
    │   └── output_parser.py       # PharmCAT JSON → structured model
    │                               # + gene categories + action symbols
    ├── preprocessing/
    │   └── preprocess_vcf.py      # clean_vcf + reference-fill
    ├── screening/
    │   ├── pgx_positions.py       # Parse phenotype.json for PGx positions
    │   └── screen_pharmacogenes.py  # Coverage screening report
    ├── reports/
    │   ├── patient_report.py      # Plain-language report generator
    │   └── clinician_report.py    # Technical report generator
    └── templates/
        ├── patient_report.html    # Jinja2 template (colorblind-safe)
        └── clinician_report.html  # Jinja2 template (colorblind-safe)
```

---

## Design system (reports)

Both reports share one visual theme with different levels of detail. Each
result is shown with a **shape**, a **label**, and a **color** — so it remains
readable for colorblind users and when printed in grayscale.

| Symbol | Meaning | Triggers |
| :---: | --- | --- |
| ▲ | **Action** | Poor metabolizer, ultrarapid metabolizer, HLA risk allele |
| ◆ | **Review** | Intermediate / rapid metabolizer, decreased or increased function |
| ✓ | **Normal** | Normal metabolizer, normal function |
| — | **No Data** | Gene not callable from the supplied VCF |

---

## Troubleshooting

- **`PharmCAT JAR not found`** → run `python download_pharmcat.py`.
- **`java: command not found`** → install Java 17+.
- **`Number of FORMAT entries does not match number of sample entries`**
  (common on GIAB VCFs) → the `clean_vcf` step handles this automatically in
  the default pipeline. Don’t use `--skip-preprocess` on such VCFs.
- **`HLA-A`, `HLA-B`, or `MT-RNR1` show as "No Data"** → expected.
  HLA requires dedicated HLA typing; MT-RNR1 is mitochondrial and absent from
  most whole-genome VCFs.
- **PDF not produced** → install WeasyPrint (`pip install weasyprint`) and its
  system dependencies, or just use the HTML output.

---

## License

PharmCAT is distributed under the Mozilla Public License 2.0 and is not
included in this repository — the download helper fetches it directly from
PharmGKB’s GitHub releases. All other code in this standalone package is
yours to use.
