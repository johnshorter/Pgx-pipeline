# PGx Reporter

End-to-end pharmacogenomics reporting pipeline. Takes a VCF file and
produces both a **patient-friendly** and a **clinician/researcher** report
(HTML + PDF) using [PharmCAT](https://github.com/PharmGKB/PharmCAT) and
CPIC/DPWG guidelines.

Developed as a bachelor thesis project.

## Pipeline

```
VCF -> validation -> (optional preprocess) -> PharmCAT -> patient + clinician reports
```

See `docs/pipeline_flowchart.png` for a visual overview.

## Requirements

- Python 3.10+
- Java 17+ (tested with Eclipse Adoptium)
- [PharmCAT JAR](https://github.com/PharmGKB/PharmCAT/releases) placed at
  `lib/pharmcat.jar` (not tracked in git)

Install Python dependencies:

```
python -m venv venv
venv\Scripts\activate     # Windows
pip install -r requirements.txt
```

## Usage

```
python generate_reports.py path/to/sample.vcf
python generate_reports.py sample.vcf.gz --output my_results
python generate_reports.py huge_wgs.vcf.gz --preprocess --java-memory 4g
```

### Useful flags

| Flag | Purpose |
|------|---------|
| `--output DIR` | Output root directory (default `./output`) |
| `--sample-id ID` | Subfolder / sample name (default: VCF basename) |
| `--preprocess` | Pre-filter VCF to PharmCAT positions (**strongly recommended for whole-genome VCFs**) |
| `--java-memory 4g` | Java `-Xmx` heap size |
| `--timeout 3600` | Max seconds to wait for PharmCAT |
| `--keep-intermediate` | Keep PharmCAT raw JSON output |
| `--skip-validation` | Skip the VCF validation step |

Run `python generate_reports.py --help` for the full list.

## Output layout

```
output/
  <sample-id>/
    patient_report.html
    patient_report.pdf
    clinician_report.html
    clinician_report.pdf
```

## Privacy

**Do not commit patient VCFs or generated reports** -- both contain
identifiable genetic data. The included `.gitignore` blocks `data/`,
`output/`, and `*.vcf*` files by default.

## Disclaimer

This project is for research and educational purposes only. It is not a
medical device and must not be used for clinical decision-making without
confirmation by a CLIA-certified laboratory and a licensed clinician.
