# Example input

[`reference_baseline.vcf`](reference_baseline.vcf) is a minimal, header-only GRCh38 VCF (no
variant rows). On its own it carries no genotypes; run it with `--ref-fill` and the pipeline
fills every pharmacogene position with the reference genotype, producing an "everything Normal"
baseline report. It exists so the pipeline can be run end-to-end with a single command:

```bash
python src/pgx_report.py examples/reference_baseline.vcf --ref-fill -o output
```

Reports are written to `output/reference_baseline/`. This is the same baseline shown in
[`output/reference_test/`](../output/reference_test/).

To see a report with real variation, run the pipeline on a genome VCF of your own, or view the
bundled [`output/HG005/`](../output/HG005/) example (a public Genome-in-a-Bottle benchmark genome).
See the [top-level README](../README.md) for full usage.
