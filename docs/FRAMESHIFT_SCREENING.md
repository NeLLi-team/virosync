# Screen for frameshifted markers

Frameshift screening can recover viral marker domains with a frameshift or
in-frame stop codon that protein prediction misses. It requires BATH and is
off by default.

## Run the example

After [installing ViroSync](getting-started.md), build BATH and run the
three-contig *Trichomonas vaginalis* G3 example. A C compiler and `make` must
be available on the host.

```bash
pixi run setup-bath
pixi run example-frameshift
cat results/example-frameshift/batch_summary.tsv
```

The contigs come from
[GCA_000002825.3](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000002825.3/).
With resource bundle v1.0.7, the `trichomonas-g3` row reports
`status=success`, `predictions=5`, and `accepted=2`. This example checks marker
rescue; it does not measure sensitivity.

## Screen your assembly

Add `--frameshift-screening` to a ViroSync run, or enable it in the config:

```yaml
phase1:
  frameshift_screening_enabled: true
```

## Read the output

Each genome writes frameshift files under
`phase1/frameshift_screening/`:

- `confirmed_frameshift_markers.tsv` lists validated domains that can seed
  regions.
- `confirmed_frameshift_proteins.faa` contains their amino-acid domain
  sequences.
- `validation/validated_marker_hits.tsv` records the evidence and validation
  result for each candidate domain.

Accepted calls are in the standard `virosync_predictions.tsv` and GFF3
outputs. An accepted per-EVE FAA can also contain a confirmed rescued domain.

## Interpret the result

The screen searches the masked assembly with the shipped `VS######` profiles.
ViroSync validates each domain against the viral marker database before it can
seed a candidate region. The normal boundary and acceptance filters apply.

A confirmed domain supports a degraded viral marker locus. It does not
reconstruct a full-length protein or coding sequence. The BATH cutoff is fixed
and has not been calibrated for each family.
