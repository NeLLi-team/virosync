# Run frameshift-sensitive VS marker screening

Frameshift screening searches the masked assembly with the shipped
`VS######` marker profiles. It can recover domains with a frameshift or
in-frame stop codon that protein prediction can miss. ViroSync validates each
domain against the viral marker database before it can seed a candidate
region. The normal boundary and acceptance filters still apply.

## Install BATH

Install a C compiler and `make` on the host before running `setup-bath`.

```bash
pixi run setup-bath
```

## Configure the screen

The screen is off by default. Enable it for one run:

```bash
pixi run virosync \
  -i assembly.fna \
  -o results/frameshift \
  --config config/orchestration.yaml \
  --frameshift-screening
```

To enable it in a config file, set:

```yaml
phase1:
  frameshift_screening_enabled: true
```

## Run the shipped example

The example uses three *Trichomonas vaginalis* G3 contigs from
[GCA_000002825.3](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000002825.3/).

```bash
pixi run example-frameshift
cat results/example-frameshift/batch_summary.tsv
```

With resource bundle v1.0.7, the row must report `status=success`,
`predictions=5`, and `accepted=2`. This example checks the rescue path. It is
not a sensitivity benchmark.

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

A confirmed rescued domain supports a degraded viral marker locus. It is not a
repaired full-length protein or coding sequence. The BATH cutoff is fixed and
has not been calibrated for each family. The screen uses only the shipped VS
profiles.
