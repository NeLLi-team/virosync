# Run frameshift-sensitive VS marker screening

Frameshift screening searches the masked assembly with the shipped
`VS######` marker profiles. It can recover domains with a frameshift or
in-frame stop codon that protein prediction can miss. ViroSync validates each
domain against the viral marker database before it can seed a candidate
region. The normal boundary and acceptance filters still apply.

The screen runs inside ViroSync with the standard Pixi environment.

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
`predictions=3`, and `accepted=2`. This example checks the rescue path. It is
not a sensitivity benchmark.

## Read the output

Each genome writes frameshift files under
`phase1/frameshift_screening/`:

- `confirmed_frameshift_markers.tsv` lists validated domains that can seed
  regions.
- `confirmed_frameshift_proteins.faa` contains their amino-acid domain
  sequences.
- `frameshift_hits.tsv` lists candidate domains and their alignment evidence.
- `frameshift_events.tsv` records candidate event positions as zero-based,
  half-open genomic intervals. An interval marks the affected codon span;
  it is not an experimentally determined mutation site.
- `validation/validated_marker_hits.tsv` records the evidence and validation
  result for each candidate domain.

Accepted calls are in the standard `virosync_predictions.tsv` and GFF3
outputs. An accepted per-EVE FAA can also contain a confirmed rescued domain.

## Interpret the result

A confirmed rescued domain supports a disrupted viral marker locus. Its
sequence is an aligned domain, not a repaired full-length protein or coding
sequence. `X` marks an uncertain residue at a frameshift, stop codon, or
ambiguous DNA. Event annotations cannot distinguish biological decay,
programmed frameshifting, introns, and assembly errors.

The screen uses only the shipped VS profiles. Its six-frame seed search can
miss loci whose conserved fragments are too weak to seed an alignment.
Native alignment scores are not calibrated probabilities. Protein HMM
E-values describe a peptide selected by the alignment procedure; they are
not calibrated genome-search E-values. Viral-reference validation and coverage
filters reduce unsupported candidates but do not establish a false-discovery
rate. See the [method](reference/frameshift-method.md) for scoring and limits.
