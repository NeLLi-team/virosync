# ViroSync

![Version](https://img.shields.io/badge/version-1.0.1-blue)
[![License: non-commercial use only](https://img.shields.io/badge/license-non--commercial-orange.svg)](LICENSE)

ViroSync finds candidate endogenous viral elements (EVEs) in assembled
eukaryotic genomes. It detects viral markers, refines boundaries, and provides
accepted regions with their evidence.

## Workflow

![ViroSync workflow](docs/virosync_workflow.png)

ViroSync runs on CPUs by default. A GPU is not required for the core workflow.

## Install ViroSync

ViroSync supports Linux x86-64 and uses [Pixi](https://pixi.sh/).

```bash
git clone https://github.com/NeLLi-team/virosync.git
cd virosync
pixi install --locked
pixi run setup-prodigal-gv
pixi run setup-virosync-resources
```

Resources download about 6 GB, use about 13 GB when installed, and need about
19 GB during setup. Allow more space for the Pixi environment and results.
Native setup builds the pinned Prodigal-GV capacity correction. The
[installation guide](docs/getting-started.md) covers its runtime directory and
direct Python or console use.

## Quick start

```bash
pixi run example
cat results/example/batch_summary.tsv
```

The `test-1` row must report `status=success`, `predictions=6`, and
`accepted=1` with resource bundle v1.0.7.

## Run your own genomes

```bash
pixi run virosync \
  -i genomes/ \
  -o results/my_run \
  --config config/orchestration.yaml \
  -w 4 \
  --threads-per-worker 16
```

`-i` accepts one uncompressed `.fna`, `.fasta`, or `.fa` file, a directory of
these files, or a text file with one FASTA path per line.

The output root contains `batch_summary.tsv` with status and counts for each
genome. Each genome directory contains:

- `phase3_synthesis/virosync_predictions.tsv`: accepted EVE calls.
- `phase3_synthesis/virosync_predictions.gff3`: accepted EVE coordinates and annotations.
- `virosync_predictions_detailed.tsv`: all candidates and their evidence.

Frameshift-sensitive marker detection is off by default. Add
`--frameshift-screening` to enable it. See the
[frameshift screening guide](docs/FRAMESHIFT_SCREENING.md) for usage and limits.

See the [documentation](https://nelli-team.github.io/virosync/) for output
details, command options, and optional analyses.

ViroSync 1.0.1 uses resource bundle v1.0.7.

ViroSync is available for non-commercial use under [LICENSE](LICENSE).
