# Get started

ViroSync requires Linux x86-64, [Git](https://git-scm.com/), and
[Pixi](https://pixi.sh/). The core workflow uses CPUs and does not need a GPU.

## Install ViroSync

```bash
git clone https://github.com/NeLLi-team/virosync.git
cd virosync
pixi install --locked
```

## Install the core resources

Resources download about 6 GB, use about 13 GB when installed, and need about
19 GB during setup. Allow more space for the Pixi environment and results.

```bash
pixi run setup-virosync-resources
```

Setup installs resource bundle v1.0.7 in `resources/virosync` by default.
See [Databases](RESOURCE_BUNDLE.md) to use another location.

Verify the install:

```bash
pixi run virosync orchestrate resources verify \
  --config config/orchestration.yaml \
  --full
```

The output must include `Version: v1.0.7`.

## Run the example

```bash
pixi run example
cat results/example/batch_summary.tsv
```

The `test-1` row must report `status=success`, `predictions=6`, and
`accepted=1`.

## Run your genome

Use one FASTA file:

```bash
pixi run virosync \
  -i genome.fna \
  -o results/my_genome \
  --config config/orchestration.yaml \
  -w 1 \
  --threads-per-worker 16
```

For several genomes, set `-i` to a directory and increase `-w`. ViroSync
accepts uncompressed `.fna`, `.fasta`, and `.fa` files. Directory input scans
only top-level files with these lowercase extensions. You can also pass a text
file with one FASTA path per line.

The output root contains `batch_summary.tsv` and one directory per genome.
Each genome has accepted calls in
`phase3_synthesis/virosync_predictions.tsv` and
`phase3_synthesis/virosync_predictions.gff3`. Use
`virosync_predictions_detailed.tsv` to inspect all scored candidates.

See [Outputs](METHODS.md#output-specification) for the other result
files and the [command-line reference](reference/cli.md) for all options.
