# ViroSync

ViroSync finds candidate endogenous viral elements (EVEs) in assembled
eukaryotic genomes. It reports their coordinates, sequences, and supporting
evidence. The core workflow runs on CPUs.

## Install and run an example

On Linux x86-64 with [Pixi](https://pixi.sh/) and Git installed:

```bash
git clone https://github.com/NeLLi-team/virosync.git
cd virosync
pixi install --locked
pixi run setup-virosync-resources
pixi run example
cat results/example/batch_summary.tsv
```

The `test-1` row must report `status=success`, `predictions=6`, and
`accepted=1` with resource bundle v1.0.7. The example uses eight CPU threads.
Allow about 19 GB for resource setup, plus space for the Pixi environment and
results.

See [Getting started](getting-started.md) for output files and input requirements.

## Run your genomes

Put one genome per uncompressed `.fna`, `.fasta`, or `.fa` file in `genomes/`,
then run from the repository directory:

```bash
pixi run virosync \
  -i genomes/ \
  -o results/my_run \
  --config config/orchestration.yaml \
  -w 1 \
  --threads-per-worker 8
```

`-i` also accepts a single FASTA file or a text file with one FASTA path per
line. See the [command-line reference](reference/cli.md) for worker counts
and resume options.

## Documentation

| Task | Page |
| --- | --- |
| Install and run the bundled example | [Getting started](getting-started.md) |
| Read the results | [Outputs and interpretation](METHODS.md) |
| Look up an option | [Command-line reference](reference/cli.md) |
| Recover frameshifted markers | [Frameshift screening](FRAMESHIFT_SCREENING.md) |
| Add domain or structural evidence | [Optional analyses](optional-features.md) |
| Install or verify databases | [Databases](RESOURCE_BUNDLE.md) |
| Inspect benchmarks | [Performance](performance.md) |

ViroSync 1.0.1 uses resource bundle v1.0.7.

ViroSync is available for non-commercial use under
[LICENSE](https://github.com/NeLLi-team/virosync/blob/main/LICENSE).
Report problems through [GitHub Issues](https://github.com/NeLLi-team/virosync/issues).
