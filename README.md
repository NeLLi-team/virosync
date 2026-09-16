# ViroSync

![Version](https://img.shields.io/badge/version-1.0.1-blue)
[![License: non-commercial use only](https://img.shields.io/badge/license-non--commercial-orange.svg)](LICENSE)

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
results. Follow [Getting started](https://nelli-team.github.io/virosync/getting-started/)
for installation details and output files.

## Run your own genomes

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
line. `batch_summary.tsv` records the status and counts for each genome.

## Documentation

- [Getting started](https://nelli-team.github.io/virosync/getting-started/): installation and the bundled example.
- [Outputs](https://nelli-team.github.io/virosync/METHODS/): result files, coordinates, and interpretation.
- [Command-line reference](https://nelli-team.github.io/virosync/reference/cli/): options, worker counts, and resume.
- [Optional analyses](https://nelli-team.github.io/virosync/optional-features/): frameshift screening, protein domains, and structural evidence.
- [Databases](https://nelli-team.github.io/virosync/RESOURCE_BUNDLE/) and [performance](https://nelli-team.github.io/virosync/performance/).

ViroSync 1.0.1 uses resource bundle v1.0.7.

ViroSync is available for non-commercial use under [LICENSE](LICENSE).
Report problems through [GitHub Issues](https://github.com/NeLLi-team/virosync/issues).
