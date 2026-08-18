# ViroSync

![Version](https://img.shields.io/badge/version-1.0.0-blue)
[![License: non-commercial use only](https://img.shields.io/badge/license-non--commercial-orange.svg)](LICENSE)

ViroSync identifies candidate endogenous viral elements (EVEs) in assembled
eukaryotic genomes. It searches for viral markers, refines candidate boundaries
with gene taxonomy and host-like sequence signals, and scores the resulting
regions.

## Workflow

![ViroSync workflow](docs/virosync_workflow.png)

ViroSync processes each genome in four phases:

1. Predict proteins after optional repeat masking.
2. Find and validate viral markers, then assemble seed regions.
3. Refine region boundaries with gene taxonomy and host-like sequence signals.
4. Score each candidate and write the accepted EVE set and reports.

## Install ViroSync

ViroSync supports Linux x86-64 and uses [Pixi](https://pixi.sh/) for its pinned
environment.

```bash
git clone https://github.com/NeLLi-team/virosync.git
cd virosync
pixi install --locked
pixi run setup-virosync-resources
```

Resource setup downloads and verifies `resources_v1_0_7_runtime.tar.gz`. See
[Get started](docs/getting-started.md) for the installation checks.

## Quick start

Run the shipped example:

```bash
pixi run virosync \
  -i example/ \
  -o results/example_16t \
  --config config/orchestration.yaml \
  -w 1 \
  --threads-per-worker 16 \
  --clean-run
```

Check the batch summary:

```bash
cat results/example_16t/batch_summary.tsv
```

The `test-1` row must have `status=success`, `predictions=6`, and `accepted=1`
with database version `v1.0.7`.

The [documentation](docs/index.md) covers input modes, command-line options,
outputs, performance tests, and optional analyses.

ViroSync is available for non-commercial use under the terms in [LICENSE](LICENSE).
