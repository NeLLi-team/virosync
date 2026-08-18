# Get started

This tutorial installs ViroSync, verifies its resources, and runs the shipped
example. It requires Linux x86-64, [Git](https://git-scm.com/), and
[Pixi](https://pixi.sh/).

## Install ViroSync

```bash
git clone https://github.com/NeLLi-team/virosync.git
cd virosync
pixi install --locked
```

The lock file pins the software environment. Do not use Conda or the system
Python for this project.

## Install the core resources

```bash
pixi run setup-virosync-resources
```

For the first interactive install, select the core-resource location and
confirm the download. For a non-interactive install, set
`VIROSYNC_DB_ROOT` before setup.

This command downloads `resources_v1_0_7_runtime.tar.gz`, verifies its archive
and manifest digests, and installs the read-only v1.0.7 resource tree under
`resources/virosync`.

Run the full resource check:

```bash
pixi run virosync orchestrate resources verify \
  --config config/orchestration.yaml \
  --full
```

The output must include:

```text
Version: v1.0.7
Authenticated payloads: 9
```

To store the resource tree outside the repository, set `VIROSYNC_DB_ROOT`
before setup:

```bash
export VIROSYNC_DB_ROOT=/data/virosync-db
pixi run setup-virosync-resources
```

Keep the variable set for later runs.

## Run the example

```bash
pixi run example
```

The task runs the genomes in `example/` and writes to `results/example/`. Check
the summary:

```bash
cat results/example/batch_summary.tsv
```

The `test-1` row must have `status=success`, `predictions=6`, and `accepted=1`.

The task passes `--clean-run`, which ignores existing outputs and starts from
scratch. To test resume, run the same input without that flag:

```bash
pixi run virosync \
  -i example/ \
  -o results/example \
  --config config/orchestration.yaml \
  -w 1 \
  --threads-per-worker 8
```

ViroSync reuses the completed genome only after it validates the run state and
all recorded output files.

## Run your genome

Point `-i` at a single FASTA file:

```bash
pixi run virosync \
  -i genome.fna \
  -o results/my_genome \
  --config config/orchestration.yaml \
  -w 1 \
  --threads-per-worker 16
```

Or at a directory of FASTA files:

```bash
pixi run virosync \
  -i genomes/ \
  -o results/my_genomes \
  --config config/orchestration.yaml \
  -w 4 \
  --threads-per-worker 16
```

ViroSync reads `.fna`, `.fasta`, and `.fa` files. It does not search
subdirectories. `-i` also accepts a text file with one FASTA path per line.
`-w` sets how many genomes run at the same time. `--threads-per-worker` sets
the threads for each of those genomes.

See the [command-line reference](reference/cli.md) for run controls and
[methods and outputs](METHODS.md) for the result files.
