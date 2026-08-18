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
pixi run virosync \
  -i example/ \
  -o results/example_16t \
  --config config/orchestration.yaml \
  -w 1 \
  --threads-per-worker 16 \
  --clean-run
```

Check the summary:

```bash
cat results/example_16t/batch_summary.tsv
```

The `test-1` row must have `status=success`, `predictions=6`, and `accepted=1`.

Run the same command without `--clean-run` to test resume. ViroSync reuses the
completed genome only after it validates the run state and all recorded output
files.

## Run your genome

```bash
pixi run virosync \
  -i genome.fna \
  -o results/my_genome \
  --config config/orchestration.yaml \
  -w 1 \
  --threads-per-worker 16
```

You can also use a directory of FASTA files or a text file with one FASTA path
per line. See the [command-line reference](reference/cli.md) for run controls
and [methods and outputs](METHODS.md) for the result files.
