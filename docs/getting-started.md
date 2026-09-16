# Getting started

ViroSync requires Linux x86-64, [Git](https://git-scm.com/), and
[Pixi](https://pixi.sh/). The core workflow uses CPUs and does not need a GPU.

## 1. Install ViroSync and its resources

Allow about 19 GB for resource setup, plus space for the Pixi environment and
results. The resources download about 6 GB and occupy about 13 GB after setup.

```bash
git clone https://github.com/NeLLi-team/virosync.git
cd virosync
pixi install --locked
pixi run setup-virosync-resources
```

Setup installs resource bundle v1.0.7 in `resources/virosync` by default.
See [Databases](RESOURCE_BUNDLE.md) to use another location.

## 2. Run the bundled example

```bash
pixi run example
cat results/example/batch_summary.tsv
```

The `test-1` row must report `status=success`, `predictions=6`, and
`accepted=1` with resource bundle v1.0.7. `predictions` counts all
candidates and `accepted` counts the reported EVEs. The example uses one
worker with eight CPU threads.

Rerunning the example replaces `results/example/test-1/`. For your own analyses, see
[Resume and restart](reference/cli.md#resume-and-restart).

## 3. Inspect the results

```bash
head -n 2 results/example/test-1/phase3_synthesis/virosync_predictions.tsv
```

This prints the header and the accepted EVE. The other main files are:

| File under `results/example/` | Contents |
| --- | --- |
| `batch_summary.tsv` | Status and candidate counts for each genome |
| `test-1/phase3_synthesis/virosync_predictions.gff3` | Accepted EVE coordinates and annotations |
| `test-1/virosync_predictions_detailed.tsv` | All six candidates and their evidence |

TSV coordinates are 0-based and half-open; GFF3 coordinates are 1-based and
inclusive. See [Outputs and interpretation](METHODS.md) for column definitions
and limits of the confidence scores.

## 4. Run your genome

Replace `genome.fna` with your assembly path and run from the repository
directory:

```bash
pixi run virosync \
  -i genome.fna \
  -o results/my_genome \
  --config config/orchestration.yaml \
  -w 1 \
  --threads-per-worker 8
```

For several genomes, set `-i` to a directory. ViroSync accepts uncompressed
`.fna`, `.fasta`, and `.fa` files. Directory input scans
only top-level files with these lowercase extensions. You can also pass a text
file with one FASTA path per line. Increase `-w` to process genomes in
parallel; budget `workers × threads-per-worker` CPU threads.

The output root contains `batch_summary.tsv` and one directory per genome.
See [Outputs and interpretation](METHODS.md) for the files written when a
genome has no candidates.
