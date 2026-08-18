# Optional analyses

Frameshift screening, TMVec2, Boltz, and InterProScan are opt-in. Each feature
needs a separate runtime or data set.

| Feature | Install or data requirement | Run option |
| --- | --- | --- |
| Frameshift screening | Pinned BATH and Easel builds | `--frameshift-screening` |
| TMVec2 | Structural Pixi environment and the pinned Lobster-24M/BFVD resource set | `--tmvec`. Add `--tmvec-gpu` to require CUDA. |
| Boltz and Foldseek | Boltz runtime and a Foldseek database prefix | `--boltz --no-skip-structural` |
| InterProScan | An external InterProScan directory with `interproscan.sh` | `--interproscan` |

## Frameshift screening

Build the pinned BATH tools:

```bash
pixi run setup-bath
```

Test them with the shipped three-contig *Trichomonas vaginalis* G3 input:

```bash
pixi run example-frameshift
cat results/example-frameshift/batch_summary.tsv
```

The row must have `status=success`, `predictions=5`, and `accepted=2`. The
[frameshift screening guide](FRAMESHIFT_SCREENING.md) describes the rescue
files and validation steps.

## TMVec2

Install the structural environment:

```bash
pixi install --locked -e structural
```

Install the pinned TMVec2 models and matching BFVD embeddings:

```bash
pixi run virosync orchestrate setup \
  --config config/orchestration.yaml \
  --no-write-config \
  --tmvec
```

The default target is `resources/virosync-optional/tmvec`. Setup downloads the
pinned BFVD bundle and the pinned Lobster-24M and TMVec2 model files. It checks
each SHA-256 before it activates the target.

Use `--tmvec-dir PATH` to select another target. The command-line path has
priority over `phase3.tmvec_database_dir` in the config. A custom
`--tmvec-url` also needs `--tmvec-resource-sha256`. An explicit TMVec2 setup
request exits with status 1 if a download or validation step fails.

ViroSync uses Lobster-24M to create 408-column residue features. TMVec2 then
creates one 512-column protein vector. The installed manifest binds those
model revisions to the BFVD embeddings. BFVD supplies the structural evidence
used by the confidence score.

The setup command downloads the Lobster-24M weights from
[`asalam91/lobster_24M`](https://huggingface.co/asalam91/lobster_24M) and the
TMVec2 weights from
[`scikit-bio/TMVec-2`](https://huggingface.co/scikit-bio/TMVec-2). The
[Lobster code](https://github.com/prescient-design/lobster) uses the Apache-2.0
license. The TMVec2 model repository uses the BSD-3-Clause license. The
Lobster-24M model repository does not state a model-weight license.

Run the model-parity and BFVD query check:

```bash
pixi run -e structural check-structural-runtime --require-tmvec
```

The check selects the CPU or CUDA reference vector from `compute.device`. It
then runs the recorded sequence through the production BFVD search path and
requires the recorded target and score.

Run ViroSync with TMVec2:

```bash
pixi run -e structural virosync \
  -i genome.fna \
  -o results/tmvec \
  --config config/orchestration.yaml \
  --device cuda \
  --tmvec \
  --tmvec-gpu
```

Plain `--tmvec` runs on the configured CPU or CUDA device. `--tmvec-gpu`
selects CUDA and requires it. An enabled TMVec2 run stops before analysis when
its runtime or resources do not pass validation.

## Boltz and Foldseek

Install the isolated Boltz runtime:

```bash
pixi install --locked --manifest-path tools/boltz_runtime/pixi.toml
```

Set `phase3.viral_structure_db` to a Foldseek database prefix and set
`phase3.boltz_use_msa_server: true`. ViroSync writes sequence-only Boltz inputs
and does not write local multiple-sequence alignments.

Check the runtime and database:

```bash
pixi run -e structural check-structural-runtime --require-boltz
```

Enable the layer with `--boltz --no-skip-structural`. The layer needs the
runtime, database, and MSA-server setting. ViroSync disables it when one of
these inputs is absent.

## InterProScan

Point `phase3.interproscan_dir` at an InterProScan installation that contains
an executable `interproscan.sh`, then check it:

```bash
pixi run -e structural check-structural-runtime --require-interproscan
```

For an InterProScan archive, pass both `--interproscan-url` and
`--interproscan-resource-sha256`. Setup checks the archive digest and preserves
the executable bit on `interproscan.sh`.

Enable annotation with `--interproscan`. ViroSync disables the layer when the
configured program is absent.

## Check all optional layers

```bash
pixi run -e structural check-structural-runtime --require-all-optional
```

This check does not download models or databases. It exits with status 1 and
lists each missing requirement when the optional setup is incomplete.
