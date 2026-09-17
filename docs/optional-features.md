# Optional analyses

The core workflow runs on CPUs. Frameshift screening uses the standard
environment; the other optional analyses need separate tools or data.

## Frameshift screening

This screen can recover viral marker domains that protein prediction misses.
Enable it with `--frameshift-screening`. See the
[frameshift screening guide](FRAMESHIFT_SCREENING.md) for an example and
the method's limits.

## TMVec2

TMVec2 adds BFVD protein-embedding evidence. It can use a CPU; a GPU is
optional.

```bash
pixi install --locked -e structural
pixi run virosync orchestrate setup \
  --config config/orchestration.yaml \
  --no-write-config \
  --tmvec
pixi run -e structural check-structural-runtime --require-tmvec
pixi run -e structural virosync \
  -i genome.fna -o results/tmvec \
  --config config/orchestration.yaml \
  --device cpu --tmvec
```

For a GPU, replace `--device cpu` with `--device cuda --tmvec-gpu`. An enabled
TMVec2 run stops if its runtime or resources fail validation.
To check CUDA with `check-structural-runtime`, first set `compute.device: cuda`
in the config.

## Boltz and Foldseek

Install the isolated Boltz runtime:

```bash
pixi install --locked --manifest-path tools/boltz_runtime/pixi.toml
```

Set the Foldseek database prefix and MSA-server mode in the config:

```yaml
phase3:
  viral_structure_db: /data/viral_structures/db
  boltz_use_msa_server: true
```

This sends query protein sequences to an online MSA service.

Check and run the analysis:

```bash
pixi run check-structural-runtime --require-boltz
pixi run virosync \
  -i genome.fna -o results/boltz \
  --config config/orchestration.yaml \
  --boltz
```

An MCP fold inferred without a query structure is not a confirmed 3D
structure. Use the Boltz and Foldseek results when query-structure evidence is
required. If the runtime or database is missing, ViroSync warns and skips this
layer.

## InterProScan

Set `phase3.interproscan_dir` to an InterProScan directory that contains an
executable `interproscan.sh`. Then check and run the layer:

```bash
pixi run check-structural-runtime --require-interproscan
pixi run virosync \
  -i genome.fna -o results/interproscan \
  --config config/orchestration.yaml \
  --interproscan
```

If the runtime is missing, ViroSync warns and skips this layer.

## Check all optional layers

```bash
pixi run -e structural check-structural-runtime --require-all-optional
```

This command checks configured tools and data. It does not install them.
