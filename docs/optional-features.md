# Optional analyses

The core workflow runs on CPUs. Enable an optional analysis with its run
option after installing the required tools and data.

| Analysis | Run option | Required tools and data |
| --- | --- | --- |
| Frameshift screening | `--frameshift-screening` | BATH and the shipped VS marker profiles. See the [worked example](FRAMESHIFT_SCREENING.md). |
| TMVec2 | `--tmvec` | The `structural` Pixi environment and the TMVec2 BFVD resource bundle. Runs on CPU or GPU. |
| Boltz and Foldseek | `--boltz` | The isolated Boltz runtime, Foldseek, a viral-structure database, and an online MSA service. |
| InterProScan | `--interproscan` | An InterProScan installation with an executable `interproscan.sh`. |

TMVec2, Boltz/Foldseek, and InterProScan add evidence to existing candidate
regions. They do not discover regions on their own.

## TMVec2 requirements

Install the `structural` environment with `pixi install --locked -e structural`.
The resource installer accepts `--tmvec` to download the configured BFVD bundle;
see [resource setup options](reference/cli.md#install-resources).

Run ViroSync in the `structural` environment. Select `--device cpu` for CPU
inference, or `--device cuda --tmvec-gpu` to require a GPU. An enabled TMVec2
analysis stops if its runtime or resources fail validation.

The runtime check uses `compute.device` from the config. Set it to `cuda` to
check CUDA support before a GPU run.

## Boltz and Foldseek requirements

Boltz has a separate Pixi manifest at `tools/boltz_runtime/pixi.toml`. Install
it with `pixi install --locked --manifest-path tools/boltz_runtime/pixi.toml`.
In the run config, set `phase3.viral_structure_db` to the Foldseek database
prefix and `phase3.boltz_use_msa_server` to `true`.

MSA-server mode sends query protein sequences to an online service. If the
runtime or database is missing, ViroSync warns and skips this analysis.

An inferred MCP fold does not confirm the protein's 3D structure. Boltz and
Foldseek provide computational query-structure evidence, not experimental
confirmation.

## InterProScan requirements

Set `phase3.interproscan_dir` to the directory containing an executable
`interproscan.sh`. If the runtime is missing, ViroSync warns and skips this
analysis.

## Runtime checks

The `check-structural-runtime` Pixi task checks configured tools and data. It
does not install tools or download models and databases. A failed required
check returns a nonzero exit status.

| Option | Check |
| --- | --- |
| `--require-tmvec` | Run a protein query against BFVD in the `structural` environment. |
| `--require-boltz` | Check the Boltz runtime, Foldseek, MSA-server setting, and database prefix. |
| `--require-interproscan` | Check that the configured `interproscan.sh` exists and is executable. |
| `--require-all-optional` | Apply all three checks in the `structural` environment. |

Without a required option, the task checks only analyses enabled in the
config. See the [command-line reference](reference/cli.md#check-optional-tools)
for the config option and command syntax.
