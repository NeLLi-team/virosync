# Command-line reference

Run ViroSync through Pixi from the repository root:

```bash
pixi run virosync [GLOBAL OPTIONS] COMMAND [COMMAND OPTIONS]
```

`pixi run virosync -i INPUT -o OUTPUT` is a shortcut for
`pixi run virosync run -i INPUT -o OUTPUT`. CLI values override YAML values only
when you supply them. Use `-h` or `--help` on any command.

<!-- cli-reference:virosync -->
## Global options

Global options must occur before the command. The bare `-i` and `-o` shortcut
also accepts them before the input option.

| Option | Behavior |
| --- | --- |
| `-v`, `--verbose` | Show diagnostic output. A command-local `-v` is also available on `run` and `orchestrate setup`. |
| `-q`, `--quiet` | Hide the banner, progress display, and final summary. Errors remain visible. |
| `--version` | Print the ViroSync version and exit. |
| `-h`, `--help` | Print help and exit. |

<!-- cli-reference:virosync-info -->
## Show system information

Print the ViroSync and database versions, Python and PyTorch versions, CUDA
state, CPU count, and installed memory.

| Option | Behavior |
| --- | --- |
| `--config PATH` | Read the database root from this config. Default: `config/orchestration.yaml`. `VIROSYNC_DB_ROOT` has priority. |

<!-- cli-reference:virosync-run -->
## Run ViroSync

Run one genome, a directory of FASTA files, or a text file with one FASTA path
per line. A directory scan uses only its top level. Relative paths in a list
file start at the current working directory. `virosync orchestrate run` calls
the same command.

### Input and execution

| Option | Behavior |
| --- | --- |
| `-i PATH`, `--input PATH` | Required input. Accepts `.fna`, `.fasta`, or `.fa`, a directory, or a list file. The path must exist. |
| `-o PATH`, `--output PATH` | Required output root. ViroSync creates one subdirectory per genome. |
| `--config PATH` | Read orchestration and pipeline defaults from this YAML file. The path must exist. |
| `--clean-run` | Start again and do not reuse completed output. Without this option, ViroSync validates run state before resume. |
| `-w N`, `--workers N` | Set the number of genome slots. Minimum: 1. Default: config value or 4. |
| `--threads-per-worker N` | Set tool threads for each genome. Minimum: 1. Default: config value or 8. |
| `--max-concurrent-genomes N` | Cap genomes in flight. Minimum: 1. It must equal `--workers` when both options are supplied. |
| `-v`, `--verbose` | Show the effective config and diagnostic logs instead of the progress display. |

### Database and tool paths

All path overrides in this table must exist.

| Option | Behavior |
| --- | --- |
| `--hmm-db PATH` | Override the viral-marker HMM database. |
| `--hmm-allowlist PATH` | Override the HMM allowlist. |
| `--marker-faa-db PATH` | Set the marker-protein FASTA build input. This has priority over `--marker-faa-dir`. |
| `--marker-faa-dir PATH` | Set the directory of marker-protein FASTA build inputs. |
| `--marker-db PATH` | Override the prebuilt marker DIAMOND database used for Phase 1 validation. `--rebuild-db` ignores it. |
| `--faa-dir PATH` | Set the protein FASTA directory used for a run-local marker database build. |
| `--gvclass-db PATH` | Override the GVClass database. |
| `--gvclass PATH` | Set the GVClass install directory and enable batch classification. Default: `VIROSYNC_GVCLASS_PATH`. The directory must contain an executable `gvclass`. A tool error omits `gvclass_results.tsv`, but the core run continues. |
| `--diamond-db PATH` | Override the Phase 3 phylogenetic DIAMOND database. |

### Marker discovery and scoring

| Option | Behavior |
| --- | --- |
| `--enable-phylogenetic`, `--disable-phylogenetic` | Enable or disable GVClass and DIAMOND phylogenetic validation. |
| `--assembly-mode MODE` | Set HMM seed handling to `default`, `fragmented`, `relaxed`, or `strict`. |
| `--high-tier-threshold FLOAT` | Set the HIGH confidence cutoff from 0 to 1. Default: config value or 0.7. |
| `--low-tier-threshold FLOAT` | Set the LOW confidence cutoff from 0 to 1. Default: config value or 0.2. It must be below the HIGH cutoff. |
| `--hmm-chunk-size N` | Set the number of predicted open reading frames per HMM chunk. Minimum: 1. |
| `--rebuild-db`, `--no-rebuild-db` | Force or avoid a run-local marker database build. Forced mode ignores `--marker-db` and needs an FAA directory plus one marker FASTA source, from CLI or config. Without a prebuilt database, ViroSync must build one even when the flag is off. |
| `--phase1-initial-window-bp N` | Set the first marker-cluster window in base pairs. Minimum: 1. Default: 10000. |
| `--phase1-initial-window-genes N` | Set the first marker-cluster window in genes. Minimum: 1. Default: 5. |
| `--phase1-min-markers-initial N` | Set the minimum marker count for an initial cluster. Minimum and default: 1. |
| `--phase1-extension-kb N` | Extend from the outer markers by this many kilobases. Minimum: 0. Default: 5. |
| `--phase1-merge-distance N` | Merge overlapping regions within this base-pair gap. Minimum: 0. Default: 1000. |
| `--frameshift-screening`, `--no-frameshift-screening` | Enable or disable BATH marker rescue. BATH must be on `PATH` when enabled. |

### Compute and optional evidence

| Option | Behavior |
| --- | --- |
| `--device DEVICE` | Select `cpu` or `cuda`. Default: config value or `cpu`. |
| `--search-backend diamond` | Select the sequence search backend. `diamond` is the only accepted value. |
| `--gpu-id N` | Select a zero-based GPU. Sets `CUDA_VISIBLE_DEVICES` and `VIROSYNC_GPU`. |
| `--skip-masking`, `--no-skip-masking` | Force masking off or enable TRF plus RepeatMasker. `--no-skip-masking` needs exactly one `execution.masking.repeatmasker_species` or `repeatmasker_library` value. |
| `--skip-structural`, `--no-skip-structural` | Skip or run Boltz and Foldseek structural homology. `--boltz` clears the skip unless you set it explicitly. |
| `--boltz`, `--no-boltz` | Enable or disable Boltz and Foldseek. Missing runtime or data disables the layer with a warning. |
| `--tmvec`, `--no-tmvec` | Enable or disable TMVec2 BFVD search. An enabled search fails before analysis when its runtime or resources do not pass validation. |
| `--tmvec-gpu`, `--no-tmvec-gpu` | Require or do not require TMVec2 on CUDA. `--tmvec-gpu` also enables TMVec2 and selects CUDA. It cannot be combined with `--device cpu` or `--no-tmvec`. |
| `--interproscan`, `--no-interproscan` | Enable or disable InterProScan. Missing runtime or data disables the layer with a warning. |
| `--use-taxonomy-ml`, `--no-taxonomy-ml` | Enable or disable Phase 2 taxonomy-boundary machine learning. |
| `--taxonomy-ml-model MODEL` | Select `logreg`, `gbdt`, or `xgboost` for taxonomy-boundary refinement. |

<!-- cli-reference:virosync-orchestrate -->
## Orchestration commands

The orchestration group contains setup, resource verification, system
information, and the `run` command.

<!-- cli-reference:virosync-orchestrate-setup -->
## Install resources

Install core resources and optional analysis resources.

An explicit TMVec2 install activates its target only after the bundle, model,
manifest, and database checks pass. A failed TMVec2 install exits with status
1. InterProScan setup remains optional and reports an unavailable archive as a
warning.

| Option | Behavior |
| --- | --- |
| `--config PATH` | Read or update this config. Default: `config/orchestration.yaml`. |
| `--db-root PATH` | Install the stable core-resource link at this path. `VIROSYNC_DB_ROOT` is also supported. |
| `--core-resource PATH_OR_URL` | Use this core archive instead of the configured source. Custom sources do not inherit the shipped digest pins. |
| `--core-version TEXT` | Set the expected core-resource version for a custom source. |
| `--core-resource-sha256 HEX` | Set the expected SHA-256 for the full core archive. |
| `--core-manifest-sha256 HEX` | Set the expected SHA-256 for `RESOURCE_MANIFEST.json` inside the archive. |
| `--tmvec`, `--no-tmvec` | Install or skip the configured TMVec2 BFVD resource set. |
| `--tmvec-url PATH_OR_URL` | Override the configured TMVec2 BFVD bundle. Use it with `--tmvec-resource-sha256`. |
| `--tmvec-resource-sha256 HEX` | Set the SHA-256 for a custom TMVec2 BFVD bundle. |
| `--tmvec-dir PATH` | Set the TMVec2 target. Priority: CLI path, `phase3.tmvec_database_dir`, then `resources/virosync-optional/tmvec` beside the core resource tree. |
| `--interproscan-url PATH_OR_URL` | Install a user-supplied InterProScan archive. Use it with `--interproscan-resource-sha256`. |
| `--interproscan-resource-sha256 HEX` | Set the SHA-256 for the InterProScan archive. |
| `--interproscan-dir PATH` | Set the InterProScan target. |
| `--boltz-db-dir PATH` | Record a Foldseek viral-structure database prefix for Boltz. |
| `--interactive-optional`, `--no-interactive-optional` | Enable or disable optional-resource prompts. Prompts appear only in an interactive terminal. Default: disabled. |
| `--force` | Reinstall resources even when the target exists. |
| `--write-config`, `--no-write-config` | Enable or disable writes of resolved paths to the config. Default: write. |
| `-v`, `--verbose` | Show source, target, config, and validation details. |

`pixi run setup-virosync-resources` calls setup with
`--no-interactive-optional --no-write-config`.

<!-- cli-reference:virosync-orchestrate-resources -->
## Resource commands

The resource group contains the `verify` command.

<!-- cli-reference:virosync-orchestrate-resources-verify -->
## Verify core resources

Check the installed core-resource identity.

| Option | Behavior |
| --- | --- |
| `--config PATH` | Read the expected version and manifest digest. Default: `config/orchestration.yaml`. |
| `--db-root PATH` | Verify this stable resource path. `VIROSYNC_DB_ROOT` is also supported. |
| `--full` | Hash all nine payloads and run semantic DIAMOND checks. Without `--full`, verification uses authenticated metadata and receipts. |

<!-- cli-reference:virosync-orchestrate-info -->
## Show orchestration information

Print the orchestration backend, ViroSync version, and input forms. This command
has no command-specific options.

## Structural preflight

`check-structural-runtime` is a Pixi task backed by
`scripts/check_structural_runtime.py`.

| Option | Behavior |
| --- | --- |
| `--config PATH` | Read optional-feature paths and states from this config. Default: `config/orchestration.yaml`. |
| `--require-tmvec` | Exit with status 1 unless the TMVec2 runtime, device-specific upstream vector, and real BFVD query pass. CPU is valid. `compute.device: cuda` also requires CUDA. |
| `--require-boltz` | Exit with status 1 unless Boltz, Foldseek, the MSA setting, and the Foldseek database pass. |
| `--require-interproscan` | Exit with status 1 unless an executable `interproscan.sh` exists in the configured directory. |
| `--require-all-optional` | Apply all three required checks. |

With no required option, the script checks only optional layers enabled in the
config. It does not download any model or database.
