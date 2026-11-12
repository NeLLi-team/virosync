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

Input can be one uncompressed genome FASTA, a directory of FASTA files, or a
text file with one FASTA path per line. A directory scan uses only its top
level and matches the lowercase extensions `.fna`, `.fasta`, and `.fa`.
Relative paths in a list file start at the current working directory.
`virosync orchestrate run` calls the same command.

### Input and execution

| Option | Behavior |
| --- | --- |
| `-i PATH`, `--input PATH` | Required input. Accepts an uncompressed `.fna`, `.fasta`, or `.fa` file, a directory, or a list file. The path must exist. |
| `-o PATH`, `--output PATH` | Required output root. ViroSync creates one subdirectory per genome. |
| `--config PATH` | Read orchestration and pipeline defaults from this YAML file. The path must exist. |
| `--clean-run` | Start again without reusing completed output. Otherwise, resume when the input, settings and saved outputs still match. |
| `--reuse-phase1-from PATH` | Start from the Phase 0 and Phase 1 results of an earlier run. `PATH` is the `-o` folder of the earlier run and must exist. The new `-o` folder must not equal it, lie inside it or contain it. Needs `--reuse-phase1-code-sha256`. See [Reuse Phase 0 and Phase 1 of an earlier run](../METHODS.md#reuse-phase-0-and-phase-1-of-an-earlier-run). |
| `--reuse-phase1-code-sha256 SHA` | Source SHA-256 of the earlier run: its `identities.code.source_sha256` in `virosync_run_state.json`. This is your statement that the Phase 0 and Phase 1 code of that source is identical to the current source, except for the host signature model. ViroSync cannot check this. Needs `--reuse-phase1-from`. |
| `--reuse-phase1-snapshot PATH` | Authenticate the Phase 1 files that Phases 2 and 3 read, for an earlier run that recorded no hash for them. `PATH` is a `sha256sum` listing made inside the `-o` folder of the earlier run. Needs `--reuse-phase1-snapshot-sha256` and `--reuse-phase1-from`. See [Reuse Phase 0 and Phase 1 of an earlier run](../METHODS.md#reuse-phase-0-and-phase-1-of-an-earlier-run). |
| `--reuse-phase1-snapshot-sha256 SHA` | SHA-256 of the file given to `--reuse-phase1-snapshot`. ViroSync refuses a listing with another hash. Needs `--reuse-phase1-snapshot`. |
| `--reuse-phase2-search` | Import the Phase 2 gene-taxonomy search of the earlier run (`phase2/boundary_diamond/genome_boundary_query.faa` and its hit tables, `boundary_diamond.tsv` or `diamond_chunk_1.tsv` to `diamond_chunk_K.tsv`) and start no new search. The snapshot supplies the SHA-256 of each file, the table set must be complete, and the gene-taxonomy database must be the one of the earlier run. Give the flag on every start of the folder. Every failed condition refuses the run. Needs `--reuse-phase1-from`. See [Reuse the Phase 2 gene-taxonomy search](../METHODS.md#reuse-the-phase-2-gene-taxonomy-search). |
| `--plots`, `--no-plots` | Enable or disable notebook and figure rendering. Default: `execution.generate_plots: false`. Scientific tables, sequence exports, evidence and audits remain required in either mode. |
| `-w N`, `--workers N` | Run up to N genomes at once. Minimum: 1. The shipped config sets 6. |
| `--threads-per-worker N` | Set tool threads for each genome. Minimum: 1. Default: `compute.threads`, or 8 if unset. ViroSync lowers the value to `compute.max_threads` (48 in the shipped config) divided by the genome limit and prints a warning. |
| `--max-concurrent-genomes N` | Set the same genome limit as `--workers`. Minimum: 1. Both values must agree if you supply both. |
| `-v`, `--verbose` | Show the effective config and diagnostic logs instead of the progress display. |

Plot-enabled runs require a valid notebook before completion. Changing `--plots`
changes the run fingerprint and can rerun scientific stages. To add plots to a
completed output without rerunning the pipeline, call
`virosync.report.generate.generate_eve_report(output_dir=Path("results/GENOME"),
genome_id="GENOME", tax_labels_path=Path("taxonomy_labels.tsv"),
host_credit_threshold=0.35)` from Python. Supply the taxonomy labels and host
threshold used for that run. This standalone renderer adds the notebook and
figures; it does not change the recorded pipeline plot mode.

Without `--clean-run`, reuse requires matching run identities and validated
checkpoints and outputs, including the repeat-candidate table. See
[Output compatibility and resume](../METHODS.md#output-compatibility-and-resume)
for the required files and schema versions.

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
| `--phase1-initial-window-bp N` | Set the first marker-cluster window and the maximum gap to an independently validated neighbor for novel-marker support, in base pairs. Minimum: 1. Default: 10000. |
| `--phase1-initial-window-genes N` | Set the first marker-cluster window in genes. Minimum: 1. Default: 5. |
| `--phase1-min-markers-initial N` | Set the minimum marker count for an initial cluster. Minimum and default: 1. |
| `--phase1-extension-kb N` | Extend from the outer markers by this many kilobases. Minimum: 0. Default: 5. |
| `--phase1-merge-distance N` | Merge overlapping regions within this base-pair gap. Minimum: 0. Default: 1000. |
| `--frameshift-screening`, `--no-frameshift-screening` | Enable or disable native codon-aware marker rescue. Off by default; uses the standard Pixi environment. |

### Compute and optional evidence

| Option | Behavior |
| --- | --- |
| `--device DEVICE` | Select `cpu` or `cuda`. Default: config value or `cpu`. |
| `--search-backend diamond` | Select the sequence search backend. `diamond` is the only accepted value. |
| `--gpu-id N` | Select a GPU by its zero-based index. |
| `--skip-masking`, `--no-skip-masking` | Force masking off or enable TRF plus RepeatMasker. `--no-skip-masking` needs exactly one `execution.masking.repeatmasker_species` or `repeatmasker_library` value. |
| `--skip-structural`, `--no-skip-structural` | Skip or run Boltz and Foldseek structural homology. `--boltz` clears the skip unless you set it explicitly. |
| `--boltz`, `--no-boltz` | Enable or disable Boltz and Foldseek. Missing runtime or data disables the layer with a warning. |
| `--tmvec`, `--no-tmvec` | Enable or disable TMVec2 BFVD search. An enabled search fails before analysis when its runtime or resources do not pass validation. |
| `--tmvec-gpu`, `--no-tmvec-gpu` | Require or do not require TMVec2 on CUDA. `--tmvec-gpu` also enables TMVec2 and selects CUDA. It cannot be combined with `--device cpu` or `--no-tmvec`. |
| `--interproscan`, `--no-interproscan` | Enable or disable InterProScan. Missing runtime or data disables the layer with a warning. |
| `--use-taxonomy-ml`, `--no-taxonomy-ml` | Enable or disable per-gene ML support for seed-preserving outward extension. |
| `--taxonomy-ml-model MODEL` | Select `logreg`, `gbdt`, or `xgboost` for that support model. Host stops and seed containment still apply. |

The core [direct-repeat assessment](../METHODS.md#direct-repeat-assessment) runs
by default and has no separate enable flag. It reads unmasked host sequence;
`--skip-masking` controls preprocessing and does not disable this assessment.

### Reference masking

Reference masking drops search hits to chosen reference genomes. It is not the
repeat masking of `--skip-masking`. By default a run does not use the proteins
of its own assembly as evidence. [Reference masking](../METHODS.md#reference-masking)
describes the mask. Changing any mask option or its config value
(`host.reference_mask_taxa`, `host.mask_own_assembly`, `host.self_reference_keys`)
starts the genome again.

| Option | Behavior |
| --- | --- |
| `--mask-taxon NAME` | Ignore every reference genome whose lineage names NAME at any rank, ignoring case. Repeat the option for several taxa. Needs the taxonomy labels file. Default: none. |
| `--mask-own-assembly`, `--no-mask-own-assembly` | Ignore or keep reference hits to the assembly named by explicit self keys, or by a `GCA_` or `GCF_` accession in the input file name when no keys are supplied. Default: on. |
| `--self-reference-key KEY` | Name a reference key of the genome's own proteomes, such as `GCA-000002455-1` or a JGI, IMG or GWH key. Repeat the option for several keys. The keys replace the accession read from the input file name. Config: `host.self_reference_keys`. |

Explicit self keys require exactly one input genome and own-assembly masking on.
Exact NCBI keys include both their hyphen and underscore accession forms; other
keys remain exact. Empty keys, protein IDs, placeholder taxa and reference-class
tags are rejected.

<!-- cli-reference:virosync-orchestrate -->
## Orchestration commands

The orchestration group contains setup, resource verification, system
information, and the `run` command.

<!-- cli-reference:virosync-orchestrate-setup -->
## Install resources

Install databases and optional analysis resources.

| Option | Behavior |
| --- | --- |
| `--config PATH` | Read or update this config. Default: `config/orchestration.yaml`. |
| `--db-root PATH` | Install the core databases here. You can also set `VIROSYNC_DB_ROOT`. |
| `--core-resource PATH_OR_URL` | Use a different core archive. Supply its version, archive SHA-256 and manifest SHA-256 with the options below. |
| `--core-version TEXT` | Set the expected core-resource version for a custom source. |
| `--core-resource-sha256 HEX` | Set the expected SHA-256 for the full core archive. |
| `--core-manifest-sha256 HEX` | Set the expected SHA-256 for `RESOURCE_MANIFEST.json` inside the archive. |
| `--tmvec`, `--no-tmvec` | Install or skip the configured TMVec2 BFVD resource set. |
| `--tmvec-url PATH_OR_URL` | Override the configured TMVec2 BFVD bundle. Use it with `--tmvec-resource-sha256`. |
| `--tmvec-resource-sha256 HEX` | Set the SHA-256 for a custom TMVec2 BFVD bundle. |
| `--tmvec-dir PATH` | Install TMVec2 here. Overrides `phase3.tmvec_database_dir`. If unset, setup chooses a directory beside the core databases. |
| `--interproscan-url PATH_OR_URL` | Install a user-supplied InterProScan archive. Use it with `--interproscan-resource-sha256`. |
| `--interproscan-resource-sha256 HEX` | Set the SHA-256 for the InterProScan archive. |
| `--interproscan-dir PATH` | Set the InterProScan target. |
| `--boltz-db-dir PATH` | Record a Foldseek viral-structure database prefix for Boltz. |
| `--interactive-optional`, `--no-interactive-optional` | Enable or disable optional-resource prompts. Prompts appear only in an interactive terminal. Default: disabled. |
| `--force` | Reinstall resources even when the target exists. |
| `--write-config`, `--no-write-config` | Enable or disable writes of resolved paths to the config. Default: write. |
| `-v`, `--verbose` | Show source, target, config, and validation details. |

<!-- cli-reference:virosync-orchestrate-resources -->
## Resource commands

The resource group contains the `verify` command.

<!-- cli-reference:virosync-orchestrate-resources-verify -->
## Verify core resources

Check the installed databases.

| Option | Behavior |
| --- | --- |
| `--config PATH` | Read the expected version and manifest digest. Default: `config/orchestration.yaml`. |
| `--db-root PATH` | Check databases at this path. You can also set `VIROSYNC_DB_ROOT`. |
| `--full` | Check file contents as well as names, sizes and version. Takes longer than the default check. |

<!-- cli-reference:virosync-orchestrate-info -->
## Show orchestration information

Print the orchestration backend, ViroSync version, and input forms. This command
has no command-specific options.

## Check optional tools

Run `pixi run -e structural check-structural-runtime` with the options below.

| Option | Behavior |
| --- | --- |
| `--config PATH` | Read optional-feature paths and states from this config. Default: `config/orchestration.yaml`. |
| `--require-tmvec` | Test TMVec2 and its database with a protein query. Runs on CPU unless the config selects CUDA. |
| `--require-boltz` | Exit with status 1 unless Boltz, Foldseek, the MSA setting, and the Foldseek database pass. |
| `--require-interproscan` | Exit with status 1 unless an executable `interproscan.sh` exists in the configured directory. |
| `--require-all-optional` | Apply all three required checks. |

With no required option, check only tools enabled in the config. A failed
required check returns a nonzero exit status. This command does not download
models or databases.


Saved Phase 1 reuse requires authenticated consumed-file closure in the old phase markers. Historical SOURCE10/SOURCE12 files whose markers lack this closure are refused; compute a fresh Phase 1 instead. A current hash or copied provenance record cannot supply missing historical authentication, unless you pass `--reuse-phase1-snapshot` and `--reuse-phase1-snapshot-sha256`. The snapshot is a `sha256sum` listing, made now inside the `-o` folder of the earlier run (`sha256sum GENOME_ID/phase1/hmm/hmm_hits_validated.tsv ...`), with one line per file and a path relative to that folder. One listing can cover many genomes. ViroSync takes a hash from it only for the Phase 1 files that Phases 2 and 3 read, and only when the phase markers record none. A listing hash that differs from a recorded hash is refused. Other entries authenticate nothing. A snapshot proves that the files did not change between the listing and the import. It does not prove that they are the files the earlier run wrote. The provenance record of each phase names the listing hash and the files that the snapshot authenticated for that phase. Files are checked again during copying. Running old runs and symlinked roots are refused. A resumed Phase 0 must have the same recorded import origin; a fresh or different-origin Phase 0 cannot receive imported Phase 1. Empty nonterminal checkpoints are refused, while legitimate authenticated terminal-zero runs remain supported.

The patched source has a new code identity. Any semantic equivalence claim must compare Phase 0/1 code against this repair source, including the host-model and no-seed refactor. Output comparison includes scientific identity, consumed Phase 0/1 artifacts and jelly-roll exports; non-JSON files stream in bounded chunks. JSON comparison still parses complete documents and can require memory proportional to document size. No runtime override, campaign credit or acceptance is granted by these options.

Only authenticated phase files are copied. Unbound intermediates are omitted;
existing downstream-consumed files without an old marker hash cause refusal, unless the snapshot options supply the hash.
Output lookup stays inside the genome run root, including fallback paths, so a
neighboring run cannot supply missing search evidence. The comparison tool requires
gene-function tables for full runs with the current resource bundle; older bundles
without those tables cannot receive a complete-comparison result. Evidence graphs
remain conditional on coherence analyses.

The output generator defaults to its own output folder as the genome root.
Nested export callers supply the genome root explicitly; fallback lookup never
leaves that root. The comparator requires enabled masked outputs to resolve
inside their own genome folders and to distinct files. Masking-off is different:
the native writer uses the original input path. A shared unmasked input is allowed
only when the native masking validator verifies that path and its bytes against
each run's recorded input identity and explicit masking-off request.
