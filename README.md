# ViroSync

![Version](https://img.shields.io/badge/version-1.0.0-blue)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

ViroSync detects candidate endogenous viral elements (EVEs) in assembled eukaryotic genomes. The active pipeline combines hallmark-marker discovery, taxonomy-guided boundary refinement, host-aware trimming, and multi-evidence confidence scoring to prioritize regions for downstream review.

## Workflow

ViroSync follows a four-phase workflow.

```text
                         ┌──────────────────────────────┐
                         │   Eukaryotic genome (FASTA)  │
                         └──────────────┬───────────────┘
                                        │
                ╔═══════════════════════▼═══════════════════════╗
                ║                   PHASE 0                     ║
                ║                Pre-processing                 ║
                ║  ─────────────────────────────────────────── ║
                ║   • optional repeat masking                   ║
                ║   • protein prediction with prodigal-gv       ║
                ╚═══════════════════════╤═══════════════════════╝
                                        │ proteins + contigs
                ╔═══════════════════════▼═══════════════════════╗
                ║                   PHASE 1                     ║
                ║         EVE Seed Region Prediction            ║
                ║  ─────────────────────────────────────────── ║
                ║   • pyhmmer scan of marker HMMs               ║
                ║   • marker hit taxonomy validation            ║
                ║   • infer host signature                      ║
                ║   • seed region assembly                      ║
                ╚═══════════════════════╤═══════════════════════╝
                                        │ seed regions
                ╔═══════════════════════▼═══════════════════════╗
                ║                   PHASE 2                     ║
                ║            Boundary Refinement                ║
                ║  ─────────────────────────────────────────── ║
                ║   • batched gene-taxonomy search              ║
                ║   • host-aware boundary trimming              ║
                ║   • gene-anchored extension/contraction       ║
                ╚═══════════════════════╤═══════════════════════╝
                                        │ refined candidates
                ╔═══════════════════════▼═══════════════════════╗
                ║                   PHASE 3                     ║
                ║      Evidence Synthesis & Scoring             ║
                ║  ─────────────────────────────────────────── ║
                ║   marker · taxonomy · compositional evidence  ║
                ║   ┌─────────────────────────────────────┐    ║
                ║   │  Optional annotation layers         │    ║
                ║   │   ◦ TMVec  (structural similarity)  │    ║
                ║   │   ◦ InterProScan  (domains)         │    ║
                ║   │   ◦ Boltz / Foldseek                │    ║
                ║   └─────────────────────────────────────┘    ║
                ║   • confidence tiering (HIGH/MEDIUM/LOW)      ║
                ║   • v2 acceptance quality gate                ║
                ╚═══════════════════════╤═══════════════════════╝
                                        │
              ┌─────────────────────────┼─────────────────────────┐
              ▼                         ▼                         ▼
   ┌────────────────────┐    ┌────────────────────┐    ┌────────────────────┐
   │ predictions.tsv    │    │   eves.fna (FASTA) │    │ summary.json /     │
   │ predictions.gff3   │    │   bed / detailed   │    │ batch_summary.tsv  │
   └────────────────────┘    └────────────────────┘    └────────────────────┘
```

1. Phase 0: optional repeat masking and protein prediction with `prodigal-gv`.
2. Phase 1: pyhmmer marker discovery plus marker-database validation to assemble seed regions.
3. Phase 2: gene-anchored boundary refinement using batched gene-taxonomy searches and host-aware trimming.
4. Phase 3: confidence scoring and report generation using marker, taxonomy, compositional, and optional structural/domain evidence.

The pipeline runs the HMM-gated workflow as its only execution path. Legacy novelty-heavy whole-proteome seeding (Path B) has been removed.
Phase 2 requires the large `gene_taxonomy_faa_db` resource; the smaller marker validation database is not accepted as a boundary-taxonomy fallback.
An opt-in prototype searches the Phase 0 proteome once and reuses the result
for both Phase 2 taxonomy consumers. It is off by default because query shape
can affect DIAMOND heuristics. Run paired prediction and boundary checks before
using it for an analysis.

## Installation and Resources

Install the project environment from the repository root:

```bash
pixi install --locked
```

Provision the pinned core resources:

```bash
pixi run setup-virosync-resources
```

Verify the installed bundle:

```bash
pixi run virosync orchestrate resources verify \
  --config config/orchestration.yaml --full
```

Expected:

```text
Version: v1.0.6
```

The public release provisions a pinned prebuilt resource bundle. Building the
full ViroSync database from raw source inputs is not part of the public setup
workflow yet.

The default resource bundle is configured in [config/orchestration.yaml](config/orchestration.yaml)
and targets `resources_v1_0_6.tar.gz`. The config pins both the complete archive
and its internal manifest by SHA-256. Version 1.0.6 packages the v1.0.5 biological
resources with authenticated payload metadata and transactional installation.
The bundle layout, marker annotation table, and resource release checks are
documented in [docs/RESOURCE_BUNDLE.md](docs/RESOURCE_BUNDLE.md).

The default install root is `resources/virosync`. Override it with either:

```bash
VIROSYNC_DB_ROOT=/path/to/virosync-db pixi run setup-virosync-resources
```

or:

```bash
pixi run virosync orchestrate setup \
  --config config/orchestration.yaml \
  --db-root /path/to/virosync-db
```

Setup validates the archive in a sibling staging directory, then switches the
stable path to a relative symlink backed by a read-only versioned directory. It
retains the previous working resource tree for recovery.

Benchmarking and manuscript assets remain private until the preprint and are
not part of the public release surface.

## Quick Start

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

Expected smoke-test result:

```bash
cat results/example_16t/batch_summary.tsv
```

The `test-1` row should report `status=success`, `predictions=6`, and
`accepted=1` with the v1.0.6 resource bundle.

Run a single genome:

```bash
pixi run virosync \
  -i genome.fasta \
  -o results/run_name \
  --config config/orchestration.yaml \
  -w 1 \
  --threads-per-worker 16
```

Run a batch directory:

```bash
pixi run virosync \
  -i data/genomes/ \
  -o results/batch_run \
  --config config/orchestration.yaml \
  -w 4 \
  --threads-per-worker 8
```

`pixi run example` remains available as a convenience wrapper around the shipped example input and repository defaults.

### Console output

Setup and pipeline runs print the ViroSync banner, software version, database
version, and one aggregate progress bar. Interactive terminals update
the bar in place. Redirected output records progress at bounded milestones.

Pass `--verbose` to the setup or run command to print configuration, phase
details, and diagnostic logs:

```bash
pixi run virosync run --verbose \
  -i genome.fasta \
  -o results/run_name \
  --config config/orchestration.yaml
```

Pass the root `--quiet` flag before the command to suppress the banner,
progress, and final summary while retaining errors:

```bash
pixi run virosync --quiet run \
  -i genome.fasta \
  -o results/run_name \
  --config config/orchestration.yaml
```

Fresh setup prompts and the download-size notice remain visible in quiet mode
because they require an explicit install-location and download decision.

Pixi can print task and cache messages before ViroSync starts. ViroSync
verbosity flags do not control those messages.

## Optional Annotation Layers

TMVec, InterProScan, and Boltz/Foldseek are optional and disabled unless explicitly enabled. The default `pixi install` is CPU-oriented and does not install the GPU/TMVec runtime or Boltz.

Behavior summary:

- TMVec runs only when `--tmvec` or config `phase3.use_tmvec_database: true` is set. It needs the optional Pixi structural/GPU environment, CUDA, the original trained TM-Vec weights, ProtT5, and precomputed TMVec embedding databases from the same model family.
- `--tmvec-gpu` makes TMVec fail fast if CUDA or embedding fails; plain `--tmvec` disables the layer with warnings when runtime/assets are unavailable.
- Boltz runs only when `--boltz` or config `phase3.use_boltz: true` is set. It uses the isolated `tools/boltz_runtime` Pixi environment plus a Foldseek database prefix. The current ViroSync Boltz integration writes sequence-only inputs, so set `phase3.boltz_use_msa_server: true` when enabling Boltz.
- InterProScan runs only when enabled and `phase3.interproscan_dir` points at an installation containing `interproscan.sh`.

Install the optional runtimes you need:

```bash
pixi install --locked -e structural
pixi install --locked --manifest-path tools/boltz_runtime/pixi.toml
```

Supported TMVec database keys are `bfvd`, `cath`, and `swissprot` for the pipeline confidence path. `bfvd` is the viral database that can increase `structural_score`; `cath` and `swissprot` are background/reference searches reported in `virosync_tmvec_proteins.tsv` as specificity context. The loader also recognizes a `pdb` embedding layout for exploratory searches, but the current EVE confidence model does not use PDB hits. The Hugging Face `scikit-bio/tmvec-2` checkpoint is not used here because it expects 408-dimensional input features; ViroSync's supported TMVec databases use 1024-dimensional ProtT5 residue embeddings.

Expected TMVec layouts under `phase3.tmvec_database_dir`:

```text
tmvec/
|-- bfvd/bfvd_embeddings.npy + bfvd_annotations.npy
|-- cath/cath_large.npy + cath_large_metadata.npy
`-- swissprot/swiss_large.npy + swiss_large_metadata.npy
```

Flat layouts are also accepted for one database at a time, such as `bfvd_embeddings.npy` next to `bfvd_annotations.npy`.
Custom BFVD embeddings generated by older ViroSync builds that logged `Using TMvec with default ProtT5 configuration` were produced by an untrained fallback model and are rejected by the TMVec loader/preflight; regenerate those embeddings with the trained `tmvec_swiss_model_large` checkpoint before using BFVD structural confidence.

Configure explicit optional-resource paths during setup:

```bash
pixi run virosync orchestrate setup \
  --config config/orchestration.yaml \
  --no-interactive-optional \
  --tmvec-dir /path/to/tmvec \
  --interproscan-dir /path/to/interproscan \
  --boltz-db-dir /path/to/bfvd
```

Validate the optional runtime after setup:

```bash
pixi run -e structural check-structural-runtime --require-tmvec
pixi run -e structural check-structural-runtime --require-boltz
pixi run -e structural check-structural-runtime --require-interproscan
```

The TMVec preflight includes a small model-compatibility smoke test. A runtime
that imports successfully but cannot produce a nonzero TMVec embedding is
treated as unavailable and does not contribute confidence.

The Boltz preflight fails when `phase3.boltz_use_msa_server` is false because
ViroSync does not yet write local MSA files into Boltz inputs.

Run with optional structural and domain annotation:

```bash
pixi run -e structural virosync \
  -i genome.fasta \
  -o results/structural_run \
  --config config/orchestration.yaml \
  --device cuda \
  --tmvec \
  --interproscan \
  --boltz \
  --no-skip-structural
```

If optional assets are missing, ViroSync prints a warning and continues with the
corresponding layer disabled.

Structural layers do not create new seed regions. They run during Phase 3 on candidates that already passed Phase 1/2, and can only annotate candidates or increase confidence through the structural component of the final score. The canonical output gate still applies after scoring.

## Output Layout

ViroSync writes one subdirectory per genome beneath the batch output root.

```text
results/<run_name>/
├── batch_summary.tsv
├── batch_report.md
└── <genome_id>/
    ├── phase0/
    ├── phase1/
    ├── phase2/
    ├── phase3/
    ├── phase0.complete.json
    ├── phase1.complete.json
    ├── phase2.complete.json
    ├── phase3.complete.json
    ├── phase3_synthesis/
    │   ├── virosync_predictions.tsv
    │   ├── virosync_predictions_detailed.tsv
    │   ├── virosync_predictions.bed
    │   ├── virosync_predictions.gff3
    │   ├── virosync_summary.json
    │   ├── evidence_profiles.json
    │   ├── interproscan_summary.tsv
    │   ├── virosync_tmvec_proteins.tsv
    │   └── virosync_jelly_roll_proteins.tsv
    ├── virosync_predictions_detailed.tsv
    ├── virosync_tsv_invariant_report.tsv
    ├── <genome_id>_eves.fna
    ├── run.log
    ├── virosync_run_complete.json
    ├── virosync_run_state.json
    └── notebooks/
        └── jupyter/eve_analysis.ipynb
```

Key output semantics:

| Path | Meaning |
|------|---------|
| `phase3_synthesis/virosync_predictions.tsv` | Canonical table of accepted predictions |
| `phase3_synthesis/virosync_predictions.bed` | Canonical 0-based half-open coordinates for accepted predictions |
| `phase3_synthesis/virosync_predictions.gff3` | Canonical GFF3 annotations for accepted predictions |
| `phase3_synthesis/virosync_predictions_detailed.tsv` | Detailed table of all Phase 3 candidates, including rejected candidates |
| `<genome_id>/virosync_predictions_detailed.tsv` | Convenience copy of the detailed table at the run root |
| `<genome_id>/<genome_id>_eves.fna` | Combined FASTA for accepted predictions only |
| `<genome_id>/virosync_tsv_invariant_report.tsv` | QA report for detailed TSV invariants |
| `<genome_id>/virosync_run_complete.json` | Human-readable final completion metadata; validated as a recorded output |
| `<genome_id>/virosync_run_state.json` | Authoritative schema-v3 run identity, status, result counts, and final artifact identities |
| `<genome_id>/phase<N>.complete.json` | Ordered phase record with dependency and artifact identities |
| `batch_summary.tsv` | Per-genome status and count summary for the batch run |

The output schema version is 4. `effective_eve_class` is one of `NCLDV`,
`MIRUS`, `PPV`, `CRESS`, `MIXED`, or `UNKNOWN`, and each accepted prediction
contributes to one class count. `PPV` contains the VP and PLV subtypes. The
detailed table writes `ppv_subtype=VP` or `ppv_subtype=PLV` only when
subtype-specific marker evidence supports one subtype and not the other.
Ambiguous PPV candidates retain `effective_eve_class=PPV` and use `.` for
`ppv_subtype`.

The detailed table groups columns by final call, candidate provenance, marker
evidence, gene taxonomy, composition and host evidence, InterProScan evidence,
and marker-set completeness. Its `taxonomy_best_hits` partition is:

```text
EUK;MITO;PLASTID;BAC;ARC;UNK;NO_HITS;NCLDV;MIRUS;PPV;CRESS;GVMAG;PHAGE
```

The four `*_top10_proteins` columns count raw top-10 prefix support. The
mutually exclusive partition assigns a viral family only when the supporting
hit reaches 25% amino-acid identity. When markers do not assign another family,
identity-qualified gene taxonomy can assign CRESS. `vp_completeness` records
VP-subtype evidence;
`ppv_completeness` combines the PPV marker sets.

The result reader maps `VP` and `PLV` class aliases to `PPV`. Batch summaries
report the parent classes only:
`ncldv`, `mirus`, `ppv`, `cress`, `mixed`, and `unknown`.

Normal runs validate schema-v3 state before reuse. The run fingerprint binds the input,
effective output-determining configuration, source code, locked runtime, enabled tools
and models, masking request, and resource manifests. Phase markers are checked in order;
each recorded artifact must retain its relative path, size, SHA-256 digest, schema, and
row count where applicable. Validation stops at the first stale phase and recomputes
that phase and everything downstream. Schema-v1/v2 outputs and unmarked partial files
are never resume evidence. `--clean-run` clears the genome output before the first
attempt; an automatic retry validates and resumes the surviving schema-v3 phase prefix.

Each successful genome run writes the executed analysis notebook `notebooks/jupyter/eve_analysis.ipynb` (rendered from the jupytext source `src/virosync/report/eve_analysis.py`).
Notebook execution also emits rendered summary figures (for example confidence, marker, and gene-category PNGs) at the genome root.

## Confidence Tiers and Quality Gates

ViroSync assigns every candidate one of three confidence tiers:

| Tier | Score | Meaning |
|------|-------|---------|
| HIGH | ≥ 0.7 | Strong multi-evidence support |
| MEDIUM | 0.2 – 0.7 | Moderate evidence; accepted for downstream review |
| LOW | < 0.2 | Weak evidence; qualified calls promoted by the v2 acceptance quality gate (see below), remainder retained in detailed TSV only |

### Priority-marker promotion

Candidates with a Major Capsid Protein (MCP) hallmark are promoted from LOW to MEDIUM
regardless of their composite score, because MCP is a family-defining marker with high
specificity.

### Non-HHG seeding quality gate (defensive)

The active pipeline seeds candidates via the HMM-gated hallmark discovery (HHG)
path only; the legacy k-mer compositional and taxonomic novelty seeding paths
have been removed. A defensive quality gate remains in Phase 3 so that any
future candidate reaching scoring without an `hhg` entry in `seed_sources` is
demoted to LOW unless it meets at least one of:

- `hallmark_count >= 3`: three or more Phase-3-validated hallmark genes
- `hallmark_count >= 1 AND has_mcp`: at least one hallmark and it is an MCP
- `non-host gene count >= 5`: five or more interior genes without host-like taxonomy

The seeding provenance for each candidate is recorded in the `seed_sources`
column of `virosync_predictions_detailed.tsv` (typical values under the
current pipeline: `hhg`, `marker_validation`).

### EVE acceptance quality gate (v2)

The pipeline applies the v2 acceptance quality gate when generating the canonical
EVE set for downstream analysis. The EVE class is resolved from
`region_classification`, falling back to the family-like column (`likely_family`
in legacy outputs, `classification` in current outputs) when the region label is
not a concrete family. A concrete single-family label always takes precedence; a
region whose seed markers span more than one viral family resolves to `MIXED` and
is gated as a viral region under the standard rule rather than discarded.

**HIGH & MEDIUM confidence:**

| EVE class | Length | Marker requirement |
|-----------|--------|-------------------|
| PPV, CRESS, MIXED | > 2 kb | has MCP OR (hallmark >= 2 with >= 1 non-ATPase) |
| NCLDV, MIRUS | > 5 kb | None |
| NCLDV, MIRUS | <= 5 kb | has MCP |

PPV (Preplasmiviricota) contains Polinton-like viruses (PLV) and virophages
(VP). The v1.0.6 database uses the parent `PPV__` taxonomy label. ViroSync
assigns a VP or PLV subtype only from unambiguous subtype-specific marker
evidence.

**LOW confidence (promoted if):**

| EVE class | Length | Marker requirement |
|---------------|--------|-------------------|
| NCLDV, MIRUS | > 5 kb | >= 2 validated hallmark markers |
| PPV, CRESS, MIXED | > 2 kb | has MCP OR (hallmark >= 2 with >= 1 non-ATPase) |

Other LOW-confidence calls are not exported. This design recovers degraded
EVEs in heavily endogenized genomes (e.g. *Chlamydomonas* sp. ICE-L) that
have NCLDV hallmarks but are classified as UNKNOWN due to host gene
accumulation around the viral core.

## Runtime Behavior

- CLI flags override YAML config values; when a flag is not set explicitly, values from [config/orchestration.yaml](config/orchestration.yaml) are used.
- A completed genome is reused only when its schema-v3 run state, ordered phase
  chain, run identity, result contract, and every recorded artifact validate.
  Use `--clean-run` to force a fresh rerun.
- Optional TMVec, InterProScan, and Boltz layers are skipped cleanly when they are not enabled or not provisioned.
- For single-genome-slot execution, `-w 1 --threads-per-worker <N>` is the most predictable mode on shared systems.

## CLI Entrypoints

| Command | Purpose |
|---------|---------|
| `pixi run virosync -i <input> -o <output>` | Run the local Python-parallel pipeline |
| `pixi run virosync orchestrate setup` | Install and validate ViroSync resources |
| `pixi run virosync info` | Print system and runtime information |
| `pixi run -e structural check-structural-runtime` | Validate optional structural/runtime dependencies |

## Documentation

| File | Purpose |
|------|---------|
| [README.md](README.md) | Operator guide and run instructions |
| [docs/METHODS.md](docs/METHODS.md) | Code-verified workflow and output specification |
| [config/orchestration.yaml](config/orchestration.yaml) | Repository-default orchestration and resource configuration |
| [CHANGELOG.md](CHANGELOG.md) | Release and unreleased change log |

## Testing

Run the tracked test suite from the repository root:

```bash
pixi install --locked
pixi run lint
pixi run python -m pytest
pixi run python -m compileall -q src tests scripts
```

Run the production release-surface guard:

```bash
pixi run check-production-ready
```

GitHub Actions run three public-readiness checks:

- `tests`: installs the locked Pixi environment, runs the same pinned Ruff task
  used locally, executes pytest and compileall, checks release surfaces, and
  verifies the CLI version.
- `production-guards`: typed-parses both shipped configurations, cross-checks
  software and resource identities, authenticates the tracked v1.0.6 manifest,
  and checks the public archive URL on schedules, manual runs, and release tags.
- `example-smoke`: runs weekly, manually from the default branch, and for
  release tags. It fully verifies provisioned resources, runs a clean example
  and unchanged resume, checks all authenticated completion artifacts and
  coordinates, and uploads path-safe summaries.

## Reproducibility Notes

- Always pass `--config config/orchestration.yaml` in scripted runs.
- Keep run outputs outside the repository when processing real datasets.
- Check `resources/virosync/DB_VERSION` to confirm the installed database bundle.
- Use [docs/METHODS.md](docs/METHODS.md) for a code-verified description of the current workflow rather than older design notes.
