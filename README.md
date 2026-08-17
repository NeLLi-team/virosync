# ViroSync

![Version](https://img.shields.io/badge/version-1.0.0-blue)
[![License: non-commercial use only](https://img.shields.io/badge/license-non--commercial-orange.svg)](LICENSE)

ViroSync detects candidate endogenous viral elements (EVEs) in assembled eukaryotic genomes. The pipeline combines hallmark-marker discovery, taxonomy-guided boundary refinement, host-aware trimming, and multi-evidence confidence scoring.

## Workflow

ViroSync processes each genome in four phases.

![ViroSync workflow with Pfam marker validation and optional frameshift-sensitive marker rescue](docs/virosync_workflow.png)

1. Phase 0: optional repeat masking and protein prediction with `prodigal-gv`.
2. Phase 1: pyhmmer marker discovery combined with Pfam domain identification, and marker-taxxonomy validation to assemble seed regions.
3. Phase 2: gene-anchored boundary refinement using batched gene-taxonomy searches and host-aware trimming.
4. Phase 3: confidence scoring and report generation using marker, taxonomy, compositional, and optional structural/domain evidence.

## Installation and Resources

Clone the repository then install it from the repository root:

```bash
pixi install --locked
```

The optional Phase 1 frameshift screen requires `bathconvert` and `bathsearch`.
Pixi does not install BATH. To use the screen, build the pinned revisions from
source as described in
[the frameshift screening guide](docs/FRAMESHIFT_SCREENING.md). Add both commands
to `PATH`, then verify them:

```bash
bathconvert -h
bathsearch -h
```

Enable the screen with `--frameshift-screening` or set
`phase1.frameshift_screening_enabled: true` in the run config.
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
Version: v1.0.7
```
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

## Quick Start

Run the standard example.

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
`accepted=1` with the v1.0.7 resource bundle.

Run the three-contig *Trichomonas vaginalis* G3 frameshift example after
installing BATH:

```bash
pixi run example-frameshift
```

This command enables `--frameshift-screening` and writes
`results/example-frameshift/batch_summary.tsv`. With the v1.0.7 resource
bundle, the `trichomonas-g3` row should report `status=success`,
`predictions=5`, and `accepted=2`.

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

`pixi run example` runs the standard example.
`pixi run example-frameshift` runs only the nested G3 fixture with the screen
enabled.

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

### EVE taxonomy class

`effective_eve_class` comes from a weighted vote over the region's genes. Every
gene with an identity-qualified viral hit votes with its top-10 DIAMOND
taxonomy:

- A hit below 25% amino-acid identity does not qualify.
- A `GVMAG` hit counts as `NCLDV`. A `PHAGE`-namespace hit counts as `PPV` when
  its resolved lineage is Preplasmiviricota, and as `PHAGE` otherwise.
- A gene whose qualified hits give one class votes for that class. A gene whose
  hits span several classes votes `VIRAL_UNKNOWN`, which carries weight but
  never wins.
  
A marker-bearing gene is searched twice, against the marker reference in
Phase 1 and again in the Phase 2b all-gene search. Weights:

| Gene | Weight |
|------|--------|
| Validated MCP marker | 5 on the marker call |
| Validated marker, all-gene search agrees | 3 |
| Validated marker, all-gene search disagrees | 2 on the marker call, 1 on the conflicting all-gene call |
| No marker | 1 on the all-gene call |

A lineage class (`NCLDV`, `MIRUS`, `PPV`, `CRESS`, or `PHAGE`) needs strictly
more than half the total weight. Half is not enough, so two genes that disagree
leave the region `VIRAL_UNKNOWN`.

A region with no viral vote at all is `VIRAL_UNKNOWN` when it carries a
validated marker and `UNKNOWN` when it carries none.

An `UNKNOWN` region has no viral evidence of any kind: no validated marker, and
no gene of its own returned a qualified viral hit. Phase 3 drops it from the
published set unless it shares an ANI cluster with a marker-bearing EVE, which
makes it a decayed copy of a real element rather than host sequence the length
rule admitted. This is the only case where the published set is smaller than the
set the acceptance gate kept.

Phase 3 then compares every accepted EVE in the genome against every other with
skani. Two EVEs join when they reach 95% average nucleotide identity over at
least 50% of either sequence, and each connected component is one cluster.
A member may donate its class only when an MCP marker's own vote is the class
that won, which is narrower than carrying an MCP. Where the donors in a cluster
agree on a single lineage class, the other members take it. Nothing propagates
when the donors disagree, or when the class they agree on is `VIRAL_UNKNOWN`.
Every pair skani reports is written to `phase3_synthesis/eve_ani_edges.tsv`,
and the detailed table records the result per region:

| Column | Meaning |
|--------|---------|
| `ani_cluster_id` | Cluster index, `.` for an EVE with no clustered relative |
| `ani_cluster_size` | Members in the cluster, `1` for a singleton |
| `ani_max_percent` | Highest ANI to another member, `.` for a singleton |
| `taxonomy_class_before_ani` | Class the region's own vote gave, filled only where propagation replaced it |
| `taxonomy_class_propagated_from` | EVE ID the class came from |

Genomes with fewer than two accepted EVEs, and genomes whose accepted regions
skani cannot sketch, get a header-only edge table and no propagation.

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

The standard pipeline seeds candidates through HMM-gated hallmark discovery
(HHG). The opt-in frameshift screen can add a rescue-only seed without an
`hhg` source. The legacy k-mer compositional and taxonomic novelty paths have
been removed. A defensive quality gate in Phase 3 demotes a candidate without
an `hhg` source to LOW unless it meets at least one of:

- `hallmark_count >= 3`: three or more Phase-3-validated hallmark genes
- `hallmark_count >= 1 AND has_mcp`: at least one hallmark and it is an MCP
- `non-host gene count >= 5`: five or more interior genes without host-like taxonomy

The `seed_sources` column records provenance. Ordinary seeds carry `hhg` and
`marker_validation`, rescue-only seeds carry `frameshift_rescue`, and a region
with both marker types carries all three.

### EVE acceptance quality gate (v2)

The pipeline applies the v2 acceptance quality gate when generating the canonical
EVE set for downstream analysis. The gate decides which candidates are accepted.
It does not set the published class. Its own class vocabulary is `NCLDV`,
`MIRUS`, `PPV`, `CRESS`, `MIXED`, and `UNKNOWN`, and the `MIXED` in the tables
below is that internal label, not an output token.

The gate resolves its class from `region_classification`, falling back to the
family-like column (`likely_family` in legacy outputs, `classification` in
current outputs) when the region label is not a concrete family. A concrete
single-family label always takes precedence; a region whose seed markers span
more than one viral family resolves to `MIXED` and is gated as a viral region
under the standard rule rather than discarded.

**HIGH & MEDIUM confidence:**

| EVE class | Length | Marker requirement |
|-----------|--------|-------------------|
| PPV, CRESS, MIXED | > 2 kb | has MCP OR (hallmark >= 2 with >= 1 non-ATPase) |
| NCLDV, MIRUS | > 5 kb | None |
| NCLDV, MIRUS | <= 5 kb | has MCP |

PPV (Preplasmiviricota) contains Polinton-like viruses (PLV) and virophages
(VP). The v1.0.7 database uses the parent `PPV__` taxonomy label. ViroSync
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
| [docs/FRAMESHIFT_SCREENING.md](docs/FRAMESHIFT_SCREENING.md) | Run and interpret frameshift-sensitive VS marker screening |
| [config/orchestration.yaml](config/orchestration.yaml) | Repository-default orchestration and resource configuration |
| [CHANGELOG.md](CHANGELOG.md) | Release and unreleased change log |

- Check `resources/virosync/DB_VERSION` to confirm the installed database bundle.
- Use [docs/METHODS.md](docs/METHODS.md) for a code-verified description of the current workflow rather than older design notes.
