# Changelog

All notable changes to ViroSync are documented in this file.

The format is based on Keep a Changelog and this project follows Semantic Versioning.

## [Unreleased]

### Added

- Opt-in Phase 1 frameshift-sensitive screening of the masked nucleotide
  assembly with BATH and the 808 shipped VS marker HMMs. ViroSync validates
  event-bearing aligned domains against the Tier-1 protein database. Confirmed,
  coverage-filtered loci can seed regions and appear in an accepted EVE FAA
  without entering the ordinary proteome. Enable the screen with
  `--frameshift-screening` or the run config.
- A tracked three-contig *Trichomonas vaginalis* G3 frameshift example, a Pixi
  task that runs it through the full pipeline, and CI checks for confirmed
  rescue markers and accepted rescued-marker FAA output.
- Compact setup and run output with software and database versions plus
  aggregate progress for resources and genome queries.
- Command-local `--verbose` flags for setup and run diagnostics.
- CRESS calls and counts in prediction, batch, report, and notebook outputs.
- Published EVE taxonomy class from a weighted vote over the element's genes:
  validated marker the all-gene search agrees with 3, validated marker alone 2,
  gene without a marker 1. A lineage needs strictly more than half the total
  weight. An MCP marker that cast a vote overrides that vote outright; MCP
  markers that disagree settle it by weight, 5 each.
- Per-genome ANI clustering of accepted EVEs at 95% identity and 50% aligned
  fraction, on by default. Members inherit the lineage class their MCP-decided
  relatives agree on.
- `phase3_synthesis/eve_ani_edges.tsv`, and detailed TSV columns
  `ani_cluster_id`, `ani_cluster_size`, `ani_max_percent`,
  `taxonomy_class_before_ani`, and `taxonomy_class_propagated_from`.
- `eve_ani_network.png` in the report notebook: one node per EVE colored by
  published class, MCP-decided nodes outlined, an edge per clustered pair,
  laid out with graphviz `sfdp`.

### Changed

- Rescue-derived and ordinary candidates remain separate through Phase 2. A
  rescue-only candidate must retain a confirmed rescued marker. Phase 3
  resolves direct ordinary/rescue overlaps without suppressing candidates that
  overlap only through another branch.
- Output schema v6 adds `canonical_selection_outcome` to the detailed TSV and
  uses v6 artifact identities for both prediction tables. The new output and
  run-fingerprint identities invalidate resume from older runs.
- Phase 3 reports the quality-gate drop count without subtracting the later
  unsupported-viral-evidence filter twice.
- Output schema v4 groups detailed TSV fields by evidence type.
- Output schema v4 to v5. `effective_eve_class` is `NCLDV`, `MIRUS`, `PPV`,
  `CRESS`, `PHAGE`, `VIRAL_UNKNOWN`, or `UNKNOWN`. `MIXED` is retired and reads
  back as `VIRAL_UNKNOWN`. Batch summaries replace the `mixed` column with
  `viral_unknown` and add `phage`.
- The v2 quality gate keeps its own class vocabulary and its own resolver, and
  still decides which candidates it accepts. The published class only labels
  them.
- Phase 3 now drops an accepted EVE that carries no validated marker, returned
  no identity-qualified viral gene hit, and shares no ANI cluster with a
  marker-bearing EVE. Nothing in such a region is viral. This is the one place
  the published set is smaller than the gate's own.
- The report notebook annotates gene-map and marker-heatmap clusters at 95%
  ANI rather than 90%, which regroups the EVEs in both figures. It reads the
  clusters from the pipeline instead of running skani itself.
- `graphviz` and `python-graphviz` are pinned in `pixi.toml`. The locked
  runtime is part of the run fingerprint, so runs made against the previous
  `pixi.lock` no longer resume; they restart from Phase 0.
- Resource setup captures downloader progress in the single ViroSync bar and
  shows redacted network errors when a download fails.
- Batch summaries report PPV as the parent class instead of separate VP and
  PLV peers.
- Detailed predictions record `ppv_subtype` only when marker evidence supports
  VP or PLV without conflicting subtype evidence.

## [1.0.0] - 2026-07-28

First public release.

ViroSync detects endogenous viral elements of giant viruses and their relatives
in assembled eukaryotic genomes. The workflow runs in four phases: optional
repeat masking and gene prediction; HMM marker discovery with a taxonomy gate
that validates markers before whole-proteome expansion; gene-anchored boundary
refinement with host-aware trimming; and multi-evidence confidence scoring that
assigns each accepted region a HIGH, MEDIUM or LOW tier.

Marker families cover Nucleocytoviricota, Mirusviricota, Preplasmiviricota
(Polinton-like viruses and virophages) and CRESS/ssDNA lineages, using the
v1.0.6 core resource bundle.

Outputs are a confidence-tiered prediction table, BED, GFF3 and FASTA of
accepted regions, per-region gene taxonomy, and a run manifest. Runs are
resumable and skip completed genomes.
