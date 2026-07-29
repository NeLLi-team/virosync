# Changelog

All notable changes to ViroSync are documented in this file.

The format is based on Keep a Changelog and this project follows Semantic Versioning.

## [Unreleased]

### Added

- Compact setup and run output with software and database versions plus
  aggregate progress for resources and genome queries.
- Command-local `--verbose` flags for setup and run diagnostics.
- CRESS calls and counts in prediction, batch, report, and notebook outputs.

### Changed

- Output schema v4 groups detailed TSV fields by evidence type.
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
