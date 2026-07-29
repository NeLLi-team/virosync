# Changelog

All notable changes to ViroSync are documented in this file.

The format is based on Keep a Changelog and this project follows Semantic Versioning.

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
