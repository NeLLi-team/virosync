# ViroSync Methods (Code-Verified, July 29, 2026)


## Workflow summary

ViroSync detects candidate endogenous viral elements (EVEs) in
assembled eukaryotic genomes with a four-phase workflow:

1.  Phase 0: optional repeat masking and protein-coding gene prediction.
2.  Phase 1: HMM marker discovery and taxonomy-gated marker validation
    to build seed regions.
3.  Phase 2: gene-anchored boundary refinement using batched taxonomy
    searches plus host-aware trimming.
4.  Phase 3: multi-evidence scoring, confidence tiering, and report
    generation.

This document reflects current repository behavior, not older design
notes.

## Scope and execution path

The pipeline runs the HMM-gated workflow as its only execution path,
using the Python genome-parallel workflow runner introduced in ViroSync
1.4. Novelty-heavy whole-proteome seeding (legacy Path B) has been
removed. Phase 2 boundaries are driven by gene extension and taxonomy
trimming; CRF-based expansion is kept in code but is not the default
boundary driver.

Primary command entrypoint:

``` bash
pixi run virosync -i <input_fasta_or_dir> -o <output_dir> --config config/orchestration.yaml
```

Default setup and run output consists of the ViroSync banner, software and
database versions, an aggregate progress bar, and a final run summary.
`--verbose` on either command enables configuration, phase, and diagnostic
output. The root `--quiet` flag suppresses the banner, progress, and final
summary but does not suppress errors. Fresh setup prompts and the download-size
notice remain visible because setup requires explicit location and download
confirmation. Pixi task and cache messages are outside the ViroSync output
controls.

## Inputs, resources, and environment

Required inputs are assembled genome FASTA files. ViroSync runs per
genome and also supports batch directory input.

Dependencies and tools are pinned in `pixi.toml` (for example: Python
3.11.14, pyhmmer 0.11.3, DIAMOND 2.2.1, and prodigal-gv 2.11.0). Core
resources are configured through `config/orchestration.yaml`; the current
default points to the SHA-256-pinned `resources_v1_0_6.tar.gz`. It contains
the v1.0.5 CRESS/ssDNA and MetaVR-derived marker/taxonomy resources plus an
authenticated schema-v1 manifest. The resource builder can also emit a bound
schema-v2 runtime artifact and source/repair artifact. The current public pin
remains schema v1. Tool citations for these components are listed in “Tool
citation map” below.

Resource setup:

``` bash
pixi install --locked
pixi run setup-virosync-resources
pixi run virosync orchestrate resources verify --full
```

The core resource bundle layout and marker annotation table checks are
documented in `docs/RESOURCE_BUNDLE.md`.

An installed core bundle is selected through its externally pinned version,
archive digest, and manifest digest. After a full installation, normal run
startup checks `DB_METADATA.json` before rehashing the large payloads. Receipt
reuse requires the exact relative versioned pointer, an immutable closed
inventory, regular single-link files, and matching device, inode, size,
modification-time, and change-time records for every manifest payload. A
current receipt returns without a payload scan. Any mismatch falls back to full
manifest validation. The receipt is a performance cache, not an independent
trust root.

The DIAMOND 2.2.1 pin was compared with 2.1.21 on six Phase 2 query sets using
the same v1.0.6 database, four warm repetitions, and position-balanced run
order. Output bytes matched for every query set. The geometric mean of the six
median 2.2.1-to-2.1.21 runtime ratios was 0.959, or about 4.1% faster.

## Authenticated run state and resume

Each genome output has one authoritative `virosync_run_state.json` with schema
version 3 and four ordered `phase<N>.complete.json` records. The run fingerprint
is canonical JSON over:

- input FASTA size and SHA-256;
- the output-determining effective configuration and coordinate/output schemas;
- installed ViroSync source and locked runtime identities;
- requested and effective execution environment, excluding thread count;
- enabled executable and immutable model identities;
- requested masking and authenticated core or optional resource identities.

Phase records form a hash chain. Each record binds the full run fingerprint,
the preceding marker digest, and the relative path, size, SHA-256, schema, and
row count of every required artifact. Phase 0 also binds the validated masking
result and the exact sequence passed to Prodigal. Phase 1 stores a lossless
`phase1/resume_state.json`; Phase 2 stores lossless refined and full resume
state plus an exact BED projection.

Resume validates the phase chain from Phase 0. At the first missing or stale
marker or artifact, that phase and every downstream phase are invalidated
through guarded no-follow deletion. Unmarked files, schema-v1/v2 state, and a
TSV header alone are never accepted. Final success is reusable only after the
canonical and detailed tables, BED, GFF3, summary, invariant report, completion
metadata, run log, notebook, class/tier counts, and input-coordinate bounds all
validate. Automatic retries set `resume=true` after the first attempt, including
when the initial command used `--clean-run`.

## Phase 0: preprocessing

Phase 0 creates the sequence context used by downstream phases.

- Gene calling: `prodigal-gv` is run in metagenomic mode. With multiple
  threads, scaffolds are split and processed in parallel, then merged
  into `proteome.fasta`.
- Optional masking: `execution.masking.backend` is `off` by default, which
  avoids removing potentially informative repetitive viral loci. Supported
  backends are `off`, `trf`, `repeatmasker`, and `trf_repeatmasker`.
  RepeatMasker backends require exactly one explicit
  `repeatmasker_species` or non-empty `repeatmasker_library`; the shipped
  configs do not assume a species. The default `strict` failure policy stops
  the run if a requested backend fails. A `fallback` policy must name either
  `off` or `trf`, and fallback results are excluded from primary benchmarks.

Phase 0 writes `phase0/masking/masking_status.json` for both disabled and
enabled masking. The status records the requested and effective backends,
verified tool versions, target or library digest, masked-base count, input and
output digests, fallback outcome, and benchmark eligibility. ViroSync verifies
that this status identifies the exact sequence passed to `prodigal-gv` and
binds its digest into the completion manifest. Outputs also include a proteome
FASTA and gene coordinates used for region assembly and boundary operations.

## Phase 1: HMM-gated seeding

Phase 1 narrows search space before expensive full taxonomy steps.

### 1) HMM scan

Predicted proteins are screened against the marker HMM database
(`models/combined.hmm` in the required v1.0.6 resource bundle) using
pyhmmer.

The HMM search wrapper accepts an optional reporting E-value cutoff.
When the cutoff is unset (`evalue_cutoff=None`), ViroSync passes
infinite pyhmmer reporting thresholds and records all reported HMM hits
for downstream validation. The June 2026 Python workflow-runner
benchmark rerun used this no-HMM-reporting-cutoff mode, followed by the
same marker-validation, boundary-refinement, and canonical output gates.

### 2) Marker validation with small taxonomy DB

Only HMM-hit proteins are searched against the marker validation
database (small Tier 1 target set). This reduces runtime versus
searching the whole proteome.

Current marker status logic:

- `validated`: at least one `NCLDV__`, `MIRUS__`, `PPV__` (Preplasmiviricota;
  legacy `VP__`/`PLV__` transitionally), or `CRESS__` top-10 marker hit with
  percent identity at least 25%.
- `validated_novel`: no Diamond hit but passes HMM-only gating
  (score/coverage/cluster criteria); in seed construction this path is
  effectively restricted to MCP-like markers.
- `supported`: partial viral support (for example GVMAG-associated)
  without full validation.
- `unvalidated`: predominantly cellular signal or below validation
  criteria.

### 3) Host signature model

From unvalidated host-like marker hits, ViroSync builds a weighted token
model of host taxonomy. This model is persisted and reused in later
host-trimming and penalty steps.

### 4) Region assembly from validated markers

Validated markers are clustered by distance (base pairs and gene count),
then iteratively extended by configurable flanks until no new markers
are captured. Overlapping regions are merged. These candidate regions
are converted into seeds for Phase 2.

## Phase 2: boundary refinement (taxonomy-first)

Phase 2 resolves seed boundaries with gene-level taxonomy context.

### 1) Seed-region taxonomy and host trimming

The default path searches genes that overlap the initial seed regions. It uses
those hits for host-aware trimming before gene extension. Set
`phase2.host_trim_enabled: false` in the config to skip the step and go straight
to gene extension.

### 2) Gene extension and boundary taxonomy

Each trimmed seed is extended by plus/minus genes (default: 5), and
overlapping extended seeds are merged. A second batched search covers all
seed, interior, flanking, and sampled control genes. Both searches use the
large gene-taxonomy database (Tier 2). The result is a per-gene taxonomy map
with top-hit and top-k support fields.

The opt-in `phase2.diamond_superset_prototype_enabled` path searches
the complete Phase 0 proteome once, before trimming, and slices that raw result
for both consumers. It requires `phase2.diamond_top_k: 10`, remains off by
default, and authenticates the raw TSV in the Phase 2 completion record.
DIAMOND uses the hit limit and query size in search heuristics, so cached
slicing is not a general proof of equivalence to the two default searches.
Evaluation must compare raw ordered hits, refined state and BED, taxonomy maps,
control statistics, host-trim tables, canonical and detailed predictions,
BED/GFF3, and summaries on paired clean runs. Compare per-region taxonomy rows
by gene ID. File order is not a downstream contract. The
prototype flag changes the run fingerprint, and authenticated resume reuses
the saved Phase 2 state without another search.

On the shipped example, paired clean DIAMOND 2.2.1 runs gave identical
ordered hits and passed 1,022 prediction, boundary, taxonomy, and provenance
checks. The prototype reduced Phase 2 from 61.3 to 40.9 seconds and total
per-genome time from 101.2 to 80.3 seconds. An authenticated resume completed
with zero reported genome time and left all 30 recorded artifacts unchanged.
The result covers one genome. The flag remains off by default.

### 3) Taxonomy-based seed refinement

Refinement can run with:

- heuristic taxonomy walk (default), or
- optional ML-assisted taxonomy refiner when enabled.

### 4) Host-aware trimming

Boundaries are trimmed using:

- host baseline fingerprints from control genes,
- per-gene top-k taxonomy signal,
- the Phase 1 host-signature model,
- local density logic for ambiguous/no-hit neighborhoods.

After trimming, boundaries are hard-constrained to the seed-specific
flanking envelope so every retained gene remains in the precomputed
taxonomy support range.

### 5) Adjacent-boundary merge

Adjacent EVEs can be merged when short inter-boundary gaps show enough
viral taxonomy support (or strong flanking viral context).

Note: current default Phase 2 path converts refined seeds directly to
`RefinedBoundary` objects and does not perform CRF-driven boundary
expansion.

## Phase 3: evidence synthesis and confidence assignment

Phase 3 combines multiple evidence streams per candidate.

Evidence includes:

- marker content and marker-family consistency,
- gene taxonomy composition (interior and flanking tracked separately),
- compositional deviation (for example KFD/GC/CUB-derived features),
- host-signature burden (fraction of host-like interior genes),
- optional InterProScan keyword/category support,
- optional TMVec structural similarity (database search),
- optional Boltz + Foldseek structural tie-breaker,
- optional phylogenetic validation (GVClass + Diamond path when
  enabled),
- coherence analysis from evidence graph features.

Structural/domain layers are candidate annotators and confidence
modifiers, not seed generators. TMVec, Boltz/Foldseek, and InterProScan
run only after Phase 1/2 have produced candidate boundaries. They can
increase confidence for those candidates, but cannot create new EVE
regions by themselves, and the final canonical output gate is still
applied after scoring.

### Confidence model

The final score is computed from weighted components plus bonuses minus
penalties:

- base weighted mixture of marker, gene taxonomy, composition, optional
  CRF contribution, InterPro, structural support, and seed-cluster
  evidence,
- additive bonuses for convergent evidence (for example marker+taxonomy
  synergy, MCP, completeness, family consistency, taxonomy divergence
  bonuses),
- penalties for high-confidence eukaryotic-only signal, elevated
  host-signature fraction, and very small non-MCP regions.

Protective rules include:

- priority-marker floors (default marker list contains MCP),
- low-confidence caps for candidates with little or no viral taxonomy
  evidence,
- optional phylogenetic rejection override behavior.

Tier mapping defaults:

- `HIGH`: score \>= 0.7
- `MEDIUM`: 0.2 \<= score \< 0.7
- `LOW`: score \< 0.2

Candidates with priority-marker evidence can be promoted from `LOW` to
`MEDIUM` under configured rules.

### Taxonomy class assignment

The published class of a region comes from a weighted vote over its genes. It
labels the region and does not feed the acceptance gate. The gate resolves its
own class from `region_classification` and the family-like columns, in a
separate vocabulary that still carries `MIXED`.

Every gene inside the refined boundary with an identity-qualified viral hit
votes with its top-10 taxonomy:

- a hit below 25% amino-acid identity does not qualify;
- reference namespaces map onto published classes: `VP` and `PLV` to `PPV`,
  `GVMAG` to `NCLDV`, and `PHAGE` to `PPV` when the target's resolved lineage
  contains Preplasmiviricota (the legacy `PHAGE__VARDNA__` records) and to
  `PHAGE` otherwise;
- one distinct class across a gene's qualified hits is that gene's vote;
  several distinct classes make the vote `VIRAL_UNKNOWN`, which carries weight
  but can never win; no qualified viral hit is no vote.

A gene without a vote is left out of the denominator. A validated marker with
no qualified hit is the HMM-only `validated_novel` case.

A marker-bearing gene is searched twice, against the marker validation database
in Phase 1 and again in the Phase 2b all-gene search. Weights:

- validated MCP marker: weight 5 on the marker call;
- validated marker, all-gene search agrees: weight 3;
- validated marker, all-gene search disagrees: weight 2 on the marker call and
  weight 1 on the conflicting all-gene call;
- no marker: weight 1 on the all-gene call.

A lineage class (`NCLDV`, `MIRUS`, `PPV`, `CRESS`, or `PHAGE`) needs strictly
more than half the total weight. Half is not enough, so two genes that disagree
leave the region `VIRAL_UNKNOWN`. A region with no viral vote at all is
`VIRAL_UNKNOWN` when it carries a validated marker and `UNKNOWN` when it
carries none.

The major capsid protein decides ahead of the weighted vote. An MCP marker that
cast a vote sets the class alone, whatever the other genes say, including when
its own top-10 spans several lineages and it votes `VIRAL_UNKNOWN`. An MCP with
no qualified viral hit casts no vote and decides nothing. MCP markers that
disagree fall through to the weighted vote at weight 5 each, which is the only
place the weights break an MCP tie.

Phase 3 drops an `UNKNOWN` region from the published set unless it shares an ANI
cluster with a marker-bearing EVE. Such a region holds no viral evidence at all,
so without a clustered relative to vouch for it the length rule admitted host
sequence. The drop runs after clustering, and it is the only step that publishes
fewer EVEs than the acceptance gate kept.

`ppv_subtype` comes from VP-specific and PLV-specific HMM markers, not from
taxonomy, and is reported only for regions published as `PPV`. The v1.0.6
database labels Preplasmiviricota references `PPV__` and holds no `VP__` or
`PLV__` records, so the top-10 hits cannot separate the two subtypes.

### ANI clustering and class propagation

Accepted regions of one genome are compared all against all with skani
(`triangle -E --medium -m 200 -s 80`, minimum aligned fraction 50), before any
class is counted or persisted. Two regions join when they reach 95% average
nucleotide identity over at least 50% of either sequence, and each connected
component is one cluster. Clusters are numbered by descending size, ties broken
by lowest member EVE ID, so a rerun of one genome reproduces the numbering.

A region may donate its class only when an MCP marker's own vote is the class
that won. That is narrower than carrying an MCP: a capsid annotation, a
structural jelly-roll call, and phylogenetic evidence also mark a region as
MCP-bearing without casting a taxonomy vote. In a cluster holding both donors
and non-donors whose donors agree on a single lineage class, every non-donor
takes that class and records the source EVE ID. Nothing propagates when the
donors disagree, or when the class they agree on is `VIRAL_UNKNOWN` or
`UNKNOWN`.
Clustering runs on the fixed accepted set and rewrites only the class. The
clustering confidence bonus stays 0.0, since scoring has already run.

Every pair skani reports is written to `phase3_synthesis/eve_ani_edges.tsv` as
`eve_a`, `eve_b`, `ani`, `af_a`, `af_b`, including pairs below the 95% ANI
threshold, so downstream readers filter from one table. Three cases give a
header-only edge table and no propagation: fewer than two accepted regions, no
`skani` binary on PATH, and regions skani cannot sketch. Any other skani
failure fails the run.

## Output specification

Per-genome outputs are split between `phase3_synthesis/` and top-level
run files.

Core synthesis files (`phase3_synthesis/`):

- `virosync_predictions.tsv` (canonical accepted predictions)
- `virosync_predictions.bed` (canonical accepted coordinates)
- `virosync_predictions.gff3` (canonical accepted annotations)
- `virosync_predictions_detailed.tsv` (all Phase 3 candidates)
- `virosync_summary.json`
- `evidence_profiles.json`
- `eve_ani_edges.tsv` (every accepted-EVE pair skani compared; header-only
  when the genome has fewer than two accepted regions)
- `virosync_tmvec_proteins.tsv` (can be header-only when TMVec is
  disabled or no hits are found)
- optional: `virosync_jelly_roll_proteins.tsv`
- `interproscan_summary.tsv` (can be header-only when InterProScan is
  disabled or no hits are found)

Top-level run files:

- `virosync_predictions_detailed.tsv` (copied summary table)
- `<genome_id>_eves.fna` (combined accepted EVE FASTA; generated when
  candidate regions are present)
- `virosync_tsv_invariant_report.tsv` (detailed table QA)
- `run.log` (timing and run summary)
- optional: `gvclass_results.tsv`

The output schema version is 5. Both prediction TSVs contain
`effective_eve_class`, with exactly one of `NCLDV`, `MIRUS`, `PPV`, `CRESS`,
`PHAGE`, `VIRAL_UNKNOWN`, or `UNKNOWN`. Canonical rows contribute to exactly
one class total. Result parsing maps the `VP` and `PLV` aliases to the parent
`PPV` class, and the `MIXED` alias to `VIRAL_UNKNOWN`.

The detailed TSV groups columns in this order:

1. identity and final calls;
2. candidate provenance;
3. marker evidence;
4. gene taxonomy;
5. composition and host evidence;
6. InterProScan evidence;
7. marker-set completeness;
8. ANI clustering.

The ANI columns are `ani_cluster_id` (`.` for a region with no clustered
relative), `ani_cluster_size` (`1` for a singleton), `ani_max_percent` (`.` for
a singleton), `taxonomy_class_before_ani`, and
`taxonomy_class_propagated_from`. The last two are filled only where a cluster
donor supplied the class.

`ppv_subtype` is `VP` or `PLV` only when subtype-specific marker evidence
supports one subtype and not the other. It is `.` for ambiguous PPV calls and
for every non-PPV parent class. The top-10 support columns are
`ncldv_top10_proteins`, `mirus_top10_proteins`, `ppv_top10_proteins`, and
`cress_top10_proteins`. The mutually exclusive `taxonomy_best_hits` partition
is ordered as:

``` text
EUK;MITO;PLASTID;BAC;ARC;UNK;NO_HITS;NCLDV;MIRUS;PPV;CRESS;GVMAG;PHAGE
```

The `*_top10_proteins` columns count raw top-10 prefix support. The disjoint
partition requires at least 25% amino-acid identity before assigning a viral
family. When markers do not assign a different concrete family,
identity-qualified gene taxonomy can assign CRESS. `vp_completeness` records
subtype evidence, while
`ppv_completeness` combines the PPV marker sets.

Batch TSV class columns are `ncldv`, `mirus`, `ppv`, `cress`, `phage`,
`viral_unknown`, and `unknown`. They are mutually exclusive and sum to
`accepted`.

Batch mode additionally writes batch summaries (`batch_summary.tsv`,
`batch_report.md`). The sibling `virosync-bench` repository uses
`batch_summary.tsv`, `virosync_predictions.tsv`, `virosync_predictions.bed`,
and `virosync_predictions.gff3` as its canonical count and coordinate
surfaces. Use `virosync_predictions_detailed.tsv` to inspect candidates
rejected by the acceptance gate.

## Reproducibility checks

The runtime freeze gate uses the same locked commands locally and in CI:

``` bash
pixi install --locked
pixi run lint
pixi run python -m pytest -q
pixi run python -m compileall -q src tests scripts
pixi run check-production-ready
pixi run virosync --version
```

The release manifest at
`release-manifests/resources_v1_0_6/RESOURCE_MANIFEST.json` is an exact copy of
the manifest inside the pinned archive. Production checks typed-parse both
shipped configs and compare their resource identities with this manifest and
the database source record. Scheduled and release-tag smoke runs additionally
perform a full resource verification, a clean shipped example, and an unchanged
resume whose schema-v3 artifact identities and counts must match the clean run.

Recommended smoke test:

``` bash
pixi run virosync \
  -i example/ \
  -o results/example_16t \
  --config config/orchestration.yaml \
  -w 1 \
  --threads-per-worker 16 \
  --clean-run
```

Quick validations:

``` bash
cat resources/virosync/DB_VERSION
cat results/example_16t/batch_summary.tsv
```

Expected operational targets for current resource pinning:

- core resource DB version: `v1.0.6`
- example batch status: `success`

The Python genome-parallel workflow runner is the benchmarked ViroSync
runtime. The harness in the sibling `virosync-bench` repository records
per-tool wall time and peak resident memory for the SynEVEs-2 and real-genome
panels. Manuscript Figure 5 summarizes those results. Its source table retains
the legacy name
`manuscript/figures/benchmark_figS1_runtime_memory_with_vr30_data.tsv`.
Runtime depends strongly on thread allocation and database caching and is
reported for completeness rather than as an optimized comparison.

Optional annotation/runtime notes:

- TMVec is optional and disabled by default. Install it with
  `pixi install --locked -e structural`, configure `phase3.tmvec_database_dir`,
  then enable it with `--tmvec` or `phase3.use_tmvec_database: true`.
- Supported TMVec confidence-path database keys are `bfvd`, `cath`, and
  `swissprot`. BFVD supplies viral structural support; CATH and
  SwissProt are reported as background/reference hits for specificity
  context.
- TMVec support requires the original trained TM-Vec weights, ProtT5,
  and precomputed embedding databases from the same model family.
  ViroSync does not use randomly initialized TMVec weights, incompatible
  `scikit-bio/tmvec-2` weights, or zero embeddings as evidence.
- Custom BFVD TMVec embeddings generated by legacy ViroSync builds that
  used the untrained fallback model are rejected during database
  loading/preflight and must be regenerated with trained TM-Vec weights
  before they can contribute viral structural support.
- InterProScan is optional and skipped when `interproscan_dir` is unset
  or missing.
- Boltz/Foldseek is optional and skipped when its runtime or structure
  DB is unavailable. Boltz uses `tools/boltz_runtime/pixi.toml`; the
  main Pixi environment does not install it by default. ViroSync
  currently writes sequence-only Boltz inputs, so
  `phase3.boltz_use_msa_server: true` is required for Boltz prediction.
- Optional paths are not pre-populated in the public default config;
  configure them with `virosync orchestrate setup` and validate with
  `pixi run -e structural check-structural-runtime --require-tmvec`
  and/or `--require-boltz`. The TMVec preflight includes a small model
  smoke test.

## Tool citation map

The following references cover the software and methods named in this
document.

- ViroSync: (NeLLi-team 2026).
- `prodigal-gv`: (Camargo et al. 2023; Hyatt et al. 2010).
- `pyhmmer` and HMM search engine context: (Larralde and Zeller 2023;
  Eddy 2011).
- DIAMOND: (Buchfink, Xie, and Huson 2015; Buchfink, Reuter, and Drost
  2021).
- TRF: (Benson 1999).
- RepeatMasker: (Smit, Hubley, and Green 2013).
- InterProScan: (Jones et al. 2014).
- TMVec: (Hamamsy et al. 2024).
- Boltz: (Wohlwend et al. 2024).
- Foldseek: (Kempen et al. 2024).
- GVClass: (Pitot, Bruna, and Schulz 2024).

The manuscript citation metadata and Quarto keys are in the sibling
`virosync-bench` repository at `manuscript/references.bib`.

## Notes on accuracy boundaries

This methods text is intentionally implementation-first. It documents
behavior verified in current code paths and config defaults. Older
references to CRF-led boundary expansion or novelty-first seeding should
be treated as legacy unless explicitly enabled and tested in a run
configuration.

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-benson1999trf" class="csl-entry">

Benson, G. 1999. “Tandem Repeats Finder: A Program to Analyze DNA
Sequences.” *Nucleic Acids Research* 27 (2): 573–80.
<https://doi.org/10.1093/nar/27.2.573>.

</div>

<div id="ref-buchfink2021diamond" class="csl-entry">

Buchfink, Benjamin, Klaus Reuter, and Hajk-Georg Drost. 2021. “Sensitive
Protein Alignments at Tree-of-Life Scale Using DIAMOND.” *Nature
Methods* 18 (4): 366–68. <https://doi.org/10.1038/s41592-021-01101-x>.

</div>

<div id="ref-buchfink2015diamond" class="csl-entry">

Buchfink, Benjamin, Chao Xie, and Daniel H. Huson. 2015. “Fast and
Sensitive Protein Alignment Using DIAMOND.” *Nature Methods* 12 (1):
59–60. <https://doi.org/10.1038/nmeth.3176>.

</div>

<div id="ref-camargo2023genomad" class="csl-entry">

Camargo, Antonio Pedro, Simon Roux, Frederik Schulz, Michal Babinski,
Yan Xu, Bin Hu, Patrick S. G. Chain, Stephen Nayfach, and Nikos C.
Kyrpides. 2023. “Identification of Mobile Genetic Elements with
geNomad.” *Nature Biotechnology* 42 (8): 1303–12.
<https://doi.org/10.1038/s41587-023-01953-y>.

</div>

<div id="ref-eddy2011hmmer" class="csl-entry">

Eddy, Sean R. 2011. “Accelerated Profile HMM Searches.” *PLoS
Computational Biology* 7 (10): e1002195.
<https://doi.org/10.1371/journal.pcbi.1002195>.

</div>

<div id="ref-hamamsy2024tmvec" class="csl-entry">

Hamamsy, Tymor, James T. Morton, Robert Blackwell, Daniel Berenberg,
Nicholas Carriero, Vladimir Gligorijevic, Charlie E. M. Strauss, Julia
Koehler Leman, Kyunghyun Cho, and Richard Bonneau. 2024. “Protein Remote
Homology Detection and Structural Alignment Using Deep Learning.”
*Nature Biotechnology* 42 (6): 975–85.
<https://doi.org/10.1038/s41587-023-01917-2>.

</div>

<div id="ref-hyatt2010prodigal" class="csl-entry">

Hyatt, Doug, Gwo-Liang Chen, Philip F. LoCascio, Miriam L. Land, Frank
W. Larimer, and Loren J. Hauser. 2010. “Prodigal: Prokaryotic Gene
Recognition and Translation Initiation Site Identification.” *BMC
Bioinformatics* 11 (1): 119. <https://doi.org/10.1186/1471-2105-11-119>.

</div>

<div id="ref-jones2014interproscan" class="csl-entry">

Jones, Philip, David Binns, Hsin-Yu Chang, Matthew Fraser, Weizhong Li,
Craig McAnulla, Hamish McWilliam, et al. 2014. “InterProScan 5:
Genome-Scale Protein Function Classification.” *Bioinformatics* 30 (9):
1236–40. <https://doi.org/10.1093/bioinformatics/btu031>.

</div>

<div id="ref-vankempen2024foldseek" class="csl-entry">

Kempen, Michel van, Stephanie S. Kim, Charlotte Tumescheit, Milot
Mirdita, Jeongjae Lee, Cameron L. M. Gilchrist, Johannes Soeding, and
Martin Steinegger. 2024. “Fast and Accurate Protein Structure Search
with Foldseek.” *Nature Biotechnology* 42 (2): 243–46.
<https://doi.org/10.1038/s41587-023-01773-0>.

</div>

<div id="ref-larralde2023pyhmmer" class="csl-entry">

Larralde, Martin, and Georg Zeller. 2023. “PyHMMER: A Python Library
Binding to HMMER for Efficient Sequence Analysis.” *Bioinformatics* 39
(5): btad214. <https://doi.org/10.1093/bioinformatics/btad214>.

</div>

<div id="ref-virosync_software_2026" class="csl-entry">

NeLLi-team. 2026. “ViroSync.” <https://github.com/NeLLi-team/virosync>.

</div>

<div id="ref-pitot2024gvclass" class="csl-entry">

Pitot, Thomas M., Tomas Bruna, and Frederik Schulz. 2024. “Conservative
Taxonomy and Quality Assessment of Giant Virus Genomes with GVClass.”
*Npj Viruses* 2 (1). <https://doi.org/10.1038/s44298-024-00069-7>.

</div>

<div id="ref-repeatmasker_open4" class="csl-entry">

Smit, A. F. A., R. Hubley, and P. Green. 2013. “RepeatMasker Open-4.0.”
<http://www.repeatmasker.org>.

</div>

<div id="ref-wohlwend2024boltz" class="csl-entry">

Wohlwend, Jeremy, Gabriele Corso, Saro Passaro, Noah Getz, Mateo Reveiz,
Ken Leidal, Wojtek Swiderski, et al. 2024. “Boltz-1 Democratizing
Biomolecular Interaction Modeling.” *bioRxiv*.
<https://doi.org/10.1101/2024.11.19.624167>.

</div>

</div>
