# Outputs and interpretation

Use the accepted-prediction table for counts and later analyses. Use the
detailed table to check why ViroSync accepted or rejected a candidate.

## Output specification

For a batch run, ViroSync writes `batch_summary.tsv` and `batch_report.md` in
the selected output directory. It also creates one subdirectory for each input
genome. The genome ID is the input FASTA name without its extension.

The main per-genome files are:

| File | Contents |
| --- | --- |
| `phase3_synthesis/virosync_predictions.tsv` | Accepted EVEs. Use this as the main result table. |
| `virosync_predictions_detailed.tsv` | All Phase 3 candidates, including rejected and overlap-suppressed candidates. The same table also remains in `phase3_synthesis/`. |
| `virosync_repeat_candidates.tsv` | Retained direct-repeat pairs for all Phase 3 candidates, linked to final EVE identifiers. The same table remains in `phase3_synthesis/`; a header-only file has no retained pairs. |
| `phase3_synthesis/virosync_predictions.bed` | Accepted EVE coordinates in BED6 format. |
| `phase3_synthesis/virosync_predictions.gff3` | Accepted EVE annotations in GFF3 format. |
| `<genome_id>_eves.fna` | Nucleotide sequences for accepted EVEs. This file is absent when there are no accepted EVEs. |
| `phase3_synthesis/virosync_summary.json` | Per-genome counts and run metadata. |
| `phase3_synthesis/eve_ani_edges.tsv` | ANI comparisons used for clustering. Can include candidates removed from the final accepted set. |
| `phase3_synthesis/gene_taxonomy/` | Per-candidate gene taxonomy tables. Written when at least one EVE is accepted. |
| `phase3_synthesis/gene_functions.tsv` | Function labels of core genes in the accepted EVEs. Columns: `eve_id`, `porf_id`, `label`, `module`, `evidence`, `evalue`; one row per labelled gene. Written when at least one EVE is accepted and the resource bundle is v1.2.0 or later. See [gene-function labels](#gene-function-labels). |
| `phase3_synthesis/gene_function_length_exclusions.tsv` | Core proteins excluded from the gene-function profile search by its length limit, one row per EVE/protein pair. Columns: `eve_id`, `porf_id`, `status`, `reason`, `length_aa`, `limit_aa`. Written with `gene_functions.tsv`; a header-only table records no length exclusions. |
| `notebooks/jupyter/eve_analysis.ipynb` | Per-genome report, including repeat-assessment summaries for accepted EVEs. Use the detailed table and repeat sidecar for rejected candidates. |
| `eve_gene_map.syntviz.json` | Report figure: the shown EVEs with their flanks, gene categories, gene functions and DIAMOND protein hits between EVE cores, as a [Syntviz](https://pypi.org/project/syntviz/) session. `syntviz view` opens it in the interactive viewer. Written when the report shows at least one EVE. |
| `eve_gene_map_syntviz.png`, `eve_gene_map_syntviz.svg` | The same figure as images. Written only when Node.js is available. |

The Syntviz report retains every shown EVE row, gene and annotation. Its
all-versus-all DIAMOND search compares available core proteins from HIGH and
MEDIUM EVEs with `--more-sensitive`, E-value ≤ 1e-5 and identity ≥ 20%. LOW
EVEs and genes without protein sequences are not searched. The report streams
the hits and retains each query protein's best hit to another EVE, ranked by
descending bitscore, ascending E-value, then stable feature IDs and alignment
coordinates. Reciprocal hits collapse to one record. These sparse hits support
ribbons and similarity ordering; they do not alter predictions or gene content.
The absence of a ribbon does not establish the absence of homology.

Gene maps outline core HMM matches that qualify for scoring or a contained
spatial rescue. Their annotations retain the model, E-value and original
DIAMOND validation status. MCP matches receive the MCP symbol; a portal or
penton match does not. The marker heatmap counts distinct proteins per
qualified model. Candidates that bypassed scoring through the optional
phylogenetic rejection path have no hallmark assessment and contribute no
qualified heatmap entries. Their existing annotations remain available.
Validated MCP annotations remain visible when another profile on the same
protein has a higher score. These plotting IDs do not change the marker count,
confidence score or taxonomy vote.
These marks do not change the gene's taxonomy colour or establish viral
integration.

When Node.js and retained hits are available, confidence tiers with 3–500 rows
are ordered by protein similarity. Larger tiers retain confidence order.
The display shows at most 20,000 links. The saved session records whether the
search ran, the number of searched EVEs and eligible core proteins, qualifying
input rows (including reciprocal hits and multiple alignments), retained hits,
and the retention policy. Explicit search coverage records only searched EVE
pairs adjacent in the final row order, in both directions, including pairs
with zero hits. A later interactive reorder can therefore show an unknown
coverage warning for a newly adjacent pair.

A run that ends before Phase 3, for example because no marker hits or seeds
remain, writes its empty `virosync_predictions.tsv`,
`virosync_predictions_detailed.tsv`, `virosync_predictions.bed`,
`virosync_predictions.gff3`, `virosync_repeat_candidates.tsv`, and summary at the
genome root. It need not have a `phase3_synthesis/` directory. Empty TSVs retain
their headers. A run that reaches Phase 3 can also have zero accepted EVEs while
retaining rejected candidates and repeat pairs in the detailed outputs.

In `batch_summary.tsv`, `predictions` counts all candidates and `accepted`
counts rows in the accepted-prediction table. Each accepted region has one
class, so the class columns sum to `accepted`.

### Coordinates

The TSV and BED files use 0-based, half-open coordinates: `[start, end)`.
Therefore, `length = end - start`. GFF3 uses 1-based, inclusive coordinates.
The GFF3 interval `start + 1` through `end` describes the same bases as its TSV
and BED row.

### Key TSV columns

Both tables use output schema 10. The repeat-evidence groups below comprise 18
columns in each table; they are separate from the 17-column
[repeat-pair sidecar](#repeat-pair-sidecar).

| Columns | Table | Interpretation |
| --- | --- | --- |
| `eve_id`, `scaffold`, `start`, `end`, `length` | Both | Region identity and final coordinates. `eve_id` is an opaque identifier; see below. |
| `confidence_tier`, `final_confidence` | Both | Rule-based evidence tier and score. |
| `effective_eve_class` | Both | Published region class. |
| `hallmark_total`, `hallmark_unique`, `mcp_gene_ids` | Both | Viral hallmark and major capsid protein evidence. |
| `gene_taxonomy_viral_interior`, `gene_taxonomy_cellular` | Accepted | Viral and cellular gene support inside the region. |
| `taxonomy_best_hits`, `ncldv_top10_proteins`, `mirus_top10_proteins`, `ppv_top10_proteins`, `cress_top10_proteins` | Detailed | Gene-level taxonomy counts. |
| `kfd`, `gc_deviation` | Both | Sequence-composition differences from the host background used by scoring. |
| `canonical_selection_outcome` | Detailed | Reason for keeping or excluding a candidate. |
| `ani_cluster_id`, `ani_cluster_size`, `ani_max_percent` | Detailed | Within-genome ANI cluster assignment. A dot marks no clustered relative where applicable. |
| `taxonomy_class_before_ani`, `taxonomy_class_propagated_from` | Detailed | Class before ANI-based transfer and the donor EVE, when transfer occurred. |
| `tir_present`, `tir_status`, `tir_candidate_count` | Both | Terminal inverted-repeat evidence and whether a pair can define the boundary. Presence alone does not mean the boundary was changed. |
| `tir_left_start`, `tir_left_end`, `tir_right_start`, `tir_right_end`, `tir_identity` | Both | Coordinates of the matched repeat segments and their identity after reverse complementation. Segments may cover only part of a long or gapped repeat. For ambiguous results, these describe the best candidate pair. |
| `tir_scan_start`, `tir_scan_end` | Both | Genomic interval searched for repeat pairs. |
| `tir_alignment_capped`, `tir_alignment_length` | Both | Whether the reported arms were shortened to 500 bp, and the full ungapped alignment length used to select their outer ends. |
| `tir_boundary_override`, `pre_tir_start`, `pre_tir_end` | Both | Whether a repeat pair controls the final boundary, and the parent coordinates before repeat refinement. Split children share these parent coordinates. |
| `tsd_sequence` | Both | Candidate 3–9 bp duplication immediately outside the retained TIR anchors, when found. Check `tsd_assessment_status` before interpreting it at a changed endpoint. |
| `repeat_evidence_id`, `repeat_input_id`, `repeat_assessed_start`, `repeat_assessed_end`, `repeat_parent_start`, `repeat_parent_end` | Both | Stable evidence and input-seed linkage, the interval actually assessed, and its parent before TIR partitioning. |
| `repeat_left_status`, `repeat_right_status`, `repeat_left_assessment`, `repeat_right_assessment` | Both | Per-end status and JSON with requested and searched host windows, clipping reasons and ambiguous-base counts. A completed search is relative to the recorded method and windows. |
| `repeat_search_parameters`, `repeat_filter_counts` | Both | JSON with the discovery method, thresholds and limits, plus aggregate counts of excluded low-complexity and overlapping pairs. |
| `direct_repeat_candidate_count`, `direct_repeat_display_id` | Both | Number of retained pairs and a display representative. The display pair has no boundary authority. |
| `tsd_assessment_status`, `tsd_anchor_start`, `tsd_anchor_end`, `tsd_assessment` | Both | Explicit short-flank assessment and retained-TIR anchors. The JSON includes per-end limits and shifted-anchor diagnostic counts. No anchor is distinct from an assessed absence of a match. |
| `recombinase_genes` | Both | Distinct protein identifiers with recombinase annotations, including nearby genes and repair recombinases. Inspect the evidence field for their locations and functions. |
| `integration_gene_evidence` | Both | JSON records with protein identifiers, genomic coordinates, location, enzyme family, annotation source, profile accession and scores. Includes DDE integrases and recombinases. A dot means no qualifying gene annotation. |
| `integration_hmm_status` | Both | `complete` means every selected predicted protein was searched, including an empty selection. `incomplete_sequence_length` means one or more selected proteins exceeded the HMM engine limit. `not_assessed` means the screen has not completed. |
| `integration_hmm_unsearched` | Both | JSON list of proteins excluded from the integration HMM search, with original identifiers, lengths, genomic coordinates, strand, EVE context, exclusion reason and limit. An empty list is `[]`; check the status to distinguish a completed screen from one not assessed. |

### Repeat-pair sidecar

`virosync_repeat_candidates.tsv` contains one row per retained pair and final
EVE candidate. Join `eve_id` and `repeat_evidence_id` to the detailed table for
the assessment windows, statuses, parameters and filter counts. A shared source
assessment can appear under more than one final EVE identifier. Rows are ordered
by `eve_id`, then `candidate_id`; row order does not rank the alternatives.

| Sidecar columns | Interpretation |
| --- | --- |
| `eve_id`, `scaffold`, `repeat_evidence_id` | Final candidate identifier, host scaffold, and source assessment identifier. |
| `candidate_id` | Deterministic identifier for the pair within its source assessment. |
| `orientation` | `direct`: the two arms have the same orientation. |
| `left_start`, `left_end`, `right_start`, `right_end` | Host coordinates of the matched arms, using 0-based, half-open intervals. |
| `identity`, `alignment_length` | Matching-base fraction from 0 to 1 and ungapped alignment length in base pairs. |
| `method` | Discovery method recorded with the pair. |
| `interpretation` | `unresolved`: sequence similarity does not establish host or viral origin. |
| `outer_start`, `outer_end` | Interval including both arms: `left_start` through `right_end`. |
| `inner_start`, `inner_end` | Interval excluding both arms: `left_end` through `right_start`. |

A header-only sidecar means there are no retained pairs. It does not mean every
endpoint was searched completely. Read `repeat_left_status` and
`repeat_right_status` in the detailed table. When no evidence record exists,
those statuses and `tsd_assessment_status` are `not_assessed`, the pair count is
zero, and the remaining evidence fields are dots.

### Output compatibility and resume

The public output and coordinate versions are independent of checkpoint and
artifact versions:

| Contract | Current version |
| --- | --- |
| Public output schema | 10 |
| Coordinate schema | 2: 0-based, half-open TSV and BED coordinates |
| Phase-1 checkpoint, `phase1/resume_state.json` | `virosync.phase1.resume_state/v3` |
| Full Phase-2 checkpoint, `phase2/resume_state.json` | `virosync.phase2.resume_state/v5` |
| Refined-boundary checkpoint, `phase2/refined_state.json` | `virosync.phase2.refined_boundaries/v6` |
| Accepted and detailed prediction artifacts | `canonical-predictions-v7` and `detailed-predictions-v7` |
| Repeat-pair artifact | `virosync.repeat_candidates/v1` |

Resume validates the input, configuration, software, environment, resources,
phase checkpoints and output artifacts. The Phase-2 checkpoints retain the
nested repeat assessments and confirmed seed intervals, without a second
candidate-coordinate pair; TSV and BED reports cannot reconstruct that state.
Incompatible checkpoint schemas cannot be reused. A changed run fingerprint
starts a new run; otherwise only the validated sequence of completed phases
can be reused.

The repeat-candidate table is required even when it has only a header. If both
root and Phase-3 copies are recorded, they must agree. Missing or altered pairs,
invalid coordinates, inconsistent counts, or broken links to detailed candidates
prevent a completed-run resume. `--clean-run` requests a new run regardless of
the saved state.

### Region identifiers

`eve_id` has the form `EVE_<scaffold>_<start>-<end>` when no other candidate
shares the coordinates. A candidate that shares its coordinates with another
candidate, for example an ordinary candidate and a frameshift-rescue candidate,
carries the suffix `-c` followed by 16 hexadecimal characters, such as
`EVE_scaffold_1_1000-5000-c3f2a9b7c1d0e4f56`. Candidates created together at
the same coordinates all carry the suffix; a candidate added later at
coordinates an existing candidate already uses carries the suffix while the
earlier candidate keeps its plain identifier. The suffix is derived from the
candidate's internal seed and is deterministic for a given input and
configuration. The same identifier appears in the TSV, BED and GFF3 outputs.
Treat `eve_id` as an opaque identifier: read coordinates from the `scaffold`,
`start` and `end` columns rather than from the identifier.

`canonical_selection_outcome` records direct retention, gate rejection, loss
of a frameshift-rescued marker, overlap selection or suppression, or lack of
viral evidence. Both `kept` and `overlap_selected` can identify accepted
candidates. Use the accepted-prediction table to determine which candidates
ViroSync accepted.

### Class meanings

| Class | Meaning |
| --- | --- |
| `NCLDV` | Nucleocytoviricota |
| `MIRUS` | Mirusviricota |
| `PPV` | Preplasmiviricota, including virophage and polinton-like virus references |
| `CRESS` | Circular Rep-encoding single-stranded DNA viruses |
| `PHAGE` | Bacteriophage-like taxonomy |
| `VIRAL_UNKNOWN` | Viral evidence is present, but no single lineage wins |
| `UNKNOWN` | No published lineage class is resolved |

`ppv_subtype` is `VP` or `PLV` only when subtype-specific marker evidence
supports that label. A dot means unresolved or not applicable. An accepted
`UNKNOWN` region requires support from an ANI cluster with a marker-bearing
EVE; this support alone does not prove viral origin.

`EUK-pEVE` is a reference class, not a published class. It marks viral-like
proteins that sit in eukaryotic genome assemblies (resource bundle v1.2.0 and
later). ViroSync counts an `EUK-pEVE` hit with at least 25% identity as viral
evidence, in marker validation, boundary refinement, gene-level taxonomy and
region scoring, and never as cellular or host evidence. The class names no
lineage: it adds no vote for `NCLDV`, `MIRUS`, `PPV` or `CRESS`, and a lineage
hit in the same top 10 always decides. A gene or region whose viral hits are all
`EUK-pEVE` is `VIRAL_UNKNOWN`, unless a structure-derived profile supplies a
class (see [Profile classes](#profile-classes-for-euk-peve-only-hits)). In the
region vote, a gene whose viral hits are all `EUK-pEVE` casts no vote when any
other gene votes for a lineage, so it never dilutes the lineage majority. It votes
`VIRAL_UNKNOWN` only when no gene votes for a lineage. Votes of `VIRAL_UNKNOWN`
that come from conflicting lineage hits do not change. The A1 seed export has no
gene taxonomy. It applies the same marker votes to the top-10 hits of the seed's
validated markers, so it publishes the class that Phase 3 gives for the same
evidence. Whenever these rules decide the class (a profile class votes, or all
qualified viral hits are `EUK-pEVE`), A1 derives the likely family with the Phase 3
function, from the region family, the marker families and the class: a region
whose viral hits are all `EUK-pEVE` without a profile class is `MIXED`, and in a
mixed region one marker family must hold at least 70% of the family-bearing
markers to name the family. Otherwise the likely family stays the family of the
seed. The `taxonomy_best_hits` string counts each
interior gene once, in the buckets `EUK`, `MITO`, `PLASTID`, `BAC`, `ARC`, `UNK`,
`NO_HITS`, `NCLDV`, `MIRUS`, `PPV`, `CRESS`, `GVMAG`, `PHAGE` and `EUK-pEVE`.

The report gene maps put a gene in the category `EUK-pEVE` when its qualifying
viral hits are all `EUK-pEVE`, or when no viral hit qualifies and its best hit is
`EUK-pEVE`. The category has a color of its own and the legend entry "EUK-pEVE
(viral-like, eukaryotic genome)". A qualifying lineage hit gives the gene the
color of that lineage instead. The category distribution shows `EUK-pEVE` genes
as a segment of their own. The report never takes an `EUK-pEVE` hit, or a
lineage rank that ends in `-pEVE`, as host.

## MCP fold classification

`phase3_synthesis/virosync_jelly_roll_proteins.tsv` reports MCP support and
fold assignment as separate fields:

- `mcp_support` states whether an MCP-like match is only a `candidate`, is
  `sequence_supported`, or is `structure_supported`. Only a supported record
  counts as MCP evidence in this annotation.
- `type` reports the inferred protein fold: `DJR`, `SJR`, `HK97`, or `UNKNOWN`.
  `UNKNOWN` means unresolved fold. It does not mean that MCP support or
  taxonomy is absent.
- `effective_eve_class` in the prediction tables reports the region-level
  taxonomy consensus from qualified viral-reference matches. Fold assignment
  does not set this class.

A supported `Mirus_MCP` marker receives an inferred `HK97` type unless
stronger optional fold evidence takes priority. This is a marker-family
inference, not structural or experimental confirmation. The `confidence`
column in the MCP table describes fold evidence; it is distinct from
`final_confidence` for an EVE.

## Capsid evidence

A hit to one of the 27 capsid-type structure-derived profiles (penton, capsid and
coat functions; see the
[profile reference](reference/capsid-models.md)) counts as capsid evidence. The
region flag `has_capsid_evidence` is true when `has_mcp` is true or any
`validated` or `validated_novel` hallmark hit of the region is such a profile. The
test reads every hit before deduplication, so a higher-scoring non-capsid profile
on the same protein does not hide the capsid hit. A hit of any other status, such
as `supported` or `unvalidated`, gives no capsid evidence, even when it outscores
the other profiles on its protein. The A1 seed export applies the same rule to the
anchors of each seed, which carry the status of their marker hit. This flag
describes the existing marker and structural evidence and supplies the capsid
shortcuts of the class-specific acceptance gate. Graduated hallmark scoring and
spatial rescue use the separate exact-model catalogue described under
[Confidence and limitations](#confidence-and-limitations). A structural-anchor
exemption is distinct from a catalogue capsid assignment. Neither the flag nor
a model-name substring adds an unlisted model to that catalogue.

A region whose viral hits are all `EUK-pEVE`, in the marker
database, the proteome database or both, and whose markers supply no profile
class, has the likely family `MIXED`. This holds whatever lineage the marker names
or the Phase 1 region class suggest, and the gate reads a likely family of `MIXED`
before the region class. The region passes the gate through the same rule as a
`VIRAL_UNKNOWN` region of mixed markers: a length above 2,000 bp, and capsid
evidence or two hallmarks with one non-ATPase hallmark. A region with any lineage
viral hit keeps its family logic.

`has_mcp` and `mcp_gene_ids` name major capsid proteins only. The class-vote
override, ANI class transfer, jelly-roll classification and Boltz input
selection still need an MCP. `evidence_profiles.json` reports both flags.

The report gene maps (overview, detailed contexts and Syntviz) draw a black
outline around each gene with capsid evidence. A gene has capsid evidence when
it has a `validated` or `validated_novel` marker hit of an MCP model or of one
of the 27 profiles, or when the jelly-roll table supports it as an MCP. A
supported MCP also has a gold circle. In the Syntviz session, MCPs are in the
style group `mcp` and the other genes with capsid evidence are in the style group
`capsid`. Both groups draw the same outline.

## Profile classes for EUK-pEVE-only hits

The 33 structure-derived profiles come from EVEs of known classes. For 32 of
them, one class holds at least 80% of the training EVEs, counted in the
`context_*` columns of `release-manifests/resources_v1_1_0/capsid_model_curation.tsv`
(`PPV`, `NCLDV`, `MIRUS`, `PHAGE`, `CRESS` and `VIRAL_UNKNOWN`, with
`VIRAL_UNKNOWN` in the denominator). The classes come from the 2026-08-08 ViroSync
screen of EVEs. Twenty-nine profiles are `PPV`, `EVE_CAPSID_NOVEL_0026` is
`MIRUS`, `0033` is `CRESS` and `0018` is `PHAGE`. `0002` has no class (NCLDV 9 of
15). `STRUCTURE_PROFILE_CLASSES` in `src/virosync/pipeline/phase3/mcp_detection.py`
holds the table, and a test derives it from the curation file.

The curation table states that an EVE context is not a taxonomic assignment of
the profile. ViroSync therefore uses the class only for a marker hit whose
qualified viral references are all `EUK-pEVE`, which names no lineage. A hit with
any lineage reference, or with no qualified viral reference, gets nothing from its
profile. For such a hit the class has three effects:

- The marker votes the class in the region vote as an ordinary marker, with the
  weights of 2 and 3. It never takes the flat weight of 5 of an MCP and never
  overrides the other genes, also for `0029`. Proteome gene records whose hits are
  only `EUK-pEVE` still cast no class.
- The class is the marker family of the hit in Phase 3 (`marker_family_hits`,
  `marker_dominant_family`) and counts as one marker of that family in the Phase 1
  region family. The likely family and the acceptance gate then use the class: a
  3 kb region whose only marker is a validated `EVE_CAPSID_NOVEL_0008` hit with
  only `EUK-pEVE` references publishes `taxonomy_class` `PPV`, the likely family
  `PPV` and takes the PPV path of the gate. A `PHAGE` class has no gate family, so
  its likely family is `UNKNOWN`.
- The A1 seed export applies the same marker votes, with the same taxonomy
  lookup, so A1 and Phase 3 publish the same class for the same evidence.

Each protein counts once. When several profiles hit one protein, the profile that
wins the per-protein deduplication of the vote gives the class: an MCP that can
vote as one first, then the higher score, then the first hit. The marker family
in Phase 3, in the Phase 1 terminal-repeat children and in A1 follows the same
winner.

Exception: Phase 1 has no `CRESS` or `PHAGE` region family for any marker, so the
`CRESS` class of `0033` and the `PHAGE` class of `0018` apply from Phase 3 on (the
vote, the marker family and the gate). In Phase 1, in its terminal-repeat
children and in the `region_classification` of the A1 export, a hit of either
profile gives the region family `UNKNOWN`. The A1 `taxonomy_class` and likely
family still follow the vote.

Hits whose references are all `EUK-pEVE` vote the profile class when the profile
has one. Without a class (`0002`, other markers such as `GVOGm0003`) they stay
`VIRAL_UNKNOWN` and use the `MIXED` rule. Two classed profiles of different
classes tie and leave the region `VIRAL_UNKNOWN` with the likely family `MIXED`.

## Gene-function labels

A resource bundle of v1.2.0 or later lets ViroSync name the genes inside each
accepted EVE. Labels follow the gene vocabulary of Bellas and Sommaruga (2026,
bioRxiv, doi:10.64898/2026.06.19.733378, Fig. 4). They are labels only. They
never change detection, boundaries, scores, tiers or any other output column.
Flanking genes get no label.

ViroSync picks one label per core gene from the first source that gives one:

1. **Curated model.** The gene's best-scoring Phase 1 HMM hit at E ≤ 1e-5 names
   the gene when it is a curated ViroSync model. Evidence: the model name. A
   renamed or relaxed OG-backed model (sources `gvclass_existing_og_to_vs` and
   `OG_backed_relaxed_low_overlap_*`) names a gene only when its top Pfam domain
   occurs in at least 50% of the sampled members of the model and belongs to the
   vocabulary. Evidence: `<model> (OG; <domain> in k/n members)`. Metadata that
   cannot be trusted gives no label from that model: no sampled members, a
   malformed count, or a count above the number of sampled members.
2. **Pfam.** The strongest hit of the gene against the Pfam families of the
   vocabulary, at an i-Evalue of 1e-5 or lower. Evidence: `Pfam:<name>`.
3. **VOG.** The strongest hit against a VOGDB group whose consensus description
   names a vocabulary gene, at a full-sequence E-value of 1e-9 or lower.
   Evidence: `VOG:<id>`.

`module` groups the labels: capsid (MCP, mCP, ATPase, Protease, Portal,
Terminase, Triplex), replication (PolB, pPolB, PolA, SF1, SF3, AEP),
integration (YR, RVE, Int), mriya (VLTF2, VLTF3, HUH), other (MTase, Tet,
SLATT, Lipase) and hallmark (RNAPL, TFIIB). `evalue` is the E-value of the hit
that gave the label.

E-values follow the `hmmscan` convention: the database size is the number of
profiles in `gene_function.hmm`, not the number of genes a run searches. A
gene gets the same E-value and label whether the run searches it alone or with
thousands of other genes. The profile files and their selection are described in the
[resource bundle](RESOURCE_BUNDLE.md#gene-function-profiles-bundle-v120-and-later).

The profile search supports proteins of at most 100,000 amino acids, measured
after terminal `*` symbols are removed. The first proteome file containing a
protein supplies its sequence; a duplicate in the rescue FASTA cannot replace
an overlength original. Longer proteins remain in the input and EVE gene
content but are excluded from this profile search. The length-exclusion table
records `status=unsearched`, `reason=sequence_length_limit`, the processed
length and the limit for each affected EVE/protein pair. This is missing
assessment, not evidence that a protein lacks a function. Supported proteins
are searched with the same profile-count normalization, and existing curated
model evidence can still label an excluded protein. Both gene-function tables
are included in the completed run's artifact hashes.

Core genes include confirmed frameshift-rescued proteins inside the EVE, which
have no gene-taxonomy record. ViroSync searches them from the rescue FASTA of
Phase 1. Labels of rescued proteins appear in `gene_functions.tsv`; the report
gene maps draw the Phase 0 genes only.

On a bundle before v1.2.0, ViroSync logs one INFO line and writes neither
gene-function table, and the report draws the gene map as before. A bundle
with only one of the two gene-function files is corrupt, and the run stops with
an error that names the missing file. With the table, the static gene map
prints each label under its core gene and drops a label that has no free line,
and the Syntviz figure uses the label as the feature name and records it, with
its evidence, as a `function` annotation.

## Workflow summary

ViroSync predicts genes, finds viral marker proteins, and validates marker
matches against viral references. It builds confirmed seed regions around
supported markers and preserves each complete seed while extending outward.
The last phase combines marker, taxonomy, composition, and optional
domain or structure evidence. Acceptance rules then select the reported
regions. Within each genome, ANI clustering can transfer a class from
MCP-bearing members that meet the transfer rules and agree on the class.

Optional InterProScan, TMVec2, and Boltz/Foldseek analyses annotate and score
existing candidates. They do not create regions without marker-seeded
candidates.

### Reference hit selection

Marker validation and host-signature estimation read the top 10 DIAMOND hits of
each marker-domain query against the marker database. They use DIAMOND's rows
as they are. The host-signature model does not create seeds.

Boundary refinement, gene-level taxonomy and the Phase 3 DIAMOND check search the
combined proteome. That database holds the same subject ID several times, several
near-identical proteins from one genome, and some proteins under both a cellular
and a viral label. Taken directly, such copies fill the top 10 and push viral
references out. For these searches ViroSync asks DIAMOND for ten times the
requested number of targets and then keeps, for each query and in DIAMOND's
order:

1. For a protein (`<genome>|<protein>`) that appears under a cellular label
   (`EUK`, `BAC`, `ARC`, `MITO`, `PLASTID`) and under a viral label, the viral
   row at the position of the first occurrence. The cellular row is dropped.
   `EUK-pEVE` is a viral label here, so an `EUK-pEVE` row replaces the `EUK` row
   of the same protein.
2. The best hit of each target genome. The genome is the subject ID up to `|`
   including its label, so `EUK__GCA-013389655-1` and `PPV__GCA-013389655-1`
   count as different genomes.
3. The first N hits, where N is the requested number of targets.

The filter rewrites the DIAMOND table in place with the same columns. The
all-against-all search among EVE core proteins in the Syntviz figure is not a
top-hit search and uses no filter.

### Seed markers

An `EUK-pEVE` reference hit with at least 25% identity validates a marker like a
viral reference. A marker with no DIAMOND match (`validated_novel`) can start a
region only when its model is an MCP or gives
[capsid evidence](#capsid-evidence).

A validated marker can start a candidate region unless its HMM model is a host
kinase family: `VS000084`, `VS000087`, `VS000094`, `VS000191` and `VS000319`
(protein kinase and ABC1 kinase models). These families share distant homologs
with NCLDV kinases, so host kinases pass marker validation. A marker from one of
these models never seeds a region. Model names match without case and without a
suffix after the first `_` or `.`. The terminal-repeat step in Phase 2 follows the
same rule. Every repeat-bounded child must also receive at least one complete
confirmed seed. A viral hit alone cannot turn an unseeded residual slice into
an EVE. A partition that cuts or loses a confirmed seed keeps the parent whole.
A blocked marker can still appear in marker counts when another marker seeded
the region that contains it. Only catalogue-qualified evidence contributes to
the graduated hallmark score.

The additive hallmark rescue examines retained Phase 1 marker records, including
`supported` and `unvalidated` records, at E <= 0.01. It does not recover HMM
models discarded by marker retention. Three independent core families on three
proteins in one compatible panel and scaffold, within a full-gene span of at
most 15,000 bp, give a provisional MEDIUM floor. A mapped MCP or CRESS capsid
alone gives a MEDIUM floor; a distinct compatible core family on a different
protein at an edge-to-edge gap of at most 5,000 bp gives a HIGH floor. These
rules do not require a DIAMOND hit and do not change its recorded status.
The search bounds candidate genes by coordinates and skips groups lacking the
required distinct families. For each left endpoint it retains the earliest
qualifying right endpoint; partner searches include long overlapping genes
without clipping their coordinates.
Phase 3 indexes retained evidence by scaffold and start coordinate, then routes
fully contained records and witness memberships to each final boundary in their
original order. Functional annotation uses the same boundary coordinates and
the retained-marker index, including rescued frameshift genes.

Each witness not wholly contained in an ordinary confirmed seed becomes a
separate confirmed seed with a stable identity and its exact full-gene hull.
Partial overlap does not trim either seed. A fully contained witness needs no
duplicate seed. Ordinary seeds keep their original identities. The normal
extension, merge and repeat-partition rules preserve every assigned seed.

### Seed-preserving extension

The complete assembled Phase 1 interval, including its padding, or the exact
full-gene hull of a hallmark-rescue witness, is the protected seed. Its original identity and coordinates remain attached to every
derived candidate. The searched flanking envelope is separate from the EVE:
finding the first or last predicted gene does not remove gene-free seed ends.

One controller starts at the seed and walks flanking genes in genomic order,
independently on each side. Identity-qualified viral hits provide positive
support. Optional taxonomy ML supplies per-gene support without using the
host-signature model as a positive feature. Ambiguous genes, genes without
hits and isolated model-matched host genes do not extend the boundary on their
own. They remain pending until connected positive support occurs within the
five-gene lookahead and the searched nucleotide envelope.
Five unsupported genes end the support search. If the envelope ends after only
one to four unsupported genes, the reported limit is the search envelope.

Three consecutive exterior genes matching the host-signature model at the
configured host-lineage-share threshold (default 0.35), with no qualifying viral
support in those genes, stop
extension before the first gene in the run. A single eukaryotic best hit does
not establish this stop. A missing model or missing taxonomy is reported as
an evidence limitation; it does not imply host origin. Overlapping genes are
checked before a supported gene can extend the boundary past another feature.
A host block crossing the seed edge can stop new growth at that edge, but
cannot cut the seed.
Phase 3 reports the retained coordinates directly. It does not export a second
candidate interval or reduction fields. Integration-gene scans use the retained
interval plus the configured Phase 1 extension flank.
Host-lineage share measures the hit weight within the host's lineage; it does
not compare a gene's hit composition with the pooled model by weighted Jaccard.
Standalone checkpoint audits must load the configured taxonomy labels before
scoring; normal fresh and resumed Phase 1 runs initialize this lookup.
The configured taxonomy weighting applies to these hit shares, including the
stored per-hit bitscores when `taxonomy_weight_mode` is `bitscore`.

The existing `phase2.host_trim_enabled` and `host_trim_score_threshold`
settings now control this exterior stop. The former window, step, host-fraction,
viral-fraction, buffer and overlap-score settings belonged to the removed
trimming stages and are no longer accepted configuration keys.

The controller distinguishes a host stop, a support gap and the search-envelope
limit. A later ML, merge, repeat or rescue operation cannot cross a host stop.
Reaching the envelope records the limit without an automatic new search.
Each independently confirmed seed beyond a stop is processed separately.
Overlapping search envelopes alone do not merge candidates. An allowed merge
retains every source seed and the tightest recorded stop on each side. The
controller checks all reachable overlapping candidates, including those
separated in the list by an independent rescue or CRESS call. It does not merge
disjoint candidate intervals. A split assigns each original seed whole to
exactly one retained child,
or leaves the parent intact.

Final gene sets, markers, composition and optional evidence are computed for
the retained coordinates before scoring. Output and resumed-state checks
require complete seed containment. Incompatible boundary checkpoints must be
recomputed rather than reconstructed from late, mutable boundary fields.
Each Phase 1 checkpoint seed must exactly match its one original seed record.
Phase 2 also checks that its complete source membership agrees with Phase 1,
both when writing new boundaries and when resuming them.
Marker family votes are deduplicated by protein identity, including when HMM
queries carry domain-coordinate suffixes. Region marker counts retain the
highest-scoring profile per protein. After merging or repeat partitioning,
the final region class and marker counts are recomputed from markers in the
retained interval for every candidate using the same marker-retention policy as
final evidence synthesis, including the selected ablation. A3 uses that same
policy at its unchanged seed coordinates. Compact CRESS regions retain their
separate classification.
The three-gene stop is a working detection rule, not a calibrated biological
cutoff; its effects require comparison with positive controls.

## Confidence and limitations

`final_confidence` is a rule-based score, not a calibrated probability that a
region is an EVE. The default tiers are `HIGH` at 0.7 or above, `MEDIUM` from
0.2 to below 0.7, and `LOW` below 0.2. Marker rules can accept a LOW candidate
without changing its reported tier. Several score terms reuse marker or
taxonomy evidence, so they are not independent observations.

The packaged `hallmark_panels.json` records exact model aliases, core and
support roles, structural anchors, capsids and excluded profiles for seven
compatible panels. Panels group evidence; they do not assign taxonomy. The
catalogue is part of the run's resource fingerprint. Exact retained HMM scores,
E-values, query envelopes and measured domain scores are carried through Phase 1
resume state. Unmeasured diagnostics are absent, including those unavailable
from frameshift evidence. HMM envelopes do not define genomic seed bounds.

Graduated core evidence requires E <= 0.01 and either `validated_novel`, a viral
top-10 reference hit at >=25% identity, or a catalogue structural-anchor
exemption. Maximum family-to-protein matching prevents repeated profiles,
families or one protein from counting twice. Strong evidence has E <= 1e-5.
Three or more independent core groups receive weight 1 each. A pair receives
weight 1 per strong group and 0.5 per weak group, using the valid assignment
with the most strong groups. A singleton receives weight 1 if strong and 0.25
if weak. The marker component is capped at 1 and combines 0.2 times the weighted
core count with a capsid term of 0.4, discounted by the applicable weak weight.
Compatible context adds 0.05 per additional independent group, up to three,
scaled by the weighted-to-unweighted core count. Function and integration hits
can supply context at E <= 1e-5, but cannot become core groups or create seeds.
Only genes wholly inside the final retained interval contribute.

The viral-gene fraction uses only interior genes with a viral reference hit at
>=25% identity. Raw NCLDV/MIRUS or PPV prefix flags remain informational and
cannot replace a zero count of identity-qualified viral genes. The seed-cluster
component describes marker support separately from this gene-taxonomy term.

A confirmed frameshift rescue and an overlapping ordinary or rescued protein
carrying the same marker model on the same strand are treated as one scoring
locus. A protein that links different model groups is left separate. Phase 3 uses
one representative for gene-taxonomy counts, marker counts and lineage votes.
The same identity defines independent catalogue cores and score-qualified
spatial witnesses, including after boundary reassessment. Original witnesses
and confirmed seed intervals remain unchanged. Evidence profiles export the
protein-to-locus mapping used for scoring.
Capsid-fold counts and their bonus also count each linked locus once, using its
supported evidence before confidence rank. Linked alternatives cannot contribute extra independent
votes to these terms. Raw marker and fold records, rescued sequences, annotation
IDs and all confirmed seed intervals remain available. Distinct models, opposite strands and
nonoverlapping intervals are not joined by this rule.
This coordinate-based heuristic cannot distinguish a fragment bridge from a
rescued read-through across real tandem copies. A1 stops before Phase 3 scoring;
its marker counts remain raw translation diagnostics and its confidence is zero.

Ordinary ORF fragmentation remains a scoring limitation. In a paired
single-deletion control, one viral gene split into two ordinary predicted ORFs:
the viral-gene fraction changed from 1/6 to 2/7 and confidence from 0.308857 to
0.433523. Rescue-locus deduplication did not change either score because the
rescued translation was absent from the taxonomy table. This control supports
marker deduplication and seed preservation, but does not establish robustness to
fragmented genes or a rescue-dependent sensitivity gain. The aligner also
reported three frameshift events for one planted deletion; event counts do not
measure independent genomic lesions.

Non-marker components and cellular penalties remain active. Small-region
penalties use the catalogue's full-weight capsid qualification. A candidate
without qualifying core evidence is capped below MEDIUM; one with fewer than
two independent core groups is capped below HIGH before floors. A qualified
capsid supplies a MEDIUM floor and a capsid with a distinct compatible core
group supplies a HIGH floor without a distance restriction. The spatial rules
add their floors only when all witness genes remain fully contained. Floors
use the configured tier thresholds, whose defaults are 0.2 and 0.7.

A retained spatial witness can pass acceptance with an explicit rescue reason,
including for an `UNKNOWN` class, and survives the no-viral-evidence ANI filter.
This does not rewrite taxonomy or DIAMOND status. Detailed evidence records
carry the original status, qualification, component values, witness members
and floor provenance. These are provisional candidate rules, not biological
confirmation of viral origin or integration. Raw HMM support does not cancel
an exterior host stop; an independently confirmed seed beyond that stop keeps
its own protection.

Reference coverage limits marker validation and lineage assignment. ViroSync
uses predicted genes to refine boundaries, so gene prediction errors and
fragmented assemblies can affect the result. Independent evidence of cellular
homology inside an EVE can lower its score. Host-signature model outputs are
diagnostics and extension-stop evidence; they do not contribute a score penalty,
restrict the viral-enrichment bonus or supply a non-host-gene-count acceptance
alternative.
Preserving a confirmed seed does not guarantee final acceptance or establish
viral integration. Scaffold order can change the background used for
composition scoring. Structural matches are computational support, not
experimental validation.

The Phase 1 hallmark HMM search skips predicted proteins longer than 100,000
amino acids. A completed run therefore does not establish that all predicted
proteins were screened for viral hallmarks.

ANI class transfer can cross a chain of pairwise links, and shared host sequence
can contribute to an alignment. Mixed insertions and capsid exchange can also
conflict with a region-level class. Inspect the detailed table, gene-taxonomy
tables, and ANI edges for important calls. A missing ANI edge does not show
that two short regions are unrelated. Short sequences may not support a
comparison.

## Integration evidence

Three repeat observations have different roles:

| Evidence | Sequence relationship | Role |
| --- | --- | --- |
| Terminal inverted repeat (TIR) | Opposite orientations, compared after reverse complementation | A qualifying, unambiguous pair can refine or split a candidate boundary. |
| Short-flank match | Exact 3–9 bp match immediately outside retained TIR anchors | Candidate target-site duplication (TSD); evidence only. |
| Long direct-repeat pair | Same orientation, at least 50 bp under the current search method | Annotation near the assessed endpoints; does not change accepted boundaries. |

### Terminal inverted repeats

Terminal inverted repeats (TIRs) are paired sequences in opposite orientations.
The repeat screen searches the original, unmasked candidate sequence and its
extended flanks, limited to the region with gene-taxonomy coverage. The flank
length is `phase1.extension_kb` (5 kb by default). A qualifying pair must enclose
at least one complete validated marker-bearing protein. Repeated hits to the
same protein provide one anchor. Unambiguous pairs with complete gene-taxonomy
coverage can propose boundaries within the controller's reachable span.
A proposal must contain its assigned confirmed seeds and respect host stops.
Incompatible repeat arms remain annotations without changing coordinates.
ViroSync recalculates composition and gene evidence for the final span.
EVE acceptance rules still apply.

Disjoint pairs can divide a merged candidate into several EVE candidates.
Leading, intervening and trailing segments remain candidates only when they
receive complete confirmed seeds, including segments without TIRs. Every
source seed must belong to exactly one retained child. Each child receives its
own classification, marker counts and
verification. The detailed table includes every retained child. The accepted
table includes children that independently pass the usual acceptance rules.
If a proposed child lacks a confirmed seed, overlaps another child or leaves
a source seed unassigned, the whole split is rejected and the parent is kept.
The parent coordinates remain in `pre_tir_start` and `pre_tir_end`.

Identical repeats on neighboring elements can also form a pair spanning both.
When local and spanning pairs compete, the scanner retains the unsplit parent
with `tir_status=ambiguous`. Nested and overlapping alternatives receive the
same treatment, even when another part of the parent has an independent pair.
A split cannot discard or bisect a confirmed seed. Gene-taxonomy evidence alone
does not create an unseeded neighboring candidate. Other predicted genes may
cross a repeat boundary.
They retain their full coordinates under the usual interval-overlap rules,
so one crossing gene can appear in both
neighboring regions. Repeat boundaries that conflict across separate parent
candidates retain the original candidates.

The search uses exact 7 bp seeds and ungapped matches of at least 50 bp at 90%
identity or better. It scores matches at +2 and mismatches at −4 and uses the
full local alignment to select matching outer bases. It reports at most 500 bp
from each outer end; these reported segments must also meet the identity
threshold. It rejects a reported left arm containing non-ACGT bases, more than
80% of one base, base entropy below 1.2 bits, or at least 80% periodic identity
for periods of 1–12 bp. The right arm is checked against the identity threshold
without a separate complexity filter. The 50 bp minimum is an empirical search
cutoff, not a biological limit. Shorter, diverged or indel-rich TIRs can be
missed. The reported matched segment can cover only part of a longer repeat;
`tir_alignment_capped=1` identifies shortened segments. The full ungapped
alignment length is recorded in `tir_alignment_length`.

ViroSync uses GenomeTools 1.6.6 to index each candidate search interval with `suffixerator`
and locate exact seven-base seeds with `repfind`. Native output contains seed
locations rather than all possible repeat pairs. A seed contributes only when
its marker-specific left and right occurrence counts produce at most 4,096
pairings. Each eligible pair is generated once. Overlapping search intervals
share base comparisons; each distinct marker-clipped interval is scored separately.
This retains alternatives needed to identify ambiguous boundaries.

The segment scorer uses mismatch and score prefix sums with an identity-indexed
range query. It selects the highest-scoring eligible segment in O(n log n)
time per group, where n is the interval length. Ties prefer higher identity,
then longer alignment, then the earlier start.

| `tir_status` | Meaning |
| --- | --- |
| `detected` | A qualifying pair defines this child's boundary. |
| `ambiguous` | Repeat evidence does not support a split with whole confirmed seeds; the parent interval is retained. |
| `shared_pair` | Distinct candidates would acquire the same interval; parent intervals are retained. |
| `taxonomy_incomplete` | The marker span or a proposed repeat-defined interval lacks required gene-taxonomy coverage. |
| `not_detected` | No qualifying pair defines this region, including a retained neighboring segment without TIRs. |
| `no_marker_anchor` | No retained marker span was available to anchor the search. |
| `not_assessed` | Repeat screening has not been applied. |

### Integration-associated genes

The default gene screen uses five bundled Pfam profiles through PyHMMER.
Both sequence and domain scores must pass the profile's curated gathering
threshold. It runs without InterProScan and covers the retained EVE interval
plus the configured flanking extension. Each annotation records
whether its gene is interior, upstream, downstream, or crosses a boundary.

The integration HMM screen searches selected proteins up to 100,000 amino
acids long. Longer proteins keep their sequences and identifiers and appear
in `integration_hmm_unsearched` for each affected EVE. Their HMM assessment
is missing. A dot in `integration_gene_evidence` does not establish absence
of an integration gene when `integration_hmm_status` is incomplete.
Independent marker and InterProScan annotations remain available for these
proteins. Screening status does not change EVE acceptance or confidence.

For a batch containing both supported and oversized proteins, the sequence
E-value search-space parameter `Z` counts all distinct selected proteins,
including those not searched. Excluded proteins contribute to this count
but have no HMM comparison.
The profile gathering thresholds still use sequence and domain scores.
`complete` describes coverage of this screen's selected predicted proteins;
it does not certify gene prediction completeness or biological absence of
an integration mechanism.

| Pfam profile | Reported function |
| --- | --- |
| PF00589.28, Phage_integrase | Tyrosine recombinase |
| PF00239.27, Resolvase | Serine recombinase |
| PF00665.33, rve | DDE integrase domain |
| PF13333.13, rve_2 | DDE integrase domain |
| PF13683.13, rve_3 | DDE integrase domain |

The software includes profile accessions, source URLs, the CC0 license and
checksums in `virosync/data/integration_profiles.json`. The run fingerprint
covers both this manifest and the profiles. Existing marker annotations and
optional InterProScan annotations retain their source. General repair/recombination
proteins such as YqaJ are labelled separately from site-specific recombinases.
Integration-gene annotations do not seed EVEs. Catalogue-mapped interior hits
can contribute compatible context to the graduated hallmark score when core
evidence is present; flanking hits cannot contribute.

A domain match predicts a protein family. It does not establish catalytic
activity, identify the enzyme that caused an insertion, or distinguish a viral
element from a transposon. Flanking host genes are labelled by location. Repeat
matches can also arise within host repeat families. Use the reported search
scope and ambiguity status when interpreting a pair.

### Direct-repeat assessment

Direct-repeat evidence is recorded by default for genome candidates. The screen
uses unmasked host sequence around the endpoints assessed in Phase 2. The
`phase1.extension_kb` setting supplies the requested radius, 5 kb by default;
the method caps it at 10,000 bases. Scaffold ends and contiguous Phase-2
taxonomy coverage clip each window. The JSON endpoint assessments record both
the requested and searched windows, reasons for limits, and ambiguous-base
counts. A pair can be present when assessment is incomplete.

The inward half of each window stops at the input interval midpoint, so a pair
spans the two endpoint neighborhoods. `midpoint_partition` records this method
constraint and permits a completed assessment when no other limit applies.

The method uses exact forward seeds and ungapped alignments. Both arms must
pass the [TIR left-arm complexity filter](#terminal-inverted-repeats). The
recorded `repeat_search_parameters` are:

| Parameter | Value |
| --- | --- |
| `method` | `forward-7mer-ungapped-score-v1` |
| `seed_bp` | 7 bp |
| `min_arm_bp` | 50 bp |
| `min_identity` | 0.90 |
| `min_seed_entropy_bits` | 1.2 bits |
| `max_seed_chain_bp` | 493 bp between consecutive seed starts |
| `max_window_radius_bp` | 10,000 bp per endpoint |
| `max_seed_pairings` | 65,536 |
| `max_diagonals` | 8,192 |
| `max_candidates` | 128 retained pairs per assessment |

Seeds with fewer left–right occurrence combinations are searched first.
Overlapping arms and low-complexity arms are excluded, with counts retained in
`repeat_filter_counts`. These thresholds have not been calibrated as universal
biological cutoffs. Diverged or indel-rich repeats can be missed. Repeat absence
does not lower confidence or exclude a degraded EVE.

| Endpoint status | Meaning |
| --- | --- |
| `completed` | The recorded method completed within its defined endpoint window, without a limiting reason. Zero pairs means no qualifying pair was retained under that search. |
| `incomplete` | A window, sequence, partner, or search limit prevented a complete assessment. Retained pairs remain available. |
| `not_assessed` | No usable assessment exists at that endpoint, including missing sequence or coverage, or an endpoint changed after assessment. |

| Recorded reasons | Interpretation |
| --- | --- |
| `missing_sequence`, `missing_coverage` | Host sequence or a contiguous taxonomy-covered window is unavailable. |
| `contig_clipped`, `taxonomy_clipped`, `window_limit` | A scaffold end, taxonomy coverage, or radius cap reduced the requested window. |
| `insufficient_window`, `insufficient_paired_window`, `partner_unavailable` | One or both endpoints lack enough searchable sequence for paired arms. |
| `ambiguous_bases`, `low_complexity` | Sequence content limits the assessment. Ambiguous bases are not bridged as ordinary mismatches. |
| `seed_pairing_limit`, `diagonal_limit`, `candidate_limit` | The corresponding cap limited the search or retained alternatives. |
| `endpoint_changed` | The final endpoint differs from the assessed endpoint. |
| `midpoint_partition` | The defined inward-window split; this reason alone permits `completed`. |

The [repeat-pair sidecar](#repeat-pair-sidecar) lists the retained alternatives. Its
outer and inner intervals record the unresolved choice between viral terminal
repeats and duplicated host sequence. They are geometric alternatives, not
statistical confidence intervals or inferred attachment sites.
The interpretation remains `unresolved`; host repeats, terminal viral repeats
and nested elements can produce the same geometry. The
`direct_repeat_display_id` selects a representative for display, not a boundary.
Direct-repeat evidence does not alter accepted coordinates, gene membership,
identifiers, classes or confidence tiers.

### Short-flank assessment

Short-flank assessment uses the exact 3–9-base comparison at retained TIR
termini. The long direct-repeat pairs do not supply these anchors.
`tsd_anchor_start` and `tsd_anchor_end` identify the assessed TIR interval.

| `tsd_assessment_status` | Meaning |
| --- | --- |
| `assessed_match` | An eligible exact short match was found at the retained TIR anchors. |
| `assessed_no_match` | The flanks were assessed and no eligible match was found. |
| `incomplete` | Clipping or ambiguous bases limit one or both flanks; inspect the per-end reasons and any retained sequence. |
| `not_assessed_no_anchor` | No retained TIR anchor pair was available. |
| `not_assessed_endpoint_changed` | A candidate endpoint changed after the source assessment. |
| `not_assessed` | No repeat-evidence record exists; anchor absence is not established. |

The `tsd_assessment` JSON retains the anchor source, matched sequence, per-end
reasons and a local shifted-anchor diagnostic. The diagnostic translates both
termini together by each nonzero integer offset from −20 to +20 bases, excluding
clipped or ambiguous controls. `null_assessed` counts the usable offsets and
`null_matches` counts those with a qualifying short match.
These correlated comparisons summarize local matching frequency; they do not
estimate independent boundary accuracy or a calibrated chance probability.
A Phase-3 readmission that changes an endpoint retains its parent
evidence and pair alternatives, marks that endpoint unassessed, clears the
display representative and records `not_assessed_endpoint_changed` for the
short-flank assessment. Its sidecar rows carry the final EVE identifier, while
the assessed and parent coordinates continue to describe the source assessment.
At unchanged retained TIR anchors, the matched sequence agrees with
`tsd_sequence`.

### Interpreting integration signatures

TIRs can identify the ends of a complete linear viral element. Host sequence
outside both ends is needed to establish an integrated context. TIRs also occur
in virophage genomes for which integration has not been demonstrated
([Yutin et al., 2015](https://doi.org/10.1186/s13062-015-0054-9)).

An exact short direct repeat immediately outside the TIRs can be a target-site
duplication (TSD). Experimentally integrated mavirus has 5–6 bp TSDs
([Fischer and Hackl, 2016](https://doi.org/10.1038/nature20593)). The
`tsd_sequence` field reports the longest exact 3–9 bp match immediately outside
the retained TIR pair. It excludes ambiguous bases and homopolymers. This is a
candidate duplication, not proof of insertion; short matches can occur by chance.

Tyrosine-recombinase integration can produce matching attachment-site cores
with different lengths and positions from DDE-associated TSDs. Absence of a TSD
does not exclude that mechanism
([Bulzu et al., 2025](https://doi.org/10.1186/s40168-025-02148-0)).

Additional evidence to inspect outside ViroSync includes reads spanning both
host–virus junctions, an orthologous empty host locus, and insertion differences
between strains. Reverse-transcriptase and RNase H genes can help identify a
nested retrotransposon as the source of an integrase. Nested retrotransposons
occur in endogenous mavirus-like elements
([Hackl et al., 2021](https://doi.org/10.7554/eLife.72674)). ViroSync does not
currently test junction reads, empty loci, attachment sites or retrotransposon
domain architecture. Short junction microhomology alone should not determine
an integration mechanism.
