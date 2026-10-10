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
| `notebooks/jupyter/eve_analysis.ipynb` | Per-genome report, including repeat-assessment summaries for accepted EVEs. Use the detailed table and repeat sidecar for rejected candidates. |

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

Both tables use output schema 9. The repeat-evidence groups below comprise 18
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
| Public output schema | 9 |
| Coordinate schema | 2: 0-based, half-open TSV and BED coordinates |
| Full Phase-2 checkpoint, `phase2/resume_state.json` | `virosync.phase2.resume_state/v4` |
| Refined-boundary checkpoint, `phase2/refined_state.json` | `virosync.phase2.refined_boundaries/v5` |
| Accepted and detailed prediction artifacts | `canonical-predictions-v7` and `detailed-predictions-v7` |
| Repeat-pair artifact | `virosync.repeat_candidates/v1` |

Resume validates the input, configuration, software, environment, resources,
phase checkpoints and output artifacts. The Phase-2 checkpoints retain the
nested repeat assessments; TSV and BED reports cannot reconstruct that state.
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

## MCP fold classification

`phase3_synthesis/virosync_jelly_roll_proteins.tsv` reports MCP support and
fold assignment as separate fields:

- `mcp_support` states whether an MCP-like match is only a `candidate`, is
  `sequence_supported`, or is `structure_supported`. Only a supported record
  counts as MCP evidence for region scoring through this analysis.
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

## Workflow summary

ViroSync predicts genes, finds viral marker proteins, and validates marker
matches against viral references. It builds candidate regions around supported
markers, then refines their boundaries with gene-level taxonomy and host
context. The last phase combines marker, taxonomy, composition, and optional
domain or structure evidence. Acceptance rules then select the reported
regions. Within each genome, ANI clustering can transfer a class from
MCP-bearing members that meet the transfer rules and agree on the class.

Optional InterProScan, TMVec2, and Boltz/Foldseek analyses annotate and score
existing candidates. They do not create regions without marker-seeded
candidates.

## Confidence and limitations

`final_confidence` is a rule-based score, not a calibrated probability that a
region is an EVE. The default tiers are `HIGH` at 0.7 or above, `MEDIUM` from
0.2 to below 0.7, and `LOW` below 0.2. Marker rules can accept a LOW candidate
without changing its reported tier. Several score terms reuse marker or
taxonomy evidence, so they are not independent observations.

Reference coverage limits marker validation and lineage assignment. ViroSync
uses predicted genes to refine boundaries, so gene prediction errors and
fragmented assemblies can affect the result. Host-like hits inside an EVE can
also lower its score. Scaffold order can change the background used for
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
coverage set boundaries before ViroSync recalculates composition and gene
evidence. These boundaries take precedence over heuristic trimming and
marker-floor widening. EVE acceptance rules still apply.

Disjoint pairs can divide a merged candidate into several EVE candidates.
Leading, intervening and trailing segments remain candidates when they contain
a validated marker or a qualified viral gene-taxonomy hit, including segments
without TIRs. Each child receives its own classification, marker counts and
verification. The detailed table includes every retained child. The accepted
table includes children that independently pass the usual acceptance rules.
The parent coordinates remain in `pre_tir_start` and `pre_tir_end`.

Identical repeats on neighboring elements can also form a pair spanning both.
When local and spanning pairs compete, the scanner retains the unsplit parent
with `tir_status=ambiguous`. Nested and overlapping alternatives receive the
same treatment, even when another part of the parent has an independent pair.
A split cannot discard or bisect a retained validated marker. Gene-taxonomy
evidence preserves a neighboring segment only when a qualified viral gene fits
wholly inside that segment. Other predicted genes may cross a repeat boundary.
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
| `ambiguous` | Repeat evidence does not support an unambiguous split with intact marker evidence; heuristic coordinates are retained. |
| `shared_pair` | Distinct candidates would acquire the same interval; heuristic coordinates are retained. |
| `taxonomy_incomplete` | The marker span or a proposed repeat-defined interval lacks required gene-taxonomy coverage. |
| `not_detected` | No qualifying pair defines this region, including a retained neighboring segment without TIRs. |
| `no_marker_anchor` | No retained marker span was available to anchor the search. |
| `not_assessed` | Repeat screening has not been applied. |

### Integration-associated genes

The default gene screen uses five bundled Pfam profiles through PyHMMER.
Both sequence and domain scores must pass the profile's curated gathering
threshold. It runs without InterProScan and covers the final EVE, its original
candidate span, and the configured flanking extension. Each annotation records
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
Integration-gene annotations do not seed EVEs or increase confidence scores.

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
