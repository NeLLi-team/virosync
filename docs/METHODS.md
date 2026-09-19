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
| `phase3_synthesis/virosync_predictions.bed` | Accepted EVE coordinates in BED6 format. |
| `phase3_synthesis/virosync_predictions.gff3` | Accepted EVE annotations in GFF3 format. |
| `<genome_id>_eves.fna` | Nucleotide sequences for accepted EVEs. This file is absent when there are no accepted EVEs. |
| `phase3_synthesis/virosync_summary.json` | Per-genome counts and run metadata. |
| `phase3_synthesis/eve_ani_edges.tsv` | ANI comparisons used for clustering. Can include candidates removed from the final accepted set. |
| `phase3_synthesis/gene_taxonomy/` | Per-candidate gene taxonomy tables. Written when at least one EVE is accepted. |

In `batch_summary.tsv`, `predictions` counts all candidates and `accepted`
counts rows in the accepted-prediction table. Each accepted region has one
class, so the class columns sum to `accepted`.

### Coordinates

The TSV and BED files use 0-based, half-open coordinates: `[start, end)`.
Therefore, `length = end - start`. GFF3 uses 1-based, inclusive coordinates.
The GFF3 interval `start + 1` through `end` describes the same bases as its TSV
and BED row.

### Key TSV columns

| Columns | Table | Interpretation |
| --- | --- | --- |
| `eve_id`, `scaffold`, `start`, `end`, `length` | Both | Region identity and final coordinates. |
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
| `tsd_sequence` | Both | Candidate target-site duplication immediately outside the repeat-defined interval, when found. |
| `recombinase_genes` | Both | Distinct protein identifiers with recombinase annotations, including nearby genes and repair recombinases. Inspect the evidence field for their locations and functions. |
| `integration_gene_evidence` | Both | JSON records with protein identifiers, genomic coordinates, location, enzyme family, annotation source, profile accession and scores. Includes DDE integrases and recombinases. A dot means no qualifying gene annotation. |
| `integration_hmm_status` | Both | `complete` means every selected predicted protein was searched, including an empty selection. `incomplete_sequence_length` means one or more selected proteins exceeded the HMM engine limit. `not_assessed` means the screen has not completed. |
| `integration_hmm_unsearched` | Both | JSON list of proteins excluded from the integration HMM search, with original identifiers, lengths, genomic coordinates, strand, EVE context, exclusion reason and limit. An empty list is `[]`; check the status to distinguish a completed screen from one not assessed. |

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
threshold. It rejects non-ACGT arms, arms with more than 80% of one
base, base entropy below 1.2 bits, and simple repeats with at least 80% periodic
identity for periods of 1–12 bp. The 50 bp minimum is an empirical search
cutoff, not a biological limit. Shorter, diverged or indel-rich TIRs can be
missed. The reported matched segment can cover only part of a longer repeat;
`tir_alignment_capped=1` identifies shortened segments. The full ungapped
alignment length is recorded in `tir_alignment_length`.

| `tir_status` | Meaning |
| --- | --- |
| `detected` | A qualifying pair defines this child's boundary. |
| `ambiguous` | Repeat evidence does not support an unambiguous split with intact marker evidence; heuristic coordinates are retained. |
| `shared_pair` | Distinct candidates would acquire the same interval; heuristic coordinates are retained. |
| `taxonomy_incomplete` | The marker span or a proposed repeat-defined interval lacks required gene-taxonomy coverage. |
| `not_detected` | No qualifying pair defines this region, including a retained neighboring segment without TIRs. |
| `no_marker_anchor` | No retained marker span was available to anchor the search. |
| `not_assessed` | Repeat screening has not been applied. |

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

### Interpreting integration signatures

TIRs can identify the ends of a complete linear viral element. Host sequence
outside both ends is needed to establish an integrated context. TIRs also occur
in virophage genomes for which integration has not been demonstrated
([Yutin et al., 2015](https://doi.org/10.1186/s13062-015-0054-9)).

An exact short direct repeat immediately outside the TIRs can be a target-site
duplication (TSD). Experimentally integrated mavirus has 5–6 bp TSDs
([Fischer and Hackl, 2016](https://doi.org/10.1038/nature20593)). The
`tsd_sequence` field reports the longest exact 3–9 bp match immediately outside
the reported pair. It excludes ambiguous bases and homopolymers. This is a
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
