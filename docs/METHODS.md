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

ANI class transfer can cross a chain of pairwise links, and shared host sequence
can contribute to an alignment. Mixed insertions and capsid exchange can also
conflict with a region-level class. Inspect the detailed table, gene-taxonomy
tables, and ANI edges for important calls. A missing ANI edge does not show
that two short regions are unrelated. Short sequences may not support a
comparison.
