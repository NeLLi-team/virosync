# Structure-derived candidate profiles

Resource bundle v1.1.0 contains 33 candidate profiles built from
structure-supported protein clusters returned by a capsid search of an EVE
catalogue. The set includes profiles with non-capsid functional evidence. The
profiles use 5,305 training proteins and belong to the 1,086-profile marker
library.

## Functional annotations

| Curated function | Profiles | Training proteins |
| --- | ---: | ---: |
| Putative penton protein | 20 | 4,518 |
| Putative adenain-like protease | 4 | 648 |
| Putative tail tube protein | 1 | 34 |
| Unresolved protease-like protein | 1 | 18 |
| Putative minor capsid protein | 2 | 28 |
| Putative coat protein | 2 | 10 |
| Putative picornavirus-like capsid domain | 1 | 11 |
| Putative capsid protein | 1 | 14 |
| Putative major capsid protein | 1 | 24 |
| Total | 33 | 5,305 |

Identifiers retain the historical form `EVE_CAPSID_NOVEL_NNNN`.
`EVE_CAPSID_NOVEL_0029` is the only profile registered as an MCP, and that
assignment remains putative. Its strongest structural reference is a
[Pleurochrysis polinton-like-virus MCP](https://rest.uniprot.org/unisave/A0A2H4YBC6?format=txt&versions=6).
Related matches reach NCLDV capsid proteins. These short proteins do not
support an exclusive virus-family assignment, so the runtime annotation
assigns no marker family by name. A detected protein can still contribute a
taxonomic vote through its independent sequence matches. A hit whose references
are only `EUK-pEVE` takes the class of the profile's training EVEs instead (see
[Profile classes](#profile-classes)).

The other 32 profiles are not MCPs. Their identifiers and source captions
do not confer MCP status. Capsid-type profiles give capsid evidence instead,
as described next.

## Capsid evidence

Twenty-seven profiles count as capsid evidence. They are the 20 penton, two
minor-capsid, two coat, one picornavirus-like capsid-domain, one capsid and
one major-capsid profiles: `0002` to `0005`, `0007` to `0014`, `0016`, `0017`,
`0019`, `0021` to `0025` and `0027` to `0033`. The curated function decides
membership; the curation table lists a penton, capsid or coat function for each.
`CAPSID_EVIDENCE_MODELS` in `src/virosync/pipeline/phase3/mcp_detection.py`
holds the list, and a test checks it against the curation table.

A hit to one of these profiles is treated as a viral capsid hit:

- Phase 1 lets it seed a candidate region when it validates on a viral or
  EUK-pEVE reference, and when it has no DIAMOND match and passes the HMM-only
  gates (`validated_novel`). Only MCP profiles seeded in that second case before.
- Phase 3 sets `has_capsid_evidence` on the region when `has_mcp` is true or
  any `validated` or `validated_novel` hit of the region is such a profile, even
  if a higher-scoring non-capsid profile on the same protein hides it from the
  marker count. A hit of another status (`supported`, `unvalidated`) gives none,
  even when it outscores the other profiles. The A1 seed export applies the same
  rule to the anchors of each seed. The flag feeds every scoring and acceptance rule that an MCP
  feeds: the MCP term of the marker score, the common bonus, the exemption from
  the small-region penalty, the confidence caps, the priority-marker floors, the
  non-HHG quality pass and the shortcuts of the v2 acceptance gate. The
  re-admit pass inherits the gate change.

`has_mcp` and `mcp_gene_ids` still name major capsid proteins only. The
class-vote override, ANI class transfer, jelly-roll classification, Boltz input
selection and the Phase 1 `high` seed confidence still need an MCP. The profile
`0029` remains the only structure-derived MCP.

The six other profiles (`0001`, `0006`, `0015`, `0020`, `0026` and `0018`) stay
viral markers. They validate and seed like any marker, but they give no capsid
evidence.

## EUK-pEVE references

Resource bundle v1.2.0 relabels the eukaryote-labelled reference hits of the 33
profiles from `EUK` to `EUK-pEVE`. These are 3,011 proteins, of which the marker
database holds 1,031 and the proteome database 1,892 after dereplication per
genome. The key of an affected genome is `EUK-pEVE__<assembly>`, and its
lineage in `taxonomy/labels.tsv` starts with `EUK-pEVE` and carries the suffix
`-pEVE` on every named rank. Viral labels on eukaryotic assemblies, such as
`PPV`, keep their class.

ViroSync counts `EUK-pEVE` as viral and never as cellular or host:

- A top-10 hit on an `EUK-pEVE` reference with at least 25% identity validates a
  marker and counts as a viral hit in boundary refinement, gene-level taxonomy
  and region scoring.
- When an `EUK__` and an `EUK-pEVE__` entry share one `<genome>|<protein>`, the
  `EUK-pEVE` row wins the twin rule.
- `EUK-pEVE` hits never enter the host-signature model or the host lineage share.
- `EUK-pEVE` carries no viral family. It sets no NCLDV, MIRUS or PPV flag, joins
  no dominant-family count and casts no family vote of its own. A marker hit
  whose references are all `EUK-pEVE` votes the class of its profile when the
  profile has one (see [Profile classes](#profile-classes)). Without a class
  (`0002`, other markers such as `GVOGm0003`) the hit stays `VIRAL_UNKNOWN`, and
  a region whose viral hits are all `EUK-pEVE`, in the marker database, the
  proteome database or both, uses the `MIXED` rule: it has the likely family
  `MIXED` and takes the mixed path of the v2 gate, whatever lineage its marker
  names or its Phase 1 region class suggest. When a lineage hit exists, the
  `EUK-pEVE` hit changes nothing: in the region vote, a gene whose viral hits are
  all `EUK-pEVE` casts no vote while any other gene votes for a lineage, apart
  from the marker vote of a profile class. The A1 seed export applies the same
  marker votes to the top-10 hits of the seed's validated markers. The
  `taxonomy_best_hits` partition counts these genes in the bucket `EUK-pEVE`; no
  output class has that name.

## Profile classes

For 32 of the 33 profiles, one class holds at least 80% of the training EVEs.
The counts are the `context_PPV`, `context_NCLDV`, `context_MIRUS`,
`context_PHAGE`, `context_CRESS` and `context_VIRAL_UNKNOWN` columns of
`release-manifests/resources_v1_1_0/capsid_model_curation.tsv`, with
`VIRAL_UNKNOWN` in the denominator. The classes come from the 2026-08-08
ViroSync screen of EVEs.

| Class | Profiles |
| --- | --- |
| `PPV` | 29: every profile except `0002`, `0018`, `0026` and `0033` |
| `MIRUS` | `0026` |
| `CRESS` | `0033` |
| `PHAGE` | `0018` |
| none | `0002` (NCLDV 9 of 15) |

`STRUCTURE_PROFILE_CLASSES` in `src/virosync/pipeline/phase3/mcp_detection.py`
holds the table, and a test derives it from the curation file.

The curation rows that say "EVE context is not a taxonomic assignment" still
hold for the profile. The class is used only for a marker hit whose qualified
viral references are all `EUK-pEVE`. A hit with a lineage reference, or with no
qualified viral reference, gets nothing from the class. For a hit that qualifies:

- The marker votes the class as an ordinary marker (weights 2 and 3), never as an
  MCP (flat weight 5, override). This holds for `0029` as well.
- The class is the marker family in Phase 3 and one marker of that family in the
  Phase 1 region family. The likely family and the gate use it, so a 3 kb region
  with only a validated `0008` hit and `EUK-pEVE` references takes the PPV path of
  the gate.
- The A1 seed export applies the same marker votes, with the same taxonomy
  lookup, to the validated markers, so A1 and Phase 3 publish the same class for
  the same evidence. A1 keeps the seed family for the likely family unless every
  qualified viral hit is `EUK-pEVE`.

Each protein counts once. When several profiles hit one protein, the one that wins
the per-protein deduplication of the vote gives the class (an MCP that can vote as
one first, then the higher score, then the first hit), in the marker family of
Phase 3, in the Phase 1 terminal-repeat children and in A1.

Exception: Phase 1 has no `CRESS` or `PHAGE` region family for any marker, so the
`CRESS` class of `0033` and the `PHAGE` class of `0018` apply from Phase 3 on (the
vote, the marker family and the gate). In Phase 1, in its terminal-repeat children
and in the `region_classification` of the A1 export, a hit of either profile gives
the region family `UNKNOWN`. The A1 `taxonomy_class` and likely family still follow
the vote.

Proteome gene records whose hits are only `EUK-pEVE` cast no class. Hits whose
references are all `EUK-pEVE` vote the profile class when the profile has one.
Without a class (`0002`, other markers such as `GVOGm0003`) they stay
`VIRAL_UNKNOWN`, and a region whose viral hits are all `EUK-pEVE` uses the `MIXED`
rule, with the likely family `MIXED`.

## Source groups

| Source EVE context | Training proteins | Profiles containing that context |
| --- | ---: | ---: |
| Preplasmiviricota (PPV) | 5,210 | 30 |
| PHAGE | 31 | 1 |
| MIRUS | 18 | 1 |
| NCLDV | 18 | 3 |
| CRESS | 14 | 1 |
| Unresolved viral context | 14 | 7 |

PPV combines polinton-like viruses, virophages, and unresolved PPV calls in
the source catalogue. Contexts can overlap within a profile; the last
column therefore does not sum to 33. These labels describe the surrounding
EVE. They do not establish a profile's taxonomic specificity or functional
assignment. ViroSync uses them only through the profile classes above.

## Curation and provenance

All 33 source profiles are retained. The source caption and curated function
are separate fields:

- `0001`, `0006`, `0015`, and `0020` retain the source caption `Product:
  capsid protein 23K`. Their curated function is putative adenain-like
  protease because the selected reference is an
  [adenovirus protease](https://www.uniprot.org/uniprotkb/D3X7C3/entry).
- `0018` retains the source caption `Capsid protein`. Its curated function is
  putative tail tube protein, supported by its strongest BFVD match and chains
  in the experimental
  [T4 baseplate structure](https://www.rcsb.org/structure/5IV5).
- `0026` retains the source caption `Major capsid protein`. Its curated
  function is unresolved protease-like protein. It matches the N-terminal
  portion of a capsid precursor and has competing experimental protease-fold
  matches. This evidence does not establish a mature-MCP function or an active
  protease.

The release
[curation table](https://github.com/NeLLi-team/virosync/blob/main/release-manifests/resources_v1_1_0/capsid_model_curation.tsv)
records each original identifier, source caption, retention decision, curated
function, source-group counts, and profile checksum.
The runtime annotation table is
`models/model_annotations_with_interpro.tsv` within the resource bundle.

Pfam annotations use the gathering thresholds from the shipped 937-profile
screening library. An empty signature records no qualifying hit among the
training members; it is not evidence that a protein lacks domains.

## Pre-extracted marker references

Pre-extraction yields 3,522 distinct reference proteins. Of these, 230 overlap
the existing marker set and 3,292 extend it, giving 3,867,592 marker proteins.
The additions comprise 2,997 EUK, 265 PPV, 12 NCLDV, 11 MIRUS, and seven
PHAGE references. These are source groups, not model-specificity assignments.

All 33 profiles search the broad reference proteome at a full-sequence
E-value cutoff of `1e-5`. Each profile searches all 46,765,750 filtered
proteins, so its E-values use the full database size. One protein exceeds
100,000 amino acids and is excluded, matching the pipeline's length limit.
No additional bit-score cutoff is applied during pre-extraction.

Qualifying hits contribute the full reference protein to `marker/marker.faa`.
Eleven source headers are normalized to existing genome keys for taxonomy
lookup; the extraction record preserves their original identifiers. Other
headers remain unchanged. A match does not change the reference's source
group or establish its function.

Existing marker proteins are preserved. Legacy organism-name suffixes are
resolved through the same taxonomy lookup used at runtime, and an existing
organism/protein pair is reused only when its sequence matches. Identical
proteins from distinct reference organisms remain separate.

One reference, `PHAGE__IMGVR_UViG_8103388301_000001_163`, has no detailed
lineage entry. It retains its source `PHAGE` label through the existing
runtime prefix fallback. Its detailed lineage remains unknown.

`marker/marker.dmnd` is rebuilt from the resulting FASTA. The runtime archive
contains that database; the source companion contains the matching FASTA.
The release
[extraction record](https://github.com/NeLLi-team/virosync/blob/main/release-manifests/resources_v1_1_0/marker_preextraction.json)
records the input checksums, thresholds, per-model counts, header mappings,
and sequence-preservation checks.

## Interpretation limits

Profile hits must pass the standard ViroSync marker-validation and region
filters.
Structural similarity supports a putative function; it does not establish
an active virion protein, a complete MCP, or a viral family.

Recovery of the training proteins tests profile loading and search behavior.
It does not measure sensitivity or false-positive rates on independent
sequences. The source collection was selected from existing EVE calls and
does not provide an independent genome-wide validation set.

The combined library retains new-profile hits for all 5,305 training proteins
after Pfam arbitration. Eight proteins from minor-capsid cluster `0027`
score best against penton profile `0028`; these functional labels overlap
in the current profile set. Both profiles give capsid evidence but are not MCPs.

In a tested source locus, `0029` produces a strong profile match, but its
viral sequence-reference hits fall below the validation identity cutoff.
ViroSync retains that hit as a candidate; it does not grant validated MCP
support.
