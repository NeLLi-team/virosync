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
assigns no marker family. A detected protein can still contribute a taxonomic
vote through its independent sequence matches.

The other 32 profiles remain non-MCP markers. Their identifiers and source
captions do not confer MCP status or its scoring bonus.

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
assignment.

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
in the current profile set. Both profiles remain non-MCP markers.

In a tested source locus, `0029` produces a strong profile match, but its
viral sequence-reference hits fall below the validation identity cutoff.
ViroSync retains that hit as a candidate; it does not grant validated MCP
support.
