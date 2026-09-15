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
