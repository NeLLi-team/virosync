# Structure-derived capsid profiles

Resource bundle v1.1.0 contains 27 candidate capsid profiles built from
structure-supported protein clusters in an EVE catalogue. The profiles use
4,605 training proteins and belong to the 1,080-profile marker library.

## Functional annotations

| Predicted function | Profiles | Training proteins |
| --- | ---: | ---: |
| Penton protein | 20 | 4,518 |
| Minor capsid protein | 2 | 28 |
| Coat protein | 2 | 10 |
| Picornavirus-like capsid domain | 1 | 11 |
| General capsid protein | 1 | 14 |
| Major capsid protein (MCP) | 1 | 24 |
| Total | 27 | 4,605 |

Identifiers retain the form `EVE_CAPSID_NOVEL_NNNN`. Only
`EVE_CAPSID_NOVEL_0029` has an MCP role. Its strongest structural reference is
a [Pleurochrysis polinton-like-virus MCP](https://rest.uniprot.org/unisave/A0A2H4YBC6?format=txt&versions=6).
Related matches reach NCLDV capsid proteins. These short proteins do not
support an exclusive virus-family assignment, so the runtime annotation
assigns no marker family. A detected protein can still contribute a
taxonomic vote through its independent sequence matches.

Penton, minor capsid, coat, and other capsid candidates provide structural
marker evidence. They do not receive MCP status or its scoring bonus.
The `CAPSID` text in an identifier is not an MCP assignment.

## Source groups

| Source EVE context | Training proteins | Profiles containing that context |
| --- | ---: | ---: |
| Preplasmiviricota (PPV) | 4,564 | 26 |
| NCLDV | 18 | 3 |
| CRESS | 14 | 1 |
| Unresolved | 9 | 5 |

PPV combines polinton-like viruses, virophages, and unresolved PPV calls in
the source catalogue. Contexts can overlap within a profile; the last
column therefore does not sum to 27. These labels describe the surrounding
EVE. They do not establish a profile's taxonomic specificity or added PHAGE
or MIRUS capsid coverage.

## Curation and provenance

Six of the 33 source profiles are excluded from this capsid set:

- `0001`, `0006`, `0015`, and `0020` match a protein labeled "capsid protein
  23K" in the source metadata. The reference is an
  [adenovirus protease](https://www.uniprot.org/uniprotkb/D3X7C3/entry).
- `0018` has stronger phage tail-tube support, including chains in the
  experimental [T4 baseplate structure](https://www.rcsb.org/structure/5IV5).
- `0026` matches the N-terminal portion of a capsid precursor and has
  competing experimental protease-fold matches. A mature-MCP function is
  unresolved.

The release's
[curation table](https://github.com/NeLLi-team/virosync/blob/main/release-manifests/resources_v1_1_0/capsid_model_curation.tsv)
records each original identifier, source caption, inclusion decision,
predicted function, source-group counts, and profile checksum.
The runtime annotation table is
`models/model_annotations_with_interpro.tsv` within the resource bundle.

Pfam annotations use the shipped 937-profile screening library and its
gathering thresholds, scanned against all training members of the selected
profiles. An empty signature records no qualifying hit in that library;
it is not evidence that a protein lacks domains.

## Interpretation limits

Profile hits must pass the standard ViroSync marker-validation and region
filters.
Structural similarity supports a putative function; it does not establish
an active virion protein, a complete MCP, or a viral family.

Recovery of the training proteins tests profile loading and search behavior.
It does not measure sensitivity or false-positive rates on independent
sequences. The source collection was selected from existing EVE calls and
does not provide an independent genome-wide validation set.

The combined library retains hits for all 4,605 selected training proteins
after Pfam arbitration. Eight proteins from minor-capsid cluster `0027`
score best against penton profile `0028`; these functional labels overlap
in the current profile set. Both profiles remain non-MCP markers.

In a tested source locus, `0029` produces a strong profile match, but its
viral sequence-reference hits fall below the validation identity cutoff.
ViroSync retains that hit as a candidate; it does not grant validated MCP
support.
