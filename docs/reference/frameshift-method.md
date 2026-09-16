# Frameshift screening method

ViroSync searches for disrupted viral marker domains in assembled DNA.
The detector uses the shipped `VS######` profile hidden Markov models and
runs with NumPy and PyHMMER in the standard environment.

Six-frame protein searches use 3,000-nucleotide chunks with 1,500-nucleotide
overlap to propose candidate windows on both DNA strands. Short chunks limit
the effect of long, stop-rich translations on HMMER's fast filters. A seed
needs only a local profile match; its alignment window extends by the full
profile length on each side.
A local codon-aware alignment then compares each window with its marker
profile. Match states consume three nucleotides for an ordinary codon or
one, two, four, or five nucleotides for a codon affected by a frameshift.
Profile insertion and deletion states accommodate ordinary amino-acid gaps.
The profile's emission and transition probabilities determine the alignment
score. Frameshift and stop penalties discourage disrupted paths.

Uncertain residues are reported as `X`; the assembly sequence is unchanged.
Each frame-disrupted path is compared with an ordinary-codon alignment of
the same DNA span. Candidate domains require protein HMM support
and viral-reference validation before they can seed a region.

## Screening thresholds

| Stage | Requirement |
| --- | --- |
| Six-frame seed | Protein HMM domain score at least 8 bits and independent domain E-value at most 10; HMMER filter thresholds are 0.1. |
| Native alignment | Profile insertion and deletion probabilities; frameshift probability penalties of 0.01 for a one-base change and 0.005 for a two-base change; stop probability penalty of 0.01. |
| Frame-disrupted candidate | At least 6 bits above an ordinary-codon alignment of the same locus. Stop-only candidates do not require this gain. |
| Protein HMM support | Sequence and independent domain E-values at most 1e-5; domain covers at least 50 percent of the profile. |
| Event support | At least eight observed amino acids on each side within the supported domain. Unknown residues do not count. |
| Viral-reference confirmation | The existing viral validation rules, including at least 25 percent amino-acid identity and 50 percent query coverage. |

Native `pid` values are blank because the profile alignment does not compute
pairwise identity. DIAMOND identity is recorded with the separate validation
evidence. Event codon spans use zero-based, half-open nucleotide coordinates.
Event protein and model positions are zero-based indices.

## Evidence and interpretation

The normalized hit table distinguishes native alignment scores from protein
HMM validation statistics. Neither quantity measures calibrated genome-wide
significance for this detector. Selecting a translated path to
fit a profile changes the statistical interpretation of a subsequent protein
HMM search. DIAMOND validation also uses that selected peptide.

The event table reports genomic codon spans and reading-frame changes.
Multiple equivalent alignments can place an event at nearby bases, especially
in repeats. An event count describes the selected alignment, not a proven
count of historical mutations.

Six-frame seeding can miss extensively fragmented or divergent loci. Resource
coverage limits which marker families can be detected. Stops are interpreted
with the standard genetic code; alternative-code genomes need separate
assessment. Introns, assembly errors, and programmed translation changes can
also produce disrupted alignments.

Native alignment windows are processed serially. The thread setting controls
the seed search and later validation tools. Runtime depends on the number
and size of seeded windows. The detector keeps the assembly and its translated
chunks in memory; bounded alignment windows do not imply bounded whole-genome
memory use.

## Relationship to BATH

[BATH](https://academic.oup.com/bioinformaticsadvances/article/4/1/vbae088/7693713)
uses codon-aware profile searches with frameshift-aware Forward algorithms
and statistical calibration. ViroSync uses a seeded local Viterbi detector
and retains separate viral-reference confirmation. It does not reproduce
BATH's statistical model. Comparable sensitivity or improved accuracy must
be demonstrated on labeled data; removing the installation dependency alone
does not establish either claim.
