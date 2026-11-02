# Scoring and acceptance

ViroSync scores each candidate with fixed rules and then applies a separate
acceptance gate. `final_confidence` is a rule-based score, not a calibrated
probability. A confidence tier and acceptance are different decisions: a LOW
candidate can be accepted, and a MEDIUM candidate can be rejected.

All evidence is computed on the final retained coordinates. Catalogue hallmarks,
support annotations and supplementary-module genes must lie wholly inside the
interval. Gene-taxonomy and marker counts use interval overlap.
[Outputs](../METHODS.md) defines the result files, and the [workflow map](https://github.com/NeLLi-team/virosync#workflow)
shows where each step runs.

## Steps from marker to accepted EVE

1. **Validate markers.** A marker HMM hit is `validated` when one of its ten best
   DIAMOND matches has the prefix `NCLDV`, `MIRUS`, `PLV`, `VP`, `PPV`, `CRESS` or
   `EUK-pEVE` at 25% identity or more. `GVMAG` and `PHAGE` matches do not validate
   a hit. A hit with no DIAMOND match can be `validated_novel` under separate HMM
   rules. `supported` and `unvalidated` hits stay diagnostic and earn no credit
   below.
2. **Find spatial witnesses.** Validated hits at E ≤ 0.01 to catalogue core models
   can form a witness (see [Spatial witnesses](#spatial-witnesses)). The widest
   witnesses outside ordinary seeds become seeds; nested or identical ones add no
   further seed.
3. **Set the interval.** Phase 2 extends each seed and sets its boundaries. Scoring
   reads the result; it does not move the coordinates.
4. **Collect context.** Interior integration genes and gene-function labels that
   the catalogue maps to a panel become support evidence at E ≤ 1e-5.
5. **Assess hallmarks.** One compatible panel is selected and gives the marker
   score, the coherence bonus and a possible floor.
6. **Assess the supplementary module.** For each candidate that reaches
   confidence scoring, the module asks whether at least three independent genes
   support a stable NCLDV or MIRUS class.
7. **Combine the components.** Weighted components, bonuses and penalties give a
   score, clipped to 0–0.99.
8. **Apply caps, then floors.** Missing hallmark evidence caps the score, and
   qualified evidence sets a minimum.
9. **Assign the tier.** HIGH is 0.7 or above, MEDIUM is 0.2 to below 0.7, and LOW
   is below 0.2.
10. **Apply the acceptance gate and selection.** The class and evidence decide
    whether a candidate is kept, and overlapping rescue and ordinary candidates are
    arbitrated.
11. **Reassess rejected neighbors.** One pass can join two adjacent rejected LOW
    candidates and score the joined interval. The gate and selection then run again.
12. **Cluster and filter.** Accepted EVEs are clustered by ANI, an agreed MCP class
    can transfer within a cluster, and candidates without viral evidence are
    removed unless a marker-bearing relative or a spatial witness supports them.

## Hallmark assessment

[`hallmark_panels.json`](https://github.com/NeLLi-team/virosync/blob/main/src/virosync/data/hallmark_panels.json)
assigns exact models to panels, families and roles. It holds 332 links in 13
panels. A panel groups compatible evidence. It does not assign taxonomy, and
the published class comes from the taxonomy vote.

Each family counts once in a panel, and each protein supplies one family. Maximum
family-to-protein matching enforces both limits. Repeated profiles of one family,
or several families on one protein, do not add cores.

The panel with the most cores wins. Ties go to more strong cores, then capsid
presence, then more compatible families, then the panel name that sorts last.

| Quantity | Rule |
| --- | --- |
| Eligible hit | Status `validated` or `validated_novel`, E ≤ 0.01, not excluded |
| Strong core | E ≤ 1e-5 |
| Weighted cores | Three or more cores: 1 each. Two cores: 1 if strong, 0.5 if weak, using the assignment with the most strong cores. One core: 1 if strong, 0.25 if weak |
| Capsid fraction | 1 for a strong capsid core; otherwise the weak weight above when a capsid core is present; otherwise 0 |
| Marker score | min(0.2 × weighted cores + 0.4 if the capsid fraction is 1, 1); a weak capsid adds 0.4 × capsid fraction |
| Coherence bonus | 0.05 × min(compatible families − 1, 3) × weighted cores / cores |

Support-only links need E ≤ 1e-5. They add compatible families inside a panel
that has a qualifying core. They never add a core, a floor or a witness.

### Spatial witnesses

A witness needs complete genes on one scaffold, in one panel:

| Rule | Evidence | Floor |
| --- | --- | --- |
| `capsid` | One capsid core | MEDIUM |
| `three_core_15kb` | Three independent core families on three proteins, in a full-gene span of at most 15,000 bp | MEDIUM |
| `capsid_partner_5kb` | A capsid core and a different core family on a different protein, at most 5,000 bp apart edge to edge | HIGH |

A witness counts only when all of its members remain inside the final interval.

### Model assignment

The table lists each family with its number of exact models. Capsid-bearing
families are in bold. `EVE_CAPSID_NOVEL_0029` sits in the MCP family of the
Nucleocytoviricota and Polinton-like virus panels with no capsid role, so it adds
a core without capsid credit.

<!-- Checked against hallmark_panels.json by tests/test_scoring_docs.py. -->
| Panel | Core families | Support-only families |
| --- | --- | --- |
| Nucleocytoviricota | A32 (5), **MCP (18)**, PolB (4), RNAPL (3), RNAPS (2), SFII (2), TFIIB (6), TopoII (1), VLTF3 (8) | A32 (2), D5_primase_domain (2), **MCP (5)**, mRNA_capping (4), PolB (2), TFIIB (2), VLTF3 (1) |
| Mirusviricota | **MCP (2)**, PolB (4), Portal (2), RNAP (5), Terminase (3), Triplex (4) | **MCP (1)**, PolB (2) |
| Virophages | ATPase (4), **MCP (9)**, Penton (27), Protease (2) | ATPase (1), **MCP (1)**, Penton (3), Protease (6), RVE (4) |
| Polinton_like_viruses | ATPase (1), **MCP (41)** | ATPase (1), **MCP (1)**, mCP (24), pPolB (1), Protease (6), YR (1) |
| Maverick_Polintons | ATPase (1), **MCP (4)** | ATPase (1), **MCP (1)**, mCP (24), pPolB (1), Protease (6), RVE (4) |
| Mriyaviruses | A32 (5), HUH (2), **MCP (1)**, VLTF2 (6), VLTF3 (7) | A32 (1), HUH (3), **MCP (1)**, VLTF3 (1) |
| CRESS | **Capsid (11)**, Rep (4) | Rep (1) |
| Unresolved_Penton | Penton (1) | none |
| Unresolved_Coat | Coat (2) | none |
| Unresolved_Picornavirus_like_capsid_domain | Picornavirus_like_capsid_domain (1) | none |
| Unresolved_Capsid_unresolved | Capsid_unresolved (1) | none |

Two more unresolved panels hold one model each, `EVE_CAPSID_NOVEL_0018` and
`EVE_CAPSID_NOVEL_0026`, and both are excluded from scoring. The fusion families
`SF1/SF3_fusion` and `PrimPol/SF3` are also excluded. Unresolved panels hold one
family, so their models cannot combine with established panels.

ViroSync has no hallmark panel for phages or retroelements. `PHAGE` is reached
through the taxonomy vote. `EUK-pEVE` is a reference class, not a published class
([class meanings](../METHODS.md#class-meanings)).

## Confidence score

```text
score = weighted components + bonuses + 0.15 × structural score
        − EUK penalty − size penalty
```

### Weighted components

The default run uses five components. Their weights are normalized to sum to 1.
A composition score of zero removes that weight from the sum.

| Component | Raw weight | Normalized | Without composition | Value |
| --- | ---: | ---: | ---: | --- |
| Marker | 0.24 | 0.293 | 0.375 | Marker score above |
| Gene taxonomy | 0.24 | 0.293 | 0.375 | min(0.7 × viral-gene fraction + 0.3 × dominant-family fraction, 1) |
| Composition | 0.18 | 0.220 | removed | 0.4 × min(KFD / 0.3, 1) + 0.3 × min(GC deviation / 0.15, 1) + 0.3 × min(codon-usage deviation / 0.5, 1) |
| InterProScan | 0.08 | 0.098 | 0.125 | Viral-domain score; 0 when InterProScan is off |
| Seed cluster | 0.08 | 0.098 | 0.125 | 1 for at least four VP/PLV, three NCLDV or three MIRUS markers; 0.5 for two of one family; otherwise 0 |

The viral-gene fraction counts interior genes with a viral reference hit at 25%
identity or more. `phase3.use_crf_in_final_score` adds a legacy boundary
component with raw weight 0.18, reduced to 0.12 for four genes or fewer and
to 0.06 for two genes or fewer. It is off by default.

### Bonuses and penalties

| Term | Condition | Value |
| --- | --- | --- |
| Domain agreement | Seed cluster 1 and InterProScan score ≥ 0.5 | +0.03 |
| Boundary and composition | Boundary confidence ≥ 0.9 and composition ≥ 0.6 | +0.03 |
| Family consistency | Two region, taxonomy, marker or domain votes for one NCLDV, MIRUS or PPV family; three or more | +0.04; +0.05 |
| Viral enrichment | Viral-gene fraction ≥ 0.15 | min(0.12, fraction × 0.4) |
| Hallmark coherence | Compatible families beyond the first | See the table above |
| Structural | Optional structural analysis | 0.15 × structural score |
| EUK penalty | Fraction of genes with three eukaryotic hits at 70% identity or more and no viral hit, × 0.20; or, above 80% high-identity eukaryotic genes and below 10% viral genes, that fraction × 0.12. The larger term applies | up to −0.12 |
| Size penalty | One, two or three genes, without full capsid credit (capsid fraction 1) | −0.15, −0.08 or −0.04 |

A family vote from taxonomy needs a dominant-family fraction of 0.25 or more, and a
marker vote needs a positive marker-complement score. The scorer also holds
taxonomy-distribution and host-divergence bonuses. The default flow does not
produce their inputs, so both stay at zero.

### Caps and floors

Caps apply first, in this order:

1. A candidate with no qualifying core and no eligible supplementary module stays
   below MEDIUM, at 0.199.
2. A candidate with fewer than two cores stays below HIGH, at 0.699.

Floors then raise the score:

| Evidence | Floor |
| --- | --- |
| Capsid core and a second core in the selected panel, at any distance | HIGH, 0.7 |
| Capsid core alone | MEDIUM, 0.2 |
| Contained spatial witness | The witness tier |

### Supplementary module

The module runs for every candidate that reaches confidence scoring, and it lifts
only the zero-core cap and the LOW gate below. It adds no score term. A candidate
is eligible when:

- The LOW and MEDIUM class resolvers agree on `NCLDV` or `MIRUS`.
- At least three nonoverlapping genes inside the interval each have a retained
  hit of that class at 25% identity or more, E ≤ 1e-5 and query coverage of 50%
  or more per alignment.
- The viral bitscore of each counted gene exceeds its best nuclear or prokaryotic
  hit, when the gene has one. Organelle hits are reported as ambiguity and are not
  comparators.

## Acceptance gate

A candidate with a contained spatial witness that supports its tier passes with
the reason `hallmark_spatial_rescue_pass`. Other candidates follow this table:

| Class | HIGH or MEDIUM | LOW |
| --- | --- | --- |
| `NCLDV`, `MIRUS` | Length above 5,000 bp, or capsid evidence | Length above 5,000 bp, and two hallmarks or an eligible module of the same class |
| `PPV`, `MIXED` | Length above 2,000 bp, and capsid evidence or two hallmarks with one non-ATPase | Same rule |
| `CRESS` | CRESS marker family and at least one hallmark | Same rule |
| Other classes | Rejected | Rejected |

`MIXED` is an internal gate class. The hallmark count here is the number of
marker-bearing proteins, one per scoring locus. It is separate from the panel
core count above.

Capsid evidence in this table is the [`has_capsid_evidence` flag](../METHODS.md#capsid-evidence).
It is separate from the capsid role in the catalogue.

## Joint rescue of rejected seeds

One nonrecursive pass after gating can join two adjacent, disjoint, rejected LOW
candidates of the same `NCLDV` or `MIRUS` class. Their searched envelopes must
overlap, every gene in the joint interval must have taxonomy, and no host stop,
host-positive run, third candidate or internal repeat boundary may separate them.
The joined interval is scored from the start. It replaces its two children only
when it passes the gate and the supplementary module counts at least one gene
inside each child's confirmed seeds.
[Bounded joint assessment](../METHODS.md#bounded-joint-assessment-of-rejected-seeds)
gives the full rules.

## Where each decision is recorded

| Decision | Output |
| --- | --- |
| Score and tier | `final_confidence`, `confidence_tier` |
| Published class | `effective_eve_class` |
| Qualifying core families | `qualifying_hallmark_genes` |
| Gate and selection outcome | `canonical_selection_outcome` in the detailed table |
| Joint rescue | `joint_seed_provenance` in the detailed table |
| Components, weights, floor, witnesses and supplementary module | `score_components` and `supplementary_viral_evidence` in `evidence_profiles.json` |
