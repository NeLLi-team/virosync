# Scoring and acceptance

## Reported capsid admission and confidence

New seed and evidence admission requires a full-sequence HMM score of at least 10 bits for a
semantic capsid, coat or Penton model. This floor also applies to validated
EUK-pEVE-backed capsids; it does not change noncapsid thresholds. Raw reported
hits remain provenance. Numeric scores read from TSV retain the same eligibility
as in-memory scores; missing, malformed or nonfinite scores do not qualify.
Once admitted, full assigned confirmed seed identities
and spans remain protected through trimming, merging and partitioning.

An eligible reported capsid supplies a MEDIUM floor only when it has validated
viral evidence, a score of at least 75 bits, or positive viral identity in its own
retained marker top hits. E-value alone does not supply this floor. Other eligible
reported capsids retain a LOW witness. A neighboring viral gene is separate
regional evidence, not a self-gene viral hit. These admission choices can remove
real weak capsid hits; no claim of zero biological loss is made.


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

1. **Validate markers.** Marker searches retain 20 matches after reference masking,
   viral/cellular protein-twin selection and reduction to one hit per labelled
   genome. A marker HMM hit is `validated` when one of these DIAMOND matches
   has the prefix `NCLDV`, `MIRUS`, `PLV`, `VP`, `PPV`, `CRESS` or
   `EUK-pEVE` at 25% identity or more. `GVMAG` and `PHAGE` matches do not validate
   a hit. A hit with no DIAMOND match can be `validated_novel` under separate HMM
   rules. Reported capsid and Penton hits can seed and receive a MEDIUM floor
   even when `supported` or `unvalidated`; their validation status stays unchanged.
   Other markers keep their validation and HMM thresholds.
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

Numerical scoring also considers two pools: `PPV` combines Virophages,
Polinton_like_viruses and Maverick_Polintons; `NCLDV_Mirus` combines
Nucleocytoviricota and Mirusviricota. These pools affect marker score, coherence
and the independent-core cap only. Spatial witnesses and automatic floors still
require evidence within an original panel. No new discovery seeds follow from a pool.

Each pool uses maximum family-to-locus matching. Equivalent MCPs within PPV,
Penton/mCP and shared PolB models count once. Mirus RPB1 models map to RNAPL;
RPB2 models map to RNAPS. NCLDV and Mirus capsid architectures remain distinct,
as do A32 and Mirus Terminase. Each model keeps its core or support-only role.

The original panel with the most cores wins. Ties go to more strong cores, then
capsid presence, then more compatible families, then the panel name that sorts
last. Only the pool containing this winning panel can enrich its numerical
evidence. An unrelated pool cannot displace it. The original panel wins an equal
evidence tie with its pool. Pool enrichment retains at least the original panel's
coherence bonus, so a new weak core cannot dilute established strong context.

For PPV subtype composition, a PLV-prefixed MCP model adds a PLV vote only when
its catalogue core roles belong solely to Polinton_like_viruses. Models shared
with Maverick_Polintons, or assigned only to Nucleocytoviricota, add no PLV vote
from their names. Qualified raw reference hits and resolved subtype-specific
taxonomy can still resolve the subtype.

| Quantity | Rule |
| --- | --- |
| Eligible hit | Status `validated` or `validated_novel`, E ≤ 0.01, not excluded |
| Strong core | E ≤ 1e-5 |
| Weighted cores | Three or more cores: 1 each. Two cores: 1 if strong, 0.5 if weak, using the assignment with the most strong cores. One core: 1 if strong, 0.25 if weak |
| Capsid fraction | 1 for a strong capsid core; otherwise the weak weight above when a capsid core is present; otherwise 0 |
| Marker score | min(0.2 × weighted cores + 0.4 if the capsid fraction is 1, 1); a weak capsid adds 0.4 × capsid fraction |
| Coherence bonus | 0.05 × min(compatible families − 1, 3) × weighted cores / cores |

Ordinary support-only credit needs E ≤ 1e-5 and a qualifying core in the same
panel. It never adds a core or a HIGH floor. The capsid exception below is separate.

The capsid exception is separate from catalogue core scoring. Every eligible reported
hit to a curated capsid, minor capsid, coat or Penton model protects its complete
gene span. A score below 10 bars new admission; eligible hits supply LOW or
MEDIUM under the reported-capsid rule above, without an E-value-only promotion. This includes semantic capsid models whose catalogue role
is support-only. Raw catalogue roles and validation remain in the evidence.
The exception alone never supplies a HIGH floor or an extra qualifying core.
The retained complete gene is sufficient for canonical acceptance, with no
minimum region length. Protease profiles such as `EVE_CAPSID_NOVEL_0001` do not
receive this exception.

### Spatial witnesses

A witness needs complete genes on one scaffold, in one panel:

| Rule | Evidence | Floor |
| --- | --- | --- |
| `reported_capsid` | Complete semantic capsid/Penton with score >=10; validated viral, score >=75 or own viral hit gives MEDIUM | LOW or MEDIUM |
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
| Integration | `tir_status` `detected`, `tsd_assessment_status` `assessed_match`, or a DDE integrase, tyrosine recombinase or serine recombinase gene inside the region; counted once | +0.05 |
| Structural | Optional structural analysis | 0.15 × structural score |
| EUK penalty | Fraction of genes with three eukaryotic hits at 70% identity or more and neither same-gene Phase1 validated viral-marker evidence nor a qualifying viral hit in the top 10 or within `phase2.viral_hit_depth` (candidate default 50), × 0.20; or, above 80% high-identity eukaryotic genes and below 10% viral genes, that fraction × 0.12. The larger term applies. Both terms exclude qualifying top-10 or deeper viral hits and same-gene Phase1 validated viral markers. This corrects the broad penalty at depth 10 too; positive evidence is unchanged. Raw top1-host diagnostics are unchanged | up to −0.12 |
| Size penalty | One, two or three genes, without full capsid credit (capsid fraction 1) | −0.15, −0.08 or −0.04 |

A family vote from taxonomy needs a dominant-family fraction of 0.25 or more, and a
marker vote needs a positive marker-complement score. The scorer also holds
taxonomy-distribution and host-divergence bonuses. The default flow does not
produce their inputs, so both stay at zero.

The integration bonus ignores ambiguous TIRs, direct repeats, GC deviation and
integrase genes in the flanks. GC deviation already counts in the composition
component. The caps below still apply, so a candidate with fewer than two cores
cannot reach HIGH through this bonus. On the 1.3.0 protist cohort, +0.05 moves
about 1,700 MEDIUM calls to HIGH, almost all of them PPV.

### Caps and floors

Caps apply first, in this order:

1. A candidate with no qualifying core and no eligible supplementary module stays
   below MEDIUM, at 0.199.
2. A candidate with fewer than two cores stays below HIGH, at 0.699.

Floors then raise the score:

| Evidence | Floor |
| --- | --- |
| Capsid core and a second core in one original panel, at any distance | HIGH, 0.7 |
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

A candidate with a contained spatial witness that supports its tier passes the
class gate with reason `hallmark_spatial_rescue_pass`. Other candidates follow
this table. The additional LOW evidence check below applies to both routes.

| Class | HIGH or MEDIUM | LOW |
| --- | --- | --- |
| `NCLDV`, `MIRUS` | Length above 5,000 bp, capsid evidence, or the short-call exception below | Length above 5,000 bp or two competitive supporting genes; plus two marker-bearing proteins or an eligible module of the same class |
| `PPV`, `MIXED` | Length above 2,000 bp, and capsid evidence or two hallmarks with one non-ATPase | Same rule |
| `CRESS` | CRESS marker family and at least one hallmark | Same rule |
| Other classes | Rejected | Rejected |

`MIXED` is an internal gate class. The hallmark count here is the number of
marker-bearing proteins, one per scoring locus. It is separate from the panel
core count above.

The current private detector experiment adds a check to every otherwise
accepted LOW candidate without a qualifying core family in
`qualifying_hallmark_genes`. It rejects a candidate supported by fewer than two
distinct genes. A qualifying core exempts the call
from this extra check; a capsid flag, model count or module flag alone does not.
The private short-call experiment also lets two competitive supporting genes
satisfy the NCLDV/MIRUS length requirement at any tier, when the candidate has
two marker-bearing proteins or an eligible module of the same class. It leaves the other
class requirements and the computed confidence unchanged. A qualifying core
alone does not grant this length exception. PPV/MIXED retain their 2,000 bp
requirement; CRESS retains its marker requirement. Search padding never supplies
protected sequence merely to meet a length cutoff.

A supporting gene lies fully inside the final interval. In its retained top ten
proteome hits after reference masking, its best viral bitscore must strictly
exceed its best cellular bitscore. Ties fail. A viral hit without a cellular hit
supports; a gene without a viral hit does not. Viral classes are NCLDV, MIRUS,
PPV (including VP and PLV), CRESS, GVMAG and EUK-pEVE. Phage-only hits do not
supply this support. A counted viral hit needs at least 25% identity. Cellular
competitors are EUK, BAC and ARC; organelle and unknown hits are neither. This
check never reads deeper hits. Existing positive marker validation is unchanged.
Repeated models on one gene
and fragments linked to one scoring locus cannot count as independent genes.

The experiment masks each background genome's whole species in both reference
searches, including pEVE-labelled references, before retaining top-ten and
top-50 host-veto evidence from the same retained Phase 2 search. A deeper native
search may change the top ten; positive support still reads only that slice.
Synthetic backgrounds determine the mask,
not their inserted viral donors. Missing or malformed gene evidence is reported
as unavailable, rather than treated as an observed zero-hit search.

The detailed table records `low_insufficient_competitive_viral_genes` or
`low_gene_evidence_unavailable` in
`normal_quality_gate_reason`. Evidence profiles report the supporting gene IDs
and assessment completeness. These checks change acceptance, not the score,
tier or coordinates. They remain a benchmark experiment pending native review.

Capsid evidence in this table is the [`has_capsid_evidence` flag](../METHODS.md#capsid-evidence).
It is separate from the capsid role in the catalogue.

## Tier rescue of NCLDV and Mirus regions

The zero-core cap and the LOW gate reject most NCLDV and Mirus regions without a
core hallmark. Cellular homologues often outrank the viral references of their
genes, so these regions can still carry several genes with a viral hit among
their top ten. After overlap selection, ViroSync applies one rescue pass to
regions whose published class (`effective_eve_class`) is `NCLDV` or `MIRUS`.

A region meets the rescue rule when:

- it is at least 5,000 bp long;
- at least two interior genes have an NCLDV hit among their top ten hits, or at
  least two have a MIRUS hit; the two counts are not added;
- `host_signature_fraction` is 0.5 or less.

A gate-rejected region that meets the rule and overlaps no kept call of any class
becomes LOW. Of overlapping rescues, the region with more viral top-ten genes
wins, then the longer region, then the smaller `eve_id`. The region becomes
MEDIUM instead when it also has an integration mark:

- `tir_status` `detected`;
- `tsd_assessment_status` `assessed_match`;
- a DDE integrase, tyrosine recombinase or serine recombinase gene inside the
  region or in its flanks;
- `gc_deviation` of 0.10 or more.

A kept LOW call that meets the rule and has a mark also becomes MEDIUM. Ambiguous
TIRs and internal direct repeats are not marks: both are common in repeat-rich
host sequence. The rescue changes the tier and acceptance, not `final_confidence`,
so a rescued MEDIUM call can have a score below 0.2.
`canonical_selection_outcome` records `tier_rescue_low`, `tier_rescue_medium`
or `tier_promoted_medium`; `normal_quality_gate_reason` keeps the gate's own
decision. The A6 gate-bypass arm skips the rescue.

## Joint assessment of adjacent seeds

One nonrecursive pass after gating can join two adjacent, disjoint candidates of
the same `NCLDV` or `MIRUS` class. Both must be rejected LOW candidates with
overlapping searched envelopes, or both must be accepted candidates separated
by at most 5,000 bp. Every gene in the joint interval must have queried taxonomy,
and an available cellular host model is required. Recorded stops, host-positive
runs, third candidates and unresolved child-bounding repeat candidates can
withhold a join. Accepted pairs also require a qualifying core hallmark in at
least one child and no repeated strict core family between them.
Eligibility and the repeat guard both require all confirmed-seed mappings, at
least one hull gene, and queried taxonomy for every hull gene. This proof permits
the DTR coverage to include the fixed hull's gene-free padding; it does not extend
the outer search windows or erase incompleteness beyond the hull.
The joined interval is scored from the start. It replaces its two children only
when it passes the gate and the supplementary module counts at least one gene
inside each child's protected seed-gene hull for rejected pairs. Search padding
does not count as that core. For accepted pairs,
each child must retain a qualifying hallmark or strict-module counted gene inside
its original interval, and the union must retain both classes and the higher
child tier. Failed unions leave the children unchanged. Independent neighbors
without distinguishing evidence can still merge.
[Bounded joint assessment](../METHODS.md#bounded-joint-assessment-of-adjacent-seeds)
gives the full rules.

## Where each decision is recorded

| Decision | Output |
| --- | --- |
| Score and tier | `final_confidence`, `confidence_tier` |
| Published class | `effective_eve_class` |
| Qualifying core families | `qualifying_hallmark_genes` |
| Normal quality-gate decision | `normal_quality_gate_reason` in the detailed table, including A6 counterfactual decisions |
| Final selection outcome and tier rescue | `canonical_selection_outcome` in the detailed table |
| Joint rescue | `joint_seed_provenance` in the detailed table |
| Components, weights, floor, witnesses and supplementary module | `score_components` and `supplementary_viral_evidence` in `evidence_profiles.json` |

### Resource-associated marker families

Marker family association uses recognized model names first, then a qualified
structure-profile class, then the bound resource annotation. A structure-source
model cannot inherit its training family from annotation text when its reference
hits do not qualify. Family provenance does not create a capsid or core-hallmark
role. Cellular and seed-blocked models retain their internal evidence and reported
association, but gain no new metadata-derived family-count credit.

Family support selects the highest-scoring ordinary, Pfam-selected hit per base
protein, with the first hit retained on score ties, then counts distinct model
names. An optional seed allowlist applies after that selection. Repeat children
use the same selected-model projection; their geometry and seed eligibility still
use all original eligible marker spans. This also corrects earlier repeat-child
counts that included multiple alternative models from one protein. Reporting
retains its per-locus count unit and resolves each hit's family once.

The annotation TSV is read beside the strictly resolved HMM resource, matching
resource fingerprinting. A TSV beside a symlink alias does not supply family
support. A missing optional TSV retains name/profile-only behavior and emits one
warning at the native run boundary.

## Confirmed cores during overlap selection

The required preservation contract applies after the normal acceptance gate and
existing unique-witness rescue. An overlap-suppressed candidate that passed those
checks must remain canonical when one of its assigned confirmed cores is absent
from the surviving candidates. Coverage requires the exact `ConfirmedSeed` tuple
(`seed_id`, scaffold, start and end), fully contained in one surviving candidate.
Coordinate overlap, a different seed identity at the same coordinates, or coverage
split across candidates does not satisfy this contract.

Restoring a candidate preserves its own coordinates, confidence, tier and assigned
cores. It does not transfer cores to another candidate, merge intervals, promote a
normal-gate reject or restore a branch whose rescued marker was excluded. The A6
gate bypass is unchanged. The selection outcome is `overlap_retained_confirmed_seed`. Candidates are
reconsidered once in the existing preference order; restored candidates add their
fully contained core identities to coverage. No subsequent pruning occurs.
