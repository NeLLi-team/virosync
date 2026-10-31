"""Two acceptance invariants that a raised confidence must never violate.

1. The v2 quality gate must be monotone in tier: a region accepted at LOW has
   to stay accepted at MEDIUM and HIGH. The PPV/PLV/VP and MIXED LOW branches
   used to be strict supersets of their own HIGH/MEDIUM branch, so a 3 kb
   region with one non-ATPase hallmark and no MCP was published at LOW and
   dropped at MEDIUM -- raising its confidence removed it from the output.

2. Only exact catalogue model evidence can earn a capsid floor.
"""

from __future__ import annotations

from types import SimpleNamespace

from virosync.pipeline.hallmark_evidence import HallmarkEvidence
from virosync.pipeline.phase3.evidence_synthesizer import (
    VerificationResult,
    assign_confidence_tier,
    calculate_eve_confidence,
)
from virosync.pipeline.phase3.output_generator import evaluate_v2_quality_gate

# Label triples per class, shaped so both the LOW resolver (likely_family /
# classification) and the HIGH/MEDIUM resolver (region_classification) land on
# the same effective class.
_CLASS_LABELS = {
    "PPV": dict(region_classification="PPV", classification="PPV", likely_family="PPV"),
    "MIXED": dict(region_classification="MIXED", classification="MIXED", likely_family="UNKNOWN"),
    "NCLDV": dict(region_classification="NCLDV", classification="NCLDV", likely_family="NCLDV"),
    "MIRUS": dict(region_classification="MIRUS", classification="MIRUS", likely_family="MIRUS"),
}

# (end, hallmark_count, hallmark_genes, has_capsid_evidence)
_EVIDENCE_VECTORS = [
    # The reported non-monotone vector: 3 kb, one non-ATPase hallmark, no MCP.
    (3000, 1, ["VP_Penton_1"], False),
    # Bare MCP, below the NCLDV/MIRUS length floor.
    (3000, 1, ["plv_mcp_1"], True),
    # Two hallmarks, one non-ATPase, above every length floor.
    (6000, 2, ["VP_Penton_1", "COG0532"], False),
    # ATPase-only support.
    (3000, 2, ["PLV_PC_054", "GVOGm0760"], False),
    # Below the 2 kb floor even with an MCP.
    (1500, 3, ["plv_mcp_1", "COG0532", "VS000001"], True),
]


def _region() -> VerificationResult:
    """A 14-gene region that matches the host everywhere but one gene.

    Scored on its evidence alone this lands at 0.0: the broad-EUK and
    host-signature penalties exceed the weighted base. Only the priority-marker
    floor can lift it, which is what makes it a probe for that floor.
    """
    result = VerificationResult(eve_id="eve", scaffold="scaf", start=0, end=12000)
    result.gene_count = 14
    result.hallmark_count = 1
    result.hallmark_diversity = 1
    # One interior viral gene: enough to clear the "no viral evidence" and
    # "<5% viral" clamps, which would otherwise mask the floor.
    result.gene_taxonomy_viral_top10 = 1
    result.gene_taxonomy_viral_interior = 1
    result.genes_with_high_pident_euk = 14
    result.host_signature_gene_count = 14
    result.host_signature_fraction = 1.0
    return result


def test_v2_gate_is_monotone_in_confidence_tier() -> None:
    """Accepted at LOW must imply accepted at MEDIUM and at HIGH."""
    for eve_class, labels in _CLASS_LABELS.items():
        for end, hallmark_count, hallmark_genes, has_capsid_evidence in _EVIDENCE_VECTORS:
            decisions = {
                tier: evaluate_v2_quality_gate(
                    SimpleNamespace(
                        confidence_tier=tier,
                        start=0,
                        end=end,
                        hallmark_count=hallmark_count,
                        hallmark_genes=hallmark_genes,
                        has_capsid_evidence=has_capsid_evidence,
                        **labels,
                    )
                )
                for tier in ("LOW", "MEDIUM", "HIGH")
            }
            if not decisions["LOW"].kept:
                continue
            context = (
                f"{eve_class} len={end} hallmark={hallmark_count} genes={hallmark_genes} capsid={has_capsid_evidence}"
            )
            for tier in ("MEDIUM", "HIGH"):
                assert decisions[tier].kept, (
                    f"{context}: accepted at LOW ({decisions['LOW'].reason}) but "
                    f"rejected at {tier} ({decisions[tier].reason}) -- raising "
                    "confidence must never drop a region"
                )


def test_mcp_substring_does_not_earn_the_catalogue_floor() -> None:
    decoy = _region()
    genuine = _region()
    for result, model in ((decoy, "ncmcp_pseudoprotein"), (genuine, "GVOGm0003")):
        result.hallmark_evidence = [
            HallmarkEvidence(
                "gene|aa1-20",
                "gene",
                result.scaffold,
                result.start,
                result.end,
                model,
                1e-8,
                70.0,
                60.0,
                1,
                20,
                "unvalidated",
                0.0,
                "retained_hmm",
            )
        ]
        result.final_confidence = calculate_eve_confidence(result, crf_confidence=0.0)
    assert decoy.final_confidence < 0.2
    assert assign_confidence_tier(decoy) == "LOW"
    assert genuine.final_confidence >= 0.2
    assert assign_confidence_tier(genuine) != "LOW"
