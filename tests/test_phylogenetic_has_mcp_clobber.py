"""Regression test: phylogenetic validation must not clobber a true has_mcp.

Past bug: `_run_phylogenetic_validation` in evidence_synthesizer.py
unconditionally assigned ``result.has_mcp = phylo_result.has_mcp``. When an
EVE's MCP was detected upstream by the HMM/Diamond path
(``result.has_mcp=True``) but GVClass happened to return ``has_mcp=False``
(different marker corpus, different region resolution), the True value was
silently overwritten to False. That flag then flows into MEDIUM-tier
promotion, LOW-filter retention, the family-consistency floor, and the v2
export-gate MCP-exception — so the bug caused real MCP-bearing EVEs to be
dropped from exports.

Fix: ``result.has_mcp = result.has_mcp or phylo_result.has_mcp`` — the flag
can be promoted but never demoted by phylogenetic validation.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.pipeline.phase3.evidence_synthesizer import (
    EvidenceSynthesizer,
    EvidenceSynthesizerConfig,
    VerificationResult,
)
from virosync.pipeline.phase3.phylogenetic_validation import (
    PhylogeneticValidationResult,
    PhylogeneticVerdict,
)


def _build_synth(tmp_path: Path) -> EvidenceSynthesizer:
    config = EvidenceSynthesizerConfig()
    config.use_phylogenetic_validation = True
    synth = EvidenceSynthesizer(
        config=config,
        genome_path=tmp_path / "fake.fna",  # property wants a truthy path
    )
    # Inject mock validator; the property's lazy init is skipped because
    # _phylogenetic_validator is not None after assignment.
    synth._phylogenetic_validator = MagicMock()
    return synth


def _phylo_result(has_mcp: bool) -> PhylogeneticValidationResult:
    return PhylogeneticValidationResult(
        eve_id="EVE_test",
        scaffold="ctg1",
        start=0,
        end=1000,
        verdict=PhylogeneticVerdict.AMBIGUOUS,
        combined_score=0.5,
        has_mcp=has_mcp,
        has_viral_markers=False,
        is_chimeric=False,
    )


def _boundary() -> RefinedBoundary:
    return RefinedBoundary(
        scaffold="ctg1",
        start=0,
        end=1000,
        seed_id="s",
        original_start=0,
        original_end=1000,
        confidence=0.5,
        posterior_probability=0.5,
    )


def _verification_result(has_mcp: bool) -> VerificationResult:
    return VerificationResult(
        eve_id="EVE_test",
        scaffold="ctg1",
        start=0,
        end=1000,
        has_mcp=has_mcp,
    )


def test_phylo_validation_does_not_clobber_true_has_mcp(tmp_path: Path) -> None:
    synth = _build_synth(tmp_path)
    synth._phylogenetic_validator.validate_eve.return_value = _phylo_result(has_mcp=False)

    result = _verification_result(has_mcp=True)
    synth._run_phylogenetic_validation(result, _boundary())

    assert result.has_mcp is True, (
        "has_mcp set upstream by HMM/Diamond must not be demoted by GVClass"
    )


def test_phylo_validation_promotes_false_has_mcp(tmp_path: Path) -> None:
    synth = _build_synth(tmp_path)
    synth._phylogenetic_validator.validate_eve.return_value = _phylo_result(has_mcp=True)

    result = _verification_result(has_mcp=False)
    synth._run_phylogenetic_validation(result, _boundary())

    assert result.has_mcp is True, (
        "GVClass-detected MCP must promote an upstream-False has_mcp to True"
    )


def test_phylo_validation_keeps_false_when_both_false(tmp_path: Path) -> None:
    synth = _build_synth(tmp_path)
    synth._phylogenetic_validator.validate_eve.return_value = _phylo_result(has_mcp=False)

    result = _verification_result(has_mcp=False)
    synth._run_phylogenetic_validation(result, _boundary())

    assert result.has_mcp is False
