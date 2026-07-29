from __future__ import annotations

from pathlib import Path

import pyhmmer

from virosync.pipeline.phase1.hhg_seeding import (
    CAPS_BITSCORE_FLOOR,
    run_hmmsearch,
)

DATA = Path(__file__).parent / "data" / "caps_floor"


def _load_caps_hmm() -> list:
    with pyhmmer.plan7.HMMFile(str(DATA / "caps.hmm")) as fh:
        return list(fh)


def test_caps_profile_has_no_ga_cutoff() -> None:
    # The capscan MCP profiles ship without a GA gathering cutoff, which is exactly
    # why the bitscore floor is applied by name (_caps_) rather than via GA.
    (hmm,) = _load_caps_hmm()
    assert hmm.name.decode().startswith("plv_mcp_caps_")
    assert not hmm.cutoffs.gathering_available()


def test_caps_floor_drops_sub_threshold_hits() -> None:
    hmms = _load_caps_hmm()
    faa = DATA / "caps.faa"

    no_floor = run_hmmsearch(faa, hmms, evalue_cutoff=1.0, enforce_ga_cutoffs=False)
    with_floor = run_hmmsearch(faa, hmms, evalue_cutoff=1.0, enforce_ga_cutoffs=True)

    # fixture has 4 sub-75 + 4 strong hits; the floor must remove the weak ones.
    assert any(h.score < CAPS_BITSCORE_FLOOR for h in no_floor)
    assert len(with_floor) < len(no_floor)
    # every surviving caps hit meets capscan's strong-hit threshold (seq AND domain).
    for h in with_floor:
        assert h.score >= CAPS_BITSCORE_FLOOR
        assert h.domain_score >= CAPS_BITSCORE_FLOOR
