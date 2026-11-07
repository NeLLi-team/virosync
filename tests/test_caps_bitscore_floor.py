from __future__ import annotations

from pathlib import Path

import pyhmmer

from virosync.pipeline.phase1.hhg_seeding import run_hmmsearch

DATA = Path(__file__).parent / "data" / "caps_floor"


def _load_caps_hmm() -> list:
    with pyhmmer.plan7.HMMFile(str(DATA / "caps.hmm")) as fh:
        return list(fh)


def test_caps_profile_has_no_ga_cutoff() -> None:
    # Capsid reporting does not require a GA gathering cutoff.
    (hmm,) = _load_caps_hmm()
    assert hmm.name.decode().startswith("plv_mcp_caps_")
    assert not hmm.cutoffs.gathering_available()


def test_capsid_reporting_retains_sub_75_hits_with_non_capsid_ga_enabled() -> None:
    hmms = _load_caps_hmm()
    faa = DATA / "caps.faa"

    no_floor = run_hmmsearch(faa, hmms, evalue_cutoff=1.0, enforce_ga_cutoffs=False)
    with_floor = run_hmmsearch(faa, hmms, evalue_cutoff=1.0, enforce_ga_cutoffs=True)

    assert any(h.score < 75 for h in no_floor)
    assert with_floor == no_floor
