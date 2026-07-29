"""Regression test for GPU memory release between genomes.

``release_gpu_memory()`` documents itself as freeing the TMVec models, the
largest GPU consumer, and the orchestrator calls it between genomes. It used to
clear only the module cache written by ``get_cached_tmvec_predictor()``, which
nothing calls, so the cache was always ``None`` and the release freed nothing.
The predictor actually in use comes from ``get_tmvec_predictor()`` and is held
by ``TMVecDatabaseSearch``, so a batch run kept every model resident.
"""

from __future__ import annotations

import virosync.pipeline.phase3.tmvec_predictor as tmvec


class _FakePredictor:
    available = True

    def __init__(self) -> None:
        self.released = False

    def release(self) -> None:
        self.released = True


def test_release_frees_predictors_from_the_live_factory(monkeypatch) -> None:
    monkeypatch.setattr(tmvec, "TMvecPredictor", lambda **_: _FakePredictor())
    tmvec._live_predictors.clear()

    predictor = tmvec.get_tmvec_predictor(device="cpu")
    assert predictor in tmvec._live_predictors
    # The old module cache stays empty: nothing populates it.
    assert tmvec._cached_predictor is None

    tmvec.release_tmvec_predictor()
    assert predictor.released, "the predictor actually in use must be released"


def test_release_survives_a_predictor_that_raises(monkeypatch) -> None:
    class _Stuck(_FakePredictor):
        def release(self) -> None:
            raise RuntimeError("device busy")

    good = _FakePredictor()
    tmvec._live_predictors.clear()
    tmvec._live_predictors.add(_Stuck())
    tmvec._live_predictors.add(good)

    tmvec.release_tmvec_predictor()

    assert good.released, "one stuck predictor must not block the others"
