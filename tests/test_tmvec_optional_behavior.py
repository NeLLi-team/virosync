from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import torch

from virosync.pipeline.phase3 import tmvec_predictor
from virosync.pipeline.phase3.tmvec_database import TMVecDatabaseSearch
from virosync.pipeline.phase3.tmvec_predictor import TMvecPredictor
from virosync.utils.database_manager import ViroSyncDatabaseManager
from virosync.utils.resource_manifest import RUNTIME_RESOURCE_FILES


def _write_tmvec_pair(
    root: Path,
    embeddings_name: str,
    metadata_name: str,
    rows: int = 2,
    dimensions: int = 512,
    object_metadata: bool = False,
) -> None:
    np.save(root / embeddings_name, np.ones((rows, dimensions), dtype=np.float32))
    if object_metadata:
        metadata = np.array(
            [{"id": f"protein-{index}"} for index in range(rows)],
            dtype=object,
        )
    else:
        metadata = np.array([f"protein-{index}" for index in range(rows)])
    np.save(root / metadata_name, metadata, allow_pickle=object_metadata)


def test_tmvec_database_search_requires_explicit_portable_root() -> None:
    with pytest.raises(ValueError, match="database_root is required"):
        TMVecDatabaseSearch()


def test_tmvec_predictor_failure_raises_when_gpu_required(
    monkeypatch,
    tmp_path: Path,
) -> None:
    def _fail_predictor(*args, **kwargs):
        raise RuntimeError("no CUDA device")

    monkeypatch.setattr(
        "virosync.pipeline.phase3.tmvec_database.get_tmvec_predictor",
        _fail_predictor,
    )
    searcher = TMVecDatabaseSearch(
        database_root=tmp_path,
        require_gpu=True,
    )

    with pytest.raises(RuntimeError, match="no CUDA device"):
        _ = searcher.predictor


def test_tmvec_predictor_failure_disables_when_gpu_not_required(
    monkeypatch,
    tmp_path: Path,
) -> None:
    def _fail_predictor(*args, **kwargs):
        raise RuntimeError("no CUDA device")

    monkeypatch.setattr(
        "virosync.pipeline.phase3.tmvec_database.get_tmvec_predictor",
        _fail_predictor,
    )
    searcher = TMVecDatabaseSearch(
        database_root=tmp_path,
        require_gpu=False,
    )

    assert searcher.predictor is None


def test_tmvec_predictor_failure_raises_after_preflight_enablement(
    monkeypatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        "virosync.pipeline.phase3.tmvec_database.get_tmvec_predictor",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("predictor failed")),
    )
    searcher = TMVecDatabaseSearch(
        database_root=tmp_path,
        device="cpu",
        fail_on_unavailable=True,
    )

    with pytest.raises(RuntimeError, match="predictor failed"):
        _ = searcher.predictor


def test_tmvec_batch_embedding_failure_raises_when_gpu_required(
    monkeypatch,
    tmp_path: Path,
) -> None:
    class FailingPredictor:
        available = True

        def embed_batch(self, sequences):
            raise RuntimeError("embedding failed")

    searcher = TMVecDatabaseSearch(
        database_root=tmp_path,
        require_gpu=True,
    )
    searcher._predictor = FailingPredictor()

    with pytest.raises(RuntimeError, match="embedding failed"):
        searcher.search_batch([("p1", "M" * 50)])


def test_tmvec_batch_embedding_failure_raises_after_preflight_enablement(
    tmp_path: Path,
) -> None:
    class FailingPredictor:
        available = True

        def embed_batch(self, sequences):
            raise RuntimeError("embedding failed")

    searcher = TMVecDatabaseSearch(
        database_root=tmp_path,
        device="cpu",
        fail_on_unavailable=True,
    )
    searcher._predictor = FailingPredictor()

    with pytest.raises(RuntimeError, match="embedding failed"):
        searcher.search_batch([("p1", "M" * 50)])


def test_tmvec_batch_database_load_failure_raises_after_preflight(
    tmp_path: Path,
) -> None:
    class Predictor:
        available = True

        def embed_batch(self, sequences):
            return np.ones((len(sequences), 2), dtype=np.float32)

    searcher = TMVecDatabaseSearch(
        database_root=tmp_path,
        device="cpu",
        fail_on_unavailable=True,
    )
    searcher._predictor = Predictor()
    searcher._load_db = lambda name: (_ for _ in ()).throw(OSError("DB load failed"))

    with pytest.raises(OSError, match="DB load failed"):
        searcher.search_batch([("p1", "M" * 50)], databases=["bfvd"])


def test_tmvec_batch_dimension_mismatch_raises_after_preflight(
    tmp_path: Path,
) -> None:
    class Predictor:
        available = True

        def embed_batch(self, sequences):
            return np.ones((len(sequences), 2), dtype=np.float32)

    searcher = TMVecDatabaseSearch(
        database_root=tmp_path,
        device="cpu",
        fail_on_unavailable=True,
    )
    searcher._predictor = Predictor()
    searcher._load_db = lambda name: {
        "embeddings": np.ones((2, 3), dtype=np.float32),
        "ids": ["db1", "db2"],
        "annotations": None,
    }

    with pytest.raises(ValueError, match="matmul"):
        searcher.search_batch([("p1", "M" * 50)], databases=["bfvd"])


def test_tmvec_sequence_search_failure_raises_after_preflight(
    tmp_path: Path,
) -> None:
    class Predictor:
        available = True

        def search_database(self, *args, **kwargs):
            raise RuntimeError("search failed")

    searcher = TMVecDatabaseSearch(
        database_root=tmp_path,
        device="cpu",
        fail_on_unavailable=True,
    )
    searcher._predictor = Predictor()
    searcher._load_db = lambda name: {
        "embeddings": np.ones((1, 2), dtype=np.float32),
        "ids": ["db1"],
        "annotations": None,
    }

    with pytest.raises(RuntimeError, match="search failed"):
        searcher.search_sequence("M" * 50, databases=["bfvd"])


def test_tmvec_zero_embedding_raises_when_gpu_required(tmp_path: Path) -> None:
    class ZeroPredictor:
        available = True

        def embed_batch(self, sequences):
            return np.zeros((len(sequences), 512), dtype=np.float32)

    searcher = TMVecDatabaseSearch(
        database_root=tmp_path,
        require_gpu=True,
    )
    searcher._predictor = ZeroPredictor()

    with pytest.raises(RuntimeError, match="invalid or zero embeddings"):
        searcher.search_batch([("p1", "M" * 50)])


def test_tmvec_zero_embedding_disables_hits_when_gpu_not_required(
    tmp_path: Path,
) -> None:
    class ZeroPredictor:
        available = True

        def embed_batch(self, sequences):
            return np.zeros((len(sequences), 512), dtype=np.float32)

    searcher = TMVecDatabaseSearch(
        database_root=tmp_path,
        require_gpu=False,
    )
    searcher._predictor = ZeroPredictor()

    def _unexpected_load_db(name: str):
        raise AssertionError("zero embeddings should not search databases")

    searcher._load_db = _unexpected_load_db

    assert searcher.search_batch([("p1", "M" * 50)], databases=["bfvd"]) == {
        "p1": {"bfvd": None}
    }


def test_tmvec_loads_trained_figshare_checkpoint_strictly(
    monkeypatch,
    tmp_path: Path,
) -> None:
    config_path = tmp_path / "tm_vec_swiss_model_large_params.json"
    config_path.write_text(
        """{
  "d_model": 1024,
  "nhead": 4,
  "num_layers": 4,
  "dim_feedforward": 2048,
  "out_dim": 512,
  "dropout": 0.1,
  "activation": "relu"
}"""
    )
    ckpt_path = tmp_path / "tm_vec_swiss_model_large.ckpt"
    torch.save(
        {"state_dict": {"encoder.weight": torch.ones(1), "mlp.weight": torch.ones(1)}},
        ckpt_path,
    )

    def _fake_download(self, url, filename, expected_md5=None):
        return config_path if filename.endswith("_params.json") else ckpt_path

    class FakeTMvecModel:
        def __init__(self, config):
            self.config = config
            self.loaded = None

        def load_state_dict(self, state_dict, strict):
            self.loaded = (state_dict, strict)

        def to(self, device):
            self.device = device
            return self

        def eval(self):
            self.evaluated = True

    monkeypatch.setattr(TMvecPredictor, "_download_model", _fake_download)
    monkeypatch.setattr(tmvec_predictor, "TMvecModel", FakeTMvecModel)

    predictor = TMvecPredictor(device="cpu")
    predictor._load_tmvec()

    state_dict, strict = predictor._tmvec_model.loaded
    assert predictor._config.d_model == 1024
    assert strict is True
    assert "projection.weight" in state_dict
    assert "mlp.weight" not in state_dict


def test_tmvec_rejects_incompatible_huggingface_tmvec2() -> None:
    predictor = TMvecPredictor(device="cpu", model_name="scikit-bio/tmvec-2")

    with pytest.raises(RuntimeError, match="expects 408-dimensional"):
        predictor._load_tmvec()


def test_phase3_honors_explicit_cpu_when_cuda_is_available(monkeypatch) -> None:
    from virosync.orchestration._flows.single_genome import phase3 as phase3_flow

    monkeypatch.setattr(torch.cuda, "is_available", lambda: True)

    assert phase3_flow._resolve_tmvec_device("cpu") == "cpu"


def test_phase3_rejects_unavailable_explicit_cuda(monkeypatch) -> None:
    from virosync.orchestration._flows.single_genome import phase3 as phase3_flow

    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)

    with pytest.raises(RuntimeError, match="CUDA.*not available"):
        phase3_flow._resolve_tmvec_device("cuda")


def test_legacy_untrained_bfvd_embeddings_raise_when_gpu_required(
    tmp_path: Path,
) -> None:
    root = tmp_path / "bfvd_root"
    emb_dir = root / "tmvec_embeddings"
    emb_dir.mkdir(parents=True)
    _write_tmvec_pair(
        emb_dir,
        "bfvd_embeddings.npy",
        "bfvd_annotations.npy",
        object_metadata=True,
    )
    (emb_dir / "bfvd_embeddings.log").write_text(
        "Using TMvec with default ProtT5 configuration\n"
    )

    searcher = TMVecDatabaseSearch(database_root=root, require_gpu=True)

    with pytest.raises(RuntimeError, match="legacy untrained/random-weight"):
        searcher._load_db("bfvd")


def test_legacy_untrained_bfvd_embeddings_disable_when_gpu_not_required(
    tmp_path: Path,
) -> None:
    root = tmp_path / "bfvd_root"
    emb_dir = root / "tmvec_embeddings"
    emb_dir.mkdir(parents=True)
    _write_tmvec_pair(
        emb_dir,
        "bfvd_embeddings.npy",
        "bfvd_annotations.npy",
        object_metadata=True,
    )
    (emb_dir / "stdout.log").write_text("Using TMvec with random initialization\n")

    searcher = TMVecDatabaseSearch(database_root=root, require_gpu=False)

    assert searcher._load_db("bfvd") is None


def test_missing_tmvec_files_accepts_supported_flat_and_nested_layouts(
    tmp_path: Path,
) -> None:
    flat_bfvd = tmp_path / "flat_bfvd"
    flat_bfvd.mkdir()
    _write_tmvec_pair(
        flat_bfvd,
        "bfvd_embeddings.npy",
        "bfvd_annotations.npy",
        object_metadata=True,
    )

    nested_bfvd = tmp_path / "nested_bfvd"
    (nested_bfvd / "tmvec_embeddings").mkdir(parents=True)
    _write_tmvec_pair(
        nested_bfvd / "tmvec_embeddings",
        "bfvd_embeddings.npy",
        "bfvd_annotations.npy",
        object_metadata=True,
    )

    flat_pdb = tmp_path / "flat_pdb"
    flat_pdb.mkdir()
    _write_tmvec_pair(flat_pdb, "embeddings.npy", "metadata.npy")

    assert ViroSyncDatabaseManager.missing_tmvec_files(flat_bfvd, ["bfvd"]) == []
    assert ViroSyncDatabaseManager.missing_tmvec_files(nested_bfvd, ["bfvd"]) == []
    assert ViroSyncDatabaseManager.missing_tmvec_files(flat_pdb, ["pdb"]) == []


@pytest.mark.parametrize("invalid_kind", ["empty", "directory"])
def test_missing_tmvec_files_rejects_non_regular_or_empty_assets(
    tmp_path: Path,
    invalid_kind: str,
) -> None:
    root = tmp_path / invalid_kind
    root.mkdir()
    invalid = root / "bfvd_embeddings.npy"
    if invalid_kind == "directory":
        invalid.mkdir()
    else:
        invalid.write_bytes(b"")
    np.save(
        root / "bfvd_annotations.npy",
        np.array([{"id": "valid"}], dtype=object),
        allow_pickle=True,
    )

    missing = ViroSyncDatabaseManager.missing_tmvec_files(root, ["bfvd"])

    assert len(missing) == 1
    assert "non-empty regular file" in missing[0]


@pytest.mark.parametrize(
    ("case", "message"),
    [
        ("corrupt_header", "valid NPY"),
        ("one_dimensional", "2-D"),
        ("empty_shape", "non-empty"),
        ("nonnumeric", "numeric"),
        ("wrong_width", "512"),
        ("row_mismatch", "row count"),
        ("truncated", "truncated"),
    ],
)
def test_tmvec_preflight_rejects_invalid_npy_assets(
    tmp_path: Path,
    case: str,
    message: str,
) -> None:
    root = tmp_path / case
    root.mkdir()
    _write_tmvec_pair(
        root,
        "bfvd_embeddings.npy",
        "bfvd_annotations.npy",
        object_metadata=True,
    )
    embeddings = root / "bfvd_embeddings.npy"
    if case == "corrupt_header":
        embeddings.write_bytes(b"not-npy")
    elif case == "one_dimensional":
        np.save(embeddings, np.ones(3, dtype=np.float32))
    elif case == "empty_shape":
        np.save(embeddings, np.ones((0, 3), dtype=np.float32))
    elif case == "nonnumeric":
        np.save(embeddings, np.array([["not", "numeric"]]))
    elif case == "wrong_width":
        np.save(embeddings, np.ones((2, 3), dtype=np.float32))
    elif case == "row_mismatch":
        np.save(
            root / "bfvd_annotations.npy",
            np.array([{"id": "only-one"}], dtype=object),
            allow_pickle=True,
        )
    elif case == "truncated":
        embeddings.write_bytes(embeddings.read_bytes()[:-1])

    missing = ViroSyncDatabaseManager.missing_tmvec_files(root, ["bfvd"])

    assert len(missing) == 1
    assert message in missing[0]


def test_tmvec_preflight_does_not_load_object_metadata_payload(
    tmp_path: Path,
    monkeypatch,
) -> None:
    _write_tmvec_pair(
        tmp_path,
        "bfvd_embeddings.npy",
        "bfvd_annotations.npy",
        object_metadata=True,
    )
    monkeypatch.setattr(np, "load", lambda *args, **kwargs: pytest.fail("payload loaded"))

    assert ViroSyncDatabaseManager.missing_tmvec_files(tmp_path, ["bfvd"]) == []


def test_tmvec_preflight_rejects_object_metadata_truncated_by_one_byte(
    tmp_path: Path,
) -> None:
    _write_tmvec_pair(
        tmp_path,
        "bfvd_embeddings.npy",
        "bfvd_annotations.npy",
        object_metadata=True,
    )
    metadata = tmp_path / "bfvd_annotations.npy"
    metadata.write_bytes(metadata.read_bytes()[:-1])

    missing = ViroSyncDatabaseManager.missing_tmvec_files(tmp_path, ["bfvd"])

    assert len(missing) == 1
    assert "pickle payload" in missing[0]
    assert "STOP" in missing[0]


def test_tmvec_preflight_rejects_object_metadata_trailing_data(
    tmp_path: Path,
) -> None:
    _write_tmvec_pair(
        tmp_path,
        "bfvd_embeddings.npy",
        "bfvd_annotations.npy",
        object_metadata=True,
    )
    metadata = tmp_path / "bfvd_annotations.npy"
    metadata.write_bytes(metadata.read_bytes() + b"unexpected")

    missing = ViroSyncDatabaseManager.missing_tmvec_files(tmp_path, ["bfvd"])

    assert len(missing) == 1
    assert "unexpected trailing data" in missing[0]


def _write_core_resource_files(root: Path, hmm_name: str) -> None:
    for rel_path in [
        "DB_VERSION",
        "DATABASE_README.txt",
        f"models/{hmm_name}",
        f"models/{hmm_name}.h3f",
        f"models/{hmm_name}.h3i",
        f"models/{hmm_name}.h3m",
        f"models/{hmm_name}.h3p",
        "models/model_annotations_with_interpro.tsv",
        "models/og_marker_name_map.tsv",
        "marker/marker.dmnd",
        "marker/marker.faa",
        "genomes/combined_proteome.dmnd",
        "taxonomy/labels.tsv",
    ]:
        path = root / rel_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("synthetic resource payload\n")


def test_core_resource_checks_prefer_combined_hmm(tmp_path: Path) -> None:
    root = tmp_path / "virosync"
    _write_core_resource_files(root, "combined.hmm")

    assert ViroSyncDatabaseManager.default_paths(root)["hmm_db"] == root / "models" / "combined.hmm"
    assert ViroSyncDatabaseManager.default_paths(root)["marker_faa_db"] == root / "marker" / "marker.faa"
    assert ViroSyncDatabaseManager._check_missing_files(root) == []
    assert ViroSyncDatabaseManager.required_files_for_path(root) == list(
        RUNTIME_RESOURCE_FILES
    )


def test_core_resource_checks_reject_combined_ga_only_bundle(tmp_path: Path) -> None:
    root = tmp_path / "virosync"
    _write_core_resource_files(root, "combined_ga.hmm")

    assert ViroSyncDatabaseManager.default_paths(root)["hmm_db"] == root / "models" / "combined.hmm"
    assert ViroSyncDatabaseManager.required_files_for_path(root) == list(
        RUNTIME_RESOURCE_FILES
    )
    assert ViroSyncDatabaseManager._check_missing_files(root) == [
        "models/combined.hmm"
    ]


def test_core_resource_checks_report_missing_hmm_group(tmp_path: Path) -> None:
    root = tmp_path / "virosync"
    for rel_path in ViroSyncDatabaseManager.REQUIRED_FILES:
        path = root / rel_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("synthetic resource payload\n")

    missing = ViroSyncDatabaseManager._check_missing_files(root)

    assert missing == ["models/combined.hmm"]


def test_runtime_resource_checks_do_not_require_source_repair_payloads(
    tmp_path: Path,
) -> None:
    root = tmp_path / "virosync"
    _write_core_resource_files(root, "combined.hmm")
    (root / "marker/marker.faa").unlink()
    for suffix in (".h3f", ".h3i", ".h3m", ".h3p"):
        (root / f"models/combined.hmm{suffix}").unlink()

    assert ViroSyncDatabaseManager._check_missing_files(root) == []
    assert ViroSyncDatabaseManager.default_paths(root)["marker_faa_db"] is None
