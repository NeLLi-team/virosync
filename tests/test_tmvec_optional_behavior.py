from __future__ import annotations

import json
import hashlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import torch

from virosync.pipeline.phase3 import tmvec_predictor
from virosync.pipeline.phase3.tmvec_database import TMVecDatabaseSearch
from virosync.pipeline.phase3.tmvec_predictor import TMvecPredictor
from virosync.utils import database_manager
from virosync.utils.database_manager import ViroSyncDatabaseManager
from virosync.utils.resource_installer import ResourceInstallError
from virosync.utils.resource_manifest import LEGACY_RUNTIME_RESOURCE_FILES


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_tmvec2_sequence_preparation_removes_prodigal_terminal_stop() -> None:
    assert TMvecPredictor._prepare_sequence("  MAKX*\n") == "MAKX"


@pytest.mark.parametrize("sequence", ["*", " \n*", "MA*K"])
def test_tmvec2_sequence_preparation_rejects_invalid_stops(sequence: str) -> None:
    with pytest.raises(ValueError, match="stop"):
        TMvecPredictor._prepare_sequence(sequence)


def test_tmvec2_runtime_architecture_matches_manifest_contract() -> None:
    config = tmvec_predictor.TMvecConfig()
    assert database_manager.TMVEC2_ARCHITECTURE == {
        "base_embedding_dim": config.d_model,
        "output_dim": config.out_dim,
        "nhead": config.nhead,
        "num_layers": config.num_layers,
        "dim_feedforward": config.dim_feedforward,
        "transformer_activation": config.activation,
        "projection_hidden_dim": config.projection_hidden_dim,
        "projection_activation": "relu",
        "dropout": config.dropout,
        "max_sequence_length": config.max_length,
    }


def _write_tmvec2_bundle(root: Path, monkeypatch) -> dict:
    base_dir = root / "models" / "lobster_24M"
    head_dir = root / "models" / "tmvec-2"
    bfvd_dir = root / "bfvd"
    reference_dir = root / "reference"
    for path in (base_dir, head_dir, bfvd_dir, reference_dir):
        path.mkdir(parents=True, exist_ok=True)

    base_names = [
        "config.json",
        "pytorch_model.bin",
        "special_tokens_map.json",
        "tokenizer_config.json",
        "vocab.txt",
    ]
    head_names = ["params.json", "tmvec-2.ckpt"]
    for name in base_names:
        (base_dir / name).write_bytes(f"base-{name}".encode())
    for name in head_names:
        (head_dir / name).write_bytes(f"head-{name}".encode())
    monkeypatch.setattr(
        database_manager,
        "TMVEC2_BASE_WEIGHT_SHA256",
        _sha256(base_dir / "pytorch_model.bin"),
    )
    monkeypatch.setattr(
        database_manager,
        "TMVEC2_HEAD_WEIGHT_SHA256",
        _sha256(head_dir / "tmvec-2.ckpt"),
    )

    embeddings = bfvd_dir / "bfvd_embeddings.npy"
    metadata = bfvd_dir / "bfvd_annotations.tsv"
    cpu_reference = reference_dir / "smoke_query.cpu.reference_embedding.npy"
    cuda_reference = reference_dir / "smoke_query.cuda.reference_embedding.npy"
    np.save(embeddings, np.ones((2, 512), dtype=np.float32))
    metadata.write_text(
        "id\tprotein_name\nBFVD-1\tprotein one\nBFVD-2\tprotein two\n",
        encoding="utf-8",
    )
    np.save(cpu_reference, np.ones(512, dtype=np.float32))
    np.save(cuda_reference, np.ones(512, dtype=np.float32))

    manifest = {
        "schema_version": 1,
        "bundle_version": "v1.0.0",
        "model": {
            "family": "tmvec2",
            "architecture": dict(database_manager.TMVEC2_ARCHITECTURE),
            "base": {
                "id": database_manager.TMVEC2_BASE_MODEL_ID,
                "revision": database_manager.TMVEC2_BASE_MODEL_REVISION,
                "files": [
                    {
                        "path": f"models/lobster_24M/{name}",
                        "sha256": _sha256(base_dir / name),
                    }
                    for name in base_names
                ],
            },
            "head": {
                "id": database_manager.TMVEC2_HEAD_MODEL_ID,
                "revision": database_manager.TMVEC2_HEAD_MODEL_REVISION,
                "files": [
                    {
                        "path": f"models/tmvec-2/{name}",
                        "sha256": _sha256(head_dir / name),
                    }
                    for name in head_names
                ],
            },
        },
        "databases": {
            "bfvd": {
                "attribution": {
                    "name": "BFVD",
                    "creator": "Kim, Rachel Seongeun",
                    "source_url": "https://bfvd.steineggerlab.workers.dev/",
                    "doi": "10.5281/zenodo.13993145",
                    "license": "CC BY 4.0",
                    "license_url": "https://creativecommons.org/licenses/by/4.0/",
                    "changes": (
                        "Converted BFVD protein sequences to "
                        "Lobster-24M/TMVec2 embeddings."
                    ),
                },
                "embeddings": {
                    "path": "bfvd/bfvd_embeddings.npy",
                    "sha256": _sha256(embeddings),
                    "rows": 2,
                },
                "metadata": {
                    "path": "bfvd/bfvd_annotations.tsv",
                    "sha256": _sha256(metadata),
                    "rows": 2,
                },
            }
        },
        "smoke_query": {
            "id": "bfvd-smoke",
            "sequence": "M" * 60,
            "database": "bfvd",
            "expected_target_id": "BFVD-1",
            "expected_score": 1.0,
            "score_tolerance": 0.01,
            "reference_embeddings": {
                "cpu": {
                    "path": "reference/smoke_query.cpu.reference_embedding.npy",
                    "sha256": _sha256(cpu_reference),
                    "dimensions": 512,
                    "atol": 1e-5,
                    "rtol": 1e-5,
                },
                "cuda": {
                    "path": "reference/smoke_query.cuda.reference_embedding.npy",
                    "sha256": _sha256(cuda_reference),
                    "dimensions": 512,
                    "atol": 1e-5,
                    "rtol": 1e-5,
                },
            },
        },
    }
    (root / "TMVEC_MANIFEST.json").write_text(
        json.dumps(manifest),
        encoding="utf-8",
    )
    return manifest


def test_tmvec_database_search_requires_explicit_portable_root() -> None:
    with pytest.raises(ValueError, match="database_root is required"):
        TMVecDatabaseSearch()


def test_tmvec_predictor_failure_raises_when_gpu_required(
    monkeypatch,
    tmp_path: Path,
) -> None:
    _write_tmvec2_bundle(tmp_path, monkeypatch)

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
    _write_tmvec2_bundle(tmp_path, monkeypatch)

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
    _write_tmvec2_bundle(tmp_path, monkeypatch)

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


def test_tmvec2_loads_local_checkpoint_strictly(
    monkeypatch,
    tmp_path: Path,
) -> None:
    model_dir = tmp_path / "models" / "tmvec-2"
    model_dir.mkdir(parents=True)
    config = tmvec_predictor.TMvecConfig()
    (model_dir / "params.json").write_text(
        json.dumps(
            {
                name: getattr(config, name)
                for name in config.__dataclass_fields__
            }
        ),
        encoding="utf-8",
    )
    torch.save(
        {"state_dict": {"projection.0.weight": torch.ones(1)}},
        model_dir / "tmvec-2.ckpt",
    )

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
            return self

    monkeypatch.setattr(tmvec_predictor, "TMvecModel", FakeTMvecModel)

    predictor = TMvecPredictor(device="cpu", model_root=tmp_path / "models")
    predictor._load_tmvec2()

    state_dict, strict = predictor._tmvec_model.loaded
    assert predictor._config == tmvec_predictor.TMvecConfig()
    assert strict is True
    assert list(state_dict) == ["projection.0.weight"]


def test_tmvec2_uses_lobster_attention_mask_for_padding() -> None:
    class FakeTokenizer:
        def __init__(self) -> None:
            self.call = None

        def __call__(self, sequences, **kwargs):
            self.call = (sequences, kwargs)
            return {
                "input_ids": torch.tensor([[0, 5, 6, 2]]),
                "attention_mask": torch.tensor([[1, 1, 1, 0]]),
            }

    class FakeLobsterBase:
        def __call__(self, **kwargs):
            assert kwargs["output_hidden_states"] is True
            features = torch.ones((1, 4, 408), dtype=torch.float32)
            return SimpleNamespace(hidden_states=(features,))

    class FakeHead:
        def __init__(self) -> None:
            self.padding_mask = None

        def __call__(self, features, padding_mask):
            assert features.shape == (1, 4, 408)
            self.padding_mask = padding_mask
            return torch.ones((1, 512), dtype=torch.float32)

    predictor = TMvecPredictor(device="cpu")
    predictor._initialized = True
    predictor._config = tmvec_predictor.TMvecConfig()
    predictor._tokenizer = FakeTokenizer()
    predictor._lobster_model = SimpleNamespace(model=FakeLobsterBase())
    predictor._tmvec_model = FakeHead()

    result = predictor.embed("A C D E F")

    assert result.shape == (512,)
    sequences, kwargs = predictor._tokenizer.call
    assert sequences == ["ACDEF"]
    assert kwargs == {
        "padding": True,
        "truncation": True,
        "max_length": 512,
        "return_tensors": "pt",
    }
    assert predictor._tmvec_model.padding_mask.tolist() == [
        [False, False, False, True]
    ]


def test_tmvec_database_passes_local_model_root(monkeypatch, tmp_path: Path) -> None:
    captured = {}
    events = []

    class Predictor:
        available = True

    def _predictor(**kwargs):
        events.append("predictor")
        captured.update(kwargs)
        return Predictor()

    def _manifest(cls, root, **kwargs):
        events.append("manifest")
        return {
            "databases": {
                "bfvd": {
                    "embeddings": {"path": "bfvd/bfvd_embeddings.npy"},
                    "metadata": {"path": "bfvd/bfvd_annotations.tsv"},
                }
            }
        }

    monkeypatch.setattr(
        "virosync.pipeline.phase3.tmvec_database.get_tmvec_predictor",
        _predictor,
    )
    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "load_tmvec_manifest",
        classmethod(_manifest),
    )
    searcher = TMVecDatabaseSearch(database_root=tmp_path, device="cpu")

    assert searcher.predictor is not None
    assert events == ["manifest", "predictor"]
    assert captured["model_root"] == tmp_path / "models"


def test_tmvec_database_loads_hash_bound_tsv_metadata(
    monkeypatch,
    tmp_path: Path,
) -> None:
    bfvd = tmp_path / "bfvd"
    bfvd.mkdir()
    np.save(bfvd / "bfvd_embeddings.npy", np.ones((2, 512), dtype=np.float32))
    (bfvd / "bfvd_annotations.tsv").write_text(
        "id\tprotein_name\nBFVD-1\tprotein one\nBFVD-2\tprotein two\n",
        encoding="utf-8",
    )
    manifest = {
        "databases": {
            "bfvd": {
                "embeddings": {"path": "bfvd/bfvd_embeddings.npy"},
                "metadata": {"path": "bfvd/bfvd_annotations.tsv"},
            }
        }
    }
    calls = []

    def _manifest(cls, root, **kwargs):
        calls.append((root, kwargs))
        return manifest

    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "load_tmvec_manifest",
        classmethod(_manifest),
    )
    searcher = TMVecDatabaseSearch(database_root=tmp_path, device="cpu")

    database = searcher._load_db("bfvd")

    assert database["ids"] == ["BFVD-1", "BFVD-2"]
    assert database["annotations"][0]["protein_name"] == "protein one"
    assert calls == [
        (
            tmp_path,
            {"verify_hashes": True, "databases": ["bfvd"]},
        )
    ]


def test_tmvec2_rejects_unbuilt_database_keys(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Unsupported TMVec2 database"):
        TMVecDatabaseSearch(database_root=tmp_path, databases=["cath"])



def test_phase3_honors_explicit_cpu_when_cuda_is_available(monkeypatch) -> None:
    from virosync.orchestration._flows.single_genome import phase3 as phase3_flow

    monkeypatch.setattr(torch.cuda, "is_available", lambda: True)

    assert phase3_flow._resolve_tmvec_device("cpu") == "cpu"


def test_phase3_rejects_unavailable_explicit_cuda(monkeypatch) -> None:
    from virosync.orchestration._flows.single_genome import phase3 as phase3_flow

    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)

    with pytest.raises(RuntimeError, match="CUDA.*not available"):
        phase3_flow._resolve_tmvec_device("cuda")


def test_tmvec_manifest_accepts_hash_bound_tmvec2_bundle(
    tmp_path: Path,
    monkeypatch,
) -> None:
    expected = _write_tmvec2_bundle(tmp_path, monkeypatch)

    actual = ViroSyncDatabaseManager.load_tmvec_manifest(
        tmp_path,
        verify_hashes=True,
        databases=["bfvd"],
    )

    assert actual == expected
    assert ViroSyncDatabaseManager.missing_tmvec_files(tmp_path, ["bfvd"]) == []


def test_tmvec_manifest_rejects_manifestless_legacy_arrays(tmp_path: Path) -> None:
    bfvd = tmp_path / "bfvd"
    bfvd.mkdir()
    np.save(bfvd / "bfvd_embeddings.npy", np.ones((1, 512), dtype=np.float32))
    (bfvd / "bfvd_annotations.tsv").write_text("id\nBFVD-1\n", encoding="utf-8")

    issues = ViroSyncDatabaseManager.missing_tmvec_files(tmp_path, ["bfvd"])

    assert len(issues) == 1
    assert "TMVEC_MANIFEST.json" in issues[0]


def test_tmvec_manifest_rejects_tampered_database(
    tmp_path: Path,
    monkeypatch,
) -> None:
    _write_tmvec2_bundle(tmp_path, monkeypatch)
    embeddings = tmp_path / "bfvd" / "bfvd_embeddings.npy"
    embeddings.write_bytes(embeddings.read_bytes()[:-1] + b"x")

    with pytest.raises(RuntimeError, match="checksum mismatch"):
        ViroSyncDatabaseManager.load_tmvec_manifest(
            tmp_path,
            verify_hashes=True,
            databases=["bfvd"],
        )


def test_tmvec_manifest_rejects_pickle_metadata(
    tmp_path: Path,
    monkeypatch,
) -> None:
    manifest = _write_tmvec2_bundle(tmp_path, monkeypatch)
    unsafe = tmp_path / "bfvd" / "bfvd_annotations.npy"
    np.save(unsafe, np.array([{"id": "BFVD-1"}], dtype=object), allow_pickle=True)
    manifest["databases"]["bfvd"]["metadata"] = {
        "path": "bfvd/bfvd_annotations.npy",
        "sha256": _sha256(unsafe),
        "rows": 1,
    }
    (tmp_path / "TMVEC_MANIFEST.json").write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(RuntimeError, match="must use TSV or JSONL"):
        ViroSyncDatabaseManager.load_tmvec_manifest(tmp_path, databases=["bfvd"])


def test_tmvec_manifest_rejects_duplicate_metadata_ids(
    tmp_path: Path,
    monkeypatch,
) -> None:
    manifest = _write_tmvec2_bundle(tmp_path, monkeypatch)
    metadata = tmp_path / "bfvd" / "bfvd_annotations.tsv"
    metadata.write_text(
        "id\tprotein_name\nBFVD-1\tprotein one\nBFVD-1\tduplicate\n",
        encoding="utf-8",
    )
    manifest["databases"]["bfvd"]["metadata"]["sha256"] = _sha256(metadata)
    (tmp_path / "TMVEC_MANIFEST.json").write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(RuntimeError, match="duplicate id 'BFVD-1'"):
        ViroSyncDatabaseManager.load_tmvec_manifest(tmp_path, databases=["bfvd"])


def test_tmvec_manifest_rejects_member_path_escape(
    tmp_path: Path,
    monkeypatch,
) -> None:
    manifest = _write_tmvec2_bundle(tmp_path, monkeypatch)
    manifest["databases"]["bfvd"]["metadata"]["path"] = "../outside.tsv"
    (tmp_path / "TMVEC_MANIFEST.json").write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(RuntimeError, match="stay inside"):
        ViroSyncDatabaseManager.load_tmvec_manifest(tmp_path, databases=["bfvd"])


def test_tmvec_manifest_rejects_wrong_architecture(
    tmp_path: Path,
    monkeypatch,
) -> None:
    manifest = _write_tmvec2_bundle(tmp_path, monkeypatch)
    manifest["model"]["architecture"]["base_embedding_dim"] = 1024
    (tmp_path / "TMVEC_MANIFEST.json").write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(RuntimeError, match="architecture mismatch"):
        ViroSyncDatabaseManager.load_tmvec_manifest(tmp_path, databases=["bfvd"])


def test_tmvec_setup_downloads_hash_pinned_model_files(
    tmp_path: Path,
    monkeypatch,
) -> None:
    manifest = _write_tmvec2_bundle(tmp_path, monkeypatch)
    file_payloads = {}
    for section_name in ("base", "head"):
        for item in manifest["model"][section_name]["files"]:
            path = tmp_path / item["path"]
            file_payloads[path.name] = path.read_bytes()
            path.unlink()
    downloads = []

    def fake_download(cls, source, target, **_kwargs):
        downloads.append(source)
        target.write_bytes(file_payloads[target.name])

    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "_copy_or_download_archive",
        classmethod(fake_download),
    )

    ViroSyncDatabaseManager._download_tmvec_models(tmp_path)

    assert len(downloads) == 7
    assert all("/resolve/" in source and "?download=true" in source for source in downloads)
    ViroSyncDatabaseManager.load_tmvec_manifest(
        tmp_path,
        verify_hashes=True,
        databases=["bfvd"],
    )


def test_optional_archive_checksum_fails_before_extract(tmp_path: Path) -> None:
    archive = tmp_path / "tmvec.tar.gz"
    archive.write_bytes(b"not the pinned archive")

    with pytest.raises(ResourceInstallError, match="checksum mismatch"):
        ViroSyncDatabaseManager._install_archive(
            tmp_path / "target",
            str(archive),
            archive_sha256="0" * 64,
        )


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


def test_manifestless_core_resource_checks_use_legacy_combined_hmm(
    tmp_path: Path,
) -> None:
    root = tmp_path / "virosync"
    _write_core_resource_files(root, "combined.hmm")

    assert ViroSyncDatabaseManager.default_paths(root)["hmm_db"] == root / "models" / "combined.hmm"
    assert ViroSyncDatabaseManager.default_paths(root)["marker_faa_db"] == root / "marker" / "marker.faa"
    assert ViroSyncDatabaseManager._check_missing_files(root) == []
    assert ViroSyncDatabaseManager.required_files_for_path(root) == list(
        LEGACY_RUNTIME_RESOURCE_FILES
    )


def test_core_resource_checks_reject_combined_ga_only_bundle(tmp_path: Path) -> None:
    root = tmp_path / "virosync"
    _write_core_resource_files(root, "combined_ga.hmm")

    assert ViroSyncDatabaseManager.default_paths(root)["hmm_db"] == root / "models" / "combined.hmm"
    assert ViroSyncDatabaseManager.required_files_for_path(root) == list(
        LEGACY_RUNTIME_RESOURCE_FILES
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
