from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import click
import pytest

from virosync.config import (
    ApplicationConfig,
    ConfigError,
    MaskingBackend,
    MaskingConfig,
    MaskingFailurePolicy,
)
from virosync.orchestration._flows.single_genome.manifest import (
    _write_completion_manifest,
)
from virosync.orchestration._flows.single_genome.resume import (
    _manifest_is_stale,
    _valid_completion_manifest,
)
from virosync.orchestration.cli import _build_pipeline_config
import virosync.orchestration._flows.single_genome.orchestrator as orchestrator
import virosync.pipeline.phase0.masking as masking
import virosync.utils.provenance as provenance


def _fasta(path: Path, sequence: str = "ACGTACGT", seq_id: str = "demo") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">{seq_id}\n{sequence}\n")
    return path


def _reseal_status_payload(payload: dict) -> dict:
    semantic = dict(payload)
    semantic.pop("result_fingerprint", None)
    payload["result_fingerprint"] = hashlib.sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return payload


def _species_config(
    backend: MaskingBackend = MaskingBackend.TRF_REPEATMASKER,
    *,
    policy: MaskingFailurePolicy = MaskingFailurePolicy.STRICT,
    fallback: MaskingBackend | None = None,
) -> MaskingConfig:
    return MaskingConfig(
        backend=backend,
        failure_policy=policy,
        fallback_backend=fallback,
        repeatmasker_species="chlorophyta",
    )


def test_masking_config_validates_target_xor_and_fallback() -> None:
    assert not MaskingConfig(backend=MaskingBackend.OFF).validate()
    assert not MaskingConfig(backend=MaskingBackend.TRF).validate()
    assert not MaskingConfig(
        backend=MaskingBackend.OFF,
        repeatmasker_species="latent-target",
    ).validate()
    assert MaskingConfig(backend=MaskingBackend.REPEATMASKER).validate()
    assert MaskingConfig(
        backend=MaskingBackend.REPEATMASKER,
        repeatmasker_species="x",
        repeatmasker_library=Path("library.fa"),
    ).validate()
    assert MaskingConfig(
        backend=MaskingBackend.TRF,
        failure_policy=MaskingFailurePolicy.FALLBACK,
        fallback_backend=MaskingBackend.REPEATMASKER,
    ).validate()
    assert not MaskingConfig(
        backend=MaskingBackend.TRF_REPEATMASKER,
        failure_policy=MaskingFailurePolicy.FALLBACK,
        fallback_backend=MaskingBackend.TRF,
        repeatmasker_species="x",
    ).validate()


def test_shipped_configs_are_off_without_an_implicit_repeatmasker_target() -> None:
    for path in (
        Path("config/orchestration.yaml"),
        Path("config/orchestration_archaeal.yaml"),
    ):
        config = ApplicationConfig.from_yaml(path).pipeline.execution.masking
        assert config.backend is MaskingBackend.OFF
        assert config.repeatmasker_species is None
        assert config.repeatmasker_library is None


def test_skip_compatibility_normalizes_fallback_request_to_strict_off() -> None:
    requested = _species_config(
        policy=MaskingFailurePolicy.FALLBACK,
        fallback=MaskingBackend.TRF,
    )
    effective = requested.with_backend(MaskingBackend.OFF)
    assert effective.backend is MaskingBackend.OFF
    assert effective.failure_policy is MaskingFailurePolicy.STRICT
    assert effective.fallback_backend is None
    assert effective.repeatmasker_species == requested.repeatmasker_species
    assert not effective.validate()


def test_cli_no_skip_requires_target_but_trf_does_not() -> None:
    default = ApplicationConfig.from_dict({"schema_version": 1})
    with pytest.raises(click.ClickException, match="repeatmasker_species"):
        _build_pipeline_config(
            yaml_config=default,
            clean_run=False,
            skip_masking=False,
        )

    latent = ApplicationConfig.from_dict(
        {
            "schema_version": 1,
            "execution": {
                "masking": {
                    "backend": "off",
                    "failure_policy": "strict",
                    "fallback_backend": None,
                    "repeatmasker_species": "chlorophyta",
                    "repeatmasker_library": None,
                }
            },
        }
    )
    enabled = _build_pipeline_config(
        yaml_config=latent,
        clean_run=False,
        skip_masking=False,
    )
    assert enabled.execution.masking.backend is MaskingBackend.TRF_REPEATMASKER

    trf = ApplicationConfig.from_dict(
        {
            "schema_version": 1,
            "execution": {"masking": {"backend": "trf"}},
        }
    )
    assert not trf.pipeline.execution.masking.validate()


def test_repeatmasker_argv_uses_exactly_one_target_and_wires_engine(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    calls: list[list[str]] = []

    def _run(command, **_kwargs):
        calls.append(list(command))
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(masking.subprocess, "run", _run)

    species_dir = tmp_path / "species"
    _fasta(species_dir / "input.fna.masked")
    masking.run_repeatmasker(
        input_fasta,
        species_dir,
        species="chlorophyta",
        engine="rmblast",
    )
    species_argv = calls.pop()
    assert species_argv[:5] == [
        "RepeatMasker",
        "-e",
        "rmblast",
        "-species",
        "chlorophyta",
    ]
    assert "-lib" not in species_argv

    library = tmp_path / "custom.fa"
    library.write_text(">repeat\nACGT\n")
    library_dir = tmp_path / "library"
    _fasta(library_dir / "input.fna.masked")
    masking.run_repeatmasker(input_fasta, library_dir, library=library)
    library_argv = calls.pop()
    assert library_argv[3:5] == ["-lib", str(library)]
    assert "-species" not in library_argv

    with pytest.raises(ConfigError, match="exactly one"):
        masking.run_repeatmasker(input_fasta, tmp_path / "neither")
    with pytest.raises(ConfigError, match="exactly one"):
        masking.run_repeatmasker(
            input_fasta,
            tmp_path / "both",
            species="x",
            library=library,
        )
    assert calls == []


def test_backend_version_probe_uses_verified_dash_v(
    monkeypatch,
) -> None:
    calls: list[list[str]] = []

    def _run(command, **_kwargs):
        calls.append(list(command))
        text = (
            "Tandem Repeats Finder, Version 4.10.0-rc.2"
            if command[0] == "trf"
            else "RepeatMasker version 4.2.2"
        )
        return SimpleNamespace(returncode=0, stdout=text + "\n", stderr="")

    monkeypatch.setattr(masking.subprocess, "run", _run)
    assert "4.10.0" in masking._probe_backend_version(MaskingBackend.TRF)
    assert "4.2.2" in masking._probe_backend_version(MaskingBackend.REPEATMASKER)
    assert calls == [["trf", "-v"], ["RepeatMasker", "-v"]]


def test_apply_mask_counts_only_new_ambiguous_bases(tmp_path: Path) -> None:
    input_fasta = _fasta(tmp_path / "input.fna", "ANnT")
    output_fasta = tmp_path / "masked.fna"
    count = masking.apply_mask(
        input_fasta,
        output_fasta,
        [masking.MaskedRegion("demo", 0, 4, "trf")],
    )
    assert count == 2
    assert "NNNN" in output_fasta.read_text()


def test_off_mode_spawns_no_tools_and_writes_valid_status(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    monkeypatch.setattr(
        masking.subprocess,
        "run",
        lambda *_args, **_kwargs: pytest.fail("off mode spawned a tool"),
    )
    result = masking.mask_genome_pipeline(
        input_fasta,
        tmp_path / "masking",
        config=MaskingConfig(),
    )
    assert result.status == "off"
    assert result.output_path == input_fasta
    assert result.masked_bases == 0
    assert result.benchmark_eligible
    assert masking.load_masking_result(
        result.status_path,
        expected_config=MaskingConfig(),
        expected_input=input_fasta,
    ) == result


def test_self_fingerprinted_illegal_off_state_is_ineligible_and_rejected(
    tmp_path: Path,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    result = masking.mask_genome_pipeline(
        input_fasta,
        tmp_path / "masking",
        config=MaskingConfig(),
    )
    illegal_result = replace(
        result,
        requested_backend=MaskingBackend.TRF,
        effective_backend=MaskingBackend.TRF,
    )
    assert not illegal_result.benchmark_eligible

    payload = json.loads(result.status_path.read_text())
    payload["requested_backend"] = "trf"
    payload["effective_backend"] = "trf"
    payload["benchmark_eligible"] = False
    result.status_path.write_text(
        json.dumps(_reseal_status_payload(payload), indent=2, sort_keys=True) + "\n"
    )
    with pytest.raises(ValueError, match="invalid masking state.*off status"):
        masking.load_masking_result(result.status_path)


def test_phase0_rejects_valid_off_result_for_requested_trf(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    output_dir = tmp_path / "result"
    off_result = masking.mask_genome_pipeline(
        input_fasta,
        output_dir / "phase0" / "masking",
        config=MaskingConfig(),
    )

    def _call_task(task, **_kwargs):
        if task is orchestrator.mask_genome_task:
            return off_result
        pytest.fail("Phase 0 continued to Prodigal after a masking request mismatch")

    monkeypatch.setattr(orchestrator, "call_task", _call_task)
    with pytest.raises(ValueError, match="masking status mismatch.*requested backend"):
        orchestrator._run_phase0_subflow(
            genome_path=input_fasta,
            output_dir=output_dir,
            genome_id="demo",
            threads=1,
            skip_masking=None,
            resume=False,
            logger=SimpleNamespace(info=lambda *_args, **_kwargs: None),
            masking=MaskingConfig(backend=MaskingBackend.TRF),
        )


def test_phase0_accepts_semantically_equal_relative_masking_result(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.chdir(tmp_path)
    input_fasta = _fasta(Path("input.fna"))
    output_dir = Path("result")

    def _call_task(task, **kwargs):
        if task is orchestrator.mask_genome_task:
            return masking.mask_genome_pipeline(
                kwargs["genome_path"],
                kwargs["output_dir"] / "masking",
                config=kwargs["masking"],
            )
        if task is orchestrator.generate_proteome_task:
            proteome = kwargs["output_dir"] / "proteome.fasta"
            proteome.write_text("")
            return proteome, 0
        raise AssertionError(f"unexpected task: {task}")

    monkeypatch.setattr(orchestrator, "call_task", _call_task)
    phase0 = orchestrator._run_phase0_subflow(
        genome_path=input_fasta,
        output_dir=output_dir,
        genome_id="demo",
        threads=1,
        skip_masking=None,
        resume=False,
        logger=SimpleNamespace(info=lambda *_args, **_kwargs: None),
        masking=MaskingConfig(),
    )
    assert phase0["masking_result"].status == "off"


def test_relative_custom_library_identity_round_trips(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.chdir(tmp_path)
    input_fasta = _fasta(Path("input.fna"))
    library = Path("repeats.fa")
    library.write_text(">repeat\nACGT\n")
    config = MaskingConfig(
        backend=MaskingBackend.REPEATMASKER,
        repeatmasker_library=library,
    )
    monkeypatch.setattr(
        masking,
        "_run_backend_attempt",
        lambda **_kwargs: ([], "RepeatMasker version"),
    )
    result = masking.mask_genome_pipeline(
        input_fasta,
        Path("masking"),
        config=config,
    )
    loaded = masking.load_masking_result(
        result.status_path,
        repeat_regions=result.repeat_regions,
        expected_config=config,
        expected_input=input_fasta,
    )
    assert loaded.to_status_payload() == result.to_status_payload()


def test_combined_success_records_versions_and_is_benchmark_eligible(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")

    def _attempt(*, backend, **_kwargs):
        region = masking.MaskedRegion(
            "demo",
            0 if backend is MaskingBackend.TRF else 4,
            2 if backend is MaskingBackend.TRF else 6,
            backend.value,
        )
        return [region], f"{backend.value} version"

    monkeypatch.setattr(masking, "_run_backend_attempt", _attempt)
    result = masking.mask_genome_pipeline(
        input_fasta,
        tmp_path / "masking",
        config=_species_config(),
    )
    assert result.status == "success"
    assert result.masked_bases == 4
    assert dict(result.backend_versions) == {
        "repeatmasker": "repeatmasker version",
        "trf": "trf version",
    }
    assert result.benchmark_eligible


def test_atomic_mask_and_status_writes_replace_symlinks_without_touching_victims(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    output_dir = tmp_path / "masking"
    output_dir.mkdir()
    final_victim = tmp_path / "final-victim.txt"
    temporary_victim = tmp_path / "temporary-victim.txt"
    status_victim = tmp_path / "status-victim.txt"
    final_victim.write_text("final sentinel\n")
    temporary_victim.write_text("temporary sentinel\n")
    status_victim.write_text("status sentinel\n")
    (output_dir / "genome.masked.fna").symlink_to(final_victim)
    (output_dir / "genome.masked.fna.tmp").symlink_to(temporary_victim)
    (output_dir / "masking_status.json").symlink_to(status_victim)
    monkeypatch.setattr(
        masking,
        "_run_backend_attempt",
        lambda **_kwargs: (
            [masking.MaskedRegion("demo", 0, 2, "trf")],
            "trf version",
        ),
    )

    result = masking.mask_genome_pipeline(
        input_fasta,
        output_dir,
        config=MaskingConfig(backend=MaskingBackend.TRF),
    )
    assert final_victim.read_text() == "final sentinel\n"
    assert temporary_victim.read_text() == "temporary sentinel\n"
    assert status_victim.read_text() == "status sentinel\n"
    assert not result.output_path.is_symlink()
    assert result.output_path.is_file()
    assert not result.status_path.is_symlink()
    assert result.status_path.is_file()


def test_success_with_configured_fallback_records_no_selected_fallback(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    config = MaskingConfig(
        backend=MaskingBackend.TRF,
        failure_policy=MaskingFailurePolicy.FALLBACK,
        fallback_backend=MaskingBackend.OFF,
    )
    monkeypatch.setattr(
        masking,
        "_run_backend_attempt",
        lambda **_kwargs: (
            [masking.MaskedRegion("demo", 0, 2, "trf")],
            "trf version",
        ),
    )
    result = masking.mask_genome_pipeline(
        input_fasta,
        tmp_path / "masking",
        config=config,
    )
    assert result.status == "success"
    assert result.configured_fallback_backend is MaskingBackend.OFF
    assert result.fallback_backend is None
    loaded = masking.load_masking_result(
        result.status_path,
        repeat_regions=result.repeat_regions,
        expected_config=config,
        expected_input=input_fasta,
    )
    assert loaded == result


def test_strict_backend_failure_is_fatal_and_records_failed_status(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")

    def _fail(**_kwargs):
        raise masking.MaskingBackendError(MaskingBackend.TRF, "synthetic failure")

    monkeypatch.setattr(masking, "_run_backend_attempt", _fail)
    output_dir = tmp_path / "masking"
    with pytest.raises(masking.MaskingBackendError, match="synthetic failure"):
        masking.mask_genome_pipeline(
            input_fasta,
            output_dir,
            config=MaskingConfig(backend=MaskingBackend.TRF),
        )
    payload = json.loads((output_dir / "masking_status.json").read_text())
    assert payload["status"] == "failed"
    assert payload["benchmark_eligible"] is False
    assert payload["output_path"] is None


def test_programming_error_escapes_without_fallback_or_success_status(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    calls = 0

    def _bug(**_kwargs):
        nonlocal calls
        calls += 1
        raise TypeError("programming bug")

    monkeypatch.setattr(masking, "_run_backend_attempt", _bug)
    output_dir = tmp_path / "masking"
    with pytest.raises(TypeError, match="programming bug"):
        masking.mask_genome_pipeline(
            input_fasta,
            output_dir,
            config=MaskingConfig(
                backend=MaskingBackend.TRF,
                failure_policy=MaskingFailurePolicy.FALLBACK,
                fallback_backend=MaskingBackend.OFF,
            ),
        )
    assert calls == 1
    assert not (output_dir / "masking_status.json").exists()


def test_duplicate_input_ids_fail_before_backend_and_cannot_fallback_off(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = tmp_path / "duplicate.fna"
    input_fasta.write_text(">dup\nACGT\n>dup\nTGCA\n")
    calls = 0

    def _attempt(**_kwargs):
        nonlocal calls
        calls += 1
        raise AssertionError("backend must not run")

    monkeypatch.setattr(masking, "_run_backend_attempt", _attempt)
    output_dir = tmp_path / "masking"
    with pytest.raises(ValueError, match="duplicate record IDs.*dup"):
        masking.mask_genome_pipeline(
            input_fasta,
            output_dir,
            config=MaskingConfig(
                backend=MaskingBackend.TRF,
                failure_policy=MaskingFailurePolicy.FALLBACK,
                fallback_backend=MaskingBackend.OFF,
            ),
        )
    assert calls == 0
    assert not (output_dir / "masking_status.json").exists()


def test_subprocess_permission_error_uses_configured_off_fallback(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    monkeypatch.setattr(
        masking.subprocess,
        "run",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            PermissionError("execution denied")
        ),
    )
    result = masking.mask_genome_pipeline(
        input_fasta,
        tmp_path / "masking",
        config=MaskingConfig(
            backend=MaskingBackend.TRF,
            failure_policy=MaskingFailurePolicy.FALLBACK,
            fallback_backend=MaskingBackend.OFF,
        ),
    )
    assert result.status == "fallback"
    assert result.effective_backend is MaskingBackend.OFF
    assert result.fallback_backend is MaskingBackend.OFF
    assert "execution denied" in result.fallback_reason


@pytest.mark.parametrize(
    "parser_error",
    [OSError("unreadable output"), ValueError("malformed output")],
)
def test_backend_output_parse_failures_use_configured_off_fallback(
    tmp_path: Path,
    monkeypatch,
    parser_error: Exception,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    trf_dat = tmp_path / "synthetic.dat"
    trf_dat.write_text("Sequence: demo\n")
    monkeypatch.setattr(masking, "run_trf", lambda *_args, **_kwargs: trf_dat)

    def _parse(_path):
        raise parser_error

    monkeypatch.setattr(masking, "parse_trf_output", _parse)
    result = masking.mask_genome_pipeline(
        input_fasta,
        tmp_path / "masking",
        config=MaskingConfig(
            backend=MaskingBackend.TRF,
            failure_policy=MaskingFailurePolicy.FALLBACK,
            fallback_backend=MaskingBackend.OFF,
        ),
    )
    assert result.status == "fallback"
    assert result.effective_backend is MaskingBackend.OFF
    assert str(parser_error) in result.fallback_reason


def test_truncated_trf_header_is_a_typed_backend_failure(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    trf_dat = tmp_path / "truncated.dat"
    trf_dat.write_text("Sequence:\n")
    monkeypatch.setattr(masking, "run_trf", lambda *_args, **_kwargs: trf_dat)
    result = masking.mask_genome_pipeline(
        input_fasta,
        tmp_path / "masking",
        config=MaskingConfig(
            backend=MaskingBackend.TRF,
            failure_policy=MaskingFailurePolicy.FALLBACK,
            fallback_backend=MaskingBackend.OFF,
        ),
    )
    assert result.status == "fallback"
    assert "malformed TRF sequence header" in result.fallback_reason


@pytest.mark.parametrize(
    ("contents", "reason"),
    [
        ("", "empty"),
        ("Parameters: 2 7 7 80 10 50 500\n", "headers"),
        ("Sequence: demo\n1 not-a-coordinate\n", "coordinates"),
    ],
)
def test_invalid_trf_artifacts_cannot_be_reported_as_success(
    tmp_path: Path,
    monkeypatch,
    contents: str,
    reason: str,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    trf_dat = tmp_path / "invalid.dat"
    trf_dat.write_text(contents)
    monkeypatch.setattr(masking, "run_trf", lambda *_args, **_kwargs: trf_dat)
    result = masking.mask_genome_pipeline(
        input_fasta,
        tmp_path / "masking",
        config=MaskingConfig(
            backend=MaskingBackend.TRF,
            failure_policy=MaskingFailurePolicy.FALLBACK,
            fallback_backend=MaskingBackend.OFF,
        ),
    )
    assert result.status == "fallback"
    assert reason.lower() in result.fallback_reason.lower()


def test_header_only_trf_artifact_is_a_valid_no_repeat_success(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    trf_dat = tmp_path / "no-repeats.dat"
    trf_dat.write_text("Sequence: demo\nParameters: 2 7 7 80 10 50 500\n")
    monkeypatch.setattr(masking, "run_trf", lambda *_args, **_kwargs: trf_dat)
    monkeypatch.setattr(masking, "_probe_backend_version", lambda _backend: "v1")
    regions, version = masking._run_backend_attempt(
        backend=MaskingBackend.TRF,
        input_fasta=input_fasta,
        attempt_dir=tmp_path / "attempt",
        config=MaskingConfig(backend=MaskingBackend.TRF),
        threads=1,
    )
    assert regions == []
    assert version == "v1"


def test_non_utf8_repeatmasker_out_is_a_typed_backend_failure(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    output_dir = tmp_path / "repeatmasker"
    output_dir.mkdir()
    (output_dir / "input.fna.out").write_bytes(b"\xff\xfe")
    monkeypatch.setattr(
        masking.subprocess,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(
            returncode=0,
            stdout="",
            stderr="",
        ),
    )
    with pytest.raises(masking.MaskingBackendError, match="cannot read"):
        masking.run_repeatmasker(
            input_fasta,
            output_dir,
            species="chlorophyta",
        )


def test_combined_failure_uses_named_trf_fallback_and_is_ineligible(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    calls: list[MaskingBackend] = []

    def _attempt(*, backend, attempt_dir, **_kwargs):
        calls.append(backend)
        if backend is MaskingBackend.REPEATMASKER:
            raise masking.MaskingBackendError(backend, "rm unavailable")
        return [masking.MaskedRegion("demo", 0, 2, "trf")], "trf version"

    monkeypatch.setattr(masking, "_run_backend_attempt", _attempt)
    result = masking.mask_genome_pipeline(
        input_fasta,
        tmp_path / "masking",
        config=_species_config(
            policy=MaskingFailurePolicy.FALLBACK,
            fallback=MaskingBackend.TRF,
        ),
    )
    assert calls == [
        MaskingBackend.TRF,
        MaskingBackend.REPEATMASKER,
        MaskingBackend.TRF,
    ]
    assert result.status == "fallback"
    assert result.effective_backend is MaskingBackend.TRF
    assert result.configured_fallback_backend is MaskingBackend.TRF
    assert result.fallback_backend is MaskingBackend.TRF
    assert not result.benchmark_eligible
    assert masking.load_masking_result(
        result.status_path,
        repeat_regions=result.repeat_regions,
        expected_config=_species_config(
            policy=MaskingFailurePolicy.FALLBACK,
            fallback=MaskingBackend.TRF,
        ),
        expected_input=input_fasta,
    ) == result


def test_failed_named_fallback_is_fatal(tmp_path: Path, monkeypatch) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")

    def _fail(*, backend, **_kwargs):
        raise masking.MaskingBackendError(backend, "still failed")

    monkeypatch.setattr(masking, "_run_backend_attempt", _fail)
    with pytest.raises(masking.MaskingBackendError, match="still failed"):
        masking.mask_genome_pipeline(
            input_fasta,
            tmp_path / "masking",
            config=_species_config(
                backend=MaskingBackend.REPEATMASKER,
                policy=MaskingFailurePolicy.FALLBACK,
                fallback=MaskingBackend.TRF,
            ),
        )


def test_attempt_directory_is_cleaned_before_backend_run(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    attempt = tmp_path / "attempt"
    attempt.mkdir()
    (attempt / "stale.dat").write_text("stale")

    def _run_trf(_input, output):
        assert not (output / "stale.dat").exists()
        path = output / "fresh.dat"
        path.write_text("Sequence: demo\n")
        return path

    monkeypatch.setattr(masking, "run_trf", _run_trf)
    monkeypatch.setattr(
        masking,
        "parse_trf_output",
        lambda _path: [masking.MaskedRegion("demo", 0, 2, "trf")],
    )
    monkeypatch.setattr(masking, "_probe_backend_version", lambda _backend: "v1")
    regions, version = masking._run_backend_attempt(
        backend=MaskingBackend.TRF,
        input_fasta=input_fasta,
        attempt_dir=attempt,
        config=MaskingConfig(backend=MaskingBackend.TRF),
        threads=1,
    )
    assert len(regions) == 1
    assert version == "v1"


@pytest.mark.parametrize(
    "region",
    [
        masking.MaskedRegion("missing", 0, 1, "trf"),
        masking.MaskedRegion("demo", -1, 1, "trf"),
        masking.MaskedRegion("demo", 2, 2, "trf"),
        masking.MaskedRegion("demo", 0, 9, "trf"),
    ],
)
def test_invalid_parser_regions_are_rejected(
    tmp_path: Path,
    monkeypatch,
    region: masking.MaskedRegion,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    trf_dat = tmp_path / "x.dat"
    trf_dat.write_text("Sequence: demo\n")
    monkeypatch.setattr(masking, "run_trf", lambda *_args: trf_dat)
    monkeypatch.setattr(masking, "parse_trf_output", lambda _path: [region])
    monkeypatch.setattr(masking, "_probe_backend_version", lambda _backend: "v1")
    with pytest.raises(masking.MaskingBackendError, match="repeat|interval|unknown"):
        masking._run_backend_attempt(
            backend=MaskingBackend.TRF,
            input_fasta=input_fasta,
            attempt_dir=tmp_path / "attempt",
            config=MaskingConfig(backend=MaskingBackend.TRF),
            threads=1,
        )


def test_repeatmasker_output_must_preserve_fasta_shape(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    wrong = _fasta(tmp_path / "wrong.masked", "ACGT", seq_id="renamed")
    monkeypatch.setattr(masking, "run_repeatmasker", lambda *_args, **_kwargs: wrong)
    with pytest.raises(masking.MaskingBackendError, match="IDs/order/lengths"):
        masking._run_backend_attempt(
            backend=MaskingBackend.REPEATMASKER,
            input_fasta=input_fasta,
            attempt_dir=tmp_path / "attempt",
            config=_species_config(MaskingBackend.REPEATMASKER),
            threads=1,
        )


def test_repeatmasker_output_must_preserve_sequence_content(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna", "ACGT")
    corrupted = _fasta(tmp_path / "corrupted.masked", "ACGA")
    monkeypatch.setattr(
        masking,
        "run_repeatmasker",
        lambda *_args, **_kwargs: corrupted,
    )
    with pytest.raises(masking.MaskingBackendError, match="changed sequence content"):
        masking._run_backend_attempt(
            backend=MaskingBackend.REPEATMASKER,
            input_fasta=input_fasta,
            attempt_dir=tmp_path / "attempt",
            config=_species_config(MaskingBackend.REPEATMASKER),
            threads=1,
        )


def test_repeatmasker_parser_does_not_claim_preexisting_lowercase(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna", "aCGTacgt")
    masked = _fasta(tmp_path / "output.masked", "acgtacgt")
    monkeypatch.setattr(
        masking,
        "run_repeatmasker",
        lambda *_args, **_kwargs: masked,
    )
    monkeypatch.setattr(masking, "_probe_backend_version", lambda _backend: "v1")
    regions, _version = masking._run_backend_attempt(
        backend=MaskingBackend.REPEATMASKER,
        input_fasta=input_fasta,
        attempt_dir=tmp_path / "attempt",
        config=_species_config(MaskingBackend.REPEATMASKER),
        threads=1,
    )
    assert [(region.start, region.end) for region in regions] == [(1, 4)]


def test_status_and_output_mutations_are_rejected(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    monkeypatch.setattr(
        masking,
        "_run_backend_attempt",
        lambda **_kwargs: (
            [masking.MaskedRegion("demo", 0, 2, "trf")],
            "trf version",
        ),
    )
    result = masking.mask_genome_pipeline(
        input_fasta,
        tmp_path / "masking",
        config=MaskingConfig(backend=MaskingBackend.TRF),
    )
    original_status = result.status_path.read_text()
    result.status_path.write_text(original_status.replace('"masked_bases": 2', '"masked_bases": 3'))
    with pytest.raises(ValueError, match="fingerprint"):
        masking.load_masking_result(result.status_path)
    result.status_path.write_text(original_status)
    result.output_path.write_text(">demo\nNNNNNNNN\n")
    with pytest.raises(ValueError, match="output SHA256"):
        masking.load_masking_result(result.status_path)


def test_status_write_is_atomic_on_serialization_failure(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    output_dir = tmp_path / "masking"

    def _explode(*_args, **_kwargs):
        raise RuntimeError("serialization failure")

    monkeypatch.setattr(masking.json, "dump", _explode)
    with pytest.raises(RuntimeError, match="serialization failure"):
        masking.mask_genome_pipeline(
            input_fasta,
            output_dir,
            config=MaskingConfig(),
        )
    assert not (output_dir / "masking_status.json").exists()


def test_completion_manifest_fails_closed_and_detects_status_mutation(
    tmp_path: Path,
) -> None:
    (tmp_path / "run.log").write_text("stale success")
    with pytest.raises(ValueError, match="without a valid"):
        _write_completion_manifest(
            tmp_path,
            genome_id="demo",
            status="success",
            fingerprint="requested",
        )
    assert not (tmp_path / "run.log").exists()

    input_fasta = _fasta(tmp_path / "input.fna")
    result = masking.mask_genome_pipeline(
        input_fasta,
        tmp_path / "phase0" / "masking",
        config=MaskingConfig(),
    )
    result.status_path.write_text("{broken")
    (tmp_path / "run.log").write_text("stale success")
    with pytest.raises(ValueError, match="without a valid"):
        _write_completion_manifest(
            tmp_path,
            genome_id="demo",
            status="success",
            fingerprint="requested",
        )
    assert not (tmp_path / "run.log").exists()
    assert not (tmp_path / "virosync_run_complete.json").exists()

    result = masking.mask_genome_pipeline(
        input_fasta,
        tmp_path / "phase0" / "masking",
        config=MaskingConfig(),
    )
    manifest = _write_completion_manifest(
        tmp_path,
        genome_id="demo",
        status="success",
        fingerprint="requested",
    )
    assert _valid_completion_manifest(manifest, expected_fingerprint="requested")
    result.status_path.write_text(result.status_path.read_text() + " ")
    assert not _valid_completion_manifest(manifest, expected_fingerprint="requested")
    assert _manifest_is_stale(tmp_path, expected_fingerprint="requested")


def test_enabled_resume_rejects_current_input_content_mutation(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    monkeypatch.setattr(
        masking,
        "_run_backend_attempt",
        lambda **_kwargs: (
            [masking.MaskedRegion("demo", 0, 2, "trf")],
            "trf version",
        ),
    )
    masking.mask_genome_pipeline(
        input_fasta,
        tmp_path / "phase0" / "masking",
        config=MaskingConfig(backend=MaskingBackend.TRF),
    )
    manifest = _write_completion_manifest(
        tmp_path,
        genome_id="demo",
        status="success",
        fingerprint="requested",
    )
    assert _valid_completion_manifest(
        manifest,
        expected_fingerprint="requested",
        expected_input=input_fasta,
    )

    _fasta(input_fasta, "TGCATGCA")

    assert not _valid_completion_manifest(
        manifest,
        expected_fingerprint="requested",
        expected_input=input_fasta,
    )
    assert _manifest_is_stale(
        tmp_path,
        expected_fingerprint="requested",
        expected_input=input_fasta,
    )


def test_provenance_reuses_validated_masking_status_versions(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    config = MaskingConfig(backend=MaskingBackend.TRF)
    monkeypatch.setattr(
        masking,
        "_run_backend_attempt",
        lambda **_kwargs: (
            [masking.MaskedRegion("demo", 0, 2, "trf")],
            "Tandem Repeats Finder, Version TEST",
        ),
    )
    result = masking.mask_genome_pipeline(
        input_fasta,
        tmp_path / "run" / "phase0" / "masking",
        config=config,
    )
    monkeypatch.setattr(provenance, "capture_tool_versions", lambda: {"virosync": "test"})
    provenance.write_provenance(
        tmp_path / "run",
        {
            "masking": config,
            "masking_status_path": result.status_path,
            "masking_status_sha256": result.status_sha256,
        },
        input_genome=input_fasta,
    )
    payload = json.loads((tmp_path / "run" / "provenance.json").read_text())
    status_payload = json.loads(result.status_path.read_text())
    assert payload["masking_status"]["backend_versions"] == status_payload["backend_versions"]
    assert payload["masking_status"]["sha256"] == result.status_sha256
    assert payload["masking_status"]["result_fingerprint"] == status_payload["result_fingerprint"]
    assert payload["masking_status"]["requested_backend"] == "trf"
    assert payload["masking_status"]["effective_backend"] == "trf"

    with pytest.raises(ValueError, match="request mismatch"):
        provenance.write_provenance(
            tmp_path / "other",
            {
                "masking": MaskingConfig(),
                "masking_status_path": result.status_path,
                "masking_status_sha256": result.status_sha256,
            },
            input_genome=input_fasta,
        )


def test_deprecated_apply_mask_keyword_remains_safe(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = _fasta(tmp_path / "input.fna")
    monkeypatch.setattr(
        masking,
        "identify_repeats",
        lambda **_kwargs: [masking.MaskedRegion("demo", 0, 2, "trf")],
    )
    result = masking.mask_genome_pipeline(
        input_fasta,
        tmp_path / "masking",
        apply_mask=True,
    )
    assert result.masked_bases == 2
    assert result.output_path.is_file()
    assert result.legacy_adapter
    assert not result.benchmark_eligible
    with pytest.raises(ValueError, match="legacy adapter.*not reusable"):
        masking.load_masking_result(result.status_path)
