"""Tests for the reference-mask settings: config, CLI, run fingerprint and provenance."""

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

import virosync.orchestration.cli as orchestration_cli
from virosync.config import ApplicationConfig, ConfigError, FeatureResolution, PipelineConfig
from virosync.orchestration._flows.single_genome.manifest import _compute_config_fingerprint
from virosync.orchestration.cli import _build_pipeline_config, orchestrate
from virosync.pipeline.reference_mask import ReferenceMask
from virosync.utils import provenance


@pytest.fixture(autouse=True)
def _skip_native_tool_preflight(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(orchestration_cli, "_preflight_native_tools", lambda: None)


def test_reference_mask_defaults_mask_the_own_assembly_only() -> None:
    config = PipelineConfig()

    assert config.host.reference_mask_taxa == ()
    assert config.host.mask_own_assembly is True
    assert config.host.self_reference_keys == ()


def test_reference_mask_settings_round_trip_through_the_application_config() -> None:
    application = ApplicationConfig.from_dict(
        {"host": {"reference_mask_taxa": ["Phytophthora", "Collodictyon triciliatum"], "mask_own_assembly": False}}
    )

    assert application.pipeline.host.reference_mask_taxa == ("Phytophthora", "Collodictyon triciliatum")
    assert application.pipeline.host.mask_own_assembly is False
    assert ApplicationConfig.from_dict(application.to_dict()) == application


@pytest.mark.parametrize(
    ("host", "message"),
    [
        ({"reference_mask_taxa": "Phytophthora"}, "must be a list"),
        ({"reference_mask_taxa": [3]}, "must be a string"),
        ({"mask_own_assembly": "no"}, "must be a boolean"),
    ],
)
def test_reference_mask_settings_are_type_checked(host: dict, message: str) -> None:
    with pytest.raises(ConfigError, match=message):
        PipelineConfig.from_dict({"host": host})


def test_reference_mask_settings_reach_the_flow_kwargs_and_overrides() -> None:
    config = PipelineConfig().with_overrides(reference_mask_taxa=("Phytophthora",), mask_own_assembly=False)

    flow_kwargs = config.to_flow_kwargs()

    assert flow_kwargs["reference_mask_taxa"] == ("Phytophthora",)
    assert flow_kwargs["mask_own_assembly"] is False


@pytest.mark.parametrize(
    "overrides",
    [{"reference_mask_taxa": ("Phytophthora",)}, {"mask_own_assembly": False}, {"self_reference_keys": ("JGI-1",)}],
    ids=["taxa", "own-assembly", "self-keys"],
)
def test_changing_the_mask_changes_the_run_fingerprint(overrides: dict) -> None:
    baseline = PipelineConfig()

    changed = baseline.with_overrides(**overrides)

    assert _compute_config_fingerprint(baseline.to_flow_kwargs()) != _compute_config_fingerprint(
        changed.to_flow_kwargs()
    )


def test_taxa_need_the_taxonomy_labels_file(tmp_path: Path) -> None:
    labels = tmp_path / "labels.tsv"
    labels.write_text("EUK__G1\tEUK|a|b|c|d|e|f\n", encoding="utf-8")
    config = PipelineConfig().with_overrides(reference_mask_taxa=("Phytophthora",))
    problem = "host.reference_mask_taxa requires databases.taxonomy_labels_file"

    assert problem in config.validate()
    assert problem not in config.with_overrides(taxonomy_labels_file=labels).validate()
    assert problem not in PipelineConfig().validate()


def test_cli_values_replace_the_yaml_mask_and_absent_values_keep_it() -> None:
    application = ApplicationConfig.from_dict({"host": {"reference_mask_taxa": ["From yaml"]}})

    kept = _build_pipeline_config(yaml_config=application, clean_run=False)
    replaced = _build_pipeline_config(
        yaml_config=application,
        clean_run=False,
        mask_taxa=("Phytophthora", "Collodictyon"),
        mask_own_assembly=False,
    )

    assert kept.host.reference_mask_taxa == ("From yaml",)
    assert kept.host.mask_own_assembly is True
    assert replaced.host.reference_mask_taxa == ("Phytophthora", "Collodictyon")
    assert replaced.host.mask_own_assembly is False


def _run_with_stubbed_runner(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *options: str) -> PipelineConfig:
    received: dict[str, object] = {}
    input_fasta = tmp_path / "input.fna"
    input_fasta.write_text(">scaffold_1\nACGT\n", encoding="utf-8")
    monkeypatch.setattr(orchestration_cli, "_resolve_pipeline_resources", lambda config, orchestration, path: config)
    monkeypatch.setattr(
        orchestration_cli,
        "_resolve_optional_features",
        lambda config: (config, {name: FeatureResolution(False, False, False) for name in ("boltz", "tmvec")}),
    )
    monkeypatch.setattr(orchestration_cli, "_validate_runtime_config", lambda config: None)

    def fake_runner(**kwargs: object) -> list[dict]:
        received.update(kwargs)
        return [
            {
                "success": True,
                "benchmark_eligible": True,
                "legacy_resume": False,
                "accepted": 0,
                "predictions": 0,
                "elapsed_sec": 0.0,
                "genome_id": "input",
            }
        ]

    monkeypatch.setattr(orchestration_cli, "run_batch_python", fake_runner)

    result = CliRunner().invoke(
        orchestrate,
        ["run", "-i", str(input_fasta), "-o", str(tmp_path / "results"), *options],
    )

    assert result.exit_code == 0, result.output
    return received["config"]


def test_run_command_collects_repeated_mask_taxa_and_the_own_assembly_switch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = _run_with_stubbed_runner(
        tmp_path,
        monkeypatch,
        "--mask-taxon",
        "Phytophthora ramorum",
        "--mask-taxon",
        "Collodictyon",
        "--no-mask-own-assembly",
    )

    assert config.host.reference_mask_taxa == ("Phytophthora ramorum", "Collodictyon")
    assert config.host.mask_own_assembly is False


def test_run_command_collects_repeated_self_reference_keys(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    config = _run_with_stubbed_runner(
        tmp_path, monkeypatch, "--self-reference-key", "JGI-12345", "--self-reference-key", "EUK__IMG-678"
    )

    assert config.host.self_reference_keys == ("JGI-12345", "EUK__IMG-678")


def test_run_command_masks_only_the_own_assembly_by_default(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    config = _run_with_stubbed_runner(tmp_path, monkeypatch)

    assert config.host.reference_mask_taxa == ()
    assert config.host.mask_own_assembly is True


def test_provenance_records_the_reference_mask(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(provenance, "capture_tool_versions", lambda: {"virosync": "test"})
    monkeypatch.setattr(provenance, "capture_gene_caller_identity", lambda version: {"tool": "prodigal-gv"})
    mask = ReferenceMask(
        genome_keys=frozenset({"GCA-052324655-1", "GCA_052324655.1", "EP00001"}),
        taxa=("Collodictyon",),
        own_assembly_keys=("GCA-052324655-1", "GCA_052324655.1"),
    )

    provenance.write_provenance(tmp_path / "run", {"reference_mask": mask.to_record()})

    document = json.loads((tmp_path / "run" / "provenance.json").read_text(encoding="utf-8"))
    assert document["configuration"]["reference_mask"] == {
        "n_genome_keys": 3,
        "taxa": ["Collodictyon"],
        "own_assembly_keys": ["GCA-052324655-1", "GCA_052324655.1"],
    }


@pytest.mark.parametrize(
    "host",
    [
        {"self_reference_keys": [""]},
        {"self_reference_keys": ["EUK__"]},
        {"self_reference_keys": [" "]},
        {"self_reference_keys": ["foo bar"]},
        {"self_reference_keys": ["foo|protein"]},
        {"self_reference_keys": ["JGI-1"], "mask_own_assembly": False},
        {"reference_mask_taxa": ["NA"]},
        {"reference_mask_taxa": ["EUK"]},
        {"reference_mask_taxa": ["unclassified"]},
        {"reference_mask_taxa": [" "]},
    ],
)
def test_invalid_reference_mask_values_fail_at_config_boundary(host: dict) -> None:
    with pytest.raises(ConfigError):
        PipelineConfig.from_dict({"host": host})


@pytest.mark.parametrize("source", ["cli", "yaml"])
def test_cli_rejects_global_self_keys_for_multiple_genomes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, source: str
) -> None:
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "one.fna").write_text(">one\nACGT\n")
    (genomes / "two.fna").write_text(">two\nACGT\n")
    options = ["--self-reference-key", "GCA-000002455-1"]
    if source == "yaml":
        config = tmp_path / "mask.yaml"
        config.write_text("host:\n  self_reference_keys: [GCA-000002455-1]\n")
        options = ["--config", str(config)]
    result = CliRunner().invoke(orchestrate, ["run", "-i", str(genomes), "-o", str(tmp_path / "out"), *options])
    assert result.exit_code != 0
    assert "self_reference_keys requires exactly one input genome" in result.output
    assert not (tmp_path / "out").exists()


def test_python_batch_rejects_global_self_keys_before_output(tmp_path: Path) -> None:
    from virosync.orchestration.python_runner import run_batch_python

    with pytest.raises(ConfigError, match="self_reference_keys requires exactly one input genome"):
        run_batch_python(
            genome_paths=[tmp_path / "a.fna", tmp_path / "b.fna"],
            output_base_dir=tmp_path / "out",
            config=PipelineConfig().with_overrides(self_reference_keys=("GCA-000002455-1",)),
            max_concurrent_genomes=1,
        )
    assert not (tmp_path / "out").exists()
