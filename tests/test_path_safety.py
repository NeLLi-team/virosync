from __future__ import annotations

from pathlib import Path

import pytest

from virosync.orchestration._flows.single_genome import orchestrator
from virosync.utils.path_safety import (
    require_strict_child,
    safe_filename_component,
    safe_filename_components,
    validate_path_component,
)


@pytest.mark.parametrize("value", ["", ".", "..", "a/b", "a\\b", "line\nbreak"])
def test_validate_path_component_rejects_structural_values(value: str) -> None:
    with pytest.raises((TypeError, ValueError)):
        validate_path_component(value, "test ID")


def test_require_strict_child_rejects_root_parent_and_outside(tmp_path: Path) -> None:
    root = tmp_path / "root"

    assert require_strict_child(root, root / "sample") == (root / "sample").resolve()
    for candidate in [root, root.parent, tmp_path / "elsewhere"]:
        with pytest.raises(ValueError, match="strict child"):
            require_strict_child(root, candidate)


def test_require_strict_child_rejects_symlink_escape(tmp_path: Path) -> None:
    root = tmp_path / "root"
    outside = tmp_path / "outside"
    root.mkdir()
    outside.mkdir()
    (root / "sample").symlink_to(outside, target_is_directory=True)

    with pytest.raises(ValueError, match="strict child"):
        require_strict_child(root, root / "sample")


def test_safe_filename_components_are_contained_distinct_and_portable(
    tmp_path: Path,
) -> None:
    raw_ids = [
        "EVE_ordinary-1.2",
        "EVE_NODE/1",
        "EVE_..",
        "EVE with space",
        "EVE_λ",
        "EVE_control\n",
        ".",
        "",
    ]
    components = [safe_filename_component(raw_id) for raw_id in raw_ids]

    assert components[0] == raw_ids[0]
    assert len(set(components)) == len(raw_ids)
    for component in components:
        assert component not in {"", ".", ".."}
        assert component.isascii()
        assert "/" not in component
        assert "\\" not in component
        assert not any(character.isspace() for character in component)
        candidate = require_strict_child(tmp_path, tmp_path / component)
        assert candidate.parent == tmp_path.resolve()


def test_safe_filename_component_map_rejects_duplicates_and_separates_aliases() -> None:
    raw_ids = ["a/b", "a|b", "../a", "a\\b", "a\nb"]

    components = safe_filename_components(raw_ids, label="EVE ID")

    assert list(components) == raw_ids
    assert len(set(components.values())) == len(raw_ids)
    with pytest.raises(ValueError, match="duplicate EVE IDs.*indices 0, 1"):
        safe_filename_components(["duplicate", "duplicate"], label="EVE ID")


def test_guarded_remove_never_calls_rmtree_for_unsafe_targets(
    tmp_path: Path,
    monkeypatch,
) -> None:
    output_root = tmp_path / "results"
    output_root.mkdir()
    sentinel = tmp_path / "sentinel.txt"
    sentinel.write_bytes(b"parent sentinel\n")
    outside = tmp_path / "outside"
    outside.mkdir()
    symlink = output_root / "sample"
    symlink.symlink_to(outside, target_is_directory=True)
    calls: list[Path] = []
    monkeypatch.setattr(orchestrator.shutil, "rmtree", lambda path: calls.append(Path(path)))

    unsafe_targets = [
        (output_root, "sample"),
        (output_root.parent, "sample"),
        (output_root / ".." / "sample", "sample"),
        (symlink, "sample"),
        (output_root / "..", ".."),
    ]
    for target, genome_id in unsafe_targets:
        with pytest.raises(ValueError):
            orchestrator._remove_output_dir(target, genome_id)

    assert calls == []
    assert sentinel.read_bytes() == b"parent sentinel\n"


@pytest.mark.parametrize("resume", [False, True])
def test_direct_single_genome_flow_rejects_unsafe_target_in_all_modes(
    tmp_path: Path,
    monkeypatch,
    resume: bool,
) -> None:
    genome_path = tmp_path / "...fna"
    genome_path.write_text(">scaffold\nACGT\n")
    sentinel = tmp_path / "sentinel.txt"
    sentinel.write_bytes(b"parent sentinel\n")
    calls: list[Path] = []
    monkeypatch.setattr(orchestrator.shutil, "rmtree", lambda path: calls.append(Path(path)))

    with pytest.raises(ValueError, match="genome ID"):
        orchestrator.single_genome_flow(
            genome_path=genome_path,
            output_dir=tmp_path / "results" / "..",
            genome_id="..",
            device="cpu",
            resume=resume,
        )

    assert calls == []
    assert sentinel.read_bytes() == b"parent sentinel\n"


def test_guarded_remove_passes_only_resolved_matching_child_to_rmtree(
    tmp_path: Path,
    monkeypatch,
) -> None:
    output_dir = tmp_path / "results" / "sample"
    output_dir.mkdir(parents=True)
    calls: list[Path] = []
    monkeypatch.setattr(orchestrator.shutil, "rmtree", lambda path: calls.append(Path(path)))

    orchestrator._remove_output_dir(output_dir, "sample")

    assert calls == [output_dir.resolve()]


def test_guarded_remove_real_filesystem_preserves_parent_and_sibling_sentinels(
    tmp_path: Path,
) -> None:
    output_root = tmp_path / "results"
    output_dir = output_root / "sample"
    sibling_dir = output_root / "other"
    output_dir.mkdir(parents=True)
    sibling_dir.mkdir()
    (output_dir / "stale.txt").write_bytes(b"stale output\n")
    sibling_sentinel = sibling_dir / "sentinel.txt"
    sibling_sentinel.write_bytes(b"sibling sentinel\n")
    parent_sentinel = tmp_path / "sentinel.txt"
    parent_sentinel.write_bytes(b"parent sentinel\n")

    orchestrator._remove_output_dir(output_dir, "sample")

    assert output_dir.exists() is False
    assert sibling_sentinel.read_bytes() == b"sibling sentinel\n"
    assert parent_sentinel.read_bytes() == b"parent sentinel\n"


def test_target_validation_preserves_relative_output_path(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.chdir(tmp_path)

    output_dir = orchestrator._validate_clean_run_target(
        Path("results/sample"),
        "sample",
    )

    assert output_dir == Path("results/sample")
    assert output_dir.is_absolute() is False
