from __future__ import annotations

from pathlib import Path
import re
import tomllib


ROOT = Path(__file__).resolve().parents[1]
WORKFLOWS = (
    ROOT / ".github/workflows/tests.yml",
    ROOT / ".github/workflows/production-guards.yml",
    ROOT / ".github/workflows/example-smoke.yml",
)


def _text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_ruff_version_and_local_ci_command_are_locked() -> None:
    with (ROOT / "pixi.toml").open("rb") as handle:
        manifest = tomllib.load(handle)
    assert manifest["feature"]["dev"]["dependencies"]["ruff"] == "==0.15.21"
    assert manifest["tasks"]["lint"] == "ruff check src tests scripts"
    lock = _text(ROOT / "pixi.lock")
    assert lock.count("/ruff-0.15.21-") >= 2
    tests_workflow = _text(ROOT / ".github/workflows/tests.yml")
    assert "run: pixi run lint" in tests_workflow
    assert "ruff-action" not in tests_workflow


def test_ci_pixi_binary_url_and_checksum_are_exact() -> None:
    installer = _text(ROOT / "scripts/ci/install_pixi.sh")
    assert (
        "https://github.com/prefix-dev/pixi/releases/download/v0.72.0/"
        "pixi-x86_64-unknown-linux-musl"
    ) in installer
    assert (
        "6304fe3178f3036e2c95151bbb318592fae5c31a77f5a6f4319bb023a479d4b9"
        in installer
    )
    assert "sha256sum" in installer


def test_workflow_permissions_actions_and_installs_are_pinned() -> None:
    checkout_count = 0
    persist_count = 0
    for path in WORKFLOWS:
        workflow = _text(path)
        assert re.search(r"(?m)^permissions:\n\s+contents: read$", workflow)
        references = re.findall(
            r"(?m)^\s*-?\s*uses:\s*[^@\s]+@([^\s#]+)", workflow
        )
        assert references
        assert all(re.fullmatch(r"[0-9a-f]{40}", ref) for ref in references)
        checkout_count += workflow.count("uses: actions/checkout@")
        persist_count += workflow.count("persist-credentials: false")
        assert "pixi.sh/install.sh" not in workflow
        assert "scripts/ci/install_pixi.sh" in workflow
        for command in re.findall(
            r"(?m)^\s*run:\s*(pixi install[^\n]*)", workflow
        ):
            assert "--locked" in command
        if "runs-on: ubuntu" in workflow:
            assert "runs-on: ubuntu-24.04" in workflow
            assert "ubuntu-latest" not in workflow
    assert checkout_count == persist_count


def test_release_smoke_runs_full_clean_and_unchanged_resume_checks() -> None:
    workflow = _text(ROOT / ".github/workflows/example-smoke.yml")
    assert "schedule:" in workflow and "workflow_dispatch:" in workflow
    assert "tags:" in workflow and '"v*"' in workflow
    assert "resources verify" in workflow and "--full" in workflow
    assert workflow.count("pixi run virosync") >= 4
    assert workflow.count("--clean-run") == 2
    assert "--write-snapshot" in workflow
    assert "--compare-snapshot" in workflow
    assert "--require-resume" in workflow
    assert workflow.count("--expect-predictions 6") == 2
    assert workflow.count("--expect-predictions 5") == 1
    assert workflow.count("--expect-accepted 1") == 2
    assert workflow.count("--expect-accepted 2") == 1
    assert workflow.count("--expect-canonical-eve-id") == 2
    assert workflow.count("--expect-detailed-eve-id") == 5
    assert "EVE_DS113495.1_18305-37386" not in workflow
    assert "EVE_DS113200.1_129184-151689" in workflow
    assert "EVE_DS113495.1_58468-89417" in workflow
    assert workflow.count("check_coordinate_outputs.py") == 3
    assert "results/ci_example/*/run.log" not in workflow
    assert "example/frameshift/trichomonas-g3.fna" in workflow
    assert "--frameshift-screening" in workflow
    assert workflow.count("--threads-per-worker 8") == 1
    assert "awk 'END { exit NR < 2 }'" in workflow
    assert "'^>.*_VSR'" in workflow
    assert "7842ebd58b96591b4b60863ee5c33e49eb79eccc" in workflow
    assert "0f4e71832d6ba1e4c65039ba4b4663c546a041fa" in workflow
    assert '>> "$GITHUB_PATH"' in workflow


def test_production_guard_checks_public_resource_size() -> None:
    workflow = _text(ROOT / ".github/workflows/production-guards.yml")
    assert "EXPECTED_BYTES=5877324818" in workflow
    assert 'test "$actual_bytes" = "$EXPECTED_BYTES"' in workflow


def test_frameshift_example_task_is_explicit_opt_in() -> None:
    pixi = _text(ROOT / "pixi.toml")
    task = next(
        line
        for line in pixi.splitlines()
        if line.startswith("example-frameshift = ")
    )
    assert "example/frameshift/trichomonas-g3.fna" in task
    assert "results/example-frameshift/" in task
    assert "--frameshift-screening" in task
    assert "--clean-run" in task
    assert 'depends-on = ["setup-databases"]' in task
