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
    assert workflow.count("pixi run virosync") >= 3
    assert workflow.count("--clean-run") == 1
    assert "--write-snapshot" in workflow
    assert "--compare-snapshot" in workflow
    assert "--require-resume" in workflow
    assert workflow.count("--expect-predictions 6") == 2
    assert workflow.count("--expect-accepted 1") == 2
    assert workflow.count("check_coordinate_outputs.py") == 2
    assert "results/ci_example/*/run.log" not in workflow
