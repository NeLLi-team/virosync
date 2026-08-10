"""Graphviz settings and runtime probe for the required ANI report."""

from __future__ import annotations

import subprocess


ANI_NETWORK_ENGINE = "sfdp"
ANI_NETWORK_FORMAT = "png"
ANI_NETWORK_RENDER_ATTR = {
    "overlap": "scale",
    "splines": "line",
    "size": "40,40",
    "dpi": "96",
}
_PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def graphviz_runtime_error(*, timeout_seconds: float = 30) -> str | None:
    """Return a diagnostic when Graphviz cannot render the required format."""

    attributes = ", ".join(
        f'{name}="{value}"' for name, value in ANI_NETWORK_RENDER_ATTR.items()
    )
    source = f"graph G {{ graph [{attributes}]; a -- b; }}\n".encode()
    # python-graphviz 0.21 invokes this same dot -K/-T command for Graph.pipe().
    command = ["dot", f"-K{ANI_NETWORK_ENGINE}", f"-T{ANI_NETWORK_FORMAT}"]
    try:
        result = subprocess.run(
            command,
            input=source,
            check=False,
            capture_output=True,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired:
        return f"Graphviz sfdp PNG render timed out after {timeout_seconds:g} seconds"
    except OSError as exc:
        return f"Graphviz runtime is unavailable: {exc}"
    if result.returncode != 0 or not result.stdout.startswith(_PNG_SIGNATURE):
        detail = result.stderr.decode(errors="replace").strip()
        return "Graphviz cannot render the required sfdp PNG report" + (
            f": {detail}" if detail else ""
        )
    return None
