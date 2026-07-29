"""Filesystem-path guards for identifiers derived from biological inputs."""

from __future__ import annotations

import hashlib
import re
import unicodedata
from collections.abc import Iterable
from pathlib import Path
from urllib.parse import quote


_PORTABLE_FILENAME_COMPONENT = re.compile(r"[A-Za-z0-9._-]+\Z")
_MAX_FILENAME_COMPONENT_LENGTH = 180


def validate_path_component(value: str, label: str = "path component") -> str:
    """Return *value* after rejecting components that can alter path structure."""
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a string")
    if not value:
        raise ValueError(f"{label} must not be empty")
    if value in {".", ".."}:
        raise ValueError(f"{label} must not be {value!r}")
    if "/" in value or "\\" in value:
        raise ValueError(f"{label} must be one path component: {value!r}")
    if any(unicodedata.category(character) == "Cc" for character in value):
        raise ValueError(f"{label} must not contain control characters: {value!r}")
    return value


def require_strict_child(root: Path, candidate: Path) -> Path:
    """Resolve and return *candidate* only when it is strictly below *root*."""
    try:
        resolved_root = Path(root).resolve(strict=False)
        resolved_candidate = Path(candidate).resolve(strict=False)
    except (OSError, RuntimeError) as exc:
        raise ValueError(f"could not resolve path safely: {candidate}") from exc

    if resolved_candidate == resolved_root or resolved_root not in resolved_candidate.parents:
        raise ValueError(
            f"path must be a strict child of {resolved_root}: {resolved_candidate}"
        )
    return resolved_candidate


def safe_filename_component(value: str) -> str:
    """Return a portable, collision-resistant filename component for a raw ID.

    Ordinary ASCII components are preserved. Values that require encoding retain a
    readable percent-encoded prefix and gain a digest, while callers keep the raw ID
    in biological records and manifests.
    """
    if not isinstance(value, str):
        raise TypeError("filename component must be a string")

    encoded_bytes = value.encode("utf-8", errors="surrogatepass")
    is_portable = (
        bool(value)
        and value not in {".", ".."}
        and _PORTABLE_FILENAME_COMPONENT.fullmatch(value) is not None
        and len(encoded_bytes) <= _MAX_FILENAME_COMPONENT_LENGTH
    )
    if is_portable:
        return value

    encoded = quote(value, safe="._-", errors="surrogatepass").replace("~", "%7E")
    digest = hashlib.sha256(encoded_bytes).hexdigest()[:16]
    prefix_limit = _MAX_FILENAME_COMPONENT_LENGTH - len(digest) - 3
    encoded = encoded[:prefix_limit] or "empty"

    # The leading percent marker keeps transformed values in a namespace that can
    # never be returned unchanged for a raw portable identifier.
    return f"%{encoded}--{digest}"


def safe_filename_components(
    values: Iterable[str],
    label: str = "identifier",
) -> dict[str, str]:
    """Return distinct safe components, rejecting duplicate or colliding values."""
    raw_values = list(values)
    indices_by_value: dict[str, list[int]] = {}
    for index, value in enumerate(raw_values):
        indices_by_value.setdefault(value, []).append(index)

    duplicates = {
        value: indices
        for value, indices in indices_by_value.items()
        if len(indices) > 1
    }
    if duplicates:
        details = "; ".join(
            f"{value!r} at indices {', '.join(str(index) for index in indices)}"
            for value, indices in sorted(duplicates.items())
        )
        raise ValueError(f"duplicate {label}s: {details}")

    components: dict[str, str] = {}
    raw_value_by_component: dict[str, str] = {}
    for value in raw_values:
        component = safe_filename_component(value)
        conflicting_value = raw_value_by_component.get(component)
        if conflicting_value is not None:
            raise ValueError(
                f"{label} filename encoding collision: "
                f"{conflicting_value!r} and {value!r}"
            )
        components[value] = component
        raw_value_by_component[component] = value
    return components
