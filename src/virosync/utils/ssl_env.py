"""
Helpers for SSL-related environment variables.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Mapping, MutableMapping

SSL_ENV_KEYS = ("SSL_CERT_FILE", "CURL_CA_BUNDLE", "REQUESTS_CA_BUNDLE")


def _is_existing_file(path_value: str) -> bool:
    return Path(path_value).expanduser().is_file()


def clear_stale_ssl_env_vars(env: MutableMapping[str, str] | None = None) -> dict[str, str]:
    target_env = os.environ if env is None else env
    removed: dict[str, str] = {}
    for key in SSL_ENV_KEYS:
        value = target_env.get(key)
        if value and not _is_existing_file(value):
            removed[key] = value
            target_env.pop(key, None)
    return removed


def sanitized_ssl_env(base_env: Mapping[str, str] | None = None) -> dict[str, str]:
    env = dict(os.environ if base_env is None else base_env)
    clear_stale_ssl_env_vars(env)
    return env
