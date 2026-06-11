#!/usr/bin/env python3
"""Validate the governed benchmark matrix without external dependencies."""

from __future__ import annotations

import sys
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - CI uses Python 3.11+.
    sys.stderr.write("validate_benchmark_matrix: Python 3.11+ is required for tomllib\n")
    sys.exit(2)


ROOT = Path(__file__).resolve().parents[1]
MATRIX = ROOT / "tools" / "benchmark_matrix.toml"
ALLOWED_TIERS = {"pr", "main", "nightly", "manual"}
REQUIRED_ROW_KEYS = {
    "id",
    "technology",
    "organism",
    "license",
    "provenance",
    "depth_class",
    "ci_tier",
    "fixtures",
}
SCRIPT_KEYS = ("pr_scripts", "main_scripts", "nightly_scripts", "manual_scripts")
LIST_KEYS = (*SCRIPT_KEYS, "fixtures", "optional_tools", "artifacts")


def fail(message: str) -> None:
    sys.stderr.write(f"validate_benchmark_matrix: {message}\n")
    sys.exit(1)


def require_rel_path(row_id: str, field: str, value: str, *, must_exist: bool) -> None:
    if Path(value).is_absolute():
        fail(f"{row_id}: {field} must be repo-relative, got absolute path {value!r}")
    if ".." in Path(value).parts:
        fail(f"{row_id}: {field} must not escape repo root: {value!r}")
    if must_exist and not (ROOT / value).exists():
        fail(f"{row_id}: {field} path does not exist: {value}")


def main() -> None:
    if not MATRIX.is_file():
        fail(f"missing {MATRIX.relative_to(ROOT)}")

    data = tomllib.loads(MATRIX.read_text())
    if data.get("schema_version") != 1:
        fail("schema_version must be 1")

    rows = data.get("rows")
    if not isinstance(rows, list) or not rows:
        fail("expected at least one [[rows]] entry")

    seen_ids: set[str] = set()
    for i, row in enumerate(rows, start=1):
        if not isinstance(row, dict):
            fail(f"row {i}: expected table")

        missing = sorted(REQUIRED_ROW_KEYS - row.keys())
        row_id = str(row.get("id", f"row {i}"))
        if missing:
            fail(f"{row_id}: missing required keys: {', '.join(missing)}")

        if row_id in seen_ids:
            fail(f"duplicate row id {row_id!r}")
        seen_ids.add(row_id)

        if row["ci_tier"] not in ALLOWED_TIERS:
            fail(f"{row_id}: ci_tier must be one of {sorted(ALLOWED_TIERS)}")

        for key in REQUIRED_ROW_KEYS - {"fixtures"}:
            if not isinstance(row[key], str) or not row[key].strip():
                fail(f"{row_id}: {key} must be a non-empty string")

        for key in LIST_KEYS:
            if key in row and (
                not isinstance(row[key], list)
                or not all(isinstance(item, str) and item.strip() for item in row[key])
            ):
                fail(f"{row_id}: {key} must be a list of non-empty strings")

        script_paths = [script for key in SCRIPT_KEYS for script in row.get(key, [])]
        if not script_paths:
            fail(f"{row_id}: at least one *_scripts list is required")
        if row["ci_tier"] == "pr" and not row.get("pr_scripts"):
            fail(f"{row_id}: ci_tier=pr requires pr_scripts")

        for script in script_paths:
            require_rel_path(row_id, "script", script, must_exist=True)
            if not script.startswith("scripts/"):
                fail(f"{row_id}: script path must live under scripts/: {script}")

        for fixture in row["fixtures"]:
            require_rel_path(row_id, "fixture", fixture, must_exist=True)

        for artifact in row.get("artifacts", []):
            require_rel_path(row_id, "artifact", artifact, must_exist=False)

        digest_manifest = row.get("digest_manifest")
        if digest_manifest is not None:
            if not isinstance(digest_manifest, str) or not digest_manifest.strip():
                fail(f"{row_id}: digest_manifest must be a non-empty string")
            require_rel_path(row_id, "digest_manifest", digest_manifest, must_exist=True)
            if not isinstance(row.get("manifest_table"), str) or not row["manifest_table"].strip():
                fail(f"{row_id}: digest_manifest requires manifest_table")

    print(f"validate_benchmark_matrix: OK ({len(rows)} rows)")


if __name__ == "__main__":
    main()
