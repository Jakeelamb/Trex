#!/usr/bin/env python3
"""Check that docs/CAPABILITIES.md names current CLI flags and benchmark scripts."""

from __future__ import annotations

import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "trex-cli" / "src" / "main.rs"
DOC = ROOT / "docs" / "CAPABILITIES.md"
SCRIPTS = ROOT / "scripts"


def fail(message: str) -> None:
    sys.stderr.write(f"validate_capabilities_doc: {message}\n")
    sys.exit(1)


def extract_assemble_flags() -> set[str]:
    lines = CLI.read_text().splitlines()
    flags: set[str] = set()
    pending_arg: str | None = None
    in_assemble = False
    brace_depth = 0

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("Assemble {"):
            in_assemble = True
            brace_depth = 1
            continue
        if not in_assemble:
            continue

        brace_depth += stripped.count("{") - stripped.count("}")
        if brace_depth <= 0:
            break

        if stripped.startswith("#[arg("):
            pending_arg = stripped
            explicit = re.search(r'long\s*=\s*"([^"]+)"', stripped)
            if explicit:
                flags.add(f"--{explicit.group(1)}")
            continue

        if pending_arg and not stripped.startswith("///") and ":" in stripped:
            field = stripped.split(":", 1)[0].strip()
            if field and "long" in pending_arg and not re.search(r'long\s*=', pending_arg):
                flags.add(f"--{field.replace('_', '-')}")
            pending_arg = None

    if not flags:
        fail("no Assemble CLI flags found")
    return flags


def main() -> None:
    if not DOC.is_file():
        fail("missing docs/CAPABILITIES.md")

    text = DOC.read_text()
    missing_flags = sorted(flag for flag in extract_assemble_flags() if flag not in text)
    if missing_flags:
        fail("docs/CAPABILITIES.md missing CLI flags: " + ", ".join(missing_flags))

    script_names = sorted(path.name for path in SCRIPTS.glob("*.sh"))
    missing_scripts = [name for name in script_names if name not in text]
    if missing_scripts:
        fail("docs/CAPABILITIES.md missing scripts: " + ", ".join(missing_scripts))

    required_phrases = [
        "Phase-1 default",
        "Phase-2 Illumina --diploid",
        "Future / deferred",
        "tools/benchmark_matrix.toml",
    ]
    missing_phrases = [phrase for phrase in required_phrases if phrase not in text]
    if missing_phrases:
        fail("docs/CAPABILITIES.md missing required sections: " + ", ".join(missing_phrases))

    print(
        f"validate_capabilities_doc: OK ({len(script_names)} scripts, {len(extract_assemble_flags())} flags)"
    )


if __name__ == "__main__":
    main()
