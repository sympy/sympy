#!/usr/bin/env python
"""Keep the legacy requirements files in sync with pyproject.toml."""

from __future__ import annotations

import argparse
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError as exc:
    raise SystemExit("Regenerating requirements requires Python 3.11+.") from exc


ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = {
    "dev-tools": ROOT / "requirements-dev.txt",
    "docs": ROOT / "doc" / "requirements.txt",
}
HEADER = (
    "# Generated from pyproject.toml by bin/generate_requirements.py.\n"
    "# Do not edit manually.\n"
)


def expand_group(groups, name, stack=()):
    """Expand PEP 735 include-group references, retaining declaration order."""
    if name in stack:
        raise ValueError("Cyclic dependency groups: " + " -> ".join((*stack, name)))
    if name not in groups:
        raise ValueError(f"Unknown dependency group: {name}")
    requirements = []
    for item in groups[name]:
        if isinstance(item, str):
            requirements.append(item)
        elif isinstance(item, dict) and set(item) == {"include-group"}:
            requirements.extend(
                expand_group(groups, item["include-group"], (*stack, name))
            )
        else:
            raise ValueError(f"Unsupported entry in dependency group {name}: {item!r}")
    return requirements


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check", action="store_true",
        help="Fail if committed requirements files differ; do not write files.",
    )
    args = parser.parse_args()

    with (ROOT / "pyproject.toml").open("rb") as f:
        groups = tomllib.load(f)["dependency-groups"]

    stale = []
    for group, path in OUTPUTS.items():
        expected = HEADER + "".join(f"{r}\n" for r in expand_group(groups, group))
        if not path.is_file() or path.read_text(encoding="utf-8") != expected:
            stale.append(path.relative_to(ROOT))
            if not args.check:
                path.write_text(expected, encoding="utf-8", newline="\n")

    if stale:
        if args.check:
            print("Out-of-date generated requirements:")
            for path in stale:
                print(f"  {path}")
            print("Run: python bin/generate_requirements.py")
            return 1
        for path in stale:
            print(f"Updated {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
