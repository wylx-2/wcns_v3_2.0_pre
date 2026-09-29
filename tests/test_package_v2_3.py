#!/usr/bin/env python3
"""Exercise v2.3 package path and integrity rejection rules."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import sys
import tempfile
from pathlib import Path


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def expect_failure(action, label: str) -> None:
    try:
        action()
    except RuntimeError:
        return
    raise AssertionError(f"expected failure was not raised: {label}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packager", type=Path, required=True)
    args = parser.parse_args()
    specification = importlib.util.spec_from_file_location("package_v2_3", args.packager)
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot import v2.3 packager")
    package = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = package
    specification.loader.exec_module(package)

    assert package.forbidden(Path("cases/input.wcns"))
    assert package.forbidden(Path("tests/test.cpp"))
    assert package.forbidden(Path("result/field.cgns"))
    assert not package.forbidden(Path("src/app/wcns_run.cpp"))
    assert not package.safe_relative(Path("../escape"))

    with tempfile.TemporaryDirectory(prefix="wcns-package-test-") as temporary:
        root = Path(temporary) / package.PACKAGE_NAME
        root.mkdir()
        (root / "CMakeLists.txt").write_text(
            'project(wcns VERSION 2.3.0 LANGUAGES CXX)\n'
            'set(WCNS_PROGRAM_VERSION "2.3")\n',
            encoding="utf-8",
        )
        (root / "README.md").write_text("WCNS v2.3\n", encoding="utf-8")
        (root / "WCNS_SOURCE_REVISION").write_text("a" * 40 + "\n", encoding="ascii")
        content = (root / "CMakeLists.txt", root / "README.md", root / "WCNS_SOURCE_REVISION")
        lines = [f"{digest(path)}  {path.name}" for path in content]
        (root / "PACKAGE_CONTENTS.sha256").write_text(
            "\n".join(lines) + "\n", encoding="utf-8"
        )
        result = package.verify_directory(root)
        assert result["revision"] == "a" * 40
        assert result["files"] == 4

        (root / "README.md").write_text("tampered\n", encoding="utf-8")
        expect_failure(lambda: package.verify_directory(root), "tampered payload")

        (root / "README.md").write_text("WCNS v2.3\n", encoding="utf-8")
        (root / "guide.md").write_text("[missing](missing.md)\n", encoding="utf-8")
        content = content + (root / "guide.md",)
        lines = [f"{digest(path)}  {path.name}" for path in content]
        (root / "PACKAGE_CONTENTS.sha256").write_text(
            "\n".join(lines) + "\n", encoding="utf-8"
        )
        expect_failure(lambda: package.verify_directory(root), "broken markdown link")

        unsafe = root / "unsafe.sha256"
        unsafe.write_text(f"{'0' * 64}  ../escape\n", encoding="utf-8")
        expect_failure(lambda: package.parse_manifest(unsafe), "manifest path traversal")

    print("v2.3 package safety rules verified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
