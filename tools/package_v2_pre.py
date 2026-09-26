#!/usr/bin/env python3
"""Create or verify the case-free WCNS v2.0_pre source directory."""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath


REPOSITORY = Path(__file__).resolve().parents[1]
PACKAGE_NAME = "WCNS_v2.0_pre"
PROJECT_VERSION = "2.0.0"
PROGRAM_VERSION = "2.0_pre"

DIRECTORY_PREFIXES = (
    Path("include"),
    Path("src"),
    Path("third_party/cgns"),
)

TOOL_SOURCES = (
    Path("tools/compare_metric_profiles.cpp"),
    Path("tools/convert_tmr_flatplate_p3d_to_cgns.cpp"),
    Path("tools/convert_tmr_p2d_to_cgns.cpp"),
    Path("tools/extract_rans_profile.cpp"),
    Path("tools/generate_couette_validation_cgns.cpp"),
    Path("tools/generate_release_cgns.cpp"),
    Path("tools/inspect_structured_mesh.cpp"),
    Path("tools/validate_release_case.cpp"),
)

EXPLICIT_MAPPINGS = {
    Path("CMakeLists.txt"): Path("CMakeLists.txt"),
    Path("LICENSE.md"): Path("LICENSE.md"),
    Path("THIRD_PARTY_NOTICES.md"): Path("THIRD_PARTY_NOTICES.md"),
    Path("算法补充.md"): Path("算法补充.md"),
    Path("docs/v2.0-pre-package-readme.md"): Path("README.md"),
    Path("docs/user-manual.md"): Path("docs/user-manual.md"),
    Path("docs/runtime-guide.md"): Path("docs/runtime-guide.md"),
    Path("docs/developer-guide.md"): Path("docs/developer-guide.md"),
    Path("docs/config-reference-2.md"): Path("docs/config-reference-2.md"),
    Path("docs/known-limitations.md"): Path("docs/known-limitations.md"),
    Path("docs/release-validation.md"): Path("docs/release-validation.md"),
    Path("docs/release-notes-2.0-pre.md"): Path("docs/release-notes-2.0-pre.md"),
    Path("docs/linux-server-guide.md"): Path("docs/linux-server-guide.md"),
    Path("docs/v2.0-pre-code-review.md"): Path("docs/v2.0-pre-code-review.md"),
    Path("docs/v2.0-pre-capability-matrix.md"): Path(
        "docs/v2.0-pre-capability-matrix.md"
    ),
    Path("docs/v2.0-pre-validation.md"): Path("docs/v2.0-pre-validation.md"),
}

FORBIDDEN_PARTS = {
    ".git",
    "cases",
    "examples",
    "tests",
    "output",
    "tmp",
    "__pycache__",
}
FORBIDDEN_SUFFIXES = {".cgns", ".dat", ".plt", ".vtk", ".vtu"}


@dataclass(frozen=True)
class PayloadEntry:
    source: Path
    destination: Path


def run_git(*arguments: str) -> str:
    result = subprocess.run(
        ["git", "-c", "core.quotepath=false", *arguments],
        cwd=REPOSITORY,
        text=True,
        encoding="utf-8",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"git {' '.join(arguments)} failed: {result.stderr.strip()}")
    return result.stdout.strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_relative(path: Path) -> bool:
    pure = PurePosixPath(path.as_posix())
    return not pure.is_absolute() and ".." not in pure.parts and bool(pure.parts)


def forbidden(path: Path) -> bool:
    lowered = {part.lower() for part in path.parts}
    return bool(lowered & FORBIDDEN_PARTS) or any(
        part.lower().startswith("build-") for part in path.parts
    ) or path.suffix.lower() in FORBIDDEN_SUFFIXES


def tracked_payload() -> list[PayloadEntry]:
    tracked = {
        Path(line)
        for line in run_git("ls-tree", "-r", "--name-only", "HEAD").splitlines()
        if line
    }
    mappings = dict(EXPLICIT_MAPPINGS)
    for path in tracked:
        if any(path == prefix or prefix in path.parents for prefix in DIRECTORY_PREFIXES):
            mappings[path] = path
    for path in TOOL_SOURCES:
        mappings[path] = path

    missing = sorted(path for path in mappings if path not in tracked)
    if missing:
        raise RuntimeError(
            "v2.0_pre payload contains files not committed in HEAD: "
            + ", ".join(path.as_posix() for path in missing)
        )
    absent = sorted(path for path in mappings if not (REPOSITORY / path).is_file())
    if absent:
        raise RuntimeError(
            "v2.0_pre payload is missing from the worktree: "
            + ", ".join(path.as_posix() for path in absent)
        )

    destinations: set[Path] = set()
    entries: list[PayloadEntry] = []
    for source, destination in mappings.items():
        if not safe_relative(destination) or forbidden(destination):
            raise RuntimeError(f"forbidden v2.0_pre destination: {destination.as_posix()}")
        if destination in destinations:
            raise RuntimeError(f"duplicate v2.0_pre destination: {destination.as_posix()}")
        destinations.add(destination)
        entries.append(PayloadEntry(source, destination))
    return sorted(entries, key=lambda item: item.destination.as_posix())


def require_clean_payload(entries: list[PayloadEntry]) -> None:
    command = ["git", "diff", "--quiet", "HEAD", "--"]
    command.extend(entry.source.as_posix() for entry in entries)
    result = subprocess.run(command, cwd=REPOSITORY, check=False)
    if result.returncode == 1:
        raise RuntimeError("selected v2.0_pre files differ from HEAD; commit them before packaging")
    if result.returncode != 0:
        raise RuntimeError("unable to verify v2.0_pre payload against HEAD")


def require_version(cmake_text: str) -> None:
    project = re.search(
        r"project\s*\(\s*wcns\s+VERSION\s+([0-9]+\.[0-9]+\.[0-9]+)",
        cmake_text,
        flags=re.IGNORECASE,
    )
    program = re.search(
        r'set\s*\(\s*WCNS_PROGRAM_VERSION\s+"([^"]+)"',
        cmake_text,
        flags=re.IGNORECASE,
    )
    if project is None or project.group(1) != PROJECT_VERSION:
        observed = "missing" if project is None else project.group(1)
        raise RuntimeError(f"CMake project version is {observed}; expected {PROJECT_VERSION}")
    if program is None or program.group(1) != PROGRAM_VERSION:
        observed = "missing" if program is None else program.group(1)
        raise RuntimeError(f"program version is {observed}; expected {PROGRAM_VERSION}")


def write_text(path: Path, value: str) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(value)


def parse_manifest(path: Path) -> dict[Path, str]:
    result: dict[Path, str] = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if match is None:
            raise RuntimeError(f"invalid PACKAGE_CONTENTS.sha256 line {line_number}")
        relative = Path(PurePosixPath(match.group(2)))
        if not safe_relative(relative) or forbidden(relative) or relative in result:
            raise RuntimeError(f"unsafe or duplicate manifest path on line {line_number}")
        result[relative] = match.group(1)
    return result


def verify_markdown_links(directory: Path) -> None:
    link_pattern = re.compile(r"\[[^\]]+\]\(([^)]+)\)")
    for document in directory.rglob("*.md"):
        for target in link_pattern.findall(document.read_text(encoding="utf-8")):
            target = target.strip()
            if not target or target.startswith("#") or "://" in target or target.startswith("mailto:"):
                continue
            relative_text = target.split("#", 1)[0]
            relative = Path(PurePosixPath(relative_text))
            resolved = (document.parent / relative).resolve()
            try:
                resolved.relative_to(directory)
            except ValueError as error:
                raise RuntimeError(
                    f"markdown link escapes package: {document.relative_to(directory)} -> {target}"
                ) from error
            if not resolved.is_file():
                raise RuntimeError(
                    f"broken package markdown link: {document.relative_to(directory)} -> {target}"
                )


def verify_directory(directory: Path) -> dict[str, object]:
    directory = directory.expanduser().resolve()
    if not directory.is_dir() or directory.name != PACKAGE_NAME:
        raise RuntimeError(f"expected a {PACKAGE_NAME} directory: {directory}")

    for path in directory.rglob("*"):
        relative = path.relative_to(directory)
        if path.is_symlink():
            raise RuntimeError(f"symbolic links are forbidden: {relative.as_posix()}")
        if forbidden(relative):
            raise RuntimeError(f"forbidden package entry: {relative.as_posix()}")

    manifest_path = directory / "PACKAGE_CONTENTS.sha256"
    revision_path = directory / "WCNS_SOURCE_REVISION"
    if not manifest_path.is_file() or not revision_path.is_file():
        raise RuntimeError("package lacks integrity metadata")
    manifest = parse_manifest(manifest_path)
    actual = {
        path.relative_to(directory)
        for path in directory.rglob("*")
        if path.is_file() and path != manifest_path
    }
    if set(manifest) != actual:
        missing = sorted(path.as_posix() for path in actual - set(manifest))
        extra = sorted(path.as_posix() for path in set(manifest) - actual)
        raise RuntimeError(f"content manifest mismatch; missing={missing}, extra={extra}")
    for relative, expected in manifest.items():
        observed = sha256(directory / relative)
        if observed != expected:
            raise RuntimeError(f"content hash mismatch: {relative.as_posix()}")

    revision = revision_path.read_text(encoding="ascii").strip()
    if not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise RuntimeError("invalid WCNS_SOURCE_REVISION")
    require_version((directory / "CMakeLists.txt").read_text(encoding="utf-8"))
    readme = (directory / "README.md").read_text(encoding="utf-8")
    if PROGRAM_VERSION not in readme:
        raise RuntimeError("package README does not identify v2.0_pre")
    verify_markdown_links(directory)
    return {
        "directory": str(directory),
        "revision": revision,
        "files": len(actual) + 1,
    }


def create_directory(destination: Path) -> Path:
    destination = destination.expanduser().resolve()
    if destination.name != PACKAGE_NAME:
        raise RuntimeError(f"destination directory must be named {PACKAGE_NAME}")
    if destination.exists():
        raise RuntimeError(f"refusing to overwrite existing destination: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)

    entries = tracked_payload()
    require_clean_payload(entries)
    require_version((REPOSITORY / "CMakeLists.txt").read_text(encoding="utf-8"))
    revision = run_git("rev-parse", "HEAD")

    with tempfile.TemporaryDirectory(prefix="wcns-v2-pre-", dir=destination.parent) as temporary:
        root = Path(temporary) / PACKAGE_NAME
        root.mkdir()
        for entry in entries:
            target = root / entry.destination
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(REPOSITORY / entry.source, target)
        write_text(root / "WCNS_SOURCE_REVISION", revision + "\n")
        files = sorted(
            (path for path in root.rglob("*") if path.is_file()),
            key=lambda path: path.relative_to(root).as_posix(),
        )
        lines = [
            f"{sha256(path)}  {path.relative_to(root).as_posix()}"
            for path in files
        ]
        write_text(root / "PACKAGE_CONTENTS.sha256", "\n".join(lines) + "\n")
        verify_directory(root)
        shutil.move(str(root), str(destination))

    result = verify_directory(destination)
    print(f"created: {result['directory']}")
    print(f"revision: {result['revision']}")
    print(f"files: {result['files']}")
    return destination


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--verify", type=Path, metavar="DIRECTORY")
    parser.add_argument(
        "--destination",
        type=Path,
        default=REPOSITORY.parent / PACKAGE_NAME,
        help=f"new directory to create; basename must be {PACKAGE_NAME}",
    )
    args = parser.parse_args()
    if args.verify is not None:
        result = verify_directory(args.verify)
        print(f"verified: {result['directory']}")
        print(f"revision: {result['revision']}")
        print(f"files: {result['files']}")
        return 0
    create_directory(args.destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
