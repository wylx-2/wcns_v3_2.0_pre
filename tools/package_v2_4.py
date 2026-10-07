"""Package the reviewed v2.4 worktree snapshot and a separate periodic-hill case.

Never overwrites an existing destination. Does not stage, commit, or include unrelated
untracked work. Revision identifies the parent; SHA256 identifies the actual payload.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import package_v2_3 as base

ROOT=Path(__file__).resolve().parents[1]
NEW_FILES=(
    'include/wcns/physics/periodic_hill.hpp',
    'include/wcns/runtime/periodic_hill.hpp',
    'src/runtime/periodic_hill.cpp',
    'tools/generate_periodic_hill_cgns.cpp',
    'docs/periodic-hill-v2.4.md',
    'docs/release-notes-2.4.md',
    'docs/v2.4-validation.md',
)
CASE_FILES=('README.md','smoke.wcns','production.wcns','production_restart.wcns',
            'grids/hill_88x48x24.cgns','grids/hill_mesh_xy.png','mesh-inspection.txt')
def manifest(root):
    paths=sorted(p for p in root.rglob('*') if p.is_file() and p.name!='PACKAGE_CONTENTS.sha256')
    (root/'PACKAGE_CONTENTS.sha256').write_text(''.join(
        hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(root).as_posix()+'\n' for p in paths),encoding='utf-8')
def verify(root):
    listed=set()
    for line in (root/'PACKAGE_CONTENTS.sha256').read_text(encoding='utf-8').splitlines():
        expected,relative=line.split('  ',1);p=(root/relative).resolve();p.relative_to(root.resolve())
        assert hashlib.sha256(p.read_bytes()).hexdigest()==expected,relative
        listed.add(p)
    actual={p.resolve() for p in root.rglob('*') if p.is_file() and p.name!='PACKAGE_CONTENTS.sha256'}
    assert actual==listed,'unexpected or missing files'
    return len(listed)
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-parent',type=Path)
    parser.add_argument('--verify',type=Path)
    a=parser.parse_args()
    if a.verify:
        print(json.dumps({'files_verified':verify(a.verify)},indent=2));return
    if a.output_parent is None: parser.error('--output-parent or --verify required')
    out=a.output_parent.resolve();source=out/'WCNS_v2.4';case=out/'WCNS_v2.4_periodic_hill_case'
    for p in (source,case,Path(str(source)+'.zip'),Path(str(case)+'.zip')):
        if p.exists(): raise RuntimeError(f'refusing to overwrite {p}')
    tracked=subprocess.check_output(['git','ls-files','-z'],cwd=ROOT).decode().split('\0')
    mapping=dict(base.EXPLICIT_MAPPINGS)
    mapping.pop(Path('docs/v2.3-package-readme.md'))
    mapping[Path('docs/v2.4-package-readme.md')]=Path('README.md')
    for name in tracked:
        p=Path(name)
        if name and any(p==prefix or prefix in p.parents for prefix in base.DIRECTORY_PREFIXES): mapping[p]=p
    for p in base.TOOL_SOURCES: mapping[p]=p
    for name in NEW_FILES: mapping[Path(name)]=Path(name)
    for p in mapping:
        if not (ROOT/p).is_file(): raise RuntimeError(f'missing source {p}')
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()+'-v2.4-worktree'
    for src,dst in mapping.items():
        target=source/dst;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/src,target)
    (source/'WCNS_SOURCE_REVISION').write_text(revision+'\n')
    base.verify_markdown_links(source)
    case_root=ROOT/'cases/manual/case08_periodic_hill'
    for name in CASE_FILES:
        target=case/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(case_root/name,target)
    shutil.copytree(case_root/'validation',case/'validation')
    for root in (source,case):
        manifest(root);count=verify(root)
        archive=shutil.make_archive(str(root),'zip',root_dir=out,base_dir=root.name)
        print(json.dumps({'directory':str(root),'files_verified':count,'archive':archive,'source_revision':revision}))
if __name__=='__main__':main()
