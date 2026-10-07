"""Stage v2.6 source and HIT cases, then finalize only after independent validation.

Never overwrites previous releases. Stage copies only named source trees/artifacts.
The finalizer updates checksum manifests and refuses existing ZIP archives.
"""
import argparse,hashlib,json,shutil,subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
NAMES=['WCNS_v2.6','WCNS_v2.6_HIT_decay_case','WCNS_v2.6_HIT_forced_case']

def manifest(root):
    files=sorted(p for p in root.rglob('*') if p.is_file() and p.name!='PACKAGE_CONTENTS.sha256')
    (root/'PACKAGE_CONTENTS.sha256').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(root).as_posix()+'\n' for p in files),encoding='utf-8')

def verify(root):
    files=set()
    for line in (root/'PACKAGE_CONTENTS.sha256').read_text(encoding='utf-8').splitlines():
        digest,relative=line.split('  ',1);path=(root/relative).resolve();path.relative_to(root.resolve())
        if hashlib.sha256(path.read_bytes()).hexdigest()!=digest:raise ValueError('checksum mismatch: '+relative)
        files.add(path)
    actual={p.resolve() for p in root.rglob('*') if p.is_file() and p.name!='PACKAGE_CONTENTS.sha256'}
    if files!=actual:raise ValueError('unexpected or missing files')
    return len(files)

def copy(src,dest):
    dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dest)

def stage(out):
    import package_v2_3 as base
    from package_v2_5 import NEW_FILES
    for name in NAMES:
        if (out/name).exists() or (out/(name+'.zip')).exists():raise ValueError('refusing existing destination: '+name)
    source=out/NAMES[0]
    mapping=dict(base.EXPLICIT_MAPPINGS);mapping.pop(Path('docs/v2.3-package-readme.md'))
    mapping[Path('docs/v2.6-package-readme.md')]=Path('README.md')
    for directory in ['include','src','third_party/cgns','third_party/fftw','tests']:
        for path in (ROOT/directory).rglob('*'):
            if path.is_file() and '__pycache__' not in path.parts and path.suffix not in ['.pyc','.exe']:
                relative=path.relative_to(ROOT);mapping[relative]=relative
    for p in base.TOOL_SOURCES:mapping[p]=p
    for name in NEW_FILES:mapping[Path(name)]=Path(name)
    for name in ['generate_hit_cgns.cpp','generate_test_cgns.cpp','analyze_hit.py','prepare_hit_reference.py','prepare_hit_restart.py','package_v2_6.py']:
        p=Path('tools')/name;mapping[p]=p
    for path in (ROOT/'docs/hit-v2.6').rglob('*'):
        if path.is_file() and 'qa' not in path.parts and path.suffix in ['.md','.docx','.png']:
            relative=path.relative_to(ROOT);mapping[relative]=relative
    for original,dest in mapping.items():copy(ROOT/original,source/dest)
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()+'-v2.6-worktree'
    (source/'WCNS_SOURCE_REVISION').write_text(revision+'\n')
    base.verify_markdown_links(source)
    for directory,name in zip(['case11_hit_decay','case12_hit_forced'],NAMES[1:]):
        origin=ROOT/'cases/manual'/directory;dest=out/name
        for path in origin.rglob('*'):
            relative=path.relative_to(origin)
            if not path.is_file() or any(x in relative.parts for x in ['output','qa','__pycache__']):continue
            if len(relative.parts)==1 or relative.parts[0] in ['grids','reference','report','tools','validation']:
                copy(path,dest/relative)
        base.verify_markdown_links(dest)
    print(json.dumps({'staged':[str(out/n) for n in NAMES],'revision':revision},indent=2))

def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--stage',type=Path);g.add_argument('--finalize',type=Path);g.add_argument('--verify',type=Path);a=p.parse_args()
    if a.verify:print(json.dumps({'verified':verify(a.verify)}));return
    if a.stage:stage(a.stage.resolve());return
    out=a.finalize.resolve()
    for name in NAMES:
        if not (out/name/'package-validation.txt').is_file():raise ValueError('missing package validation record')
        if (out/(name+'.zip')).exists():raise ValueError('refusing to overwrite archive')
    summary=[]
    for name in NAMES:
        root=out/name;manifest(root);count=verify(root)
        archive=Path(shutil.make_archive(str(root),'zip',root_dir=out,base_dir=name))
        summary.append({'directory':str(root),'files':count,'zip':str(archive),'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest()})
    (out/'WCNS_v2.6_delivery.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
