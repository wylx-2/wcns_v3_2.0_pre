"""Package the reviewed v2.5 worktree and two independent chapter-5 case bundles.
Refuses existing destinations; includes only explicit relevant source/artifact paths.
No staging, commit, production execution, or unrelated user output collection.
"""
import argparse,json,shutil,subprocess
from pathlib import Path
import package_v2_3 as base
from package_v2_4 import manifest,verify,NEW_FILES as V24_FILES
ROOT=Path(__file__).resolve().parents[1]
NEW_FILES=(*V24_FILES,
 'include/wcns/physics/chapter5.hpp','include/wcns/physics/ramp_inlet_table.inc',
 'include/wcns/runtime/chapter5.hpp','src/runtime/chapter5.cpp',
 'tools/generate_chapter5_cgns.cpp','tools/audit_chapter5_metrics.cpp',
 'tools/prepare_chapter5_meshes.py','tools/generate_ramp_inlet.py',
 'tools/analyze_chapter5.py','tools/digitize_chapter5_reference.py',
 'tools/plot_chapter5_meshes.py','tools/compare_chapter5_statistics.py',
 'docs/chapter5-v2.5-research.md','docs/chapter5-v2.5-implementation.md',
 'docs/v2.5-validation.md','docs/chapter5-source-provenance.json',
 'docs/chapter5-assets/sd7003_mesh.png','docs/chapter5-assets/ramp_mesh.png')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-parent',type=Path);p.add_argument('--verify',type=Path);a=p.parse_args()
 if a.verify:print(json.dumps({'files_verified':verify(a.verify)}));return
 if not a.output_parent:p.error('--output-parent or --verify required')
 out=a.output_parent.resolve();source=out/'WCNS_v2.5'
 cases=[('case09_sd7003',out/'WCNS_v2.5_SD7003_case'),('case10_compression_ramp',out/'WCNS_v2.5_compression_ramp_case')]
 roots=[source]+[dest for _,dest in cases]
 for root in roots:
  for path in [root,Path(str(root)+'.zip')]:
   if path.exists():raise RuntimeError(f'refusing to overwrite {path}')
 tracked=subprocess.check_output(['git','ls-files','-z'],cwd=ROOT).decode().split('\0')
 mapping=dict(base.EXPLICIT_MAPPINGS);mapping.pop(Path('docs/v2.3-package-readme.md'))
 mapping[Path('docs/v2.5-package-readme.md')]=Path('README.md')
 mapping[Path('docs/v2.5-package-validation.txt')]=Path('package-validation.txt')
 for name in tracked:
  path=Path(name)
  if name and any(path==prefix or prefix in path.parents for prefix in base.DIRECTORY_PREFIXES):mapping[path]=path
 for path in base.TOOL_SOURCES:mapping[path]=path
 for name in NEW_FILES:mapping[Path(name)]=Path(name)
 for original,target in mapping.items():
  dest=source/target;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/original,dest)
 revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()+'-v2.5-worktree'
 (source/'WCNS_SOURCE_REVISION').write_text(revision+'\n');base.verify_markdown_links(source)
 for directory,dest in cases:
  origin=ROOT/'cases/manual'/directory;dest.mkdir(parents=True)
  for name in ['README.md','production.wcns','production_restart.wcns','smoke.wcns','package-validation.txt']:shutil.copyfile(origin/name,dest/name)
  for name in ['grids','reference','report','tools']:shutil.copytree(origin/name,dest/name)
  validation=dest/'validation';validation.mkdir()
  allowed=['mesh-generation.log','metrics.log','field-checks.log','serial-mpi4.json','repartition-restart.json',
           'final-serial.log','final-mpi4.log','final-part.log','final-resume2.log','unit-results.log',
           'postprocess-test.log','production-config-check.log','release-regressions.log','mesh-regeneration.json']
  for name in allowed:
   if (origin/'validation'/name).is_file():shutil.copyfile(origin/'validation'/name,validation/name)
  shutil.copytree(origin/'validation/configs',validation/'configs')
  base.verify_markdown_links(dest)
 for root in roots:
  manifest(root);count=verify(root);archive=shutil.make_archive(str(root),'zip',root_dir=out,base_dir=root.name)
  print(json.dumps({'directory':str(root),'files_verified':count,'archive':archive,'source_revision':revision}))
if __name__=='__main__':main()
