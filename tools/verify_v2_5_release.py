"""Verify v2.5 release inputs; this does not launch production computations."""
import argparse,json,re
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--repository',type=Path,default=Path(__file__).resolve().parents[1]);r=p.parse_args().repository
s=(r/'CMakeLists.txt').read_text(encoding='utf8')
assert re.search(r'project\(wcns VERSION 2\.5\.0',s)
assert 'set(WCNS_PROGRAM_VERSION "2.5"' in s
for file in ['include/wcns/runtime/chapter5.hpp','include/wcns/physics/chapter5.hpp',
             'include/wcns/physics/ramp_inlet_table.inc','src/runtime/chapter5.cpp',
             'tools/analyze_chapter5.py','tools/prepare_chapter5_meshes.py',
             'tools/generate_chapter5_cgns.cpp','docs/chapter5-v2.5-research.md','docs/v2.5-validation.md']:
    assert (r/file).is_file(),file
for directory,kind,stem,ncells in [('case09_sd7003','sd7003','sd7003',1720320),('case10_compression_ramp','compression_ramp','ramp',8847360)]:
    c=r/'cases/manual'/directory
    assert json.loads((c/f'grids/{stem}_coarse.json').read_text())['cells']==ncells
    assert (c/f'grids/{stem}_coarse.cgns').stat().st_size>ncells*20
    for name in ['production.wcns','production_restart.wcns','smoke.wcns']:
        s=(c/name).read_text()
        for item in ['turbulence.model = none',f'initial.type = {kind}',f'benchmark.type = {kind}']:
            assert item in s,(name,item)
    assert 'run.max_steps = 3' in (c/'smoke.wcns').read_text()
    assert (c/'reference/Cp_digitized.csv').is_file()
print('v2.5 source/case/version contract: PASS (no production run)')
