"""Verify v2.4 version and the required case/release inputs."""
import argparse
import re
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--repository',type=Path,default=Path(__file__).resolve().parents[1]);a=p.parse_args()
r=a.repository
s=(r/'CMakeLists.txt').read_text(encoding='utf-8')
assert re.search(r'project\(wcns VERSION 2\.4\.0',s)
assert 'set(WCNS_PROGRAM_VERSION "2.4"' in s
for f in ('include/wcns/runtime/periodic_hill.hpp','src/runtime/periodic_hill.cpp',
          'tools/generate_periodic_hill_cgns.cpp','tools/package_v2_4.py','docs/periodic-hill-v2.4.md',
          'docs/release-notes-2.4.md','docs/v2.4-validation.md'):
    assert (r/f).is_file(),f
case=r/'cases/manual/case08_periodic_hill'
assert (case/'grids/hill_88x48x24.cgns').stat().st_size>2000000
for name in ('smoke.wcns','production.wcns','production_restart.wcns'):
    text=(case/name).read_text()
    for fragment in ('turbulence.model = none','hill.enabled = true','initial.type = periodic_hill'):
        assert fragment in text,(name,fragment)
assert 'run.max_steps = 3' in (case/'smoke.wcns').read_text()
print('v2.4 release inputs and bounded smoke contract: PASS')
