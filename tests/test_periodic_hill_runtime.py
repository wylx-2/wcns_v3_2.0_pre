"""Bounded (<=3 steps/run) serial/MPI/restart acceptance for the hill module."""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import tempfile
import time

p=argparse.ArgumentParser()
p.add_argument('--runner',type=Path,required=True)
p.add_argument('--generator',type=Path,required=True)
p.add_argument('--config',type=Path,required=True)
p.add_argument('--mpi-runner',type=Path)
p.add_argument('--mpiexec',type=Path)
p.add_argument('--work',type=Path,required=True)
a=p.parse_args()
a.work.mkdir(parents=True,exist_ok=True)
work=Path(tempfile.mkdtemp(prefix='hill-',dir=a.work)).resolve()
evidence={'directory':str(work),'runs':{}}
def command(cmd,log,expected=0):
    start=time.monotonic()
    r=subprocess.run([str(x) for x in cmd],capture_output=True,text=True,timeout=180)
    (work/log).write_text(r.stdout+r.stderr,encoding='utf-8')
    if r.returncode!=expected:
        raise RuntimeError(f'{cmd}: return {r.returncode}; see {work/log}')
    return time.monotonic()-start,r.stdout+r.stderr
mesh=work/'tiny.cgns'
command([a.generator.resolve(),mesh,32,16,12],'mesh.log')
base=a.config.read_text(encoding='utf-8')
def replace(s,key,value):
    lines=s.splitlines(); found=False
    for n,line in enumerate(lines):
        if line.split('=',1)[0].strip()==key: lines[n]=f'{key} = {value}';found=True
    if not found: lines.append(f'{key} = {value}')
    return '\n'.join(lines)+'\n'
def run(name,steps=3,ranks=1,restart=None,extra=None):
    out=work/name
    s=base
    for k,v in {'mesh.path':mesh.as_posix(),'output.directory':out.as_posix(),
                'run.max_steps':steps,'run.t_end':.01}.items(): s=replace(s,k,v)
    if restart: s=replace(s,'restart.path',restart.as_posix())
    for k,v in (extra or {}).items(): s=replace(s,k,v)
    conf=work/(name+'.wcns');conf.write_text(s,encoding='utf-8')
    cmd=[a.runner.resolve(),'--config',conf]
    if ranks>1: cmd=[a.mpiexec.resolve(),'-n',ranks,a.mpi_runner.resolve(),'--config',conf]
    duration,log=command(cmd,name+'.log',expected=2)
    if 'reason=maximum_steps' not in log or 'numerical_failure' in log: raise RuntimeError(log[-2000:])
    evidence['runs'][name]={'ranks':ranks,'steps':steps,'seconds':duration}
    return out
def checkpoint(out): return next(out.glob('*.checkpoint.step*.cgns'))
def rows(out,suffix):
    with next(out.glob('*_'+suffix+'.csv')).open() as f:
        return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
def compare(left,right):
    for table in ('hill_mean_xy','hill_profiles','hill_wall'):
        l,r=rows(left,table),rows(right,table)
        assert len(l)==len(r)>0
        error=0
        for x,y in zip(l,r):
            for k in x:
                assert math.isfinite(x[k]) and math.isfinite(y[k])
                error=max(error,abs(x[k]-y[k]))
                assert math.isclose(x[k],y[k],rel_tol=2e-10,abs_tol=2e-10),(table,k,x[k],y[k])
        evidence.setdefault('max_abs_difference',{})[left.name+' vs '+right.name+' '+table]=error
    l,r=rows(left,'hill_history')[-1],rows(right,'hill_history')[-1]
    for k in l: assert math.isclose(l[k],r[k],rel_tol=2e-10,abs_tol=2e-10),(k,l[k],r[k])
full=run('serial-full')
part=run('serial-part',2)
resumed=run('serial-resume',restart=checkpoint(part))
compare(full,resumed)
if a.mpi_runner and a.mpiexec:
    parallel=run('mpi4-full',ranks=4);compare(full,parallel)
    resumed=run('mpi2-resume',ranks=2,restart=checkpoint(part));compare(full,resumed)
    resumed=run('serial-from-mpi',restart=checkpoint(parallel),steps=3);compare(full,resumed)
# A changed window must not silently mix two definitions of the average.
bad=replace((work/'serial-resume.wcns').read_text(),'hill.statistics.start',.0001)
bad=replace(bad,'output.directory',(work/'bad').as_posix())
(work/'bad.wcns').write_text(bad)
r=subprocess.run([str(a.runner.resolve()),'--config',str(work/'bad.wcns'),'--dry-run'],capture_output=True,text=True,timeout=60)
assert r.returncode!=0 and 'signature differs' in r.stdout+r.stderr
(work/'negative.log').write_text(r.stdout+r.stderr)
evidence['strict_restart_mismatch_rejected']=True
(work/'evidence.json').write_text(json.dumps(evidence,indent=2),encoding='utf-8')
print(json.dumps(evidence,indent=2))
