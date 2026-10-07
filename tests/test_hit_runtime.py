"""Bounded HIT acceptance: 3-step MPI/restart/forcing and exact spectral identities."""
import argparse,csv,json,math,subprocess,tempfile,time
from pathlib import Path

def replace(text,key,value):
    lines=text.splitlines();found=False
    for i,line in enumerate(lines):
        if line.split('=',1)[0].strip()==key:lines[i]=f'{key} = {value}';found=True
    if not found:lines.append(f'{key} = {value}')
    return '\n'.join(lines)+'\n'

def main():
    p=argparse.ArgumentParser();p.add_argument('--runner',type=Path,required=True)
    p.add_argument('--repository',type=Path,required=True);p.add_argument('--work',type=Path,required=True)
    p.add_argument('--mpiexec',type=Path);p.add_argument('--generator',type=Path)
    a=p.parse_args();a.work.mkdir(parents=True,exist_ok=True)
    work=Path(tempfile.mkdtemp(prefix='hit-',dir=a.work)).resolve();evidence={'work':str(work),'runs':{},'comparisons':{}}
    def rows(out,name):
        with next(out.glob('*_hit_'+name+'.csv')).open() as f:return list(csv.DictReader(f))
    def numeric_compare(left,right,label):
        error=0
        for name in ['history','spectrum','spectrum_1d','correlation','means','spectrum_mean']:
            x,y=rows(left,name),rows(right,name)
            if name in ['history','spectrum','spectrum_1d','correlation']:
                last=max(float(v['step']) for v in x);x=[v for v in x if float(v['step'])==last];y=[v for v in y if float(v['step'])==last]
            assert len(x)==len(y)>0,(name,len(x),len(y))
            for r,s in zip(x,y):
                for k in r:
                    if k=='quantity':assert r[k]==s[k];continue
                    u,v=float(r[k]),float(s[k])
                    if math.isnan(u) and math.isnan(v):continue
                    scaled=abs(u-v)/max(1,abs(u),abs(v));error=max(error,scaled)
                    assert scaled<3e-8,(label,name,k,u,v,scaled)
        # Compare every final exported cell value, including spatial coordinates.
        def field(out):
            file=sorted(out.glob('*.field.step*.dat'))[-1]
            values=[]
            for line in file.read_text().splitlines():
                if line.startswith(('TITLE','VARIABLES','ZONE')):continue
                values.extend(float(x) for x in line.split())
            return values
        x,y=field(left),field(right);assert len(x)==len(y)>0
        for u,v in zip(x,y):
            scaled=abs(u-v)/max(1,abs(u),abs(v));error=max(error,scaled)
            assert scaled<3e-8,(label,'field',u,v)
        evidence['comparisons'][label]=error
    def run(case,name,ranks=1,steps=3,restart=None,extra=None,expected=2):
        root=a.repository/'cases/manual'/('case11_hit_decay' if case=='decay' else 'case12_hit_forced')
        text=(root/'smoke.wcns').read_text();out=work/name
        config={'mesh.path':(root/'grids/hit16.cgns').resolve().as_posix(),
                'hit.spectrum_file':(root/'reference'/('cbc_spectrum_42_nondimensional.dat' if case=='decay' else 'jhtdb_spectrum.dat')).resolve().as_posix(),
                'output.directory':out.as_posix(),'run.max_steps':steps,'output.field.format':'tecplot'}
        if restart:config['restart.path']=restart.as_posix()
        config.update(extra or {})
        for k,v in config.items():text=replace(text,k,v)
        path=work/(name+'.wcns');path.write_text(text)
        cmd=[str(a.runner.resolve()),'--config',str(path)]
        if ranks>1:cmd=[str(a.mpiexec),'-n',str(ranks)]+cmd
        start=time.monotonic();r=subprocess.run(cmd,capture_output=True,text=True,timeout=150)
        (work/(name+'.log')).write_text(r.stdout+r.stderr)
        assert r.returncode==expected,(name,r.returncode,r.stdout[-500:],r.stderr[-500:])
        evidence['runs'][name]={'ranks':ranks,'steps':steps,'seconds':time.monotonic()-start}
        if expected!=2:return out
        history=rows(out,'history')
        for row in history:
            assert abs(float(row['parseval_error']))<1e-11
            assert min(float(row['rho_mean']),float(row['T_mean']))>0
        assert float(history[0]['divergence_rms'])<1e-10 if not restart else True
        if case=='decay':assert float(history[-1]['K'])<float(history[0]['K']) if not restart else True
        if case=='forced' and (extra or {}).get('hit.forcing','jhtdb_shells')=='jhtdb_shells':
            for row in history:assert abs(float(row['K_forced_band'])-.43)<2e-11
        if (extra or {}).get('hit.forcing')=='constant_power':
            for row in history:assert abs(float(row['forcing_power'])-.103)<1e-11
        return out
    for case in ['decay','forced']:
        full=run(case,case+'-serial');part=run(case,case+'-part',steps=1)
        checkpoint=next(part.glob('*.checkpoint.step*.cgns'))
        resumed=run(case,case+'-resume',restart=checkpoint);numeric_compare(full,resumed,case+' serial restart')
        if a.mpiexec:
            parallel=run(case,case+'-mpi4',ranks=4);numeric_compare(full,parallel,case+' MPI4')
            resumed=run(case,case+'-resume2',ranks=2,restart=checkpoint);numeric_compare(full,resumed,case+' repartition MPI2')
        run(case,case+'-bad-n',extra={'hit.n':32},expected=1)
        run(case,case+'-bad-seed',restart=checkpoint,extra={'hit.seed':17},expected=1)
    run('forced','constant-power',extra={'hit.forcing':'constant_power'})
    run('forced','constant-band',extra={'hit.forcing':'constant_band_energy'})
    run('forced','bad-forced-cutoff',extra={'hit.cutoff':2},expected=1)
    if a.generator and a.mpiexec:
        mesh=work/'one-zone.cgns'
        subprocess.run([str(a.generator.resolve()),str(mesh),'16','1','6.283185307179586'],check=True,timeout=30,capture_output=True)
        extra={'mesh.path':mesh.as_posix()}
        full=run('forced','one-zone-serial',extra=extra)
        parallel=run('forced','one-zone-split4',ranks=4,extra=extra)
        numeric_compare(full,parallel,'one-zone auto split MPI4')
    (work/'evidence.json').write_text(json.dumps(evidence,indent=2));print(json.dumps(evidence,indent=2))

if __name__=='__main__':main()
