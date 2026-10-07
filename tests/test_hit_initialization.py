"""Bounded 16^3 HIT initialization/preparation acceptance, including MPI/restarts.

Uses NumPy only as an independent DFT oracle. Never runs the production preparation.
"""
import argparse, csv, json, math, re, subprocess, time
from pathlib import Path
import numpy as np

def replace(text, key, value):
    if value=='':return '\n'.join(line for line in text.splitlines() if line.split('=',1)[0].strip()!=key)+'\n'
    lines=text.splitlines(); found=False
    for i,line in enumerate(lines):
        if line.split('=',1)[0].strip()==key:
            lines[i]=f'{key} = {value}'; found=True
    if not found: lines.append(f'{key} = {value}')
    return '\n'.join(lines)+'\n'

def rows(path):
    with path.open() as f:return list(csv.DictReader(f))

def history(out): return rows(next(out.glob('*_hit_history.csv')))

def field(out, first=False):
    paths=sorted(out.glob('*.field.step*.dat')); path=paths[0] if first else paths[-1]
    zones=re.split(r'^ZONE .*\n',path.read_text(),flags=re.M)[1:]
    data=np.concatenate([np.fromstring(z,sep=' ').reshape(9,-1).T for z in zones])
    # Metric coordinates carry roundoff: sort logical cells, not raw floats.
    index=np.rint((data[:,:3]-data[:,:3].min(axis=0))/(data[:,:3].max(axis=0)-data[:,:3].min(axis=0))*15).astype(int)
    order=np.lexsort((index[:,2],index[:,1],index[:,0])); data=data[order]
    assert data.shape==(16**3,9)
    return data.reshape(16,16,16,9)

def fft(data):return np.fft.fftn(data[...,4:7],axes=(0,1,2))/16**3

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runner',type=Path,required=True);p.add_argument('--repository',type=Path,required=True)
    p.add_argument('--work',type=Path,required=True);p.add_argument('--mpiexec',type=Path)
    a=p.parse_args();a.work.mkdir(parents=True,exist_ok=False);a.work=a.work.resolve()
    evidence={'runs':{},'checks':{}}
    def run(kind,name,ranks=1,steps=2,restart=None,extra=None,expected=2):
        root=a.repository.resolve()/'cases/manual'/('case11_hit_decay' if kind=='prepared' else 'case12_hit_forced')
        text=(root/'smoke.wcns').read_text();out=a.work/name
        settings={'case.name':'hit-init-test','mesh.path':(root/'grids/hit16.cgns').as_posix(),
            'output.directory':out.as_posix(),'output.field.format':'tecplot','output.field.write_initial':'true',
            'output.checkpoint.write_initial':'true','output.checkpoint.every_steps':1,
            'run.max_steps':steps,'run.t_end':1,'run.max_wall_time':90,'hit.type':'decay',
            'hit.forcing':'constant_power','hit.thermostat':'false','hit.sample_every_steps':1,
            'hit.sample_times':'0','hit.statistics.start':0,'hit.statistics.end':1}
        if kind=='prepared':
            settings.update({'hit.spectrum_file':(root/'reference/cbc_spectrum_42_nondimensional.dat').as_posix(),
                'hit.initialization':'shell_spectrum','hit.preparation.time':.001,'hit.preparation.max_steps':20})
        else:
            settings.update({'hit.spectrum_file':'','hit.initialization':'analytic_random_phase',
                'hit.spectrum_amplitude':.00013,'hit.peak_wave':8})
        if restart:settings['restart.path']=restart.as_posix()
        settings.update(extra or {})
        for k,v in settings.items():text=replace(text,k,v)
        path=a.work/(name+'.wcns');path.write_text(text)
        cmd=[str(a.runner.resolve()),'--config',str(path)]
        if ranks>1:cmd=[str(a.mpiexec),'-n',str(ranks)]+cmd
        start=time.monotonic();r=subprocess.run(cmd,capture_output=True,text=True,timeout=120)
        (a.work/(name+'.log')).write_text(r.stdout+r.stderr)
        assert r.returncode==expected,(name,r.returncode,r.stdout[-800:],r.stderr[-800:])
        evidence['runs'][name]={'ranks':ranks,'return_code':r.returncode,'seconds':time.monotonic()-start}
        return out
    def compare(x,y,label):
        error=float(np.max(np.abs(field(x)-field(y))/np.maximum(1,np.abs(field(x)))))
        assert error<2e-10,(label,error)
        for k,v in history(x)[-1].items():
            left,right=float(v),float(history(y)[-1][k])
            if math.isnan(left) and math.isnan(right):continue
            assert abs(left-right)<2e-8*max(1,abs(left)),(label,k,left,right)
        evidence['checks'][label]=error
    analytic=run('analytic','analytic')
    initial=history(analytic)[0];K0=3*.00013/64*math.sqrt(2*math.pi)*8**5
    assert abs(float(initial['K'])-K0)<2e-13
    assert abs(float(initial['Mt'])-math.sqrt(2*K0)*.1)<2e-13
    assert max(abs(float(initial[k])) for k in ['u_mean','v_mean','w_mean','rho_rms','T_rms','divergence_rms'])<2e-12
    spec=fft(field(analytic,True));m=np.stack(np.meshgrid(*[np.fft.fftfreq(16)*16]*3,indexing='ij'),axis=-1)
    assert np.max(np.abs(np.sum(m*spec,axis=-1)))<1e-13
    assert abs(.5*np.sum(np.abs(spec)**2)-K0)<2e-13
    assert np.max(np.abs(spec[np.linalg.norm(m,axis=-1)>16/3]))<1e-14
    evidence['checks']['analytic_K0']=K0
    different=run('analytic','different-seed',extra={'hit.seed':1024})
    assert np.max(np.abs(field(analytic,True)-field(different,True)))>.1
    part=run('analytic','analytic-part',steps=1)
    resumed=run('analytic','analytic-resume',restart=sorted(part.glob('*.checkpoint.step*.cgns'))[-1])
    compare(analytic,resumed,'analytic strict restart')
    prepared=run('prepared','prepared')
    prep_history=history(prepared/'preparation');formal=history(prepared)
    assert abs(float(prep_history[-1]['time'])-.001)<1e-15
    assert all(float(r['statistics_weight'])==0 for r in prep_history)
    assert float(formal[0]['time'])==0 and float(formal[0]['step'])==0 and float(formal[0]['statistics_weight'])==0
    assert float(formal[0]['divergence_rms'])<2e-11
    audit=rows(next((prepared/'preparation').glob('*_rematch.csv')))
    assert max(abs(float(r['E_after'])-float(r['E_target_discrete'])) for r in audit)<1e-14
    # Independent per-mode (NOT shell) target, including author's low-k log extrapolation.
    target=np.loadtxt(a.repository/'cases/manual/case11_hit_decay/reference/cbc_spectrum_42_nondimensional.dat')
    km=np.linalg.norm(m,axis=-1);active=(km>0)&(km<=16/3);k=km[active]*2*np.pi
    logE=np.interp(np.log(k),np.log(target[:,0]),np.log(target[:,1]))
    low=k<target[0,0];slope=np.log(target[1,1]/target[0,1])/np.log(target[1,0]/target[0,0])
    logE[low]=np.log(target[0,1])+slope*np.log(k[low]/target[0,0])
    desired=2*np.exp(logE)*(2*np.pi)/(4*np.pi*km[active]**2)
    prepared_fft=fft(field(prepared,True));actual=np.sum(np.abs(prepared_fft[active])**2,axis=-1)
    assert np.max(np.abs(actual-desired)/desired)<1e-11
    evidence['checks']['per_mode_target_relative_error']=float(np.max(np.abs(actual-desired)/desired))
    before=field(prepared/'preparation');after=field(prepared,True)
    assert np.max(np.abs(before[...,[3,7,8]]-after[...,[3,7,8]]))<2e-12
    oldfft=fft(before);old=oldfft[active]
    transverse=old-m[active]*np.sum(m[active]*old,axis=-1)[:,None]/km[active,None]**2
    # Positive real multiplier on each transverse mode preserves phase correlations.
    ratio=np.sum(np.conj(transverse)*prepared_fft[active],axis=-1)/np.sum(abs(transverse)**2,axis=-1)
    assert np.max(abs(ratio.imag))<2e-12 and min(ratio.real)>0
    assert np.max(abs(prepared_fft[active]-ratio[:,None]*transverse))<1e-13
    evidence['checks']['rematch_preserves_thermodynamics_and_transverse_phases']=True
    # A zero-step formal checkpoint proves phase/time reset survives restart.
    cp0=sorted(prepared.glob('*.checkpoint.step*.cgns'))[0]
    resume0=run('prepared','prepared-resume0',restart=cp0)
    compare(prepared,resume0,'restart after rematch does not repeat preparation')
    assert not (resume0/'preparation').exists()
    interrupted=run('prepared','prep-interrupted',extra={'hit.preparation.max_steps':1})
    assert not list(interrupted.glob('*_hit_history.csv'))
    cp=sorted((interrupted/'preparation').glob('*.checkpoint.step*.cgns'))[-1]
    resumed=run('prepared','prep-resume',restart=cp)
    compare(prepared,resumed,'restart during preparation')
    # A checkpoint exactly at the preparation boundary is still marked pending.
    endcp=sorted((prepared/'preparation').glob('*.checkpoint.step*.cgns'))[-1]
    resumed=run('prepared','prep-end-resume',restart=endcp)
    compare(prepared,resumed,'restart at preparation endpoint')
    if a.mpiexec:
        mpi=run('analytic','analytic-mpi4',ranks=4);compare(analytic,mpi,'analytic MPI4')
        mpi=run('prepared','prepared-mpi4',ranks=4);compare(prepared,mpi,'prepared MPI4')
        mpi=run('prepared','prep-resume-mpi2',ranks=2,restart=cp);compare(prepared,mpi,'preparation restart repartition MPI2')
    for name,kind,settings,restart in [
        ('bad-mode','analytic',{'hit.initialization':'typo'},None),
        ('conflicting-energy','analytic',{'hit.initial_energy':1},None),
        ('negative-A','analytic',{'hit.spectrum_amplitude':-1},None),
        ('analytic-prep','analytic',{'hit.preparation.time':.01},None),
        ('forced-prep','prepared',{'hit.type':'forced'},None),
        ('changed-duration','prepared',{'hit.preparation.time':.0003},cp0),
        ('algorithm-change','prepared',{'restart.mode':'algorithm_change'},cp0),
        ('changed-A','analytic',{'hit.spectrum_amplitude':.001},sorted(part.glob('*.checkpoint.step*.cgns'))[-1]),
    ]:run(kind,name,extra=settings,restart=restart,expected=1)
    (a.work/'evidence.json').write_text(json.dumps(evidence,indent=2)+'\n')
    print(json.dumps(evidence,indent=2))

if __name__=='__main__':main()
