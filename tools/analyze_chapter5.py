"""Postprocess WCNS v2.5 accepted-time/span moments and fixed-probe histories.
Requires numpy, scipy and matplotlib. All quantities are nondimensional.
Usage: python analyze_chapter5.py --case sd7003 --prefix output/.../case --output analysis
Use --start to discard startup from history spectra (moments have their own runtime window).
"""
from pathlib import Path
import argparse,csv,json
import numpy as np
from scipy.signal import welch,coherence
from scipy.integrate import cumulative_trapezoid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as tri

def read(path):
    a=np.genfromtxt(path,delimiter=',',names=True,encoding='utf8')
    return np.atleast_1d(a)

def csvout(path,names,rows):
    with Path(path).open('w',newline='',encoding='utf8') as f:
        w=csv.writer(f);w.writerow(names);w.writerows(rows)

def crossings(x,y):
    return [float(x[i]-y[i]*(x[i+1]-x[i])/(y[i+1]-y[i])) for i in range(len(x)-1)
            if y[i]*y[i+1]<0]

def uniform_history(t,v):
    order=np.argsort(t);t=np.asarray(t)[order];v=np.asarray(v)[order]
    unique=np.r_[True,np.diff(t)>1e-14];t=t[unique];v=v[unique]
    if len(t)<128: return None
    gaps=np.diff(t);dt=float(np.max(gaps))
    if np.max(gaps)>2*np.median(gaps): raise ValueError('history gap exceeds twice the median; split the record into contiguous windows')
    grid=np.arange(t[0],t[-1]+.01*dt,dt)
    return grid,np.interp(grid,t,v),dt

def spectrum(t,v):
    data=uniform_history(t,v)
    if data is None:return None
    grid,values,dt=data
    f,p=welch(values,fs=1/dt,nperseg=min(2048,len(values)//4),noverlap=None,
              window='hann',detrend='constant',scaling='density')
    return f,p,dict(dt=dt,duration=float(grid[-1]-grid[0]),nyquist=.5/dt,
                   frequency_spacing=float(f[1]-f[0]),method='linear resampling at maximum observed dt; Hann Welch 50% overlap')

def profiles(mean,wall,case,destination,span_cells=None):
    # Explicit structured triangles prevent interpolation through the airfoil.
    use=wall[wall['ny']>0] if case=='sd7003' else wall
    use=np.sort(use,order='x')
    fields=['rho','u','v','w','p','T','uu','vv','ww','uv','uw','vw',
            'favre_u','favre_v','favre_w','favre_uu','favre_vv','favre_ww','favre_uv','favre_uw','favre_vw']
    points=[];values={f:[] for f in fields};cols=[]
    for w in use:
        rows=mean[(mean['zone']==w['zone']) & (mean['i']==w['i'])];rows=np.sort(rows,order='j')
        ids=[]
        for j in range(len(rows)+1):
            ids.append(len(points))
            if j==0:
                points.append([w['x'],w['y']])
                for f in fields:values[f].append(w['rho_wall'] if f=='rho' else w['Twall'] if f=='T' else w['p'] if f=='p' else 0)
            else:
                points.append([rows[j-1]['x'],rows[j-1]['y']])
                for f in fields:values[f].append(rows[j-1][f])
        cols.append(ids)
    triangles=[]
    for a,b in zip(cols[:-1],cols[1:]):
        if len(a)!=len(b):raise ValueError('profile postprocessor requires equal radial counts across source zones')
        for j in range(len(a)-1):triangles.extend([[a[j],b[j],b[j+1]],[a[j],b[j+1],a[j+1]]])
    p=np.array(points);triangles=np.array(triangles)
    e1=p[triangles[:,1]]-p[triangles[:,0]];e2=p[triangles[:,2]]-p[triangles[:,0]]
    cross=e1[:,0]*e2[:,1]-e1[:,1]*e2[:,0]
    triangles[cross<0]=triangles[cross<0][:,[0,2,1]]
    mesh=tri.Triangulation(p[:,0],p[:,1],triangles)
    interpolators={f:tri.LinearTriInterpolator(mesh,np.asarray(v)) for f,v in values.items()}
    targets=[.1,.5,.95] if case=='sd7003' else [70,80,90,100,110,120,130]
    edge=.2 if case=='sd7003' else 5.
    distances=np.r_[0,np.geomspace(1e-5 if case=='sd7003' else 1e-4,edge,300)]
    out=[];boundary=[]
    for target in targets:
        wd={f:float(np.interp(target,use['x'],use[f])) for f in use.dtype.names}
        nx,ny=wd['nx'],wd['ny'];norm=np.hypot(nx,ny);nx/=norm;ny/=norm
        tx,ty=ny,-nx
        x=target+distances*nx;y=wd['y']+distances*ny
        data={f:np.ma.filled(interp(x,y),np.nan) for f,interp in interpolators.items()}
        for f in fields:data[f][0]=wd['rho_wall'] if f=='rho' else wd['Twall'] if f=='T' else wd['p'] if f=='p' else 0
        ut=data['u']*tx+data['v']*ty;un=data['u']*nx+data['v']*ny
        uv=data['uu']*tx*nx+data['vv']*ty*ny+data['uv']*(tx*ny+ty*nx)
        good=np.isfinite(ut)
        if not np.all(good):
            # Stop at the first gap; never bridge a mesh hole silently.
            bad=np.flatnonzero(~good);end=bad[0]
        else:end=len(ut)
        if end<10:raise ValueError(f'normal ray x={target} has too few valid points')
        n=distances[:end];ut=ut[:end];rho=data['rho'][:end];ue=ut[-1];rhoe=rho[-1]
        u_tau=np.sqrt(abs(wd['tau_t'])/wd['rho_wall'])
        yp=n*u_tau*wd['rho_wall']/wd['mu_wall_over_Re']
        vd=np.r_[0,np.cumsum(.5*(np.sqrt(rho[1:]/wd['rho_wall'])+np.sqrt(rho[:-1]/wd['rho_wall']))*np.diff(ut))]/max(u_tau,1e-30)
        delta=[]
        for threshold in [.99,.995]:
            found=np.flatnonzero(ut>=threshold*ue)
            j=int(found[0]) if len(found) else len(ut)-1
            value=n[j] if j==0 else np.interp(threshold*ue,ut[j-1:j+1],n[j-1:j+1]);delta.append(float(value))
        flux=rho*ut/(rhoe*ue)
        displacement=float(np.trapezoid(1-flux,n));momentum=float(np.trapezoid(flux*(1-ut/ue),n))
        boundary.append([target,n[-1],ue,delta[0],delta[1],displacement,momentum,displacement/momentum if momentum else np.nan,
                         u_tau,wd['ds']*u_tau*wd['rho_wall']/wd['mu_wall_over_Re'],wd['yplus_mean'],
                         ((.2 if case=='sd7003' else 6)/span_cells*u_tau*wd['rho_wall']/wd['mu_wall_over_Re']) if span_cells else np.nan])
        for j in range(end):out.append([target,n[j],x[j],y[j],ut[j],un[j],yp[j],vd[j],uv[j]]+[data[f][j] for f in fields])
    csvout(destination/'profiles.csv',['wall_x','n','x','y','u_t','u_n','yplus','u_vd_plus','reynolds_tn']+fields,out)
    csvout(destination/'boundary_layer.csv',['wall_x','edge_distance_used','edge_velocity','delta99','delta995','delta_star','theta','H','u_tau_from_mean_tau','dxplus','first_cell_yplus_mean','dzplus'],boundary)
    return dict(profile_stations=targets,edge_distance_requested=edge,
                profile_definition='linear interpolation on explicit structured triangles along a straight wall normal; edge quantities at last valid requested distance',
                boundary_layer_integrals='compressible displacement and momentum thickness integrated to edge_distance_used')

def analyze(prefix,case,destination,start,span_cells=None,reference=None):
    destination.mkdir(parents=True,exist_ok=True)
    name=lambda suffix:Path(str(prefix)+'_benchmark_'+suffix+'.csv')
    mean,wall=read(name('mean_xy')),read(name('wall'))
    if not len(mean) or not len(wall):raise ValueError('no accumulated statistics yet')
    meta=dict(case=case,weight=float(mean['weight'][0]),samples=int(mean['samples'][0]),
              status='diagnostics only; convergence must be checked from independent time windows',
              averaging='accepted physical dt and uniform-span integrals; Reynolds and Favre covariances formed before interpolation')
    wall=np.sort(wall,order='x')
    select=(wall['ny']>0) & (wall['x']>.05) & (wall['x']<.98) if case=='sd7003' else (wall['x']>75)&(wall['x']<140)
    roots=crossings(wall['x'][select],wall['Cf_t'][select]);meta['mean_Cf_zero_crossings']=roots
    meta.update(profiles(mean,wall,case,destination,span_cells));meta['span_cells']=span_cells
    if case=='ramp':
        # q is normalized by rho_inf U_inf^3, positive into the wall.
        cp=1/(.4*2.25**2);Taw=1+np.sqrt(.72)*.2*2.25**2
        csvout(destination/'wall_heat.csv',['x','Cq_kinetic','St_free_recovery','q_into_wall','q_rms'],
               zip(wall['x'],2*wall['q_into_wall'],wall['q_into_wall']/(cp*(Taw-wall['Twall'])),wall['q_into_wall'],wall['q_rms']))
        meta['heat_definition']='Cq=q/(0.5 rho_inf U_inf^3); St=q/[rho_inf U_inf cp (Taw_inf-Tw)], Taw_inf/Tinf=1+sqrt(Pr)*(gamma-1)*M^2/2; use q itself near zero denominator'
    # Transition diagnostic: sustained threshold over 3 adjacent surface columns.
    if case=='sd7003':
        rows=[]
        for w in wall[select]:
            m=mean[(mean['zone']==w['zone'])&(mean['i']==w['i'])]
            n=(m['x']-w['x'])*w['nx']+(m['y']-w['y'])*w['ny']
            values=-m['uv'][(n>=0)&(n<=.1)]
            rows.append([w['x'],float(np.max(values)) if len(values) else np.nan])
        csvout(destination/'transition_indicator.csv',['x','maximum_minus_uv_within_n_0p1'],rows)
        a=np.asarray(rows);hits=a[:,1]>.001;idx=np.flatnonzero(hits[:-2]&hits[1:-1]&hits[2:])
        meta['transition_x_threshold_0p001_3columns']=float(a[idx[0],0]) if len(idx) else None
        meta['transition_note']='diagnostic threshold, not a universal definition; confirm sustained rise and sensitivity to threshold'
    hist=read(name('history'));hist=hist[hist['time']>=start]
    probes=read(name('probes'));probes=probes[probes['time']>=start]
    spectral=[]
    for probe in np.unique(probes['probe']):
        p=probes[probes['probe']==probe]
        for variable in ['u','v','p']:
            result=spectrum(p['time'],p[variable])
            if result:
                freq,power,info=result;spectral.append(dict(probe=int(probe),quantity=variable,**info))
                csvout(destination/f'probe{int(probe)}_{variable}_psd.csv',['f_Lref_Uref','PSD','f_PSD'],zip(freq,power,freq*power))
    for key in ['Cd_or_Cx','Cl_or_Cy','Cmz','xs','xr','wall_pressure_gradient_x','q_into_wall_max']:
        if key not in hist.dtype.names or not np.all(np.isfinite(hist[key])):continue
        result=spectrum(hist['time'],hist[key])
        if result:
            freq,power,info=result;spectral.append(dict(quantity=key,**info))
            csvout(destination/(key+'_psd.csv'),['f_Lref_Uref','PSD','f_PSD'],zip(freq,power,freq*power))
    meta['spectra']=spectral;meta['spectrum_start']=start
    if len(hist)>=128 and all(k in hist.dtype.names for k in ['xs','xr']):
        length=hist['xr']-hist['xs'];finite=length[np.isfinite(length)]
        meta['mean_instantaneous_bubble_length_for_St']=float(np.mean(finite)) if len(finite) else None
        meta['St_conversion']='multiply f_Lref_Uref by a separately converged mean separation length / Lref; do not compare f delta0/Uinf directly with f Lsep/Uinf'
    # Pressure coherence between simultaneous fixed probes, using the same resampling.
    ids=np.unique(probes['probe']);coherences=[]
    for ia in range(len(ids)):
        for ib in range(ia+1,len(ids)):
            a=probes[probes['probe']==ids[ia]];b=probes[probes['probe']==ids[ib]]
            if len(a)!=len(b) or not np.array_equal(a['time'],b['time']):raise ValueError('probe histories are not simultaneous')
            ua,ub=uniform_history(a['time'],a['p']),uniform_history(b['time'],b['p'])
            if ua is None or ub is None or np.std(ua[1])==0 or np.std(ub[1])==0:continue
            f,c=coherence(ua[1],ub[1],fs=1/ua[2],nperseg=min(2048,len(ua[0])//4),window='hann')
            name=f'pressure_coherence_{int(ids[ia])}_{int(ids[ib])}.csv';csvout(destination/name,['f_Lref_Uref','magnitude_squared_coherence'],zip(f,c));coherences.append(name)
    meta['pressure_coherence']=coherences
    if not spectral:meta['spectra_status']='not computed: fewer than 128 samples (a smoke test cannot establish spectra)'
    # Eight disjoint time blocks; coarse convergence evidence, not independence guarantee.
    if len(hist)>=128:
        edges=np.linspace(hist['time'][0],hist['time'][-1],9);blocks=[]
        for a,b in zip(edges[:-1],edges[1:]):
            t=np.r_[a,hist['time'][(hist['time']>a)&(hist['time']<b)],b]
            blocks.append([a,b]+[float(np.trapezoid(np.interp(t,hist['time'],hist[key]),t)/(b-a)) for key in ['Cd_or_Cx','Cl_or_Cy','Cmz']])
        csvout(destination/'time_blocks.csv',['start','end','Cd_or_Cx','Cl_or_Cy','Cmz'],blocks)
        meta['block_standard_errors']=list(np.std(np.asarray(blocks)[:,2:],axis=0,ddof=1)/np.sqrt(8))
        meta['block_warning']='block length must exceed integral correlation time; these estimates alone do not prove convergence'
    fig,axs=plt.subplots(2,1,figsize=(8,7),sharex=True)
    if case=='sd7003':
        for label,mask in [('upper',wall['ny']>0),('lower',wall['ny']<0)]:
            axs[0].plot(wall['x'][mask],wall['Cp'][mask],label=label)
        axs[0].invert_yaxis();axs[0].legend();mask=wall['ny']>0
    else:mask=np.ones(len(wall),dtype=bool)
    if case!='sd7003':axs[0].plot(wall['x'],wall['Cp']);axs[0].set_xlim(70,140)
    axs[1].plot(wall['x'][mask],wall['Cf_t' if case=='sd7003' else 'Cfx'][mask]);axs[1].axhline(0,color='gray',lw=.7)
    axs[0].set_ylabel('Cp');axs[1].set_ylabel('Cf_t' if case=='sd7003' else 'Cfx');axs[1].set_xlabel('x/c' if case=='sd7003' else 'x/delta0')
    fig.suptitle('Accumulated statistics (convergence not implied)');fig.tight_layout();fig.savefig(destination/'wall_statistics.png',dpi=150);plt.close(fig)
    if reference:
        compare_reference(wall,case,destination,reference)
    (destination/'analysis.json').write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n')
    return meta

def compare_reference(wall,case,destination,reference):
    """Compare digitized points without confusing interpolated values with raw data."""
    fig,axs=plt.subplots(1,2,figsize=(10,4));results=[]
    for ax,stem,column in zip(axs,['Cp','Cf' if case=='sd7003' else 'Cfx'],['Cp','Cf_t' if case=='sd7003' else 'Cfx']):
        r=read(reference/(stem+'_digitized.csv'));xname=r.dtype.names[0];yname=r.dtype.names[1]
        branches=[('upper',wall['ny']>0),('lower',wall['ny']<0)] if case=='sd7003' and stem=='Cp' else [('surface',wall['ny']>0)]
        candidates=[]
        for label,mask in branches:
            w=np.sort(wall[mask],order='x');ax.plot(w['x'],w[column],label=label)
            candidates.append(np.interp(r[xname],w['x'],w[column],left=np.nan,right=np.nan))
        pred=candidates[0] if len(candidates)==1 else np.where(np.abs(candidates[0]-r[yname])<np.abs(candidates[1]-r[yname]),candidates[0],candidates[1])
        ax.errorbar(r[xname],r[yname],yerr=(r['y_upper']-r['y_lower'])/2,fmt='ks',ms=3,label='digitized reference');ax.set_xlabel(xname);ax.set_ylabel(stem);ax.legend()
        csvout(destination/(stem+'_reference_errors.csv'),[xname,'reference','computed_interpolated','difference','reference_reading_halfwidth'],zip(r[xname],r[yname],pred,pred-r[yname],(r['y_upper']-r['y_lower'])/2))
        good=np.isfinite(pred);results.append(dict(quantity=stem,points=int(np.sum(good)),rmse=float(np.sqrt(np.mean((pred[good]-r[yname][good])**2))),
            note='digitized reference, not raw exact data; SD Cp uses nearest branch and must be visually checked near LE/TE'))
    if case=='sd7003':axs[0].invert_yaxis()
    fig.tight_layout();fig.savefig(destination/'reference_comparison.png',dpi=160);plt.close(fig)
    (destination/'reference_errors.json').write_text(json.dumps(results,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',choices=['sd7003','ramp'],required=True);p.add_argument('--prefix',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--start',type=float,default=0)
    p.add_argument('--span-cells',type=int);p.add_argument('--reference',type=Path)
    a=p.parse_args()
    if a.span_cells is not None and a.span_cells<1:p.error('--span-cells must be positive')
    print(json.dumps(analyze(a.prefix,a.case,a.output,a.start,a.span_cells,a.reference),indent=2))
