"""Reproducible conforming multi-block extrusions, in c or delta0 units.
Run with numpy/scipy; then wcns_generate_chapter5_cgns plane.bin mesh.cgns.
The original source coordinate file must be present; no silent network access.
"""
from pathlib import Path
import argparse, struct, json
import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator, CubicHermiteSpline
from scipy.ndimage import gaussian_filter1d

ROOT=Path(__file__).resolve().parents[1]

def geometric(n,total,first):
    lo,hi=1.,2.
    for _ in range(100):
        r=(lo+hi)/2
        if first*np.expm1(n*np.log(r))/np.expm1(np.log(r))>total: hi=r
        else: lo=r
    d=first*r**np.arange(n);return np.r_[0,np.cumsum(d)]*total/sum(d)

def sd7003(nx,ny,coordinates=None):
    raw=np.loadtxt(coordinates or ROOT/'cases/manual/case09_sd7003/reference/sd7003-uiuc.dat',skiprows=1)
    le=np.argmin(raw[:,0]);up=raw[:le+1][::-1];low=raw[le:]
    u=PchipInterpolator(up[:,0],up[:,1]);l=PchipInterpolator(low[:,0],low[:,1]);r=.0004
    # Keep the chordwise TE position. Blend only the final 2% of the sharp
    # UIUC section into a radius-.0004 semicircle, matching position/slope.
    cx,cy=1-r,0.
    cut=.98
    us=CubicHermiteSpline([cut,cx],[u(cut),r],[u.derivative()(cut),0])
    ls=CubicHermiteSpline([cut,cx],[l(cut),-r],[l.derivative()(cut),0])
    xx=np.linspace(cut,cx,80)
    up_points=np.vstack([up[up[:,0]<cut],np.column_stack([xx,us(xx)])])
    lo_points=np.vstack([low[low[:,0]<cut],np.column_stack([xx,ls(xx)])])
    xu=xl=cx
    au,al=np.pi/2,-np.pi/2
    # Clockwise curve starting at rightmost TE: lower arc, lower, LE, upper, upper arc.
    lower=lo_points[lo_points[:,0]<xl][::-1]
    upper=up_points[(up_points[:,0]<xu) & (up_points[:,0]>up_points[0,0])]
    aa=np.r_[np.linspace(0,al,30,endpoint=False)]
    bb=np.linspace(au,0,30)
    curve=np.vstack([np.column_stack([cx+r*np.cos(aa),cy+r*np.sin(aa)]),[xl,-r],lower,upper,[xu,r],np.column_stack([cx+r*np.cos(bb),cy+r*np.sin(bb)])])
    keep=np.r_[True,np.linalg.norm(np.diff(curve,axis=0),axis=1)>1e-10];curve=curve[keep]
    curve[-1]=curve[0]
    s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(curve,axis=0),axis=1))]
    spline=CubicSpline(s,curve,bc_type='periodic')
    ss=np.linspace(0,s[-1],20001);q=spline(ss);d=spline(ss,1);dd=spline(ss,2)
    curvature=np.abs(d[:,0]*dd[:,1]-d[:,1]*dd[:,0])/np.linalg.norm(d,axis=1)**3
    dist_te=np.minimum(ss,s[-1]-ss)
    sle=ss[np.argmin(q[:,0])]
    density=1+100*np.exp(-(dist_te/.007)**2)+10*np.exp(-((ss-sle)/.02)**2)
    cum=np.r_[0,np.cumsum(.5*(density[1:]+density[:-1])*np.diff(ss))]
    nodes=np.interp(np.linspace(0,cum[-1],nx+1),cum,ss)
    wall=spline(nodes)
    nn=np.column_stack([-d[:,1],d[:,0]])/np.linalg.norm(d,axis=1)[:,None]
    nn=gaussian_filter1d(nn,.0007/(ss[1]-ss[0]),axis=0,mode='wrap')
    normal=np.column_stack([np.interp(nodes,ss,nn[:,a]) for a in range(2)])
    normal/=np.linalg.norm(normal,axis=1)[:,None]
    # Smoothly open normal rays into a star-shaped radial map; after that use
    # geometric radial stretching out to a circular farfield of radius 30c.
    phi=-2*np.pi*nodes/s[-1]
    theta=np.unwrap(np.arctan2(4*np.sin(phi),np.cos(phi)))
    radial=np.column_stack([np.cos(theta),np.sin(theta)])
    heights=geometric(ny,30,3e-5 if ny>=56 else .001)
    xy=[]
    for h in heights:
        a=.002*(-np.expm1(-h/.002))
        xy.append(wall+a*normal+(h-a)*radial)
    xy=np.array(xy);xy[:,-1]=xy[:,0]
    return xy,dict(trailing_edge_radius_input=r,te_circle_center=[cx,cy],rounding_blend_start=cut,
                   chord_realized=float(wall[:,0].max()-wall[:,0].min()),first_height=float(heights[1]),
                   farfield='radial offset 30c from wall (approximately 30c radius)',mapping='normal displacement 0.002*(1-exp(-h/0.002)); remaining displacement radial')

def ramp(nx,ny):
    # Exact corner at source-zone boundary: 8 of 12 streamwise blocks upstream.
    # PCHIP avoids discontinuous spacing at the trip interval and buffer joins.
    s=np.array([0,48,64,128,512,704,768])/768
    x=PchipInterpolator(s,[0,2.33,2.67,10,100,140,190])(np.linspace(0,1,nx+1))
    angle=np.deg2rad(24);wall=np.maximum(x-100,0)*np.tan(angle)
    a=np.clip((x-70)/60,0,1);normal_angle=angle*a*a*(3-2*a)
    heights=geometric(ny,40,.004 if ny>60 else .002)
    xy=[]
    for h in heights:
        # Height buffer downstream; near-wall spacing unchanged to first order.
        distance=h+np.maximum(x-100,0)*np.tan(np.deg2rad(16))*(h/40)**3
        xy.append(np.column_stack([x-distance*np.sin(normal_angle),wall+distance*np.cos(normal_angle)]))
    return np.array(xy),dict(corner_x=100,angle_degrees=24,wall_x_max=190,upstream_height=40,
        first_height=float(heights[1]),trip_x=[2.33,2.67],trip_y=[0,.01],
        mapping='near-wall normal angle blended 0 to 24 deg over x=70..130; outer downstream height buffer')

def write(kind,small,case_dir=None):
    nb=8 if kind=='sd7003' else 12
    nx,ny,nz=((56,96,40) if kind=='sd7003' else (64,160,72)) if not small else ((56,96,8) if kind=='sd7003' else (8,24,8))
    span=.2 if kind=='sd7003' else 6
    case=case_dir or ROOT/'cases/manual'/('case09_sd7003' if kind=='sd7003' else 'case10_compression_ramp')
    xy,info=sd7003(nb*nx,ny,case/'reference/sd7003-uiuc.dat') if kind=='sd7003' else ramp(nb*nx,ny)
    # All four bilinear corner determinants, not only signed average area.
    di=np.diff(xy,axis=1);dj=np.diff(xy,axis=0)
    dets=[]
    for a in range(2):
        for b in range(2):
            p=di[b:ny+b];q=dj[:,a:nx*nb+a]
            dets.append(p[:,:,0]*q[:,:,1]-p[:,:,1]*q[:,:,0])
    mindet=float(np.min(dets))
    if mindet<=0:
        idx=np.unravel_index(np.argmin(dets),np.shape(dets));raise RuntimeError(f'{kind}: inverted cell {idx}: {mindet}')
    grid=case/'grids';grid.mkdir(parents=True,exist_ok=True)
    suffix='smoke' if small else 'coarse';path=grid/(kind+'_'+suffix+'.bin')
    with path.open('wb') as f:
        f.write(struct.pack('<4s5Id',b'XY25',0 if kind=='sd7003' else 1,nb,nx,ny,nz,span))
        for b in range(nb):
            local=xy[:,b*nx:(b+1)*nx+1,:]
            for component in [0,1]:f.write(np.array(local[:,:,component],dtype='<f8').tobytes())
    info.update(zones=nb,cells_per_zone=[nx,ny,nz],cells=nb*nx*ny*nz,span=span,min_corner_volume=mindet*span/nz)
    (grid/(kind+'_'+suffix+'.json')).write_text(json.dumps(info,indent=2)+'\n')
    np.savez_compressed(grid/(kind+'_'+suffix+'_xy.npz'),xy=xy)
    print(path,info)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case',choices=['sd7003','ramp']);p.add_argument('--smoke',action='store_true')
    p.add_argument('--case-dir',type=Path,help='case root with reference/ and grids/ (default: development repository)')
    args=p.parse_args();write(args.case,args.smoke,args.case_dir)
