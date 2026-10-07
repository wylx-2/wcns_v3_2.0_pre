"""Plot the saved planar structured meshes without loading their 3D CGNS copies."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
for folder,stem,panels in [
 ('case09_sd7003','sd7003',[((-1,2),(-1,1),'Near field'),((-.02,.12),(-.04,.08),'Leading edge'),((.98,1.004),(-.012,.012),'Rounded trailing edge')]),
 ('case10_compression_ramp','ramp',[((0,190),(0,100),'Full domain'),((2.25,2.75),(0,.025),'Trip strip'),((75,125),(0,18),'Compression corner')])]:
    grid=ROOT/'cases/manual'/folder/'grids';xy=np.load(grid/(stem+'_coarse_xy.npz'))['xy']
    fig,axs=plt.subplots(1,3,figsize=(14,4.3));nb=8 if stem=='sd7003' else 12;stride=(xy.shape[1]-1)//nb
    for ax,(xr,yr,title) in zip(axs,panels):
        for j in range(xy.shape[0]):ax.plot(xy[j,:,0],xy[j,:,1],color='.65',lw=.35)
        for i in range(xy.shape[1]):ax.plot(xy[:,i,0],xy[:,i,1],color='.65',lw=.35)
        for i in range(0,xy.shape[1],stride):ax.plot(xy[:,i,0],xy[:,i,1],color='tab:blue',lw=1)
        ax.plot(xy[0,:,0],xy[0,:,1],'k-',lw=1);ax.set_xlim(xr);ax.set_ylim(yr);ax.set_title(title);ax.set_xlabel('x / reference length');ax.set_ylabel('y / reference length')
        if title!='Trip strip':ax.set_aspect('equal',adjustable='box')
    fig.suptitle(stem+' production structured mesh (blue: block interfaces)');fig.tight_layout();fig.savefig(grid/'mesh_xy.png',dpi=160);plt.close(fig)
