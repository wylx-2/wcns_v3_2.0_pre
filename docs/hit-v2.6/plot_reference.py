from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent/'figures'
fig,axes=plt.subplots(1,2,figsize=(10,3.6),layout='constrained')
for station in [42,98,171]:
 a=np.loadtxt(ROOT/f'cases/manual/case11_hit_decay/reference/cbc_spectrum_{station}_nondimensional.dat')
 axes[0].loglog(a[:,0],a[:,1],'-o',ms=3,label=f'x/M = {station}')
axes[0].set(xlabel='Nondimensional angular wavenumber k*',ylabel='Radial spectrum E*(k*)',title='CBC experimental spectra')
a=np.loadtxt(ROOT/'cases/manual/case12_hit_forced/reference/jhtdb_spectrum.dat')
axes[1].loglog(a[:,0],a[:,1],lw=1.5,label='Published DNS spectrum')
axes[1].set(xlabel='Integer-shell wavenumber k',ylabel='Radial spectrum E(k)',title='JHTDB forced isotropic DNS')
for ax,dk in zip(axes,[2*np.pi,1]):
 for n in [32,64,128]:
  ax.axvline(dk*n/3,color='0.6',ls=':',lw=.8)
  ax.text(dk*n/3,1e-5,f'N={n}',rotation=90,fontsize=8,ha='right')
 ax.legend(fontsize=8);ax.grid(alpha=.2,which='both')
fig.savefig(OUT/'reference_spectra.png',dpi=200)
