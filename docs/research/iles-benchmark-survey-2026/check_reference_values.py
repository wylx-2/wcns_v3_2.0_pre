"""Recompute the report's derived numeric anchors without running WCNS."""
from pathlib import Path
import json, hashlib
import numpy as np
R=Path(__file__).resolve().parent
S=R/'sources'
tgv=np.loadtxt(S/'spectral_Re1600_512.gdiag')
peak=tgv[np.argmax(tgv[:,2])]
def numeric(name):
 rows=[]
 for line in (S/name).read_text().splitlines():
  try: a=[float(x) for x in line.split()]
  except ValueError:continue
  if len(a)>=2:rows.append(a)
 return rows
def zero(a,b):return a[0]-a[1]*(b[0]-a[0])/(b[1]-a[1])
h=numeric('hump_cf.dat');b=numeric('bfs_cf.dat')
h1=next(r for r in h if r[0]==1.09);h2=next(r for r in h if r[0]==1.15)
b1=next(r for r in b if r[0]==5.882);b2=next(r for r in b if r[0]==7.09)
summary={
 'TGV_peak_sample': dict(zip(['t','K','epsilon','enstrophy'],peak.tolist())),
 'TGV_time_interval':[float(tgv[0,0]),float(tgv[-1,0])],
 'LM180_Cf_from_header':2*.0637309**2,
 'hump_reattachment_linear_interpolation':zero(h1,h2),
 'bfs_reattachment_linear_interpolation_not_published_value':zero(b1,b2),
 'CBC_table4_Re_lambda':[71.6,65.3,60.7],
 'CBC_epsilon_m2_s3':[.474,.0633,.0174],
 'note':'Derived values are not new DNS or CFD results; original published precision and uncertainty apply.'}
assert abs(peak[0]-8.97)<1e-10
assert tgv[-1,0]==19.99
assert abs(summary['LM180_Cf_from_header']-.00812325522962)<1e-14
(R/'reference_anchors.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=True,indent=2))
