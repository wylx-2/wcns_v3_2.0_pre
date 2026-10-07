"""Independent signals test: PSD units, energy and nonuniform sampling/gap guard."""
import importlib.util
from pathlib import Path
import numpy as np
spec=importlib.util.spec_from_file_location('post',Path(__file__).resolve().parents[1]/'tools/analyze_chapter5.py')
post=importlib.util.module_from_spec(spec);spec.loader.exec_module(post)
t=np.arange(0,40,.002);u=3+2*np.sin(2*np.pi*10*t)
f,p,info=post.spectrum(t,u)
assert abs(f[np.argmax(p)]-10)<info['frequency_spacing']
assert abs(np.trapezoid(p,f)-2)<.01 # variance of amplitude-2 sine
tn=t+1e-5*np.sin(t);result=post.spectrum(tn,3+2*np.sin(2*np.pi*10*tn))
assert abs(result[0][np.argmax(result[1])]-10)<result[2]['frequency_spacing']
assert post.spectrum(t[:3],u[:3]) is None
try:post.uniform_history(np.r_[t[:200],t[400:600]],u[:400])
except ValueError:pass
else:raise AssertionError('unreported gap')
assert post.crossings(np.array([0,1,2]),np.array([1,-1,1]))==[.5,1.5]
print('postprocess PSD frequency, variance, nonuniform sampling, gap guard: PASS')
