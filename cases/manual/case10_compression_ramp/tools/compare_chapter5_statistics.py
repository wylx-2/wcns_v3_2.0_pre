"""Compare numeric accepted-time moment CSVs from equivalent short computations."""
import argparse,json
from pathlib import Path
import numpy as np

def compare(a,b):
    result={}
    for suffix in ['mean_xy','wall']:
        x=np.genfromtxt(str(a)+'_benchmark_'+suffix+'.csv',delimiter=',',names=True)
        y=np.genfromtxt(str(b)+'_benchmark_'+suffix+'.csv',delimiter=',',names=True)
        assert x.shape==y.shape and x.dtype.names==y.dtype.names
        errors={}
        for key in x.dtype.names:
            assert np.all(np.isfinite(x[key])) and np.all(np.isfinite(y[key])),key
            error=float(np.max(abs(x[key]-y[key]))/max(1.,np.max(abs(x[key]))))
            assert error<2e-9,(suffix,key,error)
            errors[key]=error
        result[suffix]=dict(max_scaled_absolute_error=max(errors.values()),columns=errors,
                            rows=len(x),weight=float(x['weight'][0]))
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('a',type=Path);p.add_argument('b',type=Path);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();r=compare(args.a,args.b);args.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v['max_scaled_absolute_error'] for k,v in r.items()}))
