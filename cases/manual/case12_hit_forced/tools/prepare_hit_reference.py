"""Reproduce HIT input spectra from bundled original numeric tables (Python stdlib)."""
import argparse,csv,hashlib,json,math
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--case',choices=['decay','forced'],required=True)
    p.add_argument('--reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True);manifest=[]
    def save(name,text,source):
        path=a.output/name;path.write_text(text,encoding='utf-8')
        manifest.append({'source':source.name,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'file':name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    if a.case=='decay':
        L=.5588;U=math.sqrt(1.5)*.222
        for station in [42,98,171]:
            source=a.reference/f'cbc_spectrum_{station}_original.csv'
            with source.open() as f:rows=list(csv.reader(f))[1:]
            out='# CBC Table 3; k*=100*k_cm^-1*Lref; E*=E_cm^3/s^2*1e-6/(Uref^2*Lref)\n'
            out+=''.join(f'{float(r[0])*100*L:.17g} {float(r[1])*1e-6/(U*U*L):.17g}\n' for r in rows)
            save(f'cbc_spectrum_{station}_nondimensional.dat',out,source)
    else:
        source=a.reference/'spectrum.txt';rows=[]
        for line in source.read_text().splitlines():
            try:r=[float(v) for v in line.split()]
            except ValueError:continue
            if len(r)==2:rows.append(r)
        if not rows or any(r[0]<=0 or r[1]<0 for r in rows):raise ValueError('invalid JHTDB spectrum')
        save('jhtdb_spectrum.dat','# JHTDB original shell energy; box length 2*pi, dk=1\n'+''.join(f'{k:.17g} {e:.17g}\n' for k,e in rows),source)
    (a.output/'conversion_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))

if __name__=='__main__':main()
