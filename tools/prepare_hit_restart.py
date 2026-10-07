"""Write a restart config beside its base, preserving strict HIT spectrum identity."""
import argparse
from pathlib import Path

def replace(text,key,value):
    rows=text.splitlines(); found=False
    for i,row in enumerate(rows):
        if row.split('=',1)[0].strip()==key: rows[i]=f'{key} = {value}'; found=True
    if not found: rows.append(f'{key} = {value}')
    return '\n'.join(rows)+'\n'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('base',type=Path);p.add_argument('checkpoint',type=Path)
    p.add_argument('--name',default='resume.wcns');p.add_argument('--output',required=True)
    a=p.parse_args();base=a.base.resolve(); checkpoint=a.checkpoint.resolve(strict=True)
    if Path(a.name).name!=a.name or not a.name.endswith('.wcns'): p.error('--name must be a .wcns filename without directories')
    dest=base.parent/a.name
    text=base.read_text(encoding='utf-8')
    for k,v in {'restart.path':checkpoint.as_posix(),'restart.mode':'strict','output.directory':a.output,'output.allow_existing':'false'}.items():
        text=replace(text,k,v)
    with dest.open('x',encoding='utf-8') as f:f.write(text)
    print(dest)

if __name__=='__main__':main()
