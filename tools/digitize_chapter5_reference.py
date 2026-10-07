"""Extract visible black reference squares, never the red thesis result.
Input: directory containing the original chapter-5 figures supplied by the user.
Output: calibrated approximate CSVs, source-image hashes, and annotated QA images.
Only clearly detected markers are retained; omitted/occluded markers are not invented.
Requires numpy, scipy, Pillow. Pixel errors are conservative digitization estimates,
not the uncertainty of the underlying reference simulation or experiment.
"""
from pathlib import Path
import argparse,csv,hashlib,json,shutil
import numpy as np
from PIL import Image,ImageDraw
from scipy.ndimage import uniform_filter,maximum_filter

ROOT=Path(__file__).resolve().parents[1]
# Axes were independently read from the labelled major ticks (native pixels).
# Each specification: file, case, output stem, square width, ROI, excluded legend,
# x pixel/value pair, y pixel/value pair, logarithmic x, column names.
SPECS=[
 ('图5.24_1.jpeg','sd7003','Cp',26,(210,85,1380,1260),(650,70,1375,290),(226,1368,0,1),(747,999,0,.5),False,('x_over_c','Cp')),
 ('图5.24_2.jpeg','sd7003','Cf',26,(244,80,1350,1240),None,(234,1352,0,1),(1014,724,0,.01),False,('x_over_c','Cf')),
 ('图5.25_1.jpg','sd7003','velocity_x0p1',45,(400,224,2300,2170),(400,200,1700,650),(529.03125,1751.0625,0,1),(2176.453125,221.203125,0,.1),False,('u_t_over_Uinf','n_over_c')),
 ('图5.25_2.jpg','sd7003','velocity_x0p5',45,(400,224,2300,2170),None,(529.03125,1751.0625,0,1),(2176.453125,221.203125,0,.1),False,('u_t_over_Uinf','n_over_c')),
 ('图5.25_3.jpg','sd7003','velocity_x0p95',45,(400,224,2300,2170),None,(529.03125,1751.0625,0,1),(2176.453125,221.203125,0,.1),False,('u_t_over_Uinf','n_over_c')),
 ('图5.33_1.jpeg','ramp','van_driest_x80',31,(405,140,1280,1260),(405,140,930,300),(188,1296,-1,3),(1274,315,0,20),True,('yplus','u_vd_plus')),
 ('图5.34_1.jpeg','ramp','Cp',31,(260,140,1345,1230),(260,130,930,285),(444,1159,80,120),(1201,321,0,.6),False,('x_over_delta0','Cp')),
 ('图5.34_2.jpeg','ramp','Cfx',31,(260,250,1345,970),None,(444,1159,80,120),(785.5,198.5,0,.003),False,('x_over_delta0','Cfx')),
]

def extract(im,size,roi,exclude):
 a=np.max(np.asarray(im),axis=2)<110;inner=25 if size==45 else size-8
 outer=uniform_filter(a.astype(float),size=size,mode='constant')*size**2
 inside=uniform_filter(a.astype(float),size=inner,mode='constant')*inner**2
 score=(outer-inside)/(size**2-inner**2)-.5*inside/inner**2
 peaks=(score==maximum_filter(score,size=max(9,int(size*.6)))) & (score>(.43 if size==45 else .58))
 yy,xx=np.indices(a.shape);x0,y0,x1,y1=roi
 peaks&=(xx>x0)&(xx<x1)&(yy>y0)&(yy<y1)
 if exclude:
  x0,y0,x1,y1=exclude;peaks&=~((xx>x0)&(xx<x1)&(yy>y0)&(yy<y1))
 y,x=np.where(peaks);found=[]
 for j in np.argsort(score[y,x])[::-1]:
  p=np.array([x[j],y[j]])
  if all(np.linalg.norm(p-q[:2])>size*.8 for q in found):found.append(np.r_[p,score[y[j],x[j]]])
 return sorted(found,key=lambda q:(q[0],q[1]))

def main(source):
 summaries=[]
 for name,case,stem,size,roi,ex,xc,yc,logx,columns in SPECS:
  dest=ROOT/'cases/manual'/('case09_sd7003' if case=='sd7003' else 'case10_compression_ramp')/'reference'
  dest.mkdir(parents=True,exist_ok=True);im=Image.open(source/name).convert('RGB');pts=extract(im,size,roi,ex)
  # Visual QA identified axis/tick intersections masquerading as partial squares.
  # Retain genuine near-wall markers; do not fill missing ones with interpolated data.
  pts=[p for p in pts if not (size==45 and p[1]>2145 and p[2]<.65)
       and not (stem=='van_driest_x80' and p[1]>1248)]
  source_name='source_'+stem+Path(name).suffix;shutil.copyfile(source/name,dest/source_name)
  def val(p,c):return c[2]+(p-c[0])*(c[3]-c[2])/(c[1]-c[0])
  # +/- 4 native pixels includes axis calibration and imperfect marker centre.
  eps=4.;rows=[];draw=ImageDraw.Draw(im)
  for number,(x,y,score) in enumerate(pts):
   vx,vy=val(x,xc),val(y,yc);dx=eps*abs((xc[3]-xc[2])/(xc[1]-xc[0]));dy=eps*abs((yc[3]-yc[2])/(yc[1]-yc[0]))
   xmin,xmax=(10**(vx-dx),10**(vx+dx)) if logx else (vx-dx,vx+dx)
   if logx:vx=10**vx
   rows.append([vx,vy,xmin,xmax,vy-dy,vy+dy,int(x),int(y),score])
   draw.ellipse((x-6,y-6,x+6,y+6),outline='lime',width=3);draw.text((x+7,y-8),str(number),fill='blue')
  with (dest/(stem+'_digitized.csv')).open('w',newline='') as f:
   w=csv.writer(f);w.writerow([*columns,'x_lower','x_upper','y_lower','y_upper','pixel_x','pixel_y','marker_score']);w.writerows(rows)
  im.save(dest/(stem+'_digitization_qa.png'))
  meta=dict(status='approximate digitization of visible black reference squares; not original numerical data',
    document='test_4.docx, chapter 5',figure=name,reference='Galbraith & Visbal 2010' if case=='sd7003' else 'Porter & Poggie 2019 as reproduced in the supplied thesis',
    source_file=source_name,source_sha256=hashlib.sha256((source/name).read_bytes()).hexdigest(),
    native_pixels=im.size,x_calibration=list(xc),y_calibration=list(yc),x_log10=logx,
    pixel_error_estimate=eps,uncertainty_note='coordinate/marker reading estimate only, not a rigorous bound on source-data uncertainty',
    marker_count=len(rows),selection='hollow-square black border score > 0.58 (profile figures 0.43), nonmaximum suppression; visual-QA axis/tick exclusions; incomplete visible-marker subset; no curve interpolation')
  (dest/(stem+'_provenance.json')).write_text(json.dumps(meta,ensure_ascii=False,indent=2)+'\n',encoding='utf8');summaries.append([case,stem,len(rows)])
 print(json.dumps(summaries))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('source',type=Path);main(p.parse_args().source)
