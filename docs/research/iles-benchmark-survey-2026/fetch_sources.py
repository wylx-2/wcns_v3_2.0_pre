"""Download small public research sources, preserving provenance and hashes."""
from pathlib import Path
import concurrent.futures, hashlib, json, urllib.request
ROOT = Path(__file__).resolve().parent / 'sources'
SOURCES = {
 'tgv_spec.pdf': 'https://cfd.ku.edu/hiocfd/case_c3.5.pdf',
 'tgv_data.tgz': 'https://cfd.ku.edu/hiocfd/data.tgz',
 'LM_Channel_0180_mean_prof.dat': 'https://turbulence.oden.utexas.edu/channel2015/data/LM_Channel_0180_mean_prof.dat',
 'LM_Channel_0180_vel_fluc_prof.dat': 'https://turbulence.oden.utexas.edu/channel2015/data/LM_Channel_0180_vel_fluc_prof.dat',
 'cbc1971.pdf': 'https://courses.washington.edu/mengr544/handouts/comtebellot-corrsin-jfm-71.pdf',
 'bfs_cf.dat': 'https://tmbwg.github.io/turbmodels/Backstep_validation/cf.exp.dat',
 'bfs_cp.dat': 'https://tmbwg.github.io/turbmodels/Backstep_validation/cp.expnew.dat',
 'hump_cf.dat': 'https://tmbwg.github.io/turbmodels/Nasahump_validation/noflow_cf.exp.dat',
 'hump_cp.dat': 'https://tmbwg.github.io/turbmodels/Nasahump_validation/noflow_cp.exp.dat',
 'jhtdb_isotropic.pdf': 'https://turbulence.pha.jhu.edu/docs/README-isotropic.pdf',
 'kth_vel1410.prof': 'https://www.mech.kth.se/~pschlatt/DATA/vel_1410_dns.prof',
 'sbli2023.pdf': 'https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C6721938F7F8F64F486C2047C6E985A1/S0022112022010382a.pdf/unsteadiness_characterisation_of_shock_waveturbulent_boundarylayer_interaction_at_moderate_reynolds_number.pdf',
 'jet2018.pdf': 'https://atowne.com/wp-content/uploads/2018/07/BresEtAl_2018_JFM_jetLES.pdf',
 'apg2017.pdf': 'https://torroja.dmt.upm.es/pubs/2017/KitsiosEtal_JFM17.pdf',
 'shockturb_fullStats.tar.gz': 'https://drive.google.com/uc?export=download&id=1-rf97iz_rISUhT4jmqleWT4HnvEJtVle',
}
def fetch(item):
 name, url = item
 out = {'file': name, 'url': url, 'retrieved': '2026-10-03'}
 try:
  if (ROOT/name).exists():
   data=(ROOT/name).read_bytes()
   return dict(out,bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),status='downloaded')
  req = urllib.request.Request(url, headers={'User-Agent':'Mozilla/5.0'})
  with urllib.request.urlopen(req, timeout=40) as r: data = r.read(40_000_001)
  if len(data)>40_000_000: raise ValueError('Unexpectedly large response')
  if name.endswith('.pdf') and not data.startswith(b'%PDF'): raise ValueError('Not PDF')
  if name.endswith(('.tgz','.gz')) and not data.startswith(b'\x1f\x8b'): raise ValueError('Not gzip')
  if name.endswith(('.dat','.prof')) and b'<html' in data.lower(): raise ValueError('HTML instead of data')
  (ROOT/name).write_bytes(data)
  out.update(bytes=len(data), sha256=hashlib.sha256(data).hexdigest(), status='downloaded')
 except Exception as e: out.update(status='unavailable', error=str(e))
 return out
if __name__ == '__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
  manifest = list(pool.map(fetch, SOURCES.items()))
 (ROOT/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
 for row in manifest: print(row['file'],row['status'],row.get('bytes',row.get('error')))
