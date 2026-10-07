"""Compare WCNS HIT statistics to bundled CBC/JHTDB data; no third-party Python dependencies.

python tools/analyze_hit.py --case decay --reference reference --output output/run --result comparison
For restart histories pass --history segment1/*_hit_history.csv (repeat as needed).
Numerical errors are reported, never interpreted as automatic physical validation.
"""
import argparse, bisect, csv, json, math, statistics
from pathlib import Path

def csv_rows(path):
    with path.open(encoding='utf-8-sig') as f:
        return [{k: (v if k == 'quantity' else float(v)) for k, v in r.items()} for r in csv.DictReader(f)]

def table(path):
    result = []
    for line in path.read_text().splitlines():
        try: values = [float(x) for x in line.split()]
        except ValueError: continue
        if len(values) >= 2: result.append(values)
    return result

def interpolate(data, x):
    """Same log-log interpolation used for experimental initial spectra."""
    if not data[0][0] <= x <= data[-1][0]: return None
    i = bisect.bisect_left([p[0] for p in data], x)
    if i == 0: return data[0][1]
    a, b = data[i-1], data[i]
    t = math.log(x / a[0]) / math.log(b[0] / a[0])
    return math.exp(math.log(a[1]) * (1-t) + math.log(b[1]) * t) if min(a[1], b[1]) > 0 else a[1]*(1-t)+b[1]*t

def integrate(rows, key, begin, end):
    if end <= begin or rows[0]['time'] > begin or rows[-1]['time'] < end: return None
    value = 0.
    for a, b in zip(rows, rows[1:]):
        lo, hi = max(begin, a['time']), min(end, b['time'])
        if hi <= lo: continue
        slope = (b[key]-a[key])/(b['time']-a['time'])
        value += (hi-lo)*(a[key]+slope*((lo+hi)/2-a['time']))
    return value/(end-begin)

def write_csv(path, rows):
    if not rows: return
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case', choices=['decay','forced'], required=True)
    p.add_argument('--reference', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--result', type=Path, required=True); p.add_argument('--history', type=Path, action='append', default=[])
    p.add_argument('--minimum-duration', type=float, default=20.)
    a = p.parse_args()
    if a.minimum_duration <= 0: p.error('--minimum-duration must be positive')
    a.result.mkdir(parents=True, exist_ok=False)
    def find(suffix):
        files = list(a.output.glob('*_hit_'+suffix+'.csv'))
        if len(files) != 1: raise ValueError(f'expected one {suffix} file in {a.output}')
        return files[0]
    records = csv_rows(find('history'))
    for path in a.history: records += csv_rows(path)
    by_time = {}
    for row in records:
        t = row['time']
        if t in by_time and abs(by_time[t]['K']-row['K']) > 1e-10*max(1, abs(row['K'])):
            raise ValueError('inconsistent overlapping histories')
        by_time[t] = row
    history = [by_time[t] for t in sorted(by_time)]
    report = {'case': a.case, 'physical_validation': 'not_automatically_assessed', 'history_range': [history[0]['time'], history[-1]['time']],
              'notes': ['Relative errors are comparisons, not pass/fail tolerances.', 'Coarse resolved viscous dissipation and derived Re_lambda omit unresolved gradients.']}
    if a.case == 'decay':
        reference = json.loads((a.reference/'cbc_table4.json').read_text())
        spectral = csv_rows(find('spectrum')); station_rows = []; errors = []
        for ref in reference['stations']:
            available = [r for r in history if abs(r['time']-ref['time']) < 1e-9]
            if not available:
                station_rows.append({'station':ref['station'], 'time':ref['time'], 'status':'not_reached', 'K':None, 'K_full_reference':ref['K'], 'K_relative_error':None}); continue
            row = available[-1]
            station_rows.append({'station':ref['station'], 'time':row['time'], 'status':'available', 'K':row['K'], 'K_full_reference':ref['K'], 'K_relative_error':(row['K']-ref['K'])/ref['K']})
            data = table(a.reference/f"cbc_spectrum_{ref['station']}_nondimensional.dat")
            frame = [s for s in spectral if abs(s['time']-row['time']) < 1e-9 and s['k'] > 0]
            for s in frame:
                target = interpolate(data, s['k'])
                if target is not None:
                    errors.append({'station':ref['station'], 'shell':s['shell'], 'k':s['k'], 'E':s['E'], 'E_reference_interpolated':target, 'relative_error':(s['E']-target)/target if target else None})
        report['stations'] = station_rows
        report['notes'] += ['K_full_reference is unfiltered experimental energy. It is not the exact target for a truncated coarse field.', 'Spectral errors are listed over overlapping k only; inspect low modes separately from the cutoff.']
        write_csv(a.result/'station_comparison.csv', station_rows); write_csv(a.result/'spectrum_comparison.csv', errors)
    else:
        means = {r['quantity']:r for r in csv_rows(find('means'))}
        weight = means['K']['weight']; mean = means['K']['time_mean']
        raw = table(a.reference/'ener_Re_time.txt'); unique = {}; duplicates = []
        for r in raw:
            if r[0] in unique:
                if unique[r[0]] != r: raise ValueError('conflicting reference duplicate')
                duplicates.append(r[0])
            unique[r[0]] = r
        dns = [{'time':r[0], 'K':r[1], 'Re_lambda':r[2]} for r in sorted(unique.values())]
        dns_mean = integrate(dns,'K',dns[0]['time'],dns[-1]['time'])
        report.update({'statistics_duration':weight, 'sufficient_duration_screen':weight >= a.minimum_duration,
                       'DNS_table_rows':len(raw), 'DNS_unique_rows':len(dns), 'DNS_duplicate_times':duplicates,
                       'DNS_table_time_range':[dns[0]['time'],dns[-1]['time']], 'DNS_table_K_time_mean':dns_mean,
                       'coarse_K_time_mean':mean if math.isfinite(mean) else None,
                       'K_relative_error':(mean-dns_mean)/dns_mean if math.isfinite(mean) else None})
        metadata = list(a.output.glob('*_hit_metadata.txt'))
        if len(metadata) != 1: raise ValueError('expected HIT metadata for averaging-window provenance')
        signature = next(line for line in metadata[0].read_text().splitlines() if line.startswith('hit_v1;')).split(';')
        finish = min(history[-1]['time'], float(signature[15])); begin = finish-weight
        if abs(begin-history[0]['time']) < 1e-12: begin = history[0]['time']
        report['accumulated_window'] = [begin, finish]
        block_means = [integrate(history,'K',begin+i*weight/5,begin+(i+1)*weight/5) for i in range(5)] if weight>0 else []
        if block_means and all(v is not None for v in block_means):
            report['five_block_K_means'] = block_means
            if weight >= a.minimum_duration:
                report['K_block_standard_error_estimate'] = statistics.stdev(block_means)/math.sqrt(5)
            else: report['uncertainty_status'] = 'duration_too_short_for_error_estimate'
            report['block_independence'] = 'not_established_by_this_script'
            report['block_duration'] = weight/5
        else: report['uncertainty_status'] = 'need_complete_history_for_accumulated_window'
        ref_s = dict((r[0],r[1]) for r in table(a.reference/'spectrum.txt')); spectrum = []
        for r in csv_rows(find('spectrum_mean')):
            if r['k'] in ref_s and math.isfinite(r['E_time_mean']):
                target = ref_s[r['k']]
                spectrum.append({'shell':r['shell'], 'k':r['k'], 'E_mean':r['E_time_mean'], 'E_DNS':target, 'relative_error':(r['E_time_mean']-target)/target})
        write_csv(a.result/'spectrum_comparison.csv',spectrum)
        report['notes'] += ['Duration screen alone does not prove stationarity or independent blocks; inspect block means and autocorrelation.', 'Independent random phases do not reproduce the DNS time trajectory. Compare distributions and means.', 'Public radial spectrum and scalar time table may use different averaging windows.']
    (a.result/'comparison.json').write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(report,ensure_ascii=True,indent=2,allow_nan=False))

if __name__ == '__main__': main()
