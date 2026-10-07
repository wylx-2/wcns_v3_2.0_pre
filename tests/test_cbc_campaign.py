"""Linux orchestration tests using a fake MPI launcher, no fluid time integration.

Usage: python3 tests/test_cbc_campaign.py /absolute/staged/package /absolute/test/workdir
"""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

PACKAGE = Path(sys.argv.pop(1)).resolve()
WORK = Path(sys.argv.pop(1)).resolve()
WORK.mkdir(parents=True, exist_ok=True)

FAKE = r'''
import os, sys
from pathlib import Path
config = Path(sys.argv[sys.argv.index('--config') + 1])
v = dict(line.split(' = ', 1) for line in config.read_text().splitlines() if ' = ' in line)
assert v['algorithm.reconstruction'] == 'mdcd_linear'
assert v['algorithm.mdcd.diss'] == '0.001'
assert v['run.cfl'] == '0.3'
assert v['time.integrator'] == 'ssprk3'
if v['algorithm.riemann'] != 'roe_all_speed':
    assert not any(k.startswith('algorithm.roe_all_speed.') for k in v)
if '--dry-run' in sys.argv: sys.exit(0)
case = v['case.name']; out = Path(v['output.directory']); out.mkdir(parents=True)
mode = os.environ.get('CBC_FAKE_MODE', 'complete')
prep = out / 'preparation'; prep.mkdir()
(prep / (case + '-preparation.checkpoint.latest.cgns')).write_text('prep checkpoint')
if mode == 'pending-prep':
    root, prefix, phase = prep, case + '-preparation', 'preparation'
else:
    root, prefix, phase = out, case, 'production'
    (out / (case + '.checkpoint.latest.cgns')).write_text('production checkpoint')
reason = 'wall_time_checkpoint' if mode.startswith('pending') else 'physical_time_reached'
t = v['run.t_end'] if mode != 'false-complete' else '0.0000001'
(root / (prefix + '.manifest.r2.txt')).write_text('stop_reason=' + reason + '\ntime=' + t + '\n')
(root / (case + '_hit_metadata.txt')).write_text('phase=' + phase + '\n')
stations = v['hit.sample_times'].split(',')
if mode == 'missing-station': stations = stations[:-1]
(root / (case + '_hit_history.csv')).write_text('step,time\n' + ''.join(str(i) + ',' + t + '\n' for i,t in enumerate(stations)))
sys.exit(2 if mode.startswith('pending') else 3 if mode == 'failure' else 0)
'''


class CampaignTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=WORK)
        self.root = Path(self.temp.name)
        shutil.copytree(PACKAGE / 'scripts', self.root / 'scripts', ignore=shutil.ignore_patterns('__pycache__'))
        shutil.copytree(PACKAGE / 'cases', self.root / 'cases')
        shutil.copyfile(PACKAGE / 'campaign.json', self.root / 'campaign.json')
        (self.root / 'inputs/grids').mkdir(parents=True)
        (self.root / 'inputs/reference').mkdir()
        for n in [16, 32, 64, 128]:
            (self.root / f'inputs/grids/hit{n}.cgns').write_text('fake test grid')
        (self.root / 'inputs/reference/cbc_spectrum_42_nondimensional.dat').write_text('1 1\n')
        (self.root / 'build').mkdir()
        (self.root / 'build/wcns_run').write_text('fake binary for launcher test')
        self.fake = self.root / 'fake.py'
        self.fake.write_text(FAKE)

    def tearDown(self):
        self.temp.cleanup()

    def invoke(self, *args, mode='complete'):
        env = dict(os.environ, CBC_FAKE_MODE=mode)
        return subprocess.run([sys.executable, str(self.root / 'scripts/run_campaign.py'), '--np', '2',
                               '--launcher', sys.executable, '--launcher-arg', str(self.fake), *args],
                              env=env, text=True, capture_output=True, timeout=30)

    def state(self, folder='results'):
        return json.loads((self.root / folder / 'campaign_state.json').read_text())

    def test_nine_complete_and_skip(self):
        run = self.invoke()
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        self.assertEqual(len(self.state()['cases']), 9)
        self.assertTrue(all(x['status'] == 'complete' for x in self.state()['cases'].values()))
        self.assertEqual(self.invoke().returncode, 0)
        self.assertTrue(all(len(x['attempts']) == 1 for x in self.state()['cases'].values()))

    def test_pending_production_resume_prefers_production(self):
        self.assertEqual(self.invoke(mode='pending-production').returncode, 2)
        first = self.state()['cases']['cbc_n032_roe']
        self.assertEqual(first['status'], 'pending')
        self.assertEqual(self.invoke().returncode, 0)
        resumed = self.state()['cases']['cbc_n032_roe']['attempts'][1]
        self.assertTrue(resumed['restart_from'].endswith('/output/cbc_n032_roe.checkpoint.latest.cgns'))

    def test_pending_preparation_resume(self):
        self.assertEqual(self.invoke(mode='pending-prep').returncode, 2)
        self.assertEqual(self.invoke().returncode, 0)
        resumed = self.state()['cases']['cbc_n032_roe']['attempts'][1]
        self.assertIn('/preparation/', resumed['restart_from'])

    def test_return_zero_before_target_is_not_complete(self):
        self.assertEqual(self.invoke(mode='false-complete').returncode, 1)
        self.assertTrue(all(x['status'] == 'failed' for x in self.state()['cases'].values()))

    def test_missing_station_is_not_complete(self):
        self.assertEqual(self.invoke(mode='missing-station').returncode, 1)
        self.assertTrue(all(x['status'] == 'failed' for x in self.state()['cases'].values()))

    def test_failure_continues_other_cases_and_retry_is_fresh(self):
        self.assertEqual(self.invoke(mode='failure').returncode, 1)
        self.assertTrue(all(len(x['attempts']) == 1 for x in self.state()['cases'].values()))
        self.assertEqual(self.invoke('--retry-failed').returncode, 0)
        self.assertTrue(all(x['attempts'][1]['restart_from'] == '' for x in self.state()['cases'].values()))

    def test_changed_input_rejected(self):
        self.assertEqual(self.invoke().returncode, 0)
        path = self.root / 'cases/cbc_n032_roe.wcns'
        path.write_text(path.read_text() + '\n# changed\n')
        run = self.invoke()
        self.assertNotEqual(run.returncode, 0)
        self.assertIn('Input/binary/package path changed', run.stderr)

    def test_dry_run_separate_status(self):
        self.assertEqual(self.invoke('--dry-run').returncode, 0)
        self.assertTrue(all(x['status'] == 'validated' for x in self.state('results_dry-run')['cases'].values()))
        self.assertFalse((self.root / 'results').exists())

    def test_budget_exhaustion_does_not_start_solver(self):
        self.assertEqual(self.invoke('--job-seconds', '1', '--checkpoint-margin', '2').returncode, 2)
        self.assertTrue(all(x['status'] == 'not_started' for x in self.state()['cases'].values()))


if __name__ == '__main__':
    unittest.main(verbosity=2)
