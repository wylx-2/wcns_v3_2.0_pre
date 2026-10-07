import hashlib
import importlib.util
import json
import lzma
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('restore_case_grids', Path(__file__).resolve().parents[1] / 'tools/restore_case_grids.py')
restore = importlib.util.module_from_spec(spec)
spec.loader.exec_module(restore)


class GridArchiveTests(unittest.TestCase):
    def setUp(self):
        work = Path(__file__).resolve().parents[1] / 'tmp/grid-archive-tests'
        work.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=work)
        self.root = Path(self.temp.name)
        folder = self.root / 'cases/grid-archives'
        folder.mkdir(parents=True)
        self.data = b'grid fixture\x00' * 50
        self.archive = folder / 'fixture.cgns.xz'
        self.archive.write_bytes(lzma.compress(self.data))
        self.item = {'path': 'cases/manual/sample/grids/input.cgns', 'bytes': len(self.data),
                     'sha256': hashlib.sha256(self.data).hexdigest(),
                     'archive': 'cases/grid-archives/fixture.cgns.xz',
                     'archive_sha256': restore.sha256(self.archive)}
        self.manifest = folder / 'manifest.json'
        self.write_manifest()

    def write_manifest(self):
        self.manifest.write_text(json.dumps({'schema_version': 1, 'grids': [self.item]}))

    def tearDown(self):
        self.temp.cleanup()

    def test_restore_and_verify_idempotent(self):
        restore.restore(self.root)
        target = self.root / self.item['path']
        self.assertEqual(target.read_bytes(), self.data)
        stamp = target.stat().st_mtime_ns
        restore.restore(self.root, verify_only=True)
        restore.restore(self.root)
        self.assertEqual(stamp, target.stat().st_mtime_ns)

    def test_existing_different_grid_is_preserved(self):
        target = self.root / self.item['path']
        target.parent.mkdir(parents=True)
        target.write_bytes(b'user data')
        with self.assertRaises(ValueError):
            restore.restore(self.root)
        self.assertEqual(target.read_bytes(), b'user data')

    def test_corrupt_archive_rejected_before_write(self):
        self.archive.write_bytes(b'corrupt')
        with self.assertRaises(ValueError):
            restore.restore(self.root)
        self.assertFalse((self.root / self.item['path']).exists())

    def test_destination_escape_rejected(self):
        self.item['path'] = '../escape.cgns'
        self.write_manifest()
        with self.assertRaises(ValueError):
            restore.restore(self.root)

    def test_verify_only_does_not_create_missing_grid(self):
        with self.assertRaises(FileNotFoundError):
            restore.restore(self.root, verify_only=True)
        self.assertFalse((self.root / self.item['path']).exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
