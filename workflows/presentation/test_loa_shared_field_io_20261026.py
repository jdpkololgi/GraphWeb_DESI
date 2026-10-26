"""Small CPU-only synthetic storage tests; no scientific data or checkpoints."""
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np
import torch

import loa_shared_field_io_20261026 as io
import loa_shared_sixfit_20261026 as s
from loa_shared_production_20261026 import checkpoint


class FieldIOTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.fg = s.Geometry(((0, 0, 0), (8, 0, 0)), (12, 8, 8), (0, 0, 0), (4, 8, 8))
        self.cg = s.Geometry(((0, 0, 0),), (3, 2, 2), (0, 0, 0), (3, 2, 2))
        self.coverage = io.coverage_for(self.fg)
        self.support = (self.coverage > 0).astype('f4') * .5
        self.chart = {'coarse_logrho': {'mean': [.1], 'std': [.3]}, 'fine_residual': {'mean': [0.], 'std': [.5]}}
        generator = torch.Generator().manual_seed(932)
        self.z = torch.randn((2, 1, 3, 2, 2), generator=generator)
        self.u = s.project(torch.randn((2, 1, 12, 8, 8), generator=generator))

    def write(self):
        path = self.root / 'field.h5'
        return io.write_field(path, self.z, self.u, self.coverage, self.support,
            self.cg, self.fg, self.chart, np.zeros(3), 'NGC', 'synthetic')

    def test_roundtrip_and_no_overwrite(self):
        record = self.write()
        self.assertTrue(record['stored_field_verified'])
        self.assertLess(record['stored_float32_mass_max_relative_error'], 2e-6)
        self.assertEqual(record['sha256'], io.sha(record['path']))
        with h5py.File(record['path']) as handle:
            self.assertTrue(np.isnan(handle['rho'][:, 4:8]).all())
            np.testing.assert_array_equal(handle['support'][:], self.support)
        with self.assertRaises(FileExistsError):
            self.write()

    def test_stored_mass_corruption_detected(self):
        record = self.write()
        with h5py.File(record['path'], 'r+') as handle:
            handle['rho'][0, 0, 0, 0] *= 4
        with self.assertRaisesRegex(ValueError, 'mass gate'):
            io.verify_field(record['path'], self.chart, np.zeros(3), [0, 0, 0], self.coverage, 2)

    def test_per_draw_mask_corruption_detected(self):
        record = self.write()
        with h5py.File(record['path'], 'r+') as handle:
            handle['rho'][1, 4, 0, 0] = 1
        with self.assertRaisesRegex(ValueError, 'mask mismatch'):
            io.verify_field(record['path'], self.chart, np.zeros(3), [0, 0, 0], self.coverage, 2)

    def test_parent_coverage_and_hole(self):
        points = np.array([[1.2, 1.3, 1.4], [9.1, 3.2, 4.1]]) * io.CELL
        self.assertEqual(io.check_coverage(self.coverage, points, np.zeros(3), 2)['parent_rows'], 2)
        with self.assertRaisesRegex(ValueError, 'uncovered'):
            io.check_coverage(self.coverage, np.array([[4., 2., 2.]]) * io.CELL, np.zeros(3), 1)
        with self.assertRaisesRegex(ValueError, 'count mismatch'):
            io.check_coverage(self.coverage, points, np.zeros(3), 1)
        self.coverage[0, 0, 0] = 0
        with self.assertRaisesRegex(ValueError, 'Partial mass'):
            io.check_coverage(self.coverage, points, np.zeros(3), 2)

    def test_support_overlap_not_filled(self):
        target = np.zeros((4, 4, 4), dtype='f4')
        seen = np.zeros_like(target, dtype=bool)
        crop = (slice(None),) * 3
        source = np.full_like(target, .5)
        source[0] = 0
        io.merge_support(target, seen, source, crop)
        io.merge_support(target, seen, source, crop)
        np.testing.assert_array_equal(target, source)
        with self.assertRaisesRegex(ValueError, 'Inconsistent support'):
            io.merge_support(target, seen, np.ones_like(source), crop)

    def test_atomic_metadata_and_checkpoint(self):
        path = self.root / 'receipt.json'
        io.atomic_json(path, {'step': 1})
        with self.assertRaises(FileExistsError):
            io.atomic_json(path, {'step': 2})
        self.assertEqual(s.read(path), {'step': 1})
        io.atomic_json(path, {'step': 2}, replace=True)
        self.assertEqual(s.read(path), {'step': 2})
        state = torch.arange(16).reshape(1, 1, 4, 2, 2)
        checkpoint(self.root / 'state.pt', state, 'test', 8, 'abc')
        checkpoint(self.root / 'state.pt', state + 1, 'test', 16, 'abc')
        restored = torch.load(self.root / 'state.pt', weights_only=False)
        self.assertEqual(restored['heun_step'], 16)
        self.assertFalse(restored['automatic_resume'])
        torch.testing.assert_close(restored['state'], state + 1)


if __name__ == '__main__':
    torch.set_num_threads(1)
    unittest.main()
