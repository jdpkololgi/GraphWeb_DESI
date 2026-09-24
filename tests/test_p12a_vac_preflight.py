"""Small fixtures only: no mock, DESI, checkpoint or protected payload access."""
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "workflows/catalog/p12a_vac_preflight.py"
spec = importlib.util.spec_from_file_location("vac_preflight", SCRIPT)
preflight = importlib.util.module_from_spec(spec)
spec.loader.exec_module(preflight)


class PreflightTests(unittest.TestCase):
    def test_hash_detects_same_size_corruption(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "checkpoint"
            path.write_bytes(b"abc")
            record = {"path": "checkpoint", "bytes": 3,
                      "sha256": hashlib.sha256(b"abc").hexdigest()}
            self.assertTrue(preflight.verify_record(record, Path(root))["verified"])
            path.write_bytes(b"abd")
            self.assertEqual(preflight.verify_record(record, Path(root))["status"], "hash_mismatch")

    def test_missing_is_not_verified(self):
        with tempfile.TemporaryDirectory() as root:
            record = {"path": "missing", "bytes": 3, "sha256": "a" * 64}
            self.assertFalse(preflight.verify_record(record, Path(root))["verified"])

    def test_large_artifact_not_read(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "large"
            path.write_bytes(b"abc")
            record = {"path": "large", "bytes": preflight.MAX_BYTES + 1, "sha256": "a" * 64}
            self.assertEqual(preflight.verify_record(record, Path(root))["status"], "over_read_limit")

    def test_historical_blind_never_becomes_fresh(self):
        p10 = {"model_phase_contract": {"training": ["ph000"],
               "validation_and_selection": ["ph006"], "sealed_blind_test": ["ph001"]}}
        e2e = {"phase_roles": {"ph007": "train", "ph014": "confirmation"}}
        rows = {row["phase"]: row for row in preflight.phase_ledger(p10, e2e)}
        self.assertIn("opened_once", rows["ph001"]["p12_role"])
        self.assertFalse(rows["ph001"]["fresh_blind_eligible"])
        self.assertTrue(rows["ph014"]["reserved_confirmation"])
        self.assertFalse(rows["ph014"]["replication_candidate_pending_exposure_audit"])
        self.assertTrue(rows["ph007"]["replication_candidate_pending_exposure_audit"])
        self.assertFalse(any(row["payload_access_authorized_by_this_report"] for row in rows.values()))

    def test_catalogue_inventory_never_asserts_live_verification(self):
        rows = list(preflight.catalogue_records({"full": [{"path": "/not/mounted.fits",
                       "columns": ["TARGETID", "Z"], "sha256": "a" * 64}]}))
        self.assertEqual(len(rows), 1)
        self.assertFalse(rows[0]["live_verified"])


if __name__ == "__main__":
    unittest.main()
