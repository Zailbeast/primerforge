"""The database layer: the cache, backups, and recovery after a restart."""
from __future__ import annotations

import sqlite3
import time
import unittest

from primerforge import config, pipeline, store

from . import DATA_DIR


class StoreTestCase(unittest.TestCase):
    """Each test starts from an empty database in the throwaway data directory."""

    def setUp(self):
        store.init()
        with store.db() as conn:
            for table in ("pairs", "variants", "runs", "cache"):
                conn.execute(f"DELETE FROM {table}")

    def _run(self, **kw) -> str:
        return store.create_run(label=kw.get("label", "test"), assay="sanger",
                                input_raw="rs334", params={"assembly": "GRCh38"},
                                n_input=kw.get("n_input", 1), specificity=False)


class Cache(StoreTestCase):
    def test_round_trip(self):
        store.cache_put("k", {"a": 1})
        self.assertEqual(store.cache_get("k"), {"a": 1})

    def test_expired_entry_is_not_returned(self):
        store.cache_put("k", {"a": 1})
        self.assertIsNone(store.cache_get("k", max_age=-1))

    def test_prune_removes_only_expired_rows(self):
        store.cache_put("fresh", 1)
        store.cache_put("stale", 2)
        with store.db() as conn:                       # age one row past every max_age
            conn.execute("UPDATE cache SET created_at=? WHERE key='stale'",
                         (time.time() - store.CACHE_MAX_AGE - 10,))
        self.assertEqual(store.prune_cache(), 1)
        self.assertEqual(store.cache_get("fresh"), 1)
        with store.db() as conn:
            self.assertEqual(conn.execute("SELECT COUNT(*) FROM cache").fetchone()[0], 1)


class Backup(StoreTestCase):
    def test_backup_is_a_readable_copy(self):
        run_id = self._run(label="backed up")
        target = store.backup(DATA_DIR / "backups" / "test.sqlite")
        self.assertTrue(target.exists())
        conn = sqlite3.connect(target)
        try:
            row = conn.execute("SELECT label FROM runs WHERE id=?", (run_id,)).fetchone()
        finally:
            conn.close()
        self.assertEqual(row[0], "backed up")

    def test_backup_creates_its_directory(self):
        target = store.backup(DATA_DIR / "nested" / "deeper" / "test.sqlite")
        self.assertTrue(target.parent.is_dir())


class RestartRecovery(StoreTestCase):
    """A design run lives on a worker thread, so a restart must not leave it 'running'."""

    def test_interrupted_run_is_closed_with_what_it_produced(self):
        run_id = self._run(n_input=3)
        store.add_variant(run_id, 0, "rs334", "ok")
        store.add_variant(run_id, 1, "rs1801133", "failed", error="nope")
        # The process dies here: no finish_run was ever called.
        self.assertEqual(pipeline.recover(), 1)

        run = store.get_run(run_id)
        self.assertEqual(run["status"], "partial")
        self.assertEqual((run["n_ok"], run["n_failed"]), (1, 1))
        self.assertIn("restarted", run["error"])
        self.assertIsNotNone(run["finished_at"])

    def test_run_with_nothing_finished_is_marked_failed(self):
        run_id = self._run()
        pipeline.recover()
        self.assertEqual(store.get_run(run_id)["status"], "failed")

    def test_finished_runs_are_left_alone(self):
        run_id = self._run()
        store.finish_run(run_id, "ok", 1, 0)
        self.assertEqual(pipeline.recover(), 0)
        run = store.get_run(run_id)
        self.assertEqual(run["status"], "ok")
        self.assertIsNone(run["error"])

    def test_recovery_is_safe_to_run_twice(self):
        self._run()
        self.assertEqual(pipeline.recover(), 1)
        self.assertEqual(pipeline.recover(), 0)


class DataDirectory(unittest.TestCase):
    def test_tests_never_touch_the_real_data_directory(self):
        # Guards the environment set up in tests/__init__.py: without it these
        # tests would write runs into the developer's own history.
        self.assertEqual(config.DATA_DIR, DATA_DIR)


if __name__ == "__main__":
    unittest.main()
