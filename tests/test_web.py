"""The web layer: health, error handling and the body-size limit.

The app is imported with logins off (the desktop default), and startup is
marked as already done so no test reaches for a genome, a BLAST database or
the network.
"""
from __future__ import annotations

import logging
import unittest

import app as webapp

webapp._started = True                                   # noqa: SLF001


@webapp.app.route("/__raises__")
def _raises():                                           # pragma: no cover - test fixture
    raise RuntimeError("deliberate test failure")


class WebTestCase(unittest.TestCase):
    def setUp(self):
        self.client = webapp.app.test_client()


class Health(WebTestCase):
    def test_health_is_public_and_describes_the_server(self):
        response = self.client.get("/api/health")
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertTrue(payload["ok"])                   # what the installer waits for
        self.assertFalse(payload["auth"])                # logins off in the desktop default

    def test_health_needs_no_session(self):
        self.assertIn("api_health", webapp.PUBLIC_ENDPOINTS)


class NotFound(WebTestCase):
    def test_api_paths_get_json(self):
        response = self.client.get("/api/no-such-thing")
        self.assertEqual(response.status_code, 404)
        self.assertIn("error", response.get_json())

    def test_pages_get_a_readable_page(self):
        response = self.client.get("/no-such-page")
        self.assertEqual(response.status_code, 404)
        self.assertEqual(response.mimetype, "text/html")
        self.assertIn("Not found", response.get_data(as_text=True))

    def test_a_client_asking_for_json_gets_json(self):
        response = self.client.get("/no-such-page", headers={"Accept": "application/json"})
        self.assertIn("error", response.get_json())


class TooLarge(WebTestCase):
    def test_oversized_body_is_refused_before_it_is_read(self):
        original = webapp.app.config["MAX_CONTENT_LENGTH"]
        webapp.app.config["MAX_CONTENT_LENGTH"] = 512    # so the test sends 1 KB, not 32 MB
        try:
            response = self.client.post("/api/blast/parse", data="A" * 1024,
                                        content_type="application/json")
        finally:
            webapp.app.config["MAX_CONTENT_LENGTH"] = original
        self.assertEqual(response.status_code, 413)
        self.assertIn("error", response.get_json())

    def test_the_limit_is_generous_enough_for_a_full_blast_form(self):
        from primerforge import blastconf
        biggest = blastconf.MAX_SEQUENCE_LENGTH * blastconf.MAX_NUM_SEQUENCES
        self.assertGreater(webapp.MAX_BODY_BYTES, biggest)


class UnhandledErrors(WebTestCase):
    def test_traceback_is_logged_and_the_user_sees_a_message(self):
        logging.disable(logging.CRITICAL)                # the handler logs the traceback
        try:
            response = self.client.get("/__raises__")
        finally:
            logging.disable(logging.NOTSET)
        self.assertEqual(response.status_code, 500)
        body = response.get_data(as_text=True)
        self.assertIn("Something went wrong", body)
        self.assertNotIn("deliberate test failure", body)   # no traceback to the browser

    def test_http_errors_keep_their_own_status(self):
        # The catch-all must not turn a 404 into a 500.
        self.assertEqual(self.client.get("/api/no-such-thing").status_code, 404)


class Cookies(unittest.TestCase):
    def test_session_cookie_is_hardened(self):
        config = webapp.app.config
        self.assertTrue(config["SESSION_COOKIE_HTTPONLY"])
        self.assertEqual(config["SESSION_COOKIE_SAMESITE"], "Lax")


if __name__ == "__main__":
    unittest.main()
