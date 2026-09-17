"""HTTP with retries.

Every reference lookup goes over the public internet, where a single dropped
connection should never cost the user a variant. Transient failures (timeouts,
connection resets, 429s and 5xx) are retried with backoff; genuine client
errors are raised immediately so the user sees the real message.
"""
from __future__ import annotations

import json
import random
import time
import urllib.error
import urllib.request

from .config import HTTP_TIMEOUT, USER_AGENT

RETRY_STATUS = {429, 500, 502, 503, 504}
MAX_ATTEMPTS = 4


class HttpError(RuntimeError):
    """A request that failed in a way retrying will not fix."""

    def __init__(self, message: str, status: int | None = None, body: str = ""):
        super().__init__(message)
        self.status = status
        self.body = body


def fetch(url: str, *, data: bytes | None = None, method: str | None = None,
          timeout: int = HTTP_TIMEOUT, accept: str = "application/json",
          attempts: int = MAX_ATTEMPTS) -> bytes:
    headers = {"Accept": accept, "User-Agent": USER_AGENT}
    if data is not None:
        headers["Content-Type"] = "application/json"
    last: Exception | None = None

    for attempt in range(attempts):
        req = urllib.request.Request(url, data=data, headers=headers, method=method)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as exc:
            body = ""
            try:
                body = exc.read().decode()[:400]
            except Exception:                                  # noqa: BLE001
                pass
            if exc.code not in RETRY_STATUS:
                detail = ""
                try:
                    detail = json.loads(body).get("error", "")
                except Exception:                              # noqa: BLE001
                    pass
                raise HttpError(detail or f"HTTP {exc.code}", exc.code, body) from exc
            last = exc
            # Honour Retry-After when the server sends one (Ensembl rate limiting).
            wait = exc.headers.get("Retry-After") if exc.headers else None
            delay = float(wait) if wait and str(wait).replace(".", "").isdigit() \
                else _backoff(attempt)
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as exc:
            last = exc
            delay = _backoff(attempt)

        if attempt < attempts - 1:
            time.sleep(min(delay, 12.0))

    reason = getattr(last, "reason", last)
    raise HttpError(f"Network request failed after {attempts} attempts ({reason}).")


def _backoff(attempt: int) -> float:
    return (0.6 * (2 ** attempt)) + random.uniform(0, 0.4)


def fetch_json(url: str, **kw) -> object:
    raw = fetch(url, **kw)
    try:
        return json.loads(raw.decode())
    except ValueError as exc:
        raise HttpError(f"Server returned malformed JSON: {exc}") from exc


def post_json(url: str, payload: dict, **kw) -> object:
    return fetch_json(url, data=json.dumps(payload).encode(), method="POST", **kw)
