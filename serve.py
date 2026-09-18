"""Run PrimerForge as a shared server (what the Linux service starts).

    PRIMERFORGE_AUTH=1 PRIMERFORGE_HOST=0.0.0.0 PRIMERFORGE_PORT=8080 python serve.py

Uses Waitress, a production WSGI server, in a single process: BLAST job workers
and the BLAT gfServers live in this process, so it must not be forked into
several workers. Configuration comes from the environment (see config.py); the
installer keeps it in /etc/primerforge/primerforge.env.
"""
from __future__ import annotations

import logging
import signal
import sys

from waitress import serve
from werkzeug.middleware.proxy_fix import ProxyFix

import app as webapp
from primerforge import auth, config


def main() -> int:
    # systemd captures stdout and stderr into the journal, so a plain stream
    # handler is all that is needed for 'journalctl -u primerforge' to show
    # tracebacks from the app with a timestamp against them.
    logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    webapp._startup()                                       # noqa: SLF001
    application = webapp.app
    if config.TRUST_PROXY:
        # Behind Caddy/nginx on this machine: trust one hop for client address and scheme.
        application.wsgi_app = ProxyFix(application.wsgi_app, x_for=1, x_proto=1, x_host=1)
    if config.AUTH_ENABLED and not any(u["is_admin"] for u in auth.list_users()):
        print("Warning: no administrator account exists. Create one with "
              "'primerforge users add NAME --admin'.", file=sys.stderr, flush=True)
    # systemd stops services with SIGTERM; exit normally so gfServers are stopped cleanly.
    signal.signal(signal.SIGTERM, lambda *_: sys.exit(0))
    print(f"PrimerForge listening on {config.HOST}:{config.PORT} "
          f"(logins {'on' if config.AUTH_ENABLED else 'OFF'}, data in {config.DATA_DIR})",
          flush=True)
    serve(application, host=config.HOST, port=config.PORT, threads=config.THREADS,
          ident="PrimerForge", channel_timeout=300,
          max_request_body_size=webapp.MAX_BODY_BYTES)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
