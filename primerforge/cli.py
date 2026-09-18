"""Command-line administration for a server install.

    python -m primerforge.cli users list
    python -m primerforge.cli users add alice [--admin] [--name "Alice Smith"]
    python -m primerforge.cli users passwd alice
    python -m primerforge.cli users role alice admin|user
    python -m primerforge.cli users disable|enable|delete alice
    python -m primerforge.cli users has-admin          # exit status 0 if one exists
    python -m primerforge.cli setup-data [--species homo_sapiens,mus_musculus] [--no-genome]
    python -m primerforge.cli backup [--dir /srv/backups] [--keep 14]
    python -m primerforge.cli status

On the server the `primerforge` command runs these as the service account with
the service's environment. New and reset passwords are generated and printed
once; the user has to choose their own at first sign-in. Pass --password-stdin
to set one instead (e.g. from a script).
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import subprocess
import sys
import time
from pathlib import Path

from . import auth, config, species, store

ROOT = Path(__file__).resolve().parent.parent


def _password(args) -> tuple[str, bool]:
    """(password, generated?)"""
    if getattr(args, "password_stdin", False):
        # Tolerate Windows line endings and a byte-order mark from the piping shell.
        value = sys.stdin.readline().lstrip("﻿").rstrip("\r\n")
        return value, False
    return auth.generate_password(), True


def _user_or_exit(name: str) -> dict:
    user = auth.find_user(name)
    if not user:
        print(f"No user called '{name}'.", file=sys.stderr)
        raise SystemExit(1)
    return user


def cmd_users(args) -> int:
    store.init()
    auth.init()
    try:
        if args.action == "has-admin":
            return 0 if any(u["is_admin"] and u["active"] for u in auth.list_users()) else 1
        if args.action == "list":
            users = auth.list_users()
            if not users:
                print("No users yet. Add one with: users add NAME --admin")
            for u in users:
                last = dt.datetime.fromtimestamp(u["last_login"]).strftime("%Y-%m-%d %H:%M") \
                    if u["last_login"] else "never"
                flags = ("" if u["active"] else " (disabled)") + \
                    (" (must change password)" if u["must_change"] else "")
                print(f"{u['username']:<24} {u['role']:<6} last sign-in {last}{flags}")
            return 0
        if args.action == "add":
            password, generated = _password(args)
            user = auth.create_user(args.username, password,
                                    role="admin" if args.admin else "user",
                                    full_name=args.name or "", must_change=generated)
            print(f"Created {user['role']} '{user['username']}'.")
            if generated:
                print(f"Temporary password: {password}")
                print("They will be asked to choose a new one when they first sign in.")
            return 0
        user = _user_or_exit(args.username)
        if args.action == "passwd":
            password, generated = _password(args)
            auth.set_password(user["id"], password, must_change=generated)
            print(f"Password for '{user['username']}' changed.")
            if generated:
                print(f"Temporary password: {password}")
        elif args.action == "role":
            auth.set_role(user["id"], args.role)
            print(f"'{user['username']}' is now {args.role}.")
        elif args.action in ("enable", "disable"):
            auth.set_active(user["id"], args.action == "enable")
            print(f"'{user['username']}' {args.action}d.")
        elif args.action == "delete":
            auth.delete_user(user["id"])
            print(f"Deleted '{user['username']}'.")
    except auth.AuthError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


def cmd_setup_data(args) -> int:
    """Everything the app needs, resumable: BLAST+, the GRCh38 genome for primer
    design, then each species' BLAST/BLAT databases and annotation."""
    failures = []
    steps = []
    if not args.no_genome:
        steps.append(("BLAST+ and the GRCh38 genome for primer design",
                      [sys.executable, str(ROOT / "tools" / "setup_genome.py")]))
    for name in [s.strip() for s in args.species.split(",") if s.strip()]:
        cmd = [sys.executable, str(ROOT / "tools" / "setup_species.py"), name]
        if args.components:
            cmd += ["--components", args.components]
        steps.append((f"BLAST/BLAT databases for {name}", cmd))
    for label, cmd in steps:
        # Every step is resumable, so a dropped connection just means trying again.
        for attempt in range(1, args.attempts + 1):
            print(f"\n==> {label}" + (f" (attempt {attempt})" if attempt > 1 else ""), flush=True)
            if subprocess.call(cmd, cwd=str(ROOT)) == 0:
                break
            if attempt < args.attempts:
                print(f"Step failed; retrying in {args.retry_wait} s", flush=True)
                time.sleep(args.retry_wait)
        else:
            failures.append(label)
    if failures:
        print("\nFinished with errors in: " + "; ".join(failures), file=sys.stderr, flush=True)
        return 1
    print("\nAll data installed.", flush=True)
    return 0


def cmd_backup(args) -> int:
    """A consistent copy of runs, tickets and accounts, safe while the server runs.

    Only the database is copied. Genomes and BLAST databases are large and can
    always be downloaded again with `setup-data`, so backing them up would cost
    tens of gigabytes to protect nothing that is not reproducible.
    """
    store.init()
    directory = Path(args.dir) if args.dir else config.DATA_DIR / "backups"
    stamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    target = store.backup(directory / f"primerforge-{stamp}.sqlite")
    removed = store.prune_cache()
    kept = sorted(directory.glob("primerforge-*.sqlite"))
    for old in kept[:-args.keep] if args.keep > 0 else []:
        old.unlink(missing_ok=True)
    size = target.stat().st_size / (1024 * 1024)
    print(f"Backed up to {target} ({size:.1f} MB)")
    if args.keep > 0 and len(kept) > args.keep:
        print(f"Removed {len(kept) - args.keep} older backup(s); keeping the newest {args.keep}.")
    if removed:
        print(f"Pruned {removed} expired cache row(s).")
    return 0


def cmd_status(_args) -> int:
    sys.path.insert(0, str(ROOT))
    from tools.setup_genome import read_status as genome_status
    from tools.setup_species import read_status as species_status
    print(f"Data directory:   {config.DATA_DIR}")
    print(f"Primer design:    genome {'ready' if config.genome_ready() else 'not installed'}, "
          f"BLAST {'ready' if config.blast_ready() else 'not ready'}")
    g = genome_status()
    if g.get("stage") not in (None, "idle"):
        print(f"  last setup step: {g.get('message')}")
    for s in species.all_species():
        comps = s["components"]
        done = [c for c, ok in comps.items() if ok]
        missing = [c for c, ok in comps.items() if not ok]
        print(f"{s['display_name']} ({s['name']}): {len(done)}/{len(comps)} components"
              + (f"; missing {', '.join(missing)}" if missing else ""))
        st = species_status(s["name"])
        if st and not st.get("done"):
            print(f"  installing: {st.get('component')}: {st.get('message')}")
        elif st and st.get("error"):
            print(f"  last install: {st['error']}")
    if config.AUTH_ENABLED:
        store.init()
        users = auth.list_users()
        print(f"Users:            {len(users)} ({sum(u['is_admin'] for u in users)} admin)")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="primerforge", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)

    users = sub.add_parser("users", help="manage accounts")
    us = users.add_subparsers(dest="action", required=True)
    us.add_parser("list")
    us.add_parser("has-admin")
    add = us.add_parser("add")
    add.add_argument("username")
    add.add_argument("--admin", action="store_true")
    add.add_argument("--name", default="")
    add.add_argument("--password-stdin", action="store_true")
    pw = us.add_parser("passwd")
    pw.add_argument("username")
    pw.add_argument("--password-stdin", action="store_true")
    role = us.add_parser("role")
    role.add_argument("username")
    role.add_argument("role", choices=auth.ROLES)
    for action in ("enable", "disable", "delete"):
        us.add_parser(action).add_argument("username")
    users.set_defaults(func=cmd_users)

    data = sub.add_parser("setup-data", help="download genomes and build databases")
    data.add_argument("--species", default="homo_sapiens",
                      help="comma-separated Ensembl names, e.g. homo_sapiens,mus_musculus")
    data.add_argument("--components", default="",
                      help="limit species components, e.g. genome,blastdb,blat")
    data.add_argument("--no-genome", action="store_true",
                      help="skip BLAST+ and the GRCh38 primer-design genome")
    data.add_argument("--attempts", type=int, default=4, help="tries per step")
    data.add_argument("--retry-wait", type=int, default=120, help="seconds between tries")
    data.set_defaults(func=cmd_setup_data)

    backup = sub.add_parser("backup", help="copy the database (runs, tickets, accounts)")
    backup.add_argument("--dir", default="", help="where to write it (default: DATA/backups)")
    backup.add_argument("--keep", type=int, default=14,
                        help="how many backups to keep; 0 keeps all")
    backup.set_defaults(func=cmd_backup)

    sub.add_parser("status", help="what is installed").set_defaults(func=cmd_status)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
