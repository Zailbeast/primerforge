#!/usr/bin/env bash
# PrimerForge server management, installed as /usr/local/bin/primerforge by install.sh.
set -euo pipefail

APP_DIR="@APP_DIR@"
CONF_DIR="@CONF_DIR@"
SERVICE_USER="@USER@"
ENV_FILE="$CONF_DIR/primerforge.env"

usage() {
  cat <<EOF
Usage: sudo primerforge COMMAND

  status                  Services, installed data and accounts
  logs                    Follow the web service log (Ctrl+C to stop)
  restart | stop | start  Control the web service
  data-progress           Follow the genome/database download (Ctrl+C to stop following)
  setup-data              Download and build the data for the configured species
                          (\$PRIMERFORGE_SPECIES in $ENV_FILE), in the background
  setup-data --species mus_musculus [--components genome,blastdb,blat]
                          Install other species now, in this terminal
  users list              Accounts
  users add NAME [--admin] [--name "Full Name"]
                          New account; prints a temporary password
  users passwd NAME       New temporary password (e.g. a forgotten one)
  users role NAME admin|user
  users disable|enable|delete NAME
  backup [--keep N]       Back up runs, tickets and accounts now (a nightly
                          backup already runs; see primerforge-backup.timer)
  version                 Installed version
  uninstall [--purge]     Remove PrimerForge (--purge also deletes data, settings and accounts)
EOF
}

need_root() {
  if [[ $EUID -ne 0 ]]; then
    echo "Run this with sudo: sudo primerforge $*" >&2
    exit 1
  fi
}

cli() {
  # shellcheck source=/dev/null
  (set -a; . "$ENV_FILE"; set +a
   cd "$APP_DIR/app"
   exec runuser -u "$SERVICE_USER" -- "$APP_DIR/venv/bin/python" -m primerforge.cli "$@")
}

unit_state() {
  local unit="$1" file="/etc/systemd/system/$1"
  [[ "$unit" == *.timer || "$unit" == *.service ]] || file="$file.service"
  if [[ -f "$file" ]]; then
    printf '  %-22s %s\n' "$unit" "$(systemctl is-active "$unit" 2>/dev/null || true)"
  fi
}

command="${1:-help}"
[[ $# -gt 0 ]] && shift
case "$command" in
  status)
    need_root "$command"
    echo "Services:"
    unit_state primerforge
    unit_state primerforge-proxy
    unit_state primerforge-setup
    unit_state primerforge-backup.timer
    echo
    cli status
    ;;
  logs)
    need_root "$command"
    exec journalctl -u primerforge -n 100 -f
    ;;
  restart|stop|start)
    need_root "$command"
    systemctl "$command" primerforge
    systemctl --no-pager --lines=0 status primerforge || true
    ;;
  data-progress)
    need_root "$command"
    if ! systemctl is-active --quiet primerforge-setup; then
      echo "The data setup is not running. Last messages:"
      journalctl -u primerforge-setup -n 25 --no-pager
      exit 0
    fi
    exec journalctl -u primerforge-setup -n 25 -f
    ;;
  setup-data)
    need_root "$command"
    if [[ $# -eq 0 ]]; then
      systemctl start --no-block primerforge-setup
      echo "Data setup running in the background. Follow it with: sudo primerforge data-progress"
    else
      cli setup-data "$@"
      systemctl try-restart primerforge
    fi
    ;;
  users)
    need_root "$command"
    cli users "$@"
    ;;
  backup)
    need_root "$command"
    cli backup "$@"
    ;;
  version)
    sed -n 's/^__version__ = "\(.*\)"/PrimerForge \1/p' "$APP_DIR/app/primerforge/__init__.py"
    ;;
  uninstall)
    need_root "$command"
    exec bash "$APP_DIR/installation/uninstall.sh" "$@"
    ;;
  help|-h|--help)
    usage
    ;;
  *)
    echo "Unknown command: $command" >&2
    usage >&2
    exit 2
    ;;
esac
