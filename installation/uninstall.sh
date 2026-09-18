#!/usr/bin/env bash
# Remove PrimerForge from this server.
#
#   sudo primerforge uninstall            # keep data, settings and accounts
#   sudo primerforge uninstall --purge    # delete everything
set -euo pipefail

PURGE=0
case "${1:-}" in
  --purge) PURGE=1 ;;
  "") ;;
  *) echo "Usage: $0 [--purge]" >&2; exit 2 ;;
esac
[[ $EUID -eq 0 ]] || { echo "Run this as root, e.g. with sudo." >&2; exit 1; }

APP_DIR=/opt/primerforge
CONF_DIR=/etc/primerforge
DATA_DIR="$(sed -n 's/^PRIMERFORGE_DATA=//p' "$CONF_DIR/primerforge.env" 2>/dev/null | head -1)"
DATA_DIR="${DATA_DIR:-/var/lib/primerforge}"

if (( PURGE )); then
  echo "This deletes PrimerForge, all genome data in $DATA_DIR, settings and user accounts."
  if [[ -t 0 ]]; then
    read -r -p "Type 'delete' to continue: " answer
    [[ "$answer" == "delete" ]] || { echo "Cancelled."; exit 1; }
  fi
fi

systemctl disable --now primerforge-backup.timer >/dev/null 2>&1 || true
rm -f /etc/systemd/system/primerforge-backup.timer
for unit in primerforge-setup primerforge-proxy primerforge-backup primerforge; do
  systemctl disable --now "$unit" >/dev/null 2>&1 || true
  rm -f "/etc/systemd/system/$unit.service"
done
systemctl daemon-reload
rm -f /usr/local/bin/primerforge
rm -rf "$APP_DIR"
echo "Removed the PrimerForge services and $APP_DIR."

if (( PURGE )); then
  rm -rf "$DATA_DIR" "$CONF_DIR" /var/lib/primerforge-caddy
  userdel primerforge >/dev/null 2>&1 || true
  userdel primerforge-caddy >/dev/null 2>&1 || true
  echo "Deleted $DATA_DIR, $CONF_DIR and the service accounts."
else
  echo "Kept genome data, accounts and history in $DATA_DIR and settings in $CONF_DIR."
  echo "Installing again reuses them; use --purge to delete them."
fi
echo "Any firewall ports opened for PrimerForge were left open."
