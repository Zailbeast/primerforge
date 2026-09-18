#!/usr/bin/env bash
# PrimerForge server installer. Run it from the installation folder:
#
#   sudo bash install.sh                                # HTTP on port 8080
#   sudo bash install.sh --domain primers.mylab.org     # HTTPS with a Let's Encrypt certificate
#   sudo bash install.sh --self-signed                  # HTTPS with a self-signed certificate
#
# The folder should first be filled on a computer with internet access by
# 'python installation/prepare.py', which adds the application (app/) and the
# dependencies (dependencies/: Python packages, NCBI BLAST+, UCSC BLAT, Caddy). The
# installer uses those bundled files and downloads only what is missing.
#
# It installs the app in /opt/primerforge, a systemd service with logins, an
# administrator account and (optionally) a Caddy HTTPS proxy, then starts
# downloading the genome data in the background (about an hour on a fast
# connection; 11 GB for human). Re-running it upgrades the app and keeps data,
# users and settings.
set -Eeuo pipefail

APP_DIR=/opt/primerforge
DATA_DIR=/var/lib/primerforge
CONF_DIR=/etc/primerforge
SERVICE_USER=primerforge
PORT=8080
DOMAIN=""
EMAIL=""
SELF_SIGNED=0
SPECIES=homo_sapiens
WITH_DATA=1
ADMIN_NAME="admin"
OPEN_FIREWALL=1
CADDY_FALLBACK_VERSION=2.11.4

SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEPS="$SRC_DIR/dependencies"

usage() {
  cat <<EOF
Usage: sudo bash install.sh [options]

  --domain NAME       Serve https://NAME with a free Let's Encrypt certificate.
                      NAME must point at this server and ports 80/443 must be reachable.
  --email ADDRESS     Contact address for Let's Encrypt expiry notices (with --domain).
  --self-signed       Serve HTTPS with a self-signed certificate (for internal networks
                      without a public domain; browsers show a warning the first time).
  --port N            Port for the app (default $PORT). Without HTTPS this is the port
                      people open; with HTTPS it only listens on 127.0.0.1.
  --data-dir DIR      Where genomes and databases go (default $DATA_DIR). Human uses
                      11 GB (about 15 GB while installing); 25 GB free is required.
  --species LIST      Species to install, comma-separated Ensembl names
                      (default $SPECIES), e.g. homo_sapiens,mus_musculus
  --no-data           Install the app only; download data later with
                      'sudo primerforge setup-data' or from Settings.
  --admin NAME        Name of the first administrator account (default $ADMIN_NAME).
  --no-firewall       Do not open ports in ufw/firewalld.
  -h, --help          Show this help.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --domain) DOMAIN="${2:?--domain needs a name}"; shift 2 ;;
    --email) EMAIL="${2:?--email needs an address}"; shift 2 ;;
    --self-signed) SELF_SIGNED=1; shift ;;
    --port) PORT="${2:?--port needs a number}"; shift 2 ;;
    --data-dir) DATA_DIR="${2:?--data-dir needs a directory}"; shift 2 ;;
    --species) SPECIES="${2:?--species needs a list}"; shift 2 ;;
    --no-data) WITH_DATA=0; shift ;;
    --admin) ADMIN_NAME="${2:?--admin needs a name}"; shift 2 ;;
    --no-firewall) OPEN_FIREWALL=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------
if [[ -t 1 ]]; then B=$'\e[1m'; G=$'\e[32m'; Y=$'\e[33m'; R=$'\e[31m'; N=$'\e[0m'; else B=""; G=""; Y=""; R=""; N=""; fi
step() { echo; echo "${B}==> $*${N}"; }
ok()   { echo "  ${G}✓${N} $*"; }
warn() { echo "  ${Y}!${N} $*"; }
die()  { echo; echo "${R}${B}Installation stopped:${N} $*" >&2; exit 1; }
trap 'die "command failed on line $LINENO: $BASH_COMMAND"' ERR

host_name() { hostname -f 2>/dev/null || hostname 2>/dev/null || cat /proc/sys/kernel/hostname; }
host_ips() {
  hostname -I 2>/dev/null && return 0
  ip -o addr show scope global 2>/dev/null | awk '{ split($4, a, "/"); printf "%s ", a[1] }'
}
first_ip() { local ips; ips="$(host_ips)"; ips="${ips# }"; echo "${ips%% *}"; }

# bundled SUBDIR PATTERN: the newest bundled dependency matching PATTERN, if there is one.
bundled() {
  local f last=""
  for f in "$DEPS/$1"/$2; do [[ -f "$f" ]] && last="$f"; done
  [[ -n "$last" ]] && echo "$last"
}

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------
step "Checking this server"
[[ $EUID -eq 0 ]] || die "run this as root: sudo bash install.sh"
[[ "$PORT" =~ ^[0-9]+$ && "$PORT" -ge 1024 && "$PORT" -le 65535 ]] || die "--port must be a number from 1024 to 65535."
[[ -z "$DOMAIN" || "$SELF_SIGNED" -eq 0 ]] || die "choose either --domain or --self-signed, not both."
[[ "$DATA_DIR" == /* ]] || die "--data-dir must be an absolute path."
[[ "$SPECIES" =~ ^[a-z0-9_,]+$ ]] || die "--species takes Ensembl names such as homo_sapiens,mus_musculus."
[[ -d /run/systemd/system ]] || die "this installer needs systemd (most Linux servers have it)."

# The application: app/ inside this folder (prepared), or the project this folder sits in.
if [[ -f "$SRC_DIR/app/serve.py" && -d "$SRC_DIR/app/primerforge" ]]; then
  APP_SRC="$SRC_DIR/app"
elif [[ -f "$SRC_DIR/../serve.py" && -d "$SRC_DIR/../primerforge" ]]; then
  APP_SRC="$(cd "$SRC_DIR/.." && pwd)"
else
  die "the application files are missing. On a computer with internet access run 'python installation/prepare.py', then copy the whole installation folder here."
fi

# shellcheck source=/dev/null
if [[ -r /etc/os-release ]]; then . /etc/os-release; fi
ok "${PRETTY_NAME:-Linux} ($(uname -m))"

ARCH="$(uname -m)"
[[ "$ARCH" == "x86_64" ]] || die "UCSC publishes BLAT for x86_64 Linux only; this server is $ARCH."

MEM_GB=$(( $(awk '/^MemTotal:/ {print $2}' /proc/meminfo) / 1024 / 1024 ))
if (( MEM_GB < 8 )); then warn "${MEM_GB} GB of memory; 16 GB is recommended (BLAT keeps the genome index in memory)."
else ok "${MEM_GB} GB of memory"; fi

probe="$DATA_DIR"; while [[ ! -d "$probe" ]]; do probe="$(dirname "$probe")"; done
FREE_GB=$(( $(df -Pk "$probe" | awk 'NR==2 {print $4}') / 1024 / 1024 ))
if (( WITH_DATA )) && (( FREE_GB < 25 )); then
  die "only ${FREE_GB} GB free for $DATA_DIR; the human data uses 11 GB (about 15 GB while installing), so keep at least 25 GB free. Free space, choose another --data-dir, or use --no-data."
elif (( FREE_GB < 40 )); then warn "${FREE_GB} GB free for $DATA_DIR; 40 GB is recommended if you will add more species."
else ok "${FREE_GB} GB free for $DATA_DIR"; fi

if [[ -f "$DEPS/SHA256SUMS" ]]; then
  (cd "$DEPS" && sed 's/\r$//' SHA256SUMS | sha256sum --quiet -c -) \
    || die "some bundled dependencies are damaged or incomplete (listed above). Copy the installation folder to the server again."
  ok "bundled dependencies verified ($(grep -c . "$DEPS/SHA256SUMS") files)"
else
  warn "no bundled dependencies in $DEPS; they will be downloaded (run prepare.py to bundle them)"
fi

# ---------------------------------------------------------------------------
# System packages
# ---------------------------------------------------------------------------
python_ok() { "$1" -c 'import sys, venv, ensurepip; sys.exit(0 if sys.version_info >= (3, 9) else 1)' 2>/dev/null; }
find_python() {
  for c in python3.14 python3.13 python3.12 python3.11 python3.10 python3.9 python3; do
    if command -v "$c" >/dev/null && python_ok "$(command -v "$c")"; then command -v "$c"; return 0; fi
  done
  return 1
}
# Everything the installer and BLAST/BLAT need from the operating system.
runtime_ready() {
  local tool
  [[ -n "$(find_python || true)" ]] || return 1
  for tool in tar gzip curl pgrep sha256sum runuser; do command -v "$tool" >/dev/null || return 1; done
  ldconfig -p 2>/dev/null | grep -q 'libgomp\.so\.1' || return 1
  [[ -d /etc/ssl/certs || -d /etc/pki/tls/certs ]]
}

step "Checking system packages"
PKG=""
if command -v apt-get >/dev/null; then PKG=apt
elif command -v dnf >/dev/null; then PKG=dnf
elif command -v yum >/dev/null; then PKG=yum
elif command -v zypper >/dev/null; then PKG=zypper; fi

if runtime_ready; then
  ok "Python 3.9+, OpenMP runtime and tools already present"
else
  [[ -n "$PKG" ]] || die "missing system packages and no supported package manager (apt, dnf, yum or zypper). Install Python 3.9+ with venv, libgomp, tar, gzip, curl and procps."
  # Only ask for curl when there is none: RHEL/Rocky ship curl-minimal, which the full
  # curl package conflicts with.
  CURL_PKG="curl"; command -v curl >/dev/null && CURL_PKG=""
  case "$PKG" in
    apt)
      export DEBIAN_FRONTEND=noninteractive
      apt-get update -qq || warn "apt-get update failed; trying the package lists already on this server"
      apt-get install -y -qq python3 python3-venv python3-pip ca-certificates $CURL_PKG tar gzip procps libgomp1 >/dev/null
      ;;
    dnf|yum)
      $PKG install -y -q python3 python3-pip ca-certificates $CURL_PKG tar gzip procps-ng libgomp >/dev/null
      ;;
    zypper)
      zypper --non-interactive --quiet install python3 python3-pip ca-certificates $CURL_PKG tar gzip procps libgomp1 >/dev/null
      ;;
  esac
  ok "installed python3, the OpenMP runtime and tools with $PKG"
fi

PY="$(find_python || true)"
if [[ -z "$PY" ]]; then
  # Older releases (e.g. RHEL/Rocky 8, openSUSE Leap) ship Python 3.6 as python3 but offer newer ones.
  warn "the default Python is older than 3.9; installing a newer one"
  case "$PKG" in
    dnf|yum) for v in python3.12 python3.11 python39; do $PKG install -y -q "$v" >/dev/null 2>&1 && break; done ;;
    zypper) for v in python312 python311; do zypper --non-interactive --quiet install "$v" >/dev/null 2>&1 && break; done ;;
    apt) apt-get install -y -qq python3.11-venv >/dev/null 2>&1 || apt-get install -y -qq python3.10-venv >/dev/null 2>&1 || true ;;
  esac
  PY="$(find_python || true)"
  [[ -n "$PY" ]] || die "PrimerForge needs Python 3.9 or newer, and none could be installed automatically."
fi
PY_VERSION="$("$PY" -c 'import sys; print("%d.%d" % sys.version_info[:2])')"
ok "Python $PY_VERSION at $PY"

# ---------------------------------------------------------------------------
# Service account and directories
# ---------------------------------------------------------------------------
step "Creating the service account and directories"
if ! id "$SERVICE_USER" >/dev/null 2>&1; then
  NOLOGIN="$(command -v nologin || echo /bin/false)"
  useradd --system --user-group --home-dir "$DATA_DIR" --no-create-home --shell "$NOLOGIN" "$SERVICE_USER"
  ok "created system user $SERVICE_USER"
else
  ok "system user $SERVICE_USER exists"
fi
install -d -m 0755 "$APP_DIR"
install -d -m 0750 -o "$SERVICE_USER" -g "$SERVICE_USER" "$DATA_DIR"
install -d -m 0755 "$CONF_DIR"          # secrets inside keep their own permissions
ok "app $APP_DIR, data $DATA_DIR, configuration $CONF_DIR"

# ---------------------------------------------------------------------------
# Application
# ---------------------------------------------------------------------------
step "Installing PrimerForge"
if systemctl is-active --quiet primerforge 2>/dev/null; then
  systemctl stop primerforge
  ok "stopped the running version"
fi

# swap_in STAGE TARGET: replace TARGET with the freshly staged copy.
swap_in() {
  rm -rf "$2.old"
  if [[ -d "$2" ]]; then mv "$2" "$2.old"; fi
  mv "$1" "$2"
  rm -rf "$2.old"
}

STAGE="$APP_DIR/app.new"
rm -rf "$STAGE" && mkdir -p "$STAGE"
tar -C "$APP_SRC" --exclude=./data --exclude=./dist --exclude=./.git --exclude=./installation \
    --exclude='__pycache__' --exclude='*.pyc' -cf - . | tar -C "$STAGE" -xf -
chown -R root:root "$STAGE"
chmod -R u=rwX,go=rX "$STAGE"
swap_in "$STAGE" "$APP_DIR/app"
VERSION="$(sed -n 's/^__version__ = "\(.*\)"/\1/p' "$APP_DIR/app/primerforge/__init__.py")"
ok "PrimerForge ${VERSION:-} in $APP_DIR/app"

# The installer itself is kept for 'primerforge uninstall' and future reference.
ISTAGE="$APP_DIR/installation.new"
rm -rf "$ISTAGE" && install -d "$ISTAGE/services"
install -m 0755 "$SRC_DIR/install.sh" "$SRC_DIR/uninstall.sh" "$SRC_DIR/primerforge.sh" "$ISTAGE/"
install -m 0644 "$SRC_DIR"/services/*.service "$SRC_DIR"/services/*.timer "$ISTAGE/services/"
sed -i 's/\r$//' "$ISTAGE"/*.sh "$ISTAGE"/services/*   # in case Windows line endings crept in
swap_in "$ISTAGE" "$APP_DIR/installation"

# A venv made by a different Python (e.g. after an OS upgrade) is rebuilt.
if [[ -x "$APP_DIR/venv/bin/python" ]] && ! "$APP_DIR/venv/bin/python" -c 'import sys' >/dev/null 2>&1; then
  rm -rf "$APP_DIR/venv"
fi
if [[ ! -x "$APP_DIR/venv/bin/python" ]] || \
   [[ "$("$APP_DIR/venv/bin/python" -c 'import sys; print("%d.%d" % sys.version_info[:2])')" != "$PY_VERSION" ]]; then
  rm -rf "$APP_DIR/venv"
  "$PY" -m venv "$APP_DIR/venv"
fi
PIP_LOG=/tmp/primerforge-pip.log
: > "$PIP_LOG"
PIP_DONE=0
if [[ -n "$(find "$DEPS/python" -maxdepth 1 -name '*.whl' -print -quit 2>/dev/null)" ]]; then
  if "$APP_DIR/venv/bin/python" -m pip install --disable-pip-version-check -q --no-index \
       --find-links "$DEPS/python" -r "$APP_DIR/app/requirements.txt" >>"$PIP_LOG" 2>&1; then
    PIP_DONE=1
    ok "Python packages installed from the bundled dependencies"
  else
    warn "the bundled Python packages do not cover Python $PY_VERSION; downloading them instead"
  fi
fi
if (( PIP_DONE == 0 )); then
  "$APP_DIR/venv/bin/python" -m pip install --disable-pip-version-check -q --upgrade pip wheel >>"$PIP_LOG" 2>&1 || true
  if ! "$APP_DIR/venv/bin/python" -m pip install --disable-pip-version-check -q \
         -r "$APP_DIR/app/requirements.txt" >>"$PIP_LOG" 2>&1; then
    # No prebuilt wheel for this platform: install a compiler and build primer3-py.
    warn "building Python packages from source (installing a compiler first)"
    case "$PKG" in
      apt) apt-get install -y -qq build-essential python3-dev >/dev/null ;;
      dnf|yum) $PKG install -y -q gcc gcc-c++ make python3-devel >/dev/null ;;
      zypper) zypper --non-interactive --quiet install gcc gcc-c++ make python3-devel >/dev/null ;;
      *) die "could not install the Python packages and no package manager to add a compiler; see $PIP_LOG" ;;
    esac
    "$APP_DIR/venv/bin/python" -m pip install --disable-pip-version-check -q \
      -r "$APP_DIR/app/requirements.txt" >>"$PIP_LOG" 2>&1 \
      || die "could not install the Python packages; see $PIP_LOG"
  fi
  ok "Python packages downloaded and installed"
fi
"$APP_DIR/venv/bin/python" -m compileall -q "$APP_DIR/app" >/dev/null || true
ok "Python environment ready ($APP_DIR/venv)"

# ---------------------------------------------------------------------------
# BLAST+ and BLAT programs
# ---------------------------------------------------------------------------
step "Installing BLAST+ and BLAT"
install -d -m 0750 -o "$SERVICE_USER" -g "$SERVICE_USER" "$DATA_DIR/blast" "$DATA_DIR/blat_bin"
if [[ -n "$(find "$DATA_DIR/blast" -type f -name blastn -print -quit 2>/dev/null)" ]]; then
  ok "BLAST+ already installed"
elif BLAST_TGZ="$(bundled blast 'ncbi-blast-*-x64-linux.tar.gz')"; then
  tar -C "$DATA_DIR/blast" -xzf "$BLAST_TGZ"
  chown -R "$SERVICE_USER":"$SERVICE_USER" "$DATA_DIR/blast"
  ok "BLAST+ $(basename "$BLAST_TGZ" | sed 's/ncbi-blast-\(.*\)+-x64-linux.tar.gz/\1/') from the bundled dependencies"
else
  ok "BLAST+ will be downloaded with the genome data"
fi
BLAT_COUNT=0
for program in gfServer gfClient blat faToTwoBit twoBitInfo; do
  if [[ -f "$DEPS/blat/$program" ]]; then
    install -m 0755 -o "$SERVICE_USER" -g "$SERVICE_USER" "$DEPS/blat/$program" "$DATA_DIR/blat_bin/$program"
    BLAT_COUNT=$((BLAT_COUNT + 1))
  fi
done
if (( BLAT_COUNT == 5 )); then ok "BLAT programs from the bundled dependencies"
else ok "BLAT programs will be downloaded with the genome data"; fi

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
step "Writing configuration"
ENV_FILE="$CONF_DIR/primerforge.env"
SECRET=""
if [[ -f "$ENV_FILE" ]]; then
  SECRET="$(sed -n 's/^PRIMERFORGE_SECRET_KEY=//p' "$ENV_FILE" | head -1)"
fi
[[ -n "$SECRET" ]] || SECRET="$("$PY" -c 'import secrets; print(secrets.token_hex(32))')"

USE_PROXY=0
if [[ -n "$DOMAIN" || "$SELF_SIGNED" -eq 1 ]]; then USE_PROXY=1; fi
if (( USE_PROXY )); then HOST=127.0.0.1; else HOST=0.0.0.0; fi

umask 027
cat > "$ENV_FILE" <<EOF
# PrimerForge service settings. Restart after editing: sudo systemctl restart primerforge
PRIMERFORGE_DATA=$DATA_DIR
PRIMERFORGE_AUTH=1
PRIMERFORGE_HOST=$HOST
PRIMERFORGE_PORT=$PORT
PRIMERFORGE_TRUST_PROXY=$USE_PROXY
PRIMERFORGE_SECURE_COOKIES=$USE_PROXY
PRIMERFORGE_THREADS=24
PRIMERFORGE_SPECIES=$SPECIES
PRIMERFORGE_SECRET_KEY=$SECRET
EOF
umask 022
chown root:"$SERVICE_USER" "$ENV_FILE"
chmod 0640 "$ENV_FILE"
ok "$ENV_FILE"

# Run a management command as the service account with the service's settings.
as_service() {
  # shellcheck source=/dev/null
  (set -a; . "$ENV_FILE"; set +a; cd "$APP_DIR/app" && \
   runuser -u "$SERVICE_USER" -- "$APP_DIR/venv/bin/python" -m primerforge.cli "$@")
}

# ---------------------------------------------------------------------------
# systemd services and the primerforge command
# ---------------------------------------------------------------------------
step "Installing the services"
PROTECT_HOME=true
case "$DATA_DIR" in /home/*|/root/*|/run/user/*) PROTECT_HOME=false ;; esac
SYSTEMCTL="$(command -v systemctl)"
render() {
  sed -e "s#@APP_DIR@#$APP_DIR#g" -e "s#@DATA_DIR@#$DATA_DIR#g" -e "s#@CONF_DIR@#$CONF_DIR#g" \
      -e "s#@USER@#$SERVICE_USER#g" -e "s#@PORT@#$PORT#g" -e "s#@PROTECT_HOME@#$PROTECT_HOME#g" \
      -e "s#@SYSTEMCTL@#$SYSTEMCTL#g" "$1"
}
render "$APP_DIR/installation/services/primerforge.service" > /etc/systemd/system/primerforge.service
render "$APP_DIR/installation/services/primerforge-setup.service" > /etc/systemd/system/primerforge-setup.service
render "$APP_DIR/installation/services/primerforge-backup.service" > /etc/systemd/system/primerforge-backup.service
render "$APP_DIR/installation/services/primerforge-backup.timer" > /etc/systemd/system/primerforge-backup.timer
render "$APP_DIR/installation/primerforge.sh" > /usr/local/bin/primerforge
chmod 0755 /usr/local/bin/primerforge
systemctl daemon-reload
systemctl enable -q --now primerforge-backup.timer
ok "primerforge.service, primerforge-setup.service, a nightly database backup and the 'primerforge' command"

# ---------------------------------------------------------------------------
# Administrator account
# ---------------------------------------------------------------------------
step "Setting up accounts"
chown -R "$SERVICE_USER":"$SERVICE_USER" "$DATA_DIR"
ADMIN_PASSWORD=""
if as_service users has-admin; then
  ok "an administrator account already exists (kept)"
else
  ADMIN_PASSWORD="$("$PY" -c 'import secrets, string; a = string.ascii_letters + string.digits; print("".join(secrets.choice(a) for _ in range(20)))')"
  printf '%s\n' "$ADMIN_PASSWORD" | as_service users add "$ADMIN_NAME" --admin --password-stdin >/dev/null
  umask 077
  printf 'PrimerForge administrator\nusername: %s\npassword: %s\n' "$ADMIN_NAME" "$ADMIN_PASSWORD" > "$CONF_DIR/initial-admin-password"
  umask 022
  chmod 0600 "$CONF_DIR/initial-admin-password"
  ok "created administrator '$ADMIN_NAME' (password saved in $CONF_DIR/initial-admin-password)"
fi

# ---------------------------------------------------------------------------
# HTTPS proxy (Caddy)
# ---------------------------------------------------------------------------
if (( USE_PROXY )); then
  step "Setting up HTTPS"
  CADDY_BIN="$APP_DIR/caddy/caddy"
  if [[ ! -x "$CADDY_BIN" ]]; then
    TMP="$(mktemp -d)"
    if CADDY_TGZ="$(bundled caddy 'caddy_*_linux_amd64.tar.gz')"; then
      CADDY_VERSION="$(basename "$CADDY_TGZ" | sed 's/caddy_\(.*\)_linux_amd64.tar.gz/\1/')"
      cp "$CADDY_TGZ" "$TMP/caddy.tar.gz"
      SOURCE="from the bundled dependencies"
    else
      CADDY_VERSION="$(curl -fsSL --max-time 20 https://api.github.com/repos/caddyserver/caddy/releases/latest 2>/dev/null \
                       | sed -n 's/.*"tag_name": *"v\([^"]*\)".*/\1/p' | head -1 || true)"
      CADDY_VERSION="${CADDY_VERSION:-$CADDY_FALLBACK_VERSION}"
      curl -fsSL --retry 3 -o "$TMP/caddy.tar.gz" \
        "https://github.com/caddyserver/caddy/releases/download/v${CADDY_VERSION}/caddy_${CADDY_VERSION}_linux_amd64.tar.gz" \
        || die "could not download Caddy ${CADDY_VERSION} from GitHub (bundle it with prepare.py)."
      SOURCE="downloaded"
    fi
    tar -C "$TMP" -xzf "$TMP/caddy.tar.gz" caddy
    install -D -m 0755 "$TMP/caddy" "$CADDY_BIN"
    rm -rf "$TMP"
    ok "Caddy $CADDY_VERSION $SOURCE"
  else
    ok "Caddy already installed"
  fi
  if ! id primerforge-caddy >/dev/null 2>&1; then
    useradd --system --user-group --home-dir /var/lib/primerforge-caddy --no-create-home \
            --shell "$(command -v nologin || echo /bin/false)" primerforge-caddy
  fi
  install -d -m 0750 -o primerforge-caddy -g primerforge-caddy /var/lib/primerforge-caddy

  GLOBAL=""
  if (( SELF_SIGNED )); then
    SITES="https://localhost"
    for name in "$(host_name)" "$(hostname 2>/dev/null || true)"; do
      if [[ -n "$name" && "$SITES" != *"https://$name,"* && "$SITES" != *"https://$name" ]]; then
        SITES="$SITES, https://$name"
      fi
    done
    for ip in $(host_ips); do
      if [[ "$ip" == *:* ]]; then SITES="$SITES, https://[$ip]"; else SITES="$SITES, https://$ip"; fi
    done
    TLS_LINE="    tls internal"
    PUBLIC_URL="https://$(first_ip)"
  else
    SITES="$DOMAIN"
    TLS_LINE=""
    if [[ -n "$EMAIL" ]]; then GLOBAL="{
    email $EMAIL
}"; fi
    PUBLIC_URL="https://$DOMAIN"
  fi
  cat > "$CONF_DIR/Caddyfile" <<EOF
# HTTPS in front of PrimerForge. Apply changes with: sudo systemctl restart primerforge-proxy
$GLOBAL

$SITES {
$TLS_LINE
    encode gzip
    request_body {
        max_size 512MB
    }
    reverse_proxy 127.0.0.1:$PORT
}
EOF
  chmod 0644 "$CONF_DIR/Caddyfile"
  "$CADDY_BIN" fmt --overwrite "$CONF_DIR/Caddyfile" >/dev/null 2>&1 || true
  "$CADDY_BIN" validate --config "$CONF_DIR/Caddyfile" --adapter caddyfile >/dev/null 2>&1 \
    || die "the generated Caddyfile is not valid; see $CONF_DIR/Caddyfile"

  render "$APP_DIR/installation/services/primerforge-proxy.service" > /etc/systemd/system/primerforge-proxy.service
  systemctl daemon-reload
  if command -v ss >/dev/null && ss -ltnH '( sport = :80 or sport = :443 )' 2>/dev/null | grep -q . \
     && ! systemctl is-active --quiet primerforge-proxy; then
    warn "another web server already uses port 80 or 443, so the HTTPS proxy was not started."
    warn "Point that server at http://127.0.0.1:$PORT, or stop it and run: sudo systemctl enable --now primerforge-proxy"
  else
    systemctl enable -q primerforge-proxy
    systemctl restart primerforge-proxy
    ok "HTTPS proxy running (primerforge-proxy.service)"
  fi
else
  # Switching from HTTPS back to plain HTTP: retire the proxy.
  if [[ -f /etc/systemd/system/primerforge-proxy.service ]]; then
    systemctl disable -q --now primerforge-proxy 2>/dev/null || true
    rm -f /etc/systemd/system/primerforge-proxy.service
    systemctl daemon-reload
  fi
  PUBLIC_URL="http://$(first_ip):$PORT"
fi

# ---------------------------------------------------------------------------
# Firewall
# ---------------------------------------------------------------------------
if (( OPEN_FIREWALL )); then
  if (( USE_PROXY )); then PORTS="80 443"; else PORTS="$PORT"; fi
  if command -v ufw >/dev/null && ufw status 2>/dev/null | grep -q "Status: active"; then
    for p in $PORTS; do ufw allow "$p/tcp" >/dev/null; done
    ok "opened port(s) $PORTS in ufw"
  elif command -v firewall-cmd >/dev/null && firewall-cmd --state >/dev/null 2>&1; then
    for p in $PORTS; do firewall-cmd -q --permanent --add-port="$p/tcp"; done
    firewall-cmd -q --reload
    ok "opened port(s) $PORTS in firewalld"
  fi
fi

# ---------------------------------------------------------------------------
# Start
# ---------------------------------------------------------------------------
step "Starting PrimerForge"
systemctl enable -q primerforge
systemctl restart primerforge
for _ in $(seq 1 60); do
  if curl -fsS --max-time 3 "http://127.0.0.1:$PORT/api/health" >/dev/null 2>&1; then break; fi
  sleep 1
done
if ! curl -fsS --max-time 3 "http://127.0.0.1:$PORT/api/health" >/dev/null 2>&1; then
  journalctl -u primerforge -n 30 --no-pager >&2 || true
  die "the service did not answer on port $PORT (log above)."
fi
ok "running (primerforge.service)"

if (( WITH_DATA )); then
  if systemctl is-active --quiet primerforge-setup; then
    ok "data download already in progress"
  else
    systemctl start --no-block primerforge-setup
    ok "started downloading data in the background (primerforge-setup.service)"
  fi
fi

trap - ERR
cat <<EOF

${G}${B}PrimerForge is installed.${N}

  Open:      ${B}$PUBLIC_URL${N}
EOF
if [[ -n "$ADMIN_PASSWORD" ]]; then
cat <<EOF
  Sign in:   ${B}$ADMIN_NAME${N} / ${B}$ADMIN_PASSWORD${N}
             (also saved in $CONF_DIR/initial-admin-password; change it under your name at top right)
EOF
fi
if (( SELF_SIGNED )); then
  echo "  Note:      the certificate is self-signed, so browsers ask you to confirm the first time."
fi
cat <<EOF

  Add people from the Users page, or: sudo primerforge users add NAME
EOF
if (( WITH_DATA )); then
cat <<EOF

  Genome data is downloading in the background: about an hour on a fast
  connection, several hours on a slow one.
  Primer design works right away using Ensembl's web services, and BLAST/BLAT
  become available as each database finishes.
    Progress:  sudo primerforge data-progress     (or Settings in the web app)
EOF
else
cat <<EOF

  No genome data was downloaded (--no-data). Start it with: sudo primerforge setup-data
EOF
fi
cat <<EOF

  Status:    sudo primerforge status
  Logs:      sudo primerforge logs
  Upgrade:   copy a newer installation folder and run its install.sh again (data and users are kept)
  Remove:    sudo primerforge uninstall
EOF
