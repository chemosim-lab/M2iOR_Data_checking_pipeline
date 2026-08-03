#!/usr/bin/env bash
# Installs the blastp binary (NCBI BLAST+) into .venv/bin/, so
# scripts/find_protein_mutations.py can shell out to it for identity/mutation
# calculations that match blast.ncbi.nlm.nih.gov exactly. No distro package
# is used (Fedora ships none) - this fetches NCBI's official prebuilt
# binaries directly, which needs no sudo.
#
# Re-run this after any `uv sync` that recreates .venv from scratch, since
# .venv/bin isn't version-controlled.
set -euo pipefail

BLAST_VERSION="2.17.0+"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_BIN="$REPO_ROOT/.venv/bin"
TARBALL="ncbi-blast-${BLAST_VERSION}-x64-linux.tar.gz"
URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/${TARBALL}"

if [ ! -d "$VENV_BIN" ]; then
    echo "error: $VENV_BIN not found - run 'uv sync' first" >&2
    exit 1
fi

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

echo "Downloading $URL ..."
curl -sL --fail -o "$WORKDIR/$TARBALL" "$URL"
curl -sL --fail -o "$WORKDIR/$TARBALL.md5" "$URL.md5"

(cd "$WORKDIR" && md5sum -c "$TARBALL.md5")

tar -xzf "$WORKDIR/$TARBALL" -C "$WORKDIR" "ncbi-blast-${BLAST_VERSION}/bin/blastp"
install -m 0755 "$WORKDIR/ncbi-blast-${BLAST_VERSION}/bin/blastp" "$VENV_BIN/blastp"

echo "Installed: $("$VENV_BIN/blastp" -version | head -1)"
