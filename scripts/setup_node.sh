#!/usr/bin/env bash
# Turn a fresh Ubuntu sub-VPS into a fleet node (GATK + scoring). Normally run by main, not by hand:
#   python -m evaluator.fleet.control provision <host>     (copies the repo, runs this script, pushes the self-test)
# Main then starts the node agent itself over ssh whenever it mines with --fleet.
set -euo pipefail
cd "$(dirname "$0")/.."
REPO="$(pwd)"
DATA_DIR="${DATA_DIR:-/minos/data}"

PY="${PY_VERSION:-3.10}"          # main's Python: every node runs the same interpreter and pinned libraries
export DEBIAN_FRONTEND=noninteractive
# a fresh VPS often runs unattended-upgrades for minutes after boot: wait for its apt lock instead of failing
APT="apt-get -o DPkg::Lock::Timeout=900"
$APT update -q
$APT install -y -q curl git rsync ca-certificates
# Docker: Ubuntu's docker.io, unless the image already ships one (e.g. vast.ai VMs come with Docker CE, whose
# containerd.io conflicts with docker.io's containerd)
command -v docker >/dev/null || $APT install -y -q docker.io
$APT install -y -q openjdk-17-jdk-headless || $APT install -y -q default-jdk-headless   # javac --release 17
systemctl enable --now docker

# upstream repo (scoring formula, the validator's templates): `control provision/sync` copies main's exact version;
# cloned only when run by hand on a machine that has none
[ -d /minos/minos_subnet ] || git clone -q https://github.com/minos-protocol/minos_subnet.git /minos/minos_subnet

# python env: uv installs Python $PY (whatever the OS ships) and the exact pinned requirements
command -v uv >/dev/null || curl -LsSf https://astral.sh/uv/install.sh | env UV_INSTALL_DIR=/usr/local/bin sh
uv venv -q --python "$PY" --allow-existing .venv
uv pip install -q --python .venv/bin/python -r requirements.txt

# images: the validator's GATK, and the two images upstream's scoring runs - hap.py at its pinned digest and bcftools
# (pulled now: pulling on first use would stall, or on a host that cannot reach the registry fail, the first scorings)
docker pull -q broadinstitute/gatk:4.5.0.0
HAPPY=$(grep -oP 'HAPPY_DOCKER_IMAGE = "\K[^"]+' /minos/minos_subnet/utils/scoring.py)
BCFTOOLS=$(grep -oP 'BCFTOOLS_DOCKER_IMAGE = "\K[^"]+' /minos/minos_subnet/utils/scoring.py)
docker pull -q "$HAPPY"
docker pull -q "$BCFTOOLS"

# the GATK batching patch, compiled against the exact jar inside the validator image
mkdir -p /opt/gatk
if [ ! -f /opt/gatk/gatk-package-4.5.0.0-local.jar ]; then
  cid=$(docker create broadinstitute/gatk:4.5.0.0)
  docker cp "$cid:/gatk/gatk-package-4.5.0.0-local.jar" /opt/gatk/ && docker rm "$cid" >/dev/null
fi
bash evaluator/gatk_patch/build.sh

# reference genomes (+ SDF for hap.py) from the platform's endpoint: chr14-22 (every round so far); a node fetches
# any other chromosome on demand. Only what is still missing is downloaded (`control genomes` can push main's copies)
DATA_DIR="$DATA_DIR" bash scripts/fetch_genome.sh

echo "node ready. workers: $(.venv/bin/python -c 'from common.settings import default_workers; print(default_workers())')"
echo "done: main starts the node agent over ssh (python -m main ... --fleet)"
