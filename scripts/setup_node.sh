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
apt-get update -q
apt-get install -y -q docker.io curl git rsync ca-certificates
apt-get install -y -q openjdk-17-jdk-headless || apt-get install -y -q default-jdk-headless   # javac --release 17
systemctl enable --now docker

# upstream repo: scoring formula and the validator's templates
[ -d /minos/minos_subnet ] || git clone -q https://github.com/minos-protocol/minos_subnet.git /minos/minos_subnet

# python env: uv installs Python $PY (whatever the OS ships) and the exact pinned requirements
command -v uv >/dev/null || curl -LsSf https://astral.sh/uv/install.sh | env UV_INSTALL_DIR=/usr/local/bin sh
uv venv -q --python "$PY" --allow-existing .venv
uv pip install -q --python .venv/bin/python -r requirements.txt

# images: the validator's GATK, and hap.py at the exact digest upstream pins (pulled now: pulling it on first use
# would stall the node's first scoring for minutes)
docker pull -q broadinstitute/gatk:4.5.0.0
HAPPY=$(grep -oP 'HAPPY_DOCKER_IMAGE = "\K[^"]+' /minos/minos_subnet/utils/scoring.py)
docker pull -q "$HAPPY"

# the GATK batching patch, compiled against the exact jar inside the validator image
mkdir -p /opt/gatk
if [ ! -f /opt/gatk/gatk-package-4.5.0.0-local.jar ]; then
  cid=$(docker create broadinstitute/gatk:4.5.0.0)
  docker cp "$cid:/gatk/gatk-package-4.5.0.0-local.jar" /opt/gatk/ && docker rm "$cid" >/dev/null
fi
bash evaluator/gatk_patch/build.sh

# reference genomes (+ SDF for hap.py): `control provision` pushes main's copies first; this only downloads what is
# still missing (fetch_genome.sh skips files already present)
DATA_DIR="$DATA_DIR" bash scripts/fetch_genome.sh

echo "node ready. workers: $(.venv/bin/python -c 'from common.settings import default_workers; print(default_workers())')"
echo "done: main starts the node agent over ssh (python -m main ... --fleet)"
