#!/usr/bin/env bash
# Turn a fresh Ubuntu sub-VPS into a fleet node (GATK + scoring). Normally run by main, not by hand:
#   python -m evaluator.fleet.control provision <host>     (copies the repo, runs this script, pushes the self-test)
# Main then starts the node agent itself over ssh whenever it mines with --fleet.
set -euo pipefail
cd "$(dirname "$0")/.."
REPO="$(pwd)"
DATA_DIR="${DATA_DIR:-/minos/data}"

apt-get update -q
DEBIAN_FRONTEND=noninteractive apt-get install -y -q docker.io openjdk-17-jdk-headless python3-venv python3-pip \
  bcftools curl git rsync
systemctl enable --now docker

# upstream repo: scoring formula and the validator's templates
[ -d /minos/minos_subnet ] || git clone -q https://github.com/minos-protocol/minos_subnet.git /minos/minos_subnet

# python env
python3 -m venv .venv && .venv/bin/pip install -q --upgrade pip && .venv/bin/pip install -q -r requirements.txt

# images: the validator's GATK, samtools; hap.py (digest pinned upstream) is pulled on first use
docker pull -q broadinstitute/gatk:4.5.0.0
docker pull -q quay.io/biocontainers/samtools:1.20--h50ea8bc_0

# the GATK batching patch, compiled against the exact jar inside the validator image
mkdir -p /opt/gatk
if [ ! -f /opt/gatk/gatk-package-4.5.0.0-local.jar ]; then
  cid=$(docker create broadinstitute/gatk:4.5.0.0)
  docker cp "$cid:/gatk/gatk-package-4.5.0.0-local.jar" /opt/gatk/ && docker rm "$cid" >/dev/null
fi
bash evaluator/gatk_patch/build.sh

# reference genomes (+ SDF for hap.py)
DATA_DIR="$DATA_DIR" bash scripts/fetch_genome.sh

echo "node ready. workers: $(.venv/bin/python -c 'from common.settings import default_workers; print(default_workers())')"
echo "done: main starts the node agent over ssh (python -m main ... --fleet)"
