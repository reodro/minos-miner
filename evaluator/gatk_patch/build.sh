#!/usr/bin/env bash
# Build the patched HaplotypeCallerEngine against the exact GATK 4.5.0.0 jar shipped in broadinstitute/gatk:4.5.0.0.
# Output: build/minos-patch.jar, placed BEFORE the GATK jar on the classpath at run time.
set -euo pipefail
cd "$(dirname "$0")"
JAR="${GATK_JAR:-/opt/gatk/gatk-package-4.5.0.0-local.jar}"
rm -rf build && mkdir -p build/classes
javac --release 17 -nowarn -cp "$JAR" -d build/classes $(find src -name '*.java')
jar cf build/minos-patch.jar -C build/classes .
echo "built $(ls build/classes/org/broadinstitute/hellbender/tools/walkers/haplotypecaller | wc -l) class files -> build/minos-patch.jar"
