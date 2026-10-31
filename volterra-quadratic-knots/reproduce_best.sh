#!/bin/bash
# Run from this folder. Up to about 95 GB of RAM. About 40 minutes on 32 vCPUs.
# Set PYTHON to change the interpreter. Needs python-flint 0.9.0, numpy, scipy, mpmath.
set -e
cd "$(dirname "$0")"
PYTHON="${PYTHON:-python3}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
C42_REPRO_DIR="logs/reproduce_$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$C42_REPRO_DIR"
nice -n 10 "$PYTHON" verifier_versions/verify_volterra_7134ee6184e563f2.py certificates/C42_compact_best.json --N 14155776 --bits 224 --threads 3 --weighted-residual --verbose > "$C42_REPRO_DIR/primary.json" 2> "$C42_REPRO_DIR/primary.progress"
nice -n 10 "$PYTHON" verifier_versions/check_bakeas_local_3d07e55ffcc2e49a.py certificates/C42_compact_best.json 294912 "$C42_REPRO_DIR/public.json" > "$C42_REPRO_DIR/public.log" 2>&1
nice -n 10 "$PYTHON" - "$C42_REPRO_DIR/public.json" <<'PY'
import json,sys
sys.exit(0 if json.load(open(sys.argv[1]))['verdict']=='PASS' else 1)
PY
printf '%s\n' "$C42_REPRO_DIR"
