#!/bin/sh
set -eu
cd "$(dirname "$0")"
PYTHON=${PYTHON:-python3}
mkdir -p reproduction-output
"$PYTHON" verify_primary.py certificates/c42_256.json --dps 100 --terms 240 --output reproduction-output/interval.json
"$PYTHON" independent/check_bounds.py certificates/c42_256.json --bits 768 --terms 400 --out reproduction-output/arb.json
"$PYTHON" rejection_controls.py > reproduction-output/rejection_controls.log
cat reproduction-output/rejection_controls.log
printf '%s\n' 'PASS: both C42 checkers and all eight rejection controls.'
