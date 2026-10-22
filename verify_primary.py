#!/usr/bin/env python3
"""C42-only command wrapper for the unchanged public interval checker."""
import argparse, hashlib, json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "sources/roehrig/scripts"))
from certify_turan import certify, iv

def main():
    p = argparse.ArgumentParser()
    p.add_argument("certificate", type=Path)
    p.add_argument("--dps", type=int, default=100)
    p.add_argument("--terms", type=int, default=240)
    p.add_argument("--output", type=Path)
    a = p.parse_args()
    if a.dps < 40 or a.terms < 40:
        p.error("at least 40 decimal places and terms required")
    iv.dps = a.dps
    raw = a.certificate.read_bytes()
    try:
        result = certify(json.loads(raw), a.terms)
        result.update(verdict="PASS", certificate_sha256=hashlib.sha256(raw).hexdigest())
        code = 0
    except Exception as exc:
        result = {"verdict": "FAIL", "reason": str(exc)}
        code = 1
    output = json.dumps(result, indent=2) + "\n"
    if a.output: a.output.write_text(output)
    print(output, end="")
    return code

if __name__ == "__main__": sys.exit(main())
