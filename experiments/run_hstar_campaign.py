#!/usr/bin/env python3
"""Corrected model-cancellation audit; injected bumps do not prove arithmetic coercivity."""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.hstar_anticancellation import (  # noqa: E402
    hstar_anticancellation_report,
    write_campaign_artifacts,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', type=Path, default=ROOT / 'results' / 'hstar_corrected')
    parser.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 4) - 2))
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('--workers must be positive')
    out, workers = args.out_dir, args.workers
    paths = write_campaign_artifacts(out, max_workers=workers, make_plot=True)
    print("PATHS", json.dumps(paths, indent=2))
    rep = hstar_anticancellation_report(max_workers=workers)
    print("status", rep["candidate_status"])
    print("analytic", json.dumps(rep.get("analytic_envelope"), indent=2, default=str)[:800])
    print("bridge", json.dumps(rep.get("arithmetic_bridge"), indent=2, default=str)[:800])
    print("rh_claimed", rep["rh_claimed"], "b_res_solved", rep["b_res_solved"])
    print((out / "SUMMARY.txt").read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
