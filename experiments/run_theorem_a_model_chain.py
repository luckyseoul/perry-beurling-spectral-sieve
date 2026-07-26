#!/usr/bin/env python3
"""Multi-T Conditional Theorem A model-chain snapshot (not arithmetic A, not RH)."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pbss.theorem_a_chain import multi_T_model_chain, package_status  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=str, default="results/theorem_a_model_chain")
    ap.add_argument("--T-list", type=str, default="12,18,24,36,48")
    ap.add_argument("--n-zeros", type=int, default=10)
    ap.add_argument("--degree", type=int, default=4)
    ap.add_argument("--n-points", type=int, default=4096)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    print("NOT AN UNCONDITIONAL PROOF OF RH", flush=True)
    T_values = [float(x) for x in args.T_list.split(",") if x.strip()]
    rows = multi_T_model_chain(
        T_values,
        degree=args.degree,
        n_zeros=args.n_zeros,
        n_points=args.n_points,
    )
    summary = {
        "campaign": "theorem_a_model_chain",
        "package_status": package_status(),
        "n_rows": len(rows),
        "T_values": T_values,
        "rows": rows,
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
        "note": (
            "Model chain exercise for Conditional Theorem A package. "
            "Full arithmetic A open. RH open."
        ),
    }
    (out / "model_chain.json").write_text(json.dumps(summary, indent=2))
    lines = [
        f"n_rows={len(rows)}",
        f"package={package_status()}",
    ]
    for r in rows:
        e = r["empirical"]
        lines.append(
            f"T={r['T']:.1f} R_cl={e['R_d_cl']:.6g} R_cl_w={e['R_d_cl_weighted']:.6g} "
            f"R_ef={e['R_d_ef']:.6g} ok={r['proved_model_decay_ok']}"
        )
    lines.append("NOT AN UNCONDITIONAL PROOF OF RH")
    (out / "model_chain.txt").write_text("\n".join(lines) + "\n")
    (out / "MODEL_CHAIN_COMPLETE").write_text(
        json.dumps(
            {
                "status": "completed",
                "n_rows": len(rows),
                "full_arithmetic_A": "open",
                "rh": "open",
                "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps({"status": "completed", "n_rows": len(rows)}, indent=2))


if __name__ == "__main__":
    main()
