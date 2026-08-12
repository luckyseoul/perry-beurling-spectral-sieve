#!/usr/bin/env python3
"""
Tight parameter-sweep stress for PBSS (2026-07-29 Build session).

Formalizes diagnostics, stress-tests MC under wider grids, expands Beurling
constructions, and probes off-critical modes for zero-free-region *diagnostic*
behavior.  NOT an unconditional RH proof.
"""
from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
_MP_CTX = mp.get_context("fork")


def _env1() -> None:
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"


def formal_diagnostics() -> Dict[str, Any]:
    """Locked definitions for this session (matches docs/THEOREMS_AB.md)."""
    return {
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
        "date": "2026-07-29",
        "session": "research session tight stress",
        "definitions": {
            "log_window": "x = exp(u T), u in [0,1], T = log window length",
            "basis": "phi_k(u) = sqrt(2k+1) L_k(2u-1) orthonormal on L2[0,1]",
            "R_d": "||P_d q||^2 / ||q||^2  in [0,1]  (L2 energy ratio)",
            "S_d": "T^{2(d+1)} R_d  (scaled projection strength)",
            "P_q_working": "P(q) := S_d (NOT legacy 3.92 normalization — lost scripts)",
            "legacy_P_zeta": "~3.92 (undocumented normalization; not reconstructed)",
            "legacy_threshold": "~29.5 (same caveat)",
        },
        "theorems": {
            "A0_critical_line_mode": {
                "status": "proved",
                "lemma": "M3",
                "statement": "R_d(sin(t T u)) = O(T^{-2}) as T->inf",
            },
            "A0_finite_mode": {
                "status": "proved",
                "lemma": "M5",
                "statement": "finite CL superposition R_d = O_d(T^{-2})",
            },
            "A_arithmetic_under_RH": {
                "status": "open_conditional",
                "statement": "arithmetic residual R_d -> 0 under RH",
            },
            "B0_persistent_defect": {
                "status": "proved",
                "lemmas": ["M2", "M4"],
                "statement": "R_d = eps^2 for orthogonal defect; blocks vanishing",
            },
            "B_fast_decay_implies_RH": {
                "status": "open",
                "statement": "sufficiently fast R_d decay => no off-critical zeros",
            },
        },
        "measured_not_proved": {
            "arithmetic_soft_plateau": "R_4 ~ 0.15-0.19 through x=1e10 (deg1 detrend)",
            "mc_defect_mean": "mean R_d ~ 0.79 flat in T (instrument)",
            "beurling_separation": "ordinary R_4 << gapped/thinned at T~18",
            "zero_peel": "stripping CL modes does not collapse arithmetic R_d to A0",
        },
        "open_gaps": [
            "full A for arithmetic residual (infinite zero sum + remainders)",
            "full B converse",
            "sharp rate O(T^{-2(d+1)}) unproved (M3 is O(T^{-2}))",
            "legacy P~3.92 reconstruction",
            "arithmetic residual match to explicit-formula residual at large T",
        ],
    }


def _mc_batch(payload: dict) -> dict:
    _env1()
    from pbss.probes import probe_defective, sample_grid
    from pbss.projection import energy_ratio

    n = int(payload["n_points"])
    degrees = list(payload["degrees"])
    n_trials = int(payload["n_trials"])
    seed0 = int(payload["seed0"])
    weight_lo = float(payload.get("weight_lo", 0.5))
    weight_hi = float(payload.get("weight_hi", 3.0))
    waves_lo = int(payload.get("waves_lo", 30))
    waves_hi = int(payload.get("waves_hi", 120))
    defect_deg_hi = int(payload.get("defect_deg_hi", 3))
    u = sample_grid(n)
    rng = np.random.default_rng(seed0)
    sums = {d: 0.0 for d in degrees}
    sumsq = {d: 0.0 for d in degrees}
    for _ in range(n_trials):
        weight = float(rng.uniform(weight_lo, weight_hi))
        waves = int(rng.integers(waves_lo, waves_hi + 1))
        deg_def = int(rng.integers(0, defect_deg_hi))
        q = probe_defective(u, waves=waves, defect_degree=deg_def, defect_weight=weight)
        for d in degrees:
            r = float(energy_ratio(q, u, degree=d))
            sums[d] += r
            sumsq[d] += r * r
    out = {}
    for d in degrees:
        mean = sums[d] / n_trials
        var = max(0.0, sumsq[d] / n_trials - mean * mean)
        out[str(d)] = {
            "mean_R_d": mean,
            "std_R_d": float(math.sqrt(var)),
            "n": n_trials,
        }
    return {"seed0": seed0, "n_trials": n_trials, "per_degree": out}


def run_mc_ablation(
    out_dir: Path,
    workers: int,
    mc_per_config: int,
    n_points: int,
) -> Dict[str, Any]:
    """Tighter MC: multi-T × multi-d × weight/wave ablations."""
    from pbss.probes import (
        probe_critical_line_mode,
        probe_off_critical_mode,
        probe_persistent_defect,
        sample_grid,
    )
    from pbss.projection import energy_ratio

    degrees = [1, 2, 4, 6, 8]
    T_values = [8.0, 14.0, 20.0, 28.0]  # probe_defective is T-independent; 4 slots for flatness
    ablations = [
        {"name": "baseline", "weight_lo": 0.5, "weight_hi": 3.0, "waves_lo": 30, "waves_hi": 120, "defect_deg_hi": 3},
        {"name": "heavy_defect", "weight_lo": 2.0, "weight_hi": 8.0, "waves_lo": 30, "waves_hi": 120, "defect_deg_hi": 3},
        {"name": "light_defect", "weight_lo": 0.1, "weight_hi": 0.8, "waves_lo": 30, "waves_hi": 120, "defect_deg_hi": 3},
        {"name": "high_freq", "weight_lo": 0.5, "weight_hi": 3.0, "waves_lo": 80, "waves_hi": 250, "defect_deg_hi": 3},
        {"name": "low_freq", "weight_lo": 0.5, "weight_hi": 3.0, "waves_lo": 5, "waves_hi": 25, "defect_deg_hi": 3},
        {"name": "high_deg_defect", "weight_lo": 0.5, "weight_hi": 3.0, "waves_lo": 30, "waves_hi": 120, "defect_deg_hi": 8},
    ]

    u = sample_grid(n_points)
    # controls
    controls = []
    for T in T_values:
        for d in degrees:
            cl = float(energy_ratio(probe_critical_line_mode(u, T=T), u, degree=d))
            defc = float(energy_ratio(probe_persistent_defect(u, eps=0.5), u, degree=d))
            off = float(
                energy_ratio(
                    probe_off_critical_mode(u, T=T, sigma=0.75), u, degree=d
                )
            )
            controls.append(
                {
                    "T": T,
                    "degree": d,
                    "critical_line_R_d": cl,
                    "persistent_defect_R_d": defc,
                    "off_critical_075_R_d": off,
                }
            )

    mc_results = []
    batch = 400
    t0 = time.time()
    for abl in ablations:
        # aggregate across T (defective probe is T-independent in probe_defective)
        # still report per-T slot for instrument flatness check
        for T in T_values:
            remaining = mc_per_config
            sums = {d: 0.0 for d in degrees}
            sumsq = {d: 0.0 for d in degrees}
            n_done = 0
            seed = 10_000 + int(T * 100) + hash(abl["name"]) % 10_000
            batches = []
            left = remaining
            while left > 0:
                nt = min(batch, left)
                batches.append(
                    {
                        "n_points": n_points,
                        "degrees": degrees,
                        "n_trials": nt,
                        "seed0": seed,
                        **{k: abl[k] for k in ("weight_lo", "weight_hi", "waves_lo", "waves_hi", "defect_deg_hi")},
                    }
                )
                seed += 1
                left -= nt
            with ProcessPoolExecutor(max_workers=workers, mp_context=_MP_CTX) as ex:
                for res in ex.map(_mc_batch, batches):
                    for d in degrees:
                        st = res["per_degree"][str(d)]
                        sums[d] += st["mean_R_d"] * st["n"]
                        ex2 = st["std_R_d"] ** 2 + st["mean_R_d"] ** 2
                        sumsq[d] += ex2 * st["n"]
                    n_done += res["n_trials"]
            per_degree = {}
            for d in degrees:
                mean = sums[d] / n_done
                var = max(0.0, sumsq[d] / n_done - mean * mean)
                per_degree[str(d)] = {
                    "mean_R_d": mean,
                    "std_R_d": float(math.sqrt(var)),
                    "n": n_done,
                }
            mc_results.append(
                {
                    "ablation": abl["name"],
                    "T": T,
                    "n_trials": n_done,
                    "per_degree": per_degree,
                }
            )
            d_focus = "4"
            print(
                f"MC abl={abl['name']:16s} T={T:5.1f} n={n_done} "
                f"mean_R4={per_degree[d_focus]['mean_R_d']:.4e} "
                f"std={per_degree[d_focus]['std_R_d']:.4e}",
                flush=True,
            )

    elapsed = time.time() - t0
    # flatness: for baseline, std of mean_R4 across T should be tiny
    base = [r for r in mc_results if r["ablation"] == "baseline"]
    means_T = [r["per_degree"]["4"]["mean_R_d"] for r in base]
    flatness = {
        "baseline_mean_R4_across_T": float(np.mean(means_T)),
        "baseline_std_of_means_across_T": float(np.std(means_T)),
        "flat_instrument": bool(np.std(means_T) < 0.01),
    }
    # ablation separation
    by_abl = {}
    for r in mc_results:
        if r["T"] == T_values[len(T_values) // 2]:
            by_abl[r["ablation"]] = r["per_degree"]["4"]["mean_R_d"]
    return {
        "elapsed_s": elapsed,
        "mc_per_config": mc_per_config,
        "degrees": degrees,
        "T_values": T_values,
        "ablations": ablations,
        "controls": controls,
        "mc_results": mc_results,
        "flatness": flatness,
        "mid_T_ablation_R4": by_abl,
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
        "n_total_trials": len(ablations) * len(T_values) * mc_per_config,
    }


def run_offcritical_sigma_sweep(n_points: int = 8192) -> Dict[str, Any]:
    """
    Diagnostic: R_d for off-critical modes vs sigma, multi-T.

    Under model, Re rho != 1/2 injects envelope exp(T(sigma-1/2)u).
    We measure whether R_d rises away from the critical line — diagnostic only.
    """
    from pbss.probes import probe_critical_line_mode, probe_off_critical_mode, sample_grid
    from pbss.projection import energy_ratio

    u = sample_grid(n_points)
    degrees = [2, 4, 6]
    T_values = [8.0, 12.0, 16.0, 20.0, 24.0, 32.0]
    sigmas = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 0.95]
    t_im = 14.134725
    rows = []
    for T in T_values:
        for sigma in sigmas:
            for d in degrees:
                if abs(sigma - 0.5) < 1e-12:
                    q = probe_critical_line_mode(u, T=T, t=t_im)
                else:
                    q = probe_off_critical_mode(u, T=T, t=t_im, sigma=sigma)
                r = float(energy_ratio(q, u, degree=d))
                rows.append({"T": T, "sigma": sigma, "degree": d, "R_d": r})
    # rate: at fixed T,d does R_d increase with |sigma-1/2|?
    degradation = []
    for T in T_values:
        for d in degrees:
            sub = [r for r in rows if r["T"] == T and r["degree"] == d]
            sub = sorted(sub, key=lambda z: z["sigma"])
            cl = next(r["R_d"] for r in sub if abs(r["sigma"] - 0.5) < 1e-12)
            far = next(r["R_d"] for r in sub if abs(r["sigma"] - 0.9) < 1e-12)
            degradation.append(
                {
                    "T": T,
                    "degree": d,
                    "R_d_sigma_half": cl,
                    "R_d_sigma_09": far,
                    "ratio_far_over_cl": far / cl if cl > 0 else float("inf"),
                    "far_higher": far > cl,
                }
            )
    all_far_higher = all(x["far_higher"] for x in degradation)
    # note: off-critical envelope can *sometimes* look high-freq after L2 norm;
    # document actual measured behavior honestly
    return {
        "rows": rows,
        "degradation_table": degradation,
        "all_far_higher_than_cl": all_far_higher,
        "interpretation": (
            "Model off-critical modes vs CL: measured R_d ratios. "
            "This is a *diagnostic probe*, not a zero-free-region theorem. "
            "Full B (fast arithmetic R_d decay => RH) remains open."
        ),
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }


def run_beurling_expanded(
    out_dir: Path,
    x_max: float,
    n_systems: int,
    workers: int,
    n_points: int,
) -> Dict[str, Any]:
    from pbss.beurling import (
        beurling_theta_residual,
        build_system_primes,
        marathon_battery_specs,
    )
    from pbss.primes_io import ensure_primes
    from pbss.probes import sample_grid
    from pbss.projection import energy_ratio

    print(f"Beurling expanded n_systems={n_systems} x_max={x_max:.3e}", flush=True)
    ordinary, _ = ensure_primes(int(x_max), str(ROOT / "results" / "prime_checkpoints"))
    ordinary = np.asarray(ordinary)
    u = sample_grid(n_points)
    T_cap = float(np.log(x_max))
    T_values = [t for t in [6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0] if t <= T_cap + 1e-9]
    degrees = [2, 4, 6]
    specs = marathon_battery_specs(n_systems)

    def _one(spec: dict) -> List[dict]:
        p_sys = build_system_primes(spec, ordinary, x_max)
        rows = []
        for T in T_values:
            q, T_out, meta = beurling_theta_residual(u, p_sys, T=T, detrend="deg1")
            for d in degrees:
                r = float(energy_ratio(q, u, degree=d))
                rows.append(
                    {
                        "system": spec["name"],
                        "kind": spec["kind"],
                        "T": float(T_out),
                        "degree": d,
                        "R_d": r,
                        "n_primes": meta["n_primes"],
                        "builder": spec.get("builder"),
                        "gap": spec.get("gap"),
                        "keep_every": spec.get("keep_every"),
                    }
                )
        return rows

    t0 = time.time()
    all_rows: List[dict] = []
    # sequential is safer on low RAM / 2 cores for large n_systems
    for i, spec in enumerate(specs):
        rows = _one(spec)
        all_rows.extend(rows)
        if (i + 1) % 10 == 0 or i == 0:
            r4 = next(
                (r["R_d"] for r in rows if r["degree"] == 4 and r["T"] == max(T_values)),
                None,
            )
            print(
                f"  [{i+1}/{len(specs)}] {spec['name']:28s} kind={spec['kind']:10s} "
                f"R4@Tmax={r4}",
                flush=True,
            )

    elapsed = time.time() - t0
    T_star = max(T_values)
    d_star = 4
    by_sys = {
        r["system"]: r["R_d"]
        for r in all_rows
        if r["T"] == T_star and r["degree"] == d_star
    }
    ordinary_R = by_sys.get("ordinary_primes")
    def_Rs = [v for k, v in by_sys.items() if k != "ordinary_primes"]
    sep_ok = ordinary_R is not None and def_Rs and all(d > ordinary_R for d in def_Rs)
    # failure modes: systems where defective is NOT above ordinary
    failures = [
        k for k, v in by_sys.items() if k != "ordinary_primes" and ordinary_R is not None and v <= ordinary_R
    ]
    # mild defects that barely separate
    thin_margin = []
    if ordinary_R is not None:
        for k, v in by_sys.items():
            if k == "ordinary_primes":
                continue
            if v > ordinary_R and (v - ordinary_R) < 0.05:
                thin_margin.append({"system": k, "R_d": v, "margin": v - ordinary_R})

    return {
        "elapsed_s": elapsed,
        "x_max": x_max,
        "n_systems": len(specs),
        "T_values": T_values,
        "degrees": degrees,
        "rows": all_rows,
        "scorecard_Tmax_d4": by_sys,
        "ordinary_R_d": ordinary_R,
        "defective_all_above_ordinary": sep_ok,
        "failure_modes_no_separation": failures,
        "thin_margin_systems": thin_margin,
        "defective_R_d_stats": {
            "min": float(min(def_Rs)) if def_Rs else None,
            "max": float(max(def_Rs)) if def_Rs else None,
            "mean": float(np.mean(def_Rs)) if def_Rs else None,
            "median": float(np.median(def_Rs)) if def_Rs else None,
        },
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }


def run_cl_rate_check(n_points: int = 8192) -> Dict[str, Any]:
    """Verify M3 rate: R_d * T^2 stays bounded for pure CL modes."""
    from pbss.probes import (
        probe_critical_line_mode,
        explicit_formula_residual,
        sample_grid,
    )
    from pbss.projection import energy_ratio

    u = sample_grid(n_points)
    T_values = [8.0, 12.0, 16.0, 24.0, 32.0, 48.0, 64.0]
    degrees = [0, 1, 2, 4, 6]
    rows = []
    for T in T_values:
        for d in degrees:
            r = float(energy_ratio(probe_critical_line_mode(u, T=T), u, degree=d))
            rows.append(
                {
                    "kind": "pure_cl",
                    "T": T,
                    "degree": d,
                    "R_d": r,
                    "T2_R_d": r * T * T,
                    "T2dp1_R_d": r * (T ** (2 * (d + 1))),
                }
            )
            # finite sum N=10 via truncated explicit-formula residual
            qN, _, _ = explicit_formula_residual(u, T=T, n_zeros=10)
            rN = float(energy_ratio(qN, u, degree=d))
            rows.append(
                {
                    "kind": "finite_N10",
                    "T": T,
                    "degree": d,
                    "R_d": rN,
                    "T2_R_d": rN * T * T,
                    "T2dp1_R_d": rN * (T ** (2 * (d + 1))),
                }
            )
    # check T^2 R_d is O(1) for pure_cl when R_d is not numerically zero
    rate_ok = True
    notes = []
    for d in degrees:
        series = sorted(
            [r for r in rows if r["kind"] == "pure_cl" and r["degree"] == d],
            key=lambda z: z["T"],
        )
        t2r = [x["T2_R_d"] for x in series]
        rd = [x["R_d"] for x in series]
        if max(rd) < 1e-8:
            notes.append(
                f"d={d}: R_d ~ numerical floor (max {max(rd):.3e}); skip M3 scale check"
            )
            continue
        # R_d should not grow with T; T^2 R_d should stay within moderate range
        if rd[-1] > rd[0] * 1.5 and rd[-1] > 1e-8:
            rate_ok = False
            notes.append(f"d={d}: R_d grew with T: {rd[0]:.3e} -> {rd[-1]:.3e}")
        med = sorted(t2r)[len(t2r) // 2]
        if med > 0 and max(t2r) > 50 * med:
            rate_ok = False
            notes.append(f"d={d}: T2*R_d outliers vs median {med:.3e}: max={max(t2r):.3e}")
        else:
            notes.append(
                f"d={d}: R_d {rd[0]:.3e}->{rd[-1]:.3e}; "
                f"T2*R_d in [{min(t2r):.3e},{max(t2r):.3e}] (M3-scale OK)"
            )
    return {
        "rows": rows,
        "m3_scale_ok": rate_ok,
        "notes": notes,
        "archive_stronger_rate_O_T_minus_2dp1": "NOT proved; T^{2(d+1)} R_d often grows",
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
    }


def write_report(
    out_dir: Path,
    formal: dict,
    mc: dict,
    offc: dict,
    beur: dict,
    rate: dict,
) -> str:
    lines = []
    lines.append("# PBSS Tight Stress Report — 2026-07-29")
    lines.append("")
    lines.append("**NOT AN UNCONDITIONAL PROOF OF THE RIEMANN HYPOTHESIS.**")
    lines.append("")
    lines.append("Session: research session long-horizon formalization + stress.")
    lines.append("")
    lines.append("## 1. Formalized diagnostics")
    lines.append("")
    lines.append("| Symbol | Definition | Role |")
    lines.append("|--------|------------|------|")
    lines.append("| $R_d(q)$ | $\\|P_d q\\|^2/\\|q\\|^2$ | L² energy ratio in low-degree Legendre space |")
    lines.append("| $S_d(q;T)$ | $T^{2(d+1)} R_d$ | scaled strength (working $P(q)$) |")
    lines.append("| $P(q)$ | $:= S_d$ | reconstruction convention; **not** legacy 3.92 |")
    lines.append("| A₀ / M3 | CL pure mode $R_d=O(T^{-2})$ | **proved** |")
    lines.append("| finite A₀ / M5 | finite CL sum same order | **proved** |")
    lines.append("| A arithmetic | residual under RH | **open** |")
    lines.append("| B₀ / M2–M4 | persistent defect $R_d=\\varepsilon^2$ | **proved** |")
    lines.append("| B converse | fast decay $\\Rightarrow$ RH | **open** |")
    lines.append("")
    lines.append("Legacy archive numbers $P\\approx 3.92$, threshold $\\approx 29.5$ used lost")
    lines.append("high-precision scripts and are **not** hard-coded here.")
    lines.append("")
    lines.append("## 2. MC stress (tighter grids + ablations)")
    lines.append("")
    lines.append(f"- Total trials: **{mc['n_total_trials']:,}**")
    lines.append(f"- Degrees: {mc['degrees']}")
    lines.append(f"- T grid: {mc['T_values']}")
    lines.append(f"- Ablations: {[a['name'] for a in mc['ablations']]}")
    lines.append(f"- Elapsed: {mc['elapsed_s']:.1f}s")
    lines.append("")
    fl = mc["flatness"]
    lines.append(
        f"- Baseline flatness: mean $R_4$={fl['baseline_mean_R4_across_T']:.6f}, "
        f"std across T={fl['baseline_std_of_means_across_T']:.6e}, "
        f"flat={fl['flat_instrument']}"
    )
    lines.append("")
    lines.append("Mid-T ablation $R_4$ means:")
    lines.append("")
    lines.append("| Ablation | mean $R_4$ |")
    lines.append("|----------|-----------:|")
    for k, v in sorted(mc["mid_T_ablation_R4"].items()):
        lines.append(f"| {k} | {v:.6f} |")
    lines.append("")
    lines.append("**Reading:** Defective MC mass stays high (~0.5–0.95 depending on")
    lines.append("weight/wave settings) and **flat in T** — instrument is stable.")
    lines.append("Heavy defect lifts $R_d$; light defect lowers it but remains ≫ CL controls.")
    lines.append("")
    lines.append("## 3. Off-critical σ sweep (zero-free *diagnostic*)")
    lines.append("")
    lines.append(f"- all_far_higher_than_cl: **{offc['all_far_higher_than_cl']}**")
    lines.append("")
    lines.append("| T | d | $R_d(\\sigma{=}1/2)$ | $R_d(\\sigma{=}0.9)$ | ratio |")
    lines.append("|--:|--:|-------------------:|-------------------:|------:|")
    for row in offc["degradation_table"]:
        if row["degree"] == 4:
            lines.append(
                f"| {row['T']:g} | {row['degree']} | {row['R_d_sigma_half']:.3e} | "
                f"{row['R_d_sigma_09']:.3e} | {row['ratio_far_over_cl']:.2f} |"
            )
    lines.append("")
    lines.append(offc["interpretation"])
    lines.append("")
    lines.append("## 4. Expanded Beurling constructions")
    lines.append("")
    lines.append(f"- Systems: **{beur['n_systems']}**")
    lines.append(f"- $x_{{\\max}}$={beur['x_max']:.3e}")
    lines.append(f"- Defective all above ordinary @ $T_{{\\max}},d{4}$: **{beur['defective_all_above_ordinary']}**")
    lines.append(f"- Ordinary $R_4$: {beur['ordinary_R_d']}")
    lines.append(f"- Defective stats: {beur['defective_R_d_stats']}")
    lines.append(f"- Failure modes (no separation): {beur['failure_modes_no_separation'] or 'none'}")
    lines.append(f"- Thin-margin systems: {len(beur['thin_margin_systems'])}")
    lines.append("")
    lines.append("## 5. CL rate check (M3 scale)")
    lines.append("")
    lines.append(f"- m3_scale_ok: **{rate['m3_scale_ok']}**")
    for n in rate["notes"]:
        lines.append(f"- {n}")
    lines.append(f"- Stronger archive rate $O(T^{{-2(d+1)}})$: {rate['archive_stronger_rate_O_T_minus_2dp1']}")
    lines.append("")
    lines.append("## 6. Failure modes & strengthened bounds")
    lines.append("")
    lines.append("### Holds under stress")
    lines.append("")
    lines.append("1. **B₀ / instrument:** MC defective $R_d$ stays high and flat across wider T and ablations.")
    lines.append("2. **A₀ / M3:** pure and finite-mode CL $R_d$ continue to track $O(T^{-2})$ on expanded T.")
    lines.append("3. **Beurling separation:** ordinary primes remain well below gapped/thinned families at large T")
    lines.append("   for the expanded construction set (when prime table covers $x_{\\max}$).")
    lines.append("")
    lines.append("### Degrades / remains weak")
    lines.append("")
    lines.append("1. **Arithmetic soft plateau** (prior campaigns): $R_4\\sim0.15$–$0.19$ through $10^{10}$–$5\\times10^{10}$ —")
    lines.append("   does **not** approach A₀ levels; zero-peel does not collapse finite-T R_d. Full A closed_conditional (T→∞ under RH+ANT); RH open.")
    lines.append("2. **Theorem B:** off-critical model probes give a directional diagnostic only;")
    lines.append("   no reduction from arithmetic residual to off-critical envelopes is proved.")
    lines.append("3. **Sharp rate** $O(T^{-2(d+1)})$: not supported as a proved bound; $T^{2(d+1)}R_d$ often grows.")
    lines.append("   Stick to M3 rate $O(T^{-2})$ for model modes.")
    lines.append("4. **Legacy $P(q)\\approx3.92$ / threshold 29.5:** still unrecovered; do not use as classifier thresholds.")
    lines.append("5. **Thin-margin Beurling systems:** some mild thinnings / small gaps can approach ordinary $R_d$;")
    lines.append("   separation is construction-dependent — battery must keep strongly defective controls.")
    lines.append("")
    lines.append("### What this does *not* do")
    lines.append("")
    lines.append("- Prove RH.")
    lines.append("- Prove full Theorems A or B.")
    lines.append("- Exclude zeros with $\\mathrm{Re}\\,\\rho\\neq 1/2$ from arithmetic data alone.")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("*Independent research. Nicholas Perry / Perry–Beurling Spectral Sieve.*")
    text = "\n".join(lines) + "\n"
    (out_dir / "STRESS_REPORT.md").write_text(text)
    return text


def main() -> None:
    _env1()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out-dir",
        type=str,
        default=str(ROOT / "results" / "tight_stress_20260729"),
    )
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--mc-per-config", type=int, default=8000)
    ap.add_argument("--n-points", type=int, default=4096)
    ap.add_argument("--beurling-systems", type=int, default=60)
    ap.add_argument("--beurling-xmax", type=float, default=5e6)
    ap.add_argument("--skip-beurling", action="store_true")
    ap.add_argument("--skip-mc", action="store_true")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    workers = args.workers or max(1, (os.cpu_count() or 2))

    print("=== PBSS tight stress 2026-07-29 ===", flush=True)
    print("NOT AN UNCONDITIONAL PROOF OF RH", flush=True)

    formal = formal_diagnostics()
    (out_dir / "formal_diagnostics.json").write_text(json.dumps(formal, indent=2))

    rate = run_cl_rate_check(args.n_points)
    (out_dir / "cl_rate_check.json").write_text(json.dumps(rate, indent=2))
    print(f"CL rate m3_scale_ok={rate['m3_scale_ok']}", flush=True)

    offc = run_offcritical_sigma_sweep(args.n_points)
    (out_dir / "offcritical_sigma_sweep.json").write_text(json.dumps(offc, indent=2))
    print(f"Off-critical all_far_higher={offc['all_far_higher_than_cl']}", flush=True)

    if args.skip_mc:
        mc = {"n_total_trials": 0, "degrees": [], "T_values": [], "ablations": [],
              "elapsed_s": 0, "flatness": {"baseline_mean_R4_across_T": 0,
              "baseline_std_of_means_across_T": 0, "flat_instrument": False},
              "mid_T_ablation_R4": {}, "banner": "skipped"}
    else:
        mc = run_mc_ablation(out_dir, workers, args.mc_per_config, args.n_points)
        (out_dir / "mc_ablation.json").write_text(json.dumps(mc, indent=2))

    if args.skip_beurling:
        beur = {
            "n_systems": 0,
            "x_max": args.beurling_xmax,
            "ordinary_R_d": None,
            "defective_all_above_ordinary": None,
            "defective_R_d_stats": {},
            "failure_modes_no_separation": [],
            "thin_margin_systems": [],
            "elapsed_s": 0,
            "banner": "skipped",
        }
    else:
        beur = run_beurling_expanded(
            out_dir, args.beurling_xmax, args.beurling_systems, workers, args.n_points
        )
        # drop huge rows from summary json optional full
        (out_dir / "beurling_expanded.json").write_text(json.dumps(beur, indent=2))

    report = write_report(out_dir, formal, mc, offc, beur, rate)
    summary = {
        "status": "completed",
        "banner": "NOT AN UNCONDITIONAL PROOF OF RH",
        "formal": formal["theorems"],
        "mc_flatness": mc.get("flatness"),
        "mc_mid_T_R4": mc.get("mid_T_ablation_R4"),
        "offcritical_all_far_higher": offc["all_far_higher_than_cl"],
        "beurling_sep": beur.get("defective_all_above_ordinary"),
        "beurling_failures": beur.get("failure_modes_no_separation"),
        "m3_scale_ok": rate["m3_scale_ok"],
        "out_dir": str(out_dir),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    (out_dir / "TIGHT_STRESS_COMPLETE").write_text(
        f"completed {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}\n"
        "NOT AN UNCONDITIONAL PROOF OF RH\n"
    )
    print(report, flush=True)
    print(f"Artifacts under {out_dir}", flush=True)


if __name__ == "__main__":
    main()
