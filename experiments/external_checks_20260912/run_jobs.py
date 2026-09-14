#!/usr/bin/env python3
"""Run a slice of model_family cases in parallel (ProcessPoolExecutor) on this node."""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model_family import run_case  # noqa: E402


def main():
    jobs_file = sys.argv[1]
    out_dir = sys.argv[2]
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    os.makedirs(out_dir, exist_ok=True)
    with open(jobs_file) as fh:
        jobs = json.load(fh)
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(run_case, job): job for job in jobs}
        for fut in as_completed(futs):
            job = futs[fut]
            try:
                res = fut.result()
            except Exception as exc:  # noqa: BLE001
                print(f"FAIL {job['name']} T={job['Ts']}: {type(exc).__name__}: {exc}", flush=True)
                continue
            name = f"{job['name']}_T{job['Ts'][0]:g}"
            path = os.path.join(out_dir, name + ".json")
            with open(path, "w") as fh:
                json.dump({"results": [res], "elapsed_s": time.time() - t0}, fh)
            v = res["per_T"][str(res["case"]["Ts"][0])] if str(res["case"]["Ts"][0]) in res["per_T"] else list(res["per_T"].values())[0]
            print(f"DONE {name} N={v['N']} n*E={v['nstar_energy']} log10Rd_d2={v['log10_Rd_d2']:.3f} "
                  f"t={time.time()-t0:.1f}s", flush=True)
    print(f"ALL DONE {len(jobs)} jobs in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
