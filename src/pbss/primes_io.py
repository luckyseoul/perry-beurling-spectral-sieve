"""
On-disk prime checkpoints for large-x segmented sieves.

Format: NumPy .npy int64 array (optionally memory-mapped for fork-shared workers).
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import List, Optional, Tuple, Union

import numpy as np

from .probes import primes_upto


def primes_checkpoint_path(root: Union[str, Path], x_max: int) -> Path:
    root = Path(root)
    return root / f"primes_le_{int(x_max)}.npy"


def meta_path(root: Union[str, Path], x_max: int) -> Path:
    return Path(root) / f"primes_le_{int(x_max)}.meta.json"


def save_primes_checkpoint(
    primes: np.ndarray,
    root: Union[str, Path],
    x_max: int,
    extra: Optional[dict] = None,
) -> Path:
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    path = primes_checkpoint_path(root, x_max)
    arr = np.asarray(primes, dtype=np.int64)
    np.save(path, arr)
    # np.save adds .npy if missing; normalize
    if not str(path).endswith(".npy"):
        path = Path(str(path) + ".npy")
    meta = {
        "x_max": int(x_max),
        "n_primes": int(arr.size),
        "dtype": "int64",
        "path": str(path),
        "created_unix": time.time(),
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    if extra:
        meta.update(extra)
    meta_path(root, x_max).write_text(json.dumps(meta, indent=2))
    return path


def load_primes_checkpoint(
    root: Union[str, Path],
    x_max: int,
    mmap: bool = True,
) -> Tuple[np.ndarray, dict]:
    root = Path(root)
    path = primes_checkpoint_path(root, x_max)
    if not path.exists():
        # np.save may have written path.npy
        alt = Path(str(path) + ".npy") if not str(path).endswith(".npy") else path
        if alt.exists():
            path = alt
        else:
            raise FileNotFoundError(path)
    primes = np.load(path, mmap_mode="r" if mmap else None)
    meta_file = meta_path(root, x_max)
    meta = json.loads(meta_file.read_text()) if meta_file.exists() else {}
    return primes, meta


# Fork workers inherit this read-only marking table (primes ≤ √new_max).
_WORKER_SIEVE_BASE: Optional[np.ndarray] = None


def _sieve_segment_job(payload: dict) -> np.ndarray:
    """Sieves one [low, high] interval using global _WORKER_SIEVE_BASE."""
    import os

    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"
    low = int(payload["low"])
    high = int(payload["high"])
    base_i = _WORKER_SIEVE_BASE
    if base_i is None:
        raise RuntimeError("sieve base not set for worker")
    mark = np.ones(high - low + 1, dtype=bool)
    for p in base_i:
        p = int(p)
        if p * p > high:
            break
        start = ((low + p - 1) // p) * p
        if start < p * p:
            start = p * p
        if start > high:
            continue
        mark[start - low : high - low + 1 : p] = False
    return (np.nonzero(mark)[0] + low).astype(np.int64)


def _extend_primes_from_base(
    base_primes: np.ndarray,
    old_max: int,
    new_max: int,
    segment_size: int = 50_000_000,
    workers: int = 0,
) -> np.ndarray:
    """
    Segmented sieve from old_max+1 .. new_max using base primes ≤ √new_max.

    Segments are sieved **in parallel** (ProcessPool, fork) — do not run
    multi-hour range extensions single-threaded on multi-core hosts.

    Returns only the *new* primes in (old_max, new_max]; caller concatenates.
    """
    import multiprocessing as mp
    import os
    from concurrent.futures import ProcessPoolExecutor

    global _WORKER_SIEVE_BASE

    old_max = int(old_max)
    new_max = int(new_max)
    if new_max <= old_max:
        return np.array([], dtype=np.int64)
    r = int(new_max**0.5) + 1
    base = np.asarray(base_primes)
    hi = int(np.searchsorted(base, r, side="right"))
    base_i = base[:hi].astype(np.int64, copy=False)
    if base_i.size == 0 or int(base_i[-1]) < r - 1:
        from .probes import primes_upto

        base_i = primes_upto(r).astype(np.int64)

    seg = int(segment_size)
    ranges: list[dict] = []
    low = old_max + 1
    while low <= new_max:
        high = min(low + seg - 1, new_max)
        ranges.append({"low": low, "high": high})
        low = high + 1
    if not ranges:
        return np.array([], dtype=np.int64)

    n_workers = int(workers) if workers and workers > 0 else max(
        1, (os.cpu_count() or 4) - 2
    )
    n_workers = min(n_workers, len(ranges))
    _WORKER_SIEVE_BASE = base_i
    ctx = mp.get_context("fork")
    # Preserve order: map returns in input order
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as ex:
        parts = list(ex.map(_sieve_segment_job, ranges, chunksize=1))
    _WORKER_SIEVE_BASE = None
    nonempty = [p for p in parts if p.size]
    if not nonempty:
        return np.array([], dtype=np.int64)
    return np.concatenate(nonempty)


def find_largest_checkpoint(
    checkpoint_dir: Union[str, Path],
) -> Optional[Tuple[int, Path]]:
    """Return (x_max, path) for largest primes_le_*.npy present, else None."""
    checkpoint_dir = Path(checkpoint_dir)
    if not checkpoint_dir.is_dir():
        return None
    best = None
    for p in checkpoint_dir.glob("primes_le_*.npy"):
        try:
            # primes_le_10000000000.npy
            stem = p.stem  # primes_le_10000000000
            x = int(stem.split("_")[-1])
        except ValueError:
            continue
        if best is None or x > best[0]:
            best = (x, p)
    return best


def ensure_primes(
    x_max: int,
    checkpoint_dir: Union[str, Path],
    *,
    force_resieve: bool = False,
    segment_size: int = 20_000_000,
    extend_from_existing: bool = True,
    workers: int = 0,
) -> Tuple[np.ndarray, dict]:
    """
    Load primes ≤ x_max from checkpoint, or sieve + save.

    If a smaller checkpoint exists and ``extend_from_existing``, extend via
    **parallel** segmented sieve from that table (ProcessPool over segments).
    """
    import os

    checkpoint_dir = Path(checkpoint_dir)
    x_max = int(x_max)
    path = primes_checkpoint_path(checkpoint_dir, x_max)
    alt = Path(str(path) + ".npy") if not str(path).endswith(".npy") else path
    if not force_resieve and (path.exists() or alt.exists()):
        primes, meta = load_primes_checkpoint(checkpoint_dir, x_max, mmap=True)
        meta["loaded_from_checkpoint"] = True
        return primes, meta

    n_workers = int(workers) if workers and workers > 0 else max(
        1, (os.cpu_count() or 4) - 2
    )
    t0 = time.time()
    method = "full_sieve"
    base_x = None
    if extend_from_existing and not force_resieve:
        found = find_largest_checkpoint(checkpoint_dir)
        if found is not None and found[0] < x_max and found[0] >= 1000:
            base_x, _ = found
            base, _ = load_primes_checkpoint(checkpoint_dir, base_x, mmap=False)
            # copy base into RAM for concat
            base_arr = np.asarray(base, dtype=np.int64).copy()
            new_part = _extend_primes_from_base(
                base_arr,
                base_x,
                x_max,
                segment_size=max(segment_size, 20_000_000),
                workers=n_workers,
            )
            primes = (
                np.concatenate([base_arr, new_part])
                if new_part.size
                else base_arr
            )
            method = f"extend_from_{base_x}_parallel_w{n_workers}"
        else:
            primes = primes_upto(x_max, segment_size=segment_size)
    else:
        primes = primes_upto(x_max, segment_size=segment_size)
    elapsed = time.time() - t0
    path = save_primes_checkpoint(
        primes,
        checkpoint_dir,
        x_max,
        extra={
            "sieve_seconds": elapsed,
            "segment_size": segment_size,
            "method": method,
            "extended_from": base_x,
            "workers": n_workers,
        },
    )
    # free dense array before memmap reload
    del primes
    primes_mm, meta = load_primes_checkpoint(checkpoint_dir, x_max, mmap=True)
    meta["loaded_from_checkpoint"] = False
    meta["sieve_seconds"] = elapsed
    meta["saved_path"] = str(path)
    meta["method"] = method
    return primes_mm, meta
