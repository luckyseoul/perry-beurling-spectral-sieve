"""
On-disk prime checkpoints for large-x segmented sieves.

Format: NumPy .npy int64 array (optionally memory-mapped for fork-shared workers).
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Optional, Tuple, Union

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


def ensure_primes(
    x_max: int,
    checkpoint_dir: Union[str, Path],
    *,
    force_resieve: bool = False,
    segment_size: int = 20_000_000,
) -> Tuple[np.ndarray, dict]:
    """
    Load primes ≤ x_max from checkpoint, or sieve + save.

    Uses segmented sieve for large n (see probes.primes_upto).
    """
    checkpoint_dir = Path(checkpoint_dir)
    x_max = int(x_max)
    path = primes_checkpoint_path(checkpoint_dir, x_max)
    alt = Path(str(path) + ".npy") if not str(path).endswith(".npy") else path
    if not force_resieve and (path.exists() or alt.exists()):
        primes, meta = load_primes_checkpoint(checkpoint_dir, x_max, mmap=True)
        meta["loaded_from_checkpoint"] = True
        return primes, meta

    t0 = time.time()
    primes = primes_upto(x_max, segment_size=segment_size)
    elapsed = time.time() - t0
    path = save_primes_checkpoint(
        primes,
        checkpoint_dir,
        x_max,
        extra={"sieve_seconds": elapsed, "segment_size": segment_size},
    )
    # reload as memmap for shared use
    primes_mm, meta = load_primes_checkpoint(checkpoint_dir, x_max, mmap=True)
    meta["loaded_from_checkpoint"] = False
    meta["sieve_seconds"] = elapsed
    meta["saved_path"] = str(path)
    return primes_mm, meta
