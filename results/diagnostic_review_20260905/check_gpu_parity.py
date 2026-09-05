import json
import os
import sys
import time
from pathlib import Path
# Match this process's CUDA headers to its installed CUDA 12.8 NVRTC.
os.environ['CUDA_PATH'] = '/usr/local/cuda-12.8'
import numpy as np
sys.path.insert(0, '/home/nick/perry-beurling-spectral-sieve/src')
from pbss.projection import energy_ratio, project_coefficients
from pbss.projection_backend import energy_ratio_auto, project_coefficients_auto, cupy_available
assert cupy_available(), 'CuPy device is required for this check'
import cupy as cp
rows = []
for n, d, grid_power, weighted in [(8192, 4, 1.0, False), (16384, 8, 1.7, True), (65536, 12, 1.0, True)]:
    u = np.linspace(0., 1., n) ** grid_power
    q = np.sin(2*np.pi*48*u) + 0.7*(u-0.5) + 0.04*np.cos(350*u)
    w = None if not weighted else 0.1 + u*u
    cpu = energy_ratio(q, u, d, weights=w)
    a_cpu = project_coefficients(q, u, d, weights=w)
    cp.cuda.Stream.null.synchronize()
    start = time.perf_counter()
    event_start, event_end = cp.cuda.Event(), cp.cuda.Event()
    event_start.record()
    gpu, backend = energy_ratio_auto(q,u,d,prefer_gpu=True,weights=w)
    a_gpu, coefficient_backend = project_coefficients_auto(q,u,d,weights=w,prefer_gpu=True)
    event_end.record()
    event_end.synchronize()
    wall = time.perf_counter()-start
    if backend != 'cupy':
        from pbss.projection import _fit_projection, _prepare_projection_inputs
        _fit_projection(*_prepare_projection_inputs(q,u,d,w), xp=cp)
    assert backend == coefficient_backend == 'cupy', (backend, coefficient_backend)
    assert np.isclose(cpu, gpu, rtol=1e-11, atol=1e-13), (cpu,gpu)
    assert np.allclose(a_cpu,a_gpu,rtol=1e-10,atol=1e-12)
    rows.append(dict(n_points=n, degree=d, grid_power=grid_power, custom_weights=weighted,
                     cpu_Rd=cpu, gpu_Rd=gpu, absolute_error=abs(cpu-gpu), backend=backend,
                     wall_seconds=wall, cuda_event_elapsed_ms=cp.cuda.get_elapsed_time(event_start,event_end)))
report = dict(device=cp.cuda.runtime.getDeviceProperties(0)['name'].decode(), cupy_version=cp.__version__,
              environment={'CUDA_PATH': os.environ['CUDA_PATH'], 'CUPY_ACCELERATORS': os.environ.get('CUPY_ACCELERATORS', '(default)')},
              cuda_runtime_version=cp.cuda.runtime.runtimeGetVersion(), nvrtc_version=list(cp.cuda.nvrtc.getVersion()),
              environment_reference='https://docs.cupy.dev/en/stable/reference/environment.html',
              note='CUDA events span host orchestration and GPU operations; elapsed time is not a kernel-only benchmark.', rows=rows)
path = Path('/home/nick/perry-beurling-spectral-sieve/results/diagnostic_review_20260905/gpu_parity.json')
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps(report,indent=2))
