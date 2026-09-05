"""Regression for a threshold crossed inside the allowed refinement tolerance."""
import importlib.util
from pathlib import Path


def test_refinement_threshold_crossing_cannot_pass(monkeypatch):
    path = Path(__file__).resolve().parents[1] / "experiments" / "run_diagnostic.py"
    spec = importlib.util.spec_from_file_location("diagnostic_experiment", path)
    experiment = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(experiment)
    coarse = [0.2995, 0.2995, 0.19, 0.8]
    fine = [0.3, 0.3, 0.19, 0.794]
    rows = []
    for name, a, b in zip(experiment.PROBES, coarse, fine):
        tolerance = 0.001 + 0.01 * b
        assert abs(a - b) <= tolerance
        rows.append({"name": name, "energy_ratio": a,
                     "refinement": {"energy_ratio": b, "stable": True}})

    class CompletedPool:
        def __init__(self, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def map(self, function, jobs):
            return rows

    monkeypatch.setattr(experiment, "ProcessPoolExecutor", CompletedPool)
    report = experiment.run()
    assert report["separation"]["status"] == "passed"
    assert report["refined_separation"]["status"] == "failed"
    assert report["grid_refinement"]["stable"]
    assert report["status"] == "inconclusive"
