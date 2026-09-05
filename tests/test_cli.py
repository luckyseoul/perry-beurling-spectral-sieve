"""Regression checks for user input, reproducible reports, and truthful status."""
from __future__ import annotations

import json

import numpy as np
import pytest

from pbss.cli import _write_json, main


def test_diagnose_infers_csv_grid_and_reproduces_report(tmp_path, capsys):
    samples = tmp_path / "samples.csv"
    np.savetxt(samples, np.linspace(1.0, 2.0, 37), delimiter=",")
    output = tmp_path / "reports" / "diagnostic.json"
    command = ["diagnose", "-i", str(samples), "--json-out", str(output)]

    assert main(command) == 0
    first = output.read_text()
    report = json.loads(first)
    assert report["mode"] == "input"
    assert report["n_points"] == report["rows"][0]["n_points"] == 37
    assert report["rows"][0]["Rd"] == pytest.approx(1.0)
    assert report["rows"][0]["l2_norm_sq"] > 0
    assert report["projection_method"] == "weighted_qr"
    assert report["grid"]["source"] == "uniform_linspace"
    assert len(report["input"]["sha256_float64_le"]) == 64
    assert "assessment" not in report
    assert main(command) == 0
    assert output.read_text() == first
    assert "n_points=37" in capsys.readouterr().out


@pytest.mark.parametrize("command", ["project", "diagnose"])
def test_commands_project_polynomial_on_custom_grid(tmp_path, command):
    grid = np.array([0.0, 0.03, 0.2, 0.6, 0.73, 0.88, 1.0])
    samples = tmp_path / "samples.npy"
    abscissae = tmp_path / "u.npy"
    output = tmp_path / "result.json"
    np.save(samples, 1.0 + 3.0 * grid - grid**2)
    np.save(abscissae, grid)
    assert main([
        command, "-i", str(samples), "--u", str(abscissae), "-d", "2",
        "--T", "3", "--json-out", str(output),
    ]) == 0
    report = json.loads(output.read_text())
    ratio = report["energy_ratio_Rd"] if command == "project" else report["rows"][0]["Rd"]
    assert ratio == pytest.approx(1.0)
    assert report["grid"]["source"] == str(abscissae)
    assert report["grid"]["n_points"] == grid.size
    if command == "project":
        assert report["coefficient_meaning"] == "weighted_least_squares_polynomial_coefficients"
        assert report["energy"] == pytest.approx(report["l2_norm_sq"])


@pytest.mark.parametrize("values, message", [
    (np.ones((3, 3)), "one-dimensional"),
    (np.array([1.0 + 2.0j, 3.0j]), "real numbers"),
    (np.array([1.0, np.nan]), "finite"),
    (np.array([np.inf, 1.0]), "finite"),
    (np.array([1.0]), "at least 2"),
    (np.zeros(8), "zero"),
])
def test_invalid_sample_arrays_fail_without_traceback(tmp_path, capsys, values, message):
    samples = tmp_path / "samples.npy"
    np.save(samples, values)
    with pytest.raises(SystemExit) as caught:
        main(["project", "-i", str(samples)])
    assert caught.value.code == 2
    error = capsys.readouterr().err
    assert message in error
    assert "Traceback" not in error


def test_missing_input_reports_path_without_traceback(tmp_path, capsys):
    missing = tmp_path / "missing.csv"
    with pytest.raises(SystemExit) as caught:
        main(["diagnose", "-i", str(missing)])
    assert caught.value.code == 2
    error = capsys.readouterr().err
    assert str(missing) in error
    assert "Traceback" not in error


def test_empty_numpy_file_fails_cleanly(tmp_path, capsys):
    samples = tmp_path / "empty.npy"
    samples.touch()
    with pytest.raises(SystemExit) as caught:
        main(["project", "-i", str(samples)])
    assert caught.value.code == 2
    error = capsys.readouterr().err
    assert str(samples) in error
    assert "Traceback" not in error


@pytest.mark.parametrize("grid", [np.linspace(1, 0, 8), np.linspace(0, 1, 7)])
def test_invalid_custom_grid_fails_without_report(tmp_path, grid):
    samples = tmp_path / "samples.npy"
    abscissae = tmp_path / "u.npy"
    output = tmp_path / "report.json"
    np.save(samples, np.ones(8))
    np.save(abscissae, grid)
    with pytest.raises(SystemExit) as caught:
        main([
            "diagnose", "-i", str(samples), "--u", str(abscissae),
            "--json-out", str(output),
        ])
    assert caught.value.code == 2
    assert not output.exists()


@pytest.mark.parametrize("arguments, message", [
    (["diagnose", "--demo", "-i", "q.npy"], "not allowed with"),
    (["diagnose", "-i", "q.npy", "--n-points", "8"], "only to the demo"),
    (["diagnose", "-i", "q.npy", "--with-primes"], "only to the demo"),
    (["diagnose", "--u", "u.npy"], "requires --input"),
    (["diagnose", "--x-max", "100"], "requires --with-primes"),
    (["sensitivity", "--confirm-53", "--k", "4"], "fixed ensemble"),
    (["sensitivity", "--confirm-53", "--n-per-class", "300"], "fixed ensemble"),
    (["sensitivity", "--confirm-53", "--sigma", "6"], "fixed ensemble"),
    (["diagnose", "--T", "nan"], "finite"),
    (["diagnose", "--T", "0"], "positive"),
    (["diagnose", "--degree", "-1"], "at least 0"),
    (["diagnose", "--n-points", "1"], "at least 2"),
    (["scorecard", "--x-max", "inf"], "finite"),
    (["scorecard", "--x-max", "2"], "at least 3"),
    (["sensitivity", "--noise", "-1"], "nonnegative"),
    (["sensitivity", "--n-per-class", "1"], "at least 2"),
])
def test_invalid_and_ignored_flags_are_errors(capsys, arguments, message):
    with pytest.raises(SystemExit) as caught:
        main(arguments)
    assert caught.value.code == 2
    assert message in capsys.readouterr().err


def test_input_output_collision_does_not_destroy_samples(tmp_path):
    samples = tmp_path / "samples.npy"
    np.save(samples, np.ones(8))
    before = samples.read_bytes()
    with pytest.raises(SystemExit) as caught:
        main(["project", "-i", str(samples), "--json-out", str(samples)])
    assert caught.value.code == 2
    assert samples.read_bytes() == before


def test_nonfinite_report_does_not_replace_existing_output(tmp_path):
    output = tmp_path / "report.json"
    output.write_text('{"previous": true}\n')
    with pytest.raises(ValueError, match="JSON compliant"):
        _write_json(str(output), {"result": float("nan")})
    assert json.loads(output.read_text()) == {"previous": True}


@pytest.mark.parametrize("arguments, status, exit_code", [
    (["--n-points", "1024"], "passed", 0),
    (["--n-points", "32"], "inconclusive", 2),
    (["--degree", "0"], "inconclusive", 2),
])
def test_demo_exit_status_matches_assessment(tmp_path, arguments, status, exit_code):
    output = tmp_path / "demo.json"
    assert main(["diagnose", "--demo", *arguments, "--json-out", str(output)]) == exit_code
    report = json.loads(output.read_text())
    assert report["assessment"]["status"] == status
    assert report["assessment"]["classifier_ok"] is (exit_code == 0)
    assert report["demo_parameters"]["waves"] == 48


def test_scorecard_partial_failure_is_retained_and_nonzero(tmp_path, monkeypatch):
    import pbss.beurling as beurling

    original = beurling.build_system_primes

    def build(spec, ordinary_primes, *, x_max):
        if spec["name"] == "gapped_gap3":
            raise ValueError("deliberate system failure")
        return original(spec, ordinary_primes, x_max=x_max)

    monkeypatch.setattr(beurling, "build_system_primes", build)
    output = tmp_path / "scorecard.json"
    assert main([
        "scorecard", "--x-max", "1000", "--n-points", "65",
        "--json-out", str(output),
    ]) == 2
    report = json.loads(output.read_text())
    assert report["status"] == "failed"
    assert [row["status"] for row in report["rows"]] == ["ok", "error", "ok"]
    assert report["rows"][1]["error"] == "deliberate system failure"
    assert report["degree"] == 4
    assert report["n_points"] == 65


def test_scorecard_runs_all_real_systems(tmp_path):
    output = tmp_path / "scorecard.json"
    assert main([
        "scorecard", "--x-max", "1000", "--n-points", "65",
        "--json-out", str(output),
    ]) == 0
    report = json.loads(output.read_text())
    assert report["status"] == "ok"
    assert [row["system"] for row in report["rows"]] == [
        "ordinary_primes", "gapped_gap3", "thinned_every3",
    ]
    assert [row["n_primes"] for row in report["rows"]] == [168, 333, 56]
    for row in report["rows"]:
        assert row["status"] == "ok"
        assert 0 <= row["Rd"] <= 1
        assert np.isfinite(row["Sd"])
        assert row["residual"]["detrend"] == "deg1"


def test_unwritable_report_path_is_a_clean_error(tmp_path, capsys):
    samples = tmp_path / "samples.npy"
    np.save(samples, np.ones(8))
    with pytest.raises(SystemExit) as caught:
        main(["project", "-i", str(samples), "--json-out", str(tmp_path)])
    assert caught.value.code == 2
    assert "Traceback" not in capsys.readouterr().err


def test_moments_cli_retains_signed_absolute_integrals_and_strict_json(tmp_path, capsys):
    output = tmp_path / 'moments' / 'result.json'
    # Before the first prime, Q(y)=-exp(y/2). This gives an independent
    # closed-form degree-zero moment and full continuous norm.
    T = 0.4
    assert main(['moments', '--T', str(T), '-d', '2', '--json-out', str(output)]) == 0
    report = json.loads(output.read_text())
    assert report['command'] == 'moments'
    assert report['metadata']['n_primes'] == 0
    assert float(report['coefficients_decimal'][0]) == pytest.approx(-2 * np.expm1(T / 2) / T)
    assert float(report['l2_norm_sq_decimal']) == pytest.approx(np.expm1(T) / T)
    assert report['coefficient_signs'][0] == -1
    assert [float(v) for v in report['affine_detrended_coefficients_decimal'][:2]] == [0.0, 0.0]
    assert report['metadata']['certified_interval_result'] is False
    assert report['metadata']['B_ABS_growth_bound_proved'] is False
    assert report['metadata']['RH_proved'] is False
    text = capsys.readouterr().out
    assert 'c_2=' in text and 'l2_norm_sq=' in text
    assert 'no all-window growth bound' in text


@pytest.mark.parametrize('arguments', [
    ['--T', 'nan'], ['--T', '-1'], ['--T', '100'],
    ['--degree', '33'], ['--max-prime-limit', '1'],
])
def test_moments_cli_rejects_unresolved_or_out_of_scope_requests(arguments, tmp_path, capsys):
    output = tmp_path / 'not_created.json'
    with pytest.raises(SystemExit) as caught:
        main(['moments', *arguments, '--json-out', str(output)])
    assert caught.value.code == 2
    assert not output.exists()
    assert 'Traceback' not in capsys.readouterr().err
