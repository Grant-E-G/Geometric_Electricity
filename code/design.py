"""Reproducible LC experiment: topology, loaded simulation, KiCad and analysis.

Python is used for SciPy and KiCad's native pcbnew API; see README for commands.
All intermediate reports and manufacturing previews live in ignored build/.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import socket
import subprocess
import uuid
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/geometric-electricity-mpl")
import numpy as np
from numpy.typing import NDArray
from scipy import linalg, signal
from scipy.optimize import (
    least_squares,
    linear_sum_assignment,
    minimize,
    minimize_scalar,
)
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "build"
HW = ROOT / "hardware"
FloatArray = NDArray[np.float64]
L = 1e-6
C = 330e-12
FP_C = "Capacitor_SMD:C_0603_1608Metric"
FP_L = "Inductor_SMD:L_0805_2012Metric"
FP_T = "TestPoint:TestPoint_Pad_D1.0mm"


def uid(name: str) -> str:
    return str(uuid.uuid5(uuid.NAMESPACE_URL, "geometric-electricity/" + name))


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n")


def graphs() -> list[dict[str, Any]]:
    """Construct induced regular disks, retaining metric and PCB coordinates."""
    e = sorted(
        (
            complex(a + b / 2, b * math.sqrt(3) / 2)
            for a in range(-6, 7)
            for b in range(-6, 7)
            if a * a + a * b + b * b <= 21
        ),
        key=lambda z: (round(abs(z), 7), math.atan2(z.imag, z.real)),
    )
    cosh_edge = math.cos(2 * math.pi / 7) / (1 - math.cos(2 * math.pi / 7))
    edge = math.acosh(cosh_edge)
    radius = math.tanh(edge / 2)
    h = [0j]
    depth = [0]
    for i in range(29):  # center + first two graph shells
        z = h[i]
        # A known neighbor fixes orientation; at center choose positive real.
        neighbor = (
            radius
            if i == 0
            else next(w for w in h[:i] if abs(hyperdistance(z, w) - edge) < 1e-7)
        )
        local = (neighbor - z) / (1 - z.conjugate() * neighbor)
        for k in range(7):
            wlocal = local * np.exp(2j * math.pi * k / 7)
            w = complex((wlocal + z) / (1 + z.conjugate() * wlocal))
            if not any(abs(w - old) < 1e-7 for old in h):
                h.append(w)
                depth.append(depth[i] + 1)
    h = sorted(h, key=lambda z: (round(abs(z), 7), math.atan2(z.imag, z.real)))
    result = []
    for name, q, points, target, center in [
        ("E", 6, e, 1.0, 38.0),
        ("H", 7, h, edge, 112.0),
    ]:
        n = len(points)
        assert n == 85, (name, n)
        metric = np.array(
            [
                [abs(a - b) if q == 6 else hyperdistance(a, b) for b in points]
                for a in points
            ]
        )
        adjacency = np.abs(metric - target) < 1e-7
        edges = [(i, j) for i in range(n) for j in range(i + 1, n) if adjacency[i, j]]
        shells = [0] + [99] * (n - 1)
        for _ in range(n):
            for a, b in edges:
                shells[b] = min(shells[b], shells[a] + 1)
                shells[a] = min(shells[a], shells[b] + 1)
        missing = q - adjacency.sum(axis=1)
        assert min(missing) >= 0 and max(shells) < 99
        # Spread hyperbolic shells instead of cramming the Poincare boundary.
        physical = []
        for i, z in enumerate(points):
            r = (
                abs(z) / max(map(abs, points))
                if q == 6
                else [0, 0.25, 0.59, 1][shells[i]]
            )
            theta = math.atan2(z.imag, z.real)
            physical.append(
                [center + 29 * r * math.cos(theta), 65 + 48 * r * math.sin(theta)]
            )
        # Sources: center, representative intermediate, boundary, rotated peer.
        intermediate = next(i for i in range(n) if shells[i] == 2)
        rotated = points[intermediate] * np.exp(2j * np.pi / q)
        peer = int(np.argmin([abs(z - rotated) for z in points]))
        sources = [0, intermediate, int(np.argmax(np.real(points))), peer]
        result.append(
            {
                "name": name,
                "coordination": q,
                "nodes": n,
                "edges": edges,
                "metric_xy": [[z.real, z.imag] for z in points],
                "pcb_xy_mm": physical,
                "shells": shells,
                "missing": missing.tolist(),
                "sources": sources,
                "adjacency": [np.flatnonzero(row).tolist() for row in adjacency],
                "edge_distance": target,
            }
        )
    return result


def hyperdistance(a: complex, b: complex) -> float:
    return math.acosh(
        max(1.0, 1 + 2 * abs(a - b) ** 2 / ((1 - abs(a) ** 2) * (1 - abs(b) ** 2)))
    )


def capacitance(
    g: dict[str, Any],
    edge_values: FloatArray | None = None,
    boundary_values: FloatArray | None = None,
    stray_pf: float = 0,
) -> FloatArray:
    values = np.full(len(g["edges"]), C) if edge_values is None else edge_values
    boundary = (
        C * np.array(g["missing"]) if boundary_values is None else boundary_values
    )
    matrix = np.diag(boundary + stray_pf * 1e-12)
    for (a, b), c in zip(g["edges"], values):
        matrix[a, a] += c
        matrix[b, b] += c
        matrix[a, b] -= c
        matrix[b, a] -= c
    return matrix


def eigenpairs(cm: FloatArray, inductance: FloatArray) -> tuple[FloatArray, FloatArray]:
    # cm v = lambda diag(1/Li) v; lambda=1/omega^2. Ascending lambda gives ascending q and descending electrical frequency.
    values, vectors = linalg.eigh(cm, np.diag(1 / inductance))
    frequencies = 1 / (2 * np.pi * np.sqrt(values))
    vectors /= np.linalg.norm(vectors, axis=0)
    return frequencies, vectors


def groups(freq: FloatArray) -> list[list[int]]:
    result: list[list[int]] = [[0]]
    for i in range(1, len(freq)):
        if abs(freq[i] - freq[i - 1]) < 1e-6 * freq[i]:
            result[-1].append(i)
        else:
            result.append([i])
    return result


def response(
    g: dict[str, Any],
    frequencies: FloatArray,
    source: int,
    probes: list[int],
    probe_pf: float,
    quality: float,
    sense: float = 1000,
) -> NDArray[np.complex128]:
    cm = capacitance(g, stray_pf=2)
    for node in probes:
        cm[node, node] += probe_pf * 1e-12
    resistance = 2 * np.pi * 8e6 * L / quality
    result = []
    for f in frequencies:
        w = 2 * np.pi * f
        y = 1j * w * cm.astype(complex)
        y[np.diag_indices(85)] += 1 / (resistance + 1j * w * L)
        for node in probes:
            y[node, node] += 1 / 10e6
        y[source, source] += 1 / (sense + 50)
        current = np.zeros(85, dtype=complex)
        current[source] = 1 / (sense + 50)
        result.append(linalg.solve(y, current))
    return np.array(result)


def select_defects(g: dict[str, Any], vectors: FloatArray) -> list[int]:
    # Perturbation leverage for long-wavelength modes, distributed over shells.
    sensitivity = np.array(
        [np.sum((vectors[a, :12] - vectors[b, :12]) ** 2) for a, b in g["edges"]]
    )
    picked = []
    used: set[int] = set()
    for shell in range(1, max(g["shells"]) + 1):
        candidates = sorted(
            (
                k
                for k, (a, b) in enumerate(g["edges"])
                if max(g["shells"][a], g["shells"][b]) == shell
            ),
            key=lambda k: -sensitivity[k],
        )
        for k in candidates:
            a, b = g["edges"][k]
            if a not in used and b not in used:
                picked.append(k)
                used.update((a, b))
                break
    for k in np.argsort(-sensitivity):
        a, b = g["edges"][int(k)]
        if a not in used and b not in used:
            picked.append(int(k))
            used.update((a, b))
        if len(picked) == 6:
            break
    return picked


def defect_sweep(g: dict[str, Any], nominal_f: FloatArray) -> list[dict[str, Any]]:
    """Predict every reversible sparse control before assembly."""
    cases: list[tuple[str, FloatArray, FloatArray]] = []
    for index in g["defects"]:
        for factor in [0.0, 0.5, 2.0]:
            edge = np.full(len(g["edges"]), C)
            edge[index] *= factor
            cases.append(
                (
                    f"edge {g['edges'][index]} C×{factor:g}",
                    edge,
                    C * np.array(g["missing"]),
                )
            )
    boundary_node = next(i for i in g["sources"] if g["missing"][i] > 0)
    boundary = C * np.array(g["missing"])
    boundary[boundary_node] -= C
    cases.append(
        (
            f"boundary {boundary_node}: remove one 330pF",
            np.full(len(g["edges"]), C),
            boundary,
        )
    )
    results = []
    for label, changed_edge, changed_boundary in cases:
        f, v = eigenpairs(
            capacitance(g, changed_edge, changed_boundary), np.full(85, L)
        )
        results.append(
            {
                "control": label,
                "first12_MHz": (f[:12] / 1e6).tolist(),
                "max_first12_sorted_frequency_shift": float(
                    np.max(abs(f[:12] / nominal_f[:12] - 1))
                ),
                "first12_ipr": np.sum(v[:, :12] ** 4, axis=0).tolist(),
                "interpretation": "Sorted spectrum and IPR prediction; mode identity requires spatial comparison.",
            }
        )
    return results


def simulate(gs: list[dict[str, Any]], trials: int) -> dict[str, Any]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    summary: dict[str, Any] = {"seed": 20261005, "trials": trials, "geometries": {}}
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for row, g in enumerate(gs):
        cm = capacitance(g)
        freq, modes = eigenpairs(cm, np.full(85, L))
        qvals = 1 / ((2 * np.pi * freq) ** 2 * L * C)
        g["defects"] = select_defects(g, modes)
        clusters = groups(freq)
        statistics = []
        for tolerance, quality, cap_tolerance in [
            (0.1, 45, 0.05),
            (0.1, 36, 0.05),
            (0.1, 45, 0.02),
            (0.05, 50, 0.02),
            (0.05, 30, 0.02),
            (0.02, 50, 0.02),
            (0.05, 100, 0.02),
        ]:
            rng = np.random.default_rng(20261005 + row)
            shifts, overlaps, resolved = [], [], []
            for _ in range(trials):
                li = L * (1 + rng.uniform(-tolerance, tolerance, 85))
                ce = C * (
                    1 + rng.uniform(-cap_tolerance, cap_tolerance, len(g["edges"]))
                )
                cb = (
                    C
                    * np.array(g["missing"])
                    * (1 + rng.uniform(-cap_tolerance, cap_tolerance, 85))
                )
                f, v = eigenpairs(capacitance(g, ce, cb, 2), li)
                # Assign all modes, then compare degenerate subspaces using principal angles.
                a, b = linear_sum_assignment(-((modes.T @ v) ** 2))
                mapped = dict(zip(a.tolist(), b.tolist()))
                shifts.append(np.max(np.abs(f[:12] / freq[:12] - 1)))
                angles = []
                for cluster in clusters[:6]:
                    basis, _ = np.linalg.qr(v[:, [mapped[k] for k in cluster]])
                    angles.append(
                        float(linalg.svdvals(modes[:, cluster].T @ basis).min() ** 2)
                    )
                overlaps.append(min(angles))
                # Separate distinct low-q clusters, not symmetry partners.
                cf = np.array([np.mean(f[c]) for c in clusters[:6]])
                # Fixed series resistance: Q scales with frequency; no guarantee from single-point Q.
                linewidth = 8e6 / quality
                resolved.append(bool(np.all(np.abs(np.diff(cf)) > linewidth)))
            statistics.append(
                {
                    "L_tolerance": tolerance,
                    "C_tolerance": cap_tolerance,
                    "Q_at_8MHz": quality,
                    "max_first12_frequency_shift_p95": float(np.quantile(shifts, 0.95)),
                    "first6_subspace_overlap_p05": float(np.quantile(overlaps, 0.05)),
                    "first6_clusters_gap_over_linewidth_pass_fraction": float(
                        np.mean(resolved)
                    ),
                }
            )
        loads = []
        source = g["sources"][1]
        scan = g["sources"][2]
        for stray in [0, 2, 5, 10]:
            for probe in [0, 10, 15, 50]:
                loaded = capacitance(g, stray_pf=stray)
                loaded[source, source] += probe * 1e-12
                loaded[scan, scan] += probe * 1e-12
                f, _ = eigenpairs(loaded, np.full(85, L))
                loads.append(
                    {
                        "stray_pf_per_node": stray,
                        "probe_pf": probe,
                        "first12_max_shift": float(
                            np.max(np.abs(f[:12] / freq[:12] - 1))
                        ),
                    }
                )
        c_choices = []
        for c_pf in [220, 330, 470, 680, 1000]:
            trial_cm = cm * (c_pf / 330) + np.eye(85) * 2e-12
            trial_cm[source, source] += 10e-12
            trial_cm[scan, scan] += 10e-12
            f, _ = eigenpairs(trial_cm, np.full(85, L))
            c_choices.append(
                {
                    "C_pf": c_pf,
                    "min_MHz": float(f.min() / 1e6),
                    "max_MHz": float(f.max() / 1e6),
                }
            )
        sweep = np.linspace(1e6, 18e6, 1001, dtype=np.float64)
        voltage = response(g, sweep, source, [source, scan], 10, 36)
        np.savez_compressed(
            OUT / f'{g["name"]}_simulation.npz',
            frequency_hz=freq,
            q=qvals,
            vectors=modes,
            sweep_hz=sweep,
            loaded_voltage=voltage,
        )
        # Exact damped modal current-pulse response, source loading omitted here.
        t = np.arange(0, 6e-6, 2e-9)
        injection = 0.1e-3 * np.exp(-0.5 * ((t - 0.3e-6) / 30e-9) ** 2)
        resistance = 2 * np.pi * 8e6 * L / 36
        eigen_c, v = linalg.eigh(capacitance(g, stray_pf=2))
        pulse_source = g["sources"][2]
        responses = []
        for k, c in enumerate(eigen_c):
            system = signal.TransferFunction(
                [L, resistance], [L * c, resistance * c, 1]
            )
            _, y, _ = signal.lsim(system, injection * v[pulse_source, k], t)
            responses.append(y)
        wave = np.array(responses).T @ v.T
        np.savez_compressed(
            OUT / f'{g["name"]}_pulse.npz',
            time_s=t,
            current_a=injection,
            voltage_v=wave,
        )
        # Derived finite-graph quantities; complete synthetic eigenbasis.
        times = np.logspace(-3, 2, 150)
        heat = np.exp(-np.outer(times, qvals)).sum(axis=1)
        zeta = np.array([np.sum(qvals ** (-s)) for s in [1.0, 2.0, 3.0]])
        green = (modes / qvals) @ modes.T
        np.savez_compressed(
            OUT / f'{g["name"]}_derived.npz',
            t=times,
            heat_trace=heat,
            zeta_s=[1.0, 2.0, 3.0],
            zeta=zeta,
            green_Q_inverse=green,
        )
        xy = np.array(g["pcb_xy_mm"])
        for a, b in g["edges"]:
            axes[row, 0].plot(xy[[a, b], 0], xy[[a, b], 1], color="0.8", lw=0.6)
        axes[row, 0].scatter(xy[:, 0], xy[:, 1], c=g["shells"], s=16)
        axes[row, 0].set_title(f'{g["name"]}: 85 nodes, {len(g["edges"])} edges')
        axes[row, 0].set_aspect("equal")
        axes[row, 1].plot(sweep / 1e6, abs(voltage[:, source]), label="source")
        axes[row, 1].plot(sweep / 1e6, abs(voltage[:, scan]), label="boundary")
        axes[row, 1].set(
            xlabel="MHz",
            ylabel="V / 1 V generator",
            title="1 kΩ injection, Q(8 MHz)=36, 10 pF probes",
        )
        axes[row, 1].legend()
        axes[row, 2].imshow(
            wave.T, aspect="auto", extent=[0, t[-1] * 1e6, 84, 0], cmap="RdBu_r"
        )
        axes[row, 2].set(
            xlabel="µs",
            ylabel="node ID",
            title="Ideal current pulse; series loss, no source resistor",
        )
        summary["geometries"][g["name"]] = {
            "edges": len(g["edges"]),
            "shells": [g["shells"].count(k) for k in range(max(g["shells"]) + 1)],
            "boundary_nodes": sum(x > 0 for x in g["missing"]),
            "missing_couplings": sum(g["missing"]),
            "min_MHz": float(freq.min() / 1e6),
            "max_MHz": float(freq.max() / 1e6),
            "first12_MHz": (freq[:12] / 1e6).tolist(),
            "first6_clusters_MHz": [float(freq[c].mean() / 1e6) for c in clusters[:6]],
            "monte_carlo": statistics,
            "loading": loads,
            "capacitor_choices": c_choices,
            "sources": g["sources"],
            "defect_edges": [g["edges"][k] for k in g["defects"]],
            "defect_sweep": defect_sweep(g, freq),
        }
    fig.tight_layout()
    fig.savefig(OUT / "simulation.png", dpi=160)
    plt.close(fig)
    write_json(OUT / "simulation.json", summary)
    return summary


def network_admittance(
    g: dict[str, Any],
    hz: float,
    *,
    quality: float = 36,
    inductance: FloatArray | None = None,
    srf_mhz: float | None = None,
    cap_esr: float = 0,
    cap_esl_nh: float = 0,
    trace_nh: float = 0,
    ground_nh: float = 0,
    stray_pf: float | FloatArray = 2,
    mutual: tuple[int, int, float] | None = None,
) -> NDArray[np.complex128]:
    """Passive lumped sensitivity model, not field-extracted PCB parasitics."""
    li = np.full(85, L) if inductance is None else inductance
    w = 2 * np.pi * hz
    branches = np.diag(
        2 * np.pi * 8e6 * L / quality + 1j * w * (li + ground_nh * 1e-9)
    ).astype(complex)
    if mutual is not None:
        a, b, k = mutual
        branches[a, b] = branches[b, a] = 1j * w * k * np.sqrt(li[a] * li[b])
    y = linalg.inv(branches)
    shunt = np.broadcast_to(np.asarray(stray_pf), (85,)) * 1e-12
    if srf_mhz is not None:
        shunt = shunt + 1 / ((2 * np.pi * srf_mhz * 1e6) ** 2 * li)
    yc = 1 / (cap_esr + 1j * w * (cap_esl_nh + trace_nh) * 1e-9 + 1 / (1j * w * C))
    for a, b in g["edges"]:
        y[a, a] += yc
        y[b, b] += yc
        y[a, b] -= yc
        y[b, a] -= yc
    y[np.diag_indices(85)] += 1j * w * shunt + yc * np.array(g["missing"])
    return y


def source_visibility(
    g: dict[str, Any], modes: FloatArray, clusters: list[list[int]]
) -> dict[str, Any]:
    """A degenerate cluster needs independent sources, not merely a visible peak."""
    result = {}
    for label, sources in [("center_only", [0]), ("four_selected", g["sources"])]:
        entries = []
        for cluster in clusters:
            singular = linalg.svdvals(modes[np.ix_(sources, cluster)])
            rank = int(np.sum(singular > 1e-8))
            entries.append(
                {
                    "multiplicity": len(cluster),
                    "source_rank": rank,
                    "full_rank": rank == len(cluster),
                    "minimum_source_singular_value": (
                        float(singular[-1]) if rank == len(cluster) else 0.0
                    ),
                }
            )
        result[label] = {
            "first6_clusters_full_rank": all(e["full_rank"] for e in entries[:6]),
            "full_rank_clusters": sum(e["full_rank"] for e in entries),
            "clusters": entries,
        }
    return result


def noise_deembedding(g: dict[str, Any]) -> tuple[list[dict[str, Any]], float]:
    """Simulate separate diagonal captures, the scan and calibration uncertainty."""
    source = g["sources"][1]
    noise = np.random.default_rng(20261006 + ord(g["name"]))
    records = []
    exact_error = 0.0
    for hz in [3e6, 4.5e6, 6e6, 8e6, 10e6, 15e6]:
        z = linalg.inv(network_admittance(g, hz))
        yp = 1 / 10e6 + 2j * np.pi * hz * 10e-12
        diagonal = np.diag(z)
        loaded = z[:, source] / (1 + yp * diagonal)
        recovered = loaded * (1 + yp * diagonal)
        exact_error = max(
            exact_error,
            float(
                np.linalg.norm(recovered - z[:, source]) / np.linalg.norm(z[:, source])
            ),
        )
        vs = z[source, source] - yp * z[source, :] * z[:, source] / (1 + yp * diagonal)
        # Z here has no source probe. Add that shunt and the generator/sense load.
        current = 0.1 / (1050 + vs + 1050 * yp * vs)
        source_v, scan_v = current * vs, current * loaded
        pre_v = source_v + 1000 * current * (1 + yp * vs)
        diagonal_i = 0.1 / (1050 + diagonal + 1050 * yp * diagonal)
        diagonal_v = diagonal_i * diagonal
        diagonal_pre = diagonal_v + 1000 * diagonal_i * (1 + yp * diagonal)
        for sample_noise_uv, gain_error, phase_deg, cp_error in [
            (0, 0, 0, 0),
            (100, 0, 0, 0),
            (500, 0, 0, 0),
            (2000, 0, 0, 0),
            (100, 0.005, 0.5, 2),
        ]:
            errors = []
            for _ in range(40):
                # Twenty-cycle, 1 GS/s acquisition. White sample noise is a sensitivity
                # parameter, not a claim about the Rigol. I/Q fitting averages samples.
                samples = int(20 * 1e9 / hz)
                sigma = sample_noise_uv * 1e-6 * np.sqrt(2 / samples)
                gain = (1 + noise.normal(0, gain_error, 3)) * np.exp(
                    1j * np.deg2rad(noise.normal(0, phase_deg, 3))
                )

                def capture(
                    values: NDArray[np.complex128],
                    channel: int,
                    gains: NDArray[np.complex128] = gain,
                    deviation: float = float(sigma),
                ) -> NDArray[np.complex128]:
                    return values * gains[channel] + deviation * (
                        noise.normal(size=85) + 1j * noise.normal(size=85)
                    )

                vin, vsource, vscan = (
                    capture(pre_v, 0),
                    capture(source_v, 1),
                    capture(scan_v, 2),
                )
                assumed_y = (
                    1 / 10e6
                    + 2j * np.pi * hz * (10 + noise.normal(0, cp_error)) * 1e-12
                )
                inferred_i = (vin - vsource) / 1000 - assumed_y * vsource
                diagonal_capture = capture(diagonal_v, 1)
                zi = diagonal_capture / (
                    (capture(diagonal_pre, 0) - diagonal_capture) / 1000
                    - assumed_y * diagonal_capture
                )
                transfer = vscan / inferred_i * (1 + assumed_y * zi)
                errors.append(
                    float(
                        np.linalg.norm(transfer - z[:, source])
                        / np.linalg.norm(z[:, source])
                    )
                )
            records.append(
                {
                    "frequency_mhz": hz / 1e6,
                    "white_sample_noise_uv_rms": sample_noise_uv,
                    "channel_gain_sigma": gain_error,
                    "channel_phase_sigma_deg": phase_deg,
                    "probe_capacitance_sigma_pf": cp_error,
                    "transfer_relative_error_p95": float(np.quantile(errors, 0.95)),
                    "scan_nodes_above_10sigma_phasor": int(
                        np.sum(
                            abs(scan_v)
                            > 10
                            * max(
                                sample_noise_uv
                                * 1e-6
                                * np.sqrt(2 / int(20 * 1e9 / hz)),
                                1e-15,
                            )
                        )
                    ),
                    "deembedding_gain_max": float(np.max(abs(1 + yp * diagonal))),
                }
            )
    assert exact_error < 1e-12
    return records, exact_error


def derived_validation(
    g: dict[str, Any], q: FloatArray, v: FloatArray
) -> dict[str, Any]:
    """Independent matrix functions, partial-tail coverage and kernel errors."""
    operator = capacitance(g) / C
    green = (v / q) @ v.T
    error = float(np.linalg.norm(green - linalg.inv(operator)) / np.linalg.norm(green))
    kernel_error = max(
        float(
            np.linalg.norm(
                (v * np.exp(-time * q)) @ v.T - linalg.expm(-time * operator)
            )
        )
        for time in [0.01, 0.1, 1, 10]
    )
    assert error < 1e-10 and kernel_error < 1e-10
    entries = []
    for count in [1, 6, 12, 24, 42, 85]:
        partial_green = (v[:, :count] / q[:count]) @ v[:, :count].T
        for time in [0.1, 1, 5, 10]:
            trace = np.exp(-time * q).sum()
            partial = np.exp(-time * q[:count]).sum()
            bound = (85 - count) * np.exp(-time * q[count - 1])
            assert partial <= trace + 1e-12 and trace <= partial + bound + 1e-12
            entries.append(
                {
                    "modes": count,
                    "time": time,
                    "heat_trace_captured_fraction": float(partial / trace),
                    "heat_tail_bound": float(bound),
                    "green_relative_frobenius_error": float(
                        np.linalg.norm(green - partial_green) / np.linalg.norm(green)
                    ),
                    "zeta_captured_fraction_s123": [
                        float(np.sum(q[:count] ** (-s)) / np.sum(q ** (-s)))
                        for s in [1.0, 2.0, 3.0]
                    ],
                }
            )
    return {
        "inverse_relative_error": error,
        "expm_absolute_error": kernel_error,
        "partial_spectra": entries,
    }


def localization_stress(g: dict[str, Any], trials: int) -> list[dict[str, Any]]:
    freq, modes = eigenpairs(capacitance(g), np.full(85, L))
    clusters = groups(freq)[:6]
    result = []
    for lt, ct, correlation in [
        (0.1, 0.05, 0),
        (0.05, 0.05, 0),
        (0.2, 0.05, 0),
        (0.2, 0.1, 0),
        (0.4, 0.1, 0),
        (0.1, 0.05, 0.1),
    ]:
        rng = np.random.default_rng(20261007 + ord(g["name"]))
        overlap, shift, participation = [], [], []
        for _ in range(trials):
            li = (
                L
                * (1 + rng.uniform(-lt, lt, 85))
                * (1 + rng.uniform(-correlation, correlation))
            )
            edge = C * (1 + rng.uniform(-ct, ct, len(g["edges"])))
            # Sample individual boundary capacitors; their bank average differs from
            # fully correlated bank variation in the original conservative sweep.
            bank = np.array(
                [np.sum(1 + rng.uniform(-ct, ct, n)) * C for n in g["missing"]]
            )
            f, volts = eigenpairs(capacitance(g, edge, bank, 2), li)
            transformed = volts / np.sqrt(li[:, None])
            transformed /= np.linalg.norm(transformed, axis=0)
            a, b = linear_sum_assignment(-abs(modes.T @ transformed) ** 2)
            mapping = dict(zip(a.tolist(), b.tolist()))
            angle, pr = [], []
            for cluster in clusters:
                basis, _ = np.linalg.qr(transformed[:, [mapping[i] for i in cluster]])
                angle.append(
                    float(linalg.svdvals(modes[:, cluster].T @ basis).min() ** 2)
                )
                density = np.sum(basis**2, axis=1) / len(cluster)
                pr.append(float(1 / np.sum(density**2)))
            overlap.append(min(angle))
            participation.append(pr)
            shift.append(float(np.max(abs(f[:12] / freq[:12] - 1))))
        result.append(
            {
                "L_tolerance": lt,
                "C_tolerance": ct,
                "common_L_lot_tolerance": correlation,
                "first6_self_adjoint_subspace_overlap_p05": float(
                    np.quantile(overlap, 0.05)
                ),
                "first12_sorted_frequency_shift_p95": float(np.quantile(shift, 0.95)),
                "cluster_density_participation_p05": np.quantile(
                    participation, 0.05, axis=0
                ).tolist(),
            }
        )
    return result


def loaded_pulse(
    g: dict[str, Any], source: int, scan: int, sigma_ns: float, time: FloatArray
) -> tuple[FloatArray, FloatArray]:
    """Independent continuous-time state model including both probes and source load."""
    cm = capacitance(g, stray_pf=2)
    conductance = np.zeros(85)
    for node in [source, scan]:
        cm[node, node] += 10e-12
        conductance[node] += 1 / 10e6
    conductance[source] += 1 / 1050
    ci = linalg.inv(cm)
    resistance = 2 * np.pi * 8e6 * L / 36
    a = np.block(
        [
            [-ci * conductance[None, :], -ci],
            [np.eye(85) / L, -np.eye(85) * resistance / L],
        ]
    )
    b = np.zeros((170, 1))
    b[:85, 0] = ci[:, source] / 1050
    c = np.eye(170)
    drive = 0.1 * np.exp(-0.5 * ((time - 0.5e-6) / (sigma_ns * 1e-9)) ** 2)
    _, states, _ = signal.lsim((a, b, c, np.zeros((170, 1))), drive, time)
    return states[:, :85], states[:, 85:]


def pulse_tests(g: dict[str, Any]) -> dict[str, Any]:
    name = g["name"]
    time = np.arange(0, 6e-6, 2e-9, dtype=np.float64)
    source, scan = g["sources"][2], g["sources"][1]
    # Match the separately stamped ngspice waveform, whose center is 0.3 us.
    shifted_time = np.arange(0, 6e-6, 2e-9, dtype=np.float64)
    wave, _ = loaded_pulse(g, source, scan, 30, shifted_time)
    spice = np.load(OUT / f"{name}_loaded_pulse.npz")
    aligned = np.column_stack(
        [
            np.interp(time - 0.2e-6, spice["time_s"], spice["voltage_v"][:, i], left=0)
            for i in range(85)
        ]
    )
    error = float(np.max(abs(wave - aligned)) / np.max(abs(wave)))
    assert error < 0.01, error
    points = [complex(*z) for z in g["metric_xy"]]
    records = []
    waves = []
    for src in [source, g["sources"][3]]:
        distance = np.array(
            [
                abs(z - points[src]) if name == "E" else hyperdistance(z, points[src])
                for z in points
            ]
        )
        for width in [30.0, 80.0, 150.0]:
            signal_wave, currents = loaded_pulse(g, src, scan, width, time)
            voltage_weight = signal_wave**2
            spreading = np.sqrt(
                (voltage_weight @ distance**2)
                / np.maximum(voltage_weight.sum(axis=1), 1e-30)
            )
            # Include edge, boundary, stray and both active probe capacitances.
            energy_capacitance = capacitance(g, stray_pf=2)
            for probe_node in [src, scan]:
                energy_capacitance[probe_node, probe_node] += 10e-12
            energy = 0.5 * L * np.sum(currents**2, axis=1) + 0.5 * np.einsum(
                "ti,ij,tj->t", signal_wave, energy_capacitance, signal_wave
            )
            env = abs(signal.hilbert(signal_wave, axis=0))
            threshold = np.maximum(0.2 * env.max(axis=0), 0.0005)
            above = (env > threshold[None, :]) & (time[:, None] >= 0.5e-6)
            first = np.argmax(above, axis=0)
            seen = above.any(axis=0)
            arrival = np.where(seen, time[first], np.nan)
            records.append(
                {
                    "source": src,
                    "pulse_sigma_ns": width,
                    "nodes_above_0_5mv_envelope": int(np.sum(seen)),
                    "peak_voltage_v": float(np.max(abs(signal_wave))),
                    "voltage_weighted_metric_radius_at_1us": float(
                        spreading[np.searchsorted(time, 1e-6)]
                    ),
                    "peak_electrical_energy_j": float(energy.max()),
                    "threshold_arrival_us": [
                        float(x * 1e6) if np.isfinite(x) else None for x in arrival
                    ],
                }
            )
            waves.append(signal_wave)
    # Physically scanned data changes the scan-probe position on every acquisition.
    fixed = waves[0]
    scanned = np.zeros_like(fixed)
    for node in range(85):
        moved, _ = loaded_pulse(g, source, node, 30, time)
        scanned[:, node] = moved[:, node]
    scan_error = float(np.linalg.norm(scanned - fixed) / np.linalg.norm(fixed))
    np.savez_compressed(
        OUT / f"{name}_pulse_tests.npz",
        time_s=time,
        fixed_probe_voltage_v=fixed,
        moving_probe_voltage_v=scanned,
        source=source,
        scan=scan,
    )
    return {
        "independent_ngspice_peak_error": error,
        "moving_probe_relative_waveform_change": scan_error,
        "cases": records,
        "note": "Threshold arrivals and voltage-weighted radius are descriptive finite-disk observables, not causal speeds or an unreflected continuum geodesic front.",
    }


def pole_recovery_tests(
    g: dict[str, Any], include_neighbor: bool = False, quality: float = 36
) -> list[dict[str, Any]]:
    """Exercise actual-network responses, held-out frequencies and spatial residues."""
    cm = capacitance(g, stray_pf=2)
    cv, v = linalg.eigh(cm)
    undamped = 1 / (2 * np.pi * np.sqrt(L * cv))
    clusters = groups(undamped)
    chosen = clusters[4:7] if include_neighbor else clusters[4:6]
    lower = min(float(undamped[c].mean()) for c in chosen)
    upper = max(float(undamped[c].mean()) for c in chosen)
    hz = np.linspace(lower - 0.04e6, upper + 0.04e6, 181)
    resistance = 2 * np.pi * 8e6 * L / quality
    # Full current-normalized response matrix for the selected four sources.
    impedance = np.array(
        [
            (v * (1 / (2j * np.pi * f * cv + 1 / (resistance + 2j * np.pi * f * L))))
            @ v[g["sources"], :].T
            for f in hz
        ]
    )
    expected = np.sort(
        [
            np.sqrt(
                float(undamped[c].mean()) ** 2 - (resistance / (4 * np.pi * L)) ** 2
            )
            for c in chosen
        ]
    )
    reports = []
    rng = np.random.default_rng(20261008 + ord(g["name"]))
    # A limited fit uses a smooth background for all out-of-window modes: residuals
    # and withheld points expose this approximation even with zero measurement noise.
    for relative_noise in [0.0, 0.001, 0.01, 0.05]:
        scale = float(np.sqrt(np.mean(abs(impedance) ** 2)))
        data = impedance + relative_noise * scale * (
            rng.normal(size=impedance.shape) + 1j * rng.normal(size=impedance.shape)
        ) / np.sqrt(2)
        train = np.arange(len(hz)) % 4 != 0
        path = OUT / f'{g["name"]}_recovery_input.npz'
        np.savez_compressed(
            path,
            frequency_hz=hz[train],
            impedance_ohm=data[train].reshape(np.sum(train), -1),
        )
        outcomes = []
        for perturb in [-0.02, 0.02]:
            output = OUT / f'{g["name"]}_recovery_fit.npz'
            fitted = fit_poles(path, (expected / 1e6 + perturb).tolist(), output)
            saved = np.load(output)
            freq = saved["frequency_hz"]
            gamma = saved["damping_hz"]
            residue = saved["residue_ohm_per_s"]
            s = 2j * np.pi * hz[:, None, None]
            poles = 2 * np.pi * (-gamma + 1j * freq)
            prediction = np.sum(
                residue[None] / (s - poles[None, :, None])
                + residue.conj()[None] / (s - poles.conj()[None, :, None]),
                axis=1,
            )
            background = saved["background_coefficients"]
            local_hz = (hz / 1e6 - float(saved["background_center_mhz"])) / float(
                saved["background_scale_mhz"]
            )
            prediction += (
                np.column_stack([local_hz**power for power in range(len(background))])
                @ background
            )
            holdout = float(
                np.linalg.norm(
                    prediction[~train] - data[~train].reshape(np.sum(~train), -1)
                )
                / np.linalg.norm(data[~train])
            )
            spatial = []
            for index in range(len(chosen)):
                cluster = chosen[-1 - index]
                matrix = residue[index].reshape(85, len(g["sources"]))
                basis, _, _ = linalg.svd(
                    np.column_stack([matrix.real, matrix.imag]), full_matrices=False
                )
                spatial.append(
                    float(
                        linalg.svdvals(v[:, cluster].T @ basis[:, : len(cluster)]).min()
                        ** 2
                    )
                )
            outcomes.append(
                {
                    **fitted,
                    "heldout_relative_rms": holdout,
                    "frequency_relative_error_max": float(
                        np.max(abs(freq / expected - 1))
                    ),
                    "spatial_cluster_subspace_overlap_min": min(spatial),
                    "seed_offset_mhz": perturb,
                }
            )
        reports.append(
            {
                "relative_complex_noise_rms": relative_noise,
                "Q_at_8mhz": quality,
                "fitted_clusters": len(chosen),
                "frequency_window_mhz": [float(hz.min() / 1e6), float(hz.max() / 1e6)],
                "expected_damped_frequency_hz": expected.tolist(),
                "fits": outcomes,
            }
        )
    return reports


def nonideal_spice_check(g: dict[str, Any]) -> float:
    """Independently stamp ESR/ESL/SRF and mutual coupling in ngspice."""
    name = g["name"]
    source, scan = g["sources"][1:3]
    pair = tuple(g["edges"][0])
    r = 2 * np.pi * 8e6 * L / 36
    lines = [
        f"* {name}: nonideal independent branch stamps",
        "Vsource vin 0 AC 1",
        "Rgen vin pre 50",
        f"Rin pre n{source} 1000",
    ]
    for i in range(85):
        lines.extend(
            [
                f"LI{i} n{i} lr{i} {L}",
                f"LG{i} lr{i} rr{i} 10n",
                f"RI{i} rr{i} 0 {r}",
                f"CSRF{i} n{i} 0 {1/((2*np.pi*95e6)**2*L)}",
                f"CSTRAY{i} n{i} 0 2p",
            ]
        )
        for j in range(g["missing"][i]):
            suffix = f"B{i}_{j}"
            lines.extend(
                [
                    f"R{suffix} n{i} b{suffix} .2",
                    f"L{suffix} b{suffix} c{suffix} 21n",
                    f"C{suffix} c{suffix} 0 330p",
                ]
            )
    lines.append(f"KPAIR LI{pair[0]} LI{pair[1]} .1")
    for k, (a, b) in enumerate(g["edges"]):
        lines.extend(
            [f"RE{k} n{a} re{k} .2", f"LE{k} re{k} ce{k} 21n", f"CE{k} ce{k} n{b} 330p"]
        )
    for tag, node in [("s", source), ("m", scan)]:
        lines.extend([f"CP{tag} n{node} 0 10p", f"RP{tag} n{node} 0 10meg"])
    output = OUT / f"{name}_nonideal_spice.txt"
    deck = OUT / f"{name}_nonideal_spice.cir"
    lines.extend(
        [
            ".control",
            "set wr_singlescale",
            "set wr_vecnames",
            "ac lin 41 2meg 17meg",
            f"wrdata {output} v(n{source}) v(n{scan})",
            "quit",
            ".endc",
            ".end",
        ]
    )
    deck.write_text("\n".join(lines) + "\n")
    subprocess.run(["ngspice", "-b", str(deck)], check=True, capture_output=True)
    data = np.loadtxt(output, skiprows=1)
    actual = data[:, 1::2] + 1j * data[:, 2::2]
    expected = []
    for hz in data[:, 0]:
        y = network_admittance(
            g,
            hz,
            srf_mhz=95,
            cap_esr=0.2,
            cap_esl_nh=1,
            trace_nh=20,
            ground_nh=10,
            mutual=(*pair, 0.1),
        )
        yp = 1 / 10e6 + 2j * np.pi * hz * 10e-12
        y[source, source] += yp + 1 / 1050
        y[scan, scan] += yp
        drive = np.zeros(85, dtype=complex)
        drive[source] = 1 / 1050
        expected.append(linalg.solve(y, drive)[[source, scan]])
    error = float(np.linalg.norm(actual - np.array(expected)) / np.linalg.norm(actual))
    assert error < 1e-6, error
    return error


def fault_tests(g: dict[str, Any]) -> list[dict[str, Any]]:
    """Response signatures of reversible defects and representative assembly faults."""
    _, v = eigenpairs(capacitance(g), np.full(85, L))
    selected = select_defects(g, v)
    sources = g["sources"]
    sweep = np.linspace(2e6, 17e6, 121)
    results = []
    cases = [(f"edge_{index}_open", ("edge", index, 0.0)) for index in selected]
    cases += [
        (f"edge_{selected[0]}_short_1ohm", ("short", selected[0], 0.0)),
        ("boundary_one_cap_removed", ("boundary", sources[2], 0.0)),
        ("inductor_source_open", ("inductor_open", sources[1], 0.0)),
        ("inductor_source_short_1ohm", ("inductor_short", sources[1], 0.0)),
    ]
    baseline = []
    for hz in sweep:
        z = linalg.inv(network_admittance(g, float(hz)))
        baseline.append(z)
    for label, (kind, index, _) in cases:
        signatures = []
        for src in sources:
            original = []
            changed = []
            for hz, z in zip(sweep, baseline):
                y = network_admittance(g, float(hz))
                w = 2 * np.pi * hz
                if kind in ["edge", "short"]:
                    a, b = g["edges"][index]
                    delta = -1j * w * C + (1.0 if kind == "short" else 0.0)
                    y[a, a] += delta
                    y[b, b] += delta
                    y[a, b] -= delta
                    y[b, a] -= delta
                elif kind == "boundary":
                    y[index, index] -= 1j * w * C
                else:
                    y[index, index] -= 1 / (2 * np.pi * 8e6 * L / 36 + 1j * w * L)
                    y[index, index] += 1.0 if kind == "inductor_short" else 0.0
                zs = z[:, src]
                zc = linalg.solve(y, np.eye(85)[:, src])
                # Two 10 pF probes: source and one common boundary monitor.
                monitor = sources[2]
                yp = 1 / 10e6 + 1j * w * 10e-12

                def source_loaded_column(
                    matrix: NDArray[np.complex128],
                    column: NDArray[np.complex128],
                    probe_y: complex = yp,
                    monitor_node: int = monitor,
                    source_node: int = src,
                ) -> NDArray[np.complex128]:
                    scanloaded = column - probe_y * matrix[:, monitor_node] * column[
                        monitor_node
                    ] / (1 + probe_y * matrix[monitor_node, monitor_node])
                    return (
                        0.1
                        * scanloaded
                        / (1050 + (1 + 1050 * probe_y) * scanloaded[source_node])
                    )

                original.append(source_loaded_column(z, zs)[[src, monitor]])
                zchange = linalg.inv(y)
                changed.append(source_loaded_column(zchange, zc)[[src, monitor]])
            old = np.array(original)
            new = np.array(changed)
            difference = new - old
            signatures.append(
                {
                    "source": src,
                    "two_channel_relative_response_change": float(
                        np.linalg.norm(difference) / np.linalg.norm(old)
                    ),
                    "largest_voltage_change_mv": float(np.max(abs(difference)) * 1000),
                }
            )
        results.append({"fault": label, "signatures": signatures})
    return results


def simulation_tests(gs: list[dict[str, Any]], trials: int) -> dict[str, Any]:
    """Complete the design's computational measurement tests without altering CAD."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    report: dict[str, Any] = {
        "seed": 20261006,
        "trials": trials,
        "geometries": {},
        "assumptions": "Noise/calibration/ESR/ESL/trace/mutual values are sensitivity parameters, not measured equipment or extracted PCB models.",
    }
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for row, g in enumerate(gs):
        print(
            f"{g['name']}: visibility, parasitics and measurement reconstruction",
            flush=True,
        )
        freq, modes = eigenpairs(capacitance(g), np.full(85, L))
        q = 1 / ((2 * np.pi * freq) ** 2 * L * C)
        source = g["sources"][1]
        scan = g["sources"][2]
        sweep = np.linspace(2e6, 17e6, 301, dtype=np.float64)
        pair = tuple(g["edges"][0])
        cases: list[tuple[str, dict[str, Any]]] = [
            ("reference", {}),
            ("SRF95", {"srf_mhz": 95}),
            ("TDK_SRF120", {"srf_mhz": 120}),
            ("FH_Q28_SRF75", {"quality": 28, "srf_mhz": 75}),
            ("SRF50", {"srf_mhz": 50}),
            ("Q20", {"quality": 20, "srf_mhz": 50}),
            ("Q10", {"quality": 10, "srf_mhz": 35}),
            (
                "ESR_ESL_trace",
                {
                    "srf_mhz": 95,
                    "cap_esr": 0.2,
                    "cap_esl_nh": 1.0,
                    "trace_nh": 20.0,
                    "ground_nh": 10.0,
                },
            ),
            ("mutual_positive", {"mutual": (*pair, 0.1)}),
            ("mutual_negative", {"mutual": (*pair, -0.1)}),
            (
                "uneven_stray",
                {"stray_pf": np.random.default_rng(44 + row).uniform(0, 10, 85)},
            ),
        ]
        responses = []
        parasitics = []
        for label, kwargs in cases:
            voltages = []
            for hz in sweep:
                y = network_admittance(g, float(hz), **kwargs)
                assert np.linalg.eigvalsh(y.real).min() > -1e-9
                yp = 1 / 10e6 + 2j * np.pi * hz * 10e-12
                y[source, source] += yp + 1 / 1050
                y[scan, scan] += yp
                current = np.zeros(85, dtype=complex)
                current[source] = 1 / 1050
                voltages.append(linalg.solve(y, current))
            response_v = np.array(voltages)
            responses.append(response_v)
            baseline = responses[0]

            def source_amplitude(
                hz: float,
                parameters: dict[str, Any] = kwargs,
                geometry: dict[str, Any] = g,
                source_node: int = source,
                scan_node: int = scan,
            ) -> float:
                y = network_admittance(geometry, hz, **parameters)
                probe_y = 1 / 10e6 + 2j * np.pi * hz * 10e-12
                y[source_node, source_node] += probe_y + 1 / 1050
                y[scan_node, scan_node] += probe_y
                drive = np.zeros(85, dtype=complex)
                drive[source_node] = 1 / 1050
                return float(abs(linalg.solve(y, drive)[source_node]))

            window_indices = np.flatnonzero(
                (sweep >= freq[0] * 0.9) & (sweep <= freq[0] * 1.02)
            )
            best = int(
                window_indices[np.argmax(abs(response_v[window_indices, source]))]
            )
            peak = minimize_scalar(
                lambda hz: -source_amplitude(float(hz)),
                bounds=(
                    max(freq[0] * 0.9, sweep[best] - 50e3),
                    min(freq[0] * 1.02, sweep[best] + 50e3),
                ),
                method="bounded",
                options={"xatol": 100},
            )
            parasitics.append(
                {
                    "scenario": label,
                    "highest_loaded_peak_mhz": float(peak.x / 1e6),
                    "relative_response_rms_change": float(
                        np.linalg.norm(response_v - baseline) / np.linalg.norm(baseline)
                    ),
                    "peak_source_v_per_generator_v": float(
                        np.max(abs(response_v[:, source]))
                    ),
                    "peak_scan_v_per_generator_v": float(
                        np.max(abs(response_v[:, scan]))
                    ),
                }
            )
            axes[row, 0].plot(sweep / 1e6, abs(response_v[:, source]), label=label)
        expected = response(g, sweep, source, [source, scan], 10, 36)
        baseline_error = float(
            np.linalg.norm(responses[0] - expected) / np.linalg.norm(expected)
        )
        assert baseline_error < 1e-10
        measurement, exact = noise_deembedding(g)
        print(
            f"{g['name']}: disorder, pole recovery and 85-node moving-probe pulses",
            flush=True,
        )
        localization = localization_stress(g, trials)
        fits = pole_recovery_tests(g, g["name"] == "H")
        cheap_fits = pole_recovery_tests(g, g["name"] == "H", quality=20)
        pulses = pulse_tests(g)
        derivation = derived_validation(g, q, modes)
        visibility = source_visibility(g, modes, groups(freq))
        report["geometries"][g["name"]] = {
            "source_visibility": visibility,
            "parasitic_scenarios": parasitics,
            "baseline_matrix_relative_error": baseline_error,
            "nonideal_ngspice_relative_error": nonideal_spice_check(g),
            "assembly_faults": fault_tests(g),
            "exact_scan_deembedding_error": exact,
            "measurement_noise": measurement,
            "disorder_localization": localization,
            "actual_network_pole_recovery": fits,
            "Q20_pole_recovery": cheap_fits,
            "pulses": pulses,
            "finite_graph_functions": derivation,
        }
        np.savez_compressed(
            OUT / f'{g["name"]}_nonideal_ac.npz',
            frequency_hz=sweep,
            scenarios=[label for label, _ in cases],
            voltage_v=np.array(responses),
        )
        axes[row, 0].set(
            title=f'{g["name"]}: nonideal source responses',
            xlabel="MHz",
            ylabel="V / V generator",
        )
        axes[row, 0].legend(fontsize=6)
        for noise in [100, 500, 2000]:
            subset = [
                entry
                for entry in measurement
                if entry["white_sample_noise_uv_rms"] == noise
                and entry["channel_gain_sigma"] == 0
            ]
            axes[row, 1].semilogy(
                [e["frequency_mhz"] for e in subset],
                [max(e["transfer_relative_error_p95"], 1e-12) for e in subset],
                label=f"{noise} µV white sample noise",
            )
        axes[row, 1].set(
            title="Noisy current reconstruction + scan deembedding",
            xlabel="MHz",
            ylabel="95% map relative error",
        )
        axes[row, 1].legend(fontsize=7)
        for count in [6, 12, 24, 42]:
            subset = [
                entry
                for entry in derivation["partial_spectra"]
                if entry["modes"] == count
            ]
            axes[row, 2].plot(
                [e["time"] for e in subset],
                [e["heat_trace_captured_fraction"] for e in subset],
                label=f"{count} modes",
            )
        axes[row, 2].set(
            title="Finite heat trace recovered from low modes",
            xlabel="dimensionless time",
            ylabel="captured fraction",
        )
        axes[row, 2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "test_simulations.png", dpi=160)
    plt.close(fig)
    write_json(OUT / "test_simulations.json", report)
    return report


def components(gs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result = []
    counts = {"L": 0, "C": 0, "TP": 0, "R": 0, "J": 0}

    def add(
        kind: str,
        value: str,
        nets: list[str],
        footprint: str,
        xy: list[float],
        angle: float = 0,
        role: str = "",
        dnp: bool = False,
    ) -> None:
        counts[kind] += 1
        ref = f"{kind}{counts[kind]}"
        result.append(
            {
                "ref": ref,
                "value": value,
                "nets": nets,
                "footprint": footprint,
                "xy": xy,
                "angle": angle,
                "role": role,
                "dnp": dnp,
                "uuid": uid(ref),
            }
        )

    for g in gs:
        name = g["name"]
        xy = np.array(g["pcb_xy_mm"])
        center = xy[0]
        for i, pos in enumerate(xy):
            delta = pos - center
            radial = delta / np.linalg.norm(delta) if i else np.array([1.0, 0.0])
            tangent = np.array([-radial[1], radial[0]])
            angle = -math.degrees(math.atan2(radial[1], radial[0]))
            net = f"{name}{i:02}"
            add(
                "L",
                "1u",
                [net, "GND"],
                FP_L,
                (pos + radial * 2.0).tolist(),
                angle,
                role=f"{net} inductor; Abracon AIML-0805-1R0K-T 10% Q>=45 at 10MHz",
            )
            add(
                "TP",
                net,
                [net],
                FP_T,
                (pos - radial * 1.3).tolist(),
                role=f"node {i}, shell {g['shells'][i]}",
            )
            add(
                "TP",
                "GND",
                ["GND"],
                FP_T,
                (pos + radial * 4.8).tolist(),
                role=f"{net} ground spring",
            )
            for m in range(g["missing"][i]):
                add(
                    "C",
                    "330p",
                    [net, "GND"],
                    FP_C,
                    (
                        pos
                        + radial * (1.3 + 2.6 * (m // 2))
                        + tangent * (3.0 if m % 2 == 0 else -3.0)
                    ).tolist(),
                    angle,
                    role=f"{net} boundary coupling {m+1}/{g['missing'][i]}",
                )
        for k, (a, b) in enumerate(g["edges"]):
            delta = xy[b] - xy[a]
            angle = -math.degrees(math.atan2(delta[1], delta[0]))
            pos = 0.5 * (xy[a] + xy[b])
            role = f"{name} edge {a}-{b}"
            if k in g["defects"]:
                role += f" defect D{g['defects'].index(k)+1}: replace 165p/open or parallel 330p"
            add(
                "C",
                "330p",
                [f"{name}{a:02}", f"{name}{b:02}"],
                FP_C,
                pos.tolist(),
                angle,
                role=role,
            )
            if k in g["defects"]:
                normal = np.array([-delta[1], delta[0]]) / np.linalg.norm(delta)
                add(
                    "C",
                    "330p",
                    [f"{name}{a:02}", f"{name}{b:02}"],
                    FP_C,
                    (pos + normal * 2).tolist(),
                    angle,
                    role=role + " parallel",
                    dnp=True,
                )
        # Compact manual selector: isolated node pads; short removable patch from drive output.
        x = 20 if name == "E" else 94
        add(
            "J",
            "SOURCE",
            [f"{name}_IN", "GND"],
            "Connector_PinHeader_2.54mm:PinHeader_1x02_P2.54mm_Vertical",
            [x, 7],
        )
        add(
            "R",
            "1k",
            [f"{name}_IN", f"{name}_DRIVE"],
            "Resistor_SMD:R_0603_1608Metric",
            [x + 6, 7],
            role="injection/sense 1%, low parasitic",
        )
        add(
            "TP",
            f"{name}_DRIVE",
            [f"{name}_DRIVE"],
            FP_T,
            [x + 9, 7],
            role="patch to exactly one selected node",
        )
        add(
            "TP",
            f"{name}_IN",
            [f"{name}_IN"],
            FP_T,
            [x + 3, 7],
            role="CH1: actual voltage before sense resistor",
        )
    return result


def optimize_placement(parts: list[dict[str, Any]]) -> None:
    """Pack oriented courtyard rectangles using separating-axis gradients."""
    anchors = np.array([p["xy"] for p in parts], dtype=float)
    for p in parts:
        p["angle"] = round(p["angle"] / 45) * 45
    angle = np.deg2rad([-p["angle"] for p in parts])
    u = np.column_stack([np.cos(angle), np.sin(angle)])
    v = np.column_stack([-np.sin(angle), np.cos(angle)])
    half = np.array(
        [
            (
                [1.75, 1.03]
                if p["ref"].startswith("L")
                else (
                    [1.85, 3.15]
                    if p["ref"].startswith("J")
                    else [1.08, 1.08] if p["ref"].startswith("TP") else [1.52, 0.82]
                )
            )
            for p in parts
        ]
    )
    offset = np.array(
        [[0, 1.27] if p["ref"].startswith("J") else [0, 0] for p in parts]
    )
    extent = abs(u) * half[:, 0, None] + abs(v) * half[:, 1, None]

    def objective(flat: FloatArray) -> tuple[float, FloatArray]:
        xy = flat.reshape(-1, 2)
        delta = xy - anchors
        gradient = 0.0002 * delta
        cost = 0.0001 * np.sum(delta**2)
        pairs = np.array(sorted(cKDTree(xy + offset).query_pairs(9)), dtype=int)
        if len(pairs):
            a, b = pairs.T
            d = xy[a] + offset[a] - xy[b] - offset[b]
            axes = np.stack([u[a], v[a], u[b], v[b]], axis=1)
            projection = np.einsum("pi,pki->pk", d, axes)
            support_a = (
                abs(np.einsum("pi,pki->pk", u[a], axes)) * half[a, 0, None]
                + abs(np.einsum("pi,pki->pk", v[a], axes)) * half[a, 1, None]
            )
            support_b = (
                abs(np.einsum("pi,pki->pk", u[b], axes)) * half[b, 0, None]
                + abs(np.einsum("pi,pki->pk", v[b], axes)) * half[b, 1, None]
            )
            penetration = support_a + support_b + 0.12 - abs(projection)
            selected = np.argmin(penetration, axis=1)
            rows = np.arange(len(pairs))
            overlap = np.maximum(0, penetration[rows, selected])
            cost += np.sum(overlap**2)
            force = (
                -2
                * overlap[:, None]
                * np.where(projection[rows, selected, None] >= 0, 1.0, -1.0)
                * axes[rows, selected]
            )
            np.add.at(gradient, a, force)
            np.add.at(gradient, b, -force)
        return float(cost), gradient.ravel()

    cache = OUT / "placement_cache.json"
    key = hashlib.sha256(
        json.dumps(
            {
                "version": 4,
                "anchors": anchors.tolist(),
                "half": half.tolist(),
                "angle": angle.tolist(),
            },
            sort_keys=True,
        ).encode()
    ).hexdigest()
    if cache.exists():
        saved = json.loads(cache.read_text())
        if saved["key"] == key:
            for part, xy in zip(parts, saved["xy"]):
                part["xy"] = xy
            return
    bounds = [
        (
            float(extent[i, axis] + 1 - offset[i, axis]),
            float((149 if axis == 0 else 117) - extent[i, axis] - offset[i, axis]),
        )
        for i in range(len(parts))
        for axis in range(2)
    ]
    result = minimize(
        objective,
        anchors.ravel(),
        jac=True,
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": 1800, "ftol": 1e-12, "gtol": 1e-6},
    )
    for p, xy in zip(parts, result.x.reshape(-1, 2)):
        p["xy"] = xy.tolist()
    write_json(cache, {"key": key, "xy": result.x.reshape(-1, 2).tolist()})
    write_json(
        OUT / "placement.json",
        {
            "converged": bool(result.success),
            "max_displacement_mm": float(
                np.linalg.norm(result.x.reshape(-1, 2) - anchors, axis=1).max()
            ),
            "objective": float(result.fun),
            "courtyard_model": "oriented rectangles with 0.12mm separation",
        },
    )


def spice_check(gs: list[dict[str, Any]]) -> None:
    """Compare independent ngspice AC and transient solutions to matrix/modal models."""
    results = {}
    for g in gs:
        name = g["name"]
        source, scan = g["sources"][1:3]
        base = [f"* {name}: independently stamped passive network"]
        r = 2 * np.pi * 8e6 * L / 36
        for i in range(85):
            base += [f"L{i} n{i} l{i} 1u", f"R{i} l{i} 0 {r}", f"CS{i} n{i} 0 2p"]
            if g["missing"][i]:
                base.append(f"CB{i} n{i} 0 {330*g['missing'][i]}p")
        for k, (a, b) in enumerate(g["edges"]):
            base.append(f"CE{k} n{a} n{b} 330p")
        ac_path = OUT / f"{name}_ac.txt"
        tran_path = OUT / f"{name}_tran.txt"
        ac = base + [
            "Vsrc vin 0 DC 0 AC 1",
            "Rout vin sense 50",
            f"Rsense sense n{source} 1000",
        ]
        for node in [source, scan]:
            ac += [f"CP{node} n{node} 0 10p", f"RP{node} n{node} 0 10meg"]
        ac += [
            ".control",
            "set wr_singlescale",
            "set wr_vecnames",
            "ac lin 41 1meg 18meg",
            f"wrdata {ac_path} v(n{source}) v(n{scan})",
            "quit",
            ".endc",
            ".end",
        ]
        ac_deck = OUT / f"{name}_ac.cir"
        ac_deck.write_text("\n".join(ac) + "\n")
        subprocess.run(["ngspice", "-b", str(ac_deck)], check=True, capture_output=True)
        data = np.loadtxt(ac_path, skiprows=1)
        actual = data[:, 1::2] + 1j * data[:, 2::2]
        expected = response(g, data[:, 0], source, [source, scan], 10, 36)[
            :, [source, scan]
        ]
        ac_error = float(np.max(abs(actual - expected)) / np.max(abs(expected)))
        assert ac_error < 1e-6, ac_error
        tran = base + [
            f"Bpulse 0 n{scan} I=0.0001*exp(-0.5*((time-0.0000003)/0.00000003)^2)",
            ".control",
            "set wr_singlescale",
            "set wr_vecnames",
            "tran 2n 6u 0 1n",
            f"wrdata {tran_path} v(n{source}) v(n{scan})",
            "quit",
            ".endc",
            ".end",
        ]
        tran_deck = OUT / f"{name}_tran.cir"
        tran_deck.write_text("\n".join(tran) + "\n")
        subprocess.run(
            ["ngspice", "-b", str(tran_deck)], check=True, capture_output=True
        )
        data = np.loadtxt(tran_path, skiprows=1)
        modal = np.load(OUT / f"{name}_pulse.npz")
        expected_tran = np.column_stack(
            [
                np.interp(data[:, 0], modal["time_s"], modal["voltage_v"][:, i])
                for i in [source, scan]
            ]
        )
        tran_error = float(
            np.max(abs(data[:, 1:] - expected_tran)) / np.max(abs(expected_tran))
        )
        assert tran_error < 0.02, tran_error
        results[name] = {
            "ac_relative_error": ac_error,
            "transient_peak_normalized_error": tran_error,
        }
        loaded_path = OUT / f"{name}_loaded_pulse.txt"
        loaded = base + [
            "Bpulse vin 0 V=0.1*exp(-0.5*((time-0.0000003)/0.00000003)^2)",
            "Rgen vin pre 50",
            f"Rsense pre n{scan} 1000",
            f"CPs n{scan} 0 10p",
            f"CPm n{source} 0 10p",
            f"RPs n{scan} 0 10meg",
            f"RPm n{source} 0 10meg",
            ".control",
            "set wr_singlescale",
            "set wr_vecnames",
            "tran 2n 6u 0 1n",
            f"wrdata {loaded_path} " + " ".join(f"v(n{i})" for i in range(85)),
            "quit",
            ".endc",
            ".end",
        ]
        deck = OUT / f"{name}_loaded_pulse.cir"
        deck.write_text("\n".join(loaded) + "\n")
        subprocess.run(["ngspice", "-b", str(deck)], check=True, capture_output=True)
        raw = np.loadtxt(loaded_path, skiprows=1)
        uniform = np.arange(0, 6e-6, 2e-9)
        wave = np.column_stack(
            [np.interp(uniform, raw[:, 0], raw[:, i + 1]) for i in range(85)]
        )
        np.savez_compressed(
            OUT / f"{name}_loaded_pulse.npz",
            time_s=uniform,
            voltage_v=wave,
            source=scan,
        )

    write_json(OUT / "spice_validation.json", results)


def symbol_definition(library: str, name: str) -> str:
    text = (Path("/usr/share/kicad/symbols") / f"{library}.kicad_sym").read_text()
    start = text.index(f'(symbol "{name}"')
    level = 0
    quoted = False
    escaped = False
    for end in range(start, len(text)):
        ch = text[end]
        if ch == '"' and not escaped:
            quoted = not quoted
        if not quoted:
            level += (ch == "(") - (ch == ")")
            if level == 0:
                return text[start : end + 1].replace(
                    f'(symbol "{name}"', f'(symbol "{library}:{name}"', 1
                )
        escaped = ch == "\\" and not escaped
    raise ValueError(name)


def schematic(parts: list[dict[str, Any]]) -> None:
    rootid = uid("schematic")
    definitions = [
        ("Device", "C"),
        ("Device", "L"),
        ("Device", "R"),
        ("Connector", "TestPoint"),
        ("Connector_Generic", "Conn_01x02"),
    ]
    lines = [
        f'(kicad_sch (version 20250114) (generator "eeschema") (uuid {rootid})',
        '(paper "A0")',
        '(title_block (title "Geometric Electricity: dual 85-node LC experiment") (rev "DRAFT-1"))',
        "(lib_symbols " + "\n".join(symbol_definition(*d) for d in definitions) + ")",
    ]
    for i, p in enumerate(parts):
        x = 15 + (i % 44) * 26
        y = 30 + (i // 44) * 24
        kind = (
            "TestPoint"
            if p["ref"].startswith("TP")
            else "Conn_01x02" if p["ref"].startswith("J") else p["ref"][0]
        )
        lib = (
            "Connector"
            if kind == "TestPoint"
            else "Connector_Generic" if kind == "Conn_01x02" else "Device"
        )

        def prop(
            k: str, v: str, dy: float, hide: bool = False, x: float = x, y: float = y
        ) -> str:
            return f'(property {json.dumps(k)} {json.dumps(v)} (at {x+2.54} {y+dy} 0) (effects (font (size 1 1)) (justify left) {"(hide yes)" if hide else ""}))'

        lines.append(
            f'(symbol (lib_id "{lib}:{kind}") (at {x} {y} 0) (unit 1) (in_bom {"no" if kind=="TestPoint" else "yes"}) (on_board yes) (dnp {"yes" if p["dnp"] else "no"}) (uuid {p["uuid"]}) '
            + prop("Reference", p["ref"], -2.54)
            + prop("Value", p["value"], 0)
            + prop("Footprint", p["footprint"], 0, True)
            + prop("Description", p["role"], 0, True)
            + f'(instances (project "geometric_electricity" (path "/{rootid}" (reference "{p["ref"]}") (unit 1)))))'
        )
        positions = (
            [(x, y)]
            if kind == "TestPoint"
            else (
                [(x - 5.08, y), (x - 5.08, y + 2.54)]
                if kind == "Conn_01x02"
                else [(x, y - 3.81), (x, y + 3.81)]
            )
        )
        for j, (px, py) in enumerate(positions):
            # Local net labels avoid a huge, unreadable tangle of 500 graph wires.
            lines.append(
                f'(global_label "{p["nets"][j]}" (shape passive) (at {px} {py} 0) (effects (font (size .9 .9)) (justify left bottom)) (uuid {uid(p["ref"]+"label"+str(j))}))'
            )
    lines.append(
        '(text "DRAFT: 1uH / 330pF; AIML-0805-1R0K-T; Q>=45 at 10MHz; use 10x probes. All node labels are global within this sheet." (at 15 15 0) (effects (font (size 2 2)) (justify left)) (uuid '
        + uid("title")
        + "))"
    )
    lines.append("(embedded_fonts no))")
    (HW / "geometric_electricity.kicad_sch").write_text("\n".join(lines) + "\n")


def board(parts: list[dict[str, Any]]) -> None:
    import pcbnew as pcb  # type: ignore[import-untyped]

    b = pcb.BOARD()
    b.SetCopperLayerCount(2)
    b.GetDesignSettings().SetBoardThickness(pcb.FromMM(1.6))
    nets = {}
    for name in sorted({net for p in parts for net in p["nets"]}):
        net = pcb.NETINFO_ITEM(b, name)
        b.Add(net)
        nets[name] = net
    footprint_templates: dict[str, Any] = {}
    for p in parts:
        library, name = p["footprint"].split(":")
        library_path = (
            HW / "GeometricElectricity.pretty"
            if library == "GeometricElectricity"
            else Path(f"/usr/share/kicad/footprints/{library}.pretty")
        )
        if p["footprint"] not in footprint_templates:
            template = pcb.FootprintLoad(str(library_path), name)
            if template is None:
                raise ValueError(p["footprint"])
            footprint_templates[p["footprint"]] = template
        fp = footprint_templates[p["footprint"]].Duplicate()
        fp.SetReference(p["ref"])
        fp.SetValue(p["value"])
        fp.SetFPID(pcb.LIB_ID(library, name))
        fp.SetPath(pcb.KIID_PATH(f"/{uid('schematic')}/{p['uuid']}"))
        fp.SetPosition(pcb.VECTOR2I(pcb.FromMM(p["xy"][0]), pcb.FromMM(p["xy"][1])))
        fp.SetOrientationDegrees(p["angle"])
        fp.Reference().SetVisible(False)
        fp.Value().SetVisible(False)
        if p["dnp"]:
            fp.SetAttributes(fp.GetAttributes() | pcb.FP_DNP)
        for pad in fp.Pads():
            pad.SetNet(nets[p["nets"][int(pad.GetNumber()) - 1]])
        b.Add(fp)
    for a, c in [
        ((0, 0), (150, 0)),
        ((150, 0), (150, 130)),
        ((150, 130), (0, 130)),
        ((0, 130), (0, 0)),
    ]:
        line = pcb.PCB_SHAPE()
        line.SetShape(pcb.SHAPE_T_SEGMENT)
        line.SetStart(pcb.VECTOR2I(pcb.FromMM(a[0]), pcb.FromMM(a[1])))
        line.SetEnd(pcb.VECTOR2I(pcb.FromMM(c[0]), pcb.FromMM(c[1])))
        line.SetWidth(pcb.FromMM(0.05))
        line.SetLayer(pcb.Edge_Cuts)
        b.Add(line)
    for label, x, y in [
        ("EUCLIDEAN", 38, 124),
        ("HYPERBOLIC", 112, 124),
        ("DRAFT - NOT FOR FABRICATION", 75, 119),
    ]:
        text = pcb.PCB_TEXT(b)
        text.SetText(label)
        text.SetPosition(pcb.VECTOR2I(pcb.FromMM(x), pcb.FromMM(y)))
        text.SetTextSize(pcb.VECTOR2I(pcb.FromMM(1.1), pcb.FromMM(1.1)))
        text.SetTextThickness(pcb.FromMM(0.15))
        text.SetLayer(pcb.F_SilkS)
        b.Add(text)
    # Re-read to invalidate copied footprint shape caches before geometric checks.
    pcb.SaveBoard(str(HW / "geometric_electricity.kicad_pcb"), b)
    b = pcb.LoadBoard(str(HW / "geometric_electricity.kicad_pcb"))
    # Place node IDs against actual copper and silkscreen bounding boxes.
    occupied = [pad.GetBoundingBox() for fp in b.GetFootprints() for pad in fp.Pads()]
    occupied += [
        item.GetBoundingBox()
        for fp in b.GetFootprints()
        for item in fp.GraphicalItems()
        if item.GetLayer() == pcb.F_SilkS
    ]
    occupied += [
        item.GetBoundingBox()
        for item in b.GetDrawings()
        if item.GetLayer() == pcb.F_SilkS
    ]
    geometry_map = {g["name"]: g for g in graphs()}
    for fp in sorted(b.GetFootprints(), key=lambda fp: fp.GetReference()):
        if not fp.GetReference().startswith("TP") or fp.GetValue() == "GND":
            continue
        label = pcb.PCB_TEXT(b)
        value = fp.GetValue()
        if re.fullmatch(r"[EH]\d{2}", value):
            g = geometry_map[value[0]]
            node = int(value[1:])
            value += "*" if node in g["sources"] else ""
        label.SetText(value)
        label.SetTextSize(pcb.VECTOR2I(pcb.FromMM(0.8), pcb.FromMM(0.8)))
        label.SetTextThickness(pcb.FromMM(0.12))
        label.SetLayer(pcb.F_SilkS)
        for radius in [1.8, 2.2, 2.8, 3.4, 4.2, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0]:
            found = False
            for angle in np.linspace(-np.pi / 2, 3 * np.pi / 2, 24, endpoint=False):
                offset = pcb.VECTOR2I(
                    pcb.FromMM(float(radius * np.cos(angle))),
                    pcb.FromMM(float(radius * np.sin(angle))),
                )
                label.SetPosition(fp.GetPosition() + offset)
                bounds = label.GetBoundingBox()
                bounds.Inflate(pcb.FromMM(0.22))
                within_board = (
                    bounds.GetX() > pcb.FromMM(0.4)
                    and bounds.GetRight() < pcb.FromMM(149.6)
                    and bounds.GetY() > pcb.FromMM(0.4)
                    and bounds.GetBottom() < pcb.FromMM(129.6)
                )
                if within_board and not any(bounds.Intersects(box) for box in occupied):
                    found = True
                    break
            if found:
                break
        if not found:
            raise RuntimeError(f"No readable node label position for {fp.GetValue()}")
        b.Add(label)
        occupied.append(bounds)
    # Short local return vias and an uninterrupted reference plane before routing.
    pads = [pad for fp in b.GetFootprints() for pad in fp.Pads()]

    def rectangles(selected: list[Any]) -> FloatArray:
        boxes = [pad.GetBoundingBox() for pad in selected]
        return np.array(
            [
                [box.GetX(), box.GetY(), box.GetRight(), box.GetBottom()]
                for box in boxes
            ],
            dtype=float,
        )

    signal_boxes = rectangles([pad for pad in pads if pad.GetNetname() != "GND"])
    smd_boxes = rectangles(
        [
            pad
            for pad in pads
            if pad.GetAttribute() == pcb.PAD_ATTRIB_SMD
            and not pad.GetParentFootprint().GetReference().startswith("TP")
        ]
    )

    def intersects(box: Any, boxes: FloatArray) -> bool:
        return bool(
            np.any(
                (boxes[:, 0] <= box.GetRight())
                & (boxes[:, 2] >= box.GetX())
                & (boxes[:, 1] <= box.GetBottom())
                & (boxes[:, 3] >= box.GetY())
            )
        )

    ground_vias = 0
    via_positions: list[tuple[int, int]] = []
    for pad in pads:
        if pad.GetNetname() != "GND":
            continue
        candidate = None
        if pad.GetParentFootprint().GetReference().startswith("TP"):
            candidate = pad.GetPosition()
        else:
            origin = pad.GetPosition()
            for radius in [0.9, 1.2, 1.6, 2.0, 2.5, 3.0]:
                for angle in np.linspace(0, 2 * np.pi, 24, endpoint=False):
                    offset = pcb.VECTOR2I(
                        pcb.FromMM(float(radius * np.cos(angle))),
                        pcb.FromMM(float(radius * np.sin(angle))),
                    )
                    point = origin + offset
                    box = pcb.BOX2I(
                        point - pcb.VECTOR2I(pcb.FromMM(0.65), pcb.FromMM(0.65)),
                        pcb.VECTOR2I(pcb.FromMM(1.3), pcb.FromMM(1.3)),
                    )
                    if any(
                        (point.x - x) ** 2 + (point.y - y) ** 2 < pcb.FromMM(0.65) ** 2
                        for x, y in via_positions
                    ):
                        continue
                    on_smd = intersects(box, smd_boxes) or intersects(box, signal_boxes)
                    if on_smd:
                        continue
                    path_clear = True
                    for fraction in np.linspace(0, 1, 12):
                        spot = pcb.VECTOR2I(
                            int(origin.x + fraction * offset.x),
                            int(origin.y + fraction * offset.y),
                        )
                        area = pcb.BOX2I(
                            spot - pcb.VECTOR2I(pcb.FromMM(0.45), pcb.FromMM(0.45)),
                            pcb.VECTOR2I(pcb.FromMM(0.9), pcb.FromMM(0.9)),
                        )
                        if intersects(area, signal_boxes):
                            path_clear = False
                            break
                    if path_clear:
                        candidate = point
                        break
                if candidate is not None:
                    break
        if candidate is not None and any(
            (candidate.x - x) ** 2 + (candidate.y - y) ** 2 < pcb.FromMM(0.65) ** 2
            for x, y in via_positions
        ):
            candidate = None
        if candidate is None:
            continue  # Router connects remaining return pads to the nearest grounded pad.
        via = pcb.PCB_VIA(b)
        via.SetPosition(candidate)
        via.SetWidth(pcb.FromMM(0.6))
        via.SetDrill(pcb.FromMM(0.3))
        via.SetViaType(pcb.VIATYPE_THROUGH)
        via.SetLayerPair(pcb.F_Cu, pcb.B_Cu)
        via.SetNet(nets["GND"])
        b.Add(via)
        if candidate != pad.GetPosition():
            track = pcb.PCB_TRACK(b)
            track.SetStart(pad.GetPosition())
            track.SetEnd(candidate)
            track.SetWidth(pcb.FromMM(0.25))
            track.SetLayer(pcb.F_Cu)
            track.SetNet(nets["GND"])
            b.Add(track)
        ground_vias += 1
        via_positions.append((candidate.x, candidate.y))
    zone = pcb.ZONE(b)
    zone.SetLayer(pcb.B_Cu)
    zone.SetNet(nets["GND"])
    zone.SetLocalClearance(pcb.FromMM(0.25))
    zone.SetPadConnection(pcb.ZONE_CONNECTION_FULL)
    outline = zone.Outline()
    outline.NewOutline()
    for x, y in [(1, 1), (149, 1), (149, 129), (1, 129)]:
        outline.Append(pcb.FromMM(x), pcb.FromMM(y))
    b.Add(zone)
    front_ground(b)
    # LoadBoard initializes the native rule engine needed by the zone filler.
    pcb.SaveBoard(str(HW / "geometric_electricity.kicad_pcb"), b)
    b = pcb.LoadBoard(str(HW / "geometric_electricity.kicad_pcb"))
    b.BuildConnectivity()
    filler = pcb.ZONE_FILLER(b)
    filler.Fill(b.Zones())
    write_json(
        OUT / "ground_returns.json",
        {
            "local_vias": ground_vias,
            "ground_pads": sum(pad.GetNetname() == "GND" for pad in pads),
            "plane": "B.Cu",
        },
    )
    pcb.SaveBoard(str(HW / "geometric_electricity.kicad_pcb"), b)
    gs = json.loads((HW / "graphs.json").read_text())
    ground_points = np.array(
        [[pcb.ToMM(x), pcb.ToMM(y)] for x, y in via_positions]
        + [p["xy"] for p in parts if p["ref"].startswith("TP") and p["value"] == "GND"]
    )
    for g in gs:
        distances, indices = cKDTree(ground_points).query(np.array(g["pcb_xy"]))
        g["probe_ground_xy"] = ground_points[indices].tolist()
        g["probe_ground_distance_mm"] = distances.tolist()
    write_json(HW / "graphs.json", gs)


def front_ground(board: Any) -> None:
    """Stitch the two-sided reference around signal routes with solid ground returns."""
    import pcbnew as pcb  # type: ignore[import-untyped]

    has_front = False
    for existing in board.Zones():
        existing.SetMinThickness(pcb.FromMM(0.15))
        existing.SetIslandRemovalMode(pcb.ISLAND_REMOVAL_MODE_ALWAYS)
        if existing.GetLayer() == pcb.F_Cu:
            existing.SetPadConnection(pcb.ZONE_CONNECTION_FULL)
            existing.SetLocalClearance(pcb.FromMM(0.25))
            has_front = True
    if has_front:
        return
    zone = pcb.ZONE(board)
    zone.SetLayer(pcb.F_Cu)
    zone.SetNet(board.FindNet("GND"))
    zone.SetLocalClearance(pcb.FromMM(0.25))
    zone.SetMinThickness(pcb.FromMM(0.15))
    zone.SetIslandRemovalMode(pcb.ISLAND_REMOVAL_MODE_ALWAYS)
    zone.SetPadConnection(pcb.ZONE_CONNECTION_FULL)
    zone.SetThermalReliefGap(pcb.FromMM(0.25))
    zone.SetThermalReliefSpokeWidth(pcb.FromMM(0.25))
    outline = zone.Outline()
    outline.NewOutline()
    for x, y in [(1, 1), (149, 1), (149, 129), (1, 129)]:
        outline.Append(pcb.FromMM(x), pcb.FromMM(y))
    board.Add(zone)


def routing_exchange(import_session: Path | None = None) -> None:
    """Exchange a placed board with a Specctra-compatible router; never release it."""
    import pcbnew as pcb  # type: ignore[import-untyped]

    path = HW / "geometric_electricity.kicad_pcb"
    current = pcb.LoadBoard(str(path))
    if import_session is not None:
        if not pcb.ImportSpecctraSES(current, str(import_session)):
            raise RuntimeError("Router session import failed")
        front_ground(current)
        current.BuildConnectivity()
        filler = pcb.ZONE_FILLER(current)
        filler.Fill(current.Zones())
        pcb.SaveBoard(str(path), current)
    else:
        # Route against fixed return copper, then refill the native plane after import.
        # An exported solid plane couples hundreds of return contacts in the router.
        for zone in list(current.Zones()):
            current.Remove(zone)
        dsn = OUT / "routing.dsn"
        if not pcb.ExportSpecctraDSN(current, str(dsn)):
            raise RuntimeError("Router design export failed")
        # Keep the large reference net out of signal routing; KiCad verifies it after refill.
        source = dsn.read_text()
        start = source.index("(class kicad_default ")
        end = source.index("(circuit", start)
        source = (
            source[:start] + re.sub(r"\bGND\b", "", source[start:end]) + source[end:]
        )
        source = (
            source[:start]
            + "(class GND GND (rule (width 200) (clearance 200)))\n    "
            + source[start:]
        )
        source = source.replace("(net GND)(type route)", "(net GND)(type protect)")
        dsn.write_text(source)


def budget(parts: list[dict[str, Any]]) -> dict[str, Any]:
    """Cost one populated board, including spares; fabrication entries are allowances."""
    used_l = sum(p["ref"].startswith("L") and not p["dnp"] for p in parts)
    used_c = sum(p["ref"].startswith("C") and not p["dnp"] for p in parts)
    if used_l > 200 or used_c > 800:
        raise ValueError("BOM exceeds budgeted purchasing quantities")
    allocation = {
        "200_inductors": round(200 * 0.09380, 2),
        "800_capacitors": round(800 * 0.0062, 2),
        "resistors_headers_spares": 5.0,
        "pcbway_hasl_boards_including_shipping": 72.76,
        "stencil_quote": 10.0,
        "component_shipping_allowance": 15.0,
        "consumables_allowance": 15.0,
        "tax_import_allowance": 20.0,
        "contingency": 20.0,
    }
    total = round(sum(allocation.values()), 2)
    assert total < 200
    report = {
        "currency": "USD",
        "limit_exclusive": 200,
        "populated_boards": 1,
        "used_inductors": used_l,
        "used_capacitors": used_c,
        "allocation": allocation,
        "planned_total": total,
        "remaining_to_limit": round(200 - total, 2),
        "landed_quotes_confirmed": False,
        "fabrication_quote": {
            "source": "User-provided PCBWay quote",
            "pcb_shipping_included": True,
            "lead_free_hasl_usd": 72.76,
            "enig_usd": 97.23,
            "stencil_usd": 10.0,
            "extra_stencil_shipping_confirmed": False,
        },
        "enig_total_with_current_bom_usd": round(total + 97.23 - 72.76, 2),
        "purchase_release": False,
        "price_check_date": "2026-10-05",
        "price_sources": {
            "inductors": "https://www.digikey.com/en/products/detail/abracon-llc/AIML-0805-1R0K-T/2662996",
            "capacitors": "https://www.lcsc.com/product-detail/C62784.html",
        },
        "preferred_substitution": {
            "inductor_mpn": "MLF2012A1R0JT000",
            "supplier_code": "C165812",
            "price_source": "https://www.lcsc.com/fr/product-detail/C165812.html",
            "quantity": 200,
            "unit_usd": 0.0381,
            "inductor_total_usd": 7.62,
            "proposed_total_usd": round(total - allocation["200_inductors"] + 7.62, 2),
            "enig_proposed_total_usd": round(
                total - allocation["200_inductors"] + 7.62 + 97.23 - 72.76, 2
            ),
            "part_savings_vs_original_allocation_usd": 16.58,
            "L_tolerance": 0.05,
            "Q_min_at_10mhz": 45,
            "SRF_min_mhz": 120,
            "status": "Recommended purchasing alternate; current CAD/BOM still name Abracon. Supplier stock/landed quote, final land-pattern review and RF sample qualification remain required.",
        },
    }
    write_json(OUT / "budget.json", report)
    return report


def verify(gs: list[dict[str, Any]], parts: list[dict[str, Any]]) -> dict[str, Any]:
    import pcbnew as pcb  # type: ignore[import-untyped]

    b = pcb.LoadBoard(str(HW / "geometric_electricity.kicad_pcb"))
    actual = {fp.GetReference(): fp for fp in b.GetFootprints()}
    assert set(actual) == {p["ref"] for p in parts}
    for p in parts:
        fp = actual[p["ref"]]
        assert {pad.GetNumber(): pad.GetNetname() for pad in fp.Pads()} == {
            str(i + 1): net for i, net in enumerate(p["nets"])
        }, p["ref"]
    # Independent KiCad schematic export, then compare pin nets against board.
    path = OUT / "schematic.net"
    subprocess.run(
        [
            "kicad-cli",
            "sch",
            "export",
            "netlist",
            "--format",
            "kicadxml",
            "-o",
            str(path),
            str(HW / "geometric_electricity.kicad_sch"),
        ],
        check=True,
        capture_output=True,
    )
    tree = ET.parse(path)
    sch_pins = {}
    for net in tree.findall("./nets/net"):
        for node in net.findall("node"):
            sch_pins[(node.attrib["ref"], node.attrib["pin"])] = net.attrib["name"]
    for p in parts:
        for i, net in enumerate(p["nets"]):
            assert sch_pins[(p["ref"], str(i + 1))] == net, (
                p["ref"],
                net,
                sch_pins.get((p["ref"], str(i + 1))),
            )
    graph_checks = {}
    for g in gs:
        prefix = g["name"]
        edges = set()
        boundary = np.zeros(85)
        for p in parts:
            if not p["ref"].startswith("C") or p["dnp"]:
                continue
            nets = [pad.GetNetname() for pad in actual[p["ref"]].Pads()]
            if all(re.fullmatch(prefix + r"\d{2}", net) for net in nets):
                edges.add(tuple(sorted(int(net[1:]) for net in nets)))
            elif "GND" in nets and any(
                re.fullmatch(prefix + r"\d{2}", net) for net in nets
            ):
                boundary_net = next(net for net in nets if net != "GND")
                boundary[int(boundary_net[1:])] += float(p["value"][:-1]) / 330
        assert edges == {tuple(e) for e in g["edges"]}
        assert np.array_equal(boundary, np.array(g["missing"]))
        graph_checks[prefix] = {
            "edges": len(edges),
            "boundary_couplings": int(sum(boundary)),
        }
    drc = OUT / "drc.json"
    subprocess.run(
        [
            "kicad-cli",
            "pcb",
            "drc",
            "--format",
            "json",
            "-o",
            str(drc),
            str(HW / "geometric_electricity.kicad_pcb"),
        ],
        check=True,
        capture_output=True,
    )
    erc = OUT / "erc.json"
    subprocess.run(
        [
            "kicad-cli",
            "sch",
            "erc",
            "--format",
            "json",
            "-o",
            str(erc),
            str(HW / "geometric_electricity.kicad_sch"),
        ],
        check=True,
        capture_output=True,
    )
    d = json.loads(drc.read_text())
    e = json.loads(erc.read_text())
    checks = {
        "schematic_pcb_pin_parity": True,
        "graphs": graph_checks,
        "drc_violations": len(d.get("violations", [])),
        "unconnected": len(d.get("unconnected_items", [])),
        "erc_violations": len(e.get("violations", [])),
        "fabrication_ready": False,
        "budget": budget(parts),
    }
    write_json(OUT / "verification.json", checks)
    return checks


def probe_map() -> None:
    """A printable coordinate map disambiguates IDs in the dense layout."""
    import matplotlib.pyplot as plt

    gs = json.loads((HW / "graphs.json").read_text())
    figure, axes = plt.subplots(1, 2, figsize=(16, 12), constrained_layout=True)
    with (OUT / "probe_map.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "Geometry",
                "Node",
                "Testpoint",
                "X_mm",
                "Y_mm",
                "Ground_X_mm",
                "Ground_Y_mm",
                "Ground_distance_mm",
                "Source",
            ]
        )
        for axis, g in zip(axes, gs):
            xy = np.array(g["pcb_xy"])
            ground = np.array(g["probe_ground_xy"])
            axis.scatter(ground[:, 0], ground[:, 1], s=8, c="gray", marker="+")
            axis.scatter(xy[:, 0], xy[:, 1], s=12, c="tab:blue")
            for i, point in enumerate(xy):
                axis.annotate(
                    f"{g['name']}{i:02d}",
                    point,
                    xytext=(2, 2),
                    textcoords="offset points",
                    fontsize=6,
                )
                writer.writerow(
                    [
                        g["name"],
                        i,
                        g["testpoint_references"][i],
                        *point,
                        *ground[i],
                        g["probe_ground_distance_mm"][i],
                        i in g["sources"],
                    ]
                )
            source_xy = xy[g["sources"]]
            axis.scatter(
                source_xy[:, 0],
                source_xy[:, 1],
                s=60,
                facecolors="none",
                edgecolors="tab:red",
                label="source candidates",
            )
            axis.set(
                xlim=(0 if g["name"] == "E" else 75, 75 if g["name"] == "E" else 150),
                ylim=(130, 0),
                xlabel="Board X (mm)",
                ylabel="Board Y (mm)",
                title=f"{g['name']}: blue = signal pad, gray + = nearest ground contact",
            )
            axis.set_aspect("equal")
            axis.grid(alpha=0.2)
            axis.legend()
    figure.savefig(OUT / "probe_map.pdf")
    plt.close(figure)


def refresh_probe_access() -> None:
    """Map current ground contacts and expose only selected standard return vias."""
    import pcbnew as pcb  # type: ignore[import-untyped]

    path = HW / "geometric_electricity.kicad_pcb"
    current = pcb.LoadBoard(str(path))
    silk = [
        item.GetBoundingBox()
        for item in list(current.GetDrawings())
        + [
            graphic
            for footprint in current.GetFootprints()
            for graphic in footprint.GraphicalItems()
        ]
        if item.GetLayer() == pcb.F_SilkS
    ]
    for box in silk:
        box.Inflate(pcb.FromMM(0.15))
    for via in current.GetTracks():
        if isinstance(via, pcb.PCB_VIA) and via.GetNetname() == "GND":
            via.SetFrontTentingMode(pcb.TENTING_MODE_TENTED)
    contacts = [
        via
        for via in current.GetTracks()
        if isinstance(via, pcb.PCB_VIA)
        and via.GetNetname() == "GND"
        and via.GetWidth(pcb.F_Cu) >= pcb.FromMM(0.6)
        and not any(box.Intersects(via.GetBoundingBox()) for box in silk)
    ] + [
        pad
        for footprint in current.GetFootprints()
        if footprint.GetReference().startswith("TP")
        for pad in footprint.Pads()
        if pad.GetNetname() == "GND"
    ]
    points = np.array(
        [
            [pcb.ToMM(item.GetPosition().x), pcb.ToMM(item.GetPosition().y)]
            for item in contacts
        ]
    )
    gs = json.loads((HW / "graphs.json").read_text())
    for g in gs:
        distances, indices = cKDTree(points).query(np.array(g["pcb_xy"]))
        g["probe_ground_xy"] = points[indices].tolist()
        g["probe_ground_distance_mm"] = distances.tolist()
        for index in indices:
            item = contacts[int(index)]
            if isinstance(item, pcb.PCB_VIA):
                item.SetFrontTentingMode(pcb.TENTING_MODE_NOT_TENTED)
    pcb.SaveBoard(str(path), current)
    write_json(HW / "graphs.json", gs)


def export_preview() -> None:
    refresh_probe_access()
    probe_map()
    path = str(HW / "geometric_electricity.kicad_pcb")
    for arguments in [
        [
            "pcb",
            "export",
            "gerbers",
            "-o",
            str(OUT / "manufacturing-preview") + "/",
            path,
        ],
        [
            "pcb",
            "export",
            "drill",
            "-o",
            str(OUT / "manufacturing-preview") + "/",
            path,
        ],
        [
            "pcb",
            "export",
            "pos",
            "--format",
            "csv",
            "--units",
            "mm",
            "--smd-only",
            "--exclude-dnp",
            "-o",
            str(OUT / "positions.csv"),
            path,
        ],
        [
            "pcb",
            "export",
            "pdf",
            "--mode-single",
            "--crossout-DNP-footprints-on-fab-layers",
            "-l",
            "F.Fab,F.SilkS,Edge.Cuts",
            "-o",
            str(OUT / "assembly.pdf"),
            path,
        ],
        [
            "sch",
            "export",
            "pdf",
            "-o",
            str(OUT / "schematic.pdf"),
            str(HW / "geometric_electricity.kicad_sch"),
        ],
    ]:
        subprocess.run(["kicad-cli", *arguments], check=True, capture_output=True)
    (OUT / "manufacturing-preview" / "NOT_FOR_FABRICATION.txt").write_text(
        "Design preview. Check verification.json and landed budget; release is disabled. Do not order.\n"
    )


def phasors(
    path: Path, frequency: float, sense: float, probe_pf: float, fixture_pf: float = 0
) -> dict[str, Any]:
    """CSV columns time_s,vin_v,vsource_v,vscan_v; voltages already probe-scaled."""
    data = np.genfromtxt(path, delimiter=",", names=True)
    t = data["time_s"]
    if not np.all(np.diff(t) > 0) or frequency <= 0 or sense <= 0:
        raise ValueError("Increasing time and positive frequency/resistance required")
    if (
        len(t) < 20
        or (t[-1] - t[0]) * frequency < 5
        or np.max(np.diff(t)) * frequency > 0.1
    ):
        raise ValueError("Need >=5 cycles and >=10 samples/cycle")
    a = np.column_stack(
        [
            np.cos(2 * np.pi * frequency * t),
            np.sin(2 * np.pi * frequency * t),
            np.ones_like(t),
        ]
    )

    fit_quality = {}

    def phasor(column: str) -> complex:
        fitted = linalg.lstsq(a, data[column])[0]
        amplitude = float(np.hypot(fitted[0], fitted[1]))
        residual = float(np.sqrt(np.mean((data[column] - a @ fitted) ** 2)))
        fit_quality[column] = {
            "amplitude_v_peak": amplitude,
            "residual_v_rms": residual,
            "residual_to_sine_rms": (
                residual / (amplitude / np.sqrt(2)) if amplitude else None
            ),
        }
        return complex(fitted[0], -fitted[1])

    vin, vs, vm = [phasor(k) for k in ["vin_v", "vsource_v", "vscan_v"]]
    delivered = (vin - vs) / sense
    source_probe_y = 1 / 10e6 + 2j * np.pi * frequency * probe_pf * 1e-12
    if probe_pf < 0 or fixture_pf < 0:
        raise ValueError("Capacitances must be nonnegative")
    fixture_y = 2j * np.pi * frequency * fixture_pf * 1e-12
    graph_current = delivered - (source_probe_y + fixture_y) * vs
    if abs(graph_current) < 1e-12:
        raise ValueError("Source current too small for reliable reconstruction")
    z = vs / graph_current
    transfer = vm / graph_current
    return {
        "frequency_hz": frequency,
        "source_current_a": [graph_current.real, graph_current.imag],
        "source_impedance_ohm": [z.real, z.imag],
        "transfer_impedance_ohm": [transfer.real, transfer.imag],
        "fixture_pf": fixture_pf,
        "fit_quality": fit_quality,
        "note": "Source probe and specified lumped fixture capacitance deembedded; moving scan probe still loads the graph. Calibrate channel gain/phase first.",
    }


def fit_poles(
    path: Path, initial_mhz: list[float], output: Path, background_order: int = 2
) -> dict[str, Any]:
    """Variable-projection fit of common damped poles to calibrated complex responses."""
    data = np.load(path)
    hz = np.asarray(data["frequency_hz"], dtype=float)
    measured = np.asarray(data["impedance_ohm"], dtype=complex)
    if measured.ndim not in [1, 2]:
        raise ValueError("Impedance must have frequency and optional channel axes")
    if measured.ndim == 1:
        measured = measured[:, None]
    seeds = np.sort(np.asarray(initial_mhz, dtype=float))
    if (
        hz.ndim != 1
        or measured.shape[0] != len(hz)
        or len(seeds) == 0
        or len(hz)
        < max(10 * len(seeds), 2 * (2 * len(seeds) + 2 * (background_order + 1)))
        or not np.all(np.diff(hz) > 0)
        or np.any(hz <= 0)
        or not np.all(np.diff(seeds) > 0)
        or np.any(seeds <= 0)
        or not np.all(np.isfinite(measured))
    ):
        raise ValueError(
            "Need finite frequency-ordered complex responses and distinct positive pole seeds"
        )
    mhz = hz / 1e6
    if np.any(seeds < mhz.min()) or np.any(seeds > mhz.max()):
        raise ValueError("Pole seeds must be inside the measured band")
    count = len(seeds)
    if background_order not in range(5):
        raise ValueError("Background order must be between zero and four")
    background_center = float(mhz.mean())
    background_scale = float(np.ptp(mhz))
    local_frequency = (mhz - background_center) / background_scale
    polynomial = np.column_stack(
        [local_frequency**power for power in range(background_order + 1)]
    )
    if not np.any(abs(measured) > 0):
        raise ValueError("Response has zero amplitude")
    scale = np.maximum(
        np.sqrt(np.mean(abs(measured) ** 2, axis=0)), np.max(abs(measured)) * 1e-3
    )
    target = np.concatenate([measured.real, measured.imag], axis=0)

    def project(parameters: FloatArray) -> tuple[FloatArray, FloatArray]:
        f, gamma = parameters[:count], np.exp(parameters[count:])
        s = 1j * mhz[:, None]
        poles = -gamma + 1j * f
        positive = 1 / (s - poles)
        negative = 1 / (s - poles.conj())
        # Real impulse response enforces conjugate pole/residue pairs.
        basis = np.column_stack(
            [
                positive + negative,
                1j * (positive - negative),
                polynomial,
                1j * polynomial,
            ]
        )
        matrix = np.concatenate([basis.real, basis.imag], axis=0)
        coefficients = linalg.lstsq(matrix, target)[0]
        return matrix, coefficients

    def residual(parameters: FloatArray) -> FloatArray:
        matrix, coefficients = project(parameters)
        return ((matrix @ coefficients - target) / scale).ravel()

    span = mhz.max() - mhz.min()
    if span <= 2e-4:
        raise ValueError("Frequency window is too narrow for this fit")
    initial = np.concatenate(
        [seeds, np.log(np.clip(seeds / (2 * 36), 1.1e-4, span * 0.9))]
    )
    lower = np.concatenate([np.full(count, mhz.min()), np.full(count, np.log(1e-4))])
    upper = np.concatenate([np.full(count, mhz.max()), np.full(count, np.log(span))])
    fitted = least_squares(
        residual,
        initial,
        bounds=(lower, upper),
        x_scale="jac",
        max_nfev=500,
        ftol=1e-10,
        xtol=1e-10,
        gtol=1e-10,
    )
    matrix, coefficients = project(fitted.x)
    reconstructed = matrix @ coefficients
    reconstruction = reconstructed[: len(hz)] + 1j * reconstructed[len(hz) :]
    order = np.argsort(fitted.x[:count])
    frequencies = fitted.x[:count][order] * 1e6
    damping = np.exp(fitted.x[count:])[order] * 1e6
    residues = (coefficients[:count] + 1j * coefficients[count : 2 * count])[order] * (
        2 * np.pi * 1e6
    )
    np.savez_compressed(
        output,
        frequency_hz=frequencies,
        damping_hz=damping,
        residue_ohm_per_s=residues,
        reconstructed_impedance_ohm=reconstruction,
        response_frequency_hz=hz,
        background_coefficients=coefficients[
            2 * count : 2 * count + background_order + 1
        ]
        + 1j * coefficients[2 * count + background_order + 1 :],
        background_center_mhz=background_center,
        background_scale_mhz=background_scale,
    )
    report = {
        "converged": bool(fitted.success),
        "background_order": background_order,
        "frequency_hz": frequencies.tolist(),
        "damping_hz": damping.tolist(),
        "undamped_frequency_hz": np.hypot(frequencies, damping).tolist(),
        "Q": (frequencies / (2 * damping)).tolist(),
        "relative_residual_rms": float(
            np.linalg.norm(reconstruction - measured) / np.linalg.norm(measured)
        ),
        "jacobian_condition": float(np.linalg.cond(fitted.jac)),
        "note": "Conjugate pole pairs with a local complex polynomial background on the positive-frequency window; background is not a global passive network model. Calibrated common-pole data required; convergence does not prove completeness, identifiability or physical eigenvectors. Check multiple seeds, background orders and held-out frequencies.",
    }
    write_json(output.with_suffix(".json"), report)
    return report


def derived(path: Path, output: Path, lowest_complete: bool = False) -> dict[str, Any]:
    """Finite-graph functionals from fitted (not raw peak-map) eigenpairs."""
    data = np.load(path)
    q = np.asarray(data["q"], dtype=np.float64)
    v = np.asarray(data["vectors"], dtype=np.float64)
    if (
        q.ndim != 1
        or v.shape != (85, len(q))
        or not np.all(q > 0)
        or not np.all(np.isfinite(q))
        or not np.all(np.isfinite(v))
    ):
        raise ValueError("Need positive q[m] and real orthonormal vectors[85,m]")
    if not np.allclose(v.T @ v, np.eye(len(q)), atol=0.02):
        raise ValueError(
            "Eigenvectors must be normalized/orthogonal fitted modes, not unresolved peak maps"
        )
    t = np.logspace(-3, 2, 150)
    weights = np.exp(-np.outer(t, q))
    heat_kernel = np.einsum("ik,tk,jk->tij", v, weights, v)
    green = (v / q) @ v.T
    np.savez_compressed(
        output,
        t=t,
        heat_kernel=heat_kernel,
        heat_trace=weights.sum(axis=1),
        zeta_s=[1.0, 2.0, 3.0],
        zeta=[float(np.sum(q ** (-x))) for x in [1.0, 2.0, 3.0]],
        green_Q_inverse=green,
        participation_ratio=1 / np.sum(v**4, axis=0),
        complete_spectrum=len(q) == 85,
    )
    if lowest_complete and len(q) < 85:
        np.savez_compressed(
            output.with_name(output.stem + "_bounds.npz"),
            t=t,
            heat_trace_lower=weights.sum(axis=1),
            heat_trace_upper=weights.sum(axis=1) + (85 - len(q)) * np.exp(-t * q.max()),
            zeta_s=[1.0, 2.0, 3.0],
            zeta_tail_upper=[(85 - len(q)) * q.max() ** (-x) for x in [1.0, 2.0, 3.0]],
        )
    return {
        "modes": len(q),
        "nodes": 85,
        "complete_spectrum": len(q) == 85,
        "interpretation": "Finite graph; missing modes produce partial sums/kernels, not complete trace/zeta/Green operator.",
    }


def capture_scope(host: str, output: Path) -> None:
    """Capture one frozen acquisition via SCPI ASCII (no optional VISA package).

    Channels: CH1 before resistor, CH2 source, CH3 scanning node. The generator
    and front-panel scales are set manually. Scope must already be stopped.
    """
    with socket.create_connection((host, 5555), timeout=10) as connection:
        stream = connection.makefile("rb")

        def command(text: str) -> None:
            connection.sendall((text + "\n").encode("ascii"))

        def query(text: str) -> str:
            command(text)
            value = stream.readline()
            if not value:
                raise ConnectionError("Scope closed connection")
            return value.decode("ascii").strip()

        identity = query("*IDN?")
        if "DHO924" not in identity.upper():
            raise ValueError(f"Expected DHO924/S, got {identity}")
        if query(":TRIGger:STATus?").upper() != "STOP":
            raise ValueError(
                "Stop the scope acquisition first so all channels share one record"
            )
        waves = []
        times = []
        for channel in [1, 2, 3]:
            command(f":WAVeform:SOURce CHANnel{channel}")
            command(":WAVeform:MODE RAW")
            command(":WAVeform:FORMat ASCii")
            command(":WAVeform:STARt 1")
            command(":WAVeform:STOP 10000")
            preamble = np.array(
                [float(x) for x in query(":WAVeform:PREamble?").split(",")]
            )
            wave = np.fromstring(query(":WAVeform:DATA?"), sep=",")
            if len(wave) < 20:
                raise ValueError("Scope returned too few samples")
            times.append(
                (np.arange(len(wave)) - preamble[6]) * preamble[4] + preamble[5]
            )
            waves.append(wave)
        if not all(
            len(t) == len(times[0]) and np.allclose(t, times[0], rtol=0, atol=1e-15)
            for t in times
        ):
            raise ValueError("Channel time bases disagree")
        np.savetxt(
            output,
            np.column_stack([times[0], *waves]),
            delimiter=",",
            header="time_s,vin_v,vsource_v,vscan_v",
            comments="",
        )
        write_json(
            output.with_suffix(".json"),
            {
                "scope": identity,
                "channels": [1, 2, 3],
                "mode": "RAW ASCII",
                "points": len(waves[0]),
                "probe_scale": "front panel must be 10x",
            },
        )


def propagation(path: Path, geometry: str, output: Path) -> None:
    """Envelope and phase with true metric distances, for measured/synthetic scans."""
    data = np.load(path)
    t = data["time_s"]
    wave = data["voltage_v"]
    source = int(data["source"]) if "source" in data else 0
    if (
        wave.shape != (len(t), 85)
        or len(t) < 50
        or not np.allclose(np.diff(t), np.diff(t)[0], rtol=1e-4, atol=1e-15)
    ):
        raise ValueError("Need a uniform time_s vector and voltage_v[time,85]")
    analytic = signal.hilbert(wave - wave.mean(axis=0), axis=0)
    envelope = abs(analytic)
    phase = np.angle(analytic)
    valid = envelope > np.max(envelope) * 0.05
    g = next(g for g in graphs() if g["name"] == geometry)
    points = [complex(*xy) for xy in g["metric_xy"]]
    euclidean = np.array([abs(z - points[source]) for z in points])
    metric = (
        euclidean
        if geometry == "E"
        else np.array([hyperdistance(z, points[source]) for z in points])
    )
    np.savez_compressed(
        output,
        time_s=t,
        envelope_v=envelope,
        phase_rad=phase,
        valid=valid,
        metric_distance=metric,
        embedding_distance=euclidean,
        source=source,
        note="Phase is masked at low amplitude. Compare only before reflected waves; phase by itself does not establish a causal front.",
    )


def selfcheck(gs: list[dict[str, Any]]) -> None:
    for g in gs:
        cm = capacitance(g)
        assert np.allclose(cm.sum(axis=1), C * np.array(g["missing"]), atol=1e-22)
        assert np.linalg.eigvalsh(cm).min() > 0
        f, v = eigenpairs(cm, np.full(85, L))
        q = 1 / ((2 * np.pi * f) ** 2 * L * C)
        assert np.linalg.norm(cm @ v - C * v * q) < 1e-20
        # Exact symmetry permutation closes and preserves edges.
        angle = 2 * np.pi / g["coordination"]
        points = np.array([complex(*z) for z in g["metric_xy"]])
        perm = np.argmin(
            abs(points[:, None] * np.exp(1j * angle) - points[None, :]), axis=1
        )
        assert len(set(perm.tolist())) == 85
        adjacency = np.zeros((85, 85), dtype=int)
        for i, neighbors in enumerate(g["adjacency"]):
            adjacency[i, neighbors] = 1
        assert np.array_equal(adjacency, adjacency[perm][:, perm])
    # Analytic phasor recovery checks current orientation and source-probe subtraction.
    t = np.arange(0, 20e-6, 1e-9)
    frequency = 1e6
    z = 25 + 12j
    current = 0.1e-3 + 0j
    vs = current * z
    scan = current * (3 - 4j)
    vin = vs + 1000 * (current + (1 / 10e6 + 2j * np.pi * frequency * 10e-12) * vs)
    samples = np.column_stack(
        [t, *[np.real(p * np.exp(2j * np.pi * frequency * t)) for p in [vin, vs, scan]]]
    )
    path = OUT / "synthetic_capture.csv"
    np.savetxt(
        path,
        samples,
        delimiter=",",
        header="time_s,vin_v,vsource_v,vscan_v",
        comments="",
    )
    recovered = phasors(path, frequency, 1000, 10)
    assert np.allclose(recovered["source_impedance_ohm"], [25, 12], atol=1e-9)
    assert np.allclose(recovered["transfer_impedance_ohm"], [3, -4], atol=1e-9)

    # Independently synthesize a noisy overlapping pair; recover poles and damping.
    fit_hz = np.linspace(4.1e6, 4.9e6, 201)
    true_hz = np.array([4.455e6, 4.601e6])
    damping_hz = np.array([62000.0, 64000.0])
    pole = 2 * np.pi * (-damping_hz + 1j * true_hz)
    residue = np.array([[2e6 + 1e6j, 1e6 - 0.5e6j], [-1e6 + 3e6j, 2e6 + 0.4e6j]])
    argument = 2j * np.pi * fit_hz[:, None, None]
    impedance = 2 + np.sum(
        residue[None] / (argument - pole[None, :, None])
        + residue.conj()[None] / (argument - pole.conj()[None, :, None]),
        axis=1,
    )
    noise = np.random.default_rng(17)
    impedance += noise.normal(0, 0.02, impedance.shape) + 1j * noise.normal(
        0, 0.02, impedance.shape
    )
    fit_path = OUT / "polefit_synthetic.npz"
    np.savez(fit_path, frequency_hz=fit_hz, impedance_ohm=impedance)
    result = fit_poles(fit_path, [4.42, 4.63], OUT / "polefit_recovered.npz")
    assert result["converged"]
    assert np.max(abs(np.array(result["frequency_hz"]) / true_hz - 1)) < 0.001
    assert np.max(abs(np.array(result["damping_hz"]) / damping_hz - 1)) < 0.02


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        choices=[
            "build",
            "simulate",
            "test-simulations",
            "check",
            "phasor",
            "derive",
            "fit",
            "capture",
            "propagate",
            "route-export",
            "route-import",
        ],
    )
    parser.add_argument("--trials", type=int, default=200)
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--eigenpairs", type=Path)
    parser.add_argument("--responses", type=Path)
    parser.add_argument("--initial-mhz", type=float, nargs="+")
    parser.add_argument("--background-order", type=int, default=2)
    parser.add_argument("--waveforms", type=Path)
    parser.add_argument("--geometry", choices=["E", "H"], default="H")
    parser.add_argument("--output", type=Path, default=OUT / "analysis.npz")
    parser.add_argument("--host")
    parser.add_argument("--session", type=Path)
    parser.add_argument("--lowest-modes-complete", action="store_true")
    parser.add_argument("--frequency", type=float, default=1e6)
    parser.add_argument("--sense", type=float, default=1000)
    parser.add_argument("--probe-pf", type=float, default=10)
    parser.add_argument("--fixture-pf", type=float, default=0)
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    HW.mkdir(exist_ok=True)
    if args.command in ["route-export", "route-import"]:
        if args.command == "route-import" and args.session is None:
            parser.error("--session required")
        routing_exchange(args.session if args.command == "route-import" else None)
        return
    if args.command == "propagate":
        if args.waveforms is None:
            parser.error("--waveforms required")
        propagation(args.waveforms, args.geometry, args.output)
        return
    if args.command == "fit":
        if args.responses is None or args.initial_mhz is None:
            parser.error("--responses and --initial-mhz required")
        print(
            json.dumps(
                fit_poles(
                    args.responses, args.initial_mhz, args.output, args.background_order
                ),
                indent=2,
            )
        )
        return
    if args.command == "derive":
        if args.eigenpairs is None:
            parser.error("--eigenpairs required")
        print(
            json.dumps(
                derived(args.eigenpairs, args.output, args.lowest_modes_complete),
                indent=2,
            )
        )
        return
    if args.command == "capture":
        if args.host is None or args.csv is None:
            parser.error("--host and --csv required")
        capture_scope(args.host, args.csv)
        return
    if args.command == "phasor":
        if args.csv is None:
            parser.error("--csv required")
        print(
            json.dumps(
                phasors(
                    args.csv, args.frequency, args.sense, args.probe_pf, args.fixture_pf
                ),
                indent=2,
            )
        )
        return
    if args.trials < 1:
        parser.error("--trials must be positive")
    gs = graphs()
    selfcheck(gs)
    if args.command == "test-simulations":
        simulate(gs, args.trials)
        spice_check(gs)
        simulation_tests(gs, args.trials)
        return
    if args.command in ["build", "simulate"]:
        summary = simulate(gs, args.trials)
        spice_check(gs)
        print(
            json.dumps(
                {
                    k: {
                        key: value
                        for key, value in g.items()
                        if key in ["min_MHz", "max_MHz", "edges"]
                    }
                    for k, g in summary["geometries"].items()
                },
                indent=2,
            )
        )
    else:
        for g in gs:
            _, v = eigenpairs(capacitance(g), np.full(85, L))
            g["defects"] = select_defects(g, v)
    parts = components(gs)
    if args.command == "build":
        optimize_placement(parts)
        ground_xy = np.array(
            [
                p["xy"]
                for p in parts
                if p["ref"].startswith("TP") and p["value"] == "GND"
            ]
        )
        for g in gs:
            locations = {
                p["value"]: p["xy"] for p in parts if p["ref"].startswith("TP")
            }
            g["pcb_xy"] = [locations[f"{g['name']}{i:02d}"] for i in range(85)]
            references = {
                p["value"]: p["ref"] for p in parts if p["ref"].startswith("TP")
            }
            g["testpoint_references"] = [
                references[f"{g['name']}{i:02d}"] for i in range(85)
            ]
            distances, indices = cKDTree(ground_xy).query(np.array(g["pcb_xy"]))
            g["probe_ground_xy"] = ground_xy[indices].tolist()
            g["probe_ground_distance_mm"] = distances.tolist()
        write_json(HW / "graphs.json", gs)
        schematic(parts)
        board(parts)
        (HW / "geometric_electricity.kicad_pro").write_text(
            json.dumps(
                {
                    "meta": {
                        "filename": "geometric_electricity.kicad_pro",
                        "version": 1,
                    },
                    "board": {
                        "design_settings": {
                            "rules": {"min_clearance": 0.2, "min_track_width": 0.2}
                        }
                    },
                },
                indent=2,
            )
            + "\n"
        )
        with (HW / "bom.csv").open("w") as stream:
            writer = csv.writer(stream)
            writer.writerow(["Reference", "Value", "Footprint", "DNP", "Role", "MPN"])
            writer.writerows(
                [p[k] for k in ["ref", "value", "footprint", "dnp", "role"]]
                + [
                    (
                        "AIML-0805-1R0K-T"
                        if p["ref"].startswith("L")
                        else (
                            "CC0603JRNPO9BN331"
                            if p["ref"].startswith("C")
                            else "generic; see design plan"
                        )
                    )
                ]
                for p in parts
                if not p["ref"].startswith("TP")
            )
        export_preview()
    if args.command in ["build", "check"]:
        print(json.dumps(verify(gs, parts), indent=2))


if __name__ == "__main__":
    main()
