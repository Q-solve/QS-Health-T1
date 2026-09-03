#!/usr/bin/env python3
"""Run a small Jordan-Wigner VQE on exported active-space integrals."""

import argparse
import importlib
import json
from pathlib import Path

import numpy as np
from qiskit.circuit.library import TwoLocal
from qiskit.primitives import StatevectorEstimator

VQE = importlib.import_module("qiskit_algorithms").VQE
COBYLA = importlib.import_module("qiskit_algorithms.optimizers").COBYLA
JordanWignerMapper = importlib.import_module(
    "qiskit_nature.second_q.mappers"
).JordanWignerMapper
FermionicOp = importlib.import_module("qiskit_nature.second_q.operators").FermionicOp


def fermionic_hamiltonian(h1, eri):
    n = h1.shape[0]
    terms = {}

    def add(label, value):
        if abs(value) > 1e-12:
            terms[label] = terms.get(label, 0.0) + complex(value)

    for p in range(n):
        for q in range(n):
            for spin in (0, 1):
                add(f"+_{2 * p + spin} -_{2 * q + spin}", h1[p, q])
    # Chemist's ERI: (pq|rs) -> a†p a†r a_s a_q.
    for p in range(n):
        for q in range(n):
            for r in range(n):
                for s in range(n):
                    value = 0.5 * eri[p, q, r, s]
                    for sp in (0, 1):
                        for sr in (0, 1):
                            add(
                                f"+_{2 * p + sp} +_{2 * r + sr} -_{2 * s + sr} -_{2 * q + sp}",
                                value,
                            )
    return FermionicOp(terms, num_spin_orbitals=2 * n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--maxiter", type=int, default=100)
    args = ap.parse_args()
    indir = Path(args.indir)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    h1 = np.load(indir / "active_h1.npy")
    eri = np.load(indir / "active_eri.npy")
    try:
        quantum_summary = json.loads((indir / "hamiltonian_summary.json").read_text())
        core_energy = float(quantum_summary.get("frozen_core_energy_hartree", 0.0))
    except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
        raise SystemExit(f"Invalid Hamiltonian summary: {exc}") from exc
    qubit_op = JordanWignerMapper().map(fermionic_hamiltonian(h1, eri))
    ansatz = TwoLocal(
        qubit_op.num_qubits, ["ry", "rz"], "cx", reps=2, entanglement="linear"
    )
    history = []

    def callback(eval_count, parameters, mean, metadata):
        try:
            energy = float(np.real(mean))
        except (TypeError, ValueError) as exc:
            raise RuntimeError("VQE returned a non-numeric energy") from exc
        try:
            iteration = int(eval_count)
        except (TypeError, ValueError) as exc:
            raise RuntimeError("VQE returned an invalid iteration") from exc
        history.append({"iteration": iteration, "energy_hartree": energy})

    vqe = VQE(
        StatevectorEstimator(),
        ansatz,
        COBYLA(maxiter=args.maxiter, tol=1e-8),
        callback=callback,
    )
    try:
        result = vqe.compute_minimum_eigenvalue(qubit_op)
        value = (
            result.electronic_energies[0]
            if hasattr(result, "electronic_energies")
            else result.eigenvalue
        )
        energy = float(np.real(value))
    except (RuntimeError, TypeError, ValueError) as exc:
        raise SystemExit(f"VQE failed: {exc}") from exc
    summary = {
        "method": "VQE",
        "mapping": "Jordan-Wigner",
        "ansatz": "TwoLocal(ry,rz,cx)",
        "optimizer": "COBYLA",
        "num_qubits": qubit_op.num_qubits,
        "electronic_energy_hartree": energy,
        "frozen_core_energy_hartree": core_energy,
        "energy_hartree": energy + core_energy,
        "iterations": len(history),
        "reference": "active-space Hamiltonian with frozen-core contribution",
    }
    (out / "vqe_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out / "vqe_convergence.csv").write_text(
        "iteration,energy_hartree\n"
        + "".join(f"{x['iteration']},{x['energy_hartree']:.12f}\n" for x in history)
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
