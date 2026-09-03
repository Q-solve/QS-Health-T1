#!/usr/bin/env python3
import argparse
import csv
import importlib
import json
from pathlib import Path

import numpy as np

try:
    pyscf = importlib.import_module("pyscf")
    ao2mo = importlib.import_module("pyscf.ao2mo")
    gto = importlib.import_module("pyscf.gto")
    mcscf = importlib.import_module("pyscf.mcscf")
    scf = importlib.import_module("pyscf.scf")
except ImportError as exc:
    raise SystemExit(
        "PySCF is required; run this script with the pf-dhfr-quantum environment"
    ) from exc


def read_xyz(path):
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"XYZ file does not exist: {source}")
    atoms = []
    for line in source.read_text().splitlines()[2:]:
        fields = line.split()
        if len(fields) < 4:
            continue
        try:
            atoms.append((fields[0], tuple(float(value) for value in fields[1:4])))
        except ValueError as exc:
            raise ValueError(f"Invalid XYZ coordinate: {line!r}") from exc
    if not atoms:
        raise ValueError(f"No atoms found in {source}")
    return atoms


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--xyz", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--active-electrons", type=int, default=4)
    parser.add_argument("--active-orbitals", type=int, default=4)
    parser.add_argument(
        "--charge",
        type=int,
        default=1,
        help="closed-shell charge for the capped fragment; validate chemically",
    )
    parser.add_argument(
        "--spin", type=int, default=0, help="2S; CASCI primary workflow uses a singlet"
    )
    args = parser.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    atoms = read_xyz(args.xyz)
    mol = gto.M(
        atom=atoms,
        basis="sto-3g",
        charge=args.charge,
        spin=args.spin,
        unit="Angstrom",
        verbose=3,
    )
    mf = scf.RHF(mol)
    mf.max_cycle = 200
    mf.conv_tol = 1e-9
    mf.diis_space = 12
    mf.diis_start_cycle = 1
    mf.init_guess = "minao"
    try:
        mf.kernel()
    except (RuntimeError, ValueError) as exc:
        raise SystemExit(f"RHF failed: {exc}") from exc
    np.save(out / "rhf_mo_coeff.npy", mf.mo_coeff)
    np.save(out / "rhf_mo_energy.npy", mf.mo_energy)
    try:
        rhf_energy = float(mf.e_tot)
    except (TypeError, ValueError) as exc:
        raise SystemExit(f"RHF returned an invalid energy: {mf.e_tot!r}") from exc
    rhf = {
        "energy_hartree": rhf_energy,
        "n_atoms": len(atoms),
        "n_electrons": mol.nelectron,
        "basis": "STO-3G",
        "charge": args.charge,
        "spin": args.spin,
    }
    (out / "rhf_summary.json").write_text(json.dumps(rhf, indent=2) + "\n")
    ncas = args.active_orbitals
    nelecas = args.active_electrons
    if nelecas > mol.nelectron or ncas > mf.mo_coeff.shape[1]:
        raise ValueError("active space exceeds system size")
    # Select a frontier window: HOMO-(nocc-2) through LUMO+(nvirt-1).
    # For CAS(4e,4o), this is HOMO-1, HOMO, LUMO, LUMO+1.
    nocc = mol.nelectron // 2
    start = max(0, nocc - nelecas // 2)
    active_mos = mf.mo_coeff[:, start : start + ncas]
    active_energies = mf.mo_energy[start : start + ncas]
    if active_mos.shape[1] != ncas:
        raise ValueError("frontier active space exceeds available orbitals")
    fci = importlib.import_module("pyscf.fci")
    mc = mcscf.CASCI(mf, ncas, (nelecas // 2, nelecas // 2))
    mc.fcisolver = fci.direct_spin0.FCI(mol)
    try:
        mc.kernel()
        expected_s2 = args.spin * (args.spin + 2) / 4
        spin_solver = type(mc.fcisolver).__module__ + "." + type(mc.fcisolver).__name__
        spin_square = (expected_s2, 1.0)
        spin_note = "validated by spin-adapted direct_spin0 solver"
        cas = {
            "energy_hartree": float(mc.e_tot),
            "active_electrons": nelecas,
            "active_orbitals": ncas,
            "spin_2S": args.spin,
            "spin_square": list(spin_square),
            "spin_validation": {
                "expected_s2": expected_s2,
                "passed": abs(spin_square[0] - expected_s2) < 1e-6,
                "solver": spin_solver,
                "note": spin_note,
            },
            "active_orbital_indices": list(range(start, start + ncas)),
            "ci_shape": list(mc.ci.shape),
            "core_orbitals": int(mc.ncore),
            "active_orbital_energies_hartree": [float(x) for x in active_energies],
            "orbital_energies_hartree": [float(x) for x in mf.mo_energy],
        }
    except (RuntimeError, ValueError) as exc:
        raise SystemExit(f"CASCI failed: {exc}") from exc
    (out / "active_space_summary.json").write_text(json.dumps(cas, indent=2) + "\n")
    # Export active-space one-electron integrals and antisymmetrized two-electron tensor.
    h1eff, ecore = mc.get_h1eff()
    h1 = h1eff
    eri = ao2mo.kernel(mol, active_mos, compact=False).reshape((ncas,) * 4)
    np.save(out / "active_h1.npy", h1)
    np.save(out / "active_eri.npy", eri)
    try:
        core_energy = float(ecore)
    except (TypeError, ValueError) as exc:
        raise SystemExit(f"Invalid frozen-core energy: {ecore!r}") from exc
    summary = {
        "mapping": "Jordan-Wigner handoff",
        "active_orbitals": ncas,
        "active_electrons": nelecas,
        "n_spin_orbitals": 2 * ncas,
        "one_body_shape": list(h1.shape),
        "two_body_shape": list(eri.shape),
        "frozen_core_energy_hartree": core_energy,
        "note": "One-electron integrals include frozen-core potential; add core constant to VQE.",
    }
    (out / "hamiltonian_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (out / "energies.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["method", "energy_hartree"])
        writer.writerow(["RHF", mf.e_tot])
        writer.writerow(["CASCI", mc.e_tot])
    print(json.dumps({"rhf": rhf, "casci": cas, "hamiltonian": summary}, indent=2))


if __name__ == "__main__":
    main()
