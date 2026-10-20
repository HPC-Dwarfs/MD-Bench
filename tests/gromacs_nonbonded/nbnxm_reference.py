#!/usr/bin/env python3
"""Exact (float64, O(N^2)) reference for the GROMACS NBNxM benchmark system.

Builds the same water system as gromacs/src/gromacs/nbnxm/benchmark (coordinates,
box, charges, LJ parameters, intramolecular exclusions) from MD-Bench's
src/common/gmxcoords.h and evaluates the non-bonded forces with the textbook
formulas, so both codes can be checked against an implementation that shares none
of their approximations.

Usage: nbnxm_reference.py --size N --coulomb {ewald,rf,none} [--cutoff 1.0] --out forces.txt
"""
import argparse
import re
import sys
from math import pi, sqrt
from pathlib import Path

import numpy as np
from scipy.special import erf as erf_, erfc as erfc_

# The GROMACS benchmark hands C6 and C12 of oxygen to the NBNxM atom data, which expects
# the pre-multiplied 6*C6 and 12*C12 of the force field. The benchmark therefore
# evaluates the force (C12 r^-12 - C6 r^-6) / r^2 with the numbers below, i.e. the LJ
# potential with C6/6 and C12/12. The reference reproduces what the benchmark computes.
C6_BENCH, C12_BENCH = 0.0026173456, 2.634129e-06
C6_O, C12_O = C6_BENCH / 6.0, C12_BENCH / 12.0
Q_O, Q_H = -0.8476, 0.4238
MOL = 3


def load_gmx_coords(header):
    text = Path(header).read_text()
    box = re.search(r"box1000\[\]\[3\]\s*=\s*\{(.*?)\};", text, re.S).group(1)
    box = np.array(re.findall(r"[-+0-9.eE]+", box), dtype=float).reshape(3, 3)
    body = re.search(r"coordinates1000\[\]\[3\]\s*=\s*\{(.*?)\};", text, re.S).group(1)
    xyz = np.array(re.findall(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?", body), dtype=float)
    return xyz.reshape(-1, 3), box


def build_system(header, size):
    x0, box0 = load_gmx_coords(header)
    if size & (size - 1):
        sys.exit("size must be a power of two")
    fac = [1, 1, 1]
    dim = 0
    f = size
    while f > 1:
        fac[dim] *= 2
        f //= 2
        dim = (dim + 1) % 3
    xs = []
    for ix in range(fac[0]):
        for iy in range(fac[1]):
            for iz in range(fac[2]):
                shift = np.array([ix * box0[0, 0], iy * box0[1, 1], iz * box0[2, 2]])
                xs.append(x0 + shift)
    x = np.concatenate(xs)
    box = np.array([fac[d] * box0[d, d] for d in range(3)])
    x = np.mod(x, box)
    n = len(x)
    typ = np.where(np.arange(n) % MOL == 0, 0, 1)
    q = np.where(typ == 0, Q_O, Q_H)
    return x, box, typ, q


def forces(x, box, typ, q, rc, coulomb, beta=None, epsfac=1.0):
    n = len(x)
    f = np.zeros((n, 3))
    mol = np.arange(n) // MOL
    isO = typ == 0
    k_rf2 = 1.0 / rc**3
    for i in range(n - 1):
        d = x[i] - x[i + 1 :]
        d -= box * np.round(d / box)
        r2 = (d * d).sum(1)
        m = r2 < rc * rc
        if not m.any():
            continue
        j = np.nonzero(m)[0] + i + 1
        d, r2 = d[m], r2[m]
        r = np.sqrt(r2)
        excl = mol[j] == mol[i]
        fs = np.zeros(len(j))
        if isO[i]:
            lj = isO[j] & ~excl
            rinv6 = r2[lj] ** -3
            fs[lj] += (12 * C12_O * rinv6 * rinv6 - 6 * C6_O * rinv6) / r2[lj]
        qq = epsfac * q[i] * q[j]
        if coulomb == "none":
            pass
        elif coulomb == "rf":
            fs += qq * (np.where(excl, 0.0, r**-3) - k_rf2)
        elif coulomb == "ewald":
            e = np.exp(-(beta * r) ** 2) * 2 * beta / sqrt(pi)
            fs += np.where(
                excl,
                qq * (-erf_(beta * r) / r**3 + e / r2),
                qq * (erfc_(beta * r) / r**3 + e / r2),
            )
        fv = fs[:, None] * d
        f[i] += fv.sum(0)
        np.subtract.at(f, j, fv)
    return f


def ewald_beta(rc, rtol=1e-5):
    beta, i = 5.0, 0
    while True:
        i += 1
        beta *= 2
        if erfc_(beta * rc) <= rtol:
            break
    lo, hi = 0.0, beta
    for _ in range(i + 60):
        beta = 0.5 * (lo + hi)
        if erfc_(beta * rc) > rtol:
            lo = beta
        else:
            hi = beta
    return beta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--header", default=str(Path(__file__).resolve().parents[2] / "src/common/gmxcoords.h"))
    ap.add_argument("--size", type=int, default=1)
    ap.add_argument("--coulomb", choices=["ewald", "rf", "none"], default="ewald")
    ap.add_argument("--cutoff", type=float, default=1.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    x, box, typ, q = build_system(a.header, a.size)
    beta = ewald_beta(a.cutoff) if a.coulomb == "ewald" else None
    f = forces(x, box, typ, q, a.cutoff, a.coulomb, beta)
    with open(a.out, "w") as fp:
        fp.write(f"# natoms {len(x)}\n")
        for i, v in enumerate(f):
            fp.write(f"{i} {v[0]:.9g} {v[1]:.9g} {v[2]:.9g}\n")


if __name__ == "__main__":
    main()
