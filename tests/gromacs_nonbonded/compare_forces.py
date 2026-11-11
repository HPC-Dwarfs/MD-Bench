#!/usr/bin/env python3
"""Compares two force dumps (lines "index fx fy fz", '#' comments).

Reports the errors relative to the RMS force magnitude of the reference, which is
the meaningful scale for atoms whose individual force is small by cancellation.
Exits with status 1 if --tol is given and the RMS error exceeds it.
"""
import argparse
import sys

import numpy as np


def load(path):
    return np.loadtxt(path, comments="#")[:, 1:4]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reference")
    ap.add_argument("test")
    ap.add_argument("--tol", type=float, default=None, help="max allowed rms error / rms reference force")
    ap.add_argument("--tol-max", type=float, default=None, help="max allowed max error / rms reference force")
    ap.add_argument("--label", default="")
    a = ap.parse_args()
    r, t = load(a.reference), load(a.test)
    if r.shape != t.shape:
        sys.exit(f"shape mismatch {r.shape} vs {t.shape}")
    scale = np.sqrt((r**2).sum(1).mean())
    d = np.linalg.norm(t - r, axis=1)
    rms, mx = np.sqrt((d**2).mean()) / scale, d.max() / scale
    net_r, net_t = r.sum(0), t.sum(0)
    print(
        f"{a.label:<34s} atoms={len(r)} rms|F|ref={scale:10.4f}  "
        f"rms_err={rms:.3e}  max_err={mx:.3e}  net_force_test={np.abs(net_t).max():.2e}"
    )
    bad = (a.tol is not None and rms > a.tol) or (a.tol_max is not None and mx > a.tol_max)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
