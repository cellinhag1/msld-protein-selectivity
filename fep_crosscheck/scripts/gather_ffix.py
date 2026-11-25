#!/usr/bin/env python3
"""Assemble the per-window, per-chunk FFIX energy matrices into the final
A.npy / Nk.npy that run_mbar.py / FastMBAR consumes. Run from inside the
target system's own root dir (expects ffix/output/energy_matrix_w{w}_c{c}.npy).

Usage: python gather_ffix.py --nwindows 20 --nchunks 5
"""
import os
import argparse
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--nwindows", type=int, required=True)
parser.add_argument("--nchunks", type=int, required=True)
args = parser.parse_args()

rundir = "ffix"
M = args.nwindows

parts = []
Nk = np.zeros(M, dtype=int)

for w in range(M):
    chunk_parts = []
    for c in range(args.nchunks):
        fn = f"{rundir}/output/energy_matrix_w{w}_c{c}.npy"
        if not os.path.isfile(fn):
            raise FileNotFoundError(f"Missing chunk file for window {w}: {fn}")
        chunk_parts.append(np.load(fn))
    full = np.concatenate(chunk_parts, axis=1)  # (M, n_frames_w)
    if full.shape[0] != M:
        raise ValueError(f"window {w} matrix has {full.shape[0]} rows, expected {M}")
    np.save(f"{rundir}/output/energy_matrix_w{w}.npy", full)
    parts.append(full)
    Nk[w] = full.shape[1]

A = np.concatenate(parts, axis=1)  # (M, sum_w n_frames_w)

np.save(f"{rundir}/output/A.npy", A)
np.save(f"{rundir}/output/Nk.npy", Nk)

print("Saved:", f"{rundir}/output/A.npy", f"{rundir}/output/Nk.npy")
print("A.shape =", A.shape, " sum(Nk) =", Nk.sum())
