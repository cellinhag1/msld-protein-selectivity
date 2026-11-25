#!/usr/bin/env python3
"""Generic energy-matrix evaluator for the FFIX fixed-lambda-window scan.

Re-evaluates every frame of one window's DCD (produced by run_ffix.py)
under every other window's fixed-lambda Hamiltonian, building a reduced
(beta*U) energy matrix for FastMBAR -- same approach as
thanh_compute_energies_prot_277_tri2.py, generalized to any of the
apo_patch/apo_6site/apo_7site systems and with the frame range split into
--chunks pieces (same fix used there for "too much data in one job").

Run with cwd == the system's own root directory. Output:
  ffix/output/energy_matrix_w{rank}_c{chunk_id}.npy   shape (size, n_frames_in_chunk)
"""
import os
import argparse
import numpy as np

import pycharmm
from pycharmm.charmm_file import CharmmFile
from pycharmm.lingo import charmm_script, get_energy_value

kB = 0.0019872043  # kcal/mol/K

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--rank", type=int, required=True, help="window index whose DCD to score (matches ffix_w{rank}.dcd)")
parser.add_argument("--size", type=int, required=True, help="total number of lambda windows")
parser.add_argument("--chunks", type=int, default=5, help="split this window's frames into this many chunks")
parser.add_argument("--chunk_id", type=int, default=0, help="which chunk (0-based) this task processes")
args = parser.parse_args()

rank, size = args.rank, args.size
temp = 298.15
beta = 1.0 / (kB * temp)

rundir = "ffix"
os.makedirs(f"{rundir}/output", exist_ok=True)

filename = f"ffix_w{rank}"
out_file = CharmmFile(file_name=f"{rundir}/out/{filename}_energy_c{args.chunk_id}.out", file_unit=35,
                       formatted=True, read_only=False)
charmm_script("outu 35")

charmm_script("bomlev -2")
charmm_script(f"""
set nblocks = 2
set nsites = 1
set nsubs1 = 2
set nreps = 1
set temp = {temp}
set lams1s1 = 0.0
set lams1s2 = 0.0
set cs1s1s1s2 = 0.0
set xs1s1s1s2 = 0.0
set xs1s2s1s1 = 0.0
set ss1s1s1s2 = 0.0
set ss1s2s1s1 = 0.0
""")
charmm_script('stream "prep/apo.inp"')
charmm_script("bomlev -1")
charmm_script('stream "prep/nbond.str"')
charmm_script("faster on")

########################################
# OPEN THIS WINDOW'S DCD & COUNT FRAMES
########################################
dcd_path = f"{rundir}/dcd/{filename}.dcd"
if not os.path.exists(dcd_path):
    raise FileNotFoundError(dcd_path)

from mdtraj.formats import DCDTrajectoryFile
with DCDTrajectoryFile(dcd_path) as f:
    xyz, _, _ = f.read()
    n_frames = xyz.shape[0]
del xyz
print(f"[INFO] Detected {n_frames} frames in {os.path.basename(dcd_path)}")

chunks = args.chunks
chunk_id = args.chunk_id
frames_per_chunk = (n_frames + chunks - 1) // chunks
start = chunk_id * frames_per_chunk
stop = min(n_frames, start + frames_per_chunk)

if start >= n_frames:
    print(f"[INFO] chunk {chunk_id} is empty (start {start} >= {n_frames}). Exiting.")
    out_file.close()
    raise SystemExit(0)

dcd = CharmmFile(file_name=dcd_path, file_unit=101, formatted=False, read_only=True)

###############################################
# BUILD PERTURBATION (REDUCED) ENERGY MATRIX  #
###############################################
lambda_values = np.linspace(0.0, 1.0, size)
n_to_read = stop - start
energy_matrix = np.full((size, n_to_read), np.nan, dtype=np.float64)

charmm_script("blade on")

for i_ensemble, lam in enumerate(lambda_values):
    charmm_script(f"set ffixlam1 = {lam:.8f}")
    charmm_script(f"set ffixlam2 = {1.0 - lam:.8f}")
    charmm_script('stream "prep/aa_stream/sc_blocksetup_ffix.inp"')

    charmm_script("traj query unit 101")
    charmm_script("traj firstu 101 nunit 1 nocheck")
    for _ in range(start):
        charmm_script("traj read")

    for j in range(n_to_read):
        charmm_script("traj read")
        charmm_script("image byseg sele all end")
        if (j + 1) % 200 == 0:
            charmm_script("blade off")
            charmm_script("blade on")
        charmm_script("energy blade")
        U = get_energy_value("ENER")
        energy_matrix[i_ensemble, j] = beta * U

charmm_script("blade off")

np.save(f"{rundir}/output/energy_matrix_w{rank}_c{chunk_id}.npy", energy_matrix)
print(f"[OK] wrote {rundir}/output/energy_matrix_w{rank}_c{chunk_id}.npy  shape={energy_matrix.shape}")
print(f"[INFO] frames={n_frames}, chunks={chunks}, chunk_id={chunk_id}, range=[{start},{stop})")

out_file.close()
