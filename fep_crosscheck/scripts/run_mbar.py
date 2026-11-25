#!/usr/bin/env python3
import argparse
import os
import numpy as np
from FastMBAR import FastMBAR

kB = 0.0019872043  # kcal/mol/K

def load_matrix(path_A):
    if os.path.isfile(path_A): return np.load(path_A)
    cand = [os.path.join(path_A, "A.npy")]
    for c in cand:
        if os.path.isfile(c): return np.load(c)
    raise FileNotFoundError(f"Could not find A matrix at {path_A}")

def load_Nk(path_Nk):
    if os.path.isfile(path_Nk): return np.load(path_Nk)
    cand = [os.path.join(path_Nk, "Nk.npy")]
    for c in cand:
        if os.path.isfile(c): return np.load(c)
    raise FileNotFoundError(f"Could not find Nk at {path_Nk}")

def main():
    ap = argparse.ArgumentParser(description="Run FastMBAR on pre-reduced energy matrix")
    ap.add_argument("--A", required=True, help="Path to A.npy")
    ap.add_argument("--Nk", required=True, help="Path to Nk.npy")
    ap.add_argument("--T", type=float, default=298.15, help="Temperature (K) for final conversion")
    ap.add_argument("--cuda", action="store_true", help="Use GPU")
    ap.add_argument("--bootstrap", action="store_true", help="Estimate uncertainties")
    
    # --- RESTORED ARGUMENTS ---
    ap.add_argument("--bootstrap_num_rep", type=int, default=100, help="Number of bootstrap iterations")
    ap.add_argument("--bootstrap_block_size", type=int, default=3, help="Block size for bootstrap")
    # --------------------------

    ap.add_argument("--out", default="mbar_results.npz", help="Output NPZ path")
    args = ap.parse_args()

    # 1. Load Data
    A = load_matrix(args.A)
    Nk = load_Nk(args.Nk).astype(int)
    M, N = A.shape

    print(f"Loaded A: {A.shape} (States x Frames)")
    print(f"Loaded Nk: sum={Nk.sum()}")

    # 2. Initialize FastMBAR
    # We assume A is ALREADY reduced (unitless kT) because thanh_compute_energies.py did (beta * U)
    print(f"Initializing FastMBAR (Bootstrap: {args.bootstrap}, Reps: {args.bootstrap_num_rep})...")
    
    fm = FastMBAR(
        energy=A,
        num_conf=Nk,
        cuda=args.cuda,
        bootstrap=args.bootstrap,
        bootstrap_num_rep=args.bootstrap_num_rep,     # <--- PASSED HERE
        bootstrap_block_size=args.bootstrap_block_size, # <--- PASSED HERE
        verbose=True
    )

    # 3. Extract Results (Unitless kT)
    DeltaF_kT = fm.DeltaF           # Matrix of all pairwise free energies
    DeltaF_std_kT = getattr(fm, "DeltaF_std", None)

    # 4. Convert Final Result to kcal/mol for Display
    kBT = kB * args.T
    
    total_ddG_kT = DeltaF_kT[M-1, 0] # Energy of Last State - First State
    total_ddG_kcal = total_ddG_kT * kBT
    
    if DeltaF_std_kT is not None:
        error_kcal = DeltaF_std_kT[M-1, 0] * kBT
    else:
        error_kcal = 0.0

    # 5. Print Summary
    print("\n" + "="*30)
    print(" FINAL FEP RESULTS ")
    print("="*30)
    print(f"Temperature: {args.T} K")
    print(f"Windows:     {M}")
    print(f"Total DDG:   {total_ddG_kcal:.4f} +/- {error_kcal:.4f} kcal/mol")
    print("="*30)

    # 6. Save Data
    np.savez(
        args.out,
        DeltaF_kT=DeltaF_kT,
        DeltaF_std_kT=DeltaF_std_kT,
        F_kT=fm.F,
        kBT=kBT, 
        total_ddG_kcal=total_ddG_kcal
    )
    
    # Save CSV in kcal/mol (User friendly)
    base_csv = os.path.splitext(args.out)[0] + "_DeltaF_kcal.csv"
    np.savetxt(base_csv, DeltaF_kT * kBT, delimiter=",", fmt="%.6f")
    print(f"Saved full matrix (kcal/mol) to: {base_csv}")

if __name__ == "__main__":
    main()
