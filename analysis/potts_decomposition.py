#!/usr/bin/env python
"""
Per-residue decomposition of a multi-mutation binding free energy, from a
Potts (direct-coupling) fit to multisite lambda-dynamics sampling.

Background
----------
Run MSLD with N mutation sites switching independently and you get one number:
the free energy of turning all N mutations on at once. That number does not say
which residue did the work.

A Potts model does. Fit the sampled joint occupancy to

    dG(S) = sum_i h_i(S_i) + sum_{i<j} J_ij(S_i, S_j)

and you recover a field h_i per site (that residue's own contribution) and a
coupling J_ij per pair (how much the residues depend on each other). The fit
needs only non-zero marginal and pairwise occupancy, not full joint sampling of
all 2^N states, so its cost grows quadratically rather than exponentially in
the number of sites.

Subtracting the apo leg from the complex leg gives the binding-relevant part:

    h_i(bind)  = h_i(complex)  - h_i(apo)
    J_ij(bind) = J_ij(complex) - J_ij(apo)

and the whole-panel binding free energy reconstructs as

    ddG_bind = sum_i h_i(bind) + 0.5 * sum_ij J_ij(bind)

which should agree with the directly computed value. That agreement is the
self-consistency check worth running before trusting any per-residue number.

Dropping the J terms gives the "additive" estimate, i.e. what you would get by
assuming the residues act independently. The gap between the two says how much
of the answer is cooperative.

Inputs
------
A directory per leg, each holding a Potts fit:

    <leg>/data/h.model.dat     physical fields, bias-corrected (kcal/mol)
    <leg>/data/J.model.dat     physical couplings (kcal/mol)
    <leg>/Result.txt           per-state free energies; first line is the
                               all-wild-type state, last is the all-mutant state

Use h.model.dat and J.model.dat, not h.bias.dat: the latter is the raw
maximum-likelihood fit with the sampling bias still in it.

Usage
-----
    python potts_decomposition.py --complex path/to/dca153_complex \
                                  --apo     path/to/dca153_apo \
                                  --labels  S663L,I778V,F787W,L804M \
                                  --bootstrap 2000
"""
import argparse
import os

import numpy as np

KT_298 = 0.001987 * 298.15  # kcal/mol


def read_model_fields(leg_dir: str, nsites: int):
    """Read the physical (bias-corrected) Potts fields and couplings.

    Each site contributes two rows/columns, wild-type then mutant. Only the
    mutant entries matter for a two-state-per-site panel.
    """
    h = np.loadtxt(os.path.join(leg_dir, "data", "h.model.dat"))
    J = np.loadtxt(os.path.join(leg_dir, "data", "J.model.dat"))
    mut_idx = [2 * k + 1 for k in range(nsites)]
    return h[mut_idx], J[np.ix_(mut_idx, mut_idx)]


def read_allmut_dG(leg_dir: str):
    """Free energy of the all-mutant state relative to all-wild-type.

    Result.txt lists one row per joint state: the site occupancies, then the
    free energy and its error in the last columns.
    """
    with open(os.path.join(leg_dir, "Result.txt")) as fh:
        lines = [ln for ln in fh.read().splitlines() if ln.strip()]
    first, last = lines[0].split(), lines[-1].split()
    if not all(x == "0" for x in first[:-3]):
        raise ValueError(f"{leg_dir}: first line is not the all-wild-type state")
    if not all(x == "1" for x in last[:-3]):
        raise ValueError(f"{leg_dir}: last line is not the all-mutant state")
    return float(last[-3]), float(last[-1])


def decompose(complex_dir: str, apo_dir: str, nsites: int):
    """Per-residue and pairwise binding contributions, plus consistency checks."""
    h_c, J_c = read_model_fields(complex_dir, nsites)
    h_a, J_a = read_model_fields(apo_dir, nsites)

    h_bind = h_c - h_a
    J_bind = J_c - J_a

    ddG_recon = h_bind.sum() + 0.5 * J_bind.sum()
    ddG_additive = h_bind.sum()

    dG_c, err_c = read_allmut_dG(complex_dir)
    dG_a, err_a = read_allmut_dG(apo_dir)
    ddG_true = dG_c - dG_a
    err_true = float(np.hypot(err_c, err_a))

    return dict(
        h_bind=h_bind,
        J_bind=J_bind,
        ddG_recon=ddG_recon,
        ddG_additive=ddG_additive,
        ddG_true=ddG_true,
        err_true=err_true,
        nsites=nsites,
    )


def subset_sum(J_bind, idx_a, idx_b):
    """Sum of couplings between two index sets.

    When the sets are the same, each unordered pair is counted once, matching
    the sum_{i<j} convention rather than the 0.5 * full-matrix form.
    """
    if idx_a is idx_b or list(idx_a) == list(idx_b):
        idx = list(idx_a)
        return float(sum(J_bind[idx[i], idx[j]]
                         for i in range(len(idx))
                         for j in range(i + 1, len(idx))))
    return float(sum(J_bind[a, b] for a in idx_a for b in idx_b))


def bootstrap_errors(complex_dir, apo_dir, nsites, n_boot=2000, seed=0):
    """Error bars on h_i(bind) by resampling the per-trial fits.

    Expects per-trial subdirectories (trial0, trial1, ...) alongside the pooled
    fit. Falls back to None when they are absent, in which case the point
    estimates stand without error bars rather than with invented ones.
    """
    rng = np.random.default_rng(seed)

    def trial_dirs(base):
        out = []
        for i in range(100):
            d = os.path.join(base, f"trial{i}")
            if os.path.isdir(os.path.join(d, "data")):
                out.append(d)
        return out

    c_trials, a_trials = trial_dirs(complex_dir), trial_dirs(apo_dir)
    if not c_trials or not a_trials:
        return None

    per_trial = []
    for cd, ad in zip(c_trials, a_trials):
        h_c, _ = read_model_fields(cd, nsites)
        h_a, _ = read_model_fields(ad, nsites)
        per_trial.append(h_c - h_a)
    per_trial = np.asarray(per_trial)

    n = len(per_trial)
    means = np.array([per_trial[rng.integers(0, n, n)].mean(axis=0)
                      for _ in range(n_boot)])
    return means.std(axis=0, ddof=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--complex", required=True, help="Potts fit directory, complex leg")
    ap.add_argument("--apo", required=True, help="Potts fit directory, apo leg")
    ap.add_argument("--labels", required=True,
                    help="comma-separated mutation labels in site order, e.g. S663L,I778V")
    ap.add_argument("--bootstrap", type=int, default=0,
                    help="bootstrap resamples for h_i error bars (needs per-trial fits)")
    args = ap.parse_args()

    labels = [s.strip() for s in args.labels.split(",") if s.strip()]
    res = decompose(args.complex, args.apo, len(labels))

    err = None
    if args.bootstrap:
        err = bootstrap_errors(args.complex, args.apo, len(labels), args.bootstrap)
        if err is None:
            print("No per-trial fits found; reporting point estimates without error bars.\n")

    print(f"{'residue':<10} {'h_i(bind)':>12}" + ("   +/-" if err is not None else ""))
    print("-" * (24 + (6 if err is not None else 0)))
    order = np.argsort(-np.abs(res["h_bind"]))
    for i in order:
        line = f"{labels[i]:<10} {res['h_bind'][i]:>+12.3f}"
        if err is not None:
            line += f"  {err[i]:.3f}"
        print(line)

    print()
    print(f"sum h_i(bind)  (additive, no coupling) = {res['ddG_additive']:+.3f} kcal/mol")
    print(f"h + J reconstruction                   = {res['ddG_recon']:+.3f} kcal/mol")
    print(f"directly computed ddG_bind             = {res['ddG_true']:+.3f} "
          f"+/- {res['err_true']:.3f} kcal/mol")

    gap = abs(res["ddG_recon"] - res["ddG_true"])
    if gap > 2 * res["err_true"] + 0.3:
        print(f"\nWARNING: reconstruction and direct value differ by {gap:.3f} kcal/mol, "
              "more than the error bar. Treat the per-residue split with caution.")
    else:
        print("\nReconstruction agrees with the direct value within error.")

    coop = res["ddG_recon"] - res["ddG_additive"]
    print(f"Cooperative part (coupling contribution) = {coop:+.3f} kcal/mol")


if __name__ == "__main__":
    main()
