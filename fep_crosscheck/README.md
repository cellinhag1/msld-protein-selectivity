# FEP cross-check

Classic fixed-λ-window free energy perturbation (FEP) with MBAR, set up from
the MSλD system in the parent directory so that the two methods can be
compared **like for like**.

This is an independent check on the MSλD numbers, not a replacement for them.

## Why

MSλD and conventional FEP answer the same question by different routes. MSλD
treats λ as a dynamical variable with an adaptively learned bias; FEP pins λ to
a fixed value in each of many independent windows and recombines them at the
end. When the two disagree, the useful question is whether the difference comes
from the method or from everything else that differs between two separately
built simulations: topology, box, patches, cutoffs, parameters.

This pipeline removes that ambiguity. It streams **the same prep files the
MSλD run already uses** and then overrides only the `BLOCK` command, throwing
out the learned bias and replacing it with a plain, zero-bias, hard-fixed λ.
Topology, dual-topology ligand or patched sidechain, box, and nonbonded setup
are all byte-identical to the MSλD run. The only thing that changes is how λ
is treated, so a discrepancy is attributable to the method itself.

## `ffix` is a constraint, not a bias

CHARMM's `ffix` keyword in the `BLOCK` command sets a per-block flag that
**suppresses integration of that block's λ equation of motion entirely**. λ
does not move.

This is a kinematic constraint and is categorically different from an
energetic bias, which lives in a separate force-constant/reference-value
structure and still lets λ evolve under a modified potential. Getting these
confused is easy, because both are configured through `BLOCK` and both are
loosely described as "restraining lambda". An earlier attempt here used a
harmonic restraint instead and had to chain windows through restart files to
stop them sticking at the boundaries; `ffix` needs none of that.

## Pipeline

1. **Set up and run the windows** — `scripts/run_ffix.py` streams the MSλD prep
   file, replaces `BLOCK` with the fixed-λ version, runs one window.
   `example_reference_block_setup/sc_blocksetup_ffix.inp` is a worked two-block
   template; for a ligand transformation the `Call` selections pick out the
   two end-state fragments instead of sidechain atoms, and the rest carries
   over unchanged.
2. **Extract energies** — `scripts/compute_energies_ffix.py` evaluates every
   saved frame at every λ, producing the cross-evaluation matrix MBAR needs.
3. **Gather** — `scripts/gather_ffix.py` concatenates the per-window,
   per-chunk matrices into `A.npy` (states × frames) and `Nk.npy` (frames per
   window).
4. **Estimate** — `scripts/run_mbar.py` wraps FastMBAR, returning the full
   pairwise free-energy matrix with optional bootstrap error bars.

`scripts/submit_ffix_chain.sh` chains stages 1–4 as Slurm jobs with
`--dependency=afterok`, so the whole scan runs unattended.

## Two documents worth reading

- **`METHOD_NOTES.txt`** — the method as a transferable recipe, with the
  `BLOCK` keywords spelled out and a list of what must be re-checked on a
  different cluster or a different system.
- **`BUILD_LOG.txt`** — the development log: why this approach, and every bug
  hit on the way to one clean window. It records the dead ends too, including
  the restrained-λ attempt this replaced.

## Known limits

- Built for **charge-neutral** transformations. A net-charge change needs
  PME/analytic charge corrections that this pipeline does not apply.
- Window count and length were chosen for a system whose barrier had already
  been located by a prior scan. They are not universal defaults: put windows
  where your own overlap is worst.
- Each window must start from a structure that is not already relaxed under a
  different λ's skewed potential. See `METHOD_NOTES.txt` on this; it is the
  easiest way to get a quietly wrong answer.
- The Slurm files carry cluster-specific module names and partitions. Set
  `SCRIPTDIR` and `SYSDIR`, and edit the module lines to match your site.

## Requirements

CHARMM or pyCHARMM, NumPy, MDTraj, and FastMBAR (1.4.6 tested, needs PyTorch),
in addition to the parent repository's requirements.
