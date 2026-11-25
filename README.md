# msld-protein-selectivity

Multisite λ-dynamics (MSλD) applied to the **protein** rather than the ligand:
alchemically mutate a set of binding-pocket residues from one isoform into
another, and work out how much each individual residue contributes to the
difference in ligand binding.

The motivating question is selectivity. A drug that binds its intended target
will often bind a close relative of that target too. If you can mutate the
pocket of isoform A into the pocket of isoform B *in silico* and measure what
that does to binding, you can ask which residues a selective compound would
have to exploit.

## Two pieces

**`setup/`** — the alchemical machinery: how to make N pocket residues
switch identity during a simulation.

**`analysis/`** — the Potts decomposition: how to get from one whole-panel
number to a per-residue breakdown.

## Sidechain-only patching, without CATS

Standard MSλD protein-mutation setups duplicate the whole amino acid per
substituent, backbone included, and then add a `CATS` restraint to hold the
duplicated backbone copies on top of each other. It works, but it means every
site carries several copies of a backbone that was never supposed to change,
and the restraint is one more thing to get right.

The setup here patches **only the sidechain** onto the shared, real backbone.
The backbone is never duplicated, so it stays in the environment block and no
`cats sele` selection is ever invoked. Fewer atoms, one less restraint, and the
part of the residue that is genuinely alchemical is the only part that is
treated alchemically.

`setup/patches/sc2{l,v,w,m,f}.rtf` define one patch per *target* amino acid, so
a single file is reused across every site mutating to that residue. They are
named for the destination: `sc2l` builds a leucine sidechain.

The restriction: neither the wild-type nor the mutant residue may be glycine or
proline, since both change what the backbone itself is. Those sites need
whole-residue patches instead.

| File | Role |
|---|---|
| `setup/patches/sc2*.rtf` | sidechain patch residues, one per target amino acid |
| `setup/aa_stream/sc_patch.str` | loads all patch definitions |
| `setup/aa_stream/sc_patchloop.inp` | applies the right patch at each site |
| `setup/aa_stream/sc_selectloop.inp` | builds the atom selections per substituent |
| `setup/aa_stream/sc_blocksetup.inp` | writes the CHARMM `BLOCK` command |
| `setup/aa_stream/coupleloop.inp` | ties sites into one shared λ (see below) |
| `setup/alf_info.py` | site counts, temperature, bias settings for ALF |

## Coupled or independent λ

`alchemical_definitions.coupled.inp` puts every site on **one shared λ**: the
whole pocket converts from isoform A to isoform B together, like a concerted
chimera swap. One reaction coordinate, and the end states are the two real
proteins.

`alchemical_definitions.independent.inp` gives **each site its own λ**, so the
simulation samples mixed pockets as well as the two pure ones.

The coupled form is the cleaner thermodynamic question but samples poorly: one
λ has to drag every mutation across together, and the transition becomes rare.
The independent form samples far better and, because it visits mixed states, is
what makes the per-residue analysis below possible at all. It is the setup used
for the analysis here.

## From one number to per-residue contributions

Independent sites give a combinatorial problem: N sites with two states each
span 2^N joint states, and no simulation visits all of them for N beyond a
handful. A Potts (direct-coupling) model solves this by fitting

    dG(S) = sum_i h_i(S_i) + sum_{i<j} J_ij(S_i, S_j)

a field `h_i` per site and a coupling `J_ij` per pair. The fit needs only
non-zero marginal and pairwise occupancy rather than full joint sampling, so
the number of parameters grows quadratically in the sites, not exponentially.

Subtract the apo leg from the complex leg and the fields become
binding-specific:

    h_i(bind) = h_i(complex) - h_i(apo)

`h_i(bind)` is what each residue contributes to the change in binding free
energy, which is directly comparable to an experimental single-point
mutagenesis measurement on that residue.

```bash
python analysis/potts_decomposition.py \
    --complex path/to/complex_leg_fit \
    --apo     path/to/apo_leg_fit \
    --labels  S663L,I778V,F787W,L804M \
    --bootstrap 2000
```

Three checks it prints, which matter more than the ranking:

- **Reconstruction against the direct value.** `sum h + 0.5 sum J` should
  reproduce the independently computed whole-panel number. If it does not, the
  fit is not describing the sampling and the per-residue split is not
  trustworthy. The tool warns when the gap exceeds the error bar.
- **Additive against coupled.** `sum h_i` alone is the answer you would get
  assuming residues act independently. The difference from the full `h + J`
  value is the genuinely cooperative part. These can differ by more than a
  kcal/mol, and an additive estimate that happens to look closer to experiment
  is a coincidence, not a better method.
- **Bootstrap error bars**, resampled across independent trials. A site that
  transitions less often than its neighbours gets a wider bar, which is the
  correct consequence rather than a reason to discard it.

## Requirements

CHARMM or pyCHARMM with BLaDE for the simulations, plus
[ALF](https://github.com/RyanLeeHayes/ALF) for adaptive landscape flattening
and the Potts/DCA fit. The analysis script needs only NumPy and reads the
fit output from disk.

Set `CHARMM_EXE` to your CHARMM binary; nothing else is machine-specific.

## License

MIT, see `LICENSE`. No simulation results are included in this repository.
