#!/usr/bin/env python3
"""Generic fixed-lambda-window (FFIX) FEP production driver, for the
coupled-lambda PDE5A->PDE6 apo systems (apo_patch / apo_6site / apo_7site).

Run with cwd == the system's own root directory (e.g. apo_patch/), so all
of that system's prep/ relative paths resolve exactly as they do for its
ALF flattening/production runs. Reuses that system's own prep/apo.inp
verbatim (topology, patches, psf/crd, crystal/image, nonbonded setup) --
the only thing this script does differently from ALF is override the
BLOCK command afterward with a plain, zero-bias, hard-fixed-lambda setup
(prep/aa_stream/sc_blocksetup_ffix.inp) instead of ALF's learned-bias one,
and hard-fix that window's lambda for the whole run (no ALF dynamics).

Each --rank gets lambda = linspace(0, 1, --size)[rank], run independently
(no restart-chaining between windows -- unlike the ALF flattening chain,
FFIX doesn't need it: every window is an independent equilibrium run at
one fixed composition, so there's no barrier to climb within a window).

Outputs land in ./ffix/{res,dcd,out,output}/.
"""
import os
import argparse
import numpy as np

import pycharmm
from pycharmm.charmm_file import CharmmFile
from pycharmm.lingo import charmm_script

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--rank", type=int, required=True, help="window index")
parser.add_argument("--size", type=int, required=True, help="total number of windows")
parser.add_argument("--nequil", type=int, default=150000, help="equilibration steps (default 300 ps @ dt=0.002)")
parser.add_argument("--nprod", type=int, default=5000000, help="production steps (default 10 ns @ dt=0.002)")
parser.add_argument("--nsavc", type=int, default=5000, help="steps between saved frames (default 10 ps/frame)")
parser.add_argument("--seed", type=int, default=20260, help="base random seed")
args = parser.parse_args()

rank, size = args.rank, args.size
lam1 = float(np.linspace(0.0, 1.0, size)[rank])
lam2 = 1.0 - lam1

temp = 298.15
dt = 0.002

rundir = "ffix"
for d in ("res", "dcd", "out", "output"):
    os.makedirs(f"{rundir}/{d}", exist_ok=True)

filename = f"ffix_w{rank}"

out_file = CharmmFile(file_name=f"{rundir}/out/{filename}.out", file_unit=35,
                       formatted=True, read_only=False)
charmm_script("outu 35")

charmm_script("bomlev -2")

# Dummy ALF-schema variables needed by prep/apo.inp -> sc_blocksetup.inp.
# Values never matter -- we override the BLOCK entirely right after.
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

# Reuse this system's own validated prep chain verbatim (topology, patches,
# psf/crd, crystal/image, and -- harmlessly -- ALF's own BLOCK, which we
# override below).
charmm_script('stream "prep/apo.inp"')

charmm_script("bomlev -1")
charmm_script('stream "prep/nbond.str"')
charmm_script("faster on")

# Same starting coordinates as every ALF flattening/production run.
charmm_script('read coor card resid name "prep/minimized.crd"')

# --- Override BLOCK with plain FFIX (hard-fixed, zero-bias) ---
charmm_script(f"set ffixlam1 = {lam1:.8f}")
charmm_script(f"set ffixlam2 = {lam2:.8f}")
charmm_script('stream "prep/aa_stream/sc_blocksetup_ffix.inp"')

charmm_script("nbonds nbxmod 5")

# prep/minimized.crd was minimized under ALF's initial, heavily-skewed
# lambda split (~0.999/0.001 per sc_blocksetup.inp's lambda0 formula), so
# the near-zero-weight substituent's atoms were never actually relaxed.
# Hard-fixing lambda at this window's value gives that substituent full
# weight instantly -- relax any resulting steric clash (CPU, pre-BLaDE)
# before running any dynamics, same as the validated ligand-FEP precedent.
charmm_script("mini sd nstep 250 nprint 50 step 0.005 tolenr 1e-3 tolgrd 1e-3")

charmm_script("shake fast bonh param")
charmm_script("scalar fbeta set 0.1 sele all end")

useblade = "prmc iprs 100 pref 1 prdv 100"

#################
# EQUILIBRATION #
#################
res_equil_w = CharmmFile(file_name=f"{rundir}/res/{filename}_equil.res", file_unit=22,
                          formatted=True, read_only=False)
charmm_script(f'open unit 26 write file name "{rundir}/res/{filename}_equil.lmd"')

charmm_script(f"""
calc pmass = ?natom * 0.12
dynamics cpt start -
  blade {useblade} -
  iseed {args.seed} {args.seed} {args.seed} {args.seed} -
  timestep {dt} -
  nstep {args.nequil} -
  nprint 1000 -
  iprfrq 1000 -
  nsavc 0 -
  nsavl -1 -
  isvfrq {args.nequil} -
  iunread -1 -
  iuncrd -1 -
  iunwri 22 -
  iunldm 26 -
  firstt {temp} -
  finalt {temp} -
  tstruc {temp} -
  tbath {temp} -
  pconstant -
  pmass @pmass -
  pref 1.0 -
  pgamma 20.0 -
  hoover -
  reft {temp} -
  tmass 1000 -
  ichecw 0 -
  ihtfrq 0 -
  ieqfrq 0 -
  iasors 1 -
  iasvel 1 -
  iscvel 0 -
  inbfrq 0 -
  ilbfrq 0 -
  imgfrq 0 -
  ntrfrq 0 -
  echeck -1
""")
res_equil_w.close()
charmm_script("close unit 26")

##############
# PRODUCTION #
##############
dcd_file = CharmmFile(file_name=f"{rundir}/dcd/{filename}.dcd", file_unit=21,
                       formatted=False, read_only=False)
res_equil_r = CharmmFile(file_name=f"{rundir}/res/{filename}_equil.res", file_unit=23,
                          formatted=True, read_only=True)
res_prod_w = CharmmFile(file_name=f"{rundir}/res/{filename}_prod.res", file_unit=24,
                         formatted=True, read_only=False)
charmm_script(f'open unit 25 write file name "{rundir}/res/{filename}_prod.lmd"')

charmm_script(f"""
calc pmass = ?natom * 0.12
dynamics cpt restart -
  blade {useblade} -
  iseed {args.seed + 1} {args.seed + 1} {args.seed + 1} {args.seed + 1} -
  timestep {dt} -
  nstep {args.nprod} -
  nprint 1000 -
  iprfrq 1000 -
  nsavc {args.nsavc} -
  nsavl -1 -
  isvfrq {args.nprod} -
  iunread 23 -
  iuncrd 21 -
  iunwri 24 -
  iunldm 25 -
  firstt {temp} -
  finalt {temp} -
  tstruc {temp} -
  tbath {temp} -
  pconstant -
  pmass @pmass -
  pref 1.0 -
  pgamma 20.0 -
  hoover -
  reft {temp} -
  tmass 1000 -
  ichecw 0 -
  ihtfrq 0 -
  ieqfrq 0 -
  iasors 1 -
  iasvel 1 -
  iscvel 0 -
  inbfrq 0 -
  ilbfrq 0 -
  imgfrq 0 -
  ntrfrq 0 -
  echeck -1
""")

res_equil_r.close()
res_prod_w.close()
dcd_file.close()
charmm_script("close unit 25")

charmm_script(f'write coor card name "{rundir}/output/final_w{rank}.crd"')
out_file.close()
