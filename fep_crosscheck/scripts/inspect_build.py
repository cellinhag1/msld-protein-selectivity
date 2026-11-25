#!/usr/bin/env python3
"""Build-only inspection script: streams a system's own driver .inp (same
topology/patches/ligand-merge chain ALF itself uses) with NO minimization
and NO dynamics, then writes a plain PDB so the structure can be opened in
PyMOL/VMD before committing any GPU time to it.

Run with cwd == the system's own root directory (e.g. tad_6site/).
Usage: python3 inspect_build.py --driver tad.inp --out inspect_build.pdb
"""
import argparse
import pycharmm
from pycharmm.lingo import charmm_script
import pycharmm.write as write

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--driver", required=True, help="e.g. tad.inp (relative to prep/)")
parser.add_argument("--out", default="inspect_build.pdb")
args = parser.parse_args()

charmm_script("bomlev -2")
charmm_script(f"""
set nblocks = 2
set nsites = 1
set nsubs1 = 2
set nreps = 1
set temp = 298.15
set lams1s1 = 0.0
set lams1s2 = 0.0
set cs1s1s1s2 = 0.0
set xs1s1s1s2 = 0.0
set xs1s2s1s1 = 0.0
set ss1s1s1s2 = 0.0
set ss1s2s1s1 = 0.0
""")
charmm_script(f'stream "prep/{args.driver}"')
charmm_script("bomlev -1")

write.coor_pdb(args.out)
print(f"[OK] wrote {args.out}")
